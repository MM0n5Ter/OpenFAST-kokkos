# circulation-solve (OLAF unit 4a; 4b `wake-wind` not in this card's scope)

Unit 4a of the OpenFAST-kokkos port, chosen by the host-time profile in `notes/olaf-unit4.md`: the lifting-line circulation iteration of OLAF with the geometry hoisted out of the `do while` (F10-B, class **B**), in the Kokkos **host** execution space in every build, CUDA included. 4b (wind and tower at the OLAF wind points, AeroDyn) is a separate later unit; its contract and audit rows W1-W6 are kept below unchanged.

- **Owner**: agent session `199ad883` (subagent, unit 4a port). Date: 2026-11-17.
- **Toolchain**: Kokkos 5.3.0 (C++20), gcc/gfortran 12.3.0, cmake 3.31.12, CUDA 12.3, OpenBLAS 0.3.28.

## Status (2026-11-17)

**Done.** Every gate in scope passes on CPU, and the CUDA build was checked on an A100 MIG slice. Nothing is committed and no ledger row was written (the parent does both).
- **FMA experiment → F1** (fast mode with FMA, consistent with architecture §3), by the rule fixed in advance. F0 (no FMA) stays as `OFK_CS_FMA=0`: bitwise with strict and with the pin, including every iteration count.
- **Speed, CPU, P2 config (8 threads, med Xeon 8480+):**
  - solve kernel 79.5 ms → 0.84 ms per solve (pin's `ui_quad_n1` time vs build + replays);
  - whole solve (wsamp) 80 → 1.5 ms;
  - `bar-r8-tree` 44.0 → 17.0 s (**2.59×**), `bar-r8-direct` 128.2 → 50.8 s (**2.52×**);
  - 1 thread: `bar-r8-tree` 98.0 → 48.6 s.
- **GPU (A100 MIG, ov4-cuda-a100):** `bar-r8-tree` 34.2 s with `OFK_CS=fortran` → 16.1 s.
- **Found and fixed by R6:** `ad_Sphere_OLAF` has no lifting-line control point (nCP_tot = 0), and the first build aborted there. The shim now keeps the pin's empty loop. The trees `ov4-*` predate the fix and `ov4b-*` have it; BAR results are bitwise the same on both (job 1926189).

## Files (uncommitted)

| File | What |
|---|---|
| `src/modules/aerodyn/src/FVW_Wings.f90` | **New overlay copy** of the pin's file. +96 lines, all under `#ifdef OPENFAST_KOKKOS`: module-level `bind(C)` interfaces (private), the private `ofk_cs_setup`, and the shim in `Wings_ComputeCirculationPolarData`: `ofk_cs_setup` once before the `do while` (skipped when nCP_tot = 0 or `OFK_CS=fortran`), `ofk_cs_vvar` instead of the `Vvar` block (`:424-447`, kept in the `else` branch), `ofk_cs_end` after the loop. Preprocessed without `OPENFAST_KOKKOS` it is the pin's file (A2). |
| `src/modules/aerodyn/src/FVW_Subs.f90` | +4 lines under `#ifdef`: `FVW_Kokkos_Release` (called by `FVW_End`) also frees this unit's handle (`ofk_cs_free`, same key `m%Sgmt`). |
| `src/modules/kokkos/src/circulation_common.hpp` | `Geom<T>` (raw pointers into the handle), constants, per-mode entry declarations. No arithmetic. |
| `src/modules/kokkos/src/circulation_body.hpp` | **The body**: `seg_geom` (ui_seg_11 up to its Γ statement, per segment over all targets, SIMD over targets) and `replay_block` (the Γ statement and the accumulation, per block of 8 targets, sources in the pin's order); Kokkos host `RangePolicy` or the same work items in a plain loop. |
| `src/modules/kokkos/src/circulation_{strict,fast,fastfma}.cpp` | The three mode TUs (5 lines each). |
| `src/modules/kokkos/src/circulation.{hpp,cpp}` | C API, handle registry, mode / FMA / serial-parallel selection, `OFK_CS_ITERS` iteration log, `OFK_CS_PROFILE` counters. |
| `src/modules/kokkos/ofk_kernels.cmake`, `CMakeLists.txt` | The four new sources; per-file flags of the three mode TUs (independent of `OFK_FAST_FMA`). |
| `tests/circulation-solve/` | `CMakeLists.txt` (oracle `_ref` = pin, `_rfm` = pin with FMA, port), `l0_cs.f90` (L0), `afi_count.c` (iteration counts via `ld --wrap`), `perf_cs.cpp` (P4), `r7_threads.sh`, `itertrace.c` (ptrace iteration counter for unmodified binaries), `iters_cmp.py`, `mutate.py`, `a2.sh`, `cpu.sbatch` (`iters`, `timing`, `rtest`), `fmatime.sbatch`, `p4.sbatch`, `gpu_mig.sbatch`, `final.sbatch`, this card. |

E1 (size): Fortran in scope 26 lines (`FVW_Wings.f90:424-447`) plus `ui_quad_n1`/`ui_seg_11` (≈120 lines of `FVW_BiotSavart.f90` whose arithmetic is reproduced); C++ written 826 lines (`circulation*`); Fortran shim 96 + 4 lines; tests and job scripts 1 559 lines.

## Contract (frozen, pin `12fa3fc`)

### 4a: circulation solve

```fortran
subroutine Wings_ComputeCirculationPolarData(z, z_prev, p, x, m, AFInfo, GammaScale, ErrStat, ErrMsg, iLabel)   ! FVW_Wings.f90:314
   type(FVW_ConstraintStateType), intent(inout) :: z        ! OUT  z%W(iW)%Gamma_LL(1:nSpan)                         (:488)
   type(FVW_ConstraintStateType), intent(in)    :: z_prev   ! IN   initial guess (unless m%FirstCall)                 (:366-377)
   type(FVW_ParameterType),       intent(in)    :: p        ! nWings, W(:)%nSpan, iNWStart (=2 by default, FVW.f90:103-104),
                                                            ! RegFunction, CircSolvMaxIter/ConvCrit/Relaxation, InductionAtCP
   type(FVW_ContinuousStateType), intent(in)    :: x        ! READ x%W(iW)%r_NW(1:3, 1:nSpan+1, 1:iNWStart+1)          (:430-433)
                                                            ! READ x%W(iW)%Eps_NW(1, 1:nSpan, 1:iNWStart): ROW 1 ONLY (:441)
   type(FVW_MiscVarType),         intent(inout) :: m        ! READ  W(:)%CP, Vwnd_CP, Vstr_CP, Norm, Tang, dl, Area, FirstCall
                                                            ! WRITE W(:)%Vtot_CP (per iteration, :453), alpha_LL, Vreln_LL (:582-583),
                                                            !       Vind_CP = Vtot_CP = -9999 on exit (:506-507); LL #b writes Sgmt/Part/Vind_*
   real(ReKi),                    intent(in)    :: GammaScale
   type(AFI_ParameterType),       intent(in)    :: AFInfo(:)
   integer(IntKi),                intent(in)    :: iLabel
```
- **One caller:** `Wings_ComputeCirculation` (`FVW_Wings.f90:295`), reached from `FVW_CalcConstrStateResidual` (`FVW.f90:1390`), which `FVW_UpdateStates` calls once per call (`FVW.f90:727`).
- **Calls per step:** 1 per `AD_UpdateStates`. In the OpenFAST glue that is 1 + `NumCrctn` (`FAST_Solver.f90:1255-1275`). Re-evaluation at the same `t` re-runs it (`FVW.f90:597-615`).
- **Measured** on `bar-r8-*`: 25.7 iterations per solve on average (7720 iterations over 300 solves, gprof job 1923304); the cap is `CircSolvMaxIter` = 30 (`FVW_IO.f90:47`).
- **What is shimmed:** only the `Vvar` block of the iteration (`FVW_Wings.f90:422-447`). It becomes:
  - before the `do while` (after `:399`): `ofk_cs_build(...)` builds the Γ-independent geometry cache, once per call;
  - in each iteration: `ofk_cs_vvar(GammaLastIter, Vvar)`, which replays the Γ statement of `ui_seg_11` (`FVW_BiotSavart.f90:125-126`) in the pin's order.
  - These stay Fortran, unchanged, because they decide the iteration count (audit C4):
    - `CirculationFromPolarData` (AirfoilInfo, L4);
    - relaxation, `MeanGamma` and the convergence test (`:463-468`);
    - the non-converged branch (`:471-480`);
    - the deflattening (`:483-490`);
    - the `-9999` safety writes.
- **`ui_quad_n1`** (`FVW_BiotSavart.f90:513`) keeps its interface and body: its only caller is bypassed under `OPENFAST_KOKKOS`.
- **Handle (as built):** the cache lives in a per-instance handle keyed by the address of `m%Sgmt` (as unit 2), in a mutex-guarded process-wide map with LRU cap (`OFK_CS_MAX_HANDLES`, 64) and a `Kokkos::push_finalize_hook`; freed by `FVW_End` through `FVW_Kokkos_Release` (`ofk_cs_free`). 3 reals per (segment, target) plus the ring corners: **1.49 MB for nCP = 87** (7 Kokkos HostSpace Views, all allocated on the first solve). It is **rebuilt from the Fortran arguments on every call and never trusted across calls**, so there is no residency and A7 holds by construction.
- **Modes (as built):** one body, three translation units: strict (baseline ISA, no contraction), fast F0 (x86-64-v3, no contraction) and fast F1 (x86-64-v3, contraction). The per-target order is the pin's in all three (ring-major, depth, segments 1-2, 3-4, 2-3, 4-1). The FMA policy of the fast mode was decided by the FMA experiment below (the draft card assumed "no FMA in either TU").

### 4b: host wind work at the OLAF wind points (AeroDyn, FVW pack)

> **InflowWind is an external library call** (user decision, 2026-11-17; decisions #14). 4b only changes how often and where AeroDyn calls it. **H1 and W1 are out of scope:** they are recorded only as an upstream suggestion for InflowWind's owners.

| Item | Routine (pin) | Change | Class | Interface |
|---|---|---|---|---|
| ~~H1~~ (upstream only) | `UniformField_GetVel` (`IfW_FlowField.f90:385-410`, `pure`, module-private) | Drop `use, intrinsic :: ieee_arithmetic` (`:386`): nothing from the module is used, and gfortran then saves and restores the FP environment on every call | B (bitwise: `scripts/unit4-mb` `ieeemb`, 0 differing of 5.6M points; 14.4× faster per point) | unchanged |
| H2 | `AD_UpdateStates` (`AeroDyn.f90:1942`), the two `AD_CalcWind` calls on `uInterp` (`:1996`) | Skip the `InflowWakeVel` part (`AD_CalcWind` `:2094-2121`) when `Wake_Mod = FVW`: after the loop, `m%Inflow(1)%InflowWakeVel` is the value from `:1977` (same `utimes(1)`, same `o%WakeLocationPoints`). Keep the rotor and GS parts | B | `AD_CalcWind` is private: add a private variant or an optional argument; `AD_CalcWind_Rotor`/`_GS` (public, `:65-66`) are unchanged |
| H3 | `SetInputsForFVW` from `AD_UpdateStates` (`AeroDyn.f90:2027` → `:4226-4235` `TwrInflArray` over **all** `r_wind` points, existing or not) | In `AD_UpdateStates`, only these slices of `V_wind` are ever read:<br>- `u(1)`: LL, NW, FW (`FVW.f90:645`, `:1378` through `uInterp`), plus source panels (`FVW_Subs.f90:2105`);<br>- `u(2)`: LL and panels (through `uInterp`, `FVW.f90:1378`).<br>The grid slice is read only from `AD_CalcOutput`'s own `SetInputsForFVW` (`FVW.f90:1618`). So apply `TwrInflArray` only to the read slices for these two calls. The `AD_CalcOutput` call (`:2378`) is unchanged | B | private (`SetInputsForFVW` is not public) |
| H4 | `SetRequestedWindPoints` (`FVW_Subs.f90:568`, grid block `:606-620`) | Write grid points only when `t+DTaero` is inside a grid's `[tStart − DTaero, tEnd + DTaero]` window, a superset of `bDoGrid` (`FVW.f90:1621-1625`); otherwise `(0,0,0)`, which `UniformField` skips (`z > 0` test, `IfW_FlowField.f90:182`; Grid3D `:229`). This removes 9600 of about 20.6k points from every wind call on the BAR cases. The tower loop still visits them: that is H3's job | B for every output and VTK file. `m%GridOutputs%uGrid` changes at non-output steps (misc, never written there) | unchanged |
| (later) | `AD_CalcOutput` wind + tower on the wake slice (`AeroDyn.f90:2365`, `:2378`) | Consumed only by grid outputs: needs an FVW query (deferred: Lane D / InflowWind on device) | B | — |

Callers that must keep working: `AD_UpdateStates`/`AD_CalcOutput` from `ADI_*` (driver), the OpenFAST glue (`FAST_Funcs.f90:459-466`), the linearization perturbations (`AeroDyn.f90:7629-7656`: they call `AD_CalcWind_Rotor`/`_GS` and `SetInputsForFVW` directly, so H3 must key on `tIndx`, not on the caller), and FAST.Farm/AWAE (`IfW_FlowField_GetVelAcc` has callers in 4 module directories: H1 speeds them all up, with bitwise-identical results).

## Audit list → L0 case

## Design

### What runs where
```
Wings_ComputeCirculationPolarData (Fortran, unchanged except the Vvar block)
  first-call guess, LL #b (unit 2), Vcst                     Fortran / unit 2, unchanged
  ofk_cs_setup -> ofk_cs_begin / ofk_cs_wing x nWings / ofk_cs_build      once per call (nCP_tot > 0)
  do while (not converged and iIter < CircSolvMaxIter)
     ofk_cs_vvar(GammaLastIter -> Vvar)                      replaces FVW_Wings.f90:424-447
     Vtot_CP, CirculationFromPolarData, relaxation, MeanGamma, convergence test   Fortran, unchanged (C4)
  ofk_cs_end(iIter)                                          counters / OFK_CS_ITERS log only
```
- **Cache (built once per call).** For every source segment j = 4·r + k (rings r in the pin's (iW, iSpan, iDepth) order; k = 1-2, 3-4, 2-3, 4-1) and every target t (control points flattened over wings, padded to a multiple of 8 with copies of target 0): `Kv` after the regularization select, `s = norm_a+norm_b` and `den = denominator` (with the Offset term). A segment that `ui_seg_11` leaves early (its three `<= PRECISION_UI` guards; NaN does not leave) is marked by `s = -1`: a computed `norm_a+norm_b` is never negative.
- **Replay (per iteration).** Per target, in the pin's order: `K = ((Γ·fourpi_inv)·Kv)·s / den` (left to right, as `:125`), the cross product recomputed from the same coordinates by the same expressions (same bits; recomputing halved the memory traffic and was faster than storing it), `V += marked ? 0 : K·c`. Adding +0 for marked segments is the pin's `Uind + 0`: `Vvar` starts at +0 and a round-to-nearest sum from +0 never becomes -0.
- **Parallel shape.** Build: `RangePolicy` over the 696 source segments, SIMD over targets. Replay: `RangePolicy` over blocks of 8 targets (SIMD lanes), each block summing every source in order, so the per-target order, and the result, do not depend on the thread count (R7). Both in `Kokkos::DefaultHostExecutionSpace` (OpenMP in every build here, CUDA included).
- **P4 serial fallback.** Below `OFK_CS_PAR_MIN` = 2048 (replay) and `OFK_CS_PAR_MIN_BUILD` = 512 (build) target × segment pairs the same work items run in a plain loop. Thresholds = the measured crossovers at 8 threads (job 1925426, below).
- **Variants.** `OFK_CS_MODE` (default: `OFK_MODE`) selects this unit's mode; in fast mode `OFK_CS_FMA` selects F0 (no contraction) or F1 (contraction). `OFK_CS=fortran` keeps the pin's loop (A/B). The default is set by the FMA experiment (`kDefaultFma` in `circulation.cpp`).
- **Constants** are passed from the shim (R3): `PRECISION_UI`, `fourpi_inv`, `MIN_EXP_VALUE` from `FVW_BiotSavart`, and the literal `1.25643_ReKi`.
- **nCP_tot = 0** (source panels only, `ad_Sphere_OLAF`): the shim keeps the pin's (empty) loop. Found by R6 (the first build aborted there); L0 now has the case.

## Audit list → L0 case
L0 = `tests/circulation-solve/l0_cs.f90`. Every case runs the pin (twice: ref vs ref), the pin with FMA, and six port variants (strict with the threshold / serial / parallel_for, F0 serial / parallel_for, F1) on copies of the same misc; values compared bitwise (NaN = NaN) and iteration counts from the AirfoilInfo call count.

| # | Trap | Where (pin) | What the port does | L0 case |
|---|---|---|---|---|
| C1 | `ui_quad_n1`'s directives are `!OMP`, not `!$OMP`: the pin's loop is **serial** | `FVW_BiotSavart.f90:527-528, 548-549` | P2 compares against the serial pin; the oracle is the pin unchanged | oracle `_ref` built from the pin's files |
| C2 | Γ enters `ui_seg_11` first in a left-associated product | `FVW_BiotSavart.f90:125` | `((G·fourpi_inv)·Kv)·s/den` | every case, all RegFunction 0-4 (`branches`); mutation "fourpi_inv·Kv first" caught |
| C3 | Early returns leave `Uindtmp = 0`, which the pin still adds | `:80, :85, :89` → `:534-546` | marked by `s = -1`, `V += 0` | `degenerate` 1-3, 5, 7 (CP on a ring corner, zero-length span segment, collapsed ring, NaN in `r_NW`, CP on the bound segments), `onLL` cases; mutation "guards dropped" caught |
| C4 | Discrete decisions amplify rounding | `FVW_Wings.f90:466-471` | class B in strict and F0 | iteration count compared in every solve; `knife_cases`: the two adjacent doubles of ConvCrit where the pin's count changes (bisection), 8 configurations; `iteration cap reached` |
| C5 | Ring geometry and segment order | `FVW_Wings.f90:432-435`, `FVW_BiotSavart.f90:532-546` | literal | every case; mutation "3-4 as 4-3" caught |
| C6 | FMA changes the cache and, through C4, possibly the count | `FVW_BiotSavart.f90:76-118` | strict and F0 without contraction; F1 measured | F1 and pin-with-FMA statistics in every run; the FMA experiment |
| C7 | RegParam = `Eps_NW(1, iSpan, iDepth)` for all four segments | `FVW_Wings.f90:441` | row 1 | `degenerate` 6 (rows 2-3 NaN); mutation "row 2" caught |
| C8 | Targets flattened over wings; nSpan may differ, may be 0 | `FVW_Wings.f90:437-443` | per-wing offsets | `wings` (1-3 wings, nSpan 1..13), random nWings 1-3 / nSpan 1-12, `a wing without control points`, `no control point (nCP_tot = 0)` |
| C9 | `Vind_CP`/`Vtot_CP` = -9999 on exit; `alpha_LL`/`Vreln_LL` from the last polar call | `FVW_Wings.f90:453, 506-507, 582-583` | Fortran, unchanged | compared in every case (also `Vind_LL` from LL #b) |
| C10 | FirstCall guess; GammaScale 0 / not 0 | `FVW_Wings.f90:357-377` | Fortran, unchanged | `FirstCall guess`, `GammaScale = 0`, `= 0.5`, random |
| C11 | P6: O(nCP²) cache | — | handle, sized on the first solve | `p6_sequence`: 20 solves, allocations after solve 1 = after solve 20 (7, 1.49 MB) |
| C12 | Re-evaluation, copies (A7) | `FVW.f90:597-615` | nothing trusted across calls | `a7_sequence`: re-evaluation, misc copy (new handle) in between, Eps only / r_NW only / CP changed, orphaned copy; strict and F0; mutation "cache reused" caught (only by A7) |
| new | The pin is not reproducible when a wing has no control point (uninitialised reads) | LL #b / `Vind_LL` | — | ref vs ref check per case; such cases are gated on Γ, α, \|Vrel\|, ErrStat and iterations |
| W1 | `use, intrinsic :: ieee_arithmetic` in a per-point `pure` function makes gfortran call `_gfortran_ieee_procedure_entry/exit` on every call. Module-level `USE IEEE_ARITHMETIC` in `NWTC_IO` (`NWTC_IO.f90:25`) does the same to procedures that pick it up (seen in the pin's `ui_leaves_nograd`, 0.7 s) | `IfW_FlowField.f90:386` | Remove the `use`; values unchanged | `l0_ww_ieee`: `IfW_FlowField_GetVelAcc` on 1e5 points, pin vs overlay, bitwise; timing ratio logged |
| W2 | The `:1996` loop writes `m%Inflow(1)` with time `utimes(i)` (not `BEMT_utimes(i)`) and runs i = 2 then 1: the wake slice ends equal to `:1977`'s only because the last pass is i = 1 | `AeroDyn.f90:1993-2003` | Skip only `InflowWakeVel`; keep the rotor and GS recompute | `l0_ww_wake`: `m%Inflow(1:2)` after `AD_UpdateStates`, pin vs overlay, bitwise, with uniform and Grid3D wind |
| W3 | `r_wind` is packed `[LL CPs \| NW (all nNWMax+1 rows) \| FW \| grid \| source panels]` at **maximum** size, existing or not | `FVW_Subs.f90:568-638` (`:578-580` note) | Slice boundaries from the same counts; H3 and H4 touch only their slices | `l0_ww_slices`: nGridOut 0/1/2, SrcPnl on and off, FW 0 |
| W4 | H4's predicate is evaluated at the end of `UpdateStates` for the *next* `CalcOutput` (`t+DTaero`), with `bDoGrid`'s `DTaero/2` tolerances | `FVW.f90:1621-1625` | Use a window one `DTaero` wider on both sides (superset) | `l0_ww_grid`: grid VTK files bitwise vs pin for tStart/tEnd on, just inside and just outside step boundaries |
| W5 | `TwrInflArray` has a real `!$OMP` loop: serial in every CUDA build (`OPENMP=OFF`), 4.7 s vs 1.4 s (8 threads) on bar-r8 | `AeroDyn.f90:6123-6139` | H3/H4 cut its points; no port | (timing only) |
| W6 | Linearization calls `SetInputsForFVW` with `tIndx = 1` on perturbed inputs | `AeroDyn.f90:7635, 7654` | H3 keys on `tIndx = 2` inside `AD_UpdateStates` only | `l0_ww_lin`: `AD_JacobianPInput` on an OLAF case: Jacobian columns bitwise vs pin (OLAF linearization is unsupported, `Wake_Mod` 3 cannot linearize: then N-A, recorded) |


## FMA experiment (the fast mode's FMA policy, decided by measurement)

**Question.** The unit-4 analysis argued that this unit's fast mode must be FMA-free, because ulp changes could flip a per-solve convergence test (`maxval(|ΔΓ|)/MeanGamma < ConvCrit`, `FVW_Wings.f90:467`) and move Γ by up to `ConvCrit`·|Γ|, an exception to architecture §3. Measured instead. Same source and same binary, chosen at run time:
- **F0**: fast mode with no contraction (`OFK_CS_FMA=0`; TU `circulation_fast.cpp`: `-march=x86-64-v3 -ffp-contract=off`);
- **F1**: fast mode with contraction (`OFK_CS_FMA=1`; `circulation_fastfma.cpp`: `-march=x86-64-v3 -ffp-contract=fast`);
- **natural envelope**: `openfast-pin-ffma` (the pin built with `-march=x86-64-v3`, i.e. FMA, OPENMP=OFF, aerodyn_driver only) against base₁ (`openfast-pin-fomp` at 1 thread; `ov-off`/`ov-pin` are the same code without OpenMP). In L0: the pin's four files built twice in one process, without (`_ref`) and with FMA (`_rfm`).

**How the iteration counts were read.**
- The pin binaries are unmodified. `itertrace` (ptrace) puts breakpoints on the entry of `Wings_ComputeCirculationPolarData` and on its in-loop call of `CirculationFromPolarData`. The traced outputs are bitwise the untraced base₁ and envelope runs.
- The port's counts come from `OFK_CS_ITERS`. On the same port binary, traced and logged counts agree on all 300 solves.
- In L0, `AFI_ComputeAirfoilCoefs` calls are counted (`ld --wrap`), for every side alike.
- "Isolated" runs put the other units in strict mode (`OFK_MODE=strict OFK_CS_MODE=fast`), so they show this unit's FMA alone.

| | pin (base₁) | pin with FMA | F0 | F1 |
|---|---|---|---|---|
| bar-r8-tree: solves / iterations | 300 / 7720 (10 at the cap of 30) | 0 solves differ | 0 differ (isolated, production, 1 and 8 threads) | **0 differ** (isolated, production, 1 and 8 threads) |
| bar-r8-direct: solves / iterations | 300 / 7718 | 0 differ | 0 differ | **0 differ** |
| r-test OLAF cases (job 1925871; 9 cases with a solve, 492 solves / 15 641 iterations; ad_HelicalWakeInf / HelicalWake have none; ad_Sphere has nCP_tot = 0) | per case: BAR 60/1566, BAR_RFP2 60/1566, B1n2 60/977, EllipticalWingInf 10/84, Kite 50/4036, QuadRotor 10/298, VerticalAxis 37/1110, glue EllipticalWing 5/76, glue 5MW_Land_QuadRotor 200/5928 | 0 of 287 differ (aerodyn_driver cases; no FMA pin build of `openfast`) | 0 of 492 differ (isolated) | **0 of 492 differ** (isolated and production) |
| L0 random cases (10 seeds, 4 070 cases, F0 and F1 vs the pin) | — | 0 with another count | 0 differing values | 0 with another count |
| L0 knife edge (ConvCrit = the adjacent doubles where the pin's count changes; 16 solves per seed) | — | flips 2–5 of 16 | 0 | flips 2–6 of 16 (ΔΓ 1.4e-5–3.2e-5 relative when flipped) |
| R8 bar-r8-tree (envelope max 6.39e-13), max drift / worst ratio / channels > 2× | — | (envelope) | 0 isolated; 2.08e-13 / 0.37 / 0 of 542 production | 5.63e-14 / 0.36 isolated; **4.04e-13 / 0.69 / 0 of 542** production |
| R8 bar-r8-direct (envelope max 2.60e-13) | — | (envelope) | 0 isolated; 9.70e-14 / 0.82 / 0 production | 3.70e-14 isolated; **9.70e-14 / 0.82 / 0** production |
| R6 (ctest `aerodyn_utest\|OLAF\|5MW_MRSemi`, ov4b-omp, 1 thread) | 18/18 (strict) | — | 18/18 | **18/18** |
| Kernel time per solve, 8 threads, bar-r8-tree, 5 interleaved pairs (job 1925772) | 79.5 ms (wsamp `ui_quad_n1`) | — | 0.869–0.874 ms | **0.834–0.842 ms** (−3.8%) |
| same, bar-r8-direct, 3 pairs | — | — | 0.876–0.882 ms | **0.841–0.849 ms** (−3.9%) |
| same, 1 thread, bar-r8-tree, 5 pairs | 78.9 ms | — | 4.108–4.124 ms | **3.830–3.889 ms** (−6.3%) |
| Whole run, 8 threads, bar-r8-tree (5 pairs) | 44.0–44.2 s | — | 16.93–17.36 s | 16.91–17.36 s |
| Whole run, 1 thread, bar-r8-tree (5 pairs) | 97.97–98.09 s | — | 48.69–48.83 s | 48.61–48.76 s |
| Whole run, 8 threads, bar-r8-direct (3 pairs) | 128.2–128.7 s | — | 51.28–51.35 s | 51.07–51.22 s |
| perf_cs replay per call, nCP = 87, serial / 8 threads | 3.08 ms (analysis µbench) | — | 146 / 30.1 µs | 135 / 28.7 µs |

**Run-to-run noise** (max − min over the repetitions of one variant):
- kernel time per solve: 0.4–1.5% (8 threads: F0 0.6%, F1 1.0%; 1 thread: F0 0.4%, F1 1.5%);
- whole run: 0.43 s at 8 threads (2.5%) and 0.14 s at 1 thread (0.3%).

**Decision (rule fixed in advance): F1.**
1. F1 passes R6 (18/18) and R8: 0 of 542 channels over 2× the envelope in both cases, worst ratio 0.82.
2. F1's solve time (the kernel this unit runs, per solve) is lower than F0's by 3.8–6.3%, beyond the noise of 0.4–1.5%. The F0 and F1 ranges do not overlap in any of the three series.

So `kDefaultFma = 1`, consistent with architecture §3, and no exception is recorded. Two caveats:
- **Whole-run wall time does not separate F0 and F1** beyond its noise: about 0.1 s of 17–49 s. The gain is real but small: the kernel is now about 1.5% of the run.
- **Claim 3 at L0 (architecture §4: fast vs strict ≤ max(1e-12, 1.1 × the pin's own FMA spread on the same case)) is not in the rule.** Per case it holds in 4 068 of 4 070 L0 cases (10 seeds). It fails in 2 (seeds 6 and 8, ratios 3.9 and 24): the control points lie on their own bound segment (`onLL`), the near-line geometry where the cross product cancels, one of them unregularized (RegFunction 0). There F1's contraction moves Γ by 1.2e-11 and 2.4e-11 relative, against the pin with FMA's 2.8e-12 and 5.2e-13. Neither changed an iteration count. If the parent wants claim 3 as a condition for F1, these two cases flip the outcome to F0: one constant (`kDefaultFma`) or `OFK_CS_FMA=0`.
- What the analysis feared does exist, but only at a knife edge that real runs never hit, and the Fortran built with FMA has it at the same rate. It belongs to the natural envelope, not to the port (architecture §3 updated).

## Gate results (2026-11-17)
CPU = Xeon 8480+ (`med`, or the login node for L0); GPU = A100 MIG 2g.20gb (`interactive_gpu`, EPYC 7443 host).

| ID | Result | Evidence |
|---|---|---|
| Claim 1, L0 (strict and F0 vs the pin, OPENMP=OFF flags) | **PASS** | `u4a-omp` `l0_cs`, seeds 1–10, 249 cases each at 200 random cases (447 at 400). In every case: 0 differing values, 0 different iteration counts, for strict (threshold / serial / parallel_for) and F0 (serial / parallel_for), knife edges included. Logs `opt/build/u4a-jobs/l0_cs.login.seed*.log`; ctest 5/5 (`l0_cs`, `_seed7`, `_1thread`, `_llvel4`, `_threads`). 2 cases per seed are gated on Γ, α, \|Vrel\|, ErrStat and iterations only, because the pin itself is not reproducible there (ref vs ref, uninitialised reads with a zero-span wing) |
| Claim 1, system (bitwise incl. iteration counts) | **PASS** | job 1925585: on ov4-fomp, strict and F0 (isolated) = base₁ from line 8, on bar-r8-tree and bar-r8-direct (300 solves each, counts identical). R6 trees: `.outb` of F0 isolated = the traced pin (ov-off) on the r-test cases (job 1925871) |
| Claim 2 (CUDA build) | **PASS** | jobs 1925691 and 1926054: `u4a-cuda-a100` `l0_cs` seeds 1 and 7 with the LL induction by the Fortran segment tree, strict and F0 bitwise, as on the CPU. With unit 2's GPU LL kernel (claim 2 of unit 2) Γ moves ≤ 4.5e-16 relative, with no count change. GPU bar-r8-tree: F0 = `OFK_CS=fortran` bitwise, iteration counts = the pin's, R8 0/542 (4.29e-13, worst ratio 0.75) |
| R6 | **PASS** | job 1925871 on ov4b-omp: ctest `aerodyn_utest\|OLAF\|5MW_MRSemi_DLL_WSt_WavesIrr`: strict 18/18, F0 18/18, F1 18/18. The first try (job 1925641, ov4-omp) failed `ad_Sphere_OLAF` (nCP_tot = 0 abort), fixed |
| R7 | **PASS** | L0: md5-identical dumps of every port variant at 1/4/8 threads and across processes (`l0_cs_threads`). System: F0 and F1, 1 vs 8 threads, identical; F0 two processes identical (job 1925585); GPU two processes identical (job 1925691) |
| R8 | **PASS** | job 1925585 (table above): 0/542 channels over 2× in both cases for F0, F1 and strict; no channel ever > 1e-6. GPU: 0/542 |
| P2 (8 threads, vs openfast-pin-fomp, fomp config) | **PASS** | job 1925586 (best of 3): bar-r8-tree pin 44.01 s, F0 17.06, F1 16.96 (**2.59×**); bar-r8-direct pin 128.22 s, F0 50.91, F1 50.82 (**2.52×**); 1 thread: bar-r8-tree pin 97.97, F0 48.76, F1 48.64 (2.01×); bar-r8-direct pin 767.82, F0 321.71, F1 321.45 (2.39×). Kernel per solve at 1 thread, direct: F0 4.114–4.137 ms, F1 3.856–3.902 ms |
| P4 (r-test OLAF cases ≤ 1.10× the pin path) | **PASS** | job 1925871, ov4b-omp fast vs ov-off (both Fortran OPENMP=OFF), best of 3, 1 / 8 threads. Ratios: ad_BAR_OLAF 0.40 / 0.29; _RegFunctionPart2 0.40 / 0.29; ad_B1n2 0.97 / **1.06** (0.33 vs 0.35 s, Kokkos start-up); ad_EllipticalWingInf 0.90 / 0.90; ad_HelicalWakeInf 0.61 / 0.44 (no solve: the other units); ad_Kite 0.40 / 0.11; ad_QuadRotor 0.25 / 0.18; ad_VerticalAxis 0.34 / 0.18; ad_Sphere 1.03 / 1.02 (nCP_tot = 0: no kernel, start-up only); glue HelicalWake 0.61 / 0.43; glue EllipticalWing 0.91 / 0.91; glue 5MW_Land_QuadRotor 0.35 / 0.18 (230 → 42 s). Worst 1.06 |
| P6 | **PASS** | L0: 7 allocations (1.49 MB) in the first solve, none in the next 19. Full run (job 1925772, memtool): the 7 `ofk_cs_*` HostSpace allocations happen before the 4th of 903 evaluations, none afterwards |
| A1 | **PASS** | `nm` `libaerodynlib.a`: no `_MOD_` symbol of ov-pin (787 global) is missing in ov4-off, ov4-omp, ov4-fomp. The only one added is unit 2's `fvw_subs::fvw_kokkos_release`; this unit adds none (`ofk_cs_*` are C symbols, `ofk_cs_setup` is private) |
| A2 | **PASS** | `tests/circulation-solve/a2.sh`: all four overlaid files preprocess without `OPENFAST_KOKKOS` to the pin's, line for line. Files touched: the table above, all in `src/modules/{aerodyn,kokkos}` and `tests/circulation-solve`, plus `notes/` and the skill |
| A5 | **PASS** | `tests/overlay/a5.sh ov-pin ov4-off` (`opt/build/u4a-a5.log`): 76/76 `.mod` identical; 17/19 aerodynlib objects identical, **`FVW_Wings.f90.o` identical**; FVW_BiotSavart/FVW_Subs differ only in `movl $imm` line numbers (3 and 71 pairs), as on units 1–3. `DOUBLE_PRECISION=OFF`: `ov4-sp-omp` builds, and ad_BAR_OLAF runs with a `r4` handle, finite outputs and the same 1566 iterations |
| A6 | **PASS** (code) | no non-const `static` per-instance data; handles in the guarded map; multi-instance A7 sequence |
| A7 | **PASS** | `a7_sequence` (strict and F0): re-evaluation, a misc copy with its own handle used in between, Eps only, r_NW only and CP changed, orphaned copy; 14 solves bitwise with fresh pin solves, CPU and GPU builds |
| Mutations | **8/8 caught** | `tests/circulation-solve/mutate.py` on a copy of `src/` (`opt/build/u4a-mut-src`, log `opt/build/u4a-mutations.log`): wrong regularization branch, Γ product order (F10-R1 style), cache not invalidated (caught only by A7), off-by-one ring range, segment 3-4 reversed, Eps row 2, guards dropped, Γ index off by one depth |
| E1 | above | E2/E5/E6: the parent's ledger (this session wrote no ledger row, as instructed) |

### P4 thresholds (job 1925426, perf_cs, med, best of 5, µs per call: build / replay)
| nCP (pairs) | serial | parallel_for, 2 threads | parallel_for, 8 threads |
|---|---|---|---|
| 6 (288) | 2.23 / 0.99 | 2.88 / 2.11 | 4.06 / 3.24 |
| 9 (648) | 6.12 / 2.77 | 4.98 / 2.86 | 4.53 / 4.16 |
| 15 (1 800) | 10.1 / 4.67 | 7.37 / 3.81 | 5.67 / 5.08 |
| 18 (2 592) | 17.9 / 8.24 | 11.7 / 6.88 | 6.95 / 5.83 |
| 87 (60 552, BAR) | 310 / 148 | 167 / 81.4 | 53.0 / 30.3 |
| 150 | 924 / 466 | 467 / 251 | 131 / 75.4 |

So: build parallel from 512 pairs, replay from 2048 pairs (the 8-thread crossovers; 2 threads cross earlier, so these are conservative). With iNWStart = 2, the replay goes parallel from 16 control points (16 × 128 = 2048 pairs) and the build from 8. The r-test P4 row below is measured with these defaults.

## Slurm jobs (logs in `opt/build/u4a-jobs/`, outputs in `opt/build/u4a-<what>-<job>/`)
| Job | What | Result |
|---|---|---|
| 1925426 `-p4` | perf_cs, 1/2/4/8 threads, nCP 3…300 | thresholds above |
| 1925585 `-iters` | FMA experiment iteration counts, claim 1 system, R7, R8 (ov4-fomp) | above |
| 1925586 `-timing` | P2, 1 and 8 threads, best of 3; wsamp solve buckets | above |
| 1925641 `-rtest` (cancelled) | R6 on ov4-omp: found the nCP_tot = 0 abort | fixed |
| 1925871 `-rtest` | R6 strict/F0/F1, r-test iteration counts, P4 (ov4b-omp) | above |
| 1925772 `-fmatime` | interleaved F0/F1 timing; P6 memtool; perf_cs | above |
| 1925691, 1926054 `-mig` | CUDA build on A100 MIG: L0, bar-r8-tree, Sphere (ov4-cuda-a100, ov4b-cuda-a100) | above |
| 1926189 `-final` | ov4b-fomp (final source, default F1) = the experiment's F1/F0 outputs | bar-r8-tree and bar-r8-direct: default = experiment F1, `OFK_CS_FMA=0` = experiment F0, strict = base₁, all IDENTICAL from line 8; the log line shows `fma=1` |

Trees: `opt/build/u4a-omp`, `u4a-cuda-a100` (L0), `u4a-mut` (mutations), `ov4-{fomp,omp,off,sp-omp,cuda-a100}` (first source), `ov4b-{omp,fomp,cuda-a100}` (with the nCP_tot = 0 fix; `ov4b-fomp`/`-cuda-a100` rebuilt with the F1 default; `ov4b-omp` still defaults to F0 because job 1925871 used it, so rebuild it before reuse).

## Lessons
1. **Measure the Fortran's own FMA build on the discrete decision before declaring an FMA exception.** The convergence-flip hazard is real, but only within ulps of a knife edge, and the pin built with FMA flips there just as often. Real inputs (600 BAR solves, the r-test set, 4 070 random L0 cases) never hit it.
2. **Iteration counts without editing anything:** `ld --wrap` around a cross-file call in L0; a 100-line ptrace breakpoint counter for whole binaries (no gdb or perf on this cluster). Validate the tracer against an in-code counter on the same binary, and check that traced outputs are bitwise the untraced ones.
3. **ref vs ref per case** found pin non-determinism (uninitialised reads with a zero-span wing). A single global ref-ref run would have hidden which cases are valid oracles.
4. **R6 found what L0 did not:** a geometry with no lifting line at all (`ad_Sphere_OLAF`). Shims must accept empty extents. L0 now covers it.
5. **The replay was compute- and memory-bound, not latency-bound.** Storing 3 values per pair and recomputing the cross product (same bits) was faster than storing 6. Vector width came from the block size: 4 lanes gave 16-byte vectors, 8 lanes 32-byte ones.
6. **Amdahl moved again:** the circulation solve went from 58% of `bar-r8-tree` to about 3%. The run is now wind (9.7 s) plus the wake kernels (3.6 s) plus tower (1.5 s): unit 4b's territory.
7. The rule "-j 12 total" was briefly exceeded when a chained build script did not stop after a failed step and overlapped the next chain. Chain scripts need `|| exit 1`.

## Resume
Nothing is pending. Job 1925586 was still running only its optional wsamp profiles of bar-r8-direct when this card was finished; its P2 numbers above are complete. Before the trees are reused:
- `ov4b-omp` was rebuilt with the F1 default after job 1925871;
- the `ov4-*` trees (first source, without the nCP_tot = 0 fix) should not be used for new evidence;
- the parent session writes the ledger row (E2/E5/E6) and commits.
