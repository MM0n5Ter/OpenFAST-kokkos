# circulation-solve (unit 4, DRAFT: contract and audit list only)

Unit 4 of the OpenFAST-kokkos port, chosen by the host-time profile in `notes/olaf-unit4.md` (session `199ad883`, analysis subagent, 2026-11-17). Two parts with separate gates:
- **4a** ports the lifting-line circulation iteration of OLAF with the geometry hoisted out of the `do while` (F10-B, class **B**). It runs in the Kokkos *host* execution space in every build, CUDA included (L1 fails: the work scales with nSpan², not with the wake).
- **4b** removes three kinds of redundant host work at the OLAF wind points: the IEEE entry/exit overhead, wind and tower re-evaluations whose results nobody reads, and grid points outside the grid window. All class **B**. It changes no OLAF kernel; it touches AeroDyn, InflowWind and one FVW pack.

The device-resident wake (phase-2 residency: F2–F5, F9) is **not** in this unit. The profile puts it at ≤ 0.1% of the run (`notes/olaf-unit4.md` §2.3).

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
- **Handle:** the cache (about 7 reals per target × ring × depth × segment; 2.9 MB for nCP = 87, 25 MB for nCP = 300) lives in a per-instance scratch keyed like unit 2 (`m%Sgmt`). It is sized on the first call (nCP_tot is fixed for a run), because a per-call allocation would fail P6. It is **rebuilt from the Fortran arguments on every call and never trusted across calls**. So there is no residency, and A7 holds by construction.
- **Modes:** strict = fast = class **B**: no FMA in either TU (audit C6). The kernel shape is "square" (targets nCP_tot, sources nCP_tot·iNWStart·4). The per-target serial order is the pin's (ring-major, then depth, then segments 1-2, 3-4, 2-3, 4-1). The Kokkos `RangePolicy` runs over targets in `DefaultHostExecutionSpace`, with a serial fallback below a measured `N_min` (P4).

### 4b: host wind work at the OLAF wind points (AeroDyn, InflowWind, FVW pack)

| Item | Routine (pin) | Change | Class | Interface |
|---|---|---|---|---|
| H1 | `UniformField_GetVel` (`IfW_FlowField.f90:385-410`, `pure`, module-private) | Drop `use, intrinsic :: ieee_arithmetic` (`:386`): nothing from the module is used, and gfortran then saves and restores the FP environment on every call | B (bitwise: `scripts/unit4-mb` `ieeemb`, 0 differing of 5.6M points; 14.4× faster per point) | unchanged |
| H2 | `AD_UpdateStates` (`AeroDyn.f90:1942`), the two `AD_CalcWind` calls on `uInterp` (`:1996`) | Skip the `InflowWakeVel` part (`AD_CalcWind` `:2094-2121`) when `Wake_Mod = FVW`: after the loop, `m%Inflow(1)%InflowWakeVel` is the value from `:1977` (same `utimes(1)`, same `o%WakeLocationPoints`). Keep the rotor and GS parts | B | `AD_CalcWind` is private: add a private variant or an optional argument; `AD_CalcWind_Rotor`/`_GS` (public, `:65-66`) are unchanged |
| H3 | `SetInputsForFVW` from `AD_UpdateStates` (`AeroDyn.f90:2027` → `:4226-4235` `TwrInflArray` over **all** `r_wind` points, existing or not) | In `AD_UpdateStates`, only these slices of `V_wind` are ever read:<br>- `u(1)`: LL, NW, FW (`FVW.f90:645`, `:1378` through `uInterp`), plus source panels (`FVW_Subs.f90:2105`);<br>- `u(2)`: LL and panels (through `uInterp`, `FVW.f90:1378`).<br>The grid slice is read only from `AD_CalcOutput`'s own `SetInputsForFVW` (`FVW.f90:1618`). So apply `TwrInflArray` only to the read slices for these two calls. The `AD_CalcOutput` call (`:2378`) is unchanged | B | private (`SetInputsForFVW` is not public) |
| H4 | `SetRequestedWindPoints` (`FVW_Subs.f90:568`, grid block `:606-620`) | Write grid points only when `t+DTaero` is inside a grid's `[tStart − DTaero, tEnd + DTaero]` window, a superset of `bDoGrid` (`FVW.f90:1621-1625`); otherwise `(0,0,0)`, which `UniformField` skips (`z > 0` test, `IfW_FlowField.f90:182`; Grid3D `:229`). This removes 9600 of about 20.6k points from every wind call on the BAR cases. The tower loop still visits them: that is H3's job | B for every output and VTK file. `m%GridOutputs%uGrid` changes at non-output steps (misc, never written there) | unchanged |
| (later) | `AD_CalcOutput` wind + tower on the wake slice (`AeroDyn.f90:2365`, `:2378`) | Consumed only by grid outputs: needs an FVW query (deferred: Lane D / InflowWind on device) | B | — |

Callers that must keep working: `AD_UpdateStates`/`AD_CalcOutput` from `ADI_*` (driver), the OpenFAST glue (`FAST_Funcs.f90:459-466`), the linearization perturbations (`AeroDyn.f90:7629-7656`: they call `AD_CalcWind_Rotor`/`_GS` and `SetInputsForFVW` directly, so H3 must key on `tIndx`, not on the caller), and FAST.Farm/AWAE (`IfW_FlowField_GetVelAcc` has callers in 4 module directories: H1 speeds them all up, with bitwise-identical results).

## Audit list → L0 case

| # | Trap | Where (pin) | What the port does | L0 case |
|---|---|---|---|---|
| C1 | `ui_quad_n1`'s directives are `!OMP`, not `!$OMP`: the pin's loop is **serial** | `FVW_BiotSavart.f90:527-528, 548-549` | P2 compares against the serial pin (it is the baseline); do not "fix" the directive in the oracle | `l0_cs_*` oracle built from the pin unchanged |
| C2 | Γ enters `ui_seg_11` first in a left-associated product: `Kv = SegGamma*fourpi_inv*Kv*(norm_a+norm_b)/denominator` | `FVW_BiotSavart.f90:125` | Replay `((G·fourpi_inv)·Kv)·(na+nb)/den`, then `·crossprod` per component. Precomputing `fourpi_inv·Kv` or `A·Γ` (F10-R1) is **not** B | `l0_cs_bitwise` (strict = fast = pin, every RegFunction 0–4) |
| C3 | Early returns leave `Uindtmp = 0`, which the pin still adds (`Uind = Uind + 0`) | `FVW_BiotSavart.f90:73, 80, 85, 89` → `:533-546` | Skipping the add is exact: `Vvar` starts at +0 (`FVW_Wings.f90:424`) and a round-to-nearest sum from +0 never becomes −0. Keep a mask for the three guards (`<=`, NaN passes, as in `ui_seg_11`) | `l0_cs_degenerate`: CP on a ring edge, zero-length span, coincident points, Γ = 0 and Γ = −0 rings, NaN in `r_NW` |
| C4 | **Discrete decisions amplify rounding:** `maxval(abs(DGamma))/MeanGamma < ConvCrit` and `iIter == CircSolvMaxIter` | `FVW_Wings.f90:466-471` | Class B in both modes, so the iteration count can never change. An R1 variant (F10-R1, FMA) can add or drop an iteration and then moves Γ by up to ConvCrit · \|Γ\| | `l0_cs_iters`: iteration count equal to the pin's per solve, on captured BAR states (L0b, steps 1, 10, 100, 300) and on states tuned to converge at ConvCrit · (1 ± 1e-12) |
| C5 | Ring geometry and orientation: P1 = (s, d), P2 = (s+1, d), P3 = (s+1, d+1), P4 = (s, d+1); segments 1→2, 3→4, 2→3, 4→1 | `FVW_Wings.f90:430-433`, `FVW_BiotSavart.f90:533-546` | Literal; the order is part of the per-target sum | `l0_cs_bitwise`; mutation: swap 3→4 for 4→3 (must fail) |
| C6 | FMA changes the cache (e.g. `ya*zb-za*yb` cancellation, `biotsavart-direct` card) and so, through C4, the iteration count | `FVW_BiotSavart.f90:59-102` | Both cache TUs in `OFK_STRICT_SOURCES` (no contraction, nvcc irrelevant: host space) | `l0_cs_bitwise` in the fast build |
| C7 | RegParam is `Eps_NW(1, iSpan, iDepth)` for **all four** segments of a ring (row 1 only); unit 2's `LatticeToSegments` uses other rows for shed and FW segments | `FVW_Wings.f90:441` | Read row 1; do not reuse unit 2's segment SoA | `l0_cs_eps`: rows 2 and 3 set to garbage, results unchanged |
| C8 | Targets are flattened over wings, in wing order; `nSpan` may differ per wing | `FVW_Wings.f90:437-443` | Offsets from a per-wing prefix sum | `l0_cs_wings`: 1, 2 and 3 wings with unequal nSpan |
| C9 | `Vind_CP`/`Vtot_CP` are API-visible misc (VTK `FVW_IO.f90:559-560`) and are overwritten with `-9999` on exit; `alpha_LL`/`Vreln_LL` come from the last `CirculationFromPolarData` | `FVW_Wings.f90:453, 506-507, 582-583` | Stay in Fortran, unchanged | `l0_cs_bitwise` compares `m%W(:)%{Vtot_CP,Vind_CP,alpha_LL,Vreln_LL}` and `z` |
| C10 | `m%FirstCall` guess path; `GammaScale = 0` (t ≤ 0) skips the rescale; `EqualRealNos` | `FVW_Wings.f90:357-377`, `:255-268` | Fortran, unchanged; the cache does not depend on Γ | `l0_cs_firstcall`, `l0_cs_gscale0` |
| C11 | P6: the cache is O(nCP²) and the pin allocates `Vvar`/`Vcst`/`DGamma` per call (allowed: Fortran) | `FVW_Wings.f90:342-385` | Device-side (host space) cache in the handle, sized on the first call | memtool: 0 allocations after call 3 |
| C12 | Re-evaluation and copies (A7): `UpdateStates` twice at the same `t`; glue copies `x`, `z` CURR→PRED every step (`FAST_Funcs.f90:459`, `FAST_Solver.f90:1275`, `FAST_Subs.f90:5299`) | `FVW.f90:597-615` | Nothing is kept between calls | `l0_cs_a7`: same call twice, interleaved with a call on a copy; bitwise |
| W1 | `use, intrinsic :: ieee_arithmetic` in a per-point `pure` function makes gfortran call `_gfortran_ieee_procedure_entry/exit` on every call. Module-level `USE IEEE_ARITHMETIC` in `NWTC_IO` (`NWTC_IO.f90:25`) does the same to procedures that pick it up (seen in the pin's `ui_leaves_nograd`, 0.7 s) | `IfW_FlowField.f90:386` | Remove the `use`; values unchanged | `l0_ww_ieee`: `IfW_FlowField_GetVelAcc` on 1e5 points, pin vs overlay, bitwise; timing ratio logged |
| W2 | The `:1996` loop writes `m%Inflow(1)` with time `utimes(i)` (not `BEMT_utimes(i)`) and runs i = 2 then 1: the wake slice ends equal to `:1977`'s only because the last pass is i = 1 | `AeroDyn.f90:1993-2003` | Skip only `InflowWakeVel`; keep the rotor and GS recompute | `l0_ww_wake`: `m%Inflow(1:2)` after `AD_UpdateStates`, pin vs overlay, bitwise, with uniform and Grid3D wind |
| W3 | `r_wind` is packed `[LL CPs \| NW (all nNWMax+1 rows) \| FW \| grid \| source panels]` at **maximum** size, existing or not | `FVW_Subs.f90:568-638` (`:578-580` note) | Slice boundaries from the same counts; H3 and H4 touch only their slices | `l0_ww_slices`: nGridOut 0/1/2, SrcPnl on and off, FW 0 |
| W4 | H4's predicate is evaluated at the end of `UpdateStates` for the *next* `CalcOutput` (`t+DTaero`), with `bDoGrid`'s `DTaero/2` tolerances | `FVW.f90:1621-1625` | Use a window one `DTaero` wider on both sides (superset) | `l0_ww_grid`: grid VTK files bitwise vs pin for tStart/tEnd on, just inside and just outside step boundaries |
| W5 | `TwrInflArray` has a real `!$OMP` loop: serial in every CUDA build (`OPENMP=OFF`), 4.7 s vs 1.4 s (8 threads) on bar-r8 | `AeroDyn.f90:6123-6139` | H3/H4 cut its points; no port | (timing only) |
| W6 | Linearization calls `SetInputsForFVW` with `tIndx = 1` on perturbed inputs | `AeroDyn.f90:7635, 7654` | H3 keys on `tIndx = 2` inside `AD_UpdateStates` only | `l0_ww_lin`: `AD_JacobianPInput` on an OLAF case: Jacobian columns bitwise vs pin (OLAF linearization is unsupported, `Wake_Mod` 3 cannot linearize: then N-A, recorded) |
