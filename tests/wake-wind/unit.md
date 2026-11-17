# wake-wind (OLAF unit 4b)

Unit 4b of the OpenFAST-kokkos port (`notes/olaf-unit4.md` §2.2, §3.1): remove dead host work at the OLAF wind points in AeroDyn, items H2-H4. Pure Fortran in the overlay, no Kokkos, class **B** (bitwise in every mode: the code is mode-independent host Fortran). Runtime A/B switch: `OFK_WW=fortran` runs the pin's code path in the same binary.

## Status (2026-11-17)

**Phase 1 done** (port, L0, A1/A2/A5): see "Gate results". **Phase 2** (R6, R7, outputs/VTK bitwise on the r-tests and the glue, P2, the wind re-measurement) is open: see "Resume".

## Files

| File | Change |
|---|---|
| `src/modules/aerodyn/src/AeroDyn.f90` | **new overlay** (pin copy + 161 lines under `#ifdef OPENFAST_KOKKOS`): H2 (`AD_CalcWind` optional `ofkSkipWake`), H3 (`SetInputsForFVW` optional `ofkStep`, `ofk_ww_TwrInflSlices`), H4 (`ofk_ww_zero_dead_grids`), `ofk_ww_enabled` (reads `OFK_WW` on every `AD_UpdateStates`, no module state), `ofk_ww_dump` (L0 hook, `OFK_WW_DUMP=<file>`) |
| `src/modules/aerodyn/src/FVW_Subs.f90` | adds `ofk_ww_slices` (slice bounds of the packed `r_wind`, same counts as `SetRequestedWindPoints`) and `ofk_ww_gridLive` (grid window predicate); unit 2/3/4a code untouched |
| `tests/wake-wind/l0_ww.sh`, `mkbts.py`, `mutate.py` | L0 case matrix, synthetic Grid3D wind file, mutation checks |

`FVW.f90` is not changed by this unit.

## Contract (frozen, pin `12fa3fc`)

> **InflowWind is an external library call** (user decision, 2026-11-17; decisions #14). 4b only changes how often and where AeroDyn calls it. **H1 and W1 are out of scope:** they are recorded only as an upstream suggestion for InflowWind's owners.

| Item | Routine (pin) | Change | Class | Interface |
|---|---|---|---|---|
| ~~H1~~ (upstream only) | `UniformField_GetVel` (`IfW_FlowField.f90:385-410`, `pure`, module-private) | Drop `use, intrinsic :: ieee_arithmetic` (`:386`): nothing from the module is used, and gfortran then saves and restores the FP environment on every call | B (bitwise: `scripts/unit4-mb` `ieeemb`, 0 differing of 5.6M points; 14.4× faster per point) | unchanged |
| H2 | `AD_UpdateStates` (`AeroDyn.f90:1942`), the two `AD_CalcWind` calls on `uInterp` (`:1996`) | Skip the `InflowWakeVel` part (`AD_CalcWind` `:2094-2121`) when `Wake_Mod = FVW`: after the loop, `m%Inflow(1)%InflowWakeVel` is the value from `:1977` (same `utimes(1)`, same `o%WakeLocationPoints`). Keep the rotor and GS parts | B | `AD_CalcWind` is private: optional argument `ofkSkipWake`; `AD_CalcWind_Rotor`/`_GS` (public, `:65-66`) are unchanged |
| H3 | `SetInputsForFVW` from `AD_UpdateStates` (`AeroDyn.f90:2027` → `:4226-4235` `TwrInflArray` over **all** `r_wind` points, existing or not) | Apply `TwrInflArray` only to the slices read after the call:<br>- `tIndx = 1`: LL, NW, FW (`FVW.f90:645`, `:1378` through `uInterp`), source panels (`FVW_Subs.f90:2105`), **and the grids `FVW_FinalWrite` may write** (see W7);<br>- `tIndx ≥ 2`: LL and panels (through `uInterp`, `FVW.f90:1378`).<br>Selected by an optional `ofkStep` argument passed only from `AD_UpdateStates`; the slices are a function of `tIndx` and `n`. The `AD_CalcOutput` call (`:2378`) and the linearization calls are unchanged | B | private (`SetInputsForFVW` is not public) |
| H4 | wind points of the grid outputs (`FVW_Subs.f90:606-620`) | **As built:** after `FVW_UpdateStates(t, n)`, in `AD_UpdateStates` where `OtherState%WakeLocationPoints = m%FVW%r_wind` (`AeroDyn.f90:2035`), the points of a grid are set to `(0,0,0)` unless `t+DTaero` or `(n+1)·DTaero` lies in `[tStart − 2·DTaero, tEnd + DTaero]` (superset of every reader's `bDoGrid`, W4/W7). `(0,0,0)` is skipped by the wind (`z ≤ 0`, `IfW_FlowField.f90:182`; Grid3D `:229`). `m%FVW%r_wind` is **not** changed, so FVW and the tower loops see the pin's points. *(Deviation from the planned site `SetRequestedWindPoints`: that routine does not see `t`, and passing it would need an `FVW.f90` edit; this site also leaves FVW's own state as in the pin.)* | B for every output and VTK file. `m%GridOutputs%uGrid` is only rebuilt when `bDoGrid` | unchanged |
| (later) | `AD_CalcOutput` wind + tower on the wake slice (`AeroDyn.f90:2365`, `:2378`) | Consumed only by grid outputs: needs an FVW query (deferred: Lane D / InflowWind on device) | B | — |

Callers that must keep working: `AD_UpdateStates`/`AD_CalcOutput` from `ADI_*` (driver), the OpenFAST glue (`FAST_Funcs.f90:459-466`, also with `NumCrctn ≥ 1`), the linearization perturbations (`AeroDyn.f90:7629-7656`: they call `AD_CalcWind_Rotor`/`_GS` and `SetInputsForFVW` directly without `ofkStep`, so they run the pin's code), and FAST.Farm/AWAE (unchanged: 4b does not touch `IfW_FlowField_GetVelAcc`).

## Audit list → L0 case (with the readers, `file:line` on the pin)

L0 = `tests/wake-wind/l0_ww.sh`: `aerodyn_driver` on variants of r-test `ad_BAR_OLAF` (TMax 2 s, DT = DTfvw = 0.1 s, `TwrPotent 1`, `TwrShadowOnWake True`) and `ad_QuadRotor_OLAF` (4 rotors, each with a tower). Each case runs three times in **one** binary: `OFK_WW=fortran` twice (reproducibility) and the default. Every output file must be byte-identical: `.outb`, every VTK file (wake, blades, grids) and the `OFK_WW_DUMP` file (at every `AD_UpdateStates`: `m%Inflow(1:2)%InflowWakeVel` outside the grid slices, `V_wind` of `u(1)` on LL+NW+FW+panels and of `u(2)` on LL+panels). The grid slices are checked through the grid VTK files, which are their only readers.

| # | Trap | Where (pin) | Readers enumerated (evidence) | What the port does | L0 case |
|---|---|---|---|---|---|
| W1 | (upstream only) IEEE wrapper in `UniformField_GetVel` | `IfW_FlowField.f90:386` | — | nothing (decision #14) | — |
| W2 | The `:1996` loop writes `m%Inflow(1)` with time `utimes(i)` and runs i = 2 then 1: the wake slice ends equal to `:1977`'s only because the last pass is i = 1 | `AeroDyn.f90:1993-2003` | `m%Inflow(1)%InflowWakeVel` is read only by `SetInputsForFVW` (`:4220`, called at `:2027` after the loop) and overwritten by `AD_CalcOutput`'s `AD_CalcWind` (`:2365`) before its own `SetInputsForFVW` (`:2378`). Inside the loop `SetInputs` reads `InflowInterp` (copied at `:1983`, before the loop), not `m%Inflow(1)`. The wake part depends only on `t`, `FlowField`, `o%WakeLocationPoints`, `p%FVW%MHK/WtrDpth`, `p%CompSeaSt`, `p%WaveField`, `m%WaveField_m` (`:2094-2121`), never on `u`: pass i=1 is the `:1977` i=1 call, pass i=2 is the `:1977` i=2 call (so no abort path is lost) and is overwritten unread | Skip only `InflowWakeVel` on both passes (`ofkSkipWake`); keep the rotor and GS recompute | `l0_ww_wake`: `base` (uniform), `grid3d` (synthetic periodic TurbSim, `mkbts.py`), `quad` (multi-rotor) |
| W3 | `r_wind` is packed `[LL CPs \| NW (all nNWMax+1 rows) \| FW \| grid \| source panels]` at **maximum** size, existing or not | `FVW_Subs.f90:568-638` (`:578-580` note); total `FVW.f90:280-289`; panel offset `FVW_Subs.f90:745-767` | `ofk_ww_slices` uses the counts of `SetRequestedWindPoints` and of `DistributeRequestedWind_Panels` (`FWnSpan` included, FW only when `nFWMax > 0`) | Slice boundaries from the same counts; H3 and H4 touch only their slices | `l0_ww_slices`: `g0` (nGridOut 0), `g2` (2 grids, one coarse in time), `fw0` (nFWPanels 0), `srcpnl` (sphere source panels moved to (−10, 0, 100) m, in the tower's potential flow), `dtfvw` (DTfvw = 2·DTaero) |
| W4 | H4's predicate is evaluated at the end of `UpdateStates` for the *next* `CalcOutput` (`t+DTaero`), with `bDoGrid`'s `DTaero/2` tolerances | `FVW.f90:1621-1625` | Wind at grid points comes only from `OtherState%WakeLocationPoints` through `AD_CalcWind` (`AeroDyn.f90:2094-2120`, the only reader; set at `:1072`, `:2035`). Its grid slice reaches `uGrid` only in `WriteVTKOutputs` (`FVW.f90:1618`), whose results are used only when `bDoGrid` (`:1626-1641`; `InducedVelocitiesAll_OnGrid_Calc` `FVW_Subs.f90:1171` and `WrVTK_FVW_Grid` `FVW_IO.f90:669-680`). Callers: `FVW_CalcOutput` (`:1539, 1541`) at `t+DTaero`, `FVW_FinalWrite` (W7) | Window `[tStart − 2·DTaero, tEnd + DTaero]` on `t+DTaero` **or** `(n+1)·DTaero` (superset with DTaero/2 to spare for every reader) | `l0_ww_grid`: `gedge` (6 grids: window on a step, on a half step, ±1e-7 inside/outside half steps, at t = 0, past TMax) |
| W5 | `TwrInflArray` has a real `!$OMP` loop: serial in every CUDA build (`OPENMP=OFF`) | `AeroDyn.f90:6123-6139` | point-wise: no state, `getLocalTowerProps`' status discarded, so one call per slice = the full call bitwise | H3/H4 cut its points; no port | (timing only: phase 2) |
| W6 | Linearization calls `SetInputsForFVW` with `tIndx = 1` on perturbed inputs | `AeroDyn.f90:7635, 7654` | they do not pass `ofkStep`: pin path. Also unreachable with OLAF: `AD_Init` rejects linearization unless `Wake_Mod` is 0 or 1 (`AeroDyn.f90:4894`); `dXdu` "not supported for FVW" (`:7680`) | none | `l0_ww_lin`: **N-A** (OLAF cannot linearize; and the linearization path is the pin's code by construction) |
| W7 (new) | `FVW_FinalWrite` reads the **grid slice of `u(1)`** left by the last `AD_UpdateStates` (the contract said the grid slice was read only from `AD_CalcOutput`) | `FVW.f90:497-519` (from `FVW_End` `:538` with `m%FVW_u(1)`, `AeroDyn.f90:1891`), `t = (iStep+1)·DTaero`, `iStep = n` (`FVW.f90:594`) | so H3 keeps the tower on the grids with `(n+1)·DTaero` in the window, and H4's window also covers a read at `(n+2)·DTaero` of points set at step `n` | H3: grid slices of `tIndx = 1` kept when live at `(n+1)·DTaero` | `l0_ww_grid`: `gend` (grids starting at TMax + 0.04, + 0.0499999, + 0.06; one at TMax − 0.05) |

Other readers checked and excluded: `have_nan` (`FVW_Subs.f90:443-488`, reads every `u%V_wind`) runs only `if (DEV_VERSION)` (`FVW.f90:760`), a `.False.` parameter (`FVW_Subs.f90:60`); `FVW_CalcContStateDeriv` reads `m%W%Vwnd_NW/FW` (`FVW.f90:952-964`; the `V_wind` line `:843` is commented); `CalcOutputForAD` (`:1444` commented); UA reads `W%Vwnd_LL`, not `V_wind`. `m%FVW_u` has 3 entries (`AeroDyn.f90:529`): with `InterpOrder = 2`, `tIndx = 3` is treated like 2 (read only through `uInterp`). Glue corrections repeat `AD_UpdateStates` at the same `(t, n)`: every predicate depends on `(tIndx, n, t)` and the counts only. MHK with `CompSeaSt`: `m%WaveField_m` is per-call interpolation scratch plus a one-time clamp-warning flag (`SeaSt_WaveField_Types.f90:57-62`): values unaffected; with MHK, `(0,0,0) + (0,0,WtrDpth)` is evaluated (no saving, still unread).

## Gate results (2026-11-17)

Evidence under `opt/build/u5-jobs/`. Builds: `ov5-omp` (OPENFAST_KOKKOS=ON, `aerodyn_driver`), `ov5-off` (OFF, `aerodynlib`), `ov5-sp-omp` (single precision, `aerodynlib`), baseline `ov-pin`.

| Gate | Result | Evidence (`opt/build/u5-jobs/`) |
|---|---|---|
| A1 | **PASS** | `nm` `libaerodynlib.a`, global `_MOD_` symbols: ov-pin 1063 = ov5-off 1063 (same set); ov5-omp 1066, none missing; added `fvw_subs::ofk_ww_slices`, `ofk_ww_gridlive` (and unit 2's `fvw_kokkos_release`). The AeroDyn routines are private. `a1.*.txt` |
| A2 | **PASS** | `tests/circulation-solve/a2.sh`: `AeroDyn.f90` (161 lines added) and `FVW_Subs.f90` preprocess to the pin's files without `OPENFAST_KOKKOS`. Files changed: the two overlay files, `tests/wake-wind/`, the two cards |
| A5 | **PASS** | `tests/overlay/a5.sh ov-pin ov5-off`: 76 `.mod` identical, `flags.make` and per-file flags identical; disassembly identical except source-line immediates of the runtime I/O calls (AeroDyn 3 lines, FVW_Subs 71, FVW_BiotSavart 3: the `#ifdef` lines shift line numbers; same class as unit 4a). `DOUBLE_PRECISION=OFF`: `ov5-sp-omp` `aerodynlib` builds. `a5.log`, `ov5-sp-omp.chain.log` |
| L0 W2 `l0_ww_wake` | **PASS** | `base`, `grid3d`, `quad`: 0 differing files (ref/ref and pin path/4b). `l0.verdicts`, `l0/l0_ww.log` |
| L0 W3 `l0_ww_slices` | **PASS** | `g0`, `g2`, `fw0`, `fwlive`, `srcpnl`, `dtfvw`: 0 differing files |
| L0 W4/W7 `l0_ww_grid` | **PASS** | `gedge` (8 grid VTK), `gend` (4 grid VTK, including the `FVW_FinalWrite` ones): 0 differing files |
| L0 W6 `l0_ww_lin` | **N-A** | OLAF cannot linearize (`AeroDyn.f90:4894`, `:7680`), and the linearization calls run the pin's code (no `ofkStep`) |
| Mutations | **6/6 caught** | window edge → `gedge gend base`; LL slice short → `base`; panel slice late → `srcpnl`; FW left out of the packing → `fwlive`; wake slices keyed on `tIndx = 2` → `base`; H2 skip on the wrong pass → `grid3d`. `mutations.log`, `mutations2.log`. First round missed two: with 35 NW panels the far wake starts only at 3.5 s > TMax (added `fwlive`), and steady wind is time-independent, so the H2 pass mutation needs `grid3d` |
| (timing glimpse, not P2) | — | login node, one run each, "Simulation CPU Time": `base` 1.07 → 0.68 s, `quad` 14.6 → 12.8 s |

## Lessons

- The contract's reader list missed `FVW_FinalWrite` (W7): it is an end-of-run reader of `u(1)` that a step-by-step reading of `UpdateStates`/`CalcOutput` does not show. Grep every caller of the reading routine (`WriteVTKOutputs`), not just the obvious one.
- A dump hook (`OFK_WW_DUMP`) of exactly the slices the contract says are read makes a dead-work removal checkable at every step, and catches slice mutations whose effect on `.outb` would be small.

## Resume

Phase 2 (a new agent) runs R6, R7, outputs/VTK bitwise, P2 and the wind re-measurement. Nothing is queued.

**Build trees** (`opt/build/`):
- `ov5-omp` (omp config, ON): configured; only `aerodyn_driver` and its libraries are built. Safe to reuse: run `tools/openfast-build.sh omp opt/build/ov5-omp` for every target (openfast, `aerodyn_utest`).
- `ov5-off` (OFF, `aerodynlib` only) and `ov5-sp-omp` (single precision, `aerodynlib` only): A5 evidence, safe to reuse.
- `u5-mut`, `u5-mut-repo`: mutation trees. **Unsafe**: the last build holds a mutant.
- Still to build: `ov5-fomp` (P2 config) and `ov5-cuda-a100` (GPU host is where W5's serial tower loop costs most). Chain them with `|| exit 1`, at most `-j 12` in total.

**Commands and expected results:**
1. L0 again on any new build: `tests/wake-wind/l0_ww.sh <tree>/modules/aerodyn/aerodyn_driver opt/build/u5-jobs/l0-<tree>`. Expect 11 PASS lines (about 1.5 min on the login node).
2. R6: `ctest` in `ov5-omp` with the default and with `OFK_WW=fortran`, both in strict and fast modes. Include `aerodyn_utest` and the glue OLAF cases `5MW_Land_QuadRotor_DLL_WSt_OLAF`, `HelicalWake_OLAF`, `EllipticalWing_OLAF`, `MHK_RM1_Floating_MR` and `5MW_MRSemi_DLL_WSt_WavesIrr`. One of them has `NumCrctn = 1`: check which and keep it, because it is the corrections case. MHK exercises `PosOffset`: H4 saves nothing there but must stay bitwise. Also run `py_ad_B1n2_OLAF` (Python interface). Expect the same pass/fail set as `OFK_WW=fortran`.
3. Outputs and VTK bitwise: for every r-test `ad_*OLAF` case and the glue cases above, run the default and `OFK_WW=fortran` in one binary and `cmp` every output file. The `run`/`same` functions in `l0_ww.sh` do this. Expect 0 differing files, in strict and in fast mode.
4. P2 and the wind re-measurement: `bar-r8` and `bar-perf` in `ov5-fomp` and `ov5-cuda-a100`, default vs `OFK_WW=fortran` (job `ofk-wake-wind-perf`, `med`). Expected per step: wind points go from 5 × (wake + grid) to 3 × wake outside the grid windows; tower points go to about 0.5× (r8) to 0.6× (perf) (`notes/olaf-unit4.md` §2.2). There is no H1, so each point still pays the IEEE wrapper: re-measure; do not use the §3.2 projection.
