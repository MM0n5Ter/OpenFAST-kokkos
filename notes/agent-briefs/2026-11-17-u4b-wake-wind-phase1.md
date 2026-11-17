<!-- Brief for unit 4b phase 1 (port + L0 + API gates), from TEMPLATE-port.md. Cost: tools/cost/ledger.csv. -->

You are doing **phase 1 (port + L0 + API gates)** of unit **4b `wake-wind`** of the OpenFAST-kokkos port. Repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. The repo's `CLAUDE.md`, `notes/process.md` and the `ofk-port-unit` skill apply. This is the first unit run under the cost rules in process.md §6, so its cost will be compared with earlier units (`notes/runs.md`).

**Goal.** Remove dead host work at the OLAF wind points in AeroDyn: items **H2, H3, H4** (pure Fortran, no Kokkos). This class-B change must give bitwise the same outputs (and VTK files) as the build without it, in strict and in fast mode. Phase 1 is done when the code is in the overlay, the L0 tests for W2–W4 and W6 pass, and A1/A2/A5 pass. Whole-run gates (R6, R7, outputs/VTK bitwise, P2, the wind re-measurement) are **phase 2**, which a separate agent runs from your "Resume" notes.

**Read first** (grep and read ranges; do not read whole Fortran files):
- `tests/circulation-solve/unit.md`, section "### 4b" (the contract table H2–H4, callers) and audit rows **W2–W6**. Create `tests/wake-wind/unit.md` and move the 4b contract and the W rows there (in `tests/circulation-solve/unit.md`, leave one line pointing to the new card).
- `notes/olaf-unit4.md` §2.2 (call graph, sizes, which slices are read), §3.1 (4b scope), §3.3 item 3 (the deadness proofs are the hardest gate).
- `notes/overlay.md` (how an overlay file is added), and the skill's lessons for units 3 and 4a.

**Already decided** (do not reopen):
- #14: **InflowWind is an external library call.** Do not edit, overlay or port anything in InflowWind. H1/W1 are upstream suggestions only.
- #12: no CPU↔GPU bitwise requirement. #15: unit 4a's fast mode defaults to F0. Neither changes here, because 4b is mode-independent host Fortran.
- The pin `deps/openfast` is never edited.

**Ruled out:** porting the wind or tower evaluation to C++/Kokkos (Lane D, later); changing `AD_CalcOutput`'s wind or tower call (the "(later)" row).

**Scope.**
- May change: overlay copies of `modules/aerodyn/src/AeroDyn.f90` (new overlay file) and `src/modules/aerodyn/src/FVW_Subs.f90` (already overlaid: add, don't disturb unit 2/4a code), `FVW.f90` only if unavoidable; `tests/wake-wind/` (new); notes. Put every change under `#ifdef OPENFAST_KOKKOS`, so that preprocessing without it gives the pin's file (A2). Add a runtime A/B switch, `OFK_WW=fortran`, that restores the pin's behaviour, like `OFK_CS=fortran`, so phase 2 can compare in one binary.
- Each change must be legal because **no reader exists**. Enumerate the readers in the card (W2–W4 evidence columns, `file:line` on the pin), including the glue with corrections (`NumCrctn` ≥ 1), grid outputs on and off, source panels, linearization (W6: key H3 on `tIndx`, not on the caller).
- Build trees: prefix `ov5-*` (e.g. `tools/openfast-build.sh omp opt/build/ov5-omp`) and `u5-*` for test trees. **Do not touch** `ov4b-*`, `ov-cuda-a100`, `ov-tests-iv-cuda-a100`, `u3-cuda-a100`, `u3d-*`, `fork*`, `openfast-pin*`: queued or reference jobs use them. At most `-j 12` in total; chain builds with `|| exit 1`.
- Jobs: none are needed in phase 1 (login node only). If you need one, name it `ofk-wake-wind-<what>` and use `med`.
- No commits and no ledger rows.

**Gates for phase 1:** A1, A2, A5 (`tests/overlay/a5.sh`, `tests/circulation-solve/a2.sh` as templates). L0 (login node): `l0_ww_wake` (W2), `l0_ww_slices` (W3), `l0_ww_grid` (W4), `l0_ww_lin` (W6, or N-A with the reason recorded). Each compares pin vs overlay **bitwise**: `m%Inflow(1:2)` and the `V_wind` slices that are read, uniform and Grid3D wind. Before writing a new harness, check whether a whole-driver comparison on a short r-test case covers an item more cheaply; L0 should cover the cases r-test does not (nGridOut 0/1/2, grid windows at step boundaries, SrcPnl on/off, FW = 0). Mutation check: break each predicate (window edge, slice boundary, `tIndx`) on a copy of `src/`; L0 must catch each one.

**Output discipline.** Scripts print one verdict line per check and write details to log files under `opt/build/u5-jobs/`. Read logs with `grep`. Reuse `tests/circulation-solve/*.sbatch`, `tests/induced-velocities/run_case.sh` and `r8.py` patterns rather than writing new scaffolding.

**Budget.** About **80 requests / 7 USD** (API-equivalent). At the cap, write the state into `tests/wake-wind/unit.md` "Resume" and hand back.

**Hand back:**
1. a gate table (ID, PASS/FAIL, evidence path);
2. the deadness argument per item (H2, H3, H4): the readers enumerated and why none reads the removed work;
3. the "Resume" section for phase 2: build trees ready, commands to run, expected results;
4. open issues as questions, changed files, and build trees that are safe or unsafe to reuse.
