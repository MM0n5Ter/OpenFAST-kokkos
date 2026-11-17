<!-- Brief given to subagent ae5c936b (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are porting **unit 4a, `circulation-solve`**, of the OpenFAST-kokkos port, and running an **experiment that decides the FMA policy for it**.

The repo is `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. Read these first, in this order:
1. `CLAUDE.md`
2. `notes/architecture.md`
3. `notes/criteria.md`
4. `notes/overlay.md`
5. `.claude/skills/ofk-port-unit/SKILL.md`: your procedure.
6. **`notes/olaf-unit4.md`**: the analysis and the profile. §2.1 and §3 describe 4a: what it is, its design and its gates.
7. **`tests/circulation-solve/unit.md`**: the draft card, with the Contract and Audit list. You complete it.
8. As worked examples of finished units: `tests/induced-velocities/unit.md` and `tests/treecode/unit.md`.

## The unit
The analysis recommends **F10-B**. `Wings_ComputeCirculationPolarData` (`FVW_Wings.f90:422-447` on the pin) recomputes the Γ-independent segment geometry in every iteration of the `do while`, through `ui_quad_n1`→`ui_seg_11`. That loop is serial because of the `!OMP` typo, and it averages about 26 iterations per solve. It is 58% of the `bar-r8-tree` run, 23.7 s of 40.7 s.

The fix:
- build the geometry once per solve, as a kernel in the Kokkos **host** execution space (CUDA builds have it too);
- replay only `ui_seg_11`'s Γ statement per iteration, in the pin's order;
- shim only the `Vvar` block, overlaying `src/modules/aerodyn/src/FVW_Wings.f90`;
- keep the geometry cache in a per-instance handle (P6).

Follow the card for the contract, the audit rows and the gates. The 4b part of the card (wind and tower) is **not** in your scope. InflowWind is an external library and must never be edited.

## The FMA experiment: the user wants this answered by data
The analysis argued that fast mode must be FMA-free for 4a. The reason: ulp changes can flip the per-solve convergence test, which shifts Γ by up to `ConvCrit`·|Γ| (relative 1e-3). That would be an exception to architecture §3, which keeps FMA in fast mode. **Do not assume either answer.** Build fast mode both ways and measure:
- **Variant F0:** fast without FMA. If the order is the pin's, this should be bitwise identical to strict and the pin.
- **Variant F1:** fast with FMA (and the fast-mode freedoms the architecture allows), same source. `OFK_FAST_FMA` in `src/cmake/ofk_options.cmake` may already give you the switch; otherwise add one.
- **Natural envelope:** the Fortran pin built with FMA vs without FMA, `opt/build/openfast-pin-ffma` vs `ov-off`/`ov-pin`. Check what these trees are before using them. This shows whether the Fortran *itself* flips iteration counts under FMA.

Measure on `bar-r8-tree` and `bar-r8-direct` (`opt/build/iv-cases`), plus the r-test OLAF cases:
1. **Per-solve iteration counts** of F0, F1, the pin and the pin with FMA. Count how many solves differ, by how many iterations, and when. Add an env-gated counter dump to the shim, or use whatever is cheapest and does not change results.
2. **R8:** drift of F0 and F1 against base₁ (`opt/build/iv-r8-1860700/<case>.fomp1/ad_driver.out`), using the FMA envelope (`opt/build/iv-r8ffma-1860993/<case>.ffma/ad_driver.out`) and `tests/induced-velocities/r8.py`. Report channels over 2× and the maximum drift.
3. **Speed:**
   - time per solve and whole-run wall time for F0 vs F1 vs the pin, at 1 and 8 threads, best of 3, on the `med` partition;
   - include a serial-vs-parallel threshold for P4. There are only ~87 targets per iteration.

**Decision rule, fixed in advance; apply it and report:**
- **Pick F1** (FMA, consistent with architecture §3) only if F1 passes R6 and R8 (≤ 2× envelope on every channel) **and** F1's solve time is faster than F0's by more than the run-to-run noise (state the noise).
- **Otherwise pick F0** (no FMA, bitwise for free), and record it as a declared, evidence-based exception in architecture §3.
- Either way, write the data table into the card and into a short "FMA experiment" subsection of `notes/olaf-unit4.md`.

## Gates
Everything in the card for 4a:
- claim 1, bitwise including iteration counts, for strict and for F0;
- R6;
- R7;
- R8;
- P2, against `openfast-pin-fomp` and the `fomp` config, 8 threads;
- P4;
- P6;
- A1, A2, A5 (`tests/overlay/a5.sh`) and A7;
- mutation checks: plant 4–6 meaningful bugs, e.g. wrong regularization branch, Γ order, cache not invalidated when `Eps` or the geometry changes, off-by-one in the segment range. Expect the L0 test to catch each one.

Also check correctness on CUDA builds (the kernel runs in the host space there, but the build and link must work):
- run one MIG job on `interactive_gpu` (A100 2g.20gb, no queue, `-c 4`) with an `ov4-cuda-a100` tree;
- the GPU queues for full GPUs are days long: do not submit full-GPU jobs.

## Hard rules
- Never edit `deps/openfast` and never use `fork/`.
- New build trees go in new directories, e.g. `opt/build/ov4-*` and `u4a-*`. Do **not** rebuild or reconfigure trees that queued jobs use: `ov-cuda-a100`, `ov-tests-iv-cuda-a100`, `u3-cuda-a100`, `u3d-frozen-cuda-*`, `openfast-pin-fomp`, `fork-cuda`, `iv-cuda`, `iv-omp`. Running them is fine.
- Compile on the login node with `module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0` and at most `-j 12`.
- Runs on AMD GPU nodes need the portable OpenBLAS on `LD_LIBRARY_PATH` (see the skill).
- `srun` fails from this shell; use `sbatch`. Name jobs `ofk-circulation-solve-<what>`.
- Never edit a script while a job runs it.
- Compare text outputs from line 8 on.
- The Python environment is `opt/venv-rtest`.
- No git commits, no ledger row. The parent session reviews and commits.
- If the design is wrong, or a gate cannot pass, record the evidence in the card and report it. Do not work around it.

## Deliverables
- **`tests/circulation-solve/unit.md` completed:** design, gate results with job IDs, the FMA experiment table and decision, and lessons.
- Updates to `notes/olaf-unit4.md`, plus `notes/architecture.md` and the skill if something generalizes.
- **A final report, sent ONCE**, in no more than about 35 lines:
  - gates with numbers;
  - the FMA experiment table and which variant the rule picked;
  - the speedups: solve and whole run, CPU;
  - the files changed;
  - the job IDs;
  - open issues.
- If you run long, write a resume checkpoint into the card before stopping.
