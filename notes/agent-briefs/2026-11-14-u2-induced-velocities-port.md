<!-- Brief given to subagent a76bb357 (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are doing **unit 2 of 3** in the OpenFAST → C++/Kokkos port. Repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. The three units establish the port's architecture, and the user wants us to **learn something from each unit**, so your lessons section matters as much as the code.

## Read first
1. `CLAUDE.md`.
2. **`notes/architecture.md`**, all of it. It now includes unit 1's lessons: the ODR namespace rule, the FMA finding, the CPU vectorization recipe and the build facts.
3. `notes/criteria.md`.
4. `.claude/skills/ofk-port-unit/SKILL.md`.
5. `notes/olaf-placement.md`: §1 step anatomy, §2 boundary, §4 residency, §5 fusions (F1, F5), §6 audit items 4, 5, 9, 11, 14.
6. Unit 1's card `tests/biotsavart-direct/unit.md`, and its code in `fork/modules/kokkos/src/biotsavart_direct_*`, which you build on.
7. Unit 1's integration into the fork's CMake (`git -C fork diff`).

## State
- `fork/` is the submodule `MM0n5Ter/openfast`, local branch `kokkos`, at pin `12fa3fc`, with **uncommitted** unit-1 changes. Leave them uncommitted, and keep them working.
- `deps/openfast` is the frozen pin. **Never edit it.**
- Full fork builds exist under `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/fork-{off,omp,cuda}`.
- r-test is initialized; read the skill's toolchain section for the r-test setup and the Python venv.
- **Unit 1 results:**
  - strict mode is bitwise to the Fortran;
  - fast mode is 1.07–3.95× Fortran+OpenMP on CPU and up to about 120× on H100;
  - strict-mode r-test `.outb` data is identical to the OFF build.
  - **P6 failed**, because device memory is allocated on every call. Fixing it is your job.
- **Open user decision:** FMA in fast mode (`decisions.md`). Keep the current default (`OFK_FAST_FMA=ON`) and don't spend time on it. Note anything you see that bears on it.

## Unit 2: `induced-velocities`
Routines, in `fork/modules/aerodyn/src/FVW_Subs.f90` (pin line numbers):
- `InducedVelocitiesAll_Init` (:1237);
- `InducedVelocitiesAll_Calc` (:1285);
- the matching `_End` / cleanup;
- `LiftingLineInducedVelocities` (:1488);
- `WakeInducedVelocities` (:1379), as far as it calls the above;
- the packing they rely on (`PackPanelsToSegments` / `LatticeToSegments` / `SegmentsToPart` and similar).

**The treecode stays Fortran; it is unit 3.** With the tree velocity method, the tree path still runs in Fortran. With the direct methods, everything goes through Kokkos.

Goals, in priority order:
1. **A persistent per-instance device handle (fixes P6).**
   - Device buffers (sources, targets, outputs) live in a handle created at init or lazily, grown by doubling, and freed at end, with a finalize-hook fallback.
   - **Constraint:** registry types (`FVW_Types.f90`) are never hand-edited. Pick one:
     - (a) add a field via `FVW_Registry.txt` and regenerate with the fork's `openfast-registry`, which is allowed because it is not hand-editing;
     - (b) a module-level map keyed per instance.
   - Justify the choice in the unit card against these hazards:
     - glue code copying misc/state (`FVW_CopyMisc`, `FVW_CopyContState`);
     - re-evaluation and rollback at the same `t` (`FVW.f90:597-610`);
     - several rotors or instances in one process;
     - the C binding loaded from Python;
     - `*_End` never being called.
   - Show with Kokkos Tools (memory events) or counters that there is no allocation after step 3.
2. **The direct-method paths through `InducedVelocitiesAll_Calc` and `LiftingLineInducedVelocities` on Kokkos**, using one upload of the source set per step where the Fortran repacks every call. That is fusion F5, "LL upload-once"; declare its class.
   - Keep strict and fast modes. Strict must stay bitwise against the Fortran OPENMP=OFF build.
3. **End-to-end in the fork, CPU and GPU.**
   - Run the aerodyn_driver OLAF r-test cases and `aerodyn_utest` on fork-OFF, fork-OMP (strict and fast) and fork-CUDA (fast and strict). Strict must give `.outb` identical to OFF on CPU.
   - Most r-test OLAF cases use the treecode for the wake. So also run at least one case **with a direct wake method**: copy the case to a scratch dir outside the repo and change only its input `VelocityMethod`. That way the GPU actually carries the wake.
4. **R8 (trajectory drift)** on one longer OLAF case, e.g. `ad_BAR_OLAF` or HelicalWake, lengthened if needed:
   - first the Fortran baseline at `OMP_NUM_THREADS=1` vs `N`, which gives the envelope;
   - then fork fast on CPU and on GPU vs baseline₁;
   - report per channel the max drift and the step where it first exceeds 1e-6 relative.
   - This is the first system-level evidence for chain claim 4.
5. **F1, packing fused on device**, only if time remains. Otherwise write down precisely what it would take.

## Gates (record in `tests/induced-velocities/unit.md` with PASS/FAIL and evidence)
- **L0:** chain claims 1 and 3 for the new entry points, with every `VelocityMethod` and pack flag combination (audit #9), plus sentinel tails.
- **R6, R7, R8, P6.**
- **A1, A2, A5, A7.** A7 is copy/rollback/re-evaluation on the handle.
- **P1/P3:** per-step time split into kernel, H↔D and host, for the direct-method case on CPU and GPU.

## Rules
- Don't edit `deps/openfast`.
- **No git commits or pushes** anywhere.
- Toolchain: `source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28`. Kokkos 5.3.0 installs are in `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/`. Compile everything on the login node.
- GPU runs: `sbatch -p h100,h200 --gres=gpu:1 -c 8`, short time limits, job names `ofk-induced-velocities-<what>`. The PI-wide limit is 16 GPUs. The queue can take hours: submit as early as possible, then continue with other work. If a job is still pending when you are done, report its id and log path; don't idle.
- CPU timing jobs: the `med` partition with `-c 8`.
- Count `cpu_core_h` and `gpu_h` from `sacct` for your jobs (`tools/cost/compute_hours.sh`). Don't estimate.
- Write only in `fork/` and `tests/induced-velocities/`. You may add to `tests/biotsavart-direct/` only if unit 1's tests need adapting to your changes; say so if you do. Put suggested edits to shared docs in your report instead of making them.

## Report back (concise)
1. **Files and fork diff summary.**
2. **Gate table with evidence,** the handle design decision, and the R8 numbers.
3. **Timings.**
4. **Slurm jobs and log paths.**
5. **Lessons for unit 3 (treecode) and for the architecture:** what worked, what didn't, what you'd change in `architecture.md` or the skill, and what cost the most time.
6. **Suggested exact edits to shared docs.**
7. **Wall-clock time spent.**
