<!-- Brief given to subagent a103544e (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are finishing **unit 3 (`treecode`)** of the OpenFAST-kokkos port: OpenFAST's OLAF particle treecode (`grow_tree_part`, `ui_tree_part`) ported to C++/Kokkos with strict and fast modes. Repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. Read these first, in order: `CLAUDE.md`, `notes/architecture.md`, `notes/criteria.md`, `notes/overlay.md`, `.claude/skills/ofk-port-unit/SKILL.md`, and then **`tests/treecode/unit.md`**. That card is your resume point: its "Status / resume" section lists what is done with evidence, what is half-done, and the exact next steps.

## What changed since that checkpoint
- The `fork/` layout is retired. The code now lives in the repo as an overlay on the pin: `src/modules/aerodyn/src/FVW*.f90` and `src/modules/kokkos/` (the treecode is `treecode_*.{hpp,cpp}` there). Build OpenFAST with `tools/openfast-build.sh <pin|off|omp|fomp|cuda-a100|cuda-h100> <builddir>`. Build test trees from `tests/treecode/CMakeLists.txt` (see `notes/overlay.md` §3 and `src/cmake/ofk_sources.cmake`). Wherever the card says `fork3-*`, use overlay builds instead. Never use or edit `fork/`, and never edit `deps/openfast`.
- These overlay trees already exist, built from the current code (unit 3 checkpoint included; Kokkos treecode is the default, `OFK_TC=fortran` selects the Fortran tree): `opt/build/ov-off`, `ov-omp`, `ov-fomp`, `ov-pin`, `ov-cuda-a100`, `ov-cuda-h100`, `ov-tests-iv-omp`, `ov-tests-tc-omp`. They are under `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build`.
- Overlay system level already done (job 1863342, `notes/overlay.md` §4.3–4.4): ctest `aerodyn_utest|OLAF` 17/17 on ov-off and ov-omp in strict, fast and fast with `OFK_TC=fortran`. GPU ctest 11/11. BAR wall times on 2 threads: tree case 97.9 s (off), 88.5 s (fast with `OFK_TC=fortran`), 56.2 s (fast with the Kokkos tree). Do not redo this. Use it.
- A100 timing sweep of `perf_tc` with the host-side tree build (job 1862690): `opt/build/u3-perf-a100-1862690/sweep.csv` and `prof.*.log`.

## Your scope: finish unit 3
1. **System level on CPU:**
   - claim 1 at system level: strict vs the `ov-off` text output on `opt/build/iv-cases/bar-r8-tree`;
   - R8 against base₁ `opt/build/iv-r8-1860700/bar-r8-tree.fomp1/ad_driver.out` with the FMA envelope `opt/build/iv-r8ffma-1860993/bar-r8-tree.ffma/ad_driver.out` (`tests/induced-velocities/r8.py`). Also report the OFF-vs-OMP drift, as the card notes;
   - P2 at system level: `ov-fomp` vs `openfast-pin-fomp` at 8 threads, with an `OFK_TC=fortran` A/B;
   - P6 on a full run.
   
   Use the `med` partition, which has a queue, but a short one. `tests/treecode/cpu.sbatch` probably needs its paths moved from fork3 to the overlay.
2. **Fast-mode tree build on the device.** The design is already decided in the card, "Exact next steps" item 3: level-synchronous K1–K5, no thrust sort, host build kept for host backends.
   - Implement it, then gate it with `l0_tc` R2, claim 3 on the GPU, R7 and P6 (P6 counts only allocations with bytes > 0).
   - Run correctness on `interactive_gpu` (A100 MIG 2g.20gb, no queue, 4 h limit, `-c 4`).
3. **Mutation checks** (card item 4), then A1/A2/A5/A7 on the overlay builds. For A5 see `tests/overlay/a5.sh`.
4. **GPU timing.** Full GPUs queue for a day or more:
   - `-p a100` for sm_80 trees, `-p h100,h200` for sm_90 trees;
   - PI-wide limit of 16 GPUs;
   - name jobs `ofk-treecode-<what>`.
   
   Once the device build passes, submit one timing sweep per architecture with the new binaries. Use `tests/treecode/gpu.sbatch perf`, pointed at your new trees. Submit and move on: do not wait for these jobs to start.
5. **The segment tree** (VelocityMethod 4, not the default) is **out of scope**. Leave its section in the card as deferred.

## Hard rules
- **Do not rebuild or reconfigure any tree that a queued job uses.** Job 1875922 (A100) uses `ov-cuda-a100`, `ov-tests-iv-cuda-a100`, `u3-cuda-a100` and `openfast-pin-fomp`. Job 1862732 (H100) uses `u3-cuda-h100`. Job 1860223 uses `fork-cuda`, `iv-cuda` and `iv-omp`.
  - Once you change the C++ source, build into **new** directories, e.g. `u3d-*` and `ov3-*`.
  - Note that 1875922 compiles nothing. It runs the existing binaries, so source edits are safe; only rebuilds of those trees are not.
  - Do not cancel other people's jobs, or any job other than ones you submitted.
- Compile everything, CUDA included, on the login node with `module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0` and at most `-j 12`.
  - GPU and full-OpenFAST runs on AMD nodes need `export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH`.
  - `srun` does not work from this shell. Use `sbatch`, or `sbatch --wait` for short jobs.
  - Never edit a bash script while a job is running it. Copy it first.
- Compare text outputs from line 8 on (line 7 and earlier hold a timestamp). The Python environment is `opt/venv-rtest`.
- Do **not** git commit or push. Do not write the cost ledger. The parent session reviews and commits.
- If something in the design turns out to be wrong, or a gate cannot pass, stop working around it. Record the evidence in the card and report it.

## Deliverables
- **`tests/treecode/unit.md` updated.** Replace the checkpoint "Status / resume" section with final status. Fill "Gate results" with every gate, its number, and the job or log that backs it. Extend "Lessons" with what is generalizable for later units, and keep it short.
- If you learn something that every unit needs, add it to `.claude/skills/ofk-port-unit/SKILL.md` and/or `notes/architecture.md` (§7 row for unit 3).
- **Final report, sent ONCE**, no more than about 40 lines:
  - gates, pass or fail, with numbers;
  - the speedups: CPU kernel, CPU system and GPU kernel, as far as measured;
  - the files you changed;
  - the Slurm job IDs you submitted, with what each one is for and which are still pending;
  - open issues.
  
  An earlier subagent kept re-sending its report in a loop. Send yours once and stop.
- If you run long, write a resume checkpoint into the card before you stop, in the same format as the current one.
