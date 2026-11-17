<!-- Brief given to subagent aed8806d (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are porting the first real OpenFAST unit to C++/Kokkos in the repo `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. Nothing of OLAF has been ported yet; you are writing the first unit.

## Goal
Port the OLAF **direct Biot–Savart kernels** in `deps/openfast/modules/aerodyn/src/FVW_BiotSavart.f90`:
- `ui_seg` (vortex segments)
- `ui_part_nograd` (vortex particles)

Put them behind their frozen Fortran interfaces, and prove equivalence with L0 differential tests on Kokkos-OpenMP and Kokkos-CUDA (H100). Unit name: `biotsavart-direct`.

## Read first (they define how to do this)
1. `CLAUDE.md`, then follow `.claude/skills/ofk-port-unit/SKILL.md` step by step. That skill is your procedure.
2. `notes/olaf-placement.md`:
   - §3 kernel shapes: `ui_seg` is used both "square" (wake on wake) and "skinny" (few lifting-line targets, many sources). Square is the month-1 shape; keep class B there.
   - §6 audit list: every item that touches these two routines needs an L0 case. Known traps:
     - `MINNORM = 1e-4` is a default-kind (single-precision) literal, so it equals 9.9999997473787516e-05 in double;
     - regularization `select case` branches, each with its own guard style;
     - `PRECISION_UI = epsilon/100`;
     - accumulation into `Uind`;
     - "3 × n++" over-sized arrays: touch only the first n columns;
     - 1-based `SegConnct`.
   - Run `scripts/literal-kinds.py deps/openfast/modules/aerodyn/src/FVW_BiotSavart.f90` and `scripts/kernel-inventory.py` on that file yourself.
3. `notes/criteria.md`: gates A*, R*, P*, and equivalence classes B, R1, R2.
4. `notes/build-kokkos-fortran.md` §6: boundary rules. `notes/harness.md`: oracle, `nm` check, unit-card template in §0.
5. `examples/fortran-kokkos/`: the verified template (frozen ref, `_ref`-renamed oracle via `rename_modules.py`, shim, Kokkos kernel, L0 test, C-binding smoke test). Copy its structure.

## Decisions already made (do not re-decide)
- **Kokkos 5.3.0, C++20.**
  - Installs: `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/kokkos-5.3.0-omp` (Serial+OpenMP) and `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/kokkos-5.3.0-cuda-h100` (CUDA sm_90 + OpenMP + Serial).
  - For CUDA, compile with `<install>/bin/nvcc_wrapper` and `NVCC_WRAPPER_DEFAULT_COMPILER=g++`.
  - Use `find_package(Kokkos 5.3...<6 REQUIRED)`.
  - The template's CMake still says 4.5; I am switching it in parallel. Don't edit `examples/`.
- Toolchain: `source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0`.
  - **Compile everything (CUDA included) on the login node.** Only *running* on a GPU needs a job.
- `ReKi` comes from the `OPENFAST_DOUBLE_PRECISION` macro. Target double precision; the single-precision build should still compile.

## Where files go
`deps/openfast` is the frozen baseline: **never edit it**. The fork that will hold the real port is still undecided, so use these locations:
- `modules/kokkos/src/biotsavart_direct.cpp` (+ header if needed), and a `modules/kokkos/CMakeLists.txt` building `ofk_kernels` (one library for all units; reuse the lazy-init runtime pattern from `examples/fortran-kokkos/src/ofk_runtime.*`).
- `port/modules/aerodyn/src/FVW_BiotSavart.f90`: a copy of the pin's file where only the bodies of the two routines change, to `#ifdef OPENFAST_KOKKOS` → `call ofk_...` / `#else` original body. Same module name and public interface.
- `tests/biotsavart-direct/`:
  - `unit.md` (unit card: contract, equivalence class per output, audit list with L0 case names, gate results with evidence);
  - L0 tests;
  - a CMakeLists that builds:
    - the oracle from the **pin's** file renamed `_ref`;
    - the port;
    - the tests.
- Build trees go outside the repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/biotsavart-direct-{omp,cuda}`.
  - `FVW_BiotSavart.f90` depends on NWTC_Library (and maybe FVW_Types). Build those from the pin's sources in a build tree outside the repo, e.g. configure OpenFAST's own CMake from `deps/openfast` into `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/openfast-pin` and build only the needed library targets. Link against them; do not rename shared modules.
- Another agent session may be editing shared docs at the same time. Only write the files listed above. If you think a shared file (`notes/*`, `README.md`, `CLAUDE.md`, skills, `examples/`) should change, say so in your report instead of editing it.

## Gates to run and report
- **L0 R1:**
  - ref-vs-ref = exactly 0 first;
  - then port-vs-ref on Kokkos-OpenMP (8 threads, seeds 1 and 7, a few hundred cases);
  - then on CUDA.
  - Report the worst relative error per backend and per regularization branch, and whether each is bitwise (class B) or needs R1.
- **GPU run:** submit with `sbatch -p h100 --gres=gpu:1 -c 8 -t 00:15:00 -J ofk-biotsavart-direct-l0`, logging to the build dir.
  - The PI-wide limit is 16 GPUs; check `squeue` first and don't block on `srun`.
  - If the job is still pending when everything else is done, report the job id and log path instead of waiting indefinitely.
- **API:** the shim keeps the module's public interface. Show that the pin's callers of these routines compile unchanged against the port module (e.g. compile the pin's `FVW_Subs.f90`/`FVW_Wings.f90` or whichever files call them against your port module; use `scripts/callers.sh` to find them).
- **P (quick, not the full perf gate):** one timing per backend vs the Fortran (OpenMP) at, say, 10⁴ and 10⁵ targets×sources, noting the host↔device copy cost per call.

## Don'ts
- No git commits.
- Don't loosen equivalence classes without stating why in `unit.md`.
- Don't port the treecode, `InducedVelocitiesAll_*`, or the wake residency yet. Note what you learned that affects them.

## Report back (concise)
- Files created.
- Gate table (ID, PASS/FAIL, evidence).
- Worst L0 errors per backend.
- Timings.
- Slurm job ids.
- Problems or open questions.
- Suggested edits to shared docs.
- Your rough wall-clock time spent.
