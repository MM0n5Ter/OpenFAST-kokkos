<!-- Brief given to subagent a6a51a5d (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are restructuring the OpenFAST → C++/Kokkos port in `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. The user wants a **clean composition instead of a full fork**. This is infrastructure work: no new physics, and **no change in behaviour**. Proving "no change" is the main deliverable.

## Read first
1. `CLAUDE.md` and `notes/architecture.md` (§1 describes the old fork layout, which you are replacing).
2. `notes/build-kokkos-fortran.md` §5 and the "Verified in the fork (unit 1)" subsection.
3. The unit cards `tests/biotsavart-direct/unit.md` and `tests/induced-velocities/unit.md`: they hold the build and test recipes and the evidence you must reproduce.
4. `.claude/skills/ofk-port-unit/SKILL.md`, toolchain section: AMD GPU nodes, portable OpenBLAS, shared-tree rules.

## Today's layout (to be replaced)
- `deps/openfast`: the upstream pin `12fa3fc`, frozen. Never edit it.
- `fork/`: a local submodule (`MM0n5Ter/openfast`, branch `kokkos`, never pushed, no commits), with **uncommitted** changes (`git -C fork status`, `git -C fork diff --stat`):
  - `CMakeLists.txt` (+15): `OPENFAST_KOKKOS` option, `find_package(Kokkos 5.3...<6)`, `add_subdirectory(modules/kokkos)`;
  - `modules/aerodyn/CMakeLists.txt` (+7): `aerodynlib` gets `OPENFAST_KOKKOS` and links `ofk_kernels`;
  - `modules/aerodyn/src/FVW.f90`, `FVW_BiotSavart.f90`, `FVW_Subs.f90`: shims under `#ifdef OPENFAST_KOKKOS`;
  - new `modules/kokkos/` (`ofk_kernels`, with `kokkos_compilation(GLOBAL)` for CUDA, `install(EXPORT OpenFASTLibraries)`, and the `OFK_FAST_*` options).
- Unit tests in `tests/biotsavart-direct/` and `tests/induced-velocities/` find the port through a CMake variable `OFK_FORK`.

## Target layout
- `src/modules/aerodyn/src/FVW.f90`, `FVW_BiotSavart.f90`, `FVW_Subs.f90`: **full copies of only the files we change**, at OpenFAST's own relative paths. Any file not in `src/` comes from `deps/openfast`.
- `src/modules/kokkos/`: the `ofk_kernels` library (CMakeLists and sources).
- `src/cmake/ofk_overlay.cmake`: injected into OpenFAST's **unmodified** CMake with `-DCMAKE_PROJECT_INCLUDE=<repo>/src/cmake/ofk_overlay.cmake`. OpenFAST requires 3.15 and the cluster has cmake 3.31, so `cmake_language(DEFER)` is available. It must:
  - define the `OPENFAST_KOKKOS` option (default OFF) and the `OFK_*` options;
  - when ON: `find_package(Kokkos 5.3...<6)` and `add_subdirectory(<repo>/src/modules/kokkos ${CMAKE_BINARY_DIR}/ofk_kernels)`;
  - with a deferred call at the end of the top-level directory: replace in `aerodynlib`'s `SOURCES` every `modules/<mod>/src/X.f90` that has an overlay in `src/` (discover them by globbing `src/modules/*/src/*`; don't hard-code the list), add the `OPENFAST_KOKKOS` compile definition (PRIVATE) and link `ofk_kernels`.
  - The overlay replacement should apply even with the option OFF. The `#else` branches are the original bodies, so OFF must reproduce the pin exactly; that is the A5 check.
  - Handle the install/export-set requirement that unit 1 hit.
- `tools/openfast-build.sh <config> <builddir>` (or similar): one command that configures `deps/openfast` with the overlay for the configs `off | omp | fomp | cuda-h100 | cuda-a100`. Use the exact flags from the unit cards: OpenBLAS, `DOUBLE_PRECISION`, `OPENMP` on or off, Kokkos roots, `CMAKE_CXX_COMPILER=g++` for CUDA. A plain pin build without the overlay remains the baseline.
- Unit tests: replace `OFK_FORK` with the overlay. Port sources come from `src/` where overlaid, and from `deps/openfast` otherwise. The tests must still build their `_ref` oracle from `deps/openfast`.
- **r-test:** the fork had a hardlink copy of r-test at `fork/reg_tests/r-test`. With the overlay, the build's source dir is `deps/openfast`, whose `reg_tests/r-test` is already initialized. Check that ctest's OLAF cases run from an overlay build.

## Also: implement the FMA decision (user, 2026-11-14: keep FMA)
- The chain-claim-3 bound is declared as `max(1e-12, 1.1 × the pin Fortran's own FMA-vs-no-FMA spread on the same inputs and seed)`, per regularization branch.
- Today `l0_ui_seg_strictfast` and `l0_ui_seg_fast` fail at 1e-12 (worst 2.99e-12 on seed 1, 1.1e-11 on seed 7); the Fortran's own `l0_ui_seg_ffma` spread is the same.
- Implement that bound in the unit 1 tests (and unit 2's if it applies), e.g. as a ctest fixture that runs the `_ffma` measurement first and passes the per-branch spread on. Keep it simple and explicit.
- Afterwards all unit 1 ctests should pass.

## Constraints: don't break running work
- **Unit 3 (treecode) is in progress in another agent**, which has been asked to checkpoint. It created `fork/modules/kokkos/src/treecode_*` and `tests/treecode/*`, and uses build trees `opt/build/fork3-*` and `u3-*` plus Slurm jobs named `ofk-treecode-*`.
  - Migrate everything else first.
  - Before you copy the treecode files or retarget `tests/treecode/`, wait until `tests/treecode/unit.md` contains a "Status / resume" section **and** no `ofk-treecode` build process is running on the login node (`ps -u $USER`). Poll every few minutes; don't busy-wait.
  - Then sync those files into `src/modules/kokkos/` and retarget `tests/treecode/CMakeLists.txt` and its scripts. Building and testing unit 3 is **not** your job, beyond checking that it configures.
- **Queued GPU jobs 1860007, 1860223 and 1860779** (h100/h200, start ~2026-11-16) use the existing trees `opt/build/fork-cuda`, `iv-cuda`, `biotsavart-direct-v1-cuda`, and job 1860007 reads `fork/reg_tests`.
  - Do not modify, rebuild or delete those trees.
  - **Do not delete the `fork/` directory.** Remove it from git's view instead: deinit/remove the submodule registration (`.gitmodules`, `.git/config`, the index), keep the files on disk, and add `fork/` to `.gitignore` with a comment that it is deleted after those jobs run.
- Use new build trees: `opt/build/ov-{off,omp,fomp,cuda-h100,cuda-a100}` and `opt/build/ov-tests-*`. Compile on the login node, with at most `-j 12` alongside other users.
- `deps/openfast`: never edit.
- No git commits or pushes.
- Toolchain: `source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28`. Kokkos 5.3.0 installs are in `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/` (`-omp`, `-cuda-h100`, `-cuda-a100`).
- GPU runs: `interactive_gpu` (A100 MIG, no queue). Put `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib` first on `LD_LIBRARY_PATH`. CPU runs: `med`, `-c 8`. Job names: `ofk-infra-<what>`.

## Prove "no change in behaviour" (record this in a new `notes/overlay.md`)
1. **A5:** the overlay build with OFF is identical to a plain pin build: same `.mod` md5s, and `.text` differs at most in line-number immediates, as unit 1 checked.
2. Unit 1 and unit 2 L0 suites pass on overlay-built ports, on CPU OpenMP and on CUDA A100 MIG, with the **same numbers** as recorded in their unit cards. Strict stays bitwise; the fast worst errors are the same.
3. ctest `aerodyn_utest|OLAF` on `ov-off`, `ov-omp` strict and `ov-omp` fast gives the same pass/fail as the fork builds. Strict OLAF outputs equal OFF.
4. One full-precision BAR case (from `opt/build/iv-cases`) in strict and fast on `ov-omp` gives data lines identical to the old `fork-omp` outputs for the same mode. Find them under `opt/build/iv-r8-1860700` or `iv-perf*`; if none match exactly, rerun on `fork-omp` for comparison.

## Report back (concise)
1. New layout tree, and how to build (the exact commands).
2. The four "no change" results with evidence.
3. The FMA-bound implementation and the unit 1 ctest results.
4. What remains in `fork/` and when it can be deleted.
5. Problems, and suggested edits to shared docs: `architecture.md` §1, `CLAUDE.md`, `README.md`, the skills, `build-kokkos-fortran.md`. **Don't edit those yourself**; the parent will.
6. Wall-clock time.
