---
name: ofk-port-unit
description: Port one OpenFAST Fortran unit (an API contract such as ui_seg, InducedVelocitiesAll_*, a FlowField routine) to C++/Kokkos behind its frozen Fortran interface, with L0 differential tests and the API/Regression/Performance gates. Use after ofk-analyze-unit has produced the placement note, when asked to "port X to Kokkos", "write the Kokkos kernel for X", or "add the L0 test for X".
---

# Port a unit to Kokkos behind its frozen interface

Template to copy: `examples/fortran-kokkos/`. It is verified on Kokkos 5.3.0 OpenMP: L0 206/206 bitwise, `capi_smoke` OK. On the H100: worst relative error 2.3e-18.
The rules at the Fortran↔C++ boundary are in `notes/build-kokkos-fortran.md` §6. Follow them all.
Gates: `notes/criteria.md`. Harness: `notes/harness.md`.

## Inputs you need before writing code
- The unit's placement note (`notes/<module>-placement.md`). If it is missing, run the `ofk-analyze-unit` skill first.
- The unit card `tests/<unit>/unit.md`, with:
  - the contract (routines and signatures);
  - the declared equivalence class per output (B / R1 / R2);
  - the audit list.

## Steps (one unit = one branch = one PR)

1. **Freeze the contract.** Copy the routine signatures from the pin into the unit card. Mark intent(inout) arrays that are accumulated into, over-sized arrays ("n++"), and 1-based index arrays.
2. **L0 first (it is the oracle).**
   - Build `oracle_<unit>` from the **pin's** files with `examples/fortran-kokkos/rename_modules.py` (`X` → `X_ref`). Do not rename shared modules (`NWTC_Library`, `*_Types`).
   - Run the `nm` check (harness.md §1.4).
   - Write `tests/<unit>/l0_<routine>.f90` from `examples/fortran-kokkos/tests/l0_ui_part.f90`:
     - a fixed edge case per audit item;
     - every `select case` value;
     - sentinel tails and nonzero initial outputs;
     - non-contiguous sections;
     - then randomized cases (`OFK_L0_SEED`, `OFK_L0_CASES`, `OFK_L0_RTOL`).
   - Run it ref-vs-ref first. It must give exactly 0.
   - The class-B oracle is the pin built with **OPENMP=OFF** and its exact flags (`-fdefault-real-8 -fdefault-double-8`, per-file `-fno-trapping-math`). With OpenMP on, the Fortran's `!$OMP SIMD` reductions reorder its own sums.
   - Check L0 sensitivity by mutation: break each guard or constant by hand, and confirm L0 fails.
   - To link the pin's NWTC library: `module load openblas/0.3.28`, `-DBLA_VENDOR=OpenBLAS`. Worked example: `tests/biotsavart-direct/`.
3. **Kernel.** Follow `notes/architecture.md` §2–3:
   - Put the shared body in `src/modules/kokkos/src/<unit>_body.hpp`, and instantiate it in `<unit>_strict.cpp` (FMA off) and `<unit>_fast.cpp` (FMA on).
   - Select the mode at runtime with `OFK_MODE`. Strict must be class B against the Fortran on the CPU; fast must reach the declared class against strict.
   - Wrap the body in a mode namespace (ODR, architecture.md §3). Follow the CPU vectorization recipe there, and check it with `-fopt-info-vec`.
   - Measure the Fortran's own FMA spread (`l0_*_ffma`) before declaring the fast-mode bound.
   - Template: `src/modules/kokkos/src/biotsavart_direct_*` and `tests/biotsavart-direct/`.
   - Both modes must be deterministic.

   Also:
   - `extern "C" ofk_<unit>_<routine>`; `ReKi` from `OPENFAST_DOUBLE_PRECISION`; `LayoutLeft` views of extent `(3,n)`.
   - Constants reproduced bit for bit. **Pass flag-dependent constants (default-kind literals such as `MINNORM`, `epsilon`-based ones) from the Fortran shim.** Do not hard-code them: OpenFAST's double build uses `-fdefault-real-8`, so `1e-4` is a double there.
   - Compile `ofk_kernels` with `-fno-math-errno -fno-trapping-math`. These are gfortran's defaults and do not change values.
   - nvcc rejects an extended lambda whose first capture of a variable is inside `if constexpr`. Capture before the branch.
   - One functor or template instantiation per regularization branch. Keep the guards literally.
   - Same per-target sum order for class B.
   - Pick the pattern from the placement note (square: `RangePolicy`; skinny: `TeamPolicy` reduce, class R1; tree: flat arrays and an explicit stack).
   - Persistent device buffers live in a per-instance handle (`type(c_ptr)` in a unit-internal Fortran type). No C++ globals per instance. Catch every exception.
4. **Shim.** In the overlay copy of the Fortran file (`src/modules/<mod>/src/<file>`). If the file is not there yet, copy it from `deps/openfast/modules/<mod>/src/` first. The overlay picks it up automatically (`notes/overlay.md`)., keep the module name and public interface, and change only the body:
   `#ifdef OPENFAST_KOKKOS` → `call ofk_<unit>_<routine>(...)`, `#else` → the original body.
5. **Run the gates** (`tools/gate.sh <unit> api|reg|perf` once it exists; until then the commands in harness.md):
   - **R1** on Serial/OpenMP (login node) and CUDA (`srun -p h100 --gres=gpu:1`; the cluster Kokkos CUDA build is sm_90);
   - **A1–A6**;
   - **R6** (`aerodyn_utest`, `ctest -R OLAF`);
   - **P2** (no CPU slowdown) and **P4** (small cases).
6. **Record** in the unit card: PASS/FAIL per ID with evidence, worst L0 error per backend, timings, and E1–E6.
7. **Record cost (money and time)** before ending every session, including sessions that end in failure: failed attempts are part of the unit's cost.
   - Name every Slurm job `ofk-<unit>-<what>`.
   - Then run:
     ```bash
     tools/cost/session_cost.py --session <id-prefix> --ledger tools/cost/ledger.csv --unit <unit> --phase port|debug \
         --note "<one line>" --gpu-h <from tools/cost/compute_hours.sh> [--human-min N]
     ```

## Toolchain on this cluster
- `module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0`.
- Kokkos 5.3.0 installs (shared, C++20) under `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/`:
  - `kokkos-5.3.0-omp` (Serial + OpenMP);
  - `kokkos-5.3.0-cuda-h100` (CUDA sm_90 + OpenMP; compile with its `bin/nvcc_wrapper`, `NVCC_WRAPPER_DEFAULT_COMPILER=g++`).
- Compile everything, CUDA included, on the login node. Only GPU *runs* need a job.
- OpenFAST builds: `tools/openfast-build.sh pin|off|omp|fomp|cuda-a100|cuda-h100 <builddir>` (login node). CUDA keeps `CMAKE_CXX_COMPILER=g++`; `kokkos_compilation(GLOBAL)` is in `src/modules/kokkos/ofk_kernels.cmake`. GPU partitions:
  - `-p h100,h200` (both sm_90; `cuda-h100` trees);
  - `-p a100` (sm_80; `cuda-a100` trees).
  - Both queues can wait hours to days. `sbatch --test-only` estimates are unreliable.
- r-test:
  - it is initialized in `deps/openfast`;
  - overlay builds use `deps/openfast/reg_tests/r-test` directly;
  - use the Python at `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest` (numpy, bokeh, pyOpenFAST), plus `vtk`, which `py_ad_B1n2_OLAF` needs;
  - build trees under `opt/build/` need the symlink `opt/build/reg_tests/r-test -> deps/openfast/reg_tests/r-test`, because the pyOpenFAST tests import `driver_utilities` from `<build>/../reg_tests/r-test/modules`. It already exists.
  - compare `.outb` channel data with `tests/biotsavart-direct/cmp_outb.py`.
- The Noether `pack4` install is static and non-PIC, so it cannot go into the C-binding `.so`. Do not use it for OpenFAST.
- **Unit 2 lessons on running jobs:**
  - P2 compares against the `fomp` build (Fortran `OPENMP=ON` plus Kokkos-OpenMP), and against `openfast-pin-fomp` on the same node.
  - Never run concurrent jobs with `OMP_PROC_BIND=spread`; they pin to the same cores.
  - `bc` is not installed.
  - Every new translation unit with device arithmetic must be listed as strict or fast in `src/modules/kokkos/ofk_kernels.cmake`, because nvcc uses FMA by default.
  - R8 and claim 1 at system level need full-precision text output, not `.outb`.
- **GPU nodes are AMD EPYC, and the cluster OpenBLAS is a zen4 (AVX-512) build.**
  - The A100 nodes (EPYC 7443, Zen3) die with "Illegal instruction" in `dgemm_` on any full OpenFAST run.
  - In every GPU job script, after `module load`, put this first on the library path: `export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH`. That is our build: LP64, DYNAMIC_ARCH, Haswell baseline.
  - The login node and `med` are Intel Sapphire Rapids, so a test that passes there proves nothing about the GPU nodes' CPU side.
- **Shared source tree:** `src/` is shared by every unit in flight. Never rebuild another unit's build tree, and never let a queued job's build tree be reconfigured. Each unit uses its own build dirs.
- **Fast GPU correctness checks:** `interactive_gpu` has A100 MIG slices with no queue and a 4 h limit. Use `kokkos-5.3.0-cuda-a100` and the `*-a100` build trees (`opt/build/build-a100.sh`). It is not for timing.
- CPU timing jobs: use the `med` partition, because `interactive_cpu` caps a user at 4 CPUs. The h100 queue can wait hours, so submit GPU runs early.
- Login nodes have no GPU. GPU jobs share a PI-wide limit of 16 GPUs: check `squeue`/`stree --limits` before submitting, and use `sbatch` rather than a blocking `srun` when the group is at its limit.

- **Unit 3 lessons (device algorithms, mutation checks):**
  - Write a device algorithm generic over the execution space, with team size 1 on host backends, and give it a runtime switch that forces it on the CPU (`OFK_TC_DEVBUILD=1`). Then the whole CPU L0 tests it on the login node in seconds, before any GPU job.
  - Give every device rewrite of a host algorithm a self-check switch that runs both and compares the exact parts bit for bit (integers, copies, min/max; `OFK_TC_DEVBUILD_CHECK=1`). It turns "R2 passes" into "same structure", and it is what catches the device-only mutations: R2 alone does not.
  - Determinism without atomics: a fixed pairwise tree over team scratch (fixed team size) for float sums; integer scans for every partition and offset. Stable ranks inside a team tile come from one exclusive team scan of byte-packed `uint64` octant counts (tile <= 255).
  - Run mutation checks on a **copy** of `src/` with its own test tree (`tests/treecode/` + `examples/fortran-kokkos/`, `-DOPENFAST_SOURCE=<repo>/deps/openfast`), never on `src/` itself: another build may be reading `src/` at the same time. Template: `tests/treecode/mutate.py`.
  - Slurm stores a copy of the batch script at submission, so editing it under a *queued* job is safe; helper scripts the job calls at runtime (`run_case.sh`, `r8.py`) are not.
  - When queued jobs hold your build trees, build the new code into prefixed trees (`U3=u3d tests/treecode/build_u3.sh …`, `ov3-*`) and point job scripts at them through environment variables, not by editing paths.

## Never
- Edit `deps/openfast`.
- Edit a caller to make a port fit (that fails A2).
- Loosen a declared equivalence class without a sign-off in the unit card.
- Port a unit's tests in the same step as the unit (testing.md: the oracle stays one step behind).
