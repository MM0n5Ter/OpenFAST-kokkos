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
   - Put the shared body in `fork/modules/kokkos/src/<unit>_body.hpp`, and instantiate it in `<unit>_strict.cpp` (FMA off) and `<unit>_fast.cpp` (FMA on).
   - Select the mode at runtime with `OFK_MODE`. Strict must be class B against the Fortran on the CPU; fast must reach the declared class against strict.
   - Wrap the body in a mode namespace (ODR, architecture.md §3). Follow the CPU vectorization recipe there, and check it with `-fopt-info-vec`.
   - Measure the Fortran's own FMA spread (`l0_*_ffma`) before declaring the fast-mode bound.
   - Template: `fork/modules/kokkos/src/biotsavart_direct_*` and `tests/biotsavart-direct/`.
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
4. **Shim.** In the fork's Fortran file (`fork/modules/<mod>/src/`, edited in place), keep the module name and public interface, and change only the body:
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
- CUDA fork build: keep `CMAKE_CXX_COMPILER=g++`; `kokkos_compilation(GLOBAL)` is already in `fork/modules/kokkos`. GPU jobs: use `-p h100,h200` (both sm_90).
- r-test:
  - it is initialized in `deps/openfast`;
  - the fork's copy is a hardlink copy (`cp -al deps/openfast/reg_tests/r-test/. fork/reg_tests/r-test/`, since `--reference` refuses a shallow source);
  - use the Python at `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest` (numpy, bokeh, pyOpenFAST);
  - `py_ad_*` cases assume the build tree is `<src>/build`, so they fail out of tree on baseline and port alike;
  - compare `.outb` channel data with `tests/biotsavart-direct/cmp_outb.py`.
- The Noether `pack4` install is static and non-PIC, so it cannot go into the C-binding `.so`. Do not use it for OpenFAST.
- **Unit 2 lessons on running jobs:**
  - P2 compares against a fork built with `OPENMP=ON` plus Kokkos-OpenMP.
  - Never run concurrent jobs with `OMP_PROC_BIND=spread`; they pin to the same cores.
  - `bc` is not installed.
  - Every new translation unit with device arithmetic must be listed as strict or fast in `fork/modules/kokkos/CMakeLists.txt`, because nvcc uses FMA by default.
  - R8 and claim 1 at system level need full-precision text output, not `.outb`.
- **GPU nodes are AMD EPYC, and the cluster OpenBLAS is a zen4 (AVX-512) build.**
  - The A100 nodes (EPYC 7443, Zen3) die with "Illegal instruction" in `dgemm_` on any full OpenFAST run.
  - In every GPU job script, after `module load`, put this first on the library path: `export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH`. That is our build: LP64, DYNAMIC_ARCH, Haswell baseline.
  - The login node and `med` are Intel Sapphire Rapids, so a test that passes there proves nothing about the GPU nodes' CPU side.
- **Shared fork tree:** the fork's working tree is shared by every unit in flight. Never rebuild another unit's build tree, and never let a queued job's build tree be reconfigured. Each unit uses its own build dirs.
- **Fast GPU correctness checks:** `interactive_gpu` has A100 MIG slices with no queue and a 4 h limit. Use `kokkos-5.3.0-cuda-a100` and the `*-a100` build trees (`opt/build/build-a100.sh`). It is not for timing.
- CPU timing jobs: use the `med` partition, because `interactive_cpu` caps a user at 4 CPUs. The h100 queue can wait hours, so submit GPU runs early.
- Login nodes have no GPU. GPU jobs share a PI-wide limit of 16 GPUs: check `squeue`/`stree --limits` before submitting, and use `sbatch` rather than a blocking `srun` when the group is at its limit.

## Never
- Edit `deps/openfast`.
- Edit a caller to make a port fit (that fails A2).
- Loosen a declared equivalence class without a sign-off in the unit card.
- Port a unit's tests in the same step as the unit (testing.md: the oracle stays one step behind).
