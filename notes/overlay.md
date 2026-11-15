# The overlay layout (replaces the `fork/` submodule, 2026-11-14)

The port is no longer a fork of OpenFAST. OpenFAST is built from the frozen pin `deps/openfast` (never edited), and a CMake overlay swaps in our copies of the few files we change and adds the C++/Kokkos library. This note says how it works, how to build, and the evidence that the restructuring changed no behaviour.

## 1. Layout

| Path | What |
|---|---|
| `deps/openfast` | The pin `12fa3fc`, frozen. Every file that is not in `src/` comes from here. Its `reg_tests/r-test` is initialized and used by ctest. |
| `src/modules/<mod>/src/<file>` | **Full copies of only the files we change**, at OpenFAST's own relative paths. Today: `src/modules/aerodyn/src/FVW.f90`, `FVW_BiotSavart.f90`, `FVW_Subs.f90`. Ported bodies sit under `#ifdef OPENFAST_KOKKOS`; the `#else` branches are the pin's bodies. |
| `src/modules/kokkos/` | `ofk_kernels`: `CMakeLists.txt` (stand-alone / test trees), `ofk_kernels.cmake` (the target; also `include()`d by the overlay), `src/*.{hpp,cpp}`. |
| `src/cmake/ofk_overlay.cmake` | The overlay, injected into OpenFAST's unmodified CMake with `-DCMAKE_PROJECT_INCLUDE=<repo>/src/cmake/ofk_overlay.cmake`. |
| `src/cmake/ofk_options.cmake` | `OFK_FAST_FMA`, `OFK_FAST_HOST_FLAGS`, `OFK_FAST_VECTOR_LIBM`, defined once (overlay and `ofk_kernels` both include it). |
| `src/cmake/ofk_sources.cmake` | For test trees: `ofk_port_sources()` (the port's copy of a file = `src/<path>` if overlaid, else the pin's) and `ofk_add_kernels()`. |
| `tools/openfast-build.sh` | One command per configuration (below). |
| `tests/<unit>/` | Unit cards and L0 tests. The oracle (`_ref`) is still built from `deps/openfast`; the port from `src/` where overlaid, else from `deps/openfast`. |
| `tests/overlay/` | The evidence scripts of this note: `a5.sh`, `cpu.sbatch` (`l0`, `sys`), `gpu_mig.sbatch`. |
| `fork/` | The old layout. Git-ignored and unregistered as a submodule; kept on disk only for queued jobs (§6). |

## 2. What the overlay does

`CMAKE_PROJECT_INCLUDE` runs after every `project()` call; the overlay acts once, for OpenFAST's top level, and schedules `_ofk_overlay_apply` with `cmake_language(DEFER)` at the end of OpenFAST's top-level `CMakeLists.txt`. At that point every target exists and OpenFAST's global settings are final (`add_definitions(-DOPENFAST_DOUBLE_PRECISION)`, the OpenMP flags and `link_libraries`, LAPACK link options), so `ofk_kernels` sees what it saw when the fork added it from its top-level `CMakeLists.txt`.

1. **Source overlay, always (also with `OPENFAST_KOKKOS=OFF`).** It globs `src/modules/*/src/*` (except `modules/kokkos`), walks every target of the build, and replaces each `SOURCES` entry that resolves to `deps/openfast/<same path>` with `src/<path>`. Source-file properties OpenFAST set on the upstream file are copied to the copy (today: `COMPILE_FLAGS -fno-trapping-math` on `FVW_BiotSavart.f90`). It fails the configure if an overlay file has no upstream counterpart or replaces nothing. Configure prints one `-- ofk overlay:` line per replacement.
2. **`OPENFAST_KOKKOS=ON`** (option defined by the overlay, default OFF): `find_package(Kokkos 5.3...<6)`; `include(src/modules/kokkos/ofk_kernels.cmake)` defines `ofk_kernels` in OpenFAST's top-level directory (a deferred call may not `add_subdirectory()`: "Subdirectories may not be created during deferred execution"); each target that received an overlay file (today only `aerodynlib`) gets `OPENFAST_KOKKOS` (PRIVATE) and links `ofk_kernels` (plain signature, as OpenFAST links `aerodynlib`).
3. **Export set.** `ofk_kernels.cmake` installs `ofk_kernels` into `EXPORT OpenFASTLibraries` when `CMAKE_PROJECT_NAME` is OpenFAST, with a `$<BUILD_INTERFACE:…>` include directory (unit 1's requirement: the exported `aerodynlib` links it).
4. CUDA: unchanged, `kokkos_compilation(GLOBAL)` with g++ as `CMAKE_CXX_COMPILER`.

The `aerodyn_utest_trapping_flags` test runs a separate `cmake` on `cmake/tests/aerodyn_trapping`, which reads the pin's `modules/aerodyn/CMakeLists.txt`; it does not see the overlay (and needs not: it checks the pin's per-file flags).

## 3. How to build

```bash
cd <repo>; OPT=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt
tools/openfast-build.sh off       $OPT/build/ov-off         # overlay, OPENFAST_KOKKOS=OFF (= the pin)
tools/openfast-build.sh omp       $OPT/build/ov-omp         # Kokkos-OpenMP, Fortran OPENMP=OFF (old fork-omp)
tools/openfast-build.sh fomp      $OPT/build/ov-fomp        # Kokkos-OpenMP + Fortran OPENMP=ON (P2 config)
tools/openfast-build.sh cuda-a100 $OPT/build/ov-cuda-a100   # sm_80, interactive_gpu MIG
tools/openfast-build.sh cuda-h100 $OPT/build/ov-cuda-h100   # sm_90, h100/h200
tools/openfast-build.sh pin       $OPT/build/ov-pin         # the plain pin, no overlay: the baseline
```
Every config: Release, `DOUBLE_PRECISION=ON`, `BUILD_TESTING=ON`, `BLA_VENDOR=OpenBLAS`, gfortran/g++ 12.3, modules `gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28`. Extra `-D` options pass through; `OFK_BUILD=0` configures only, `OFK_TARGETS="aerodyn_driver"` limits the build, `OFK_JOBS` (default 12). Logs: `<builddir>.cfg.log`, `<builddir>.build.log`. The raw command is:
```bash
cmake -S deps/openfast -B <build> -DCMAKE_PROJECT_INCLUDE=$PWD/src/cmake/ofk_overlay.cmake \
      -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DBUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CXX_COMPILER=g++ -DOPENMP=OFF \
      -DOPENFAST_KOKKOS=ON -DKokkos_ROOT=$OPT/kokkos-5.3.0-omp
```
Unit test trees (the `OFK_FORK` variable is gone; nothing else changed in their command lines):
```bash
cmake -S tests/biotsavart-direct  -B $OPT/build/ov-tests-bsd-omp -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_Fortran_COMPILER=gfortran -DOPENFAST_PIN_BUILD=$OPT/build/openfast-pin \
      -DKokkos_ROOT=$OPT/kokkos-5.3.0-omp -DCMAKE_CXX_COMPILER=g++
cmake -S tests/induced-velocities -B $OPT/build/ov-tests-iv-omp  (same options)
# CUDA (A100): -DKokkos_ROOT=$OPT/kokkos-5.3.0-cuda-a100 -DCMAKE_CXX_COMPILER=$OPT/kokkos-5.3.0-cuda-a100/bin/nvcc_wrapper -DOFK_L0_RTOL=1e-14
```

## 4. Evidence: no change in behaviour

The sources are unchanged: `src/` is a byte copy of the fork's working tree at the end of unit 3's checkpoint (21:27; `cmp` of the three `.f90` files and `diff -r` of `modules/kokkos/src` against `fork/` are empty). The only edits are in build files: `ofk_kernels`' CMake was split into `CMakeLists.txt` + `ofk_kernels.cmake` (paths relative to the file instead of `CMAKE_CURRENT_SOURCE_DIR`), and its three options moved to `src/cmake/ofk_options.cmake`. Note that the fork's shared files already carried unit 3's code when they were copied (`FVW_Subs.f90` +80 lines, `induced_velocities.{hpp,cpp}`, the treecode TUs), so "the current code" below includes the ported treecode (`OFK_TC=kokkos` by default).

Logs are under `opt/build/`; build trees `ov-*` (OpenFAST) and `ov-tests-*` (unit test trees).

### 4.0 Same sources, same machine code (overlay vs fork)
`ov-fomp` (overlay, Kokkos-OpenMP + Fortran OPENMP=ON) against `fork3-fomp` (unit 3's fork build, same configuration): the full disassembly (`objdump -d -r`, every code section including COMDAT template instantiations, with relocations) of **all 27 objects of `aerodynlib` (19) and `ofk_kernels` (8) is identical**. Compile flags of `ofk_kernels` (`flags.make` and per-file options) are identical to `fork-omp`'s; `aerodynlib`'s too.
- Caveat found on the way: `fork3-fomp`'s `induced_velocities.cpp.o`, `treecode_strict.cpp.o` and `treecode_fast.cpp.o` are **stale** (compiled 21:09, before unit 3's last edits of `treecode_impl.hpp` 21:14 and `treecode_common.hpp` 21:20). For the comparison those three were recompiled from `fork/` with `fork3-fomp`'s exact commands (scratch output; the tree itself was not touched). Unit 3 should not trust `fork3-fomp` binaries.

### 4.1 A5: overlay with OPENFAST_KOKKOS=OFF = the plain pin
`tests/overlay/a5.sh opt/build/ov-pin opt/build/ov-off` (log `ov-a5.log`; `ov-pin` = `tools/openfast-build.sh pin`, same flags, no overlay):
- `.mod`: **76/76 identical md5** (every module both trees built, including `fvw`, `fvw_subs`, `fvw_biotsavart`, `fvw_vortextools`, `fvw_types`, `aerodyn`).
- `aerodynlib` objects, full disassembly: **17/19 identical**; `FVW_BiotSavart.f90.o` differs in 3 of 4032 lines, `FVW_Subs.f90.o` in 71 of 21740. Every one of the 74 differing instructions is `movl $imm, off(%rsp)` (the line number passed to the Fortran runtime for `print`/`STOP`/error messages); the immediates move by +38/+77/+78 (`FVW_BiotSavart`, after each shim) and +86/+295 (`FVW_Subs`, after the module-level interfaces and the shims). Instruction counts are equal. Same result as unit 1's fork-off check.
- Compile flags: `flags.make` identical; per-file flags identical for all 19 sources (`-fno-trapping-math` on `FVW_BiotSavart` (carried over by the overlay) and `FVW_VortexTools`).

### 4.2 L0 of units 1 and 2 on overlay-built trees: the same numbers
**CPU** (job 1862969/1863010, `med`, Xeon 8480+, 8 threads; `tests/overlay/cpu.sbatch l0`; logs `ov-cpu.ofk-infra-l0.<job>.log`, `ov-l0-<job>/`):
- Unit 1 (`ov-tests-bsd-omp`): the script repeats the commands of unit 1's job 1859729 (seeds 1 and 7: claim 1 rtol 0, claim 3 rtol 1e-12, Fortran-FMA envelope, Fortran-fopenmp). **All 84 per-branch and summary lines are identical** to job 1859729's (`diff` empty). So: claim 1 bitwise (0 bit differences, every branch, both seeds); `ui_part_nograd` fast 1.05e-14; `ui_seg` fast 2.99e-12 (seed 1) and 1.105e-11 (seed 7), exactly as recorded.
- Unit 2 (`ov-tests-iv-omp`): 15/15 ctest pass; the 117 summary lines (`l0_iv:`, Basic/Part/TreePart/pack/A7, P6) are **identical** to the fork-built `u3-iv-omp` run of the same code (21:21). Values: claim 1 bitwise; claim 3 Basic 8.176e-14, Part 3.095e-15 (unit card: 8.2e-14 / 3.1e-15), TreePart 1.855e-15; P6 0 allocations.

**GPU** (job 1863100, `interactive_gpu` A100 MIG 2g.20gb, `tests/overlay/gpu_mig.sbatch`; log `ov-gpu.ofk-infra-gpu.1863100.log`):
- Unit 1 (`ov-tests-bsd-cuda-a100`): every per-branch line of job 1862373 (fork-built `bsd-cuda-a100`) appears unchanged; the only lines that differ are the two old `283 passed, 7 failed` summaries of claim 3, which now pass under the FMA bound (§5). Claim 2: LambOseen 3.148e-16, Exp 8.009e-16, all other branches bitwise.
- Unit 2 (`ov-tests-iv-cuda-a100`): 12/12 pass; the summary lines are **identical** to the fork-built `u3-iv-cuda-a100` run. Claim 2: Basic 1.626e-16, Part 4.35e-17 (card: 1.63e-16 / 4.35e-17); claim 3 8.245e-14 (card 8.25e-14).

### 4.3 System level: ctest `aerodyn_utest|OLAF`
(job pending, filled below)

### 4.4 BAR full-precision cases
(job pending, filled below)

## 5. Claim 3 with FMA: the per-branch bound in unit 1's ctest

User decision (2026-11-14): fast mode keeps FMA, and chain claim 3 is declared as **max(1e-12, 1.1 × the pin Fortran's own FMA-vs-no-FMA spread on the same inputs and seed), per regularization branch**.

Implementation (`tests/biotsavart-direct/`):
- `l0_util.f90`: `OFK_L0_ENVELOPE_OUT=<file>` writes the worst gated error per branch (with seed and case count) after a run; `OFK_L0_ENVELOPE=<file>` reads it and gates branch `r` at `max(OFK_L0_RTOL, OFK_L0_ENVELOPE_FACTOR × spread(r))` (factor default 1.1). A seed or case-count mismatch aborts (`error stop 2`). The summary prints each branch's bound and worst/bound.
- `CMakeLists.txt`: `l0_<r>_ffma` (pin Fortran built with `-march=x86-64-v3`, i.e. with FMA, vs strict = the Fortran without FMA, claim 1) is a ctest **fixture setup** (`FIXTURES_SETUP fma_<r>`) that writes `envelope/l0_<r>.ffma.txt`; `l0_<r>_strictfast` and `l0_<r>_fast` require it and read it. `OFK_L0_FMA_FACTOR` (cache, 1.1) and `OFK_L0_FAST_RTOL` (1e-12) are the two knobs. `OFK_L0_SEED=7 ctest` works: both runs see the same seed.
- Unit 2 needs no change: its claim 3 (8.2e-14) passes 1e-12 with a 12× margin, so the bound reduces to 1e-12.

Results (bound and worst/bound per branch; identical on CPU and on the A100):

| | None | Rankine | LambOseen | Vatistas | Offset |
|---|---|---|---|---|---|
| `ui_seg`, seed 1: bound | 3.29e-12 | 2.57e-12 | 1e-12 | 1e-12 | 1e-12 |
| worst/bound | 0.909 | 0.909 | 0.096 | 0.056 | 0.273 |
| `ui_seg`, seed 7: bound | 1.22e-11 | 1e-12 | 2.91e-12 | 1e-12 | 1e-12 |
| worst/bound | 0.909 | 0.170 | 0.909 | 0.080 | 0.412 |

`ui_part_nograd`: bound 1e-12 in every branch, worst/bound ≤ 0.010. A ratio of exactly 0.909 = 1/1.1 means the fast port moves by exactly as much as the Fortran does under FMA.

**Unit 1 ctest: 23/23 pass** (seed 1) and 16/16 of the claim-1/claim-3 subset at seed 7, on CPU (`ov-tests-bsd-omp`, job 1863010) and 23/23 at both seeds on the A100 MIG (`ov-tests-bsd-cuda-a100`, job 1863100). Before: 21/23 (`l0_ui_seg_strictfast`, `l0_ui_seg_fast` failed at 1e-12).

## 6. What remains in `fork/`


(filled below)
