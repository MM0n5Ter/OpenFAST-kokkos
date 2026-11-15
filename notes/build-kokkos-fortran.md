# How to build Kokkos and link it with OpenFAST's Fortran

> **Update 2026-11-14:** §2 was run with **Kokkos 5.3.0** (release tarball), C++20, gcc 12.3, cuda 12.3.0 and cmake 3.31.12.
> - The script is `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/build-kokkos-5.3.0.sh`.
> - Installs: `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/kokkos-5.3.0-{omp,cuda-h100}`. The CUDA install is `Kokkos_ARCH_HOPPER90` instead of the a100 example below. Both are compiled on the login node.
> - §3 passes on the OpenMP install (see `examples/fortran-kokkos/README.md`). The CUDA toy passes on the H100 (job 1858536: worst relative error 2.3e-18, `capi_smoke` OK).
> - An earlier pass with 4.5.1, the cluster's source tree, also passed on the H100. That install is superseded.
> - The pre-built `pack4/third_party/kokkos-install` is **static, non-PIC, and has no OpenMP**, and its CMake config points at a deleted `nvcc_wrapper`. Do not use it for OpenFAST: the C-binding `.so` needs PIC/shared.
>
> **Status (2026-11-12): the commands in this note have not been run yet.**
> - Checked: the Fortran side of `examples/fortran-kokkos` (shim, renamed oracle, L0 test, shared library loaded from Python) was compiled and run against a plain C++ stand-in for the kernel.
> - Checked by reading the Kokkos 5.3.0 source, not by building: the API names (`create_mirror_view_and_copy`, `is_finalized`, `push_finalize_hook`, `numbers::pi`), the CMake options, and the architecture names.
> - Week 1, Lane A runs §2–§4 and fixes this note (milestone M1).

## 1. Toolchain on our cluster

| Need | Use | Why |
|---|---|---|
| Fortran | `gfortran` 12.3 (`module load gcc/12.3.0`) | OpenFAST accepts only GNU, Intel or Flang (`cmake/OpenfastFortranOptions.cmake:48-52`). nvfortran is rejected. |
| C++ (host) | `g++` 12.3 | Kokkos 5 needs GCC ≥ 10.4 and C++20 |
| CUDA | `cuda/12.6.3`; `cuda/12.8.1` or newer for Blackwell | Kokkos 5 needs NVCC ≥ 12.2. CUDA 12.6 supports gcc 12 as host compiler. |
| CMake | `cmake/3.31.12` | Kokkos 5 needs ≥ 3.22. OpenFAST needs ≥ 3.15. |
| Python | `python/3.11.9` with numpy | r-test scripts, our oracle renamer, the C-binding tests |
| BLAS/LAPACK | a module if there is one, else `-DUSE_LOCAL_STATIC_LAPACK=ON` | OpenFAST needs LAPACK. That option downloads and builds a reference copy. |

Kokkos GPU architecture per Slurm partition (one Kokkos install per architecture):

| Partition | GPU | Kokkos option |
|---|---|---|
| `a100` | A100 | `-DKokkos_ARCH_AMPERE80=ON` (development) |
| `h100`, `h200` | H100/H200 | `-DKokkos_ARCH_HOPPER90=ON` (final results) |
| `l40s` | L40S | `-DKokkos_ARCH_ADA89=ON` |
| `b200`, `b300` | B200/B300 | `-DKokkos_ARCH_BLACKWELL100=ON` / `BLACKWELL103` (CUDA ≥ 12.8) |
| `interactive_gpu` | MIG slices (2g.20gb, A100-class *(verify)*) | `AMPERE80` |

Compiling works on a login node, which has no GPU. CUDA executables must run on a GPU node (`srun -p a100 --gres=gpu:1 …`; add `-A <account>` if needed).

## 2. Build Kokkos from source (two installs)

```bash
module load gcc/12.3.0 cmake/3.31.12 python/3.11.9
export OFK_PREFIX=/path/shared/by/the/team/opt          # pick once and write it in decisions.md
git clone --depth 1 -b 5.3.0 https://github.com/kokkos/kokkos.git kokkos-5.3.0

# (a) Host only: Serial + OpenMP. Used for the Kokkos-OpenMP build of OpenFAST.
cmake -S kokkos-5.3.0 -B build-kokkos-omp \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=$OFK_PREFIX/kokkos-5.3.0-omp \
  -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON \
  -DBUILD_SHARED_LIBS=ON
cmake --build build-kokkos-omp -j 16 && cmake --install build-kokkos-omp

# (b) CUDA (device) + OpenMP (host), one per GPU architecture.
module load cuda/12.6.3
export NVCC_WRAPPER_DEFAULT_COMPILER=g++                 # host compiler behind nvcc
cmake -S kokkos-5.3.0 -B build-kokkos-a100 \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=$PWD/kokkos-5.3.0/bin/nvcc_wrapper -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=$OFK_PREFIX/kokkos-5.3.0-cuda-a100 \
  -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON \
  -DBUILD_SHARED_LIBS=ON
cmake --build build-kokkos-a100 -j 16 && cmake --install build-kokkos-a100
```

Why these options:
- **`BUILD_SHARED_LIBS=ON`:** one `libkokkoscore.so` per process, even when the port sits in several libraries or inside a `dlopen`ed C binding. Two static copies of Kokkos in one process would mean two runtimes.
- **No relocatable device code** (the default). It does not mix with shared libraries, and we don't need it: each kernel lives in one `.cpp`.
- **No `Kokkos_ARCH_NATIVE`:** login and compute nodes may have different CPUs.
- **Debug variant** (optional, for L0 failures): add `-DKokkos_ENABLE_DEBUG=ON -DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON` and install to `…-debug`.
- The CUDA install's default execution space is `Cuda`, and its host space is `OpenMP`. The OpenMP install's default space is `OpenMP`. The same source runs on both.

## 3. Check the install with the toy (`examples/fortran-kokkos`)

```bash
# OpenMP build: runs on the login node
cmake -S examples/fortran-kokkos -B build-toy-omp \
  -DKokkos_ROOT=$OFK_PREFIX/kokkos-5.3.0-omp -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build-toy-omp -j
OMP_NUM_THREADS=8 ctest --test-dir build-toy-omp --output-on-failure      # l0_ui_part, capi_smoke

# CUDA build: compile here, run on a GPU node
cmake -S examples/fortran-kokkos -B build-toy-a100 \
  -DKokkos_ROOT=$OFK_PREFIX/kokkos-5.3.0-cuda-a100 \
  -DCMAKE_CXX_COMPILER=$OFK_PREFIX/kokkos-5.3.0-cuda-a100/bin/nvcc_wrapper -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build-toy-a100 -j
srun -p a100 --gres=gpu:1 -t 00:10:00 env OFK_L0_RTOL=1e-10 ctest --test-dir build-toy-a100 --output-on-failure
```

Expected results:
- `l0_ui_part: 56 passed, 0 failed`.
  - On OpenMP the worst relative error should be ≤ 1e-14: per-point sums run in the same order as the Fortran.
  - On CUDA it may be larger, because of FMA contraction and device `exp`. That is why the CUDA run sets `OFK_L0_RTOL=1e-10`.
- `capi_smoke: OK`. The shared library loads from Python, Kokkos starts inside it, and the process exits cleanly.

## 4. Build OpenFAST: the baseline and the Kokkos builds

**Baseline** (unmodified pin, the oracle):
```bash
cmake -S deps/openfast -B build-of-base \
  -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DOPENMP=ON \
  -DBUILD_TESTING=ON -DBUILD_UNIT_TESTING=ON -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
cmake --build build-of-base -j 16
git -C deps/openfast submodule update --init --depth 1 reg_tests/r-test
cd build-of-base && OMP_NUM_THREADS=1 ctest -R 'OLAF|5MW_MRSemi_DLL_WSt_WavesIrr|aerodyn_utest' -j 8 --output-on-failure
```

**Kokkos builds** (in the port fork; same flags plus):

| Build | Extra CMake flags |
|---|---|
| Kokkos-OpenMP | `-DOPENFAST_KOKKOS=ON -DKokkos_ROOT=$OFK_PREFIX/kokkos-5.3.0-omp -DCMAKE_CXX_COMPILER=g++` |
| Kokkos-CUDA | `-DOPENFAST_KOKKOS=ON -DKokkos_ROOT=$OFK_PREFIX/kokkos-5.3.0-cuda-a100 -DCMAKE_CXX_COMPILER=$OFK_PREFIX/kokkos-5.3.0-cuda-a100/bin/nvcc_wrapper` |

Use one build directory per variant: `build-of-base`, `build-of-kk-omp`, `build-of-kk-a100`.

## 5. The CMake changes in the OpenFAST fork

Three small edits. All of them are owned by Lane A.

**(1) Top-level `CMakeLists.txt`,** next to the other `option(...)` lines and before the modules are added:
```cmake
option(OPENFAST_KOKKOS "Build ported units in C++/Kokkos (OpenFAST-kokkos)" off)
if (OPENFAST_KOKKOS)
  find_package(Kokkos 5 REQUIRED)                      # -DKokkos_ROOT=<install>
  add_subdirectory(modules/kokkos)                     # defines ofk_kernels
endif()
```

**(2) New `modules/kokkos/CMakeLists.txt`:** one library for **all** ported units, so the process has one Kokkos runtime and one `ofk_runtime`:
```cmake
add_library(ofk_kernels STATIC
  src/ofk_runtime.cpp
  src/fvw_induction.cpp          # OLAF unit: direct kernels + treecode
)
target_include_directories(ofk_kernels PUBLIC src)
target_link_libraries(ofk_kernels PUBLIC Kokkos::kokkos)
set_property(TARGET ofk_kernels PROPERTY POSITION_INDEPENDENT_CODE ON)   # ends up in .so files too
```

**(3) The module that owns the unit** (`modules/aerodyn/CMakeLists.txt`):
```cmake
if (OPENFAST_KOKKOS)
  target_link_libraries(aerodynlib ofk_kernels)        # plain signature: OpenFAST uses it for aerodynlib,
                                                       # and CMake refuses to mix plain and keyword forms
  target_compile_definitions(aerodynlib PRIVATE OPENFAST_KOKKOS)
endif()
```

How the pieces link, and why this works:
- **Linker language.** CMake links an executable or shared library with the language of highest linker preference among its sources *and its static dependencies*. C++ beats Fortran. So every target that links `aerodynlib` is linked by the C++ compiler, and CMake adds the gfortran runtime itself. This covers `aerodyn_driver`, `openfast`, `aerodyn_utest` and `libaerodyn_inflow_c_binding.so`. If a target still ends up linked by gfortran (undefined `std::`/`Kokkos::` symbols), set `set_target_properties(<t> PROPERTIES LINKER_LANGUAGE CXX)`.
- **Kokkos flags never reach Fortran.** Kokkos guards its compile options with `$<COMPILE_LANGUAGE:CXX>` and its link options with `$<LINK_LANGUAGE:CXX>` (`cmake/kokkos_tribits.cmake:289-316` in Kokkos 5.3.0).
- **Precision.** `add_definitions(-DOPENFAST_DOUBLE_PRECISION)` in `OpenfastFortranOptions.cmake` applies to C++ as well. So the C++ side picks `using ReKi = double` or `float` from the same macro.
- **The Fortran files keep both bodies during the month.** OpenFAST compiles `.f90` with `-cpp`, so a ported routine can be:
  ```fortran
  subroutine ui_seg(...)            ! interface unchanged
  #ifdef OPENFAST_KOKKOS
     call ofk_fvw_ui_seg(...)       ! bind(C) shim
  #else
     ... original Fortran body ...
  #endif
  ```
  With the option off, the fork builds exactly like upstream, which also makes it upstreamable. The `#else` branch is deleted when the repo goes Fortran-free.

**Compiler choice for the CUDA build** (open decision; try in this order):
- **A, the default.** `CMAKE_CXX_COMPILER=nvcc_wrapper` for the whole OpenFAST configure. This is the standard Kokkos route. All of OpenFAST's own C++ goes through nvcc_wrapper too: the registry tool, and openfast-cpp if enabled. That is slower to compile but should work.
- **B.** Keep `g++` as the C++ compiler, and call `kokkos_compilation(TARGET ofk_kernels)` after `find_package(Kokkos)`. Kokkos's launcher then swaps in nvcc_wrapper for that target only. Check the final link: Kokkos's CUDA link options apply when the link language is C++.
- **C, a full isolation fallback.** Build `ofk_kernels` as its own CMake project into `libofk.so`, which exports only `extern "C"` functions. OpenFAST links it as an imported library and never sees Kokkos.

### Verified in the fork (unit 1, 2026-11-14)
- `option(OPENFAST_KOKKOS OFF)`, `find_package(Kokkos 5.3...<6)` and `add_subdirectory(modules/kokkos)` in the top-level `CMakeLists.txt`.
- `aerodynlib` gets `OPENFAST_KOKKOS` (PRIVATE) and links `ofk_kernels`.
- `ofk_kernels` must be in `install(EXPORT OpenFASTLibraries)`, with a `$<BUILD_INTERFACE:…>` include directory.
- **CUDA:** `kokkos_compilation(GLOBAL)`, keeping g++ as `CMAKE_CXX_COMPILER`. Option A, a global `nvcc_wrapper`, fails on `awae/src/vtk.cpp`.

## 6. Rules at the Fortran ↔ C++ boundary

Every unit follows these. The toy shows each one, and L0 tests check them.

| Topic | Rule |
|---|---|
| Interface | The Fortran module keeps its **name and public interface**. Inside, an `interface … bind(C, name='ofk_<unit>_<routine>')` declares the C++ function. C++ side: `extern "C"`. |
| Scalars | `integer(c_int), value` (IntKi = 4-byte = `c_int`). Logical: pass `integer(c_int)` 0/1, not `logical` (`c_bool` only if both sides use it). |
| Reals | `real(c_double)`. Passing `real(ReKi)` to it **fails to compile** when ReKi ≠ c_double, which is how a single-precision build is caught. For single precision, add a `c_float` interface. |
| Arrays | Explicit-shape dummies `(3,*)`. A contiguous actual is passed by address. A non-contiguous one is copied in and out by the compiler (correct, but it costs time; `-fcheck=array-temps` shows where). |
| Layout | Fortran `A(3,n)` is column-major → `Kokkos::View<double**, Kokkos::LayoutLeft>` with extents `(3,n)`, indexed `A(k,i)` 0-based. Lower bounds other than 1 (`0:N`) are passed as explicit offsets. |
| Index arrays | Values inside index arrays are **Fortran 1-based**: `SegConnct` holds point indices. Subtract 1 in C++, and say so in the audit list. |
| Extents | Many OLAF arrays are "3 × n++": larger than used. Touch only `[0,n)`. The L0 test checks the tail with sentinels. |
| Side effects | Keep `intent(inout)` meaning. OLAF **accumulates** into `Uind`. Copy it to the device and add; never overwrite. |
| Derived types | Registry types and pointer types (`T_Tree`, `T_Node`) do not cross the boundary. Flatten them into arrays. Long-lived C++ state lives behind an opaque `type(c_ptr)` handle stored in a **unit-internal** type, never in a registry `*_Types.f90`. |
| Instances | No per-instance state in C++ globals: FAST.Farm, multi-rotor and the Python bindings can create several AeroDyn/OLAF instances in one process. One handle per instance. |
| Errors | Routines with `ErrStat`/`ErrMsg` get an `int*` status and a `char*` buffer of `ErrMsgLen` (8196, `NWTC_Base.f90:36`), and the Fortran shim copies the message. Catch every C++ exception in the `extern "C"` function; nothing is thrown into Fortran. |
| Kokkos lifetime | `ofk::ensure_initialized()` on first use; `Kokkos::finalize` from `atexit`, registered after initialization so it runs before the CUDA runtime tears down. Persistent Views are released through `Kokkos::push_finalize_hook`. Never initialize inside an OpenMP parallel region. If the host program already initialized Kokkos, leave it alone. |
| Device memory | Allocate device buffers in `*_Init` (or on first use) and keep them. Grow them by doubling as the wake grows. Allocating per call is fine for the toy and too slow for OLAF. |
| Threads | Fortran OpenMP and Kokkos-OpenMP share libgomp and `OMP_NUM_THREADS`. Set `OMP_PROC_BIND=spread OMP_PLACES=threads` (otherwise Kokkos warns). Don't call a Kokkos kernel from inside a Fortran `!$OMP PARALLEL` region. |
| Reproducibility | Keep the Fortran summation order (per-point accumulators, then one add). For bitwise L0 on CPU, compile both sides without FMA contraction (`-ffp-contract=off`). On CUDA, compare with a tolerance (`--fmad=false` only for debugging). |

## 7. Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `undefined reference to std::…` / `Kokkos::…` at link | Target linked by gfortran. Set `LINKER_LANGUAGE CXX` on it, or make sure `ofk_kernels` is a static dependency. |
| `undefined reference to _gfortran_…` | Linked by C++ without the Fortran runtime. Let CMake pick the linker (§5), or add `gfortran` to the link. |
| `relocation … can not be used when making a shared object; recompile with -fPIC` | `ofk_kernels` or Kokkos built without PIC. Set `POSITION_INDEPENDENT_CODE ON` and use the shared Kokkos install. |
| `Kokkos::initialize() has already been called` / two Kokkos runtimes | Kokkos linked statically into two libraries. Use the shared Kokkos install and one `ofk_kernels`. |
| `cudaErrorNoKernelImageForDevice` / "no kernel image is available" | Kokkos built for another GPU architecture. Use the install matching the partition (§1). |
| CUDA errors at program exit | Finalize ordering. Check that `ensure_initialized` registered the atexit handler after `Kokkos::initialize`. Free persistent Views in a finalize hook. |
| `OMP_PROC_BIND environment variable not set` | A warning only. Set `OMP_PROC_BIND=spread OMP_PLACES=threads`. |
| `Type mismatch in argument` at the `bind(C)` call | ReKi is not `c_double` (single-precision build), or an integer kind differs. Fix the interface, not the caller. |
| `cmake` refuses to mix keyword and plain `target_link_libraries` | Use the plain signature for OpenFAST's own targets (§5). |
| Fortran `.mod` clash between oracle and port | Each library needs its own `Fortran_MODULE_DIRECTORY` (see the toy's CMake). |
