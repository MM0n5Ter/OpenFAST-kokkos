# Fortran → C++/Kokkos: the per-unit template

A toy with the same shape as an OLAF unit: a cut-down `ui_part_nograd` (vortex-particle Biot–Savart, exponential core). Every ported OpenFAST unit follows this pattern. The full explanation is in [`notes/build-kokkos-fortran.md`](../../notes/build-kokkos-fortran.md) and [`notes/harness.md`](../../notes/harness.md).

| File | Role | In OpenFAST |
|---|---|---|
| `ref/toy_kinds.f90` | Shared kinds, never renamed | `NWTC_Library` |
| `ref/toy_biotsavart.f90` | **Frozen** Fortran: the oracle | a file in `deps/openfast` |
| `src/toy_biotsavart.f90` | Port: **same module name and interface**, body forwards through `bind(C)` | the fork's version of the file |
| `src/toy_biotsavart_kokkos.cpp` | `extern "C"` Kokkos kernel | `modules/kokkos/src/<unit>.cpp` |
| `src/ofk_runtime.{hpp,cpp}` | Lazy `Kokkos::initialize`, `finalize` at exit | shared by all units |
| `src/toy_c_binding.f90` | `bind(C)` entry in a **shared** library | `AeroDyn_Inflow_C_Binding.f90` |
| `rename_modules.py` | Renames the frozen file's modules to `*_ref` for the oracle | reused as is |
| `tests/l0_ui_part.f90` | L0: ref vs port in one process, random inputs + edge cases (singularities, sentinel tails, accumulation, non-contiguous sections) | `tests/<unit>/l0_*.f90` |
| `tests/capi_smoke.py` | API gate: load the `.so` from Python with ctypes, the way pyOpenFAST does | `py_*_regression` cases |

Build and run (needs a Kokkos 5 install, see the build note §2):
```bash
cmake -S . -B build -DKokkos_ROOT=<kokkos install> [-DCMAKE_CXX_COMPILER=<install>/bin/nvcc_wrapper]
cmake --build build -j
ctest --test-dir build --output-on-failure          # on CUDA: OFK_L0_RTOL=1e-10
```

**Checked so far (2026-11-12):**
- The Fortran side was compiled with gfortran 12.3 (`-fcheck=all`) and linked against a plain serial C++ stand-in for `ofk_toy_ui_part`.
- The stand-in was not Kokkos, so `toy_biotsavart_kokkos.cpp`, `ofk_runtime.cpp` and `CMakeLists.txt` have **not** been built yet.
- Results:
  - L0: 56/56 and 206/206 cases (seeds 1 and 7) at relative error 0;
  - `capi_smoke`: OK;
  - the oracle exports only `__toy_biotsavart_ref_MOD_*`, so it coexists with the port's `__toy_biotsavart_MOD_*`.
