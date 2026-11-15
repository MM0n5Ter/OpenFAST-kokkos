# OpenFAST-kokkos build options, in one place. Included by the overlay (src/cmake/ofk_overlay.cmake,
# for an OpenFAST build) and by src/modules/kokkos/CMakeLists.txt (for the unit test trees, which add
# ofk_kernels directly). option()/set(CACHE) never override a value that is already in the cache.
include_guard(GLOBAL)

# Host ISA flags for the FAST translation units only: FMA contraction needs FMA instructions,
# and the baseline x86-64 ISA has none. x86-64-v3 = AVX2 + FMA (Haswell, 2013, and later).
# Set to "" to keep the compiler's default ISA (fast mode then has no FMA on the host), or to
# e.g. -march=native for the build machine. The STRICT units never get these flags.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(_ofk_fast_default "-march=x86-64-v3")
else()
  set(_ofk_fast_default "")
endif()
set(OFK_FAST_HOST_FLAGS "${_ofk_fast_default}" CACHE STRING
    "Host ISA flags for the fast-mode translation units of ofk_kernels (FMA needs e.g. -march=x86-64-v3)")

# OFK_FAST_FMA=OFF keeps the fast mode's blocked order and vectorization but compiles it without
# contraction. ui_seg's Biot-Savart segment formula is ill-conditioned near the segment line
# (cross product cancellation), and FMA alone moves it by up to ~1e-11 of sum|c| on L0 inputs,
# exactly as much as the Fortran itself moves when built with FMA (tests/biotsavart-direct/unit.md).
# Kept ON by user decision (2026-11-14); claim 3 then uses the Fortran's own FMA spread as bound.
option(OFK_FAST_FMA "fast mode: allow FMA contraction (architecture.md §3)" ON)

# Vector libm (glibc libmvec) for exp in the fast host TUs: lets the exp-based regularizations
# vectorize, at the price of evaluating exp for every pair (the Fortran skips it in the far
# field). Changes exp by a few ulp in fast mode only (src/ofk_vector_libm.hpp). Measured
# SLOWER on Xeon 8480+ (part Exp 0.99x vs 1.12x, seg LambOseen 0.74x vs 1.31x of the Fortran,
# 1 thread, 4000x4000), so OFF by default; kept as the named switch for unit 3's near field.
option(OFK_FAST_VECTOR_LIBM "fast mode, host: vectorized exp from glibc libmvec (slower here)" OFF)
