// Vector math library switch for FAST-mode host translation units (architecture.md §2: the
// math library is one of the three named differences between runs).
//
// When the TU is compiled with -DOFK_VECTOR_LIBM=1 (CMake option OFK_FAST_VECTOR_LIBM, fast TUs
// only), on x86-64 glibc with GCC, exp/expf are redeclared with a SIMD clone ("notinbranch"),
// exactly as glibc's own <bits/math-vector.h> does under -ffast-math. GCC then vectorizes loops
// that call exp with libmvec's _ZGV*_exp (glibc >= 2.22, linked through libm.so's AS_NEEDED
// libmvec). Results differ from the scalar libm by a few ulp and are deterministic.
// Never active in device code or in strict TUs.
#pragma once

#include <cmath>

#if defined(OFK_VECTOR_LIBM) && OFK_VECTOR_LIBM && defined(__GNUC__) && !defined(__clang__) && \
    !defined(__CUDACC__) && defined(__x86_64__) && defined(__GLIBC__) && !OFK_BSD_STRICT
#define OFK_VECTOR_LIBM_ACTIVE 1
extern "C" double exp(double) noexcept __attribute__((simd("notinbranch")));
extern "C" float expf(float) noexcept __attribute__((simd("notinbranch")));
#else
#define OFK_VECTOR_LIBM_ACTIVE 0
#endif
