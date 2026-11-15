// Unit biotsavart-direct: THE interaction bodies (architecture.md §2). One target-source
// contribution per call, templated on the regularization branch. Every loop shape and every
// backend calls these, in both modes; nothing else in the unit computes physics.
//
// Translation of FVW_BiotSavart.f90 (pin 12fa3fc):
//   seg_term  <- ui_seg inner loop body (select case at :168): None :176-196, Rankine
//                :215-236, LambOseen :255-277, Vatistas :293-313, Offset :331-347
//   part_term <- ui_part_nograd inner loop body (select case at :379): None :387-396,
//                Exp :407-420, Compact :432-446
// The guards are kept literally, including the Vatistas `<= PRECISION_UI -> cycle` style
// against the nested `> PRECISION_UI` style of the other segment branches: they differ for
// NaN (NaN passes the `<=`-cycle test and propagates; the nested `>` test skips it).
//
// FMA: this header is compiled twice, once per mode TU, with contraction off (strict) or on
// (fast). It must therefore live in a mode-specific namespace (OFK_BSD_MODE_NS), or the linker
// would merge the two instantiations of the same inline template (ODR/COMDAT) and one mode
// would silently run the other mode's arithmetic.
#pragma once

#include "biotsavart_direct_common.hpp"

#include "ofk_vector_libm.hpp"

#if !defined(OFK_BSD_MODE_NS) || !defined(OFK_BSD_STRICT)
#error "biotsavart_direct_body.hpp is included only from a mode TU that defines OFK_BSD_MODE_NS and OFK_BSD_STRICT"
#endif

namespace ofk_bsd {
namespace OFK_BSD_MODE_NS {

constexpr bool kStrict = OFK_BSD_STRICT;
// Math-library switch (architecture.md §2, the third named difference). In fast mode on a host
// backend built with OFK_VECTOR_LIBM, exp has a SIMD clone (glibc libmvec), so the bodies
// evaluate exp for every pair and select the Fortran's arm afterwards: the selected value is
// the same expression, only the exp implementation differs (libmvec, < 4 ulp). Strict mode
// and device backends keep the Fortran's branch around exp and the scalar library.
constexpr bool kVecLibm = !kStrict && kOnHost && (OFK_VECTOR_LIBM_ACTIVE != 0);

// Fortran MAX/MIN of two reals. gfortran >= 9 lowers them to MAX_EXPR/MIN_EXPR, which make no
// promise for NaN ("whatever is fastest"); a compare-select is the same value for every
// non-NaN input. Kokkos::fmax/fmin (IEEE fmax: NaN-ignoring) would also match on non-NaN
// inputs but cannot be vectorized by GCC without -ffinite-math-only, and the Fortran's own
// !$OMP SIMD loops rely on MAX vectorizing. NaN inputs give NaN outputs either way here
// (the mask multiplies a NaN cross product), and L0 checks that.
template <class T>
KOKKOS_FORCEINLINE_FUNCTION T fmax_f(const T a, const T b) { return (a > b) ? a : b; }
template <class T>
KOKKOS_FORCEINLINE_FUNCTION T fmin_f(const T a, const T b) { return (a < b) ? a : b; }

// One segment (xa..zb = CP - P1, CP - P2) at one control point. Returns false when the
// Fortran skips the segment (its guards); o = Kv*crossprod is meaningful only when true (the
// caller adds `ok ? o : 0`, which is bitwise the Fortran's skip, see impl.hpp).
// The guard expressions are the Fortran's, comparison by comparison; only their evaluation
// is branch-free (all quantities are computed, then the guards are combined), so the lane
// loop can be if-converted and vectorized. For the segments the Fortran keeps, every value is
// computed by the same expression in the same order as in the Fortran.
template <int Reg, class T>
KOKKOS_FORCEINLINE_FUNCTION bool seg_term(const T xa, const T ya, const T za, const T xb,
                                          const T yb, const T zb, const T gamma,
                                          const T regparam, const Consts<T>& k, T& o1, T& o2,
                                          T& o3) {
  const T norm_a = Kokkos::sqrt(xa * xa + ya * ya + za * za);
  const T norm_b = Kokkos::sqrt(xb * xb + yb * yb + zb * zb);
  T denominator = norm_a * norm_b * (norm_a * norm_b + xa * xb + ya * yb + za * zb);
  const T c1 = ya * zb - za * yb, c2 = za * xb - xa * zb, c3 = xa * yb - ya * xb;
  T norm2_orth = c1 * c1 + c2 * c2 + c3 * c3;
  const T norm2_r0 = (xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb);
  bool ok;
  if constexpr (Reg == kRegVatistas) {
    // :299, :303, :306  `if (x <= PRECISION_UI) cycle` (NaN passes and propagates)
    ok = !(denominator <= k.precision_ui) & !(norm2_orth <= k.precision_ui) &
         !(norm2_r0 <= k.precision_ui);
  } else {
    // nested `if (x > PRECISION_UI) then` (NaN is skipped), e.g. :184, :187, :189
    ok = (denominator > k.precision_ui) & (norm2_orth > k.precision_ui) &
         (norm2_r0 > k.precision_ui);
  }
  T Kv;
  if constexpr (Reg == kRegNone) {  // :192
    Kv = gamma * k.fourpi_inv * (norm_a + norm_b) / (denominator + k.mindenom);
  } else if constexpr (Reg == kRegRankine) {  // :224-232
    norm2_orth = norm2_orth / norm2_r0;
    const T r_bar2 = norm2_orth / (regparam * regparam);
    if (r_bar2 < T(1)) Kv = r_bar2; else Kv = T(1);
    Kv = gamma * k.fourpi_inv * Kv * (norm_a + norm_b) / (denominator + k.mindenom);
  } else if constexpr (Reg == kRegLambOseen) {  // :264-273
    norm2_orth = norm2_orth / norm2_r0;
    const T r_bar2 = norm2_orth / (regparam * regparam);
    const T exp_value = -k.lamb_oseen * r_bar2;
    if constexpr (kVecLibm) {
      // exp only matters when exp_value >= MIN_EXP_VALUE; clamping below that keeps libmvec off
      // its slow underflow path and leaves every used value unchanged (NaN passes through).
      const T arg = (exp_value < k.min_exp_value) ? k.min_exp_value : exp_value;
      const T e = Kokkos::exp(arg);
      Kv = (exp_value < k.min_exp_value) ? T(1) : T(1) - e;
    } else {
      if (exp_value < k.min_exp_value) Kv = T(1);
      else Kv = T(1) - Kokkos::exp(exp_value);
    }
    Kv = gamma * k.fourpi_inv * Kv * (norm_a + norm_b) / (denominator + k.mindenom);
  } else if constexpr (Reg == kRegVatistas) {  // :309-312
    norm2_orth = norm2_orth / norm2_r0;
    const T r_bar2 = norm2_orth / (regparam * regparam);
    Kv = r_bar2 / Kokkos::sqrt(T(1) + r_bar2 * r_bar2);
    Kv = gamma * k.fourpi_inv * Kv * (norm_a + norm_b) / (denominator + k.mindenom);
  } else {  // :342-343
    static_assert(Reg == kRegOffset, "unknown segment regularization");
    denominator = denominator + regparam * regparam * norm2_r0;
    Kv = gamma * k.fourpi_inv * (norm_a + norm_b) / (denominator + k.mindenom);
  }
  o1 = Kv * c1; o2 = Kv * c2; o3 = Kv * c3;
  return ok;
}

// One particle (dx..dz = CP - Part, alpha, RegParam) at one control point. Returns false when
// the Fortran `cycle`s (Exp, r2 < MINNORM2); o = C*sp is meaningful only when true.
template <int Reg, class T>
KOKKOS_FORCEINLINE_FUNCTION bool part_term(const T dx, const T dy, const T dz, const T ax,
                                           const T ay, const T az, const T regparam,
                                           const Consts<T>& k, T& o1, T& o2, T& o3) {
  const T r2 = dx * dx + dy * dy + dz * dz;
  T sp;
  if constexpr (Reg == kRegNone) {  // :386-397
    const T msk = (r2 >= k.minnorm2) ? T(1) : T(0);  // merge(1,0, r2>=MINNORM2)
    const T r2s = fmax_f(r2, k.minnorm2);
    const T rn = Kokkos::sqrt(r2s);
    const T r3 = r2s * rn;
    const T Cx = ay * dz - az * dy;
    const T Cy = az * dx - ax * dz;
    const T Cz = ax * dy - ay * dx;
    sp = msk * k.fourpi_inv / r3;
    o1 = Cx * sp; o2 = Cy * sp; o3 = Cz * sp;
    return true;
  } else if constexpr (Reg == kRegPartExp) {  // :406-421
    // :409 `if (r2 < MINNORM2) cycle` (on singularity), evaluated branch-free like seg_term
    const bool ok = !(r2 < k.minnorm2);
    const T rn = Kokkos::sqrt(r2);
    const T r3 = r2 * rn;
    const T rc3 = regparam * regparam * regparam;
    const T Cx = ay * dz - az * dy;
    const T Cy = az * dx - ax * dz;
    const T Cz = ax * dy - ay * dx;
    if constexpr (kVecLibm) {
      // exp is used only when r3 <= 8 rc3, i.e. -(r3/rc3) >= -8 (both exact in IEEE); clamp
      // there (NaN passes through) so far-field lanes do not hit libmvec's underflow path.
      const T arg0 = -(r3 / rc3);
      const T arg = (arg0 < -k.part_reg_cut3) ? -k.part_reg_cut3 : arg0;
      const T e = Kokkos::exp(arg);
      sp = (r3 > k.part_reg_cut3 * rc3) ? k.fourpi_inv / r3 : (T(1) - e) * k.fourpi_inv / r3;
    } else {
      if (r3 > k.part_reg_cut3 * rc3) sp = k.fourpi_inv / r3;  // r > 2rc: skip the exp
      else sp = (T(1) - Kokkos::exp(-(r3 / rc3))) * k.fourpi_inv / r3;
    }
    o1 = Cx * sp; o2 = Cy * sp; o3 = Cz * sp;
    return ok;
  } else {  // :431-447
    static_assert(Reg == kRegPartCompact, "unknown particle regularization");
    const T msk = (r2 >= k.minnorm2) ? T(1) : T(0);
    const T r2s = fmax_f(r2, k.minnorm2);
    const T rn = Kokkos::sqrt(r2s);
    const T r3 = r2s * rn;
    const T a = k.part_reg_c2 * regparam;
    const T rc2 = fmax_f(a * a, k.minnorm2);                      // (C2*Reg)**2, floored
    const T rc3 = rc2 * fmax_f(k.part_reg_c2 * regparam, k.minnorm);
    const T tt = fmin_f(r2s / rc2, T(1));
    const T Cx = ay * dz - az * dy;
    const T Cy = az * dx - ax * dz;
    const T Cz = ax * dy - ay * dx;
    // merge(): both arms are finite thanks to the floors above
    sp = (r2s >= rc2) ? k.fourpi_inv / r3
                      : (T(35) + tt * (T(-42) + T(15) * tt)) * k.fourpi_inv / (T(8) * rc3);
    sp = msk * sp;
    o1 = Cx * sp; o2 = Cy * sp; o3 = Cz * sp;
    return true;
  }
}

}  // namespace OFK_BSD_MODE_NS
}  // namespace ofk_bsd
