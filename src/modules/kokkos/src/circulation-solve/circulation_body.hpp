// Unit circulation-solve: THE body (architecture.md §2), split at the Gamma statement of
// ui_seg_11 (FVW_BiotSavart.f90:125 on the pin):
//   seg_geom  <- ui_seg_11 :67-123, everything that does not depend on SegGamma, per (segment, target)
//   replay    <- ui_seg_11 :125-126 and ui_quad_n1's accumulation :532-546, per target, in the pin's
//                order: rings in (iW, iSpan, iDepth) order, then segments 1-2, 3-4, 2-3, 4-1.
//                The cross product (:81) is recomputed here by the same expressions from the same
//                inputs (same bits) instead of being stored: the replay is bound by memory traffic,
//                and 3 cached values per pair instead of 6 halve it (P2).
// Both run in the Kokkos HOST execution space in every build (CUDA builds too): the work scales
// with nSpan^2, not with the wake (criteria L1 fails), and the iteration around it calls
// AirfoilInfo (L4). See tests/circulation-solve/unit.md.
//
// Every value is computed by the Fortran's expression in the Fortran's operation order; only
// the guards are evaluated branch-free (all lanes computed, early-return lanes marked), so the
// target loop vectorizes. The replay adds `masked ? 0 : Kv*crossprod`, which is the pin's
// `Uind = Uind + Uindtmp` with Uindtmp = 0 (its early returns): Vvar starts at +0 and a
// round-to-nearest sum that starts at +0 never becomes -0, so adding +0 changes nothing.
//
// Compiled once per mode TU, inside a mode namespace (ODR, architecture.md §3):
//   strict  : baseline ISA, -ffp-contract=off                (bitwise vs the pin)
//   fast    : OFK_FAST_HOST_FLAGS (x86-64-v3), -ffp-contract=off   (F0 of the FMA experiment)
//   fastfma : OFK_FAST_HOST_FLAGS (x86-64-v3), -ffp-contract=fast  (F1 of the FMA experiment)
#pragma once

#include "circulation_common.hpp"

#include <Kokkos_Core.hpp>

#include <cmath>

#if !defined(OFK_CS_NS)
#error "circulation_body.hpp is included only from a mode TU that defines OFK_CS_NS"
#endif

#define OFK_CS_PRAGMA(x) _Pragma(#x)
#if defined(_OPENMP)
#define OFK_CS_SIMD OFK_CS_PRAGMA(omp simd)
#else
#define OFK_CS_SIMD
#endif

// Host-only force-inline (these run in the host space only; KOKKOS_FORCEINLINE_FUNCTION would make them
// __host__ __device__ under nvcc, where the constexpr tables below are not visible). Under g++ this is what
// KOKKOS_FORCEINLINE_FUNCTION expands to.
#define OFK_CS_FORCEINLINE inline __attribute__((always_inline))

namespace ofk_cs {
namespace OFK_CS_NS {

using HostExec = Kokkos::DefaultHostExecutionSpace;

// ui_seg_11's segment order inside ui_quad_n1 (FVW_BiotSavart.f90:532-546): (Pa, Pb) =
// (P1,P2), (P3,P4), (P2,P3), (P4,P1), as corner indices into the ring record.
constexpr int kSegA[4] = {0, 2, 1, 3};
constexpr int kSegB[4] = {1, 3, 2, 0};

// ui_quad_n1: DP1 = CP - Pa, DP2 = CP - Pb (:532, :536, :540, :544); crossprod (:81). Shared by the
// build and the replay, so both compute the same bits.
template <class T>
OFK_CS_FORCEINLINE void seg_cross(const T xa, const T ya, const T za, const T xb, const T yb, const T zb,
                                           T& c1, T& c2, T& c3) {
  c1 = ya * zb - za * yb;
  c2 = za * xb - xa * zb;
  c3 = xa * yb - ya * xb;
}

template <class T>
struct SegEnds {
  T pax, pay, paz, pbx, pby, pbz, reg;
};
template <class T>
OFK_CS_FORCEINLINE SegEnds<T> seg_ends(const Geom<T>& g, const int j) {
  const int r = j >> 2, sg = j & 3;
  const T* P = g.ring + 13 * r;
  return {P[3 * kSegA[sg]], P[3 * kSegA[sg] + 1], P[3 * kSegA[sg] + 2],
          P[3 * kSegB[sg]], P[3 * kSegB[sg] + 1], P[3 * kSegB[sg] + 2], P[12]};
}

// The Gamma-independent part of ui_seg_11 for source segment j on every target (lanes t).
template <int Reg, class T>
inline void seg_geom(const Geom<T>& g, const int j) {
  const SegEnds<T> e = seg_ends(g, j);
  const T pax = e.pax, pay = e.pay, paz = e.paz, pbx = e.pbx, pby = e.pby, pbz = e.pbz;
  const T reg = e.reg;
  const T prec = g.k.precision_ui;
  const std::size_t o = std::size_t(j) * std::size_t(g.nTpad);
  const T* __restrict cx = g.cx;
  const T* __restrict cy = g.cy;
  const T* __restrict cz = g.cz;
  T* __restrict kvo = g.kv + o;
  T* __restrict so = g.s + o;
  T* __restrict deno = g.den + o;
  const T lamb = g.k.lamb_oseen, minexp = g.k.min_exp_value;
  auto elem = [&](const int t) {
    // ui_quad_n1: DP1 = CP - Pa, DP2 = CP - Pb (:532, :536, :540, :544)
    const T xa = cx[t] - pax, ya = cy[t] - pay, za = cz[t] - paz;
    const T xb = cx[t] - pbx, yb = cy[t] - pby, zb = cz[t] - pbz;
    const T norm_a = std::sqrt(xa * xa + ya * ya + za * za);                        // :76
    const T norm_b = std::sqrt(xb * xb + yb * yb + zb * zb);                        // :77
    T denominator = norm_a * norm_b * (norm_a * norm_b + xa * xb + ya * yb + za * zb);  // :78
    T c1, c2, c3;
    seg_cross(xa, ya, za, xb, yb, zb, c1, c2, c3);                                  // :81
    T norm2_orth = c1 * c1 + c2 * c2 + c3 * c3;                                     // :82
    const T norm2_r0 = (xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb);  // :86
    // :80, :85, :89 `if (x <= PRECISION_UI) return` (NaN does not return)
    const bool skip = (denominator <= prec) | (norm2_orth <= prec) | (norm2_r0 <= prec);
    norm2_orth = norm2_orth / norm2_r0;                                             // :93
    T Kv;
    if constexpr (Reg == kRegNone) {  // :96-97
      Kv = T(1);
    } else if constexpr (Reg == kRegRankine) {  // :98-104
      const T r_bar2 = norm2_orth / (reg * reg);
      Kv = (r_bar2 < T(1)) ? r_bar2 : T(1);
    } else if constexpr (Reg == kRegLambOseen) {  // :105-112
      const T r_bar2 = norm2_orth / (reg * reg);
      const T exp_value = -lamb * r_bar2;
      if (exp_value < minexp) Kv = T(1);
      else Kv = T(1) - std::exp(exp_value);
    } else if constexpr (Reg == kRegVatistas) {  // :113-115
      const T r_bar2 = norm2_orth / (reg * reg);
      Kv = r_bar2 / std::sqrt(T(1) + r_bar2 * r_bar2);
    } else {  // :116-118
      static_assert(Reg == kRegOffset, "unknown segment regularization");
      Kv = T(1);
      denominator = denominator + reg * reg * norm2_r0;
    }
    kvo[t] = Kv;
    so[t] = skip ? T(-1) : norm_a + norm_b;
    deno[t] = denominator;
  };
  if constexpr (Reg == kRegLambOseen) {
    // the Fortran's branch around exp is kept (scalar libm exp, as gfortran calls it)
    for (int t = 0; t < g.nTpad; ++t) elem(t);
  } else {
    OFK_CS_SIMD
    for (int t = 0; t < g.nTpad; ++t) elem(t);
  }
}

// One replay work item: targets [b*kLanes, (b+1)*kLanes), every source in the pin's order.
template <class T>
inline void replay_block(const Geom<T>& g, const T* __restrict gam, T* __restrict V, const int b) {
  T v1[kLanes], v2[kLanes], v3[kLanes];
  for (int l = 0; l < kLanes; ++l) v1[l] = v2[l] = v3[l] = T(0);  // Vvar = 0.0_ReKi (FVW_Wings.f90:424), +0
  const int t0 = b * kLanes;
  const std::size_t np = std::size_t(g.nTpad);
  const T f = g.k.fourpi_inv;
  T cx[kLanes], cy[kLanes], cz[kLanes];
  for (int l = 0; l < kLanes; ++l) { cx[l] = g.cx[t0 + l]; cy[l] = g.cy[t0 + l]; cz[l] = g.cz[t0 + l]; }
  const int D = g.iNWStart;
  for (int kc = 0, r = 0; kc < g.nT; ++kc) {  // rings in the pin's order: (iW, iSpan) = kCP, then iDepth
    const T G = gam[kc];                         // Gamm = GammaLastIter(kCP) for both depths (FVW_Wings.f90:429)
    for (int d = 0; d < D; ++d, ++r)
    for (int sg = 0; sg < 4; ++sg) {
      const int j = 4 * r + sg;
      const SegEnds<T> e = seg_ends(g, j);
      const std::size_t o = std::size_t(j) * np + std::size_t(t0);
      const T* __restrict kv = g.kv + o;
      const T* __restrict s = g.s + o;
      const T* __restrict den = g.den + o;
      OFK_CS_SIMD
      for (int l = 0; l < kLanes; ++l) {
        T c1, c2, c3;
        seg_cross(cx[l] - e.pax, cy[l] - e.pay, cz[l] - e.paz, cx[l] - e.pbx, cy[l] - e.pby, cz[l] - e.pbz, c1, c2, c3);
        // Kv = SegGamma*fourpi_inv*Kv*(norm_a+norm_b)/denominator (:125), left to right
        T K = G * f;
        K = K * kv[l];
        K = K * s[l];
        K = K / den[l];
        const bool skip = s[l] < T(0);
        // Uind(1:3,icp) = Uind(1:3,icp) + Uindtmp(1:3), Uindtmp = Kv*crossprod (:126, :534)
        v1[l] += skip ? T(0) : K * c1;
        v2[l] += skip ? T(0) : K * c2;
        v3[l] += skip ? T(0) : K * c3;
      }
    }
  }
  for (int l = 0; l < kLanes; ++l) {
    const int t = t0 + l;
    if (t < g.nT) {
      V[3 * t] = v1[l];
      V[3 * t + 1] = v2[l];
      V[3 * t + 2] = v3[l];
    }
  }
}

template <int Reg, class T>
void build_reg(const Geom<T>& g, const bool par) {
  if (par) {
    const Geom<T> gg = g;
    Kokkos::parallel_for("ofk_cs_build", Kokkos::RangePolicy<HostExec>(0, g.nSrc),
                         [=](const int j) { seg_geom<Reg>(gg, j); });
    HostExec().fence("ofk_cs_build");
  } else {
    for (int j = 0; j < g.nSrc; ++j) seg_geom<Reg>(g, j);
  }
}

template <class T>
void build_t(const Geom<T>& g, const bool par) {
  switch (g.reg) {
    case kRegNone: build_reg<kRegNone>(g, par); break;
    case kRegRankine: build_reg<kRegRankine>(g, par); break;
    case kRegLambOseen: build_reg<kRegLambOseen>(g, par); break;
    case kRegVatistas: build_reg<kRegVatistas>(g, par); break;
    case kRegOffset: build_reg<kRegOffset>(g, par); break;
    default: break;  // rejected by the entry TU (the Fortran STOPs, :119-123)
  }
}

template <class T>
void vvar_t(const Geom<T>& g, const T* gam, T* V, const bool par) {
  const int nb = g.nTpad / kLanes;
  if (par) {
    const Geom<T> gg = g;
    Kokkos::parallel_for("ofk_cs_vvar", Kokkos::RangePolicy<HostExec>(0, nb),
                         [=](const int b) { replay_block(gg, gam, V, b); });
    HostExec().fence("ofk_cs_vvar");
  } else {
    for (int b = 0; b < nb; ++b) replay_block(g, gam, V, b);
  }
}

void build(const Geom<double>& g, bool par) { build_t(g, par); }
void build(const Geom<float>& g, bool par) { build_t(g, par); }
void vvar(const Geom<double>& g, const double* gam, double* V, bool par) { vvar_t(g, gam, V, par); }
void vvar(const Geom<float>& g, const float* gam, float* V, bool par) { vvar_t(g, gam, V, par); }

}  // namespace OFK_CS_NS
}  // namespace ofk_cs
