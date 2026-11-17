// Unit biotsavart-direct: loop shapes and the per-mode dispatcher (architecture.md §2-3).
// Included ONLY by biotsavart_direct_strict.cpp and biotsavart_direct_fast.cpp, which define
//   OFK_BSD_MODE_NS     strict | fast   (namespace: keeps the two instantiations apart)
//   OFK_BSD_STRICT      1 | 0
// and are compiled with FMA contraction off | on. Both TUs see exactly this source.
//
// Mode = (summation order, FMA). Everything else is shared:
//  - layout: the sources are packed once per call into SoA columns (a pure copy, values
//    unchanged), so the inner loop reads unit-stride arrays in both modes;
//  - order:
//      strict: kLanes = 1. Per target, sources in Fortran order, one accumulator triple,
//              one add into U -> identical arithmetic to the Fortran (class B on the CPU);
//      fast  : kLanes = 8 on host backends: source j goes to partial sum (j mod 8) in
//              blocks of 8, the tail fills lanes 0..r-1, then a fixed pairwise tree. The
//              lane count is a compile-time constant, so the result does not depend on the
//              thread count or the SIMD width. kLanes = 1 on device backends (parallelism
//              comes from targets; the order is then the strict order, only FMA differs).
//              The "skinny" team shape is fast-only (team size fixed by OFK_BSD_TEAM).
//  - "ok ? o : 0" accumulation: bitwise equal to the Fortran's skip, because the
//    accumulators start at +0 and can never become -0 (x + (+0) == x for every x != -0).
#pragma once

#include "biotsavart_direct_body.hpp"

#include <cstddef>

#ifndef OFK_BSD_STRICT
#error "define OFK_BSD_STRICT (1 = strict, 0 = fast) before including biotsavart_direct_impl.hpp"
#endif

// `omp simd` on the lane loop only asserts that the lanes are independent (each lane owns its
// accumulator, there is no reduction clause), so it does not change any result.
#if !OFK_BSD_STRICT && defined(_OPENMP) && !defined(__CUDACC__)
#define OFK_BSD_SIMD _Pragma("omp simd")
#else
#define OFK_BSD_SIMD
#endif

#define OFK_BSD_STR2(x) #x
#define OFK_BSD_STR(x) OFK_BSD_STR2(x)
#define OFK_BSD_LABEL(name) ("ofk_bsd_" name "_" OFK_BSD_STR(OFK_BSD_MODE_NS))

namespace ofk_bsd {
namespace OFK_BSD_MODE_NS {

constexpr int kLanes = kStrict ? 1 : (kOnHost ? 8 : 1);
static_assert((kLanes & (kLanes - 1)) == 0, "kLanes must be a power of two");

// Fixed pairwise tree over the lanes: lane i += lane i+w for w = L/2, L/4, ..., 1.
template <int L, class T>
KOKKOS_FORCEINLINE_FUNCTION T lanes_sum(T* a) {
  for (int w = L / 2; w >= 1; w /= 2)
    for (int i = 0; i < w; ++i) a[i] = a[i] + a[i + w];
  return a[0];
}

// ------------------------------------------------------------------------------------------
// ui_seg. SoA columns: x1 y1 z1 x2 y2 z2 gamma regparam (regparam = 0 for kRegNone: not read).
constexpr int kSegCols = 8;

// Column pointers of the SoA block (unit stride, so the lane loop vectorizes).
template <class T, int N>
struct Cols {
  const T* c[N];
  Cols(const T* base, std::size_t ld) {
    for (int i = 0; i < N; ++i) c[i] = base + i * ld;
  }
};

template <int Reg, class T>
KOKKOS_FORCEINLINE_FUNCTION void seg_acc(const T cx, const T cy, const T cz,
                                         const Cols<T, kSegCols>& s, const int j,
                                         const Consts<T>& k, T& u1, T& u2, T& u3) {
  const T xa = cx - s.c[0][j], ya = cy - s.c[1][j], za = cz - s.c[2][j];
  const T xb = cx - s.c[3][j], yb = cy - s.c[4][j], zb = cz - s.c[5][j];
  T o1 = T(0), o2 = T(0), o3 = T(0);
  const bool ok = seg_term<Reg, T>(xa, ya, za, xb, yb, zb, s.c[6][j], s.c[7][j], k, o1, o2, o3);
  u1 = u1 + (ok ? o1 : T(0));
  u2 = u2 + (ok ? o2 : T(0));
  u3 = u3 + (ok ? o3 : T(0));
}

// Pack the AoS arrays of `a` into the SoA block S (rows 0..nSeg-1 of a buffer whose row
// capacity may be larger: unit induced-velocities keeps S in a persistent per-instance handle).
template <class T>
void pack_seg_into(const SegArgs<T>& a, const DV2<T>& S) {
  if (a.nSeg == 0) return;
  const auto P = a.P;
  const auto C = a.C;
  const auto G = a.G;
  const auto R = a.R;
  const int pOff = a.pOff;
  const bool hasR = R.extent(0) > 0;
  Kokkos::parallel_for(
      OFK_BSD_LABEL("ui_seg_pack"), Kokkos::RangePolicy<ExecSpace>(0, a.nSeg),
      KOKKOS_LAMBDA(const int j) {
        const int p1 = C(0, j) - pOff, p2 = C(1, j) - pOff;  // 1-based -> packed 0-based
        S(j, 0) = P(0, p1); S(j, 1) = P(1, p1); S(j, 2) = P(2, p1);
        S(j, 3) = P(0, p2); S(j, 4) = P(1, p2); S(j, 5) = P(2, p2);
        S(j, 6) = G(j);
        S(j, 7) = hasR ? R(j) : T(0);
      });
}

template <class T>
DV2<T> pack_seg(const SegArgs<T>& a) {
  DV2<T> S(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ofk_bsd_seg_soa"), a.nSeg, kSegCols);
  pack_seg_into(a, S);
  return S;
}

// Square: RangePolicy over targets, kLanes partial sums over sources.
template <int Reg, class T>
void seg_square(const SegArgs<T>& a, const DV2<T>& S) {
  const int n = a.nSeg;
  const auto CP = a.CP;
  const auto U = a.U;
  const Consts<T> k = a.k;
  const Cols<T, kSegCols> s(S.data(), S.stride(1));
  Kokkos::parallel_for(
      OFK_BSD_LABEL("ui_seg_square"), Kokkos::RangePolicy<ExecSpace>(0, a.nCP),
      KOKKOS_LAMBDA(const int icp) {
        const T cx = CP(0, icp), cy = CP(1, icp), cz = CP(2, icp);
        const auto sl = s;  // column pointers in registers, not re-read from the closure
        T u1[kLanes], u2[kLanes], u3[kLanes];
        for (int l = 0; l < kLanes; ++l) { u1[l] = T(0); u2[l] = T(0); u3[l] = T(0); }
        const int nFull = n - n % kLanes;
        for (int j0 = 0; j0 < nFull; j0 += kLanes) {
          OFK_BSD_SIMD
          for (int l = 0; l < kLanes; ++l)
            seg_acc<Reg, T>(cx, cy, cz, sl, j0 + l, k, u1[l], u2[l], u3[l]);
        }
        for (int l = 0; l < n - nFull; ++l)
          seg_acc<Reg, T>(cx, cy, cz, sl, nFull + l, k, u1[l], u2[l], u3[l]);
        const T s1 = lanes_sum<kLanes>(u1), s2 = lanes_sum<kLanes>(u2), s3 = lanes_sum<kLanes>(u3);
        U(0, icp) = U(0, icp) + s1;
        U(1, icp) = U(1, icp) + s2;
        U(2, icp) = U(2, icp) + s3;
      });
}

template <class T>
struct V3 {
  T x, y, z;
  KOKKOS_INLINE_FUNCTION V3() : x(0), y(0), z(0) {}
  KOKKOS_INLINE_FUNCTION V3& operator+=(const V3& o) {
    x += o.x; y += o.y; z += o.z;
    return *this;
  }
};

}  // namespace OFK_BSD_MODE_NS
}  // namespace ofk_bsd

namespace Kokkos {
template <class T>
struct reduction_identity<ofk_bsd::OFK_BSD_MODE_NS::V3<T>> {
  KOKKOS_FORCEINLINE_FUNCTION static ofk_bsd::OFK_BSD_MODE_NS::V3<T> sum() {
    return ofk_bsd::OFK_BSD_MODE_NS::V3<T>();
  }
};
}  // namespace Kokkos

namespace ofk_bsd {
namespace OFK_BSD_MODE_NS {

// Skinny (fast only): league over (target, chunk of sources), a team reduction per chunk with
// a FIXED team size, then an ordered sum over chunks. Deterministic for a fixed backend.
// `scratch` (optional): persistent buffer for the per-chunk partial sums, >= 3*nChunk*nCP values;
// if it is too small, a temporary is allocated (the stateless unit-1 entry points).
template <int Reg, class T>
void seg_skinny(const SegArgs<T>& a, const DV2<T>& S, const DV1<T>& scratch) {
  using Policy = Kokkos::TeamPolicy<ExecSpace>;
  using Member = typename Policy::member_type;
  const int nCP = a.nCP, nSeg = a.nSeg;
  const int chunk = a.tune.chunk;
  const int nChunk = nSeg > 0 ? (nSeg + chunk - 1) / chunk : 0;
  const int team = kOnHost ? 1 : a.tune.team;
  const auto CP = a.CP;
  const auto U = a.U;
  const Consts<T> k = a.k;
  const Cols<T, kSegCols> s(S.data(), S.stride(1));
  const std::size_t nPartial = std::size_t(nChunk) * std::size_t(nCP);
  const DV1<T> partial = scratch.extent(0) >= 3 * nPartial
                             ? scratch
                             : DV1<T>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ofk_bsd_partial"),
                                      3 * nPartial);
  if (nChunk > 0)
    Kokkos::parallel_for(
        OFK_BSD_LABEL("ui_seg_skinny"), Policy(nCP * nChunk, team),
        KOKKOS_LAMBDA(const Member& m) {
          const int lr = m.league_rank();
          const int icp = lr % nCP, ich = lr / nCP;  // neighbouring teams share a chunk
          const int b = ich * chunk, e = (b + chunk < nSeg) ? b + chunk : nSeg;
          const T cx = CP(0, icp), cy = CP(1, icp), cz = CP(2, icp);
          V3<T> sum;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(m, b, e),
              [&](const int j, V3<T>& acc) { seg_acc<Reg, T>(cx, cy, cz, s, j, k, acc.x, acc.y, acc.z); },
              Kokkos::Sum<V3<T>>(sum));
          Kokkos::single(Kokkos::PerTeam(m), [&]() {
            const std::size_t q = 3 * (std::size_t(ich) * nCP + icp);
            partial(q) = sum.x; partial(q + 1) = sum.y; partial(q + 2) = sum.z;
          });
        });
  Kokkos::parallel_for(
      OFK_BSD_LABEL("ui_seg_skinny_sum"), Kokkos::RangePolicy<ExecSpace>(0, nCP),
      KOKKOS_LAMBDA(const int icp) {
        T u1 = 0, u2 = 0, u3 = 0;
        for (int ich = 0; ich < nChunk; ++ich) {
          const std::size_t q = 3 * (std::size_t(ich) * nCP + icp);
          u1 = u1 + partial(q); u2 = u2 + partial(q + 1); u3 = u3 + partial(q + 2);
        }
        U(0, icp) = U(0, icp) + u1;
        U(1, icp) = U(1, icp) + u2;
        U(2, icp) = U(2, icp) + u3;
      });
}

inline Shape resolve_shape(const Tuning& t, int nCP, int nSrc, bool skinny_ok) {
  if (kStrict || !skinny_ok) return Shape::kSquare;  // strict: Fortran order, square only
  if (t.shape == Shape::kSquare || t.shape == Shape::kSkinny) return t.shape;
  return (!kOnHost && nCP <= t.skinny_max_cp && nSrc >= t.skinny_min_seg) ? Shape::kSkinny
                                                                         : Shape::kSquare;
}

// Evaluate a pre-packed SoA block S (rows 0..a.nSeg-1) at a.CP, accumulating into a.U.
// a.P/a.C/a.G/a.R are not read.
template <class T>
Shape run_seg_soa(const SegArgs<T>& a, const DV2<T>& S, const DV1<T>& scratch) {
  const Shape sh = resolve_shape(a.tune, a.nCP, a.nSeg, true);
#define OFK_BSD_SEG_CASE(REG)                                       \
  case REG:                                                         \
    if (sh == Shape::kSkinny) seg_skinny<REG, T>(a, S, scratch);    \
    else seg_square<REG, T>(a, S);                                  \
    break;
  switch (a.reg) {
    OFK_BSD_SEG_CASE(kRegNone)
    OFK_BSD_SEG_CASE(kRegRankine)
    OFK_BSD_SEG_CASE(kRegLambOseen)
    OFK_BSD_SEG_CASE(kRegVatistas)
    OFK_BSD_SEG_CASE(kRegOffset)
  }
#undef OFK_BSD_SEG_CASE
  return sh;
}

template <class T>
Shape run_seg(const SegArgs<T>& a) {
  return run_seg_soa(a, pack_seg(a), DV1<T>());
}

// ------------------------------------------------------------------------------------------
// ui_part_nograd. SoA columns: px py pz alpha_x alpha_y alpha_z regparam (0 for kRegNone).
constexpr int kPartCols = 7;

template <int Reg, class T>
KOKKOS_FORCEINLINE_FUNCTION void part_acc(const T cx, const T cy, const T cz,
                                          const Cols<T, kPartCols>& s, const int j,
                                          const Consts<T>& k, T& u1, T& u2, T& u3) {
  const T dx = cx - s.c[0][j], dy = cy - s.c[1][j], dz = cz - s.c[2][j];
  T o1 = T(0), o2 = T(0), o3 = T(0);
  const bool ok = part_term<Reg, T>(dx, dy, dz, s.c[3][j], s.c[4][j], s.c[5][j], s.c[6][j], k, o1, o2, o3);
  u1 = u1 + (ok ? o1 : T(0));
  u2 = u2 + (ok ? o2 : T(0));
  u3 = u3 + (ok ? o3 : T(0));
}

template <class T>
void pack_part_into(const PartArgs<T>& a, const DV2<T>& S) {
  if (a.nPart == 0) return;
  const auto P = a.P;
  const auto A = a.A;
  const auto R = a.R;
  const bool hasR = R.extent(0) > 0;
  Kokkos::parallel_for(
      OFK_BSD_LABEL("ui_part_pack"), Kokkos::RangePolicy<ExecSpace>(0, a.nPart),
      KOKKOS_LAMBDA(const int j) {
        S(j, 0) = P(0, j); S(j, 1) = P(1, j); S(j, 2) = P(2, j);
        S(j, 3) = A(0, j); S(j, 4) = A(1, j); S(j, 5) = A(2, j);
        S(j, 6) = hasR ? R(j) : T(0);
      });
}

template <class T>
DV2<T> pack_part(const PartArgs<T>& a) {
  DV2<T> S(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ofk_bsd_part_soa"), a.nPart, kPartCols);
  pack_part_into(a, S);
  return S;
}

template <int Reg, class T>
void part_square(const PartArgs<T>& a, const DV2<T>& S) {
  const int n = a.nPart;
  const auto CP = a.CP;
  const auto U = a.U;
  const Consts<T> k = a.k;
  const Cols<T, kPartCols> s(S.data(), S.stride(1));
  Kokkos::parallel_for(
      OFK_BSD_LABEL("ui_part_square"), Kokkos::RangePolicy<ExecSpace>(0, a.nCP),
      KOKKOS_LAMBDA(const int icp) {
        const T cx = CP(0, icp), cy = CP(1, icp), cz = CP(2, icp);
        const auto sl = s;  // column pointers in registers, not re-read from the closure
        T u1[kLanes], u2[kLanes], u3[kLanes];
        for (int l = 0; l < kLanes; ++l) { u1[l] = T(0); u2[l] = T(0); u3[l] = T(0); }
        const int nFull = n - n % kLanes;
        for (int j0 = 0; j0 < nFull; j0 += kLanes) {
          OFK_BSD_SIMD
          for (int l = 0; l < kLanes; ++l)
            part_acc<Reg, T>(cx, cy, cz, sl, j0 + l, k, u1[l], u2[l], u3[l]);
        }
        for (int l = 0; l < n - nFull; ++l)
          part_acc<Reg, T>(cx, cy, cz, sl, nFull + l, k, u1[l], u2[l], u3[l]);
        const T s1 = lanes_sum<kLanes>(u1), s2 = lanes_sum<kLanes>(u2), s3 = lanes_sum<kLanes>(u3);
        U(0, icp) = U(0, icp) + s1;
        U(1, icp) = U(1, icp) + s2;
        U(2, icp) = U(2, icp) + s3;
      });
}

template <class T>
Shape run_part_soa(const PartArgs<T>& a, const DV2<T>& S) {
  switch (a.reg) {
    case kRegNone: part_square<kRegNone, T>(a, S); break;
    case kRegPartExp: part_square<kRegPartExp, T>(a, S); break;
    case kRegPartCompact: part_square<kRegPartCompact, T>(a, S); break;
  }
  return Shape::kSquare;
}

template <class T>
Shape run_part(const PartArgs<T>& a) {
  return run_part_soa(a, pack_part(a));
}

template Shape run_seg_soa<double>(const SegArgs<double>&, const DV2<double>&, const DV1<double>&);
template Shape run_seg_soa<float>(const SegArgs<float>&, const DV2<float>&, const DV1<float>&);
template Shape run_part_soa<double>(const PartArgs<double>&, const DV2<double>&);
template Shape run_part_soa<float>(const PartArgs<float>&, const DV2<float>&);
template Shape run_seg<double>(const SegArgs<double>&);
template Shape run_seg<float>(const SegArgs<float>&);
template Shape run_part<double>(const PartArgs<double>&);
template Shape run_part<float>(const PartArgs<float>&);

}  // namespace OFK_BSD_MODE_NS
}  // namespace ofk_bsd
