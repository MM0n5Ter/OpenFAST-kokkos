// Unit circulation-solve (OLAF unit 4a, F10-B): types shared by the mode TUs and the entry TU.
// No arithmetic here. The physics is in circulation_body.hpp; the handle and the C API in
// circulation.cpp; the design in tests/circulation-solve/unit.md.
//
// What is computed (pin 12fa3fc): the Vvar block of Wings_ComputeCirculationPolarData
// (FVW_Wings.f90:424-447), i.e. for every lifting-line control point (target) the velocity
// induced by the first iNWStart near-wake panels of every wing, each a vortex ring of four
// segments with the circulation GammaLastIter of its spanwise station:
//   do iW / iSpan / iDepth = 1..iNWStart       (ring r, Gamma = GammaLastIter(kCP))
//     ui_quad_n1 on all targets                (FVW_BiotSavart.f90:513-552)
//       segments 1-2, 3-4, 2-3, 4-1            (ui_seg_11, FVW_BiotSavart.f90:52-128)
// Everything in ui_seg_11 before :125 is independent of Gamma. The geometry cache holds it per
// (segment, target) and is built once per solve; each iteration replays only
//   Kv = SegGamma*fourpi_inv*Kv*(norm_a+norm_b)/denominator ; Uind = Kv*crossprod
// and accumulates per target in the pin's order (ring-major, depth, segment), so the result is
// bitwise the pin's (class B) when contraction is off.
#pragma once

#include <cstddef>

namespace ofk_cs {

// Targets per replay work item (one SIMD block). The per-target order does not depend on it.
constexpr int kLanes = 8;

// FVW_BiotSavart.f90 regularization ids (:21-25)
constexpr int kRegNone = 0, kRegRankine = 1, kRegLambOseen = 2, kRegVatistas = 3, kRegOffset = 4;

// The constants of ui_seg_11, as the Fortran evaluated them (passed from the shim: they are
// ReKi parameters or ReKi literals of FVW_BiotSavart, criteria R3).
template <class T>
struct Consts {
  T precision_ui;   // PRECISION_UI = epsilon(1.0_ReKi)/100 (:12)
  T fourpi_inv;     // fourpi_inv = 0.25_ReKi/ACOS(-1.0_ReKi) (:32)
  T min_exp_value;  // MIN_EXP_VALUE = -10.0_ReKi (:14)
  T lamb_oseen;     // the literal 1.25643_ReKi (:107)
};

// One solve's geometry and cache, as raw host pointers into the handle's Views.
//   targets  cx/cy/cz[nTpad]: m%W(:)%CP flattened over wings (FVW_Wings.f90:437-443), padded to a
//            multiple of kLanes with copies of target 0 (padded lanes are computed and discarded);
//   rings    ring[13*r]: P1(3) P2(3) P3(3) P4(3) RegParam, r = kCP*iNWStart + (iDepth-1) with
//            P1 = r_NW(:,s,d), P2 = r_NW(:,s+1,d), P3 = r_NW(:,s+1,d+1), P4 = r_NW(:,s,d+1) (:432-435),
//            RegParam = Eps_NW(1,s,d) (row 1 for all four segments, :441);
//   cache    [j*nTpad + t] for source segment j = 4*r + k (k: 0 = 1-2, 1 = 3-4, 2 = 2-3, 3 = 4-1):
//            kv  = Kv after the regularization select (:95-123),
//            s   = norm_a+norm_b, or -1 where ui_seg_11 returns early (its three guards; a
//                  computed norm_a+norm_b is never negative: sqrt of a sum of squares is +0 or more,
//                  or NaN, and NaN is not < 0),
//            den = denominator (with the Offset term, :118).
//            The cross product (:89) is not stored: the replay recomputes it (same expressions).
template <class T>
struct Geom {
  int nT = 0, nTpad = 0, nRing = 0, nSrc = 0, iNWStart = 0, reg = 0;
  Consts<T> k{};
  const T* cx = nullptr;
  const T* cy = nullptr;
  const T* cz = nullptr;
  const T* ring = nullptr;
  T* kv = nullptr;
  T* s = nullptr;
  T* den = nullptr;
};

// Per-mode entry points, defined once per mode TU (circulation_{strict,fast,fastfma}.cpp) from the
// same body. par = run the Kokkos host-space parallel_for (else the same work items in a plain
// loop, P4). gam = GammaLastIter(1:nT); V = Vvar(1:3,1:nT), overwritten.
#define OFK_CS_DECLARE_MODE(NS)                                                           \
  namespace NS {                                                                          \
  void build(const Geom<double>& g, bool par);                                            \
  void build(const Geom<float>& g, bool par);                                             \
  void vvar(const Geom<double>& g, const double* gam, double* V, bool par);               \
  void vvar(const Geom<float>& g, const float* gam, float* V, bool par);                  \
  }
OFK_CS_DECLARE_MODE(strict)
OFK_CS_DECLARE_MODE(fast)
OFK_CS_DECLARE_MODE(fastfma)
#undef OFK_CS_DECLARE_MODE

}  // namespace ofk_cs
