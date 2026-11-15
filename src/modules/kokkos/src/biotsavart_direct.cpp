// Unit biotsavart-direct: C entry points (extern "C"), argument checks, host<->device views
// and the strict/fast switch. No physics here: the bodies are in biotsavart_direct_body.hpp,
// the loop shapes in biotsavart_direct_impl.hpp, compiled once per mode in
// biotsavart_direct_strict.cpp and biotsavart_direct_fast.cpp.
//   FVW_BiotSavart::ui_seg          (pin modules/aerodyn/src/FVW_BiotSavart.f90:135-357)
//   FVW_BiotSavart::ui_part_nograd  (pin modules/aerodyn/src/FVW_BiotSavart.f90:360-455)
// Contract, equivalence classes and audit list: tests/biotsavart-direct/unit.md (OpenFAST-kokkos).
//
// Rules applied (notes/build-kokkos-fortran.md §6, notes/architecture.md §2-3):
//  - SegConnct values are 1-based column indices into SegPoints; they are offset at use;
//  - only the target columns of the given range are read or written; Uind is accumulated;
//  - every exception is caught at the extern "C" boundary (the Fortran routines have no
//    ErrStat, so a failure is fatal, as the Fortran STOP is).
// Environment (read once, process-wide):
//   OFK_MODE = strict | fast (ofk_runtime.hpp)
//   OFK_BSD_SHAPE = auto | square | skinny (fast mode only; strict is always square)
//   OFK_BSD_SKINNY_MAX_CP (1024), OFK_BSD_SKINNY_MIN_SEG (4096), OFK_BSD_CHUNK (1024),
//   OFK_BSD_TEAM (128): skinny-shape thresholds and sizes
//   OFK_BSD_PROFILE = 1: fence and print alloc+H2D / kernel / D2H times per call to stderr
#include "biotsavart_direct.hpp"

#include <Kokkos_Core.hpp>

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "biotsavart_direct_common.hpp"
#include "biotsavart_direct_consts.hpp"
#include "ofk_runtime.hpp"

namespace ofk_bsd {

template <class T>
using HostU2 = Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template <class T>
using HostU1 = Kokkos::View<T*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

struct Config {
  Tuning tune;
  bool profile = false;
};

const Config& config() {
  static const Config c = [] {
    Config r;
    auto geti = [](const char* n, int d) {
      const char* v = std::getenv(n);
      return (v && *v) ? std::atoi(v) : d;
    };
    if (const char* s = std::getenv("OFK_BSD_SHAPE")) {
      if (!std::strcmp(s, "square")) r.tune.shape = Shape::kSquare;
      else if (!std::strcmp(s, "skinny")) r.tune.shape = Shape::kSkinny;
    }
    r.tune.skinny_max_cp = geti("OFK_BSD_SKINNY_MAX_CP", r.tune.skinny_max_cp);
    r.tune.skinny_min_seg = geti("OFK_BSD_SKINNY_MIN_SEG", r.tune.skinny_min_seg);
    r.tune.chunk = geti("OFK_BSD_CHUNK", r.tune.chunk);
    if (r.tune.chunk < 1) r.tune.chunk = 1;
    r.tune.team = geti("OFK_BSD_TEAM", r.tune.team);
    if (r.tune.team < 1) r.tune.team = 1;
    r.profile = geti("OFK_BSD_PROFILE", 0) != 0;
    return r;
  }();
  return c;
}

struct Profile {
  bool on;
  Kokkos::Timer timer;
  double t[3] = {0, 0, 0};
  int stage = 0;
  explicit Profile(bool o) : on(o) {}
  void mark() {  // end of a stage: alloc+H2D, kernel (incl. SoA pack), D2H
    if (!on) return;
    Kokkos::fence("ofk_bsd_profile");
    t[stage++] = timer.seconds() * 1e3;
    timer.reset();
  }
};

const Tuning& tuning() { return config().tune; }

const char* reg_name_seg(int r) {
  static const char* n[] = {"None", "Rankine", "LambOseen", "Vatistas", "Offset"};
  return (r >= 0 && r <= 4) ? n[r] : "?";
}
const char* reg_name_part(int r) {
  static const char* n[] = {"None", "Exp", "Compact"};
  return (r >= 0 && r <= 2) ? n[r] : "?";
}
const char* shape_name(Shape s) { return s == Shape::kSkinny ? "skinny" : "square"; }

template <class T>
void ui_seg(int iCPStart, int iCPEnd, int ldCP, const T* CPs, int iSegStart, int iSegEnd,
            int ldPts, int nPts, const T* SegPoints, int ldC, const int* SegConnct,
            const T* SegGamma, int Reg, const T* RegParam, int ldU, T* Uind) {
  const int nCP = iCPEnd - iCPStart + 1;
  if (nCP <= 0) return;  // Fortran: empty CP loop, nothing touched
  if (Reg < kRegNone || Reg > kRegOffset)
    throw std::invalid_argument("Unknown RegFunction for segment " + std::to_string(Reg));
  if (iCPStart < 1 || ldCP < 3 || ldU < 3 || ldPts < 3 || ldC < 2)
    throw std::invalid_argument("bad array extents");
  const int nSeg = (iSegEnd >= iSegStart) ? iSegEnd - iSegStart + 1 : 0;
  if (nSeg > 0 && iSegStart < 1) throw std::invalid_argument("iSegStart < 1");

  // Columns of SegPoints referenced by the selected segments (SegConnct is 1-based). Only
  // that block is copied to the device; out-of-range indices are an error, not UB.
  int pMin = 1, pMax = 0;
  if (nSeg > 0) {
    pMin = INT_MAX; pMax = INT_MIN;
    const int* c0 = SegConnct + std::size_t(iSegStart - 1) * ldC;
    for (int is = 0; is < nSeg; ++is)
      for (int r = 0; r < 2; ++r) {
        const int v = c0[std::size_t(is) * ldC + r];
        pMin = v < pMin ? v : pMin;
        pMax = v > pMax ? v : pMax;
      }
    if (pMin < 1 || pMax > nPts)
      throw std::out_of_range("SegConnct index outside 1..size(SegPoints,2): [" +
                              std::to_string(pMin) + "," + std::to_string(pMax) + "]");
  }
  const int nP = pMax - pMin + 1;

  const auto& cfg = config();
  const ofk::Mode mode = ofk::mode();
  Profile prof(cfg.profile);

  HostU2<const T> hCP(CPs + std::size_t(iCPStart - 1) * ldCP, ldCP, nCP);
  HostU2<const T> hP(SegPoints + std::size_t(pMin - 1) * ldPts, ldPts, nP);
  HostU2<const int> hC(SegConnct + std::size_t(nSeg > 0 ? iSegStart - 1 : 0) * ldC, ldC, nSeg);
  HostU1<const T> hG(SegGamma + (nSeg > 0 ? iSegStart - 1 : 0), nSeg);
  HostU1<const T> hR(RegParam + (nSeg > 0 ? iSegStart - 1 : 0), Reg == kRegNone ? 0 : nSeg);
  HostU2<T> hU(Uind + std::size_t(iCPStart - 1) * ldU, ldU, nCP);

  SegArgs<T> a;
  a.nCP = nCP; a.nSeg = nSeg; a.pOff = pMin; a.reg = Reg;
  a.k = make_consts<T>(T(0));  // MINNORM is not used by ui_seg
  a.tune = cfg.tune;
  // Host backends: these are the same unmanaged Views, nothing is copied.
  a.CP = Kokkos::create_mirror_view_and_copy(MemSpace{}, hCP);
  a.P = Kokkos::create_mirror_view_and_copy(MemSpace{}, hP);
  a.C = Kokkos::create_mirror_view_and_copy(MemSpace{}, hC);
  a.G = Kokkos::create_mirror_view_and_copy(MemSpace{}, hG);
  a.R = Kokkos::create_mirror_view_and_copy(MemSpace{}, hR);
  auto U = Kokkos::create_mirror_view_and_copy(MemSpace{}, hU);  // copy in: we accumulate
  a.U = U;
  prof.mark();
  const Shape sh = (mode == ofk::Mode::strict) ? strict::run_seg<T>(a) : fast::run_seg<T>(a);
  prof.mark();
  Kokkos::deep_copy(hU, U);  // fences; a no-op when U is hU
  prof.mark();
  if (prof.on)
    std::fprintf(stderr,
                 "[ofk_bsd] ui_seg %s %s %s nCP=%d nSeg=%d alloc+H2D=%.3f ms kernel=%.3f ms "
                 "D2H=%.3f ms\n",
                 ofk::mode_name(mode), reg_name_seg(Reg), shape_name(sh), nCP, nSeg, prof.t[0],
                 prof.t[1], prof.t[2]);
}

template <class T>
void ui_part_nograd(int nCPs, int ldCP, const T* CPs, int nPart, int ldP, const T* Part, int ldA,
                    const T* Alpha, int Reg, const T* RegParam, int ldU, T* UIout, T minnorm) {
  if (nCPs <= 0) return;
  if (Reg < kRegNone || Reg > kRegPartCompact)
    throw std::invalid_argument("Wrong regularization function for particles " +
                                std::to_string(Reg));
  if (ldCP < 3 || ldP < 3 || ldA < 3 || ldU < 3) throw std::invalid_argument("bad array extents");
  if (nPart < 0) nPart = 0;

  const auto& cfg = config();
  const ofk::Mode mode = ofk::mode();
  Profile prof(cfg.profile);

  HostU2<const T> hCP(CPs, ldCP, nCPs);
  HostU2<const T> hP(Part, ldP, nPart);
  HostU2<const T> hA(Alpha, ldA, nPart);
  HostU1<const T> hR(RegParam, Reg == kRegNone ? 0 : nPart);  // not read for idRegNone
  HostU2<T> hU(UIout, ldU, nCPs);

  PartArgs<T> a;
  a.nCP = nCPs; a.nPart = nPart; a.reg = Reg;
  a.k = make_consts<T>(minnorm);
  a.tune = cfg.tune;
  a.CP = Kokkos::create_mirror_view_and_copy(MemSpace{}, hCP);
  a.P = Kokkos::create_mirror_view_and_copy(MemSpace{}, hP);
  a.A = Kokkos::create_mirror_view_and_copy(MemSpace{}, hA);
  a.R = Kokkos::create_mirror_view_and_copy(MemSpace{}, hR);
  auto U = Kokkos::create_mirror_view_and_copy(MemSpace{}, hU);
  a.U = U;
  prof.mark();
  if (mode == ofk::Mode::strict) strict::run_part<T>(a); else fast::run_part<T>(a);
  prof.mark();
  Kokkos::deep_copy(hU, U);
  prof.mark();
  if (prof.on)
    std::fprintf(stderr,
                 "[ofk_bsd] ui_part_nograd %s %s nCP=%d nPart=%d alloc+H2D=%.3f ms "
                 "kernel=%.3f ms D2H=%.3f ms\n",
                 ofk::mode_name(mode), reg_name_part(Reg), nCPs, nPart, prof.t[0], prof.t[1],
                 prof.t[2]);
}

[[noreturn]] void fatal(const char* where, const char* what) {
  std::fprintf(stderr, "[ERROR] %s: %s\n", where, what);
  std::fflush(stderr);
  std::abort();
}

}  // namespace ofk_bsd

// ---------------------------------------------------------------------------------------------
#define OFK_BSD_GUARD(NAME, CALL)                       \
  try {                                                 \
    ofk::ensure_initialized();                          \
    CALL;                                               \
  } catch (const std::exception& e) {                   \
    ofk_bsd::fatal(NAME, e.what());                     \
  } catch (...) {                                       \
    ofk_bsd::fatal(NAME, "unknown C++ exception");      \
  }

extern "C" {

void ofk_biotsavart_direct_ui_seg_r8(int iCPStart, int iCPEnd, int ldCP, const double* CPs,
                                     int iSegStart, int iSegEnd, int ldPts, int nPts,
                                     const double* SegPoints, int ldConnct, const int* SegConnct,
                                     const double* SegGamma, int RegFunction,
                                     const double* RegParam, int ldU, double* Uind) {
  OFK_BSD_GUARD("ofk_biotsavart_direct_ui_seg_r8",
                ofk_bsd::ui_seg<double>(iCPStart, iCPEnd, ldCP, CPs, iSegStart, iSegEnd, ldPts,
                                        nPts, SegPoints, ldConnct, SegConnct, SegGamma,
                                        RegFunction, RegParam, ldU, Uind))
}

void ofk_biotsavart_direct_ui_seg_r4(int iCPStart, int iCPEnd, int ldCP, const float* CPs,
                                     int iSegStart, int iSegEnd, int ldPts, int nPts,
                                     const float* SegPoints, int ldConnct, const int* SegConnct,
                                     const float* SegGamma, int RegFunction,
                                     const float* RegParam, int ldU, float* Uind) {
  OFK_BSD_GUARD("ofk_biotsavart_direct_ui_seg_r4",
                ofk_bsd::ui_seg<float>(iCPStart, iCPEnd, ldCP, CPs, iSegStart, iSegEnd, ldPts,
                                       nPts, SegPoints, ldConnct, SegConnct, SegGamma, RegFunction,
                                       RegParam, ldU, Uind))
}

void ofk_biotsavart_direct_ui_part_nograd_r8(int nCPs, int ldCP, const double* CPs, int nPart,
                                             int ldP, const double* Part, int ldA,
                                             const double* Alpha, int RegFunction,
                                             const double* RegParam, int ldU, double* UIout,
                                             double minnorm) {
  OFK_BSD_GUARD("ofk_biotsavart_direct_ui_part_nograd_r8",
                ofk_bsd::ui_part_nograd<double>(nCPs, ldCP, CPs, nPart, ldP, Part, ldA, Alpha,
                                                RegFunction, RegParam, ldU, UIout, minnorm))
}

void ofk_biotsavart_direct_ui_part_nograd_r4(int nCPs, int ldCP, const float* CPs, int nPart,
                                             int ldP, const float* Part, int ldA,
                                             const float* Alpha, int RegFunction,
                                             const float* RegParam, int ldU, float* UIout,
                                             float minnorm) {
  OFK_BSD_GUARD("ofk_biotsavart_direct_ui_part_nograd_r4",
                ofk_bsd::ui_part_nograd<float>(nCPs, ldCP, CPs, nPart, ldP, Part, ldA, Alpha,
                                               RegFunction, RegParam, ldU, UIout, minnorm))
}

}  // extern "C"
