// Unit circulation-solve: per-instance handles, the C API, mode/variant selection, the serial vs
// parallel_for threshold (P4), the iteration log of the FMA experiment and counters. No physics
// here (circulation_body.hpp). C API and environment: circulation.hpp.
#include "circulation.hpp"

#include "circulation_common.hpp"
#include "ofk_runtime.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace ofk_cs {
namespace {

// FMA policy of this unit's fast mode, decided by the FMA experiment (tests/circulation-solve/unit.md, rule fixed
// in advance): F1 passed R6 and R8 and its solve kernel is faster than F0's beyond the run-to-run noise, so the
// fast mode keeps FMA as architecture.md §3 says. 1 = F1 (contraction), 0 = F0 (no contraction: bitwise to strict
// and the pin; OFK_CS_FMA=0).
constexpr int kDefaultFma = 1;

struct Config {
  bool enabled = true;
  int mode = -1;  // -1: follow OFK_MODE; else ofk::Mode
  int fma = kDefaultFma;
  long long parMin = 0, parMinBuild = 0;
  bool profile = false;
  std::size_t maxHandles = 64;
  std::string itersPath;
};

long long env_ll(const char* name, long long dflt) {
  const char* s = std::getenv(name);
  return (s && *s) ? std::atoll(s) : dflt;
}

const Config& config() {
  static Config c = [] {
    Config r;
    if (const char* s = std::getenv("OFK_CS")) r.enabled = std::strcmp(s, "fortran") != 0;
    if (const char* s = std::getenv("OFK_CS_MODE")) {
      if (!std::strcmp(s, "strict")) r.mode = static_cast<int>(ofk::Mode::strict);
      else if (!std::strcmp(s, "fast")) r.mode = static_cast<int>(ofk::Mode::fast);
    }
    if (const char* s = std::getenv("OFK_CS_FMA")) r.fma = std::atoi(s) != 0 ? 1 : 0;
    // P4 thresholds in targets x segments per solve: the serial/parallel_for crossovers at 8 threads on
    // Xeon 8480+ (job 1925426: build 288-648 pairs, replay 1800-2592 pairs; tests/circulation-solve/unit.md, P4)
    r.parMin = env_ll("OFK_CS_PAR_MIN", 2048);
    r.parMinBuild = env_ll("OFK_CS_PAR_MIN_BUILD", 512);
    if (const char* s = std::getenv("OFK_CS_PROFILE")) r.profile = std::atoi(s) != 0;
    r.maxHandles = static_cast<std::size_t>(std::max<long long>(1, env_ll("OFK_CS_MAX_HANDLES", 64)));
    if (const char* s = std::getenv("OFK_CS_ITERS")) r.itersPath = s;
    return r;
  }();
  return c;
}

// test-hook overrides (ofk_cs_set_variant)
int g_mode_override = -1, g_fma_override = -1, g_par_override = -1;

using HostMem = Kokkos::HostSpace;

struct Stats {
  long long nAlloc = 0, bytes = 0, nSolve = 0, nIter = 0, nBuild = 0, nVvar = 0, nCap = 0;
  double tBuild = 0, tVvar = 0;
  bool lastPar = false;
};

template <class T>
struct Handle {
  using V = Kokkos::View<T*, HostMem>;
  const void* key = nullptr;
  V cx, cy, cz, ring, kv, s, den;
  Geom<T> g;
  int nWings = 0, wingsSeen = 0;
  bool built = false;
  int mode = 1, fma = 0;  // the variant of the current solve (fixed at begin)
  bool parBuild = false, parVvar = false;
  Stats st;
};

template <class T>
void ensure(Handle<T>& h, typename Handle<T>::V& v, std::size_t n, const char* name) {
  if (v.extent(0) >= n) return;
  v = typename Handle<T>::V(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), n);
  h.st.nAlloc++;
  h.st.bytes += static_cast<long long>(n * sizeof(T));
}

struct Entry {
  std::unique_ptr<Handle<double>> d;
  std::unique_ptr<Handle<float>> f;
  std::uint64_t lastUse = 0;
};
struct Registry {
  std::mutex mu;
  std::unordered_map<const void*, Entry> map;
  std::uint64_t clock = 0;
  bool hook = false;
};
Registry& registry() {
  static Registry* r = new Registry;  // never destroyed (no static-destruction race with Kokkos)
  return *r;
}

template <class T>
void report_one(const Handle<T>& h, const char* why) {
  const Stats& s = h.st;
  std::fprintf(stderr,
               "[ofk_cs] %s key=%p r%d: nT=%d nSrc=%d solves=%lld iters=%lld (cap reached %lld) builds=%lld "
               "replays=%lld allocs=%lld (%.2f MB) t[s] build=%.4f replay=%.4f per solve=%.3e\n",
               why, h.key, int(sizeof(T)), h.g.nT, h.g.nSrc, s.nSolve, s.nIter, s.nCap, s.nBuild, s.nVvar, s.nAlloc,
               s.bytes / 1e6, s.tBuild, s.tVvar, s.nSolve ? (s.tBuild + s.tVvar) / double(s.nSolve) : 0.0);
}
void report_entry(const Entry& e, const char* why) {
  if (e.d) report_one(*e.d, why);
  if (e.f) report_one(*e.f, why);
}

template <class T>
std::unique_ptr<Handle<T>>& slot(Entry& e);
template <>
std::unique_ptr<Handle<double>>& slot<double>(Entry& e) { return e.d; }
template <>
std::unique_ptr<Handle<float>>& slot<float>(Entry& e) { return e.f; }

template <class T>
Handle<T>& get(const void* key) {
  Registry& r = registry();
  std::lock_guard<std::mutex> lock(r.mu);
  if (!r.hook) {
    r.hook = true;
    Kokkos::push_finalize_hook([] {
      Registry& rr = registry();
      std::lock_guard<std::mutex> l(rr.mu);
      if (config().profile)
        for (auto& kv : rr.map) report_entry(kv.second, "finalize");
      rr.map.clear();  // Views must die before Kokkos::finalize
    });
  }
  auto it = r.map.find(key);
  if (it == r.map.end()) {
    if (r.map.size() >= config().maxHandles) {
      auto lru = r.map.begin();
      for (auto i = r.map.begin(); i != r.map.end(); ++i)
        if (i->second.lastUse < lru->second.lastUse) lru = i;
      if (config().profile) report_entry(lru->second, "evicted");
      r.map.erase(lru);
    }
    it = r.map.emplace(key, Entry{}).first;
  }
  it->second.lastUse = ++r.clock;
  auto& p = slot<T>(it->second);
  if (!p) {
    p = std::make_unique<Handle<T>>();
    p->key = key;
  }
  return *p;
}

template <class T>
Handle<T>* find(const void* key) {
  Registry& r = registry();
  std::lock_guard<std::mutex> lock(r.mu);
  auto it = r.map.find(key);
  return it == r.map.end() ? nullptr : slot<T>(it->second).get();
}

[[noreturn]] void fatal(const char* what) {
  std::fprintf(stderr, "[ofk_cs] FATAL: %s\n", what);
  std::fflush(stderr);
  std::abort();
}

void announce() {
  static std::once_flag once;
  std::call_once(once, [] {
    const Config& c = config();
    const int m = c.mode >= 0 ? c.mode : static_cast<int>(ofk::mode());
    std::printf("[ofk_cs] circulation-solve: mode=%s%s fma=%d host threads=%d par_min=%lld par_min_build=%lld%s\n",
                ofk::mode_name(static_cast<ofk::Mode>(m)), c.mode >= 0 ? " (OFK_CS_MODE)" : "",
                m == static_cast<int>(ofk::Mode::fast) ? c.fma : 0,
                Kokkos::DefaultHostExecutionSpace().concurrency(), c.parMin, c.parMinBuild,
                c.itersPath.empty() ? "" : " (OFK_CS_ITERS on)");
    std::fflush(stdout);
  });
}

// --- the FMA experiment's iteration log (OFK_CS_ITERS) ---------------------------------------
std::mutex g_log_mu;
std::FILE* g_log = nullptr;
long long g_log_n = 0;
void log_iters(int nIter, int maxIter) {
  const Config& c = config();
  if (c.itersPath.empty()) return;
  std::lock_guard<std::mutex> l(g_log_mu);
  if (!g_log) {
    g_log = std::fopen(c.itersPath.c_str(), "w");
    if (!g_log) return;
    std::fprintf(g_log, "# solve iterations cap (OFK_CS_ITERS, unit circulation-solve)\n");
  }
  std::fprintf(g_log, "%lld %d %d\n", ++g_log_n, nIter, maxIter);
  std::fflush(g_log);
}

using clk = std::chrono::steady_clock;
double secs(clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); }

template <class T>
void begin_t(const void* key, int nWings, int nCPtot, int iNWStart, int reg, const T* consts4) {
  ofk::ensure_initialized();
  announce();
  if (reg < kRegNone || reg > kRegOffset) {
    std::printf("Unknown SgmtReg %d\n", reg);  // as ui_seg_11 (FVW_BiotSavart.f90:120-121), which then STOPs
    fatal("unknown RegFunction");
  }
  if (nWings < 1 || nCPtot < 1 || iNWStart < 1) fatal("ofk_cs_begin: empty lifting line or iNWStart < 1");
  Handle<T>& h = get<T>(key);
  const Config& c = config();
  h.mode = g_mode_override >= 0 ? g_mode_override : (c.mode >= 0 ? c.mode : static_cast<int>(ofk::mode()));
  h.fma = (h.mode == static_cast<int>(ofk::Mode::fast)) ? (g_fma_override >= 0 ? g_fma_override : c.fma) : 0;
  Geom<T>& g = h.g;
  g.nT = nCPtot;
  g.nTpad = (nCPtot + kLanes - 1) / kLanes * kLanes;
  g.iNWStart = iNWStart;
  g.nRing = nCPtot * iNWStart;
  g.nSrc = 4 * g.nRing;
  g.reg = reg;
  g.k.precision_ui = consts4[0];
  g.k.fourpi_inv = consts4[1];
  g.k.min_exp_value = consts4[2];
  g.k.lamb_oseen = consts4[3];
  // Capacity: nCP_tot and iNWStart are fixed for a run, so this allocates on the first solve only (P6)
  const std::size_t nTp = std::size_t(g.nTpad), nC = std::size_t(g.nSrc) * nTp;
  ensure(h, h.cx, nTp, "ofk_cs_cx");
  ensure(h, h.cy, nTp, "ofk_cs_cy");
  ensure(h, h.cz, nTp, "ofk_cs_cz");
  ensure(h, h.ring, std::size_t(13) * std::size_t(g.nRing), "ofk_cs_ring");
  ensure(h, h.kv, nC, "ofk_cs_kv");
  ensure(h, h.s, nC, "ofk_cs_s");
  ensure(h, h.den, nC, "ofk_cs_den");
  g.cx = h.cx.data(); g.cy = h.cy.data(); g.cz = h.cz.data(); g.ring = h.ring.data();
  g.kv = h.kv.data(); g.s = h.s.data(); g.den = h.den.data();
  h.nWings = nWings;
  h.wingsSeen = 0;
  h.built = false;
  const long long pairs = static_cast<long long>(g.nT) * g.nSrc;
  const bool threads = Kokkos::DefaultHostExecutionSpace().concurrency() > 1;
  h.parBuild = g_par_override >= 0 ? g_par_override == 1 : (threads && pairs >= c.parMinBuild);
  h.parVvar = g_par_override >= 0 ? g_par_override == 1 : (threads && pairs >= c.parMin);
  h.st.lastPar = h.parVvar;
}

template <class T>
void wing_t(const void* key, int kCP0, int nSpan, int ldR1, int ldR2, int ldR3, const T* rNW, int ldE1, int ldE2,
            int ldE3, const T* eNW, int ldC1, int ldC2, const T* CP) {
  Handle<T>* hp = find<T>(key);
  if (!hp) fatal("ofk_cs_wing before ofk_cs_begin");
  Handle<T>& h = *hp;
  Geom<T>& g = h.g;
  const int D = g.iNWStart;
  // The pin reads r_NW(1:3, 1:nSpan+1, 1:iNWStart+1), Eps_NW(1, 1:nSpan, 1:iNWStart), CP(1:3, 1:nSpan)
  if (nSpan == 0) { h.wingsSeen++; return; }  // a wing without control points adds no target and no ring
  if (nSpan < 0 || kCP0 < 0 || kCP0 + nSpan > g.nT || ldR1 < 3 || ldR2 < nSpan + 1 || ldR3 < D + 1 || ldE1 < 1 ||
      ldE2 < nSpan || ldE3 < D || ldC1 < 3 || ldC2 < nSpan)
    fatal("ofk_cs_wing: extents smaller than the pin reads (r_NW, Eps_NW or CP)");
  T* cx = h.cx.data();
  T* cy = h.cy.data();
  T* cz = h.cz.data();
  for (int i = 0; i < nSpan; ++i) {
    const std::size_t o = std::size_t(i) * std::size_t(ldC1);
    cx[kCP0 + i] = CP[o];
    cy[kCP0 + i] = CP[o + 1];
    cz[kCP0 + i] = CP[o + 2];
  }
  auto R = [&](int c, int s, int d) {  // r_NW(c, s, d), 1-based
    return rNW[std::size_t(c - 1) + std::size_t(ldR1) * (std::size_t(s - 1) + std::size_t(ldR2) * std::size_t(d - 1))];
  };
  T* ring = h.ring.data();
  for (int s = 1; s <= nSpan; ++s) {
    for (int d = 1; d <= D; ++d) {
      T* P = ring + 13 * (std::size_t(kCP0 + s - 1) * std::size_t(D) + std::size_t(d - 1));
      for (int c = 1; c <= 3; ++c) {
        P[c - 1] = R(c, s, d);              // P1 = r_NW(:, iSpan  , iDepth  )  (FVW_Wings.f90:432)
        P[3 + c - 1] = R(c, s + 1, d);      // P2 = r_NW(:, iSpan+1, iDepth  )  (:434)
        P[6 + c - 1] = R(c, s + 1, d + 1);  // P3 = r_NW(:, iSpan+1, iDepth+1)  (:435)
        P[9 + c - 1] = R(c, s, d + 1);      // P4 = r_NW(:, iSpan  , iDepth+1)  (:433)
      }
      // RegParam = Eps_NW(1, iSpan, iDepth) for all four segments (:441)
      P[12] = eNW[std::size_t(ldE1) * (std::size_t(s - 1) + std::size_t(ldE2) * std::size_t(d - 1))];
    }
  }
  h.wingsSeen++;
}

template <class T>
void build_t(const void* key) {
  Handle<T>* hp = find<T>(key);
  if (!hp) fatal("ofk_cs_build before ofk_cs_begin");
  Handle<T>& h = *hp;
  Geom<T>& g = h.g;
  if (h.wingsSeen != h.nWings) fatal("ofk_cs_build: not every wing was given");
  T* cx = h.cx.data();
  T* cy = h.cy.data();
  T* cz = h.cz.data();
  for (int t = g.nT; t < g.nTpad; ++t) { cx[t] = cx[0]; cy[t] = cy[0]; cz[t] = cz[0]; }  // padded lanes
  const auto t0 = clk::now();
  if (h.mode == static_cast<int>(ofk::Mode::strict)) strict::build(g, h.parBuild);
  else if (h.fma) fastfma::build(g, h.parBuild);
  else fast::build(g, h.parBuild);
  h.st.tBuild += secs(t0, clk::now());
  h.st.nBuild++;
  h.built = true;
}

template <class T>
void vvar_t(const void* key, const T* gam, T* V) {
  Handle<T>* hp = find<T>(key);
  if (!hp || !hp->built) fatal("ofk_cs_vvar without a cache built for this solve");
  Handle<T>& h = *hp;
  const auto t0 = clk::now();
  if (h.mode == static_cast<int>(ofk::Mode::strict)) strict::vvar(h.g, gam, V, h.parVvar);
  else if (h.fma) fastfma::vvar(h.g, gam, V, h.parVvar);
  else fast::vvar(h.g, gam, V, h.parVvar);
  h.st.tVvar += secs(t0, clk::now());
  h.st.nVvar++;
}

template <class T>
void end_t(const void* key, int nIter, int maxIter) {
  if (Handle<T>* h = find<T>(key)) {
    h->st.nSolve++;
    h->st.nIter += nIter;
    if (nIter >= maxIter) h->st.nCap++;
    h->built = false;  // never trusted across calls
  }
}

}  // namespace
}  // namespace ofk_cs

#define OFK_CS_GUARD(name, ...)                                        \
  try {                                                                \
    __VA_ARGS__                                                        \
  } catch (const std::exception& e) {                                  \
    std::fprintf(stderr, "[ofk_cs] %s: exception: %s\n", name, e.what()); \
    std::abort();                                                      \
  } catch (...) {                                                      \
    std::fprintf(stderr, "[ofk_cs] %s: unknown exception\n", name);    \
    std::abort();                                                      \
  }

#define OFK_CS_RK(rk, call_d, call_f)               \
  if ((rk) == 8) { call_d; }                         \
  else if ((rk) == 4) { call_f; }                    \
  else { ofk_cs::fatal("unsupported real kind"); }

extern "C" {

int ofk_cs_enabled(void) { return ofk_cs::config().enabled ? 1 : 0; }

void ofk_cs_begin(const void* key, int rk, int nWings, int nCPtot, int iNWStart, int reg, const void* precision_ui,
                  const void* fourpi_inv, const void* min_exp_value, const void* lamb_oseen) {
  OFK_CS_GUARD("ofk_cs_begin", {
    if (rk == 8) {
      const double k[4] = {*static_cast<const double*>(precision_ui), *static_cast<const double*>(fourpi_inv),
                           *static_cast<const double*>(min_exp_value), *static_cast<const double*>(lamb_oseen)};
      ofk_cs::begin_t<double>(key, nWings, nCPtot, iNWStart, reg, k);
    } else if (rk == 4) {
      const float k[4] = {*static_cast<const float*>(precision_ui), *static_cast<const float*>(fourpi_inv),
                          *static_cast<const float*>(min_exp_value), *static_cast<const float*>(lamb_oseen)};
      ofk_cs::begin_t<float>(key, nWings, nCPtot, iNWStart, reg, k);
    } else {
      ofk_cs::fatal("unsupported real kind");
    }
  })
}

void ofk_cs_wing(const void* key, int rk, int kCP0, int nSpan, int ldR1, int ldR2, int ldR3, const void* rNW,
                 int ldE1, int ldE2, int ldE3, const void* eNW, int ldC1, int ldC2, const void* CP) {
  OFK_CS_GUARD("ofk_cs_wing", {
    OFK_CS_RK(rk,
              ofk_cs::wing_t<double>(key, kCP0, nSpan, ldR1, ldR2, ldR3, static_cast<const double*>(rNW), ldE1, ldE2,
                                     ldE3, static_cast<const double*>(eNW), ldC1, ldC2,
                                     static_cast<const double*>(CP)),
              ofk_cs::wing_t<float>(key, kCP0, nSpan, ldR1, ldR2, ldR3, static_cast<const float*>(rNW), ldE1, ldE2,
                                    ldE3, static_cast<const float*>(eNW), ldC1, ldC2,
                                    static_cast<const float*>(CP)))
  })
}

void ofk_cs_build(const void* key, int rk) {
  OFK_CS_GUARD("ofk_cs_build", { OFK_CS_RK(rk, ofk_cs::build_t<double>(key), ofk_cs::build_t<float>(key)) })
}

void ofk_cs_vvar(const void* key, int rk, const void* gam, void* V) {
  OFK_CS_GUARD("ofk_cs_vvar", {
    OFK_CS_RK(rk, ofk_cs::vvar_t<double>(key, static_cast<const double*>(gam), static_cast<double*>(V)),
              ofk_cs::vvar_t<float>(key, static_cast<const float*>(gam), static_cast<float*>(V)))
  })
}

void ofk_cs_end(const void* key, int nIter, int maxIter) {
  OFK_CS_GUARD("ofk_cs_end", {
    ofk_cs::end_t<double>(key, nIter, maxIter);
    ofk_cs::end_t<float>(key, nIter, maxIter);
    ofk_cs::log_iters(nIter, maxIter);
  })
}

void ofk_cs_free(const void* key) {
  OFK_CS_GUARD("ofk_cs_free", {
    auto& r = ofk_cs::registry();
    std::lock_guard<std::mutex> lock(r.mu);
    auto it = r.map.find(key);
    if (it != r.map.end()) {
      if (ofk_cs::config().profile) ofk_cs::report_entry(it->second, "free");
      r.map.erase(it);
    }
  })
}

void ofk_cs_set_variant(int mode, int fma, int par) {
  ofk_cs::g_mode_override = (mode == 0 || mode == 1) ? mode : -1;
  ofk_cs::g_fma_override = (fma == 0 || fma == 1) ? fma : -1;
  ofk_cs::g_par_override = (par == 0 || par == 1) ? par : -1;
}

void ofk_cs_counters(const void* key, long long* out) {
  for (int i = 0; i < 8; ++i) out[i] = 0;
  OFK_CS_GUARD("ofk_cs_counters", {
    auto add = [&](const auto* h) {
      if (!h) return;
      const ofk_cs::Stats& s = h->st;
      out[0] += s.nAlloc; out[1] += s.bytes; out[2] += s.nSolve; out[3] += s.nIter;
      out[4] += s.nBuild; out[5] += s.nVvar; out[7] = s.lastPar ? 1 : 0;
    };
    add(ofk_cs::find<double>(key));
    add(ofk_cs::find<float>(key));
    auto& r = ofk_cs::registry();
    std::lock_guard<std::mutex> lock(r.mu);
    out[6] = static_cast<long long>(r.map.size());
  })
}

}  // extern "C"
