// Unit induced-velocities: per-instance handles, lattice mirror, source sets and evaluation.
// C API and the design rationale: induced_velocities.hpp. Physics: none here (the bodies are
// unit biotsavart-direct's, the packing is induced_velocities_pack.cpp).
//
// Environment (read once, process-wide):
//   OFK_MODE            strict | fast (ofk_runtime.hpp)
//   OFK_IV_SOURCES      auto (default) | packed | lattice: where the device source set comes from
//                       (auto = lattice on device backends, packed on host backends)
//   OFK_IV_PROFILE=1    fence between phases and accumulate per-phase times; print them per
//                       handle at free / finalize (P3)
//   OFK_IV_MAX_HANDLES  live-handle cap (default 64), least recently used evicted
//   OFK_IV_GAP          dirty-span merge gap in elements (default 65536)
//   OFK_TC              kokkos (default) | fortran: VelocityMethod 2 (particle treecode) in C++/
//                       Kokkos (unit treecode) or in the original Fortran (A/B switch)
//   OFK_TC_DEVBUILD     auto (default) | 0 | 1: fast-mode tree build on the device
//                       (treecode_dbuild.hpp). auto = on device backends for n >= OFK_TC_DEVBUILD_MIN
//                       (default 4096); 1 = always (also on host backends: tests the algorithm on
//                       the CPU); 0 = never (the host build, uploaded). Strict always builds on the host.
//   OFK_TC_DEVBUILD_CHECK=1  after every device build, also run the host fast build and require the
//                       same topology (node and leaf arrays, radii, maxRegParam, leaf sources) bit for bit
// plus unit 1's OFK_BSD_* shape switches, which the evaluation uses unchanged.
//
// Unit treecode (VelocityMethod 2, idVelocityTreePart) shares this handle: its source set is the
// particle set (from the Fortran-packed Part), plus a flat tree built from it (treecode_impl.hpp).
#include "induced_velocities.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "biotsavart-direct/biotsavart_direct_common.hpp"
#include "biotsavart-direct/biotsavart_direct_consts.hpp"
#include "induced_velocities_common.hpp"
#include "common/ofk_runtime.hpp"
#include "treecode/treecode_common.hpp"

namespace ofk_iv {

using ofk_bsd::PartArgs;
using ofk_bsd::SegArgs;
using StageSpace = Kokkos::SharedHostPinnedSpace;  // pinned host memory on device backends

constexpr int kMethodSeg = 1;   // idVelocityBasic
constexpr int kMethodTree = 2;  // idVelocityTreePart (unit treecode)
constexpr int kMethodPart = 3;  // idVelocityPart

// ------------------------------------------------------------------------------------------
// Configuration
struct Config {
  int sources = 0;  // 0 auto, 1 packed, 2 lattice
  bool profile = false;
  std::size_t maxHandles = 64;
  long gap = 65536;
  bool treeKokkos = true;  // OFK_TC
  int devBuild = -1;       // OFK_TC_DEVBUILD: -1 auto, 0 never, 1 always
  long devBuildMin = 4096;  // OFK_TC_DEVBUILD_MIN
  bool devBuildCheck = false;
};
const Config& config() {
  static const Config c = [] {
    Config r;
    if (const char* s = std::getenv("OFK_IV_SOURCES")) {
      if (!std::strcmp(s, "packed")) r.sources = 1;
      else if (!std::strcmp(s, "lattice")) r.sources = 2;
    }
    if (const char* s = std::getenv("OFK_IV_PROFILE")) r.profile = (*s == '1');
    if (const char* s = std::getenv("OFK_IV_MAX_HANDLES")) r.maxHandles = std::max(1, std::atoi(s));
    if (const char* s = std::getenv("OFK_IV_GAP")) r.gap = std::max(0L, std::atol(s));
    if (const char* s = std::getenv("OFK_TC")) r.treeKokkos = std::strcmp(s, "fortran") != 0;
    if (const char* s = std::getenv("OFK_TC_DEVBUILD")) r.devBuild = (*s == '0') ? 0 : (*s == '1') ? 1 : -1;
    if (const char* s = std::getenv("OFK_TC_DEVBUILD_MIN")) r.devBuildMin = std::max(2L, std::atol(s));
    if (const char* s = std::getenv("OFK_TC_DEVBUILD_CHECK")) r.devBuildCheck = (*s == '1');
    return r;
  }();
  return c;
}
std::atomic<int> g_sources_override{-1};
bool use_lattice() {
  int m = g_sources_override.load();
  if (m < 0) m = config().sources;
  return m == 2 || (m == 0 && !kOnHost);
}

// ------------------------------------------------------------------------------------------
// Statistics (always counted; times only with OFK_IV_PROFILE=1, which fences between phases)
struct Stats {
  long long nAlloc = 0, bytesAlloc = 0, nSrc = 0, nEval = 0, h2d = 0, d2h = 0, latCmp = 0, latUp = 0;
  long long latChanged = 0;  // lattice bytes that differed from the mirror (= what an upload-once scheme must move)
  double tSync = 0, tUpload = 0, tPack = 0, tCpH2D = 0, tKernel = 0, tD2H = 0, tHostAdd = 0;
  long long nTree = 0, treeNodes = 0, treeLeaf = 0, treeDepth = 0;  // last tree built
  long long nDevBuild = 0, nDevCheck = 0;
  double tBuild = 0, tTreeH2D = 0;
};

struct Timer {
  bool on;
  Kokkos::Timer t;
  explicit Timer(bool o) : on(o) {}
  void lap(double& acc) {
    if (!on) return;
    Kokkos::fence("ofk_iv_profile");
    acc += t.seconds();
    t.reset();
  }
};

struct Src {
  bool valid = false;
  int method = 0, reg = 0, regPart = 0, nSeg = 0, nPart = 0;
};

// Lattice layout of one wing in the flat mirror: 6 arrays (r, Gamma, Eps) x (NW, FW).
struct WingLay {
  long off = -1;
  long n[6] = {0, 0, 0, 0, 0, 0};       // allocated elements per array
  long aoff[6] = {0, 0, 0, 0, 0, 0};    // element offset of each array
  int dims[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  long synced[6] = {0, 0, 0, 0, 0, 0};  // leading depth slabs known identical in shadow and mirror
  long total() const { return n[0] + n[1] + n[2] + n[3] + n[4] + n[5]; }
};

template <class T>
struct Handle {
  using Stage1 = Kokkos::View<T*, StageSpace>;
  using Stage2 = Kokkos::View<T**, Kokkos::LayoutLeft, StageSpace>;
  const void* key = nullptr;
  Stats st;
  // lattice mirror (device) and its host shadow (pinned; device backends only)
  std::vector<WingLay> wings;
  DV1<T> lat;
  Stage1 shadow;
  std::vector<std::pair<long, long>> dirty;  // [lo, hi) element spans to upload
  // block table
  Kokkos::View<Block*, MemSpace> dblk;
  Kokkos::View<Block*, StageSpace> hblk;
  // sources (SoA), packed-mode host staging, targets/outputs, skinny scratch
  DV2<T> seg, part;
  Stage2 hseg, hpart;
  DV1<T> cp, u, scratch;
  Stage1 hcp, hu;
  Src src;
  long capSeg = 0, capPart = 0;  // capacity hints (ofk_iv_reserve)
  // treecode (method 2): host build buffers, device copies (device backends), validity
  ofk_tc::HostTree<T> ht;
  ofk_tc::NodeV<T> dnr, dmom;
  ofk_tc::NodeIV dni;
  DV2<T> dleaf;
  bool treeValid = false;
  int treeMode = -1;
  T treeThr = T(0);
  // fast mode, device build (treecode_dbuild.hpp): build buffers; the tree lands in dnr/dmom/dni/dleaf
  ofk_tc::DevBuild<T> db;
  bool treeDev = false;     // the current tree was built by build_device
  int tNodes = 0, tLeaf = 0;  // the current tree's node and leaf counts
};

// ------------------------------------------------------------------------------------------
// Buffers: allocate on first need, grow by doubling, never shrink (P6).
template <class V>
void ensure1(V& v, std::size_t n, const char* name, Stats& st, bool keep = false) {
  if (n == 0 || v.extent(0) >= n) return;
  const std::size_t cap = std::max<std::size_t>(n, 2 * v.extent(0));
  if (keep && v.extent(0) > 0) {
    Kokkos::resize(Kokkos::WithoutInitializing, v, cap);
  } else {
    v = V();
    v = V(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), cap);
  }
  st.nAlloc++;
  st.bytesAlloc += (long long)(cap * sizeof(typename V::value_type));
}
template <class V>
void ensure2(V& v, std::size_t n, std::size_t cols, const char* name, Stats& st) {
  if (n == 0 || (v.extent(0) >= n && v.extent(1) == cols)) return;
  const std::size_t cap = std::max<std::size_t>(n, 2 * v.extent(0));
  v = V();
  v = V(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), cap, cols);
  st.nAlloc++;
  st.bytesAlloc += (long long)(cap * cols * sizeof(typename V::value_type));
}

// ------------------------------------------------------------------------------------------
// Handle registry: process-wide map key -> handles (one per real kind). Never destroyed (no
// static-destruction order problem with Kokkos::finalize); emptied by a finalize hook.
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
  long long evicted = 0;
};
Registry& registry() {
  static Registry* r = new Registry;
  return *r;
}

template <class T>
void report_one(const Handle<T>& h, const char* why) {
  const Stats& s = h.st;
  std::fprintf(stderr,
               "[ofk_iv] %s key=%p r%d: sources=%lld evals=%lld allocs=%lld (%.1f MB) "
               "H2D=%.1f MB D2H=%.1f MB lattice compared=%.1f MB changed=%.1f MB uploaded=%.1f MB | "
               "t[s] sync=%.4f upload=%.4f pack=%.4f cpH2D=%.4f kernel=%.4f D2H=%.4f hostadd=%.4f | "
               "tree: builds=%lld (device %lld) last nodes=%lld leaf=%lld depth=%lld t[s] build=%.4f treeH2D=%.4f\n",
               why, h.key, int(sizeof(T)), s.nSrc, s.nEval, s.nAlloc, s.bytesAlloc / 1e6,
               s.h2d / 1e6, s.d2h / 1e6, s.latCmp / 1e6, s.latChanged / 1e6, s.latUp / 1e6, s.tSync, s.tUpload,
               s.tPack, s.tCpH2D, s.tKernel, s.tD2H, s.tHostAdd, s.nTree, s.nDevBuild, s.treeNodes, s.treeLeaf,
               s.treeDepth, s.tBuild, s.tTreeH2D);
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

// The handle of `key` for real kind T, created lazily. The caller owns the instance, so the
// handle is used by one thread at a time; the map itself is guarded.
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
    if (r.map.size() >= config().maxHandles) {  // evict the least recently used
      auto lru = r.map.begin();
      for (auto i = r.map.begin(); i != r.map.end(); ++i)
        if (i->second.lastUse < lru->second.lastUse) lru = i;
      if (config().profile) report_entry(lru->second, "evicted");
      if (r.evicted++ == 0)
        std::fprintf(stderr, "[ofk_iv] note: more than %zu live handles; evicting the least recently used\n",
                     config().maxHandles);
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
  if (it == r.map.end()) return nullptr;
  return slot<T>(it->second).get();
}

// ------------------------------------------------------------------------------------------
// Lattice mirror: slab-wise validation against the host shadow (device backends) or the mirror
// itself (host backends), bitwise (memcmp), so any change by anyone (glue copies, rollback,
// Map_LL_NW, FakeGroundEffect, a different x) is seen. Unchanged slabs are not uploaded.
template <class T>
void mark_dirty(Handle<T>& h, long lo, long hi) {
  if (kOnHost) return;
  if (!h.dirty.empty() && lo - h.dirty.back().second <= config().gap) {
    h.dirty.back().second = std::max(h.dirty.back().second, hi);
  } else {
    h.dirty.emplace_back(lo, hi);
  }
}

template <class T>
void sync_array(Handle<T>& h, WingLay& w, int a, const T* src, long slab, long nAct) {
  if (nAct <= 0 || slab <= 0) return;
  T* sh = (kOnHost ? h.lat.data() : h.shadow.data()) + w.aoff[a];
  const std::size_t bytes = std::size_t(slab) * sizeof(T);
  for (long d = 0; d < nAct; ++d) {
    const T* s = src + d * slab;
    T* t = sh + d * slab;
    bool same = false;
    if (d < w.synced[a]) {
      same = (std::memcmp(s, t, bytes) == 0);
      h.st.latCmp += (long long)bytes;
    }
    if (!same) {
      std::memcpy(t, s, bytes);
      h.st.latChanged += (long long)bytes;
      mark_dirty(h, w.aoff[a] + d * slab, w.aoff[a] + (d + 1) * slab);
    }
  }
  w.synced[a] = std::max(w.synced[a], nAct);
}

template <class T>
void lattice(Handle<T>& h, int iW, int nWings, int nNWact, int nSpanNWP, int nDepNWP,
             int nSpanNWG, int nDepNWG, const T* rNW, const T* gNW, const T* eNW, int nFWact,
             int nSpanFWP, int nDepFWP, int nSpanFWG, int nDepFWG, const T* rFW, const T* gFW,
             const T* eFW) {
  if (iW < 1 || iW > nWings) throw std::invalid_argument("ofk_iv_lattice: iW out of range");
  if (nNWact + 1 > nDepNWP || nNWact > nDepNWG || nFWact + 1 > nDepFWP || nFWact > nDepFWG)
    throw std::invalid_argument("ofk_iv_lattice: active depth exceeds the allocated lattice");
  Timer tm(config().profile);
  if (int(h.wings.size()) != nWings) {
    h.wings.assign(nWings, WingLay{});
    h.dirty.clear();
  }
  WingLay& w = h.wings[iW - 1];
  const int dims[10] = {nSpanNWP, nDepNWP, nSpanNWG, nDepNWG, nSpanFWP,
                        nDepFWP,  nSpanFWG, nDepFWG, 0,        0};
  if (iW > 1 && h.wings[iW - 2].off < 0)
    throw std::invalid_argument("ofk_iv_lattice: wings must be passed in order 1..nWings");
  const long off = (iW == 1) ? 0 : h.wings[iW - 2].off + h.wings[iW - 2].total();
  bool same = (w.off == off);
  for (int i = 0; i < 10; ++i) same = same && (w.dims[i] == dims[i]);
  if (!same) {  // new or changed layout: nothing of this wing is known in sync
    w = WingLay{};
    w.off = off;
    for (int i = 0; i < 10; ++i) w.dims[i] = dims[i];
    w.n[0] = 3L * nSpanNWP * nDepNWP;
    w.n[1] = 1L * nSpanNWG * nDepNWG;
    w.n[2] = 3L * nSpanNWG * nDepNWG;
    w.n[3] = 3L * nSpanFWP * nDepFWP;
    w.n[4] = 1L * nSpanFWG * nDepFWG;
    w.n[5] = 3L * nSpanFWG * nDepFWG;
    long o = off;
    for (int a = 0; a < 6; ++a) { w.aoff[a] = o; o += w.n[a]; }
    for (int k = iW; k < nWings; ++k) h.wings[k] = WingLay{};  // later wings move
  }
  const std::size_t need = std::size_t(off + w.total());
  ensure1(h.lat, need, "ofk_iv_lattice", h.st, /*keep=*/true);
  if (!kOnHost) ensure1(h.shadow, need, "ofk_iv_lattice_shadow", h.st, /*keep=*/true);
  sync_array(h, w, 0, rNW, 3L * nSpanNWP, nNWact + 1);
  sync_array(h, w, 1, gNW, 1L * nSpanNWG, nNWact);
  sync_array(h, w, 2, eNW, 3L * nSpanNWG, nNWact);
  sync_array(h, w, 3, rFW, 3L * nSpanFWP, nFWact + 1);
  sync_array(h, w, 4, gFW, 1L * nSpanFWG, nFWact);
  sync_array(h, w, 5, eFW, 3L * nSpanFWG, nFWact);
  tm.lap(h.st.tSync);
}

template <class T>
void upload_dirty(Handle<T>& h) {
  if (kOnHost) { h.dirty.clear(); return; }
  for (const auto& sp : h.dirty) {
    const auto r = Kokkos::make_pair(std::size_t(sp.first), std::size_t(sp.second));
    Kokkos::deep_copy(Kokkos::subview(h.lat, r), Kokkos::subview(h.shadow, r));
    const long long b = (long long)(sp.second - sp.first) * (long long)sizeof(T);
    h.st.latUp += b;
    h.st.h2d += b;
  }
  h.dirty.clear();
}

// ------------------------------------------------------------------------------------------
// Capacity from the Fortran side (P6). The misc arrays are allocated once at their maximum
// (FVW_InitMiscVarsPostParam: Sgmt for nNWMax/nFWMax, Part, CPs), so the first call can size
// every buffer for the whole simulation; growth by doubling remains for anything larger.
template <class T>
void reserve(Handle<T>& h, long latTot, long capSeg, long capPart, int nBlocks) {
  ensure1(h.lat, std::size_t(std::max(0L, latTot)), "ofk_iv_lattice", h.st, /*keep=*/true);
  if (!kOnHost) ensure1(h.shadow, std::size_t(std::max(0L, latTot)), "ofk_iv_lattice_shadow", h.st, /*keep=*/true);
  ensure1(h.dblk, std::size_t(std::max(1, nBlocks)), "ofk_iv_blocks", h.st);
  if (!kOnHost) ensure1(h.hblk, std::size_t(std::max(1, nBlocks)), "ofk_iv_blocks_stage", h.st);
  h.capSeg = std::max(h.capSeg, capSeg);
  h.capPart = std::max(h.capPart, capPart);
}

// ------------------------------------------------------------------------------------------
// Source sets
template <class T>
void sources_lattice(Handle<T>& h, int method, int reg, int regPart, int partPerSeg,
                     int iDepthStart, int nNW, int nFW, bool mirror, bool fwShed,
                     bool lastNWShed, int nSegF, int nPartF) {
  if (method != kMethodSeg && method != kMethodPart)
    throw std::invalid_argument("ofk_iv_sources_lattice: method must be 1 or 3");
  Timer tm(config().profile);
  h.src.valid = false;
  const int nW = int(h.wings.size());
  for (const auto& w : h.wings)
    if (w.off < 0) throw std::invalid_argument("ofk_iv_sources_lattice: lattice of some wing not passed");
  // PackPanelsToSegments (FVW_Subs.f90:913): NW blocks if nNW-iDepthStart >= 0, FW blocks if nFW > 0
  std::vector<Block> blk;
  int panels = 0, segs = 0;
  auto add = [&](const WingLay& w, bool fw, int nAct, int iDS, bool shed, bool shedLast, bool hack) {
    Block b;
    b.nSpan = fw ? w.dims[4] : w.dims[0];
    b.nDepth = nAct + 1;
    b.iDepthStart = iDS;
    b.shed = shed; b.shedLast = shedLast; b.hackEps = hack;
    b.rOff = w.aoff[fw ? 3 : 0]; b.gOff = w.aoff[fw ? 4 : 1]; b.eOff = w.aoff[fw ? 5 : 2];
    b.rLd2 = fw ? w.dims[4] : w.dims[0];
    b.gLd1 = fw ? w.dims[6] : w.dims[2];
    b.eLd2 = fw ? w.dims[6] : w.dims[2];
    if (b.gLd1 < b.nSpan - 1) throw std::invalid_argument("ofk_iv: Gamma has fewer rows than span panels");
    const int rows = b.nDepth - b.iDepthStart;
    b.nPanels = (rows > 0 && b.nSpan > 1) ? (b.nSpan - 1) * rows : 0;
    b.panelBase = panels;
    b.segBase = segs;
    b.nSeg = block_segments(b);
    panels += b.nPanels;
    segs += b.nSeg;
    blk.push_back(b);
  };
  if (nNW - iDepthStart >= 0)
    for (int iW = 0; iW < nW; ++iW) add(h.wings[iW], false, nNW, iDepthStart, true, lastNWShed, false);
  if (nFW > 0)
    for (int iW = 0; iW < nW; ++iW) add(h.wings[iW], true, nFW, 1, fwShed, fwShed, true);
  const int nSegHalf = segs;
  const int nSeg = mirror ? 2 * segs : segs;
  if (nSeg != nSegF)
    throw std::runtime_error("ofk_iv: closed-form segment count " + std::to_string(nSeg) +
                             " differs from the Fortran pack's " + std::to_string(nSegF));
  // drop empty blocks from the panel search (they have no panels)
  std::vector<Block> act;
  for (const auto& b : blk)
    if (b.nPanels > 0) act.push_back(b);
  const int nBlocks = int(act.size());
  ensure1(h.dblk, std::max(1, nBlocks), "ofk_iv_blocks", h.st);
  if constexpr (kOnHost) {
    for (int i = 0; i < nBlocks; ++i) h.dblk(i) = act[i];
  } else {
    ensure1(h.hblk, std::max(1, nBlocks), "ofk_iv_blocks_stage", h.st);
    for (int i = 0; i < nBlocks; ++i) h.hblk(i) = act[i];
    const auto r = Kokkos::make_pair(std::size_t(0), std::size_t(nBlocks));
    if (nBlocks > 0) Kokkos::deep_copy(Kokkos::subview(h.dblk, r), Kokkos::subview(h.hblk, r));
    h.st.h2d += (long long)(nBlocks * sizeof(Block));
  }
  upload_dirty(h);
  tm.lap(h.st.tUpload);
  ensure2(h.seg, std::size_t(std::max<long>(nSeg, h.capSeg)), ofk_bsd::kSegSoaCols, "ofk_iv_seg_soa", h.st);
  if (nSeg > 0)
    pack::lattice_to_segments<T>(h.lat, BlockView(h.dblk), nBlocks, panels, nSegHalf, mirror, h.seg);
  int nPart = 0;
  if (method == kMethodPart) {
    nPart = nSeg * partPerSeg;
    if (nPart != nPartF)
      throw std::runtime_error("ofk_iv: particle count " + std::to_string(nPart) +
                               " differs from the Fortran's " + std::to_string(nPartF));
    ensure2(h.part, std::size_t(std::max<long>(nPart, h.capPart)), ofk_bsd::kPartSoaCols, "ofk_iv_part_soa", h.st);
    if (nSeg > 0) pack::segments_to_particles<T>(h.seg, nSeg, partPerSeg, h.part);
  }
  tm.lap(h.st.tPack);
  h.src = Src{true, method, reg, regPart, nSeg, nPart};
  h.st.nSrc++;
}

// Upload rows 0..n-1 of the staged SoA, column by column (each column is contiguous).
template <class T, class D, class H>
void upload_rows(Handle<T>& h, const D& d, const H& s, int n, int cols) {
  const auto r = Kokkos::make_pair(std::size_t(0), std::size_t(n));
  for (int c = 0; c < cols; ++c)
    Kokkos::deep_copy(Kokkos::subview(d, r, c), Kokkos::subview(s, r, c));
  h.st.h2d += (long long)n * cols * (long long)sizeof(T);
}

template <class T>
void sources_packed(Handle<T>& h, int method, int reg, int nSeg, int ldPts, int nPts,
                    const T* Pts, int ldC, const int* C, const T* Gam, const T* Eps, int regPart,
                    int nPart, int ldP, const T* PP, int ldA, const T* PA, const T* PR) {
  if (method != kMethodSeg && method != kMethodPart && method != kMethodTree)
    throw std::invalid_argument("ofk_iv_sources_packed: method must be 1, 2 or 3");
  Timer tm(config().profile);
  h.src.valid = false;
  h.treeValid = false;
  using HostExec = Kokkos::DefaultHostExecutionSpace;
  if (nSeg < 0) nSeg = 0;
  if (nPart < 0) nPart = 0;
  if (method == kMethodSeg) {
    if (nSeg > 0 && (ldPts < 3 || ldC < 2)) throw std::invalid_argument("ofk_iv: bad extents");
    for (int j = 0; j < nSeg; ++j)
      for (int q = 0; q < 2; ++q) {
        const int v = C[std::size_t(j) * ldC + q];
        if (v < 1 || v > nPts)
          throw std::out_of_range("ofk_iv: SegConnct index " + std::to_string(v) + " outside 1.." +
                                  std::to_string(nPts));
      }
    ensure2(h.seg, std::size_t(std::max<long>(nSeg, h.capSeg)), ofk_bsd::kSegSoaCols, "ofk_iv_seg_soa", h.st);
    if (!kOnHost) ensure2(h.hseg, std::size_t(std::max<long>(nSeg, h.capSeg)), ofk_bsd::kSegSoaCols, "ofk_iv_seg_stage", h.st);
    T* S = kOnHost ? h.seg.data() : h.hseg.data();
    const std::size_t ld = kOnHost ? h.seg.extent(0) : h.hseg.extent(0);
    if (nSeg > 0)
      Kokkos::parallel_for(
          "ofk_iv_pack_soa_seg", Kokkos::RangePolicy<HostExec>(0, nSeg), [=](const int j) {
            const std::size_t p1 = std::size_t(C[std::size_t(j) * ldC] - 1) * ldPts;
            const std::size_t p2 = std::size_t(C[std::size_t(j) * ldC + 1] - 1) * ldPts;
            S[j + 0 * ld] = Pts[p1]; S[j + 1 * ld] = Pts[p1 + 1]; S[j + 2 * ld] = Pts[p1 + 2];
            S[j + 3 * ld] = Pts[p2]; S[j + 4 * ld] = Pts[p2 + 1]; S[j + 5 * ld] = Pts[p2 + 2];
            S[j + 6 * ld] = Gam[j];
            S[j + 7 * ld] = reg == 0 ? T(0) : Eps[j];  // RegParam is not read for RegNone
          });
    HostExec().fence();
    if (!kOnHost && nSeg > 0) upload_rows(h, h.seg, h.hseg, nSeg, ofk_bsd::kSegSoaCols);
  } else {
    if (nPart > 0 && (ldP < 3 || ldA < 3)) throw std::invalid_argument("ofk_iv: bad extents");
    // Tree (method 2): the particles are only read by the host build (the device gets the tree
    // and its packed leaves), so on device backends only the host staging is filled. The build
    // reads the REAL RegParam even for RegNone (maxRegParam, :738), so it is always copied.
    const bool tree = method == kMethodTree;
    if (kOnHost || !tree)
      ensure2(h.part, std::size_t(std::max<long>(nPart, h.capPart)), ofk_bsd::kPartSoaCols, "ofk_iv_part_soa", h.st);
    if (!kOnHost) ensure2(h.hpart, std::size_t(std::max<long>(nPart, h.capPart)), ofk_bsd::kPartSoaCols, "ofk_iv_part_stage", h.st);
    T* S = kOnHost ? h.part.data() : h.hpart.data();
    const std::size_t ld = kOnHost ? h.part.extent(0) : h.hpart.extent(0);
    if (nPart > 0)
      Kokkos::parallel_for(
          "ofk_iv_pack_soa_part", Kokkos::RangePolicy<HostExec>(0, nPart), [=](const int j) {
            const std::size_t p = std::size_t(j) * ldP, a = std::size_t(j) * ldA;
            S[j + 0 * ld] = PP[p]; S[j + 1 * ld] = PP[p + 1]; S[j + 2 * ld] = PP[p + 2];
            S[j + 3 * ld] = PA[a]; S[j + 4 * ld] = PA[a + 1]; S[j + 5 * ld] = PA[a + 2];
            S[j + 6 * ld] = (regPart == 0 && !tree) ? T(0) : PR[j];
          });
    HostExec().fence();
    if (!kOnHost && !tree && nPart > 0) upload_rows(h, h.part, h.hpart, nPart, ofk_bsd::kPartSoaCols);
  }
  tm.lap(h.st.tPack);
  h.src = Src{true, method, reg, regPart, method == kMethodSeg ? nSeg : 0,
              method != kMethodSeg ? nPart : 0};
  h.st.nSrc++;
}

// ------------------------------------------------------------------------------------------
// Treecode (method 2): build the flat tree from the particle source set (host build, in the
// current mode), then make it visible to the kernels: on host backends the evaluation reads the
// host build buffers directly (unmanaged views); on device backends they are copied into the
// handle's persistent device buffers (sized from the capacity hint on the first call, P6).
template <class V>
void ensureR(V& v, std::size_t n, std::size_t cols, const char* name, Stats& st) {
  if (n == 0 || (v.extent(0) >= n && v.extent(1) == cols)) return;
  const std::size_t cap = std::max<std::size_t>(n, 2 * v.extent(0));
  v = V();
  v = V(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), cap, cols);
  st.nAlloc++;
  st.bytesAlloc += (long long)(cap * cols * sizeof(typename V::value_type));
}

// Fast mode, device build: is it used for this n? (strict always builds on the host)
inline bool use_devbuild(ofk::Mode mode, int n) {
  if (mode != ofk::Mode::fast || n < 2) return false;
  const Config& c = config();
  if (c.devBuild == 0) return false;
  if (c.devBuild == 1) return true;
  return !kOnHost && long(n) >= c.devBuildMin;
}

// OFK_TC_DEVBUILD_CHECK: the device-built topology must equal the host fast build's bit for bit
// (ni, radius, maxRegParam, leaf sources); centers and moments are only compared for information.
template <class T>
void devbuild_check(Handle<T>& h, const ofk_tc::PartCols<T>& pc, const ofk_tc::TreeConsts<T>& tc) {
  ofk_tc::fast::build<T>(h.ht, pc, tc);
  const auto& t = h.ht;
  auto fail = [&](const std::string& w) {
    throw std::runtime_error("OFK_TC_DEVBUILD_CHECK: device build differs from the host build: " + w);
  };
  if (t.nNodes != h.tNodes || t.nLeaf != h.tLeaf)
    fail("nodes " + std::to_string(h.tNodes) + " vs " + std::to_string(t.nNodes) + ", leaves " +
         std::to_string(h.tLeaf) + " vs " + std::to_string(t.nLeaf));
  const std::size_t nn = std::size_t(t.nNodes);
  if (nn > 0) {
    auto ni = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(h.dni, Kokkos::make_pair(std::size_t(0), nn), Kokkos::ALL));
    auto nr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(h.dnr, Kokkos::make_pair(std::size_t(0), nn), Kokkos::ALL));
    double dc = 0;
    for (std::size_t i = 0; i < nn; ++i) {
      for (int q = 0; q < ofk_tc::kNI; ++q)
        if (ni(i, q) != t.ni[i * ofk_tc::kNI + q]) fail("ni(" + std::to_string(i) + "," + std::to_string(q) + ")");
      for (int q = 3; q < ofk_tc::kNR; ++q)
        if (nr(i, q) != t.nr[i * ofk_tc::kNR + q]) fail("radius/maxReg of node " + std::to_string(i));
      for (int q = 0; q < 3; ++q) {
        const double r = std::abs(double(nr(i, q)) - double(t.nr[i * ofk_tc::kNR + q])) / std::max(1.0, std::abs(double(nr(i, 3))));
        dc = std::max(dc, r);
      }
    }
    if (dc > 1e-10) fail("barycenter differs by " + std::to_string(dc) + " of the radius");
  }
  if (t.nLeaf > 0) {
    auto lf = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), h.dleaf);
    for (int c = 0; c < 7; ++c)
      for (int r = 0; r < t.nLeaf; ++r) {
        const T a = lf(r, c), b = t.leaf[std::size_t(r) + std::size_t(c) * t.ld];
        if (!(a == b || (a != a && b != b))) fail("leaf row " + std::to_string(r) + " col " + std::to_string(c));
      }
  }
  h.st.nDevCheck++;
}

template <class T>
void tree_build(Handle<T>& h, T thr) {
  if (!h.src.valid || h.src.method != kMethodTree)
    throw std::logic_error("ofk_iv_tree_build: no particle source set (method 2) for this instance");
  Timer tm(config().profile);
  const int n = h.src.nPart;
  ofk_tc::PartCols<T> pc;
  pc.n = n;
  {
    const T* base = kOnHost ? h.part.data() : h.hpart.data();
    const std::size_t ld = kOnHost ? h.part.extent(0) : h.hpart.extent(0);
    for (int c = 0; c < 7; ++c) pc.c[c] = base + c * ld;
  }
  ofk_tc::TreeConsts<T> tc;
  tc.degenerate = thr;
  tc.rootFactor = ofk_bsd::lit<T>(1.001, 1.001f);
  tc.eqTol = T(100) * std::numeric_limits<T>::epsilon() / T(2);
  tc.ncrit = 32;
  const ofk::Mode mode = ofk::mode();
  // capacity: nodes ~ n/4 at most in practice (leaf buckets up to 32), leaves <= n
  auto ensure_tree = [&](long nNodes, long nLeaf) {
    const std::size_t capN = std::size_t(std::max<long>({nNodes, h.capPart / 4, 64L}));
    const std::size_t capL = std::size_t(std::max<long>({nLeaf, h.capPart, long(n), 1L}));
    ensureR(h.dnr, capN, ofk_tc::kNR, "ofk_tc_nr", h.st);
    ensureR(h.dmom, capN, ofk_tc::kNM, "ofk_tc_mom", h.st);
    ensureR(h.dni, capN, ofk_tc::kNI, "ofk_tc_ni", h.st);
    ensure2(h.dleaf, capL, 7, "ofk_tc_leaf", h.st);
  };
  if (use_devbuild(mode, n)) {
    // ---- fast mode, device build: the particles go up once, the tree is built where it is read
    if constexpr (!kOnHost) {
      ensure2(h.part, std::size_t(std::max<long>(n, h.capPart)), ofk_bsd::kPartSoaCols, "ofk_iv_part_soa", h.st);
      upload_rows(h, h.part, h.hpart, n, ofk_bsd::kPartSoaCols);
      tm.lap(h.st.tTreeH2D);
    }
    const long long a0 = h.db.nAlloc, b0 = h.db.bytes;
    ofk_tc::fast::build_device<T>(h.db, h.part, n, tc, h.capPart);
    h.st.nAlloc += h.db.nAlloc - a0;
    h.st.bytesAlloc += h.db.bytes - b0;
    ensure_tree(h.db.nNodes, h.db.nLeaf);
    ofk_tc::FlatTree<T> out;
    out.nr = h.dnr; out.mom = h.dmom; out.ni = h.dni; out.leaf = h.dleaf;
    ofk_tc::fast::write_device<T>(h.db, h.part, out);
    tm.lap(h.st.tBuild);
    h.treeDev = true;
    h.tNodes = h.db.nNodes;
    h.tLeaf = h.db.nLeaf;
    h.st.nDevBuild++;
    h.st.treeDepth = h.db.depth;
    if (config().devBuildCheck) {
      ExecSpace().fence();
      devbuild_check(h, pc, tc);
    }
  } else {
    if (mode == ofk::Mode::strict) ofk_tc::strict::build<T>(h.ht, pc, tc);
    else ofk_tc::fast::build<T>(h.ht, pc, tc);
    tm.lap(h.st.tBuild);
    h.treeDev = false;
    h.tNodes = h.ht.nNodes;
    h.tLeaf = h.ht.nLeaf;
    h.st.treeDepth = h.ht.depth;
    if constexpr (!kOnHost) {
      const auto& t = h.ht;
      ensure_tree(t.nNodes, t.nLeaf);
      using H1 = Kokkos::View<const T*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
      using H1i = Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
      using D1 = Kokkos::View<T*, MemSpace, Kokkos::MemoryUnmanaged>;
      using D1i = Kokkos::View<int*, MemSpace, Kokkos::MemoryUnmanaged>;
      long long b = 0;
      if (t.nNodes > 0) {
        const std::size_t nn = std::size_t(t.nNodes);
        Kokkos::deep_copy(D1(h.dnr.data(), nn * ofk_tc::kNR), H1(t.nr.data(), nn * ofk_tc::kNR));
        Kokkos::deep_copy(D1(h.dmom.data(), nn * ofk_tc::kNM), H1(t.mom.data(), nn * ofk_tc::kNM));
        Kokkos::deep_copy(D1i(h.dni.data(), nn * ofk_tc::kNI), H1i(t.ni.data(), nn * ofk_tc::kNI));
        b += (long long)(nn * (ofk_tc::kNR + ofk_tc::kNM) * sizeof(T) + nn * ofk_tc::kNI * sizeof(int));
      }
      if (t.nLeaf > 0) {
        const std::size_t ldd = h.dleaf.extent(0);
        for (int c = 0; c < 7; ++c)
          Kokkos::deep_copy(D1(h.dleaf.data() + c * ldd, std::size_t(t.nLeaf)),
                            H1(t.leaf.data() + std::size_t(c) * t.ld, std::size_t(t.nLeaf)));
        b += (long long)(std::size_t(t.nLeaf) * 7 * sizeof(T));
      }
      h.st.h2d += b;
      tm.lap(h.st.tTreeH2D);
    }
  }
  h.treeValid = true;
  h.treeMode = int(mode);
  h.treeThr = thr;
  h.st.nTree++;
  h.st.treeNodes = h.tNodes;
  h.st.treeLeaf = h.tLeaf;
}

// The flat tree as the kernels see it.
template <class T>
ofk_tc::FlatTree<T> flat_tree(Handle<T>& h) {
  ofk_tc::FlatTree<T> f;
  const auto& t = h.ht;
  f.nNodes = h.tNodes;
  f.nLeaf = h.tLeaf;
  if (kOnHost && !h.treeDev) {
    using U = Kokkos::MemoryUnmanaged;
    f.nr = ofk_tc::NodeV<T>(Kokkos::View<T**, Kokkos::LayoutRight, MemSpace, U>(
        const_cast<T*>(t.nr.data()), std::size_t(t.nNodes), ofk_tc::kNR));
    f.mom = ofk_tc::NodeV<T>(Kokkos::View<T**, Kokkos::LayoutRight, MemSpace, U>(
        const_cast<T*>(t.mom.data()), std::size_t(t.nNodes), ofk_tc::kNM));
    f.ni = ofk_tc::NodeIV(Kokkos::View<int**, Kokkos::LayoutRight, MemSpace, U>(
        const_cast<int*>(t.ni.data()), std::size_t(t.nNodes), ofk_tc::kNI));
    f.leaf = DV2<T>(Kokkos::View<T**, Kokkos::LayoutLeft, MemSpace, U>(
        const_cast<T*>(t.leaf.data()), std::size_t(std::max(t.ld, 1)), 7));
  } else {
    f.nr = h.dnr; f.mom = h.dmom; f.ni = h.dni; f.leaf = h.dleaf;
  }
  return f;
}

// ------------------------------------------------------------------------------------------
// Evaluation
template <class T>
void eval(Handle<T>& h, int nCPs, int capCP, int ldCP, const T* CPs, int ldU, T* Uind, T minnorm,
          T bf = T(0), T dd = T(0)) {
  if (!h.src.valid) throw std::logic_error("ofk_iv_eval: no valid source set for this instance");
  if (h.src.method == kMethodTree) {
    if (!h.treeValid) throw std::logic_error("ofk_iv_eval: tree not built (ofk_iv_tree_build)");
    if (h.treeMode != int(ofk::mode())) tree_build(h, h.treeThr);  // mode switched since the build
  }
  if (nCPs <= 0) return;
  if (ldCP < 3 || ldU < 3) throw std::invalid_argument("ofk_iv_eval: bad extents");
  Kokkos::Profiling::pushRegion("ofk_iv_eval");
  Timer tm(config().profile);
  const ofk::Mode mode = ofk::mode();
  const Src& s = h.src;
  const ofk_bsd::Consts<T> k = ofk_bsd::make_consts<T>(minnorm);
  const ofk_bsd::Tuning& tune = ofk_bsd::tuning();
  DV2<const T> CP;
  DV2<T> U;
  if constexpr (kOnHost) {
    // Fortran memory directly: no copy in, no copy out, the kernel accumulates into Uind.
    CP = DV2<const T>(CPs, ldCP, nCPs);
    U = DV2<T>(Uind, ldU, nCPs);
  } else {
    const std::size_t n3 = 3 * std::size_t(nCPs);
    const std::size_t c3 = 3 * std::size_t(std::max(nCPs, capCP));  // capCP = size(CPs,2) (P6)
    ensure1(h.cp, c3, "ofk_iv_cp", h.st);
    ensure1(h.hcp, c3, "ofk_iv_cp_stage", h.st);
    ensure1(h.u, c3, "ofk_iv_u", h.st);
    ensure1(h.hu, c3, "ofk_iv_u_stage", h.st);
    T* hc = h.hcp.data();
    for (int i = 0; i < nCPs; ++i)
      for (int c = 0; c < 3; ++c) hc[3 * std::size_t(i) + c] = CPs[std::size_t(i) * ldCP + c];
    const auto r = Kokkos::make_pair(std::size_t(0), n3);
    Kokkos::deep_copy(Kokkos::subview(h.cp, r), Kokkos::subview(h.hcp, r));
    // The kernel computes U = U + s with U = +0, i.e. s exactly (s is never -0); the host then
    // adds s to Uind: the same single addition the Fortran does, so no copy-in of Uind.
    Kokkos::deep_copy(Kokkos::subview(h.u, r), T(0));
    h.st.h2d += (long long)(n3 * sizeof(T));
    CP = DV2<const T>(h.cp.data(), 3, nCPs);
    U = DV2<T>(h.u.data(), 3, nCPs);
  }
  tm.lap(h.st.tCpH2D);
  if (s.method == kMethodTree) {
    ofk_tc::EvalArgs<T> a;
    a.nCP = nCPs; a.reg = s.regPart; a.k = k; a.CP = CP; a.U = U; a.bf = bf; a.dd = dd;
    // PartRegFloorFactor (FVW_BiotSavart.f90:39): None 0, Compact PART_REG_C2, default PART_REG_NRAD
    a.ff = s.regPart == ofk_bsd::kRegNone ? T(0) : (s.regPart == ofk_bsd::kRegPartCompact ? k.part_reg_c2 : T(2));
    const ofk_tc::FlatTree<T> f = flat_tree(h);
    if (mode == ofk::Mode::strict) ofk_tc::strict::eval<T>(f, a);
    else ofk_tc::fast::eval<T>(f, a);
  } else if (s.method == kMethodSeg) {
    SegArgs<T> a;
    a.nCP = nCPs; a.nSeg = s.nSeg; a.reg = s.reg; a.k = k; a.tune = tune;
    a.CP = CP; a.U = U;
    // The skinny shape's partial sums (fast mode; skinny is chosen only for nCP <= skinny_max_cp
    // unless forced): sized for the largest source set and target count it can see (P6).
    const bool skinnyPossible = tune.shape == ofk_bsd::Shape::kSkinny ||
                                (tune.shape == ofk_bsd::Shape::kAuto && nCPs <= tune.skinny_max_cp);
    if (mode == ofk::Mode::fast && s.nSeg > 0 && skinnyPossible) {
      const std::size_t segCap = std::size_t(std::max<long>(s.nSeg, h.capSeg));
      const std::size_t nChunk = (segCap + tune.chunk - 1) / tune.chunk;
      const std::size_t cpCap = tune.shape == ofk_bsd::Shape::kSkinny
                                    ? std::size_t(std::max(nCPs, capCP))
                                    : std::size_t(std::min(std::max(nCPs, capCP), tune.skinny_max_cp));
      ensure1(h.scratch, 3 * nChunk * std::max<std::size_t>(cpCap, nCPs), "ofk_iv_scratch", h.st);
    }
    if (mode == ofk::Mode::strict) ofk_bsd::strict::run_seg_soa<T>(a, h.seg, h.scratch);
    else ofk_bsd::fast::run_seg_soa<T>(a, h.seg, h.scratch);
  } else {
    PartArgs<T> a;
    a.nCP = nCPs; a.nPart = s.nPart; a.reg = s.regPart; a.k = k; a.tune = tune;
    a.CP = CP; a.U = U;
    if (mode == ofk::Mode::strict) ofk_bsd::strict::run_part_soa<T>(a, h.part);
    else ofk_bsd::fast::run_part_soa<T>(a, h.part);
  }
  tm.lap(h.st.tKernel);
  if constexpr (!kOnHost) {
    const std::size_t n3 = 3 * std::size_t(nCPs);
    const auto r = Kokkos::make_pair(std::size_t(0), n3);
    Kokkos::deep_copy(Kokkos::subview(h.hu, r), Kokkos::subview(h.u, r));
    h.st.d2h += (long long)(n3 * sizeof(T));
    tm.lap(h.st.tD2H);
    const T* hu = h.hu.data();
    for (int i = 0; i < nCPs; ++i)
      for (int c = 0; c < 3; ++c) {
        T& o = Uind[std::size_t(i) * ldU + c];
        o = o + hu[3 * std::size_t(i) + c];
      }
    tm.lap(h.st.tHostAdd);
  } else {
    ExecSpace().fence("ofk_iv_eval");
  }
  h.st.nEval++;
  Kokkos::Profiling::popRegion();
}

[[noreturn]] void fatal(const char* where, const char* what) {
  std::fprintf(stderr, "[ERROR] %s: %s\n", where, what);
  std::fflush(stderr);
  std::abort();
}

}  // namespace ofk_iv

// ---------------------------------------------------------------------------------------------
#define OFK_IV_GUARD(NAME, CALL)                     \
  try {                                              \
    ofk::ensure_initialized();                       \
    CALL;                                            \
  } catch (const std::exception& e) {                \
    ofk_iv::fatal(NAME, e.what());                   \
  } catch (...) {                                    \
    ofk_iv::fatal(NAME, "unknown C++ exception");    \
  }

#define OFK_IV_RK(NAME, RK, BODY)                                      \
  do {                                                                 \
    if ((RK) == 8) { using T = double; BODY; }                         \
    else if ((RK) == 4) { using T = float; BODY; }                     \
    else throw std::invalid_argument(std::string(NAME) + ": rk must be 4 or 8"); \
  } while (0)

extern "C" {

int ofk_iv_use_lattice(void) { return ofk_iv::use_lattice() ? 1 : 0; }

void ofk_iv_lattice(const void* key, int rk, int iW, int nWings, int nNWact, int nSpanNWP,
                    int nDepNWP, int nSpanNWG, int nDepNWG, const void* rNW, const void* gNW,
                    const void* eNW, int nFWact, int nSpanFWP, int nDepFWP, int nSpanFWG,
                    int nDepFWG, const void* rFW, const void* gFW, const void* eFW) {
  OFK_IV_GUARD("ofk_iv_lattice", OFK_IV_RK("ofk_iv_lattice", rk, {
    ofk_iv::lattice<T>(ofk_iv::get<T>(key), iW, nWings, nNWact, nSpanNWP, nDepNWP, nSpanNWG,
                       nDepNWG, (const T*)rNW, (const T*)gNW, (const T*)eNW, nFWact, nSpanFWP,
                       nDepFWP, nSpanFWG, nDepFWG, (const T*)rFW, (const T*)gFW, (const T*)eFW);
  }))
}

void ofk_iv_sources_lattice(const void* key, int rk, int method, int reg, int regPart,
                            int partPerSeg, int iDepthStart, int nNW, int nFW, int mirror,
                            int fwShed, int lastNWShed, int nSegF, int nPartF) {
  OFK_IV_GUARD("ofk_iv_sources_lattice", OFK_IV_RK("ofk_iv_sources_lattice", rk, {
    ofk_iv::sources_lattice<T>(ofk_iv::get<T>(key), method, reg, regPart, partPerSeg, iDepthStart,
                               nNW, nFW, mirror != 0, fwShed != 0, lastNWShed != 0, nSegF, nPartF);
  }))
}

void ofk_iv_sources_packed(const void* key, int rk, int method, int reg, int nSeg, int ldPts,
                           int nPts, const void* Pts, int ldC, const int* Connct,
                           const void* Gamma, const void* Eps, int regPart, int nPart, int ldP,
                           const void* PartP, int ldA, const void* PartAlpha,
                           const void* PartReg) {
  OFK_IV_GUARD("ofk_iv_sources_packed", OFK_IV_RK("ofk_iv_sources_packed", rk, {
    ofk_iv::sources_packed<T>(ofk_iv::get<T>(key), method, reg, nSeg, ldPts, nPts, (const T*)Pts,
                              ldC, Connct, (const T*)Gamma, (const T*)Eps, regPart, nPart, ldP,
                              (const T*)PartP, ldA, (const T*)PartAlpha, (const T*)PartReg);
  }))
}

int ofk_iv_sources_valid(const void* key, int method, int nSeg, int nPart) {
  int ok = 0;
  OFK_IV_GUARD("ofk_iv_sources_valid", {
    auto chk = [&](const auto* h) {
      if (h && h->src.valid && h->src.method == method &&
          (method == ofk_iv::kMethodSeg ? h->src.nSeg == nSeg : h->src.nPart == nPart) &&
          (method != ofk_iv::kMethodTree || h->treeValid))
        ok = 1;
    };
    chk(ofk_iv::find<double>(key));
    chk(ofk_iv::find<float>(key));
  })
  return ok;
}

void ofk_iv_eval(const void* key, int rk, int nCPs, int capCP, int ldCP, const void* CPs, int ldU,
                 void* Uind, const void* minnorm) {
  OFK_IV_GUARD("ofk_iv_eval", OFK_IV_RK("ofk_iv_eval", rk, {
    ofk_iv::eval<T>(ofk_iv::get<T>(key), nCPs, capCP, ldCP, (const T*)CPs, ldU, (T*)Uind,
                    *(const T*)minnorm);
  }))
}

void ofk_iv_eval_tree(const void* key, int rk, int nCPs, int capCP, int ldCP, const void* CPs,
                      int ldU, void* Uind, const void* minnorm, const void* branchFactor,
                      const void* distanceDirect) {
  OFK_IV_GUARD("ofk_iv_eval_tree", OFK_IV_RK("ofk_iv_eval_tree", rk, {
    auto& h = ofk_iv::get<T>(key);
    if (!h.src.valid || h.src.method != ofk_iv::kMethodTree)
      throw std::logic_error("no tree source set for this instance");
    ofk_iv::eval<T>(h, nCPs, capCP, ldCP, (const T*)CPs, ldU, (T*)Uind, *(const T*)minnorm,
                    *(const T*)branchFactor, *(const T*)distanceDirect);
  }))
}

void ofk_iv_tree_build(const void* key, int rk, const void* degenerate) {
  OFK_IV_GUARD("ofk_iv_tree_build", OFK_IV_RK("ofk_iv_tree_build", rk, {
    ofk_iv::tree_build<T>(ofk_iv::get<T>(key), *(const T*)degenerate);
  }))
}

int ofk_tc_enabled(void) { return ofk_iv::config().treeKokkos ? 1 : 0; }

void ofk_iv_reserve(const void* key, int rk, long long latTot, long long capSeg, long long capPart,
                    int nBlocks) {
  OFK_IV_GUARD("ofk_iv_reserve", OFK_IV_RK("ofk_iv_reserve", rk, {
    ofk_iv::reserve<T>(ofk_iv::get<T>(key), long(latTot), long(capSeg), long(capPart), nBlocks);
  }))
}

void ofk_iv_release_sources(const void* key) {
  OFK_IV_GUARD("ofk_iv_release_sources", {
    if (auto* h = ofk_iv::find<double>(key)) { h->src.valid = false; h->treeValid = false; }
    if (auto* h = ofk_iv::find<float>(key)) { h->src.valid = false; h->treeValid = false; }
  })
}

void ofk_iv_free(const void* key) {
  OFK_IV_GUARD("ofk_iv_free", {
    auto& r = ofk_iv::registry();
    std::lock_guard<std::mutex> lock(r.mu);
    auto it = r.map.find(key);
    if (it != r.map.end()) {
      if (ofk_iv::config().profile) ofk_iv::report_entry(it->second, "free");
      r.map.erase(it);
    }
  })
}

void ofk_iv_counters(const void* key, long long* out) {
  for (int i = 0; i < 9; ++i) out[i] = 0;
  OFK_IV_GUARD("ofk_iv_counters", {
    auto add = [&](const auto* h) {
      if (!h) return;
      const ofk_iv::Stats& s = h->st;
      out[0] += s.nAlloc; out[1] += s.bytesAlloc; out[2] += s.nSrc; out[3] += s.nEval;
      out[4] += s.h2d; out[5] += s.d2h; out[6] += s.latCmp; out[7] += s.latUp;
    };
    add(ofk_iv::find<double>(key));
    add(ofk_iv::find<float>(key));
    auto& r = ofk_iv::registry();
    std::lock_guard<std::mutex> lock(r.mu);
    out[8] = (long long)r.map.size();
  })
}

void ofk_iv_set_sources_mode(int m) { ofk_iv::g_sources_override.store((m >= 0 && m <= 2) ? m : -1); }

void ofk_iv_report(void) {
  auto& r = ofk_iv::registry();
  std::lock_guard<std::mutex> lock(r.mu);
  for (auto& kv : r.map) ofk_iv::report_entry(kv.second, "report");
}

}  // extern "C"
