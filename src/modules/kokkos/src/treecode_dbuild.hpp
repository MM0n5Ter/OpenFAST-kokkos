// Unit treecode, FAST mode only: the tree build on the device (included by treecode_impl.hpp in
// treecode_fast.cpp; never in the strict TU). Used on device backends (OFK_TC_DEVBUILD=auto), and
// on host backends only on request (OFK_TC_DEVBUILD=1, for testing the algorithm on the CPU).
//
// Level-synchronous: all nodes of one depth are processed together, one BFS record per node.
// Per level (records [r0, r1)):
//   A  chunk layout: every record's particle range is cut into chunks of kCH positions
//      (exclusive scan of the chunk counts; one scalar D2H)
//   K1 per chunk (one team): weighted position sums, sum of weights, maxRegParam, min/max xyz,
//      octant of every particle (from the record's GEOMETRIC center, `>`), octant counts
//   K2 per record (one team): combine its chunks -> the Fortran's per-node decisions
//      (EqualRealNos dead test, leaf rule R <= maxReg .or. n <= NCRIT, degenerate extent),
//      the vorticity barycenter, and each chunk's octant base inside the record (integer scans)
//   K3/K4 moments about the barycenter, per chunk then per record
//   K5 children: exclusive scan of the child counts (one scalar D2H), a stable partition of every
//      split record's positions by octant (integer ranks: exact), and the child records with
//      G_child = G + (0.5*(R/2))*(+-1), R/2: the Fortran's float recursion, so the octant decisions
//      (and with min/max/counts, which are exact, the whole topology) equal the Fortran's and the
//      host build's; only the barycenters and moments round differently (sum order, FMA).
// Then a bottom-up pass (subtree size, subtree leaves), a top-down pass (preorder index, leaf
// start), and write_device: the same flat preorder/skip arrays as the host build, leaves gathered
// in the same order (a leaf node's particles in increasing index order, as node%iPart; a split
// node's singletons in octant order).
// Every floating-point sum is a per-thread strided loop plus a fixed pairwise tree over a fixed
// team size: deterministic (R7). No sort library, no atomics, no per-call allocation (P6).
#pragma once

#if !OFK_BSD_STRICT

#include <limits>

namespace ofk_tc {
namespace OFK_BSD_MODE_NS {
namespace dbuild {

constexpr int kTS = kOnHost ? 1 : 128;      // team size (power of two, <= 255: byte-packed ranks)
constexpr int kCH = kOnHost ? 2048 : 512;   // particle positions per chunk
// ri columns
constexpr int I_B = 0, I_E = 1, I_ST = 2, I_LN = 3, I_NC = 4, I_FC = 5, I_CNT = 6, I_C0 = 14,
              I_NCH = 15, I_SZ = 16, I_NLS = 17, I_PRE = 18, I_LSB = 19, kRI = 20;
constexpr int kPA = 11;  // chunk partials: sw*x sw*y sw*z sw maxReg min(x y z) max(x y z)
using Policy = Kokkos::TeamPolicy<ExecSpace>;
using Member = Policy::member_type;
template <class V>
using Scratch = Kokkos::View<V*, ExecSpace::scratch_memory_space, Kokkos::MemoryUnmanaged>;

// The record that owns chunk c of the level [r0, r1) (chunk starts strictly increase).
KOKKOS_INLINE_FUNCTION int rec_of_chunk(const NodeIV& ri, int r0, int r1, int c) {
  int lo = r0, hi = r1 - 1;
  while (lo < hi) {
    const int mid = (lo + hi + 1) >> 1;
    if (ri(mid, I_C0) <= c) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

// Fixed pairwise tree over the team: sh[r*W + w] holds thread r's partial; the result ends in
// sh[0..W). The team size is a power of two.
template <int W, class V, class Op>
KOKKOS_INLINE_FUNCTION void tree_reduce(const Member& tm, const Scratch<V>& sh, const Op& op) {
  const int r = tm.team_rank();
  tm.team_barrier();
  for (int s = tm.team_size() / 2; s > 0; s >>= 1) {
    if (r < s)
      for (int w = 0; w < W; ++w) sh(r * W + w) = op(w, sh(r * W + w), sh((r + s) * W + w));
    tm.team_barrier();
  }
}

template <class T>
struct PAOp {  // join of the chunk partials
  KOKKOS_INLINE_FUNCTION T operator()(int w, T a, T b) const {
    if (w < 4) return a + b;
    if (w == 4 || w >= 8) return fmax_f(a, b);
    return (b < a) ? b : a;
  }
};
struct IntSum {
  KOKKOS_INLINE_FUNCTION int operator()(int, int a, int b) const { return a + b; }
};
template <class T>
struct TSum {
  KOKKOS_INLINE_FUNCTION T operator()(int, T a, T b) const { return a + b; }
};

template <class V>
void grow2(V& v, std::size_t rows, std::size_t cols, const char* name, long long& nAlloc, long long& bytes,
           bool keep) {
  if (v.extent(0) >= rows && v.extent(1) == cols) return;
  const std::size_t cap = std::max<std::size_t>(rows, 2 * v.extent(0));
  if (keep && v.extent(0) > 0) {
    Kokkos::resize(v, cap, cols);
  } else {
    v = V();
    v = V(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), cap, cols);
  }
  nAlloc++;
  bytes += (long long)(cap * cols * sizeof(typename V::value_type));
}
template <class V>
void grow1(V& v, std::size_t n, const char* name, long long& nAlloc, long long& bytes) {
  if (v.extent(0) >= n) return;
  const std::size_t cap = std::max<std::size_t>(n, 2 * v.extent(0));
  v = V();
  v = V(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), cap);
  nAlloc++;
  bytes += (long long)(cap * sizeof(typename V::value_type));
}

template <class T>
void ensure_records(DevBuild<T>& d, std::size_t nR, bool keep) {
  if (d.capR >= nR) return;
  const std::size_t cap = std::max<std::size_t>(nR, 2 * d.capR);
  grow2(d.ri, cap, kRI, "ofk_tc_db_ri", d.nAlloc, d.bytes, keep);
  grow2(d.rg, cap, 4, "ofk_tc_db_rg", d.nAlloc, d.bytes, keep);
  grow2(d.rnr, cap, kNR, "ofk_tc_db_rnr", d.nAlloc, d.bytes, keep);
  grow2(d.rmom, cap, kNM, "ofk_tc_db_rmom", d.nAlloc, d.bytes, keep);
  d.capR = d.ri.extent(0);
}
template <class T>
void ensure_chunks(DevBuild<T>& d, std::size_t nC) {
  if (d.capC >= nC) return;
  const std::size_t cap = std::max<std::size_t>(nC, 2 * d.capC);
  grow2(d.cpa, cap, kPA, "ofk_tc_db_cpa", d.nAlloc, d.bytes, false);
  grow2(d.cmom, cap, kNM, "ofk_tc_db_cmom", d.nAlloc, d.bytes, false);
  grow2(d.ccnt, cap, 8, "ofk_tc_db_ccnt", d.nAlloc, d.bytes, false);
  grow2(d.cbase, cap, 8, "ofk_tc_db_cbase", d.nAlloc, d.bytes, false);
  d.capC = d.cpa.extent(0);
}

}  // namespace dbuild

template <class T>
void build_device(DevBuild<T>& d, const DV2<T>& P, const int n, const TreeConsts<T>& c, const long capHint) {
  using namespace dbuild;
  using RP = Kokkos::RangePolicy<ExecSpace>;
  if (n < 2) throw std::invalid_argument("ofk_tc build_device: n >= 2 (n <= 1 is built on the host)");
  // ---- capacity (P6): from the hint on the first call, doubling as the fallback
  const std::size_t capP = std::size_t(std::max<long>(n, capHint));
  if (d.capP < std::size_t(n)) {
    const std::size_t cp = std::max<std::size_t>(capP, 2 * d.capP);
    grow1(d.perm, cp, "ofk_tc_db_perm", d.nAlloc, d.bytes);
    grow1(d.perm2, cp, "ofk_tc_db_perm2", d.nAlloc, d.bytes);
    grow1(d.oct, cp, "ofk_tc_db_oct", d.nAlloc, d.bytes);
    d.capP = d.perm.extent(0);
  }
  ensure_records(d, std::max<std::size_t>(64, capP / 4), true);
  ensure_chunks(d, capP / kCH + d.capR + 1);
  if (d.lvl.capacity() < 128) d.lvl.reserve(128);

  // ---- root (grow_tree_part :601-606): bounding box (min/max: exact in any order)
  T mx, Mx, my, My, mz, Mz;
  Kokkos::parallel_reduce(
      OFK_TC_LABEL("db_box"), RP(0, n),
      KOKKOS_LAMBDA(const int i, T& a, T& A, T& b, T& B, T& e, T& E) {
        const T x = P(i, 0), y = P(i, 1), z = P(i, 2);
        if (x < a) a = x;
        if (x > A) A = x;
        if (y < b) b = y;
        if (y > B) B = y;
        if (z < e) e = z;
        if (z > E) E = z;
      },
      Kokkos::Min<T>(mx), Kokkos::Max<T>(Mx), Kokkos::Min<T>(my), Kokkos::Max<T>(My),
      Kokkos::Min<T>(mz), Kokkos::Max<T>(Mz));
  const T R0 = fmax_f(fmax_f(std::abs(Mx - mx), std::abs(My - my)), std::abs(Mz - mz)) * c.rootFactor;
  const T G0[3] = {(Mx + mx) / T(2), (My + my) / T(2), (Mz + mz) / T(2)};
  const auto perm = d.perm, perm2 = d.perm2, oct = d.oct;
  {
    const auto ri = d.ri;
    const auto rg = d.rg;
    const T g0 = G0[0], g1 = G0[1], g2 = G0[2], r0v = R0;
    Kokkos::parallel_for(OFK_TC_LABEL("db_root"), RP(0, n), KOKKOS_LAMBDA(const int q) {
      perm(q) = q;
      if (q == 0) {
        ri(0, I_B) = 0; ri(0, I_E) = n;
        ri(0, I_PRE) = 0; ri(0, I_LSB) = 0;
        rg(0, 0) = g0; rg(0, 1) = g1; rg(0, 2) = g2; rg(0, 3) = r0v;
      }
    });
  }
  d.lvl.clear();
  d.lvl.push_back(0);
  d.lvl.push_back(1);

  const T big = std::numeric_limits<T>::max();
  const T eqTol = c.eqTol, degen = c.degenerate;
  const int ncrit = c.ncrit;
  for (int L = 0;; ++L) {
    const int r0 = d.lvl[L], r1 = d.lvl[L + 1], nL = r1 - r0;
    const auto ri = d.ri;
    const auto rg = d.rg, rnr = d.rnr, rmom = d.rmom;
    // ---- A: chunk layout
    int nC = 0;
    Kokkos::parallel_scan(
        OFK_TC_LABEL("db_chunks"), RP(r0, r1),
        KOKKOS_LAMBDA(const int k, int& upd, const bool fin) {
          const int nch = (ri(k, I_E) - ri(k, I_B) + kCH - 1) / kCH;
          if (fin) { ri(k, I_C0) = upd; ri(k, I_NCH) = nch; }
          upd += nch;
        },
        nC);
    ensure_chunks(d, std::size_t(nC));
    const auto cpa = d.cpa, cmom = d.cmom;
    const auto ccnt = d.ccnt, cbase = d.cbase;
    // ---- K1: chunk partials and octants
    {
      const std::size_t shm = Scratch<T>::shmem_size(kTS * kPA) + Scratch<int>::shmem_size(kTS * 8);
      Kokkos::parallel_for(
          OFK_TC_LABEL("db_k1"), Policy(nC, kTS).set_scratch_size(0, Kokkos::PerTeam(shm)),
          KOKKOS_LAMBDA(const Member& tm) {
            const int ch = tm.league_rank();
            const int k = rec_of_chunk(ri, r0, r1, ch);
            const int qb = ri(k, I_B) + (ch - ri(k, I_C0)) * kCH;
            const int qe = (qb + kCH < ri(k, I_E)) ? qb + kCH : ri(k, I_E);
            const T g0 = rg(k, 0), g1 = rg(k, 1), g2 = rg(k, 2);
            T a[kPA] = {T(0), T(0), T(0), T(0), T(0), big, big, big, -big, -big, -big};
            int cn[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            for (int q = qb + tm.team_rank(); q < qe; q += tm.team_size()) {
              const int i = perm(q);
              const T x = P(i, 0), y = P(i, 1), z = P(i, 2);
              const T a1 = P(i, 3), a2 = P(i, 4), a3 = P(i, 5);
              const T w = Kokkos::sqrt(a1 * a1 + a2 * a2 + a3 * a3);
              a[0] = a[0] + w * x;
              a[1] = a[1] + w * y;
              a[2] = a[2] + w * z;
              a[3] = a[3] + w;
              a[4] = fmax_f(a[4], P(i, 6));
              a[5] = (x < a[5]) ? x : a[5];
              a[6] = (y < a[6]) ? y : a[6];
              a[7] = (z < a[7]) ? z : a[7];
              a[8] = fmax_f(a[8], x);
              a[9] = fmax_f(a[9], y);
              a[10] = fmax_f(a[10], z);
              const int o = (x > g0 ? 1 : 0) + (y > g1 ? 2 : 0) + (z > g2 ? 4 : 0);
              oct(q) = o;
              cn[o] += 1;
            }
            const Scratch<T> sh(tm.team_scratch(0), kTS * kPA);
            const Scratch<int> si(tm.team_scratch(0), kTS * 8);
            const int r = tm.team_rank();
            for (int w = 0; w < kPA; ++w) sh(r * kPA + w) = a[w];
            for (int o = 0; o < 8; ++o) si(r * 8 + o) = cn[o];
            tree_reduce<kPA>(tm, sh, PAOp<T>());
            tree_reduce<8>(tm, si, IntSum());
            if (r == 0) {
              for (int w = 0; w < kPA; ++w) cpa(ch, w) = sh(w);
              for (int o = 0; o < 8; ++o) ccnt(ch, o) = si(o);
            }
          });
    }
    // ---- K2: per record decisions (grow_tree_part_substep :732-822), chunk octant bases
    {
      const std::size_t shm = Scratch<T>::shmem_size(kTS * kPA) + Scratch<int>::shmem_size(kTS * 8) +
                              Scratch<int>::shmem_size(1);
      Kokkos::parallel_for(
          OFK_TC_LABEL("db_k2"), Policy(nL, kTS).set_scratch_size(0, Kokkos::PerTeam(shm)),
          KOKKOS_LAMBDA(const Member& tm) {
            const int k = r0 + tm.league_rank();
            const int c0 = ri(k, I_C0), nch = ri(k, I_NCH);
            T a[kPA] = {T(0), T(0), T(0), T(0), T(0), big, big, big, -big, -big, -big};
            int cn[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            const PAOp<T> op;
            for (int j = c0 + tm.team_rank(); j < c0 + nch; j += tm.team_size()) {
              for (int w = 0; w < kPA; ++w) a[w] = op(w, a[w], cpa(j, w));
              for (int o = 0; o < 8; ++o) cn[o] += ccnt(j, o);
            }
            const Scratch<T> sh(tm.team_scratch(0), kTS * kPA);
            const Scratch<int> si(tm.team_scratch(0), kTS * 8);
            const Scratch<int> sst(tm.team_scratch(0), 1);
            const int r = tm.team_rank();
            for (int w = 0; w < kPA; ++w) sh(r * kPA + w) = a[w];
            for (int o = 0; o < 8; ++o) si(r * 8 + o) = cn[o];
            tree_reduce<kPA>(tm, sh, op);
            tree_reduce<8>(tm, si, IntSum());
            if (r == 0) {
              const int np = ri(k, I_E) - ri(k, I_B);
              const T R = rg(k, 3);
              const T wTot = sh(3), maxReg = sh(4);
              int st = 0, ln = 0, nc = 0;
              // EqualRealNos(abs(wTot), 0.0_ReKi): dead node, nothing stored
              const T aa = Kokkos::abs(wTot);
              const T frac = fmax_f(Kokkos::abs(aa + T(0)), T(1));
              if (!(Kokkos::abs(aa - T(0)) <= frac * eqTol)) {
                rnr(k, 0) = sh(0) / wTot;
                rnr(k, 1) = sh(1) / wTot;
                rnr(k, 2) = sh(2) / wTot;
                rnr(k, 3) = R;
                rnr(k, 4) = maxReg;
                st = 1; ln = np;  // leaf node (early termination or degenerate extent)
                if (!(R <= maxReg || np <= ncrit)) {
                  bool degenerate = false;
                  for (int o = 0; o < 8; ++o) {
                    if (si(o) > 1 && si(o) == np) {
                      const T ext = fmax_f(fmax_f(Kokkos::abs(sh(8) - sh(5)), Kokkos::abs(sh(9) - sh(6))),
                                           Kokkos::abs(sh(10) - sh(7)));
                      if (ext < degen) degenerate = true;
                    }
                  }
                  if (!degenerate) {
                    st = 2; ln = 0;
                    for (int o = 0; o < 8; ++o) {
                      ln += (si(o) == 1) ? 1 : 0;
                      nc += (si(o) >= 2) ? 1 : 0;
                    }
                  }
                }
              }
              ri(k, I_ST) = st; ri(k, I_LN) = ln; ri(k, I_NC) = nc;
              for (int o = 0; o < 8; ++o) ri(k, I_CNT + o) = si(o);
              sst(0) = st;
            }
            tm.team_barrier();
            if (sst(0) != 2) return;
            for (int o = 0; o < 8; ++o)
              Kokkos::parallel_scan(Kokkos::TeamThreadRange(tm, nch), [&](const int j, int& upd, const bool fin) {
                const int v = ccnt(c0 + j, o);
                if (fin) cbase(c0 + j, o) = upd;
                upd += v;
              });
          });
    }
    // ---- K3: chunk moments about the barycenter (:750-772); K4: per record
    {
      const std::size_t shm = Scratch<T>::shmem_size(kTS * kNM);
      Kokkos::parallel_for(
          OFK_TC_LABEL("db_k3"), Policy(nC, kTS).set_scratch_size(0, Kokkos::PerTeam(shm)),
          KOKKOS_LAMBDA(const Member& tm) {
            const int ch = tm.league_rank();
            const int k = rec_of_chunk(ri, r0, r1, ch);
            if (ri(k, I_ST) == 0) return;
            const int qb = ri(k, I_B) + (ch - ri(k, I_C0)) * kCH;
            const int qe = (qb + kCH < ri(k, I_E)) ? qb + kCH : ri(k, I_E);
            const T v0 = rnr(k, 0), v1 = rnr(k, 1), v2 = rnr(k, 2);
            T M[kNM];
            for (int w = 0; w < kNM; ++w) M[w] = T(0);
            for (int q = qb + tm.team_rank(); q < qe; q += tm.team_size()) {
              const int i = perm(q);
              const T A[3] = {P(i, 3), P(i, 4), P(i, 5)};
              const T D[3] = {P(i, 0) - v0, P(i, 1) - v1, P(i, 2) - v2};
              for (int kk = 0; kk < 3; ++kk) {
                M[0 + kk] = M[0 + kk] + A[kk];
                M[3 + kk] = M[3 + kk] + A[kk] * D[0];
                M[6 + kk] = M[6 + kk] + A[kk] * D[1];
                M[9 + kk] = M[9 + kk] + A[kk] * D[2];
                M[12 + kk] = M[12 + kk] + A[kk] * D[0] * D[0];
                M[15 + kk] = M[15 + kk] + A[kk] * D[1] * D[0];
                M[18 + kk] = M[18 + kk] + A[kk] * D[1] * D[1];
                M[21 + kk] = M[21 + kk] + A[kk] * D[2] * D[0];
                M[24 + kk] = M[24 + kk] + A[kk] * D[2] * D[1];
                M[27 + kk] = M[27 + kk] + A[kk] * D[2] * D[2];
              }
            }
            const Scratch<T> sh(tm.team_scratch(0), kTS * kNM);
            const int r = tm.team_rank();
            for (int w = 0; w < kNM; ++w) sh(r * kNM + w) = M[w];
            tree_reduce<kNM>(tm, sh, TSum<T>());
            if (r == 0)
              for (int w = 0; w < kNM; ++w) cmom(ch, w) = sh(w);
          });
      Kokkos::parallel_for(
          OFK_TC_LABEL("db_k4"), Policy(nL, kTS).set_scratch_size(0, Kokkos::PerTeam(shm)),
          KOKKOS_LAMBDA(const Member& tm) {
            const int k = r0 + tm.league_rank();
            if (ri(k, I_ST) == 0) return;
            const int c0 = ri(k, I_C0), nch = ri(k, I_NCH);
            T M[kNM];
            for (int w = 0; w < kNM; ++w) M[w] = T(0);
            for (int j = c0 + tm.team_rank(); j < c0 + nch; j += tm.team_size())
              for (int w = 0; w < kNM; ++w) M[w] = M[w] + cmom(j, w);
            const Scratch<T> sh(tm.team_scratch(0), kTS * kNM);
            const int r = tm.team_rank();
            for (int w = 0; w < kNM; ++w) sh(r * kNM + w) = M[w];
            tree_reduce<kNM>(tm, sh, TSum<T>());
            if (r == 0)
              for (int w = 0; w < kNM; ++w) rmom(k, w) = sh(w);
          });
    }
    // ---- K5: children. Slots (scan), records grow (keep contents) if needed.
    int nNext = 0;
    Kokkos::parallel_scan(
        OFK_TC_LABEL("db_children"), RP(r0, r1),
        KOKKOS_LAMBDA(const int k, int& upd, const bool fin) {
          const int nc = ri(k, I_ST) == 2 ? ri(k, I_NC) : 0;
          if (fin) ri(k, I_FC) = r1 + upd;
          upd += nc;
        },
        nNext);
    if (nNext > 0) ensure_records(d, std::size_t(r1) + std::size_t(nNext), true);
    {
      const auto ri2 = d.ri;  // may have been reallocated
      const auto rg2 = d.rg;
      // stable partition of each split record's positions by octant (node%branches(:)%iPart)
      const std::size_t shm = Scratch<unsigned long long>::shmem_size(1);
      Kokkos::parallel_for(
          OFK_TC_LABEL("db_partition"), Policy(nC, kTS).set_scratch_size(0, Kokkos::PerTeam(shm)),
          KOKKOS_LAMBDA(const Member& tm) {
            const int ch = tm.league_rank();
            const int k = rec_of_chunk(ri2, r0, r1, ch);
            if (ri2(k, I_ST) != 2) return;
            const int b = ri2(k, I_B);
            const int qb = b + (ch - ri2(k, I_C0)) * kCH;
            const int qe = (qb + kCH < ri2(k, I_E)) ? qb + kCH : ri2(k, I_E);
            int base[8];
            int off = 0;
            for (int o = 0; o < 8; ++o) {
              base[o] = b + off + cbase(ch, o);
              off += ri2(k, I_CNT + o);
            }
            const Scratch<unsigned long long> tot(tm.team_scratch(0), 1);
            for (int t0 = qb; t0 < qe; t0 += kTS) {
              const int len = (t0 + kTS < qe) ? kTS : qe - t0;
              Kokkos::parallel_scan(Kokkos::TeamThreadRange(tm, len),
                                    [&](const int j, unsigned long long& upd, const bool fin) {
                                      const int o = oct(t0 + j);
                                      const unsigned long long v = 1ull << (8 * o);
                                      if (fin) {
                                        const int rank = int((upd >> (8 * o)) & 0xffull);
                                        perm2(base[o] + rank) = perm(t0 + j);
                                        if (j == len - 1) tot(0) = upd + v;
                                      }
                                      upd += v;
                                    });
              tm.team_barrier();
              const unsigned long long tt = tot(0);
              for (int o = 0; o < 8; ++o) base[o] += int((tt >> (8 * o)) & 0xffull);
              tm.team_barrier();
            }
          });
      Kokkos::parallel_for(
          OFK_TC_LABEL("db_copyback"), Policy(nC, kTS), KOKKOS_LAMBDA(const Member& tm) {
            const int ch = tm.league_rank();
            const int k = rec_of_chunk(ri2, r0, r1, ch);
            if (ri2(k, I_ST) != 2) return;
            const int qb = ri2(k, I_B) + (ch - ri2(k, I_C0)) * kCH;
            const int qe = (qb + kCH < ri2(k, I_E)) ? qb + kCH : ri2(k, I_E);
            for (int q = qb + tm.team_rank(); q < qe; q += tm.team_size()) perm(q) = perm2(q);
          });
      // child records in octant order (:841-852)
      if (nNext > 0)
        Kokkos::parallel_for(OFK_TC_LABEL("db_childrec"), RP(r0, r1), KOKKOS_LAMBDA(const int k) {
          if (ri2(k, I_ST) != 2) return;
          const int b = ri2(k, I_B);
          const T R = rg2(k, 3);
          const T halfSize = R / T(2);
          const T q5 = T(0.5) * halfSize;
          const T g0 = rg2(k, 0), g1 = rg2(k, 1), g2 = rg2(k, 2);
          int j = ri2(k, I_FC), off = 0;
          for (int o = 0; o < 8; ++o) {
            const int cnt = ri2(k, I_CNT + o);
            if (cnt >= 2) {
              ri2(j, I_B) = b + off;
              ri2(j, I_E) = b + off + cnt;
              rg2(j, 0) = g0 + q5 * ((o & 1) ? T(1) : T(-1));
              rg2(j, 1) = g1 + q5 * ((o & 2) ? T(1) : T(-1));
              rg2(j, 2) = g2 + q5 * ((o & 4) ? T(1) : T(-1));
              rg2(j, 3) = halfSize;
              ++j;
            }
            off += cnt;
          }
        });
    }
    if (nNext == 0) break;
    d.lvl.push_back(r1 + nNext);
  }
  const int nLev = int(d.lvl.size()) - 1;
  d.depth = nLev - 1;
  d.nRec = d.lvl.back();
  const auto ri = d.ri;
  // ---- bottom-up: subtree sizes (alive nodes) and leaves
  for (int L = nLev - 1; L >= 0; --L)
    Kokkos::parallel_for(OFK_TC_LABEL("db_up"), RP(d.lvl[L], d.lvl[L + 1]), KOKKOS_LAMBDA(const int k) {
      const int st = ri(k, I_ST);
      int sz = st > 0 ? 1 : 0, nls = st > 0 ? ri(k, I_LN) : 0;
      if (st == 2)
        for (int j = ri(k, I_FC); j < ri(k, I_FC) + ri(k, I_NC); ++j) { sz += ri(j, I_SZ); nls += ri(j, I_NLS); }
      ri(k, I_SZ) = sz;
      ri(k, I_NLS) = nls;
    });
  // ---- top-down: preorder index and leaf start (node, its own leaves, then its subtrees)
  for (int L = 0; L < nLev - 1; ++L)
    Kokkos::parallel_for(OFK_TC_LABEL("db_down"), RP(d.lvl[L], d.lvl[L + 1]), KOKKOS_LAMBDA(const int k) {
      if (ri(k, I_ST) != 2) return;
      int p = ri(k, I_PRE) + 1, l = ri(k, I_LSB) + ri(k, I_LN);
      for (int j = ri(k, I_FC); j < ri(k, I_FC) + ri(k, I_NC); ++j) {
        ri(j, I_PRE) = p;
        ri(j, I_LSB) = l;
        p += ri(j, I_SZ);
        l += ri(j, I_NLS);
      }
    });
  int root[2];
  {
    Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> h(root, 2);
    Kokkos::deep_copy(h, Kokkos::subview(ri, 0, Kokkos::make_pair(I_SZ, I_NLS + 1)));
  }
  d.nNodes = root[0];
  d.nLeaf = root[1];
}

template <class T>
void write_device(const DevBuild<T>& d, const DV2<T>& P, const FlatTree<T>& out) {
  using namespace dbuild;
  using RP = Kokkos::RangePolicy<ExecSpace>;
  if (d.nRec <= 0) return;
  const auto ri = d.ri;
  const auto rnr = d.rnr, rmom = d.rmom;
  const auto perm = d.perm;
  const auto nr = out.nr, mom = out.mom;
  const auto ni = out.ni;
  const auto leaf = out.leaf;
  Kokkos::parallel_for(OFK_TC_LABEL("db_write"), RP(0, d.nRec), KOKKOS_LAMBDA(const int k) {
    if (ri(k, I_ST) == 0) return;
    const int row = ri(k, I_PRE), ln = ri(k, I_LN);
    for (int w = 0; w < kNR; ++w) nr(row, w) = rnr(k, w);
    for (int w = 0; w < kNM; ++w) mom(row, w) = rmom(k, w);
    ni(row, 0) = ln > 0 ? ri(k, I_LSB) : -1;
    ni(row, 1) = ln;
    ni(row, 2) = row + ri(k, I_SZ);
  });
  Kokkos::parallel_for(
      OFK_TC_LABEL("db_leaves"), Policy(d.nRec, kTS), KOKKOS_LAMBDA(const Member& tm) {
        const int k = tm.league_rank();
        const int st = ri(k, I_ST);
        if (st == 0 || ri(k, I_LN) == 0) return;
        const int b = ri(k, I_B), ls = ri(k, I_LSB);
        if (st == 1) {  // node%leaves = node%iPart (increasing particle index)
          for (int q = b + tm.team_rank(); q < ri(k, I_E); q += tm.team_size()) {
            const int i = perm(q);
            for (int col = 0; col < 7; ++col) leaf(ls + q - b, col) = P(i, col);
          }
        } else if (tm.team_rank() == 0) {  // singletons, octant order (octant2leaves)
          int row = ls, off = 0;
          for (int o = 0; o < 8; ++o) {
            const int cnt = ri(k, I_CNT + o);
            if (cnt == 1) {
              const int i = perm(b + off);
              for (int col = 0; col < 7; ++col) leaf(row, col) = P(i, col);
              ++row;
            }
            off += cnt;
          }
        }
      });
}

template void build_device<double>(DevBuild<double>&, const DV2<double>&, int, const TreeConsts<double>&, long);
template void build_device<float>(DevBuild<float>&, const DV2<float>&, int, const TreeConsts<float>&, long);
template void write_device<double>(const DevBuild<double>&, const DV2<double>&, const FlatTree<double>&);
template void write_device<float>(const DevBuild<float>&, const DV2<float>&, const FlatTree<float>&);

}  // namespace OFK_BSD_MODE_NS
}  // namespace ofk_tc

#endif  // !OFK_BSD_STRICT
