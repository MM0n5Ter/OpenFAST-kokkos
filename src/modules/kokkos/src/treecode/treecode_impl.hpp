// Unit treecode: the tree build and the tree evaluation, ONE source, compiled twice:
//   treecode_strict.cpp (OFK_BSD_MODE_NS=strict, contraction off)  -> class B vs the Fortran (CPU)
//   treecode_fast.cpp   (OFK_BSD_MODE_NS=fast,   contraction on)   -> class R2 vs the direct sum
// The near-field interaction is unit biotsavart-direct's part_term (biotsavart_direct_body.hpp),
// i.e. the physics is not duplicated: ui_leaves_nograd (FVW_VortexTools.f90:1705) has the same
// loop body as ui_part_nograd (FVW_BiotSavart.f90:379), which part_term translates.
// ODR: everything here lives in ofk_tc::<mode> (and part_term in ofk_bsd::<mode>, which is
// compiled with the same flags in biotsavart_direct_<mode>.cpp), see architecture.md §3.
//
// Mode differences (all named, architecture.md §2):
//  - strict: the Fortran's arithmetic in the Fortran's order everywhere (build: per node, particles
//    in node%iPart order; the weight norm via pow(.,0.5) as gfortran emits for `**0.5_ReKi`;
//    evaluation: preorder traversal, per-bucket leaf sums added into the running Uind_tmp, the
//    multipole terms in the Fortran's statement order, r**3/5/7 with gfortran's powi chains).
//  - fast: same tree (same topology rules), FMA allowed, sqrt for the norm, and on host backends
//    the leaf sums of a target run in kLanes = 8 fixed partial sums across all of its buckets
//    (fixed order, independent of the thread count) so they vectorize.
#pragma once

#include "biotsavart-direct/biotsavart_direct_body.hpp"
#include "treecode_common.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#if !OFK_BSD_STRICT && defined(_OPENMP) && !defined(__CUDACC__)
#define OFK_TC_SIMD _Pragma("omp simd")
#else
#define OFK_TC_SIMD
#endif

#define OFK_TC_STR2(x) #x
#define OFK_TC_STR(x) OFK_TC_STR2(x)
#define OFK_TC_LABEL(name) ("ofk_tc_" name "_" OFK_TC_STR(OFK_BSD_MODE_NS))

namespace ofk_tc {
namespace OFK_BSD_MODE_NS {

using ofk_bsd::kRegNone;
using ofk_bsd::kRegPartCompact;
using ofk_bsd::kRegPartExp;
using ofk_bsd::OFK_BSD_MODE_NS::fmax_f;
using ofk_bsd::OFK_BSD_MODE_NS::kStrict;
using ofk_bsd::OFK_BSD_MODE_NS::part_term;

constexpr int kLanes = kStrict ? 1 : (kOnHost ? 8 : 1);

// ==========================================================================================
// Build (host). grow_tree_part :557 + grow_tree_part_substep :708 + compact_tree_part :644.
// The Fortran grows the first two levels and then up to 64 subtrees in OpenMP loops
// (grow_tree_part_parallel :893); every node's arithmetic depends only on its own particles, in
// node%iPart order, so the tree does NOT depend on the thread count or the task order (checked:
// L0 against the OPENMP=OFF oracle and R7). Same here: the top kSplitDepth levels are built
// serially, the subtrees below in parallel (host threads, one fragment each), and the fragments
// are spliced back in preorder; the result is the serial result, bit for bit.
// node%iPart of a node = the parent's iPart restricted to the octant, in the parent's order, i.e.
// increasing particle index. A stable counting partition of idx[b,e) by octant keeps that order.
constexpr int kSplitDepth = 2;
constexpr int kSplitMinN = 8192;  // below this, one serial fragment

template <class T>
struct Builder {
  HostTree<T>& h;  // idx / tmp (disjoint ranges per task)
  Frag<T>& f;      // output fragment
  const PartCols<T>& p;
  const TreeConsts<T>& c;
  int splitDepth;  // depth at which subtrees become tasks (-1: never)

  T px(int i) const { return p.c[0][i]; }
  T py(int i) const { return p.c[1][i]; }
  T pz(int i) const { return p.c[2][i]; }

  int new_node() {
    const int me = f.nNodes++;
    f.nr.resize(std::size_t(f.nNodes) * kNR);
    f.mom.resize(std::size_t(f.nNodes) * kNM);
    f.ni.resize(std::size_t(f.nNodes) * kNI);
    for (int q = 0; q < kNM; ++q) f.mom[std::size_t(me) * kNM + q] = T(0);  // node%Moments = 0
    f.ni[std::size_t(me) * kNI + 0] = -1;
    f.ni[std::size_t(me) * kNI + 1] = 0;
    return me;
  }
  void push_leaf(int i) {  // compact_tree_part: LeafBuf row = [Px Py Pz Ax Ay Az Rc]
    const std::size_t r = std::size_t(f.nLeaf++);
    f.leaf.resize((r + 1) * 7);
    for (int col = 0; col < 7; ++col) f.leaf[r * 7 + col] = p.c[col][i];
  }
  void leaves_all(int me, int b, int e) {  // node%leaves = node%iPart
    f.ni[std::size_t(me) * kNI + 0] = f.nLeaf;
    f.ni[std::size_t(me) * kNI + 1] = e - b;
    for (int q = b; q < e; ++q) push_leaf(h.idx[q]);
  }

  // grow_tree_part_substep for the node holding particles idx[b,e), geometric center G, radius R.
  void node(int b, int e, const T G[3], const T R, int depth) {
    if (depth == splitDepth) {  // a subtree task: placeholder slot in this (top) fragment
      Task<T> t;
      t.b = b; t.e = e; t.depth = depth; t.R = R;
      t.G[0] = G[0]; t.G[1] = G[1]; t.G[2] = G[2];
      const int me = new_node();
      f.ni[std::size_t(me) * kNI + 0] = -2;
      f.ni[std::size_t(me) * kNI + 1] = int(h.tasks.size());
      f.ni[std::size_t(me) * kNI + 2] = f.nNodes;
      h.tasks.push_back(t);
      return;
    }
    if (depth > f.depth) f.depth = depth;
    const int np = e - b;
    // --- barycenter of vorticity (:732-740)
    T VortC[3] = {T(0), T(0), T(0)};
    T wTot = T(0);
    T maxReg = T(0);
    for (int q = b; q < e; ++q) {
      const int i = h.idx[q];
      const T a1 = p.c[3][i], a2 = p.c[4][i], a3 = p.c[5][i];
      const T s = a1 * a1 + a2 * a2 + a3 * a3;
      T wLoc;
      if constexpr (kStrict) wLoc = std::pow(s, T(0.5));  // gfortran: `**0.5_ReKi` -> pow
      else wLoc = std::sqrt(s);
      VortC[0] = VortC[0] + wLoc * px(i);
      VortC[1] = VortC[1] + wLoc * py(i);
      VortC[2] = VortC[2] + wLoc * pz(i);
      wTot = wTot + wLoc;
      maxReg = fmax_f(maxReg, p.c[6][i]);
    }
    // --- EqualRealNos(abs(wTot), 0.0_ReKi) (NWTC_Num.f90:1647): dead node, nothing stored (:742-746)
    {
      const T a = std::abs(wTot);
      const T frac = fmax_f(std::abs(a + T(0)), T(1));
      if (std::abs(a - T(0)) <= frac * c.eqTol) return;
    }
    VortC[0] = VortC[0] / wTot;
    VortC[1] = VortC[1] / wTot;
    VortC[2] = VortC[2] / wTot;
    const int me = new_node();
    T* nr = &f.nr[std::size_t(me) * kNR];
    nr[0] = VortC[0]; nr[1] = VortC[1]; nr[2] = VortC[2]; nr[3] = R; nr[4] = maxReg;
    // --- moments about VortC (:750-772); M(k,m) at mom[(m-1)*3 + k-1]
    T* M = &f.mom[std::size_t(me) * kNM];
    for (int q = b; q < e; ++q) {
      const int i = h.idx[q];
      const T A[3] = {p.c[3][i], p.c[4][i], p.c[5][i]};
      const T D[3] = {px(i) - VortC[0], py(i) - VortC[1], pz(i) - VortC[2]};
      for (int k = 0; k < 3; ++k) {
        M[0 + k] = M[0 + k] + A[k];                 // M0_000
        M[3 + k] = M[3 + k] + A[k] * D[0];          // M1_100
        M[6 + k] = M[6 + k] + A[k] * D[1];          // M1_010
        M[9 + k] = M[9 + k] + A[k] * D[2];          // M1_001
        M[12 + k] = M[12 + k] + A[k] * D[0] * D[0];  // 5:  j=1,k=1
        M[15 + k] = M[15 + k] + A[k] * D[1] * D[0];  // 6:  j=2,k=1
        M[18 + k] = M[18 + k] + A[k] * D[1] * D[1];  // 7:  j=2,k=2
        M[21 + k] = M[21 + k] + A[k] * D[2] * D[0];  // 8:  j=3,k=1
        M[24 + k] = M[24 + k] + A[k] * D[2] * D[1];  // 9:  j=3,k=2
        M[27 + k] = M[27 + k] + A[k] * D[2] * D[2];  // 10: j=3,k=3
      }
    }
    // --- early leaf termination (:775): radius <= 1.0*maxRegParam .or. nPart <= NCRIT
    if (R <= maxReg || np <= c.ncrit) {
      leaves_all(me, b, e);
      f.ni[std::size_t(me) * kNI + 2] = f.nNodes;
      return;
    }
    // --- octants from the GEOMETRIC center (:782-793)
    int cnt[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int q = b; q < e; ++q) {
      const int i = h.idx[q];
      const int o = (px(i) > G[0] ? 1 : 0) + (py(i) > G[1] ? 2 : 0) + (pz(i) > G[2] ? 4 : 0);
      h.tmp[q] = o;
      cnt[o]++;
    }
    // --- all particles in one octant and a degenerate extent: all leaves (:806-822)
    for (int o = 0; o < 8; ++o) {
      if (cnt[o] > 1 && cnt[o] == np) {
        T mx = px(h.idx[b]), Mx = mx, my = py(h.idx[b]), My = my, mz = pz(h.idx[b]), Mz = mz;
        for (int q = b; q < e; ++q) {
          const int i = h.idx[q];
          Mx = fmax_f(Mx, px(i)); mx = (px(i) < mx) ? px(i) : mx;
          My = fmax_f(My, py(i)); my = (py(i) < my) ? py(i) : my;
          Mz = fmax_f(Mz, pz(i)); mz = (pz(i) < mz) ? pz(i) : mz;
        }
        const T ext = fmax_f(fmax_f(std::abs(Mx - mx), std::abs(My - my)), std::abs(Mz - mz));
        if (ext < c.degenerate) {
          leaves_all(me, b, e);
          f.ni[std::size_t(me) * kNI + 2] = f.nNodes;
          return;
        }
      }
    }
    // --- stable partition of idx[b,e) by octant (node%branches(:)%iPart, :856-869); the second
    // half of tmp (2n entries) is the scratch
    int off[9];
    off[0] = 0;
    for (int o = 0; o < 8; ++o) off[o + 1] = off[o] + cnt[o];
    {
      int pos[8];
      for (int o = 0; o < 8; ++o) pos[o] = b + off[o];
      int* scratch = &h.tmp[std::size_t(h.ld)];
      for (int q = b; q < e; ++q) scratch[pos[h.tmp[q]]++] = h.idx[q];
      for (int q = b; q < e; ++q) h.idx[q] = scratch[q];
    }
    // --- singletons are this node's leaves, numbered by octant (octant2leaves, :800)
    int nLeaves = 0;
    for (int o = 0; o < 8; ++o) nLeaves += (cnt[o] == 1) ? 1 : 0;
    if (nLeaves > 0) {
      f.ni[std::size_t(me) * kNI + 0] = f.nLeaf;
      f.ni[std::size_t(me) * kNI + 1] = nLeaves;
      for (int o = 0; o < 8; ++o)
        if (cnt[o] == 1) push_leaf(h.idx[b + off[o]]);
    }
    // --- branches in octant order (:841-852): radius halfSize, center GeomC + 0.5*halfSize*(+-1)
    const T halfSize = R / T(2);
    const T q5 = T(0.5) * halfSize;
    for (int o = 0; o < 8; ++o) {
      if (cnt[o] < 2) continue;
      const T C[3] = {G[0] + q5 * ((o & 1) ? T(1) : T(-1)), G[1] + q5 * ((o & 2) ? T(1) : T(-1)),
                      G[2] + q5 * ((o & 4) ? T(1) : T(-1))};
      node(b + off[o], b + off[o + 1], C, halfSize, depth + 1);
    }
    f.ni[std::size_t(me) * kNI + 2] = f.nNodes;
  }
};

template <class T>
void reset(Frag<T>& f) {
  f.nr.clear(); f.mom.clear(); f.leaf.clear(); f.ni.clear();
  f.nNodes = 0; f.nLeaf = 0; f.depth = 0;
}

// Splice the top fragment and the task fragments into the final preorder arrays of h.
template <class T>
void assemble(HostTree<T>& h) {
  const Frag<T>& top = h.top;
  h.map.resize(std::size_t(top.nNodes) + 1);
  int pos = 0, nLeafTot = 0;
  for (int t = 0; t < top.nNodes; ++t) {
    h.map[t] = pos;
    const int* ni = &top.ni[std::size_t(t) * kNI];
    if (ni[0] == -2) { pos += h.frags[ni[1]].nNodes; nLeafTot += h.frags[ni[1]].nLeaf; }
    else { pos += 1; nLeafTot += ni[1]; }
  }
  h.map[top.nNodes] = pos;
  h.nNodes = pos;
  h.nLeaf = nLeafTot;
  h.nr.resize(std::size_t(pos) * kNR);
  h.mom.resize(std::size_t(pos) * kNM);
  h.ni.resize(std::size_t(pos) * kNI);
  const std::size_t ld = std::size_t(h.ld);
  int lp = 0;  // next leaf row
  auto copy_leaves = [&](const Frag<T>& f, int from, int n) {
    for (int r = 0; r < n; ++r)
      for (int col = 0; col < 7; ++col) h.leaf[std::size_t(lp + r) + col * ld] = f.leaf[std::size_t(from + r) * 7 + col];
  };
  for (int t = 0; t < top.nNodes; ++t) {
    const int* ni = &top.ni[std::size_t(t) * kNI];
    const int o = h.map[t];
    if (ni[0] == -2) {
      const Frag<T>& f = h.frags[ni[1]];
      for (int j = 0; j < f.nNodes; ++j) {
        for (int q = 0; q < kNR; ++q) h.nr[std::size_t(o + j) * kNR + q] = f.nr[std::size_t(j) * kNR + q];
        for (int q = 0; q < kNM; ++q) h.mom[std::size_t(o + j) * kNM + q] = f.mom[std::size_t(j) * kNM + q];
        const int* fj = &f.ni[std::size_t(j) * kNI];
        int* hj = &h.ni[std::size_t(o + j) * kNI];
        hj[0] = fj[1] > 0 ? fj[0] + lp : -1;
        hj[1] = fj[1];
        hj[2] = fj[2] + o;
      }
      copy_leaves(f, 0, f.nLeaf);
      lp += f.nLeaf;
      if (f.depth > h.depth) h.depth = f.depth;
    } else {
      for (int q = 0; q < kNR; ++q) h.nr[std::size_t(o) * kNR + q] = top.nr[std::size_t(t) * kNR + q];
      for (int q = 0; q < kNM; ++q) h.mom[std::size_t(o) * kNM + q] = top.mom[std::size_t(t) * kNM + q];
      int* hj = &h.ni[std::size_t(o) * kNI];
      hj[1] = ni[1];
      hj[0] = ni[1] > 0 ? lp : -1;
      hj[2] = h.map[ni[2]];
      copy_leaves(top, ni[0], ni[1]);
      lp += ni[1];
    }
  }
  if (top.depth > h.depth) h.depth = top.depth;
}

template <class T>
void build(HostTree<T>& h, const PartCols<T>& p, const TreeConsts<T>& c) {
  const int n = p.n;
  h.nNodes = 0;
  h.nLeaf = 0;
  h.depth = 0;
  h.ld = n;
  h.tasks.clear();
  reset(h.top);
  if (n <= 0) return;  // :584 node%nPart = 0: nothing is ever evaluated
  if (h.leaf.size() < std::size_t(n) * 7) h.leaf.resize(std::size_t(n) * 7);
  h.idx.resize(n);
  h.tmp.resize(2 * std::size_t(n));
  for (int i = 0; i < n; ++i) h.idx[i] = i;
  if (n == 1) {  // :589-599
    Builder<T> B{h, h.top, p, c, -1};
    const int me = B.new_node();
    T* nr = &h.top.nr[0];
    nr[0] = p.c[0][0]; nr[1] = p.c[1][0]; nr[2] = p.c[2][0]; nr[3] = T(0); nr[4] = p.c[6][0];
    h.top.mom[0] = p.c[3][0]; h.top.mom[1] = p.c[4][0]; h.top.mom[2] = p.c[5][0];  // monopole only
    B.leaves_all(me, 0, 1);
    h.top.ni[2] = 1;
    assemble(h);
    return;
  }
  // :601-606 domain and root cell
  T Mx = p.c[0][0], mx = Mx, My = p.c[1][0], my = My, Mz = p.c[2][0], mz = Mz;
  for (int i = 0; i < n; ++i) {
    Mx = fmax_f(Mx, p.c[0][i]); mx = (p.c[0][i] < mx) ? p.c[0][i] : mx;
    My = fmax_f(My, p.c[1][i]); my = (p.c[1][i] < my) ? p.c[1][i] : my;
    Mz = fmax_f(Mz, p.c[2][i]); mz = (p.c[2][i] < mz) ? p.c[2][i] : mz;
  }
  const T R = fmax_f(fmax_f(std::abs(Mx - mx), std::abs(My - my)), std::abs(Mz - mz)) * c.rootFactor;
  const T G[3] = {(Mx + mx) / T(2), (My + my) / T(2), (Mz + mz) / T(2)};
  {
    Builder<T> B{h, h.top, p, c, n >= kSplitMinN ? kSplitDepth : -1};
    B.node(0, n, G, R, 0);
  }
  const int nTasks = int(h.tasks.size());
  if (int(h.frags.size()) < nTasks) h.frags.resize(nTasks);
  if (nTasks > 0) {
    Kokkos::parallel_for(
        OFK_TC_LABEL("build_subtrees"),
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>>(0, nTasks),
        [&](const int t) {
          Frag<T>& f = h.frags[t];
          reset(f);
          const Task<T>& k = h.tasks[t];
          Builder<T> B{h, f, p, c, -1};
          B.node(k.b, k.e, k.G, k.R, k.depth);
        });
    Kokkos::DefaultHostExecutionSpace().fence();
  }
  assemble(h);
}

// ==========================================================================================
// Evaluation (ui_tree_part :1512, ui_tree_part_11 :1545). One thread per target, stackless
// preorder walk, accumulate into U once per target (as `Uind(1:3,icp) = Uind + Uind_tmp`).

// The far-field branch of ui_tree_part_11 (:1572-1694), statement by statement. (x,y,z) = DeltaP.
template <class T>
KOKKOS_FORCEINLINE_FUNCTION void multipole(const T* M, const T x, const T y, const T z, const T r,
                                           const T fourpi_inv, const T c3fp, const T c75fp,
                                           T& u1, T& u2, T& u3) {
#define MM(k, m) M[((m)-1) * 3 + ((k)-1)]
  // gfortran's powi chains (tree-ssa-math-opts powi_table): r**3 = r*r^2, r**5 = r^2*r^3,
  // r**7 = r^3*r^4 (checked in the assembly of the pin's build and by L0).
  const T r2 = r * r;
  const T r3 = r * r2;
  const T r4 = r2 * r2;
  const T r5 = r2 * r3;
  const T r7 = r3 * r4;
  const T coeff3 = fourpi_inv / r3;
  // order 0
  {
    const T phi1 = MM(1, 1) * coeff3, phi2 = MM(2, 1) * coeff3, phi3 = MM(3, 1) * coeff3;
    const T U1 = phi2 * z - phi3 * y;
    const T U2 = phi3 * x - phi1 * z;
    const T U3 = phi1 * y - phi2 * x;
    u1 = u1 + U1; u2 = u2 + U2; u3 = u3 + U3;
  }
  // order 1 (coeff5 = 3*fourpi_inv/r**5, the constant product folded as gfortran does)
  const T coeff5 = c3fp / r5;
  {
    T U1 = T(0), U2 = T(0), U3 = T(0);
    T coeff = coeff5 * x;
    U1 = U1 - coeff * (MM(3, 2) * y - MM(2, 2) * z);
    U2 = U2 - coeff * (MM(1, 2) * z - MM(3, 2) * x);
    U3 = U3 - coeff * (MM(2, 2) * x - MM(1, 2) * y);
    coeff = coeff5 * y;
    U1 = U1 - coeff * (MM(3, 3) * y - MM(2, 3) * z);
    U2 = U2 - coeff * (MM(1, 3) * z - MM(3, 3) * x);
    U3 = U3 - coeff * (MM(2, 3) * x - MM(1, 3) * y);
    coeff = coeff5 * z;
    U1 = U1 - coeff * (MM(3, 4) * y - MM(2, 4) * z);
    U2 = U2 - coeff * (MM(1, 4) * z - MM(3, 4) * x);
    U3 = U3 - coeff * (MM(2, 4) * x - MM(1, 4) * y);
    U1 = U1 + coeff3 * MM(3, 3) - coeff3 * MM(2, 4);
    U2 = U2 + coeff3 * MM(1, 4) - coeff3 * MM(3, 2);
    U3 = U3 + coeff3 * MM(2, 2) - coeff3 * MM(1, 3);
    u1 = u1 + U1; u2 = u2 + U2; u3 = u3 + U3;
  }
  // order 2: the loop form the Fortran keeps (:1609-1626, "BAR_OLAF test fails ... unrolled")
  {
    T U1 = T(0), U2 = T(0), U3 = T(0);
    const T coeff = coeff5 * T(0.5);
    const T coeff7 = -(c75fp / r7);  // -7.5_ReKi*fourpi_inv/r**7 = -((7.5*fourpi_inv)/r**7)
    const T D[3] = {x, y, z};
    for (int i = 1; i <= 3; ++i) {
      const int a = 5 + 2 * (i / 2) + 3 * (i / 3);
      U1 = U1 + coeff * (y * MM(3, a) - z * MM(2, a));
      U2 = U2 + coeff * (z * MM(1, a) - x * MM(3, a));
      U3 = U3 + coeff * (x * MM(2, a) - y * MM(1, a));
      for (int j = 1; j <= i; ++j) {
        const T ieqj = (i == j) ? T(1) : T(2);
        const T coeff7ij = D[i - 1] * D[j - 1] * coeff7;
        const int b = 3 + i + j + i / 3;
        U1 = U1 + ieqj * coeff7ij * (y * MM(3, b) - z * MM(2, b));
        U2 = U2 + ieqj * coeff7ij * (z * MM(1, b) - x * MM(3, b));
        U3 = U3 + ieqj * coeff7ij * (x * MM(2, b) - y * MM(1, b));
      }
    }
    // r5 terms 2 (:1686-1691); M2_11=5 M2_21=6 M2_22=7 M2_31=8 M2_32=9 M2_33=10
    U1 = U1 + coeff5 * (y * MM(3, 7) - z * MM(2, 10));
    U2 = U2 + coeff5 * (z * MM(1, 10) - x * MM(3, 5));
    U3 = U3 + coeff5 * (x * MM(2, 5) - y * MM(1, 7));
    U1 = U1 + coeff5 * (z * MM(3, 9) + x * MM(3, 6) - x * MM(2, 8) - y * MM(2, 9));
    U2 = U2 + coeff5 * (x * MM(1, 8) + y * MM(1, 9) - y * MM(3, 6) - z * MM(3, 8));
    U3 = U3 + coeff5 * (y * MM(2, 6) + z * MM(2, 8) - z * MM(1, 9) - x * MM(1, 6));
    u1 = u1 + U1; u2 = u2 + U2; u3 = u3 + U3;
  }
#undef MM
}

template <int Reg, class T>
KOKKOS_FORCEINLINE_FUNCTION void leaf_acc(const T cx, const T cy, const T cz, const T* const* s,
                                          const int j, const Consts<T>& k, T& u1, T& u2, T& u3) {
  T o1 = T(0), o2 = T(0), o3 = T(0);
  // ui_leaves_nograd: dx = CP(1) - LeafBuf(i,1), ...
  const bool ok = part_term<Reg, T>(cx - s[0][j], cy - s[1][j], cz - s[2][j], s[3][j], s[4][j],
                                    s[5][j], s[6][j], k, o1, o2, o3);
  u1 = u1 + (ok ? o1 : T(0));
  u2 = u2 + (ok ? o2 : T(0));
  u3 = u3 + (ok ? o3 : T(0));
}

template <int L, class T>
KOKKOS_FORCEINLINE_FUNCTION T lanes_sum(T* a) {
  for (int w = L / 2; w >= 1; w /= 2)
    for (int i = 0; i < w; ++i) a[i] = a[i] + a[i + w];
  return a[0];
}

template <int Reg, class T>
void eval_reg(const FlatTree<T>& t, const EvalArgs<T>& a) {
  const int nNodes = t.nNodes;
  const auto nrv = t.nr;
  const auto momv = t.mom;
  const auto niv = t.ni;
  const auto CP = a.CP;
  const auto U = a.U;
  const Consts<T> k = a.k;
  const T bf = a.bf, dd0 = a.dd, ff = a.ff;
  const T fourpi_inv = k.fourpi_inv;
  // gfortran folds the constant products 3.0_ReKi*fourpi_inv and 7.5_ReKi*fourpi_inv
  const T c3fp = T(3) * fourpi_inv;
  const T c75fp = T(7.5) * fourpi_inv;
  const T* l0 = t.leaf.data();
  const std::size_t ld = t.leaf.stride(1);
  Kokkos::parallel_for(
      OFK_TC_LABEL("eval"), Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>>(0, a.nCP),
      KOKKOS_LAMBDA(const int icp) {
        // nvcc: no first capture inside `if constexpr` (unit 1 lesson): capture k here
        const Consts<T> kc = k;
        const T* s[7];
        for (int c = 0; c < 7; ++c) s[c] = l0 + c * ld;
        const T cx = CP(0, icp), cy = CP(1, icp), cz = CP(2, icp);
        T u1 = T(0), u2 = T(0), u3 = T(0);  // Uind_tmp
        T w1[kLanes], w2[kLanes], w3[kLanes];  // fast, host: leaf partial sums (lanes)
        for (int l = 0; l < kLanes; ++l) { w1[l] = T(0); w2[l] = T(0); w3[l] = T(0); }
        int i = 0;
        while (i < nNodes) {
          // distDirect = max(BranchFactor*radius + PartRegFloorFactor*maxRegParam, DistanceDirect)
          const T dist = fmax_f(bf * nrv(i, 3) + ff * nrv(i, 4), dd0);
          const T x = -nrv(i, 0) + cx, y = -nrv(i, 1) + cy, z = -nrv(i, 2) + cz;
          const T r = Kokkos::sqrt(x * x + y * y + z * z);
          if (r < dist) {
            const int ls = niv(i, 0), ln = niv(i, 1);
            if (ln > 0) {
              if constexpr (kLanes == 1) {
                if constexpr (kStrict) {
                  T ux = T(0), uy = T(0), uz = T(0);  // ui_leaves_nograd's own sum, then added
                  for (int j = ls; j < ls + ln; ++j) leaf_acc<Reg, T>(cx, cy, cz, s, j, kc, ux, uy, uz);
                  u1 = u1 + ux; u2 = u2 + uy; u3 = u3 + uz;
                } else {
                  for (int j = ls; j < ls + ln; ++j) leaf_acc<Reg, T>(cx, cy, cz, s, j, kc, w1[0], w2[0], w3[0]);
                }
              } else {
                const int nFull = ln - ln % kLanes;
                for (int j0 = ls; j0 < ls + nFull; j0 += kLanes) {
                  OFK_TC_SIMD
                  for (int l = 0; l < kLanes; ++l) leaf_acc<Reg, T>(cx, cy, cz, s, j0 + l, kc, w1[l], w2[l], w3[l]);
                }
                for (int l = 0; l < ln - nFull; ++l)
                  leaf_acc<Reg, T>(cx, cy, cz, s, ls + nFull + l, kc, w1[l], w2[l], w3[l]);
              }
            }
            i = i + 1;  // first branch (preorder), or skip(i) == i+1 when there is none
          } else {
            multipole<T>(&momv(i, 0), x, y, z, r, fourpi_inv, c3fp, c75fp, u1, u2, u3);
            i = niv(i, 2);
          }
        }
        if constexpr (!kStrict) {
          u1 = u1 + lanes_sum<kLanes>(w1);
          u2 = u2 + lanes_sum<kLanes>(w2);
          u3 = u3 + lanes_sum<kLanes>(w3);
        }
        U(0, icp) = U(0, icp) + u1;
        U(1, icp) = U(1, icp) + u2;
        U(2, icp) = U(2, icp) + u3;
      });
}

template <class T>
void eval(const FlatTree<T>& t, const EvalArgs<T>& a) {
  if (a.nCP <= 0) return;
  switch (a.reg) {
    case kRegNone: eval_reg<kRegNone, T>(t, a); break;
    case kRegPartExp: eval_reg<kRegPartExp, T>(t, a); break;
    case kRegPartCompact: eval_reg<kRegPartCompact, T>(t, a); break;
    default: throw std::invalid_argument("ofk_tc: wrong regularization function for particles");
  }
}

template void build<double>(HostTree<double>&, const PartCols<double>&, const TreeConsts<double>&);
template void build<float>(HostTree<float>&, const PartCols<float>&, const TreeConsts<float>&);
template void eval<double>(const FlatTree<double>&, const EvalArgs<double>&);
template void eval<float>(const FlatTree<float>&, const EvalArgs<float>&);

}  // namespace OFK_BSD_MODE_NS
}  // namespace ofk_tc
#include "treecode_dbuild.hpp"
