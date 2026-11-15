// Unit treecode: types shared by the entry TU (induced_velocities.cpp, which owns the per-instance
// handle) and the two mode TUs (treecode_strict.cpp, treecode_fast.cpp). Data only: nothing here
// does floating-point arithmetic, so it is safe to include from TUs with different FMA flags.
//
// Translation of the OLAF particle treecode (FVW_VortexTools.f90, pin 12fa3fc):
//   grow_tree_part :557 (+ grow_tree_part_parallel :893, _rec :685, _substep :708,
//   compact_tree_part :644)  ->  build (host, Fortran order)  /  build_device (fast, device)
//   ui_tree_part :1512 (+ ui_tree_part_11 :1545, ui_leaves_nograd :1705)  ->  eval
//
// The Fortran tree is a pointer tree (T_Node%branches, %leaves, %iPart) walked recursively. Here it
// is a FLAT tree in depth-first preorder, the order in which ui_tree_part_11 visits it:
//   node i, then its subtree, then skip(i) = the first node after the subtree.
// A traversal is a loop "i = open ? i+1 : skip(i)" with no stack and no recursion (GPU-legal), and
// visits the nodes in exactly the Fortran's order. Dead nodes (zero vorticity, nPart = 0 after the
// substep) contribute nothing in the Fortran (skipped, no arithmetic) and are not stored.
// The leaf sources are packed per node in preorder, as compact_tree_part packs Tree%LeafBuf.
#pragma once

#include "biotsavart_direct_common.hpp"

#include <vector>

namespace ofk_tc {

using ofk_bsd::Consts;
using ofk_bsd::DV1;
using ofk_bsd::DV2;
using ofk_bsd::ExecSpace;
using ofk_bsd::kOnHost;
using ofk_bsd::MemSpace;

// Per node, node-major (one node's values are contiguous: a target reads a whole node).
//   nr : cx cy cz radius maxRegParam   (center = vorticity barycenter, as node%center after growth)
//   mom: Moments(1:3,1:10) column-major, i.e. mom(i, (m-1)*3 + (k-1)) = node%Moments(k,m)
//   ni : leafStart (0-based row of `leaf`, -1 if none) leafN  skip
constexpr int kNR = 5;
constexpr int kNM = 30;
constexpr int kNI = 3;
template <class T>
using NodeV = Kokkos::View<T**, Kokkos::LayoutRight, MemSpace>;
using NodeIV = Kokkos::View<int**, Kokkos::LayoutRight, MemSpace>;

// The flat tree as the evaluation reads it (views may be unmanaged wrappers of host buffers on
// host backends, or the handle's persistent device buffers).
template <class T>
struct FlatTree {
  int nNodes = 0;
  int nLeaf = 0;
  NodeV<T> nr, mom;
  NodeIV ni;
  DV2<T> leaf;  // (rows >= nLeaf, 7) SoA: px py pz ax ay az regparam (unit 1's particle layout)
};

// Build-time constants (all from the Fortran side or ReKi literals, audit items R3/13).
template <class T>
struct TreeConsts {
  T degenerate;  // FVW_VortexTools.f90:814 default-kind literal 1.0e-5 (passed from the shim)
  T rootFactor;  // :595 1.001_ReKi
  T eqTol;       // NWTC EqualRealNos: 100*epsilon/2 (dead-node test, :733)
  int ncrit;     // :26 TREE_LEAF_NCRIT = 32
};

// Host build buffers (persist in the handle: clear() keeps capacity, so no per-call allocation in
// steady state). The leaf SoA has leading dimension `ld` (= number of particles).
// A build fragment: a preorder piece of the tree (the top levels, or one subtree built by one
// host thread). Leaves are row-major here (7 per row); the assembled tree's leaves are SoA.
// In the top fragment, a node with ni[0] == -2 is a placeholder for subtree task ni[1].
template <class T>
struct Frag {
  std::vector<T> nr, mom, leaf;
  std::vector<int> ni;
  int nNodes = 0, nLeaf = 0, depth = 0;
};
template <class T>
struct Task {
  int b = 0, e = 0, depth = 0;
  T G[3] = {T(0), T(0), T(0)};
  T R = T(0);
};
template <class T>
struct HostTree {
  std::vector<T> nr, mom, leaf;
  std::vector<int> ni, idx, tmp;
  int nNodes = 0, nLeaf = 0, ld = 0, depth = 0;
  Frag<T> top;
  std::vector<Frag<T>> frags;
  std::vector<Task<T>> tasks;
  std::vector<int> map;
};

// Particle input of a build: column pointers (px py pz ax ay az regparam), host memory, ld = row
// stride. regparam is the REAL RegParam (the build uses it for maxRegParam even for RegNone).
template <class T>
struct PartCols {
  const T* c[7];
  int n = 0;
};

template <class T>
struct EvalArgs {
  int nCP = 0, reg = 0;
  DV2<const T> CP;  // (3, nCP), any leading dimension
  DV2<T> U;         // accumulated into
  Consts<T> k;
  T bf = 0, dd = 0, ff = 0;  // BranchFactor, DistanceDirect, PartRegFloorFactor(RegFunction)
};

// Device build buffers of the fast mode on device backends (treecode_dbuild.hpp). Every node
// ever created is a BFS record k (levels contiguous); the flat preorder tree is written at the end.
//   bint(k, .): 0 b, 1 e (range of perm), 2 status (0 dead, 1 leaf, 2 split), 3 leafN, 4 nChild,
//               5 firstChild, 6..13 cnt[8], 14 childBase
//   bG(k, .)  : geometric center (3), radius
//   bpre(k, .): 0 subtree size (alive nodes), 1 subtree leaves, 2 preorder index, 3 leaf start
constexpr int kBI = 15;
template <class T>
struct DevBuild {
  Kokkos::View<int*, MemSpace> perm, perm2, nodeOf, lvl, oct, excl;  // per particle position
  NodeV<T> bnr, bmom, bG;
  NodeIV bint, bpre, nodeBase;
  Kokkos::View<T*, MemSpace> box;      // root box reduction result (6)
  Kokkos::View<int*, MemSpace> scal;   // device scalars (scan totals, root size/leaves)
  std::size_t capN = 0, capB = 0;
  long long nAlloc = 0, bytes = 0;
  int depth = 0, nB = 0;
};

namespace strict {
template <class T> void build(HostTree<T>& h, const PartCols<T>& p, const TreeConsts<T>& c);
template <class T> void eval(const FlatTree<T>& t, const EvalArgs<T>& a);
}  // namespace strict
namespace fast {
template <class T> void build(HostTree<T>& h, const PartCols<T>& p, const TreeConsts<T>& c);
template <class T> void eval(const FlatTree<T>& t, const EvalArgs<T>& a);
// Fast mode, device backends: build the tree on the device from the particle SoA P (n rows, 7
// cols, device). Writes the flat tree into out (views already sized: nodes <= capNodes, leaves
// >= n). Returns false if out's node capacity is too small (nothing written; out.nNodes = the
// number needed). Same topology rules as the Fortran (and the host build); see treecode_dbuild.hpp.
template <class T> bool build_device(DevBuild<T>& d, const DV2<T>& P, int n, const TreeConsts<T>& c, FlatTree<T>& out);
}  // namespace fast

}  // namespace ofk_tc
