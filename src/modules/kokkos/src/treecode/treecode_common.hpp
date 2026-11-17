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

#include "biotsavart-direct/biotsavart_direct_common.hpp"

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

// Device build buffers of the fast mode (treecode_dbuild.hpp). Level-synchronous: every node
// ever created is a BFS record k (the records of one level are contiguous, children of a record
// are contiguous and in octant order); the flat preorder tree is written at the end.
//   ri(k, .)   : see dbuild::I_* (range in perm, status, leaf count, children, octant counts,
//                chunk layout, subtree size / leaves, preorder index, leaf start)
//   rg(k, .)   : geometric center (3), radius R (the Fortran's float recursion)
//   rnr, rmom  : the node's kNR / kNM values (vorticity barycenter, R, maxRegParam; moments)
//   cpa, ccnt, cbase, cmom: per chunk (kCH particle positions of one record) partial sums,
//                octant counts, octant base of the chunk inside its record, partial moments
// Sized on the first call from the capacity hint (P6); records grow by doubling as a fallback.
template <class T>
struct DevBuild {
  Kokkos::View<int*, MemSpace> perm, perm2, oct;  // per particle position
  NodeIV ri;
  NodeV<T> rg, rnr, rmom;
  NodeV<T> cpa, cmom;
  NodeIV ccnt, cbase;
  std::vector<int> lvl;  // first record of each level, then the end
  std::size_t capP = 0, capR = 0, capC = 0;
  long long nAlloc = 0, bytes = 0;
  int nNodes = 0, nLeaf = 0, depth = 0, nRec = 0;
};

namespace strict {
template <class T> void build(HostTree<T>& h, const PartCols<T>& p, const TreeConsts<T>& c);
template <class T> void eval(const FlatTree<T>& t, const EvalArgs<T>& a);
}  // namespace strict
namespace fast {
template <class T> void build(HostTree<T>& h, const PartCols<T>& p, const TreeConsts<T>& c);
template <class T> void eval(const FlatTree<T>& t, const EvalArgs<T>& a);
// Fast mode, device backends: build the tree on the device from the particle SoA P (rows >= n,
// 7 cols, device memory), n >= 2. build_device fills d (topology, nodes, d.nNodes, d.nLeaf);
// write_device then writes the flat preorder tree into out (nr/mom/ni rows >= d.nNodes, leaf
// rows >= n). Same topology rules as the Fortran (and the host build); see treecode_dbuild.hpp.
template <class T> void build_device(DevBuild<T>& d, const DV2<T>& P, int n, const TreeConsts<T>& c, long capHint);
template <class T> void write_device(const DevBuild<T>& d, const DV2<T>& P, const FlatTree<T>& out);
}  // namespace fast

}  // namespace ofk_tc
