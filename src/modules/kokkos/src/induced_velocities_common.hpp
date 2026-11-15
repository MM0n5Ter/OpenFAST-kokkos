// Unit induced-velocities: types shared by the entry TU (induced_velocities.cpp) and the pack TU
// (induced_velocities_pack.cpp). Data only.
//
// The device-side packing (fusion F1 of notes/olaf-placement.md §5) reproduces, from the
// lattice arrays x%W(iW)%{r,Gamma,Eps}_{NW,FW}, the segments that the Fortran
//   PackPanelsToSegments (FVW_Subs.f90:913) -> LatticeToSegments (FVW_VortexTools.f90:283)
//   [+ ground mirror, FVW_Subs.f90:973-989]
// writes into Sgmt%Points/Connct/Gamma/Epsilon, but directly in the SoA layout the kernels read
// (x1 y1 z1 x2 y2 z2 gamma eps), in the Fortran's segment order: the serial iHeadC counter is
// replaced by a closed-form offset per panel (criteria F3; checked against the Fortran count in
// every call). SegmentsToPart (FVW_VortexTools.f90:396) is the second kernel: particle k of
// segment j is k = j*nPartPerSeg + iSub (closed form).
// Both are elementwise copies/arithmetic in the Fortran's operation order, compiled without
// contraction in every mode: class B.
#pragma once

#include "biotsavart_direct_common.hpp"

namespace ofk_iv {

using ofk_bsd::DV1;
using ofk_bsd::DV2;
using ofk_bsd::ExecSpace;
using ofk_bsd::kOnHost;
using ofk_bsd::MemSpace;

// One call of LatticeToSegments: a (wing, NW|FW) lattice, viewed as Fortran arrays
//   r(3, rLd2, *), Gamma(gLd1, *), Eps(3, eLd2, *)
// stored at element offsets rOff/gOff/eOff of the flat lattice buffer. Indices below are the
// Fortran 1-based (iSpan, iDepth) of LatticeToSegments.
struct Block {
  int nSpan = 0;        // size(LatticePoints,2): points along the span
  int nDepth = 0;       // size(LatticePoints,3): active points along the depth (nNW+1 | nFW+1)
  int iDepthStart = 1;  // first depth index packed
  int shed = 1, shedLast = 1, hackEps = 0;  // bShedVorticity, bShedLastVorticity, bHackEpsilon
  long rOff = 0, gOff = 0, eOff = 0;
  int rLd2 = 0, gLd1 = 0, eLd2 = 0;
  int panelBase = 0;    // first panel (0-based, over all blocks) of this block
  int nPanels = 0;      // (nSpan-1) * (nDepth-iDepthStart)
  int segBase = 0;      // first segment (0-based) of this block
  int nSeg = 0;         // segments written by this block (closed form)
};

using BlockView = Kokkos::View<const Block*, MemSpace>;

// Number of segments LatticeToSegments emits for a block (FVW_VortexTools.f90:327-392):
// per panel row: (nSpan-1) * (1 + shed) + 1 (the 2-3 segment of the last span panel);
// last row: + (nSpan-1) if shed && shedLast (the 4-3 segments).
inline int block_segments(const Block& b) {
  const int rows = b.nDepth - b.iDepthStart;
  const int nS = b.nSpan - 1;
  if (rows <= 0 || nS <= 0) return 0;
  const int perRow = nS * (b.shed ? 2 : 1) + 1;
  return rows * perRow + ((b.shed && b.shedLast) ? nS : 0);
}

namespace pack {
// Lattice -> SoA segments S(j, 0..7) for j < nSegHalf, and (mirror) the ground image at
// j + nSegHalf: z -> -z, gamma -> -gamma (FVW_Subs.f90:976-988).
template <class T>
void lattice_to_segments(const DV1<T>& lat, const BlockView& blocks, int nBlocks, int nPanels,
                         int nSegHalf, bool mirror, const DV2<T>& S);
// Segments S(0..nSeg-1) -> particles PS(j*n + iSub, 0..6), n = nPartPerSeg.
template <class T>
void segments_to_particles(const DV2<T>& S, int nSeg, int nPartPerSeg, const DV2<T>& PS);
}  // namespace pack

}  // namespace ofk_iv
