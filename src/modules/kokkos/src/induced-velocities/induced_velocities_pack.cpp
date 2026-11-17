// Unit induced-velocities: device-side packing (fusion F1), lattice -> segments -> particles.
// Compiled WITHOUT contraction in every mode (listed with the strict TUs in ../CMakeLists.txt,
// so nvcc gets --fmad=false): the pack is class B, identical to the Fortran on every backend
// that has IEEE division and sqrt (CUDA's default -prec-div/-prec-sqrt).
//
// Translation of (pin 12fa3fc):
//   LatticeToSegments  FVW_VortexTools.f90:283-394 (one Block per call)
//   mirror             FVW_Subs.f90:973-989
//   SegmentsToPart     FVW_VortexTools.f90:396-447
#include "induced_velocities_common.hpp"

namespace ofk_iv {
namespace pack {

template <class T>
void lattice_to_segments(const DV1<T>& lat, const BlockView& blocks, const int nBlocks,
                         const int nPanels, const int nSegHalf, const bool mirror,
                         const DV2<T>& S) {
  if (nPanels <= 0) return;
  Kokkos::parallel_for(
      "ofk_iv_pack_lattice", Kokkos::RangePolicy<ExecSpace>(0, nPanels),
      KOKKOS_LAMBDA(const int p) {
        int ib = 0;
        while (ib + 1 < nBlocks && p >= blocks(ib + 1).panelBase) ++ib;  // nBlocks <= 2*nWings
        const Block& b = blocks(ib);
        const int nS = b.nSpan - 1;                 // panels along the span
        const int local = p - b.panelBase;
        const int row = local / nS;
        const int iSpan = local % nS + 1;           // Fortran 1-based
        const int iDepth = b.iDepthStart + row;
        const int nDepth = b.nDepth;
        const T* r = lat.data() + b.rOff;
        const T* g = lat.data() + b.gOff;
        const T* e = lat.data() + b.eOff;
        // r(c, is, id), Gamma(is, id), Eps(c, is, id); c, is, id 1-based
        auto R = [&](int c, int is, int id) { return r[(c - 1) + 3 * ((is - 1) + long(b.rLd2) * (id - 1))]; };
        auto G = [&](int is, int id) { return g[(is - 1) + long(b.gLd1) * (id - 1)]; };
        auto E = [&](int c, int is, int id) { return e[(c - 1) + 3 * ((is - 1) + long(b.eLd2) * (id - 1))]; };

        // Closed-form iHeadC (criteria F3): panels in (iDepth, iSpan) order, each emitting
        // 1-2 (if shed), 1-4, 4-3 (last depth row, if shed and shedLast), 2-3 (last span panel).
        const int per = b.shed ? 2 : 1;
        const bool lastRow = (iDepth == nDepth - 1);
        const int extra = (b.shed && b.shedLast) ? 1 : 0;
        const int perRow = nS * per + 1;
        int j = b.segBase + row * perRow + (iSpan - 1) * (per + (lastRow ? extra : 0));

        // Point numbering of a panel: 2---3 / 1---4 (span up, depth right)
        const int s1 = iSpan, d1 = iDepth;          // point 1
        const int s2 = iSpan + 1, d2 = iDepth;      // point 2
        const int s3 = iSpan + 1, d3 = iDepth + 1;  // point 3
        const int s4 = iSpan, d4 = iDepth + 1;      // point 4

        auto put = [&](int jj, int sa, int da, int sb, int db, T gam, T eps) {
          const T xa = R(1, sa, da), ya = R(2, sa, da), za = R(3, sa, da);
          const T xb = R(1, sb, db), yb = R(2, sb, db), zb = R(3, sb, db);
          S(jj, 0) = xa; S(jj, 1) = ya; S(jj, 2) = za;
          S(jj, 3) = xb; S(jj, 4) = yb; S(jj, 5) = zb;
          S(jj, 6) = gam; S(jj, 7) = eps;
          if (mirror) {
            const int jm = jj + nSegHalf;
            S(jm, 0) = xa; S(jm, 1) = ya; S(jm, 2) = -za;
            S(jm, 3) = xb; S(jm, 4) = yb; S(jm, 5) = -zb;
            S(jm, 6) = -gam; S(jm, 7) = eps;
          }
        };

        T Gamma12, Eps12, Gamma41, Eps41;
        if (iDepth == b.iDepthStart) {
          Gamma12 = G(iSpan, iDepth);
          Eps12 = E(1, iSpan, iDepth);
        } else {
          Gamma12 = G(iSpan, iDepth) - G(iSpan, iDepth - 1);
          Eps12 = (E(1, iSpan, iDepth) + E(1, iSpan, iDepth - 1)) / T(2);
        }
        if (iSpan == 1) {
          Gamma41 = G(iSpan, iDepth);
          Eps41 = b.hackEps ? E(2, iSpan, iDepth) : E(3, iSpan, iDepth);
        } else {
          Gamma41 = G(iSpan, iDepth) - G(iSpan - 1, iDepth);
          Eps41 = (E(3, iSpan, iDepth) + E(3, iSpan - 1, iDepth)) / T(2);
        }
        if (b.shed) put(j++, s1, d1, s2, d2, Gamma12, Eps12);                // 1-2
        put(j++, s1, d1, s4, d4, -Gamma41, Eps41);                           // 1-4
        if (lastRow && b.shed && b.shedLast)
          put(j++, s4, d4, s3, d3, -G(iSpan, iDepth), E(1, iSpan, iDepth));  // 4-3
        if (iSpan == nS) put(j++, s2, d2, s3, d3, G(iSpan, iDepth), E(3, iSpan, iDepth));  // 2-3
      });
}

template <class T>
void segments_to_particles(const DV2<T>& S, const int nSeg, const int nPartPerSeg,
                           const DV2<T>& PS) {
  if (nSeg <= 0) return;
  const T n = T(nPartPerSeg);
  Kokkos::parallel_for(
      "ofk_iv_pack_particles", Kokkos::RangePolicy<ExecSpace>(0, nSeg),
      KOKKOS_LAMBDA(const int js) {
        const T P1[3] = {S(js, 0), S(js, 1), S(js, 2)};
        const T DP[3] = {S(js, 3) - P1[0], S(js, 4) - P1[1], S(js, 5) - P1[2]};
        const T SegLen = Kokkos::sqrt(DP[0] * DP[0] + DP[1] * DP[1] + DP[2] * DP[2]);
        T SegDir[3] = {T(0), T(0), T(0)};
        if (SegLen > T(0)) {
          SegDir[0] = DP[0] / SegLen; SegDir[1] = DP[1] / SegLen; SegDir[2] = DP[2] / SegLen;
        }
        const T gam = S(js, 6);
        const T PartInt[3] = {DP[0] * gam / n, DP[1] * gam / n, DP[2] * gam / n};
        const T PartEps = S(js, 7);
        const T PartLen = SegLen / n;
        for (int iSub = 0; iSub < nPartPerSeg; ++iSub) {
          const int k = js * nPartPerSeg + iSub;
          const T f = (T(0.5) + T(iSub)) * PartLen;  // (0.5_ReKi+iSubPart)*PartLen*SegDir
          PS(k, 0) = P1[0] + f * SegDir[0];
          PS(k, 1) = P1[1] + f * SegDir[1];
          PS(k, 2) = P1[2] + f * SegDir[2];
          PS(k, 3) = PartInt[0]; PS(k, 4) = PartInt[1]; PS(k, 5) = PartInt[2];
          PS(k, 6) = PartEps;
        }
      });
}

template void lattice_to_segments<double>(const DV1<double>&, const BlockView&, int, int, int, bool, const DV2<double>&);
template void lattice_to_segments<float>(const DV1<float>&, const BlockView&, int, int, int, bool, const DV2<float>&);
template void segments_to_particles<double>(const DV2<double>&, int, int, const DV2<double>&);
template void segments_to_particles<float>(const DV2<float>&, int, int, const DV2<float>&);

}  // namespace pack
}  // namespace ofk_iv
