// Unit biotsavart-direct: types shared by the entry TU (biotsavart_direct.cpp) and the two
// mode TUs (biotsavart_direct_strict.cpp, biotsavart_direct_fast.cpp). Data only: nothing in
// this header does floating-point arithmetic that FMA contraction could change, so it is safe
// to include from TUs compiled with different -ffp-contract / --fmad flags.
#pragma once

#include <Kokkos_Core.hpp>

#include <type_traits>

namespace ofk_bsd {

// FVW_BiotSavart.f90:22-28
constexpr int kRegNone = 0;
constexpr int kRegRankine = 1;
constexpr int kRegLambOseen = 2;
constexpr int kRegVatistas = 3;
constexpr int kRegOffset = 4;
constexpr int kRegPartExp = 1;
constexpr int kRegPartCompact = 2;

// Module constants (FVW_BiotSavart.f90:12-33). Built on the host by make_consts (entry TU),
// carried by value into the kernels.
template <class T>
struct Consts {
  T precision_ui;   // :12  epsilon(1.0_ReKi)/100
  T min_exp_value;  // :14  -10.0_ReKi
  T part_reg_cut3;  // :16  PART_REG_NRAD**3 = 2.0_ReKi**3
  T part_reg_c2;    // :17  1.6_ReKi
  T mindenom;       // :18  0.0_ReKi (kept: the Fortran adds it)
  T minnorm;        // :20  1e-4, DEFAULT-KIND literal: value comes from the Fortran side
  T minnorm2;       // :377 MINNORM*MINNORM
  T fourpi_inv;     // :32  0.25_ReKi/ACOS(-1.0_ReKi)
  T lamb_oseen;     // :107 1.25643_ReKi
};

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;
constexpr bool kOnHost = Kokkos::SpaceAccessibility<Kokkos::HostSpace, MemSpace>::accessible;

template <class T>
using DV2 = Kokkos::View<T**, Kokkos::LayoutLeft, MemSpace>;
template <class T>
using DV1 = Kokkos::View<T*, Kokkos::LayoutLeft, MemSpace>;

// Loop shape (architecture.md §2). The dispatcher (in the mode TU) resolves kAuto.
enum class Shape : int { kAuto = 0, kSquare = 1, kSkinny = 2 };

// Process-wide tuning switches, read once from the environment (entry TU). Not per-instance.
struct Tuning {
  Shape shape = Shape::kAuto;  // OFK_BSD_SHAPE = auto | square | skinny
  int skinny_max_cp = 1024;    // OFK_BSD_SKINNY_MAX_CP
  int skinny_min_seg = 4096;   // OFK_BSD_SKINNY_MIN_SEG
  int chunk = 1024;            // OFK_BSD_CHUNK: sources per team (skinny)
  int team = 128;              // OFK_BSD_TEAM: team size (skinny); fixed, so the order is fixed
};

// The process-wide tuning switches (OFK_BSD_*), read once (biotsavart_direct.cpp).
const Tuning& tuning();

// ui_seg on the device: arrays exactly as the Fortran passed them (AoS, 1-based SegConnct
// values, pOff = smallest referenced point index so P holds only the referenced block).
template <class T>
struct SegArgs {
  int nCP = 0, nSeg = 0, pOff = 1, reg = 0;
  DV2<const T> CP, P;
  DV2<const int> C;
  DV1<const T> G, R;  // R has extent 0 for kRegNone (not read)
  DV2<T> U;           // accumulated into
  Consts<T> k;
  Tuning tune;
};

template <class T>
struct PartArgs {
  int nCP = 0, nPart = 0, reg = 0;
  DV2<const T> CP, P, A;
  DV1<const T> R;  // extent 0 for kRegNone (not read)
  DV2<T> U;
  Consts<T> k;
  Tuning tune;
};

// Implemented once per mode (biotsavart_direct_impl.hpp, instantiated in _strict.cpp and
// _fast.cpp). Return the shape that ran.
// SoA column counts of the packed sources (biotsavart_direct_impl.hpp):
//   segments : x1 y1 z1 x2 y2 z2 gamma regparam
//   particles: px py pz alpha_x alpha_y alpha_z regparam
constexpr int kSegSoaCols = 8;
constexpr int kPartSoaCols = 7;

// run_seg / run_part: the stateless unit-1 entry (pack into a temporary SoA, then evaluate).
// run_seg_soa / run_part_soa: evaluate an already packed SoA block S of shape (>= n, cols)
// (rows 0..n-1 used; the row capacity may be larger) at a.CP, accumulating into a.U;
// `scratch` is a persistent buffer for the skinny shape's partial sums (may be empty).
// Used by unit induced-velocities, whose per-instance handle owns S, scratch, CP and U (P6).
namespace strict {
template <class T> Shape run_seg(const SegArgs<T>& a);
template <class T> Shape run_part(const PartArgs<T>& a);
template <class T> Shape run_seg_soa(const SegArgs<T>& a, const DV2<T>& S, const DV1<T>& scratch);
template <class T> Shape run_part_soa(const PartArgs<T>& a, const DV2<T>& S);
}  // namespace strict
namespace fast {
template <class T> Shape run_seg(const SegArgs<T>& a);
template <class T> Shape run_part(const PartArgs<T>& a);
template <class T> Shape run_seg_soa(const SegArgs<T>& a, const DV2<T>& S, const DV1<T>& scratch);
template <class T> Shape run_part_soa(const PartArgs<T>& a, const DV2<T>& S);
}  // namespace fast

}  // namespace ofk_bsd
