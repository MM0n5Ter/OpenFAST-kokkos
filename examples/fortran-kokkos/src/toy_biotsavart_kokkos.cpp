// C++/Kokkos implementation of toy_biotsavart::ui_part (see ../ref for the Fortran).
//
// Boundary rules applied here (see notes/build-kokkos-fortran.md):
//  - extern "C", scalars by value, arrays as raw pointers from Fortran;
//  - Fortran column-major (3,n) -> Kokkos LayoutLeft, extents (3,n), 0-based;
//  - only columns [0,nCPs) of Uind are read or written;
//  - Uind is accumulated into, not overwritten;
//  - per-CP local accumulation, then one add (same order as the Fortran);
//  - no C++ exception crosses into Fortran.
#include <Kokkos_Core.hpp>
#include <cstdio>
#include <cstdlib>
#include <exception>

#include "ofk_runtime.hpp"

namespace {
template <class T>
using HostUnmanaged =
    Kokkos::View<T, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

constexpr double kFourPiInv = 0.25 / Kokkos::numbers::pi;
constexpr double kMinNorm2 = 1.0e-8;

void ui_part(int nCPs, const double* cps, int nPart, const double* part, const double* alpha,
             const double* regParam, double* uind) {
  using Mem = Kokkos::DefaultExecutionSpace::memory_space;

  HostUnmanaged<const double**> hCP(cps, 3, nCPs);
  HostUnmanaged<const double**> hP(part, 3, nPart);
  HostUnmanaged<const double**> hA(alpha, 3, nPart);
  HostUnmanaged<const double*> hR(regParam, nPart);
  HostUnmanaged<double**> hU(uind, 3, nCPs);

  // On a host backend these return the same (unmanaged) Views and copy nothing.
  // On CUDA they allocate device Views and copy. A real unit keeps device
  // buffers alive between calls (allocate in *_Init, free in *_End).
  auto CP = Kokkos::create_mirror_view_and_copy(Mem{}, hCP);
  auto P = Kokkos::create_mirror_view_and_copy(Mem{}, hP);
  auto A = Kokkos::create_mirror_view_and_copy(Mem{}, hA);
  auto R = Kokkos::create_mirror_view_and_copy(Mem{}, hR);
  auto U = Kokkos::create_mirror_view_and_copy(Mem{}, hU);  // copy in: we accumulate

  Kokkos::parallel_for(
      "ofk_toy_ui_part", Kokkos::RangePolicy<>(0, nCPs), KOKKOS_LAMBDA(const int i) {
        double ux = 0, uy = 0, uz = 0;
        for (int p = 0; p < nPart; ++p) {
          const double dx = CP(0, i) - P(0, p);
          const double dy = CP(1, i) - P(1, p);
          const double dz = CP(2, i) - P(2, p);
          const double r2 = dx * dx + dy * dy + dz * dz;
          if (r2 < kMinNorm2) continue;
          const double r = Kokkos::sqrt(r2);
          const double q = r / R(p);
          const double sp = kFourPiInv / (r2 * r) * (1.0 - Kokkos::exp(-q * q * q));
          ux += sp * (A(1, p) * dz - A(2, p) * dy);
          uy += sp * (A(2, p) * dx - A(0, p) * dz);
          uz += sp * (A(0, p) * dy - A(1, p) * dx);
        }
        U(0, i) += ux;
        U(1, i) += uy;
        U(2, i) += uz;
      });
  Kokkos::deep_copy(hU, U);  // fences; a no-op when U is hU
}
}  // namespace

extern "C" void ofk_toy_ui_part(int nCPs, const double* cps, int nPart, const double* part,
                                const double* alpha, const double* regParam, double* uind) {
  try {
    ofk::ensure_initialized();
    ui_part(nCPs, cps, nPart, part, alpha, regParam, uind);
  } catch (const std::exception& e) {
    // This routine has no ErrStat argument in Fortran, so a failure is fatal.
    // Units whose Fortran interface has ErrStat/ErrMsg return them instead.
    std::fprintf(stderr, "ofk_toy_ui_part: %s\n", e.what());
    std::abort();
  }
}
