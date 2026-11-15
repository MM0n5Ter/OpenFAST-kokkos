// Unit biotsavart-direct: C entry points behind FVW_BiotSavart::ui_seg and
// FVW_BiotSavart::ui_part_nograd (modules/aerodyn/src/FVW_BiotSavart.f90 on the pin 12fa3fc:
// :135 and :360). Called only from the Fortran shim in modules/aerodyn/src/FVW_BiotSavart.f90.
//
// Conventions (OpenFAST-kokkos notes/build-kokkos-fortran.md §6):
//  - every 2-D array is Fortran column-major A(ld, ncol); `ld` is size(A,1) of the actual
//    argument, so (3,n), (4,n) and padded arrays all work;
//  - ui_seg: iCPStart..iCPEnd, iSegStart..iSegEnd are the Fortran 1-based ranges, and the
//    values in SegConnct are 1-based column indices into SegPoints (nPts = size(SegPoints,2));
//  - Uind / UIout are accumulated into, and only columns of the given range are written;
//  - `minnorm` is FVW_BiotSavart::MINNORM as the Fortran compiler evaluated it (it is a
//    default-kind literal, so its value depends on -fdefault-real-8; passing it keeps the
//    C++ side identical to the Fortran side under any flags);
//  - _r8 = real(c_double) (ReKi with OPENFAST_DOUBLE_PRECISION), _r4 = real(c_float).
// The arithmetic mode (strict | fast) is chosen at run time by OFK_MODE (ofk_runtime.hpp).
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void ofk_biotsavart_direct_ui_seg_r8(int iCPStart, int iCPEnd, int ldCP, const double* CPs,
                                     int iSegStart, int iSegEnd, int ldPts, int nPts,
                                     const double* SegPoints, int ldConnct, const int* SegConnct,
                                     const double* SegGamma, int RegFunction,
                                     const double* RegParam, int ldU, double* Uind);
void ofk_biotsavart_direct_ui_seg_r4(int iCPStart, int iCPEnd, int ldCP, const float* CPs,
                                     int iSegStart, int iSegEnd, int ldPts, int nPts,
                                     const float* SegPoints, int ldConnct, const int* SegConnct,
                                     const float* SegGamma, int RegFunction,
                                     const float* RegParam, int ldU, float* Uind);

void ofk_biotsavart_direct_ui_part_nograd_r8(int nCPs, int ldCP, const double* CPs, int nPart,
                                             int ldP, const double* Part, int ldA,
                                             const double* Alpha, int RegFunction,
                                             const double* RegParam, int ldU, double* UIout,
                                             double minnorm);
void ofk_biotsavart_direct_ui_part_nograd_r4(int nCPs, int ldCP, const float* CPs, int nPart,
                                             int ldP, const float* Part, int ldA,
                                             const float* Alpha, int RegFunction,
                                             const float* RegParam, int ldU, float* UIout,
                                             float minnorm);

#ifdef __cplusplus
}
#endif
