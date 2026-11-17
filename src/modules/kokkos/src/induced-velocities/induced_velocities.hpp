// Unit induced-velocities: C entry points behind the direct-method paths of
//   FVW_Subs::InducedVelocitiesAll_Init / _Calc / _End    (pin FVW_Subs.f90:1237, :1285, :1331)
//   FVW_Subs::LiftingLineInducedVelocities                 (pin FVW_Subs.f90:1488)
// and, through them, WakeInducedVelocities (:1379), InducedVelocitiesAll_OnGrid_Calc (:1133) and
// srcPnl_ExtVelocities_OnPanels (:2081). Called only from the Fortran shims in FVW_Subs.f90 / FVW.f90.
//
// Per-instance state (P6): a handle owning every device buffer (sources, targets, outputs,
// lattice mirror, scratch), grown by doubling and never reallocated per call. Handles live in a
// process-wide map keyed by `key` = the address of the instance's T_Sgmt (m%Sgmt): every entry
// point of the unit receives that object, the registry types are not edited, and a copy of the
// misc type (FVW_CopyMisc) is a new address, hence a new, independent handle. Nothing in a handle
// is trusted without validation: the lattice mirror is re-validated slab by slab (bitwise memcmp
// against a host shadow) on every call, so glue copies of x, rollback and re-evaluation at the
// same t, and key reuse after a free can never produce stale sources (tests/induced-velocities/unit.md).
//
// Conventions (notes/build-kokkos-fortran.md §6):
//  - `rk` = storage size of ReKi in bytes (8 or 4); real arrays are passed as void* of that kind;
//  - 2-D arrays are Fortran column-major with leading dimension `ld`;
//  - integer arrays are c_int; Connct values are 1-based;
//  - Uind is accumulated into, only columns 1..nCPs are touched;
//  - every exception is caught at this boundary and is fatal (as the Fortran STOP is).
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// 1 if the sources are built on the device from the lattice (OFK_IV_SOURCES=lattice, the default
// on device backends), 0 if they are copied from the Fortran-packed Sgmt/Part (=packed, the default
// on host backends, where that copy is free). The shim calls ofk_iv_lattice only when 1.
int ofk_iv_use_lattice(void);

// Lattice of wing iW (1-based, wings passed in order 1..nWings before ofk_iv_sources_lattice):
// r_NW(3,nSpanP,nDepNWP), Gamma_NW(nSpanG,nDepNWG), Eps_NW(3,nSpanG,nDepNWG), same for FW, with
// their ALLOCATED extents; nNWact/nFWact = m%nNW/m%nFW (depths 1..n+1 of r, 1..n of Gamma/Eps
// are synchronized).
void ofk_iv_lattice(const void* key, int rk, int iW, int nWings, int nNWact, int nSpanNWP,
                    int nDepNWP, int nSpanNWG, int nDepNWG, const void* rNW, const void* gNW,
                    const void* eNW, int nFWact, int nSpanFWP, int nDepFWP, int nSpanFWG,
                    int nDepFWG, const void* rFW, const void* gFW, const void* eFW);

// Build the source set on the device from the lattice (PackPanelsToSegments + mirror [+
// SegmentsToPart]). method = VelocityMethod (1 segments, 3 particles), reg = RegFunction of the
// segments, regPart = RegFunctionPart. nSegF/nPartF = the counts the Fortran pack produced
// (checked against the closed-form count: a mismatch is fatal).
void ofk_iv_sources_lattice(const void* key, int rk, int method, int reg, int regPart,
                            int partPerSeg, int iDepthStart, int nNW, int nFW, int mirror,
                            int fwShed, int lastNWShed, int nSegF, int nPartF);

// Build the source set from the Fortran-packed arrays (Sgmt and, for method 3, Part).
void ofk_iv_sources_packed(const void* key, int rk, int method, int reg, int nSeg, int ldPts,
                           int nPts, const void* Pts, int ldC, const int* Connct,
                           const void* Gamma, const void* Eps, int regPart, int nPart, int ldP,
                           const void* PartP, int ldA, const void* PartAlpha,
                           const void* PartReg);

// 1 if the handle holds a valid source set for (method, nSeg, nPart), else 0.
int ofk_iv_sources_valid(const void* key, int method, int nSeg, int nPart);

// Evaluate the current source set at CPs(:,1:nCPs), accumulating into Uind(1:3,1:nCPs).
// capCP = size(CPs,2) (capacity hint, P6). minnorm = FVW_BiotSavart::MINNORM as compiled
// (default-kind literal, audit item 1).
void ofk_iv_eval(const void* key, int rk, int nCPs, int capCP, int ldCP, const void* CPs, int ldU,
                 void* Uind, const void* minnorm);

// Capacity hints from the Fortran misc arrays, which are allocated once at their maximum:
// latTot = total elements of all wings' lattice arrays, capSeg = size(Sgmt%Gamma),
// capPart = size(Part%P,2) (0 if unallocated), nBlocks = 2*nWings. Called before each source
// build; it allocates only on the first call (or if a hint grows).
void ofk_iv_reserve(const void* key, int rk, long long latTot, long long capSeg, long long capPart,
                    int nBlocks);

// Unit treecode (VelocityMethod 2). After ofk_iv_sources_packed(method = 2, the particles of
// Part): build the flat tree (grow_tree_part). degenerate = the default-kind literal 1.0e-5 of
// grow_tree_part_substep (FVW_VortexTools.f90:814) as the Fortran compiled it, as a ReKi.
void ofk_iv_tree_build(const void* key, int rk, const void* degenerate);
// ui_tree_part: evaluate the tree at CPs(:,1:nCPs), accumulating into Uind. branchFactor,
// distanceDirect: ReKi scalars (TreeBranchFactor(iVel), DistanceDirect).
void ofk_iv_eval_tree(const void* key, int rk, int nCPs, int capCP, int ldCP, const void* CPs,
                      int ldU, void* Uind, const void* minnorm, const void* branchFactor,
                      const void* distanceDirect);
// 1 unless OFK_TC=fortran (then the shims call the original Fortran treecode).
int ofk_tc_enabled(void);

// InducedVelocitiesAll_End: the source set is no longer valid (buffers are kept).
void ofk_iv_release_sources(const void* key);
// FVW_End: free the instance's handle. Safe if there is none. Handles that are never freed
// (copies destroyed by registry code, *_End not called) are freed by a Kokkos finalize hook and
// bounded by OFK_IV_MAX_HANDLES (least recently used is evicted).
void ofk_iv_free(const void* key);

// Test and profiling hooks.
//  counters: [0] device allocations by the handle, [1] bytes allocated, [2] source builds,
//  [3] evaluations, [4] bytes H2D, [5] bytes D2H, [6] lattice bytes compared, [7] lattice bytes
//  uploaded (0 on host backends), [8] number of live handles (process). Missing handle: zeros
//  except [8].
void ofk_iv_counters(const void* key, long long* out9);
// 0 = auto (default), 1 = packed, 2 = lattice; anything else restores OFK_IV_SOURCES.
void ofk_iv_set_sources_mode(int m);
// Print the per-phase profile of every live handle to stderr (also done at free/finalize when
// OFK_IV_PROFILE=1).
void ofk_iv_report(void);

#ifdef __cplusplus
}
#endif
