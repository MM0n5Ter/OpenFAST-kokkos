// Unit circulation-solve (OLAF unit 4a, F10-B): C API called by the OPENFAST_KOKKOS shim in
// FVW_Wings.f90 (Wings_ComputeCirculationPolarData). Design and gates: tests/circulation-solve/unit.md.
//
// Per solve (one call of Wings_ComputeCirculationPolarData):
//   ofk_cs_begin  sizes, regularization and the ui_seg_11 constants (from the shim, criteria R3)
//   ofk_cs_wing   per wing: its control points (targets) and its first iNWStart near-wake rings
//   ofk_cs_build  the Gamma-independent geometry cache, once
//   ofk_cs_vvar   per iteration: Vvar(1:3,1:nCP_tot) from GammaLastIter(1:nCP_tot)
//   ofk_cs_end    iteration count (OFK_CS_ITERS log, counters)
// Arrays are the Fortran's own (by reference, column-major, no copy of Gamma or Vvar); rk = bytes
// of ReKi (8 or 4). key = the instance's m%Sgmt (as unit induced-velocities). Nothing is trusted
// across calls: every solve re-reads its inputs and rebuilds the cache (criteria A7 by
// construction). Only the Views' capacity persists (P6).
//
// Environment (read once, process-wide):
//   OFK_CS=kokkos|fortran     kokkos (default): this unit; fortran: the pin's ui_quad_n1 loop (A/B)
//   OFK_CS_MODE=strict|fast   this unit's mode (default: OFK_MODE)
//   OFK_CS_FMA=0|1            fast mode only: 0 = TU circulation_fast (no contraction, F0: bitwise),
//                             1 = TU circulation_fastfma (contraction, F1). Default 1 (FMA experiment)
//   OFK_CS_PAR_MIN=<pairs>    targets x segments from which the replay runs as a parallel_for (P4)
//   OFK_CS_PAR_MIN_BUILD=<pairs>  the same for the cache build
//   OFK_CS_ITERS=<file>       one line per solve: "<solve> <iterations> <cap>" (FMA experiment, C4)
//   OFK_CS_PROFILE=1          per-handle counters and times on stderr at free / finalize
//   OFK_CS_MAX_HANDLES=<n>    live-handle cap (default 64), least recently used evicted
#pragma once

extern "C" {
int ofk_cs_enabled(void);
void ofk_cs_begin(const void* key, int rk, int nWings, int nCPtot, int iNWStart, int reg, const void* precision_ui,
                  const void* fourpi_inv, const void* min_exp_value, const void* lamb_oseen);
// kCP0: 0-based offset of the wing's first control point in the flattened target list.
// r_NW(ldR1, ldR2, ldR3), Eps_NW(ldE1, ldE2, ldE3), CP(ldC1, ldC2): the Fortran extents.
void ofk_cs_wing(const void* key, int rk, int kCP0, int nSpan, int ldR1, int ldR2, int ldR3, const void* rNW,
                 int ldE1, int ldE2, int ldE3, const void* eNW, int ldC1, int ldC2, const void* CP);
void ofk_cs_build(const void* key, int rk);
void ofk_cs_vvar(const void* key, int rk, const void* gam, void* V);
void ofk_cs_end(const void* key, int nIter, int maxIter);
void ofk_cs_free(const void* key);
// Test hook (L0). mode: -1 = environment, 0 = strict, 1 = fast (this unit only; ofk_set_mode sets the
// process mode); fma: -1 = environment, 0/1 = F0/F1; par: -1 = threshold, 0 = serial, 1 = parallel_for.
void ofk_cs_set_variant(int mode, int fma, int par);
// out: allocations, bytes allocated, solves, iterations, builds, replays, live handles, last solve parallel (0/1)
void ofk_cs_counters(const void* key, long long* out);
}
