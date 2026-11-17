/* L0 of unit circulation-solve: count the calls of AirfoilInfo::AFI_ComputeAirfoilCoefs (linked with
 * -Wl,--wrap=__airfoilinfo_MOD_afi_computeairfoilcoefs). CirculationFromPolarData calls it once per
 * control point (FVW_Wings.f90:573 on the pin), so a solve's iteration count is
 * calls / nCP_tot - (m%FirstCall ? 1 : 0), for the oracle and the port alike, without touching either.
 * Fortran interface: AFI_ComputeAirfoilCoefs(AOA, Re, UserProp, p, AFI_interp, errStat, errMsg) with the
 * hidden length of errMsg last (gfortran ABI). */
#include <stddef.h>

long long ofk_afi_calls = 0;

void __real___airfoilinfo_MOD_afi_computeairfoilcoefs(void*, void*, void*, void*, void*, void*, char*, size_t);

void __wrap___airfoilinfo_MOD_afi_computeairfoilcoefs(void* aoa, void* re, void* up, void* p, void* out, void* es,
                                                      char* em, size_t len) {
  ++ofk_afi_calls;
  __real___airfoilinfo_MOD_afi_computeairfoilcoefs(aoa, re, up, p, out, es, em, len);
}
