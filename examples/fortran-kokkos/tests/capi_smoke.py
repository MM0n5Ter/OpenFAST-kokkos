#!/usr/bin/env python3
"""API-gate smoke test for the shared library, the way pyOpenFAST loads
libaerodyn_inflow_c_binding.so: ctypes + the bind(C) symbol, no knowledge of
the port. Checks that the library loads (Kokkos initializes lazily inside it),
gives finite results, accumulates into Uind, and that the process exits
cleanly (Kokkos finalizes at exit).

Usage: capi_smoke.py path/to/libtoy_c_binding.so
"""
import ctypes, math, random, sys

lib = ctypes.CDLL(sys.argv[1])
f = lib.Toy_C_Induced
D = ctypes.POINTER(ctypes.c_double)
f.argtypes = [ctypes.c_int, D, ctypes.c_int, D, D, D, D]
f.restype = None

random.seed(0)
ncp, npart = 64, 128
arr = lambda n, lo, hi: (ctypes.c_double * n)(*[random.uniform(lo, hi) for _ in range(n)])
cps, part, alpha = arr(3 * ncp, -5, 5), arr(3 * npart, -5, 5), arr(3 * npart, -1, 1)
reg = arr(npart, 0.05, 1.05)
u1 = (ctypes.c_double * (3 * ncp))()
f(ncp, cps, npart, part, alpha, reg, u1)
u2 = (ctypes.c_double * (3 * ncp))(*u1)
f(ncp, cps, npart, part, alpha, reg, u2)          # second call adds the same amount again
assert all(math.isfinite(x) for x in u1), "non-finite output"
assert max(abs(x) for x in u1) > 0, "all-zero output"
assert all(abs(b - 2 * a) <= 1e-12 * max(1.0, abs(a)) for a, b in zip(u1, u2)), "Uind is not accumulated"
print("capi_smoke: OK")
