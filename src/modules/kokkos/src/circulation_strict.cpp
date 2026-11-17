// Unit circulation-solve, STRICT mode: baseline ISA, FMA contraction off (-ffp-contract=off). Bitwise against the pin (OPENMP=OFF) on the host.
// Same source as the other circulation_*.cpp: only the namespace macro and the per-file flags
// (ofk_kernels.cmake) differ. Selection at run time: circulation.cpp (OFK_CS_MODE / OFK_MODE, OFK_CS_FMA).
#define OFK_CS_NS strict
#include "circulation_body.hpp"
