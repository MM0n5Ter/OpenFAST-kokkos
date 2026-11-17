// Unit circulation-solve, FAST mode, variant F1 of the FMA experiment: OFK_FAST_HOST_FLAGS (x86-64-v3), contraction ON (-ffp-contract=fast).
// Same source as the other circulation_*.cpp: only the namespace macro and the per-file flags
// (ofk_kernels.cmake) differ. Selection at run time: circulation.cpp (OFK_CS_MODE / OFK_MODE, OFK_CS_FMA).
#define OFK_CS_NS fastfma
#include "circulation_body.hpp"
