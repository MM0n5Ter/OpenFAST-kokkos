// Unit biotsavart-direct, STRICT mode: the Fortran's per-target serial order, FMA contraction
// off (this file is compiled with -ffp-contract=off, and --fmad=false for CUDA device code;
// see ../CMakeLists.txt). Bitwise against the OPENMP=OFF Fortran on Serial/OpenMP.
// Same source as biotsavart_direct_fast.cpp: only the two macros below and the flags differ.
#define OFK_BSD_MODE_NS strict
#define OFK_BSD_STRICT 1
#include "biotsavart_direct_impl.hpp"
