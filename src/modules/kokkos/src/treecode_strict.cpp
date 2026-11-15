// Unit treecode, STRICT mode: the Fortran tree, bit for bit on the CPU (same structure, preorder
// traversal, multipole sums and acceptance criterion), FMA contraction off (-ffp-contract=off,
// --fmad=false for CUDA device code; see ../CMakeLists.txt). Same source as treecode_fast.cpp.
#define OFK_BSD_MODE_NS strict
#define OFK_BSD_STRICT 1
#include "treecode_impl.hpp"
