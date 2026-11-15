// Unit treecode, FAST mode (the default): the same tree, FMA contraction on, vectorized leaf sums
// on host backends (fixed lanes), and the device build (fast only). Class R2 against the direct sum.
// Same source as treecode_strict.cpp: only the two macros below and the flags differ.
#define OFK_BSD_MODE_NS fast
#define OFK_BSD_STRICT 0
#include "treecode_impl.hpp"
