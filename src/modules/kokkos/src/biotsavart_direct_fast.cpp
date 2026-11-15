// Unit biotsavart-direct, FAST mode (the default): fixed blocked order (kLanes partial sums,
// independent of the thread count), FMA contraction on (this file is compiled with
// -ffp-contract=fast plus OFK_FAST_HOST_FLAGS, e.g. -march=x86-64-v3, so FMA instructions
// exist; CUDA device code keeps nvcc's default --fmad=true; see ../CMakeLists.txt).
// Same source as biotsavart_direct_strict.cpp: only the two macros below and the flags differ.
#define OFK_BSD_MODE_NS fast
#define OFK_BSD_STRICT 0
#include "biotsavart_direct_impl.hpp"
