// Kokkos runtime and process-wide switches shared by every translated unit (one per process).
// See notes/architecture.md §3 (modes) and §5 (lifetime) in the OpenFAST-kokkos repo.
#pragma once

namespace ofk {
// Initialize Kokkos on first use. Callers (OpenFAST drivers, the C bindings
// loaded from Python, FAST.Farm) are not edited, so nobody else can do it.
// Must be called from serial code, never inside an OpenMP parallel region.
void ensure_initialized();

// Reproducibility mode (architecture.md §3), process-wide, read once from OFK_MODE:
//   strict: the Fortran's per-target serial order, no FMA (bitwise vs the Fortran on CPU);
//   fast  : fixed blocked order (constant number of partial sums), FMA on. The default.
// The selected mode is printed once to stdout on first use.
enum class Mode : int { strict = 0, fast = 1 };
Mode mode();
const char* mode_name(Mode m);
}  // namespace ofk

extern "C" {
// Optional explicit finalize (e.g. from a unit's *_End). Safe to call twice.
void ofk_finalize();
// Mode override for tests (0 = strict, 1 = fast; anything else restores OFK_MODE).
// Not thread-safe with respect to running kernels: call it between kernel calls.
void ofk_set_mode(int m);
int ofk_get_mode();
}
