// Kokkos runtime shared by every translated unit (one per process).
#pragma once

namespace ofk {
// Initialize Kokkos on first use. Callers (OpenFAST drivers, the C bindings
// loaded from Python, FAST.Farm) are not edited, so nobody else can do it.
// Must be called from serial code, never inside an OpenMP parallel region.
void ensure_initialized();
}  // namespace ofk

extern "C" {
// Optional explicit finalize (e.g. from a unit's *_End). Safe to call twice.
void ofk_finalize();
}
