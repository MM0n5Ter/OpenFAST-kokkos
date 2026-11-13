#include "ofk_runtime.hpp"

#include <Kokkos_Core.hpp>
#include <cstdlib>
#include <mutex>

namespace {
bool g_we_initialized = false;

void finalize_if_ours() {
  // Persistent Views owned by units must be released before this point:
  // register their cleanup with Kokkos::push_finalize_hook.
  if (g_we_initialized && Kokkos::is_initialized() && !Kokkos::is_finalized())
    Kokkos::finalize();
}
}  // namespace

namespace ofk {
void ensure_initialized() {
  static std::once_flag once;
  std::call_once(once, [] {
    if (Kokkos::is_initialized() || Kokkos::is_finalized()) return;  // host owns Kokkos
    // Settings come from the environment: KOKKOS_NUM_THREADS / OMP_NUM_THREADS,
    // KOKKOS_DEVICE_ID, KOKKOS_TOOLS_LIBS, KOKKOS_DISABLE_WARNINGS.
    Kokkos::initialize();
    g_we_initialized = true;
    // atexit handlers run in reverse order of registration. This one is
    // registered after the CUDA runtime initialized, so it runs before the
    // CUDA runtime tears down.
    std::atexit(finalize_if_ours);
  });
}
}  // namespace ofk

extern "C" void ofk_finalize() { finalize_if_ours(); }
