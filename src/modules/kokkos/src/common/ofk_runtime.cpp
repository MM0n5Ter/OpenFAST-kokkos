#include "ofk_runtime.hpp"

#include <Kokkos_Core.hpp>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>

namespace {
bool g_we_initialized = false;

void finalize_if_ours() {
  // Persistent Views owned by units must be released before this point:
  // register their cleanup with Kokkos::push_finalize_hook.
  if (g_we_initialized && Kokkos::is_initialized() && !Kokkos::is_finalized())
    Kokkos::finalize();
}

// -1: not decided yet (read OFK_MODE), otherwise an ofk::Mode value.
std::atomic<int> g_mode{-1};

int mode_from_env() {
  const char* v = std::getenv("OFK_MODE");
  int m = static_cast<int>(ofk::Mode::fast);  // default (architecture.md §3)
  bool bad = false;
  if (v && *v) {
    if (!std::strcmp(v, "strict")) m = static_cast<int>(ofk::Mode::strict);
    else if (!std::strcmp(v, "fast")) m = static_cast<int>(ofk::Mode::fast);
    else bad = true;
  }
  std::printf("[ofk] OFK_MODE=%s%s%s%s\n", ofk::mode_name(static_cast<ofk::Mode>(m)),
              (v && *v) ? "" : " (default; set OFK_MODE=strict|fast)", bad ? " (ignored invalid value '" : "",
              bad ? v : "");
  if (bad) std::printf("')\n");
  std::fflush(stdout);
  return m;
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

Mode mode() {
  int m = g_mode.load(std::memory_order_acquire);
  if (m < 0) {
    static std::once_flag once;
    static int from_env = 1;
    std::call_once(once, [] { from_env = mode_from_env(); });
    int expected = -1;
    g_mode.compare_exchange_strong(expected, from_env);
    m = g_mode.load(std::memory_order_acquire);
  }
  return static_cast<Mode>(m);
}

const char* mode_name(Mode m) { return m == Mode::strict ? "strict" : "fast"; }
}  // namespace ofk

extern "C" void ofk_finalize() { finalize_if_ours(); }

extern "C" void ofk_set_mode(int m) {
  if (m == 0 || m == 1) {
    (void)ofk::mode();  // print the environment's choice once, so logs stay comparable
    g_mode.store(m, std::memory_order_release);
  } else {
    g_mode.store(-1, std::memory_order_release);
  }
}

extern "C" int ofk_get_mode() { return static_cast<int>(ofk::mode()); }
