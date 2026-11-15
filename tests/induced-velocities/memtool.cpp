// Minimal Kokkos Tools library for P6 (criteria.md): logs every Kokkos allocation and
// deallocation with its memory space, label and size, and the number of "ofk_iv_eval" profiling
// regions (one per induced-velocity evaluation of unit induced-velocities) seen before it.
// Load with KOKKOS_TOOLS_LIBS=<path>/libofk_memtool.so; OFK_MEMTOOL_OUT=<file> (default stderr).
// The summary line at finalize gives the last evaluation index after which anything was allocated.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
FILE* out = nullptr;
long long nEval = 0, nAlloc = 0, nFree = 0, lastAllocEval = -1, bytesAlloc = 0;
}  // namespace

struct SpaceHandle { char name[64]; };

extern "C" void kokkosp_init_library(const int, const uint64_t, const uint32_t, void*) {
  const char* f = std::getenv("OFK_MEMTOOL_OUT");
  out = (f && *f) ? std::fopen(f, "w") : stderr;
  if (!out) out = stderr;
  std::fprintf(out, "# ofk_memtool: event space label bytes evals_before\n");
}
extern "C" void kokkosp_finalize_library() {
  std::fprintf(out, "# summary: allocations=%lld (%.3f MB) deallocations=%lld evals=%lld last_alloc_after_eval=%lld\n",
               nAlloc, bytesAlloc / 1e6, nFree, nEval, lastAllocEval);
  if (out != stderr) std::fclose(out);
}
extern "C" void kokkosp_allocate_data(const SpaceHandle h, const char* label, const void*, const uint64_t size) {
  ++nAlloc; bytesAlloc += (long long)size; lastAllocEval = nEval;
  std::fprintf(out, "alloc %s %s %llu %lld\n", h.name, label, (unsigned long long)size, nEval);
}
extern "C" void kokkosp_deallocate_data(const SpaceHandle h, const char* label, const void*, const uint64_t size) {
  ++nFree;
  std::fprintf(out, "free %s %s %llu %lld\n", h.name, label, (unsigned long long)size, nEval);
}
extern "C" void kokkosp_push_profile_region(const char* name) {
  if (!std::strcmp(name, "ofk_iv_eval")) ++nEval;
}
extern "C" void kokkosp_pop_profile_region() {}
