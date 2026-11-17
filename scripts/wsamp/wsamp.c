// Wall-clock stack sampler for the main thread (LD_PRELOAD). Unit 4 analysis, no source changes needed.
// Env: WSAMP_OUT=<file prefix> (required to enable), WSAMP_HZ (default 1000), WSAMP_DEPTH (max 96).
// Output: <prefix>.stacks (binary: per sample int32 depth + depth*uint64 PCs), <prefix>.maps (copy of /proc/self/maps).
#define _GNU_SOURCE
#include <execinfo.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/syscall.h>
#ifndef sigev_notify_thread_id
#define sigev_notify_thread_id _sigev_un._tid
#endif
#define MAXD 96
static void **buf; static int32_t *dep; static size_t cap, n; static int depth = 64;
static timer_t tid; static int on; static char prefix[1024];
static void h(int s, siginfo_t *si, void *uc) {
  (void)s; (void)si; (void)uc;
  if (n >= cap) return;
  int d = backtrace(buf + n * MAXD, depth);
  dep[n++] = d;
}
__attribute__((constructor)) static void init(void) {
  const char *o = getenv("WSAMP_OUT"); if (!o) return;
  snprintf(prefix, sizeof prefix, "%s", o);
  const char *hz = getenv("WSAMP_HZ"); double f = hz ? atof(hz) : 1000.0;
  const char *dd = getenv("WSAMP_DEPTH"); if (dd) { depth = atoi(dd); if (depth > MAXD) depth = MAXD; }
  cap = 2000000; buf = calloc(cap * MAXD, sizeof(void*)); dep = calloc(cap, sizeof(int32_t));
  void *tmp[4]; backtrace(tmp, 4);               // load libgcc's unwinder outside the handler
  struct sigaction sa; memset(&sa, 0, sizeof sa); sa.sa_sigaction = h; sa.sa_flags = SA_SIGINFO | SA_RESTART;
  sigaction(SIGPROF, &sa, 0);
  struct sigevent se; memset(&se, 0, sizeof se);
  se.sigev_notify = SIGEV_THREAD_ID; se.sigev_signo = SIGPROF; se.sigev_notify_thread_id = (pid_t)syscall(SYS_gettid);
  if (timer_create(CLOCK_MONOTONIC, &se, &tid)) { perror("wsamp timer_create"); return; }
  long ns = (long)(1e9 / f); struct itimerspec its = {{ns / 1000000000, ns % 1000000000}, {ns / 1000000000, ns % 1000000000}};
  timer_settime(tid, 0, &its, 0); on = 1;
}
__attribute__((destructor)) static void fini(void) {
  if (!on) return; on = 0; timer_delete(tid);
  char fn[1100]; snprintf(fn, sizeof fn, "%s.stacks", prefix);
  FILE *f = fopen(fn, "wb");
  for (size_t i = 0; i < n; i++) { fwrite(&dep[i], 4, 1, f); for (int j = 0; j < dep[i]; j++) { uint64_t a = (uint64_t)buf[i * MAXD + j]; fwrite(&a, 8, 1, f); } }
  fclose(f);
  snprintf(fn, sizeof fn, "%s.maps", prefix);
  FILE *m = fopen("/proc/self/maps", "r"), *w = fopen(fn, "w"); char line[4096];
  while (fgets(line, sizeof line, m)) fputs(line, w); fclose(m); fclose(w);
  fprintf(stderr, "[wsamp] %zu samples -> %s.stacks\n", n, prefix);
}
