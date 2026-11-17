/* itertrace: per-solve iteration counts of OLAF's circulation solve in an UNMODIFIED executable (the pin,
 * the pin built with FMA, or a port build), for the FMA experiment of unit circulation-solve.
 * No gdb or perf on this cluster, so this is a minimal ptrace breakpoint counter:
 *   - an int3 at the entry of Wings_ComputeCirculationPolarData starts a new solve;
 *   - an int3 at the call instruction of CirculationFromPolarData INSIDE the do-while loop
 *     (FVW_Wings.f90:458 on the pin) counts one iteration.
 * Each hit: restore the byte, step one instruction, re-insert, continue. The traced program's arithmetic
 * is untouched (its outputs are compared bitwise against an untraced run by the job script).
 * Threads (Kokkos-OpenMP) are followed (PTRACE_O_TRACECLONE); the breakpoints are in serial code.
 *
 *   itertrace <out file> <solve entry addr hex> <loop call addr hex> -- <exe> [args...]
 * out: one line per solve "<solve> <iterations>". Addresses from objdump (non-PIE executables only).
 * Build: gcc -O2 -o itertrace itertrace.c */
#define _GNU_SOURCE
#include <errno.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ptrace.h>
#include <sys/types.h>
#include <sys/user.h>
#include <sys/wait.h>
#include <unistd.h>

static long word0[2];
static uintptr_t addr[2];

static void arm(pid_t pid, int k) {
  long w = word0[k];
  long w3 = (w & ~0xffL) | 0xcc;
  if (ptrace(PTRACE_POKETEXT, pid, (void*)addr[k], (void*)w3) < 0) { perror("poke"); exit(3); }
}

int main(int argc, char** argv) {
  if (argc < 6 || strcmp(argv[4], "--")) {
    fprintf(stderr, "usage: itertrace <out> <solve addr> <loop call addr> -- exe [args]\n");
    return 2;
  }
  FILE* out = fopen(argv[1], "w");
  if (!out) { perror("out"); return 2; }
  addr[0] = strtoull(argv[2], NULL, 16);
  addr[1] = strtoull(argv[3], NULL, 16);
  pid_t child = fork();
  if (child == 0) {
    ptrace(PTRACE_TRACEME, 0, NULL, NULL);
    execvp(argv[5], argv + 5);
    perror("execvp");
    _exit(127);
  }
  int st;
  waitpid(child, &st, 0);  /* stopped at exec */
  ptrace(PTRACE_SETOPTIONS, child, NULL, (void*)(PTRACE_O_TRACECLONE | PTRACE_O_EXITKILL));
  for (int k = 0; k < 2; ++k) {
    errno = 0;
    word0[k] = ptrace(PTRACE_PEEKTEXT, child, (void*)addr[k], NULL);
    if (errno) { perror("peek"); return 3; }
    arm(child, k);
  }
  ptrace(PTRACE_CONT, child, NULL, NULL);
  long nsolve = 0, iters = 0, total = 0;
  int rc = 0;
  for (;;) {
    pid_t pid = waitpid(-1, &st, __WALL);
    if (pid < 0) break;
    if (WIFEXITED(st) || WIFSIGNALED(st)) {
      if (pid == child) { rc = WIFEXITED(st) ? WEXITSTATUS(st) : 128 + WTERMSIG(st); break; }
      continue;
    }
    if (!WIFSTOPPED(st)) continue;
    int sig = WSTOPSIG(st);
    if (sig == SIGTRAP && (st >> 16) == PTRACE_EVENT_CLONE) { ptrace(PTRACE_CONT, pid, NULL, NULL); continue; }
    if (sig == SIGSTOP) { ptrace(PTRACE_CONT, pid, NULL, NULL); continue; }  /* new thread */
    if (sig == SIGTRAP) {
      struct user_regs_struct r;
      ptrace(PTRACE_GETREGS, pid, NULL, &r);
      int k = -1;
      for (int i = 0; i < 2; ++i) if (r.rip - 1 == addr[i]) k = i;
      if (k >= 0) {
        if (k == 0) {
          if (nsolve > 0) fprintf(out, "%ld %ld\n", nsolve, iters);
          nsolve++; iters = 0;
        } else {
          iters++; total++;
        }
        ptrace(PTRACE_POKETEXT, pid, (void*)addr[k], (void*)word0[k]);
        r.rip = addr[k];
        ptrace(PTRACE_SETREGS, pid, NULL, &r);
        ptrace(PTRACE_SINGLESTEP, pid, NULL, NULL);
        waitpid(pid, &st, __WALL);
        arm(pid, k);
        ptrace(PTRACE_CONT, pid, NULL, NULL);
        continue;
      }
      ptrace(PTRACE_CONT, pid, NULL, NULL);  /* not ours */
      continue;
    }
    ptrace(PTRACE_CONT, pid, NULL, (void*)(long)sig);  /* deliver other signals */
  }
  if (nsolve > 0) fprintf(out, "%ld %ld\n", nsolve, iters);
  fclose(out);
  fprintf(stderr, "itertrace: %ld solves, %ld iterations, exit %d\n", nsolve, total, rc);
  return rc;
}
