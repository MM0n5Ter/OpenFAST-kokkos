// P4 / FMA experiment of unit circulation-solve: time per cache build and per replay (Vvar) through the
// C API, serial loop vs Kokkos host parallel_for, for strict, F0 and F1, at lifting-line sizes from the
// r-test cases (a few control points) to nCP = 300. Synthetic rotor geometry (3 blades, iNWStart = 2,
// Vatistas by default), BAR-like spacing. Prints one CSV row per (nCP, variant, path):
//   csv,<nCP>,<variant>,<serial|parfor>,<build us>,<replay us>,<threads>
// Times: best of `reps` repetitions of the mean over `calls` calls.
//   perf_cs [reps=5] [reg=3] [nCP list...]
#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern "C" {
void ofk_cs_begin(const void*, int, int, int, int, int, const void*, const void*, const void*, const void*);
void ofk_cs_wing(const void*, int, int, int, int, int, int, const void*, int, int, int, const void*, int, int,
                 const void*);
void ofk_cs_build(const void*, int);
void ofk_cs_vvar(const void*, int, const void*, void*);
void ofk_cs_end(const void*, int, int);
void ofk_cs_set_variant(int, int, int);
void ofk_set_mode(int);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    const int reps = argc > 1 ? std::atoi(argv[1]) : 5;
    const int reg = argc > 2 ? std::atoi(argv[2]) : 3;
    std::vector<int> sizes;
    for (int i = 3; i < argc; ++i) sizes.push_back(std::atoi(argv[i]));
    if (sizes.empty()) sizes = {3, 6, 9, 12, 18, 24, 36, 48, 60, 87, 120, 150, 300};
    const int threads = Kokkos::DefaultHostExecutionSpace().concurrency();
    static int keyobj;
    const void* key = &keyobj;
    const double prec = 2.220446049250313e-16 / 100, fpi = 0.25 / std::acos(-1.0), minexp = -10.0, lamb = 1.25643;
    for (int nCP : sizes) {
      const int nW = 3, D = 2;
      std::vector<int> nS(nW, nCP / nW);
      for (int i = 0; i < nCP % nW; ++i) nS[i]++;
      std::vector<std::vector<double>> r(nW), e(nW), cp(nW);
      for (int w = 0; w < nW; ++w) {
        const int n = nS[w];
        const double th = 2 * std::acos(-1.0) * w / nW;
        r[w].assign(3 * (n + 1) * (D + 1), 0.0);
        e[w].assign(3 * n * D, 0.0);
        cp[w].assign(3 * std::max(n, 1), 0.0);
        for (int d = 0; d <= D; ++d)
          for (int s = 0; s <= n; ++s) {
            const double rr = 1.0 + 8.0 * s / std::max(n, 1);
            double* P = &r[w][3 * (s + (n + 1) * d)];
            P[0] = 0.25 + 0.3 * d; P[1] = rr * std::cos(th); P[2] = rr * std::sin(th) + 0.01 * d;
          }
        for (int s = 0; s < n; ++s) {
          for (int d = 0; d < D; ++d) e[w][3 * (s + n * d)] = 0.2;
          const double rr = 1.0 + 8.0 * (s + 0.5) / std::max(n, 1);
          cp[w][3 * s] = 0.05; cp[w][3 * s + 1] = rr * std::cos(th); cp[w][3 * s + 2] = rr * std::sin(th);
        }
      }
      std::vector<double> gam(nCP), V(3 * nCP);
      for (int i = 0; i < nCP; ++i) gam[i] = 1.0 + 0.01 * i;
      struct Var { const char* name; int mode, fma; } vars[] = {{"strict", 0, -1}, {"F0", 1, 0}, {"F1", 1, 1}};
      for (const Var& v : vars) {
        for (int par = 0; par <= 1; ++par) {
          ofk_cs_set_variant(v.mode, v.fma, par);
          const int calls = std::max(20, 2000000 / std::max(1, nCP * nCP * 8));
          double bestB = 1e30, bestV = 1e30;
          auto setup = [&] {
            ofk_cs_begin(key, 8, nW, nCP, D, reg, &prec, &fpi, &minexp, &lamb);
            int k0 = 0;
            for (int w = 0; w < nW; ++w) {
              ofk_cs_wing(key, 8, k0, nS[w], 3, nS[w] + 1, D + 1, r[w].data(), 3, nS[w], D, e[w].data(), 3,
                          std::max(nS[w], 1), cp[w].data());
              k0 += nS[w];
            }
          };
          for (int rep = 0; rep < reps; ++rep) {
            auto t0 = std::chrono::steady_clock::now();
            for (int c = 0; c < calls; ++c) { setup(); ofk_cs_build(key, 8); ofk_cs_end(key, 0, 30); }
            auto t1 = std::chrono::steady_clock::now();
            setup(); ofk_cs_build(key, 8);
            auto t2 = std::chrono::steady_clock::now();
            for (int c = 0; c < calls; ++c) ofk_cs_vvar(key, 8, gam.data(), V.data());
            auto t3 = std::chrono::steady_clock::now();
            ofk_cs_end(key, calls, 30);
            bestB = std::min(bestB, std::chrono::duration<double, std::micro>(t1 - t0).count() / calls);
            bestV = std::min(bestV, std::chrono::duration<double, std::micro>(t3 - t2).count() / calls);
          }
          std::printf("csv,%d,%s,%s,%.3f,%.3f,%d\n", nCP, v.name, par ? "parfor" : "serial", bestB, bestV, threads);
          std::fflush(stdout);
        }
      }
    }
  }
  Kokkos::finalize();
  return 0;
}
