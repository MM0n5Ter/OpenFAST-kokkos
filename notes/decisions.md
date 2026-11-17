# Decisions and story (drafted 2026-11-12)

## The story
- OpenFAST is the reference wind-turbine simulator (NREL/NLR). It is about 390k lines of hand-written Fortran plus 140k generated, and it runs on CPUs only.
- Its expensive mid-fidelity option, the **OLAF free-vortex wake**, takes days per production run (1.7–5 days, upstream PRs #3430/#3464). Its induction kernels act on 10⁵–10⁶ wake markers, which is GPU scale. Upstream says GPU OLAF "is not currently supported… NREL is not currently funded to work on that" (issue #2852).
- OpenFAST is a **library** as well as a program:
  - its glue code drives every module through one framework API;
  - AMR-Wind/ExaWind link it through `FAST_Library`;
  - pyOpenFAST loads its C bindings;
  - ROSCO plugs in through the controller DLL interface;
  - WEIS reads and writes its input and output files.

  A standalone rewrite would break all of these users. We migrate it **in place**, keep every interface, and gate every step.

## Decisions
1. **Pinned upstream: OpenFAST `dev` `12fa3fc` (2026-11-06), not the v5.0.0 release (`2895884`, 2026-03-12).**
   - Upstream rewrote exactly our target kernels on `dev` in Sep 2026 (#3430/#3457/#3464). They added a compact particle core, SIMD near-field sums and a new leaf layout, and made CPU OLAF about 3× faster. Porting v5.0.0 would compete against an old baseline.
   - r-test is pinned through OpenFAST's own submodule pointer (`reg_tests/r-test` → `dfb039d`).
   - Moving the pin is a deliberate commit that re-runs `scripts/` and updates the numbers.
2. **One-month goal: the OLAF induction unit. End state: the whole repository Fortran-free.** It is reached incrementally. The repository builds and passes its tests after every step. See [roadmap.md](roadmap.md).
3. **The unit of work is an API contract, not a directory.**
   - Every module exposes the same framework API (`*_Init / _UpdateStates / _CalcOutput / _End` with registry types). One contract template covers all of them.
   - Inside a module, smaller units have their own contracts. The first one is `InducedVelocitiesAll_Init/Calc/End` plus the direct kernels `ui_seg`/`ui_part_nograd`. It takes flat arrays, so it needs no registry types.
   - Details in [openfast-structure.md](openfast-structure.md).
4. **For each unit:**
   1. Freeze its interface as a contract.
   2. Rewrite the unit in C++/Kokkos.
   3. The Fortran module **keeps its name and public interface**. Its body becomes a thin `bind(C)` shim, so remaining Fortran callers compile unchanged. The shim is deleted when its last Fortran caller is gone.
   4. External APIs are kept permanently:
      - the C bindings (ADI, IfW, HD, SeaState, MoorDyn) and `FAST_Library`;
      - the Bladed-style controller DLL interface;
      - input files, output channels (`.out`/`.outb`) and VTK output.
   5. Registry-generated types (`*_Types.f90`) are not edited by hand. They change only by retargeting `openfast-registry` (long-term phase 1).
5. **Extend OpenFAST's existing CMake build.** We add an `OPENFAST_KOKKOS` option and one C++ library, `ofk_kernels`, that holds every ported unit, so the process has one Kokkos runtime. Never generate a new build system: ParEval-Repo finds build generation is the main reason repo-level translation fails.
   - Fortran: gfortran (gcc 12.3). OpenFAST's CMake accepts only GNU, Intel or Flang (`cmake/OpenfastFortranOptions.cmake:48-52`), so nvfortran is out.
   - C++: g++ for the Kokkos-OpenMP build, `nvcc_wrapper` for the Kokkos-CUDA build. See [build-kokkos-fortran.md](build-kokkos-fortran.md).
6. **Kokkos 5.3.0, built from source** (the cluster has no Kokkos module). It needs CMake ≥ 3.22, GCC ≥ 10.4, NVCC ≥ 12.2 and C++20.
   - **Built 2026-11-14** (gcc 12.3, cuda 12.3, cmake 3.31.12), shared libraries, C++20:
     - `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/kokkos-5.3.0-omp` (Serial + OpenMP);
     - `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/kokkos-5.3.0-cuda-h100` (CUDA sm_90 + OpenMP + Serial).
   - Recipe: `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/build-kokkos-5.3.0.sh omp|cuda-h100`. Source: the release tarball, sha256 `f1875cac…ef88`.
   - Code is C++20, and CMake uses `find_package(Kokkos 5.3...<6)`. The earlier 4.5.1 installs are superseded; do not use them.
   - CUDA code is **compiled on the login node** (nvcc from `module load cuda/12.3.0`). Only runs need a GPU job.
   - Two installs: OpenMP, and CUDA (one per GPU architecture).
   - The port initializes Kokkos lazily on first use, because its callers cannot be edited.
7. **Three gates per commit:**
   - **API:** callers compile, link and run with zero source edits, including the C bindings loaded from Python. Output channels and VTK structure are part of the API.
   - **Regression:** the unit's L0 differential tests; OpenFAST's unit tests; r-test at its own tolerances (rtol 1e-2, magnitude-scaled atol).
   - **Performance:** time per step against the number of wake markers, comparing Fortran+OpenMP, Kokkos-OpenMP and Kokkos-CUDA.
8. **Tests are not ported first.** r-test runs the binaries and compares outputs, so it doesn't depend on the source language. OpenFAST's Fortran unit tests (test-drive, `FVW_Tests.f90`) stay as API consumers of the ported unit. Before porting a unit we **add** tests. Details in [testing.md](testing.md) and [harness.md](harness.md).
9. **Four people work in four lanes.**
   - Lanes conflict only on shared files: CMake files, registry files, and module-level variables. Each shared file has exactly one owner, and all build changes go through Lane A.
   - One unit = one branch = one PR, merged only when all gates pass.
10. **We report agent cost and effort per unit, in money and time:**
    - tokens and dollars (USD API-equivalent), plus GPU/CPU hours;
    - agent wall and active hours, and human minutes;
    - E_k (from ParEval-Repo) and the number of human interventions.

    Every session appends a row to `tools/cost/ledger.csv` (criteria.md E1–E7, `tools/cost/README.md`).
12. **Reproducibility: strict and fast modes from one source; no CPU↔GPU bitwise requirement** (user, 2026-11-14).
    - Bitwise agreement with the Fortran is required only of the strict mode on the CPU, where it proves the translation is faithful.
    - The fast mode (the default) may reorder sums, use FMA and change data layout. It is held to bounded error against strict, to run-to-run determinism, and to R8 at simulation level.
    - Reasons:
      - GPU speed needs reorganized data and order;
      - a bitwise copy of the Fortran is not proof of correctness, because the Fortran's own results vary with compiler, flags and OpenMP.
    - Details: [architecture.md](architecture.md).
13. **Three units establish the architecture** (user, 2026-11-14):
    1. the direct kernels;
    2. `InducedVelocitiesAll_*` plus lifting-line induction, running end to end in the fork;
    3. the particle treecode.

    After them, a working architecture must exist that carries every high-level idea, even if it is not perfect. Each unit is ported by a fresh-context subagent.
14. **Unit 4 boundary: other modules are external library calls** (user, 2026-11-17). Unit 4b may change how often and where AeroDyn calls InflowWind. It does not edit, overlay or port InflowWind, so the IEEE fix (`IfW_FlowField.f90:386`) goes to InflowWind's owners as an upstream suggestion. The general rule: a unit's boundary stops at calls into a module it does not own. See `notes/olaf-unit4.md` §3.1.
15. **Unit 4a fast mode: no FMA by default** (user, 2026-11-17). The circulation solve's fast mode uses F0 (no contraction; bitwise with strict and the pin, iteration counts included). F1 (contraction) is an explicit opt-in, `OFK_CS_FMA=1`. The FMA experiment's rule picked F1 on a 4–6% kernel gain, but the whole run shows no measurable benefit and F1 exceeds the per-case claim-3 limit in 2 L0 cases. The limit stays unchanged, and the two cases are kept as reproducible L0 cases (`tests/circulation-solve/unit.md`). Revisit F1 if profiling later shows a meaningful end-to-end benefit; no further FMA experiments until then. An exception, for this unit only, to the open item "FMA in fast mode: keep it" below.
11. **Precision.** Target `DOUBLE_PRECISION=ON` (ReKi = 8 bytes, the default). The single-precision build must keep compiling, since upstream just added single-precision unit tests (#3511, the pinned commit). The C++ side picks its real type from the `OPENFAST_DOUBLE_PRECISION` macro that OpenFAST already defines for all languages.

## Open items / TODO
- [x] **Where does the port live? An overlay in this repo, not a fork** (user, 2026-11-14, replacing the same day's fork decision).
  - **Changed files** are full copies at OpenFAST's own relative paths under `src/modules/<mod>/src/`. **New C++** goes in `src/modules/kokkos/`.
  - `src/cmake/ofk_overlay.cmake` is injected into OpenFAST's unmodified CMake (`-DCMAKE_PROJECT_INCLUDE`). It swaps the overlay files into the module targets and links `ofk_kernels`.
  - `deps/openfast` stays the untouched reference and build source.
  - Reason: we only touch a few module files, not the whole of OpenFAST. The overlay shows exactly what changed.
  - The earlier local `fork/` submodule (never pushed) is being migrated; see `notes/overlay.md` when it exists. Its directory stays on disk until the queued GPU jobs that read it have run.
- [x] Build both Kokkos installs and run `examples/fortran-kokkos` on them. Done with 5.3.0 (decision #6). The CUDA run is on h100 rather than a100.
- [x] Compiler setup for the CUDA build: keep g++, and call `kokkos_compilation(GLOBAL)` in `fork/modules/kokkos` (unit 1). A global `nvcc_wrapper` fails on OpenFAST's `awae/src/vtk.cpp` (`fast_float.h` incomplete type).
- [x] **FMA in fast mode: keep it** (user, 2026-11-14). Where FMA moves results more than the class bound, the bound for chain claim 3 is **1.1 × the pin Fortran's own FMA-vs-no-FMA spread on the same inputs** (`l0_*_ffma`).
  - Unit 1's `ui_seg`: 1.1e-11 on random targets placed near the segment line. The Fortran built with FMA moves by the same amount.
  - Unit 2, realistic lattices: 8.2e-14, the same as the Fortran with FMA.
  - Exception: unit 4a defaults to no FMA (decision #15).
- [ ] Build a scaled OLAF performance case: scale up `ad_BAR_OLAF` or `ad_B1n2_OLAF`, or ask NLR for the OC7 Phase III case.
- [ ] Baseline r-test on our cluster at `OMP_NUM_THREADS=1`. Record known failures and timings before any change.
- [ ] Check whether FAST.Farm cases can run OLAF per turbine. If yes, FAST.Farm joins the API gate for the OLAF unit.
- [ ] Assign people to lanes.
- [x] **Kokkos version:** 5.3.0, as decision #6 says (user, 2026-11-14). It is built, and the template passes on it.
- [x] Run the toy on an h100.
  - Kokkos 4.5.1: job 1852028 passed. L0 206/206 at worst relative error 2.3e-18, `capi_smoke` OK.
  - Kokkos 5.3.0: job 1858536 **passed** (H100 NVL): L0 206/206, worst relative error 2.3e-18, `capi_smoke` OK.
- [ ] **CPU particle kernels: bitwise or vectorized?** (`tests/biotsavart-direct/unit.md`)
  - On Kokkos-OpenMP, `ui_part_nograd` is class B (bitwise against the OPENMP=OFF Fortran), but it runs at only 0.48–0.91× the Fortran built with `-fopenmp`, which vectorizes with `!$OMP SIMD reduction`. That fails P2 (≥ 0.95×).
  - Options:
    - (a) keep B and accept the gap;
    - (b) vectorize like the Fortran (class R1, about 1e-15);
    - (c) a runtime switch.
  - `ui_seg` is at 1.02×.
- [ ] Per-call device allocation in the direct kernels (P6 FAIL). The fix is a persistent handle held by `InducedVelocitiesAll_*` / `LiftingLineInducedVelocities`, which is the next unit.
- [x] OLAF placement and fusion analysis: [olaf-placement.md](olaf-placement.md). Key finding: lifting-line induction is a direct full-wake sum, 2–4× per step, the same order of cost as the wake treecode. It should reach the GPU with the wake kernel (M3), not after it.
