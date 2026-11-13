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
10. **We report agent cost and effort per unit:** tokens, dollars, E_k (from ParEval-Repo), and the number of human interventions.
11. **Precision.** Target `DOUBLE_PRECISION=ON` (ReKi = 8 bytes, the default). The single-precision build must keep compiling, since upstream just added single-precision unit tests (#3511, the pinned commit). The C++ side picks its real type from the `OPENFAST_DOUBLE_PRECISION` macro that OpenFAST already defines for all languages.

## Open items / TODO
- [ ] **Where does the port live?** Recommendation: a fork (`MM0n5Ter/openfast`, branch `kokkos`) added as a second submodule, so the work can be offered upstream as a PR. `deps/openfast` stays the untouched baseline either way.
- [ ] Build both Kokkos installs and run `examples/fortran-kokkos` on them (OpenMP on a login node, CUDA on an `a100` node). Until then, [build-kokkos-fortran.md](build-kokkos-fortran.md) is unverified.
- [ ] Choose the compiler setup for the CUDA build: `nvcc_wrapper` for the whole configure, or `kokkos_compilation(TARGET ofk_kernels)`. See [build-kokkos-fortran.md](build-kokkos-fortran.md) §5.
- [ ] Build a scaled OLAF performance case: scale up `ad_BAR_OLAF` or `ad_B1n2_OLAF`, or ask NLR for the OC7 Phase III case.
- [ ] Baseline r-test on our cluster at `OMP_NUM_THREADS=1`. Record known failures and timings before any change.
- [ ] Check whether FAST.Farm cases can run OLAF per turbine. If yes, FAST.Farm joins the API gate for the OLAF unit.
- [ ] Assign people to lanes.
