# Testing

## What OpenFAST already has
Counts from `scripts/rtest-inventory.sh`:
- **r-test** (`reg_tests/`, data in the `reg_tests/r-test` submodule): **213 CTest cases**, each running a binary on an input deck and comparing outputs against a stored baseline.
  - By registration function:
    - `of_regression` 56, `hd_regression` 26, `of_regression_linear` 20;
    - `ad_regression` 19, `md_regression` 18, `sd_regression` 16;
    - `seast_regression` 11, `of_regression_py` 8, `bd_regression` 8;
    - `ff_regression` 7, `ifw_regression` 6, `of_cpp_interface_regression` 3;
    - Python bindings: `py_*_regression` 7.
  - **16 OLAF cases:**
    - `ad_` B1n2, BAR, BAR_RegFunctionPart2, EllipticalWingInf, HelicalWakeInf, Kite, QuadRotor, Sphere and VerticalAxis (all `_OLAF`);
    - glue: `5MW_Land_QuadRotor_DLL_WSt_OLAF`, `5MW_MRSemi_DLL_WSt_WavesIrr`, `EllipticalWing_OLAF`, `HelicalWake_OLAF`;
    - Python: `EllipticalWing_OLAF_py`, `HelicalWake_OLAF_py`, `py_ad_B1n2_OLAF`.
  - Select them with `ctest -R OLAF`, plus `-R 5MW_MRSemi_DLL_WSt_WavesIrr`. **Not** with `-L olaf`, because the `ad_*_OLAF` cases are labelled `aerodyn;bem`.
  - Tolerances: `CTEST_RTEST_RTOL "2"` and `CTEST_RTEST_ATOL "1.9"` (`reg_tests/CMakeLists.txt:85-86`) are *orders of magnitude*: relative 1e-2, and an absolute floor scaled to each channel's magnitude. These are engineering tolerances, loose enough for reordered GPU reductions.
  - The OLAF cases are labelled "VALUE FOR UNIT TEST ONLY": tiny wakes, good for correctness, useless for timing.
- **Unit tests** (`unit_tests/`, framework **test-drive**, `BUILD_TESTING=ON` + `BUILD_UNIT_TESTING=ON`): `aerodyn_utest` runs `modules/aerodyn/tests/test_AD_FVW.F90`, which calls the routines in `FVW_Tests.f90`:
  - `Test_BiotSavart_Sgmt` :189, `Test_BiotSavart_Part` :288, `Test_BiotSavart_PartTree` :342 (tree vs direct);
  - `Test_BiotSavart_SrcPnl` :454, `Test_SegmentsToPart` :494, `FVW_Test_WakeInducedVelocities` :759, `Test_SrcPnl_Sphere` :869.
  - Other modules have their own `*_utest`.
- **Profiling hooks:** `OLAF_PROFILING` (`FVW_Subs.f90:61`, a compile-time `.False.`) switches on tic/toc around the induction routines.

## Our test levels
| Level | What | Runs |
|---|---|---|
| L0 | Unit differential test, in one process. The frozen Fortran (module renamed `_ref`) is the oracle and the port is under test. Inputs are randomized and include edge cases, plus inputs **captured from real r-test runs**. | Every commit |
| L1 | The unit's upstream unit tests (`aerodyn_utest`) and the r-test cases that exercise it (`ctest -R OLAF`) | Every commit |
| L2 | All 213 r-test cases at upstream tolerances | Nightly and pre-merge |
| L3 | Downstream workflows from [downstream.md](downstream.md): the Python C-binding cases, the static (Simulink-style) build links, output header and VTK structure checks | Pre-merge for units on an external path; otherwise nightly |

How to build each level: [harness.md](harness.md).

- **Baseline first.** Before porting anything, run L1 and L2 on the **unmodified** pin and record known failures and timings. r-test baselines come from upstream's reference machines, so some cases may fail on our cluster before any change. A gate fails only on regressions.
- **Determinism.** OLAF with OpenMP is not deterministic (#3401/#3402), so baselines use `OMP_NUM_THREADS=1`.
- **Validation ladder** (from FESOM2):
  - bit-for-bit on the serial CPU backend where achievable;
  - tolerance otherwise;
  - the treecode is checked against the direct sum, not bitwise against the Fortran tree.

## Do we port the tests first? No
- **r-test does not depend on the source language.** It runs binaries and compares files, and nothing in it needs porting. That independence is what makes it a valid judge of the port.
- **The oracle must not be translated by the same agent before the code it checks**, or both can share the same mistake. Tests stay at least one step behind the code they check.
- **The Fortran unit tests are API consumers.** `FVW_Tests.f90` calls `ui_seg`, `ui_part_nograd` and the tree routines through their Fortran interfaces, and it keeps doing so after the port. That exercises the API gate for free. They are ported or replaced by C++ tests only **after** their unit, before the shim is deleted.
- **Before porting a unit, we *add* tests:**
  1. L0 tests with randomized and captured inputs;
  2. a scaled performance case;
  3. output-header and VTK-structure checks.
- **`deps/openfast` is the frozen oracle.** After a unit's Fortran is deleted from the port, L0 still compares against the baseline build.
