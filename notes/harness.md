# How to build the test harness for a translated file or library

Every unit, whether one file (`FVW_BiotSavart.f90`) or a library (`InflowWind`), gets the same harness. It has one check per gate, at four levels (L0–L3, see [testing.md](testing.md)). The template is `examples/fortran-kokkos`.

> **Status:** the L0 pattern (the renamed-module oracle) was compiled and run on the toy, against a C++ stand-in for the kernel. Everything else in this note is the plan for weeks 1–2 (Lane C).

## 0. Per-unit layout

```
tests/<unit>/
  unit.md            unit card: files, contract, callers, downstream, owner, status, agent cost
  CMakeLists.txt     oracle library + L0 executables; added to the fork's build when BUILD_TESTING=ON
  l0_<routine>.f90   differential tests (one per boundary routine)
  capture/           inputs/outputs captured from r-test runs (git-ignored if large; seeds + scripts kept)
  perf/              the scaled case and its run script
```

Unit card template (`unit.md`):
```markdown
# <unit>
- Files (frozen → port): modules/aerodyn/src/FVW_BiotSavart.f90, …
- Contract: routines + signatures, side effects (accumulates Uind; touches only 1:nCPs), index conventions (SegConnct 1-based)
- Callers (scripts/callers.sh): …     Downstream (downstream.md): …
- L0: l0_ui_seg, l0_ui_part, l0_tree      L1: aerodyn_utest, ctest -R OLAF      L2: all     L3: py_ad_B1n2_OLAF, …
- Audit list: (from the agent) index offsets, layouts, reductions, constants, each with the Fortran file:line
- Owner / status / tokens / $ / E_k / interventions
```

## 1. L0a: the twin-module oracle (frozen Fortran vs port, in one process)

Purpose: compare the port against the **exact frozen Fortran**, routine by routine, on many inputs, in seconds.

How:
1. Take the unit's files **from `deps/openfast`** (never from the fork). Rename only the unit's own modules with `examples/fortran-kokkos/rename_modules.py`. For example, `FVW_BiotSavart` → `FVW_BiotSavart_ref` and `FVW_VortexTools` → `FVW_VortexTools_ref`, and every `use` of them inside those files. Fortran mangles module procedures as `__<module>_MOD_<name>`, so `__fvw_biotsavart_ref_MOD_ui_seg` and the port's `__fvw_biotsavart_MOD_ui_seg` coexist.
2. Their dependencies (`NWTC_Library`, `FVW_Types`) are **not** renamed. They are linked once, so both sides see the same kinds and types.
3. Build the renamed files into `oracle_<unit>`, with its own `Fortran_MODULE_DIRECTORY`.
4. Check that the oracle exports nothing else that could clash. This must print nothing:
   ```bash
   nm --defined-only liboracle_<unit>.a | awk '$2=="T"{print $3}' | grep -v '_ref_MOD_'
   ```
   A hit means the file has `bind(C)` names or non-module procedures. Rename those too, or move that routine to L0b.
5. The L0 driver (Fortran, so it calls both through their real Fortran interfaces) `use`s `<mod>_ref` and `<mod>` side by side. It makes inputs, calls both, and compares. Use the test-drive framework that OpenFAST already ships (`unit_tests/test-drive`) once there is more than one test.

**Inputs.** Each test has fixed edge cases and then randomized cases. Log the seed (`OFK_L0_SEED`), so any failure can be replayed.

| Edge case | Why |
|---|---|
| n = 0, n = 1 | Empty loops, off-by-one |
| Points on top of each other (r → 0) | Singularity guards (`MINNORM`) |
| Very small and very large core radii | Regularization branches; `exp` underflow |
| Every `RegFunction` / `VelocityMethod` value | Every `select case` branch |
| Arrays larger than n, with sentinel tails | "3 × n++": the port must not touch the tail |
| Nonzero initial output | Accumulate vs overwrite |
| Non-contiguous sections | Copy-in/out through the `bind(C)` interface |
| Realistic geometry (from L0b captures) | Random clouds miss wake structure (sheets, roll-up) |

Randomized inputs also stop a port from passing while computing less than the full kernel (the "reward hacking" seen by Gupta et al.).

**Comparison.** For each output: `max|port − ref| / max(max|ref|, tiny) ≤ rtol`, plus exact equality on the sentinel tails.

| Backend | rtol | Notes |
|---|---|---|
| Kokkos-Serial / OpenMP, same summation order | 1e-14 (bitwise with `-ffp-contract=off` on both sides) | Per-point sums are independent of thread count |
| Kokkos-CUDA | 1e-10 | FMA, device `exp`/`sqrt` |
| Treecode (any backend) | Not compared bitwise with the Fortran tree | Compare tree and direct sum on both sides. The port's tree error must stay within the Fortran tree's error on the same inputs (as `Test_BiotSavart_PartTree` does). |

## 2. L0b: capture and replay (inputs from real runs)

Purpose: test the port on the inputs it will really get: real wake geometry, sizes and growth from step to step. Also, a routine that is hard to call in isolation can be tested this way.

How:
1. Make a **copy** of the baseline source in the build tree (never edit `deps/openfast`). Add a probe at the boundary routine, e.g. `InducedVelocitiesAll_Calc`.
2. When `OFK_CAPTURE_DIR` is set, the probe writes the routine's inputs before the call and its outputs after it, for the first K calls (`OFK_CAPTURE_MAX`, default 20).
3. Format: one file per call, Fortran `access='stream'`. Each array is written as a record: name length and name, kind (`r8`/`i4`), rank, extents, data. Small and language-neutral, so C++ tests can read it too.
4. Run the OLAF r-test cases (`ad_BAR_OLAF`, `ad_B1n2_OLAF`, `ad_QuadRotor_OLAF`, …) and the scaled case with capture on.
5. The replay test reads each record, calls the **port** through the same Fortran interface, and compares with the captured outputs at the L0 tolerances.
6. Keep the capture scripts and seeds in git. Keep large capture files outside git and regenerate them from the pin.

## 3. L1 and L2: OpenFAST's own tests, through ctest

```bash
# L1 for the OLAF unit: unit tests + the r-test cases that run OLAF
OMP_NUM_THREADS=1 ctest --test-dir build-of-kk-omp -R 'aerodyn_utest|OLAF|5MW_MRSemi_DLL_WSt_WavesIrr' --output-on-failure
# L2: everything
OMP_NUM_THREADS=1 ctest --test-dir build-of-kk-omp -j 16
```
- Select OLAF cases by **name** (`-R OLAF`). The `ad_*_OLAF` cases carry the labels `aerodyn;bem`, so `-L olaf` misses them.
- Compare against our **baseline run**, not against "all pass". The gate is: no case that passes on `build-of-base` fails on the port.
- r-test writes per-case outputs and plots under the build tree. On failure, keep the `.out`/`.outb` pair for the unit card.

## 4. L3: downstream checks

| Check | How |
|---|---|
| Python C binding | `ctest -R 'py_ad_B1n2_OLAF\|_OLAF_py'` (pyOpenFAST loads `libaerodyn_inflow_c_binding.so` and `FAST_Library` through ctypes). Pattern from `examples/fortran-kokkos/tests/capi_smoke.py`. |
| Static build | Configure once with `BUILD_SHARED_LIBS=OFF` (the Simulink requirement) and link `openfast` and `aerodyn_driver` |
| Several instances | `ad_QuadRotor_OLAF` (multi-rotor), plus a small test that runs two AeroDyn/OLAF instances in one process |
| Output headers | Channel names and units in `.out`/`.outb` headers equal the baseline's, byte for byte |
| VTK | With `WrVTK` on, the same file names and the same array names and sizes as the baseline |

## 5. API gate: mechanical checks

1. **Scope of the diff.** `git diff --name-only <pin>..HEAD` in the fork lists only the unit's files (from the unit card) plus build files owned by Lane A. No caller is edited.
2. **Fortran interface.** The unit's modules keep their public names and argument lists. Proof: every caller compiles unchanged. Also compare the exported module symbols:
   ```bash
   nm --defined-only build-of-base/modules/aerodyn/libaerodynlib.a  | awk '$2=="T"{print $3}' | grep -i '_MOD_' | sort > base.txt
   nm --defined-only build-of-kk-omp/modules/aerodyn/libaerodynlib.a | awk '$2=="T"{print $3}' | grep -i '_MOD_' | sort > port.txt
   comm -23 base.txt port.txt          # must be empty: nothing the baseline exported has disappeared
   ```
3. **C symbols.** The same comparison with `nm -D --defined-only` on every shared library the baseline builds (`find build-of-base -name "*.so"`; the build needs `BUILD_SHARED_LIBS=ON` for most of them, `libaerodyn_inflow_c_binding.so` is always shared). New symbols are allowed only with the `ofk_` prefix, or from Kokkos.
4. **Inputs.** The r-test input decks run unchanged. No new required input keys.

## 6. Performance harness

- **Case.** A scaled OLAF case: start from `ad_BAR_OLAF` or `ad_B1n2_OLAF`. Sweep the number of wake panels (`nNWPanels`, `nFWPanels`) so wake markers range from 10³ to 10⁶. Fix the number of steps, and exclude the first steps as warm-up.
- **Matrix:**
  - Fortran with `OMP_NUM_THREADS` = 1 and = N (all cores on the node);
  - Kokkos-OpenMP with N threads;
  - Kokkos-CUDA on A100, and on H100 for the final results.
  - Each run uses `VelocityMethod` 1 (direct segments), 2 (particle treecode, the default) and 3 (direct particles).
- **Timers** (all three):
  - wall-clock of `aerodyn_driver`;
  - OLAF's own tic/toc (set `OLAF_PROFILING = .True.` at `FVW_Subs.f90:61` in the *perf copy* of each build);
  - Kokkos Tools kernel timer for the port (`KOKKOS_TOOLS_LIBS=libkp_kernel_timer.so`, from github.com/kokkos/kokkos-tools, built once).
- **Fairness.** Same compiler flags, same node type, same thread pinning (`OMP_PROC_BIND=spread OMP_PLACES=threads`). The GPU time includes host↔device copies.
- **Output.** One CSV row per run: `commit, build, backend, VelocityMethod, n_markers, steps, wall_s, s_per_step, kernel_s`. A plot script turns it into the scaling figures.

## 7. One command per gate

Lane C writes `tools/gate.sh <unit> {api|reg|perf}` in week 2, so the agent and the humans run the same checks:
- `api`: §5 checks 1–4.
- `reg`: L0a + L0b + L1 (every commit). Add `--full` for L2 + L3.
- `perf`: the matrix in §6, appended to `tests/<unit>/perf/results.csv`.

Each prints PASS/FAIL per check and writes a JSON summary, which goes into the unit card and into the agent's next prompt.

## 8. From the toy to the OLAF unit

| Toy (`examples/fortran-kokkos`) | OLAF unit |
|---|---|
| `ref/toy_biotsavart.f90` | `deps/openfast/modules/aerodyn/src/FVW_BiotSavart.f90`, `FVW_VortexTools.f90`, and the induction drivers in `FVW_Subs.f90` |
| `ref/toy_kinds.f90` (not renamed) | `NWTC_Library`, `FVW_Types` (not renamed) |
| `src/toy_biotsavart.f90` (same interface, shim body) | The fork's versions of those files, with `#ifdef OPENFAST_KOKKOS` shims |
| `src/toy_biotsavart_kokkos.cpp` | `src/modules/kokkos/src/<unit>/` (one directory per unit: direct kernels, then the treecode) |
| `src/ofk_runtime.cpp` | `modules/kokkos/src/common/ofk_runtime.cpp` (shared by all units) |
| `src/toy_c_binding.f90` + `tests/capi_smoke.py` | `AeroDyn_Inflow_C_Binding.f90` + the `py_ad_*_OLAF` r-test cases |
| `tests/l0_ui_part.f90` | `tests/fvw-induction/l0_ui_seg.f90`, `l0_ui_part.f90`, `l0_tree.f90`, plus replay of captured `InducedVelocitiesAll_Calc` calls |
