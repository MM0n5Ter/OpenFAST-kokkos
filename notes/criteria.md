# Porting criteria: what "done" means for a unit, check by check

Drafted 2026-11-14. This refines the three gates in [decisions.md](decisions.md) #7 (API, Regression, Performance) into numbered, checkable items. It adds two gates that the OLAF analysis ([olaf-placement.md](olaf-placement.md)) showed we need: **Placement** and **Fusion legality**.
- Every item has an ID. The unit card (`tests/<unit>/unit.md`) lists each ID with PASS/FAIL/N-A and evidence (a command and its output, or a file:line).
- `tools/gate.sh` (harness.md §7) prints the same IDs.

## A. API gate: callers and downstream users see no change

| ID | Check | How (mechanical) |
|---|---|---|
| A1 | Fortran module names and public interfaces are unchanged | `nm` comparison of `_MOD_` symbols, base vs port: nothing disappears (harness.md §5.2) |
| A2 | The diff touches only the unit's files plus Lane-A build files | `git diff --name-only <pin>..HEAD` ⊆ unit card file list |
| A3 | C symbols of every shared library are unchanged; new ones carry the `ofk_` prefix or come from Kokkos | harness.md §5.3 |
| A4 | Inputs, output channels (`.out`/`.outb` headers byte for byte), VTK file and array names are unchanged | harness.md §4 |
| A5 | `-DOPENFAST_KOKKOS=OFF` builds the upstream code path. `DOUBLE_PRECISION=OFF` still compiles. | Two extra configure+build jobs. The OFF build must match the baseline's `_MOD_` list exactly. |
| A6 | **Several instances** in one process (multi-rotor, FAST.Farm, Python) keep separate state: no C++ globals hold per-instance data | `ad_QuadRotor_OLAF`, plus a two-instance driver; grep the unit's C++ for `static` non-const objects |
| A7 | **Residency contract** (only once a unit owns device-resident state, olaf-placement.md §4): every Fortran entry point that reads the state sees current values, and every **copy, rollback and re-evaluation** path keeps working. In OLAF that means `FVW_CopyContState` (glue predictor-corrector, `m%x1/x2` overcycling `FVW.f90:649,680`) and `UpdateStates` called twice at the same `t` (`RollBackPreviousTimeStep`, `FVW.f90:597-610`). | The unit card lists those entry points (from `scripts/callers.sh` + grep). A test calls each one after a device-side update. Another test re-runs `UpdateStates` at the same `t` and checks the result matches a single call. |
| A8 | Kokkos version: builds against Kokkos 5.3.0 on both installs (`-omp`, `-cuda-h100`) | CMake uses `find_package(Kokkos 5.3...<6)`; code is C++20 (decision #6) |

## R. Regression gate: same answers, at the declared equivalence class

**Modes and the equivalence chain** ([architecture.md](architecture.md) §3–4). Every kernel has two modes from one source:
- **strict**: Fortran order, no FMA. It must be class **B** against the Fortran on Serial/OpenMP.
- **fast** (the default): free order and layout, FMA on.

Correctness is argued as a chain:
1. Fortran = strict CPU, bitwise.
2. strict CPU ≈ strict GPU, a few ulp.
3. strict ≈ fast, at the declared class (R1 or R2), on the same backend.
4. The whole simulation stays within R8 and r-test tolerances.

CPU↔GPU bitwise agreement is **not** required (decision #12). Both modes must be deterministic run to run (R7). Error bounds are relative to Σ|contributions| per target, so cancellation does not cause false failures.

**Equivalence classes.** Each kernel and each fusion declares one class for its **fast** mode in the unit card. Moving to a looser class needs a human sign-off recorded in the unit card.

| Class | Meaning | L0 tolerance (max\|port−ref\| / max\|ref\|) | Where it applies |
|---|---|---|---|
| **B** (bitwise) | Same arithmetic, same order | 0 on Kokkos-Serial/OpenMP with `-ffp-contract=off` on both sides; ≤ 1e-12 on CUDA (FMA, `sqrt`/`exp` ulps) | Elementwise kernels, packing, "square" N-body sums with per-target serial loops |
| **R1** (reordered) | Same math, different summation order | ≤ 1e-12 on all backends, **and bitwise run-to-run** for a fixed backend and thread count (no float atomics in sums) | Skinny LL sums with team reductions, mean-velocity reductions, split sums (F5-ii, F8) |
| **R2** (approximation) | A different valid approximation (treecode, cached far field) | Error vs the *direct* sum ≤ the Fortran tree's error vs direct on the same input × 1.1 | Treecode, multipole order, MAC variants |
| **S** (system) | r-test tolerances | rtol 1e-2, magnitude-scaled atol (testing.md) | Whole-module runs |

| ID | Check | How |
|---|---|---|
| R1 | L0 passes at the declared class on Serial/OpenMP **and** CUDA | `l0_*` executables; seeds 1 and 7, 200 cases; seed logged |
| R2 | L0 **covers every branch**: every `select case` value (`RegFunction` 0–4, `RegFunctionPart` 0–2, `VelocityMethod` 1–4), every guard in the audit list, every flag combination of packing (`bShedVorticity` × `bShedLastVorticity` × mirror) | The unit card's audit list maps each item to an L0 case name. Coverage: `gcov` on the frozen Fortran during L0 shows every branch line hit. |
| R3 | **Flag- and precision-dependent constants** are reproduced (e.g. `MINNORM=1e-4` is a default-kind literal: exactly 1e-4 under OpenFAST's double-build `-fdefault-real-8`, single otherwise) | `scripts/literal-kinds.py <unit files>` lists the default-kind real literals. Each one is passed from the Fortran shim (or taken as `ReKi`) rather than hard-coded, or it is declared as a class-A change. The oracle is built with the pin's exact flags. One L0 case per threshold constant, with inputs straddling it. |
| R4 | Sentinel tails are untouched; `Uind` is accumulated into | L0 sentinel and nonzero-initial-value cases (template `l0_ui_part`) |
| R5 | Captured inputs (L0b) from r-test OLAF runs replay at the declared class | harness.md §2 |
| R6 | `aerodyn_utest` passes; `ctest -R 'OLAF\|5MW_MRSemi_DLL_WSt_WavesIrr'` has no case that passes on the baseline and fails on the port | `OMP_NUM_THREADS=1` (testing.md) Note: identical `.outb` means agreement to about 1.5e-5 of each channel's range (16-bit scaling), not bitwise. |
| R7 | Determinism: two runs of the same case on the same backend give identical outputs; the fast mode also gives them across thread counts | `diff` of two `.outb` files. For L0, compare dumped outputs by md5 across processes and thread counts, and **reject empty dumps** (template: `tests/biotsavart-direct/r7_threads.sh`; gfortran `newunit=` once emptied them silently). |
| R8 | **Trajectory drift stays inside the baseline's own envelope.** Free-vortex wakes amplify 1e-16 differences over many steps (roll-up is chaotic), so a class-R1 kernel can pass L0 and still drift on a long run. Reference envelope: the Fortran baseline against itself with `OMP_NUM_THREADS=1` vs N, whose reductions are already reordered (#3401/#3402). | On the scaled case: per output channel, max\|port − base₁\| ≤ 2 × max\|base_N − base₁\| over the same window. Report the step where each first exceeds 1e-6 relative. Use **full-precision text output**: `.outb` stores 16-bit scaled channels and is too coarse. If base_N equals base_1 (envelope 0, as in units 1–2), the proposed envelope is the pin built with FMA (`-march=x86-64-v3`) at 1 thread; this awaits user sign-off. |

## P. Performance gate: faster where it matters, never slower where it doesn't

| ID | Check | Threshold |
|---|---|---|
| P1 | Kernel time vs Fortran+OpenMP on the same node and the same N | Reported per kernel; the unit card states the target |
| P2 | **CPU users lose nothing**: `OPENFAST_KOKKOS=ON` with the OpenMP backend vs the Fortran+OpenMP build | Time per step ≥ 0.95× Fortran speed at every N in the sweep Compare like with like: the fork built with `OPENMP=ON` plus Kokkos-OpenMP, so that Fortran code not yet ported keeps its threads (unit 2: an OPENMP=OFF fork looked 2× slower on the tree wake only because of this). |
| P3 | GPU time is broken down into kernel / H↔D transfer / host work left over (Amdahl) | All three are reported. Transfers > 25% of the GPU step time trigger a residency item (olaf-placement.md §4). |
| P4 | **Small-case overhead**: the r-test-sized cases do not get slower | ≤ 1.10× of the baseline's wall time. If the GPU loses below some N, the unit falls back to the host execution space below a measured `N_min` (runtime switch, same source). |
| P5 | Scaling sweep N_P = 10³…10⁶ with the crossover point (GPU vs CPU) | One CSV row per run (harness.md §6) |
| P6 | Memory: device memory is allocated in `*_Init` or grown by doubling, never per call in the time loop | Kokkos Tools memory-events log shows no allocations after step 3 |

## L. Placement gate (new): is this routine on the right side?

A routine goes to the **device** only if all of these hold:

| ID | Condition |
|---|---|
| L1 | Its work scales with the wake (a loop over N_P, N_S, N_part or M targets), not with nSpan/nWings: `scripts/kernel-inventory.py` column `wake` > 0 |
| L2 | It runs per step (or per DTfvw), not at init or output time only (from the call graph and the step anatomy, olaf-placement.md §1) |
| L3 | Its inputs are already on the device, or are produced by a neighbour that will also move, or the transfer is < 25% of its own measured CPU time |
| L4 | It calls nothing that must stay on the host (other modules' Fortran: AirfoilInfo, UA, InflowWind, ModMesh; I/O; LAPACK) |

Anything that fails L1–L4 stays on the host. The decision and its reason go in the module's placement note.

## F. Fusion legality gate (new): may these two steps be one kernel?

| ID | Condition |
|---|---|
| F1 | Producer → consumer over the same (or a closed-form-mapped) index space, with an O(N) intermediate that **no caller reads** (otherwise it is API-visible: keep it, or write it as a side output) |
| F2 | The fused version's equivalence class is declared, and it is no looser than the class the unit card accepts for that output |
| F3 | Serial counters (`iHead = iHead + 1`) are replaced by closed-form offsets, and an L0 case checks the offsets for every flag combination |
| F4 | No fusion across a framework API boundary (`*_UpdateStates` ↔ `*_CalcOutput`): the glue may call them any number of times in any order (predictor-corrector, Jacobians) |
| F5 | The unfused kernels stay available behind a switch (`OFK_FUSE=0`) until the fused one has passed R1–R7 twice |

## E. Cost and effort: money and time per unit (every unit, for the scale-up projection)

Measured with `tools/cost/` (method and rules in [tools/cost/README.md](../tools/cost/README.md)). One row per session phase goes into `tools/cost/ledger.csv`, and the unit card sums them.

| ID | Record | Source |
|---|---|---|
| E1 | Size: hand-written Fortran lines in scope; C++ lines written; shim lines | `git diff --stat`, `scripts/loc-by-module.sh` |
| E2 | **Money, model:** tokens by kind (input, output, cache write, cache read) and **USD (API-equivalent)**, per phase (analysis, port, debug, review), including subagents | `tools/cost/session_cost.py` |
| E3 | Human interventions: count, plus one line each in the unit card saying what the agent could not do | `human_prompts` column + unit card |
| E4 | Gate iterations: how many `tools/gate.sh` runs until all PASS, and which ID failed most | gate JSON summaries |
| E5 | **Time:** agent wall-clock and active hours per phase; human attention in minutes (self-reported); calendar days from branch creation to merge | `session_cost.py` (`wall_h`, `active_h`), ledger `human_minutes`, git |
| E6 | **Money, compute:** CPU core-hours and GPU-hours of the unit's builds, tests and performance runs (jobs named `ofk-<unit>-*`) | `tools/cost/compute_hours.sh` |
| E7 | E_k (ParEval-Repo) and the derived rates: $ per kLOC ported, active agent hours per kLOC, human minutes per unit | computed from E1–E6 at the end of each unit |

A unit is not "done" until E1–E6 are in its card. These rates drive the projection in roadmap.md (phases 2–4).
