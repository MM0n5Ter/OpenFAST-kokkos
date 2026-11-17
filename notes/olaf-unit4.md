# OLAF unit 4: where the host time goes after units 1–3, and what unit 4 is

Drafted 2026-11-17 (session `199ad883`, unit-4 analysis subagent). The analysis is profile-first; nothing was ported. This note builds on [olaf-placement.md](olaf-placement.md) and refers to its sections instead of repeating them. Code is cited as `file:line` on the pin `12fa3fc`. The profile numbers come from `scripts/wsamp/` and `scripts/unit4-mb/` (§6).

## 0. Summary

**The profile overturns the plan.** After units 1–3, two host costs dominate every BAR run, on the CPU and on the GPU:
1. **The lifting-line circulation iteration, at 16–58% of the CPU wall time (49–58% on the tree cases) and 44–58% of the GPU runs.** It is `Wings_ComputeCirculationPolarData` → `ui_quad_n1` → `ui_seg_11`, and it is serial because of the `!OMP` typo.
2. **AeroDyn's wind and tower evaluation at the OLAF wind points, at 10–30% on the CPU and 32–38% on the GPU.** `AD_CalcWind` runs 5× per step and `TwrInflArray` 3× per step. 85% of the uniform-wind time is gfortran's IEEE entry/exit overhead.

The wake time-stepping that phase-2 residency would move to the device (CCSD means, Euler update, `PropagateWake`) takes **0.05–0.12 s of 40–110 s**. The host pack of `Sgmt`/`Part` takes 0.3–0.75 s.

**Profile** (ov-fomp = the P2 config, Kokkos-OpenMP + Fortran OpenMP, 8 threads, Xeon 8480+, 300 steps). Wall-clock samples on the main thread; seconds, with the share of the run in brackets (jobs 1923289, 1923304):

| Bucket | `bar-r8-tree` | `bar-r8-direct` | `bar-perf-tree` | `bar-perf-part` | `bar-perf-direct` |
|---|---|---|---|---|---|
| **Total** | **40.7** | **74.8** | **48.4** | **110.7** | **146.2** |
| Kernels: wake (tree build + eval, or direct) | 3.65 | 37.71 | 6.97 | 68.88 | 104.87 |
| Kernels: lifting-line induction (2 per step) | 0.69 | 0.67 | 1.06 | 1.06 | 1.02 |
| **Circulation solve: `ui_quad_n1` iterations** | **23.65 (58%)** | **23.81 (32%)** | **23.65 (49%)** | **23.86 (22%)** | **23.65 (16%)** |
| Circulation solve: rest (polar lookup, loop) | 0.19 | 0.16 | 0.19 | 0.18 | 0.19 |
| **`AD_CalcWind` (InflowWind at all OLAF wind points)** | **9.74 (24%)** | **9.70 (13%)** | **12.20 (25%)** | **12.20 (11%)** | **12.18 (8%)** |
| `SetInputsForFVW` (`TwrInflArray` at all wind points) | 1.45 | 1.44 | 2.49 | 2.60 | 2.48 |
| Pack/unpack (Fortran pack + `ofk` source set), wake + LL | 0.39 | 0.30 | 0.75 | 0.75 | 0.63 |
| CCSD means/dxdt + Euler + `PropagateWake` | 0.05 | 0.05 | 0.12 | 0.08 | 0.09 |
| FVW other (panelling, maps, interp, copies) | 0.16 | 0.16 | 0.32 | 0.35 | 0.32 |
| VTK/grid output | 0.14 | 0.17 | 0.13 | 0.14 | 0.19 |
| AeroDyn other, init, driver, I/O | 0.59 | 0.64 | 0.56 | 0.55 | 0.54 |

**GPU split (Amdahl).** It uses measured GPU runs plus the host profile of `ov-omp`. `ov-omp` has Fortran `OPENMP=OFF`, as every CUDA build does. Its host time × 0.82–0.84 reproduces the measured GPU host time (§1.3).

| | `bar-r8-tree` (A100 MIG, job 1887850) | `bar-r8-direct` (H100, job 1860779) | `bar-perf-direct` (H100, job 1860779) |
|---|---|---|---|
| Measured total | 34.1 s | 37.8 s | 46.2 s |
| GPU kernels (+ tree build, + transfers) | 1.6 s | 4.4 s | 6.9 s |
| Host: circulation solve (est.) | 19.7 s | 20.2 s | 20.4 s |
| Host: wind + tower at wind points (est.) | 11.9 s | 12.1 s | 17.4 s |
| Host: everything else (est.) | 0.9 s | 1.1 s | 1.5 s |

**Recommendation: unit 4 is the host-side lifting-line circulation solve (4a) plus the wind-point host work (4b). It is not residency.**

| | Scope | Class | Expected gain (BAR, from the profile and microbenchmarks) |
|---|---|---|---|
| **4a `circulation-solve`** | F10-B (§2.1): the geometry of the circulation iteration is built once per solve (a "square" Kokkos kernel in the **host** execution space), and each iteration replays only the Γ statement of `ui_seg_11` in the pin's order | **B** (bitwise, so the iteration count cannot change) | Circulation solve 23.7 s → **2.3 s serial / about 0.4 s with 8 threads**. CPU: **1.2–2.4×** whole run (tree cases most) |
| **4b `wake-wind`** | ~~H1 IEEE fix~~ (out of scope: InflowWind is an external library, user 2026-11-17); H2 skip 2 of the 5 wake-point wind evaluations per step, whose values are recomputed identically; H3 tower influence only on the slices read; H4 grid points only inside the grid window | **B** | Wind 9.7–12.2 s → about 0.3–0.5 s; tower (serial, as on the GPU) 4.8–8.4 s → 2.5–4.9 s |
| **4a + 4b** | | B | **CPU (P2 config): r8-tree 40.7 → about 7.3 s (5.6×); perf-tree 48.4 → about 12.4 s (3.9×); r8-direct 74.8 → about 41 s (kernel-bound, 1.8×).** **GPU: r8-direct 37.8 → about 8.4 s (4.5×); r8-tree 34.1 → about 5.4 s (6.3×); perf-direct 46.2 → about 13.4 s (3.4×).** After unit 4 the GPU runs are 30–52% kernel |

- **F10 is a bigger win than residency:** 23.7 s against ≤ 0.12 s of time-stepping (about 200×), or against 0.4–0.9 s with the pack included (about 30×).
- **Phase-2 residency (F2–F5, F9) stays deferred** until InflowWind's wake-point evaluation is on the device (Lane D). Until then the wind crosses the boundary every step whoever owns `x`, and the residency gain is ≤ 2% of a GPU step after unit 4 (§2.3).

**Hardest gates:**
- 4a:
  - **claim 1 / R8 must stay bitwise, including iteration counts**. One flipped convergence test moves Γ by up to `ConvCrit`·|Γ|. So fast mode = strict mode for this kernel (no FMA), which is a declared exception to architecture.md §3;
  - **P4**: 87 targets per iteration, so the serial fallback below `N_min` must be measured;
  - **P6**: an O(nCP²) cache that must live in a handle.
- 4b: **A2/A4 across AeroDyn and InflowWind** (a 7.8k-line `AeroDyn.f90` overlay copy and another lane's module), and the **proof obligation that each skipped slice is never read** (H2–H4, card audit W2–W4).
- A7 is **not** hard for unit 4: neither part keeps state between calls. It stays the hardest gate of the later residency unit (§2.3).

## 1. Profile: method and per-step anatomy

### 1.1 Method
- **Sampler** (`scripts/wsamp/wsamp.c`, `LD_PRELOAD`): a POSIX timer on `CLOCK_MONOTONIC` delivers `SIGPROF` to the **main thread only**, at 500 Hz, and the handler records `backtrace()`.
  - So the samples measure where the main thread spends **wall** time, including time it spends inside OpenMP and Kokkos-OpenMP regions as thread 0. That is the Amdahl view.
  - It needs no rebuild: the existing binaries have symbols, and gfortran's unwind tables are on by default on x86-64.
  - Overhead: `bar-r8-tree` takes 40.80 s sampled against 40.73 s unsampled in the same job (0.2%). `bar-perf-tree` takes 50.1 s against 47.9 s in job 1887061, a different node.
- **Buckets** (`scripts/wsamp/wsamp_report.py --buckets`): the first rule in a fixed priority order that matches any frame of the stack. Kernels are split by caller (`WakeInducedVelocities` vs `LiftingLineInducedVelocities`).
- **Call counts:** gprof on a new `-pg` build (`opt/build/u4-fomp-pg`, fomp config, 1 thread, job 1923304), which gives exact counts.
- **Microbenchmarks** (`scripts/unit4-mb/`, job 1923596):
  - one circulation iteration using the pin's own `ui_seg_11`/`ui_quad_n1`, extracted verbatim and compiled with the pin's flags;
  - `UniformField_GetVel` with and without the IEEE `use`.
  - They reproduce the profile: 3.08 ms per iteration × 7720 iterations = 23.8 s, against 23.65–23.87 s sampled.

### 1.2 Calls per step (gprof, `bar-r8-tree`, 300 steps; `bar-r8-direct` and `bar-perf-tree` match within 0.1%)

| What | Calls / step | Note |
|---|---|---|
| `AD_UpdateStates` → `FVW_UpdateStates` | 1 | AeroDyn driver; the glue has 1 + `NumCrctn` |
| Circulation solves | 1 | `FVW.f90:727` → `FVW_Wings.f90:295` |
| Circulation iterations | **25.7** (7720 / 300; cap 30, `FVW_IO.f90:47`) | `CirculationFromPolarData` 7721 = 7720 + the first-call guess |
| `ui_quad_n1` | 13 433 (= 25.7 × 522) | 522 = nCP_tot(87) × iNWStart(2) × nWings(3) |
| `ui_seg_11` | 1.56 M (= 25.7 × 60 552) | 60 552 = 87 targets × 87 rings × 2 depths × 4 segments |
| `LiftingLineInducedVelocities` | 2 | LL #b, #d (olaf-placement §1); `ofk_iv` evals 903 = 300 wake + 600 LL + 3 grid |
| `WakeInducedVelocities` (tree build + eval) | 1 | 303 builds = 300 + 3 grid outputs |
| `AD_CalcWind` (→ `AD_CalcWind_Rotor` 1501) | **5** | 4 in `AD_UpdateStates` (`AeroDyn.f90:1977` × size(u) = 2, `:1996` × 2), 1 in `AD_CalcOutput` (`:2365`) |
| `UniformField_GetVel` (points) | 93 k (r8), 117 k (perf) | 5 × all `r_wind` points with z > 0: LL + NW + FW at **maximum** size (`FVW_Subs.f90:578-580`) + 9600 grid points (`:606-620`) |
| `SetInputsForFVW` → `TwrInflArray` | **3** | 2 in `AD_UpdateStates` (`:2027`), 1 in `AD_CalcOutput` (`:2378`); over **all** `r_wind` points (`:4229`) |
| `PropagateWake`, `FakeGroundEffect`, `SetRequestedWindPoints` | 1 each | |
| `Map_LL_NW`, `Map_NW_FW`, `Wings_Panelling` | 3 each | |

### 1.3 The GPU host: why "× 0.82–0.84 of `ov-omp`"
- CUDA builds compile the Fortran with `OPENMP=OFF` (`ov-cuda-*/CMakeCache.txt`). On a GPU node, therefore, `TwrInflArray`'s `!$OMP` loop (`AeroDyn.f90:6123-6139`) runs serially, while the circulation solve and the wind are serial in every build.
- `ov-omp` (Kokkos-OpenMP, Fortran `OPENMP=OFF`) is that host side on the Xeon. Its non-kernel time is 39.8 s on both r8 cases.
- The measured GPU host time is 0.84× that on the H100 node: 37.8 − 4.34 kernel − 0.1 transfers = 33.4 s for `bar-r8-direct`.
- It is 0.82× on the A100 MIG node: 34.1 − 1.05 − 0.42 build − 0.1 = 32.5 s for `bar-r8-tree`.
- The GPU rows in §0 and §3.2 scale `ov-omp`'s buckets by 0.84 (H100) or 0.82 (A100 MIG).

`ov-omp`, 8 threads, s (job 1923289):

| | circ | wind | tower | pack | step | rest | non-kernel total |
|---|---|---|---|---|---|---|---|
| `bar-r8-tree` | 24.03 | 9.71 | 4.77 | 0.32 | 0.04 | 0.90 | 39.8 |
| `bar-r8-direct` | 24.00 | 9.72 | 4.73 | 0.33 | 0.05 | 0.96 | 39.8 |
| `bar-perf-tree` | 23.88 | 12.18 | 8.40 | 0.76 | 0.12 | 1.00 | 46.3 |
| `bar-perf-part` | 23.99 | 12.19 | 8.39 | 0.73 | 0.11 | 1.00 | 46.4 |
| `bar-perf-direct` | 24.03 | 12.18 | 8.36 | 0.60 | 0.09 | 1.05 | 46.3 |

**Earlier readings this corrects.** The unit-2 card (job 1860779) put the remaining ~33 s of the H100 run in "wake stepping, lifting line, rest of AeroDyn" and named residency as the next gain. In fact:
- wake stepping is 0.05 s;
- the lifting-line *induction* is on the GPU;
- the 33 s are the circulation iteration (about 20 s) and the wind/tower at the wind points (about 12 s).

## 2. Candidates the profile points to

### 2.1 Circulation solve (4a)

**Call graph** (olaf-placement §1, refined):
```
FVW_UpdateStates (FVW.f90:567)
 └─ FVW_CalcConstrStateResidual (FVW.f90:727 → :1353) → Wings_ComputeCirculation (:1390 → FVW_Wings.f90:241)
     └─ Wings_ComputeCirculationPolarData (FVW_Wings.f90:295 → :314)
         ├─ [m%FirstCall] CirculationFromPolarData (:365)                      host, AirfoilInfo
         ├─ LiftingLineInducedVelocities(iDepthStart=iNWStart+1) (:391)        LL #b → unit 2 kernel (GPU)
         └─ do while (.not.bConverged .and. iIter<CircSolvMaxIter)  (:421)     25.7 iterations per solve (measured)
              ├─ Vvar = Σ_rings ui_quad_n1(...)  (:424-447)                    60 552 ui_seg_11 per iteration: 99% of the solve
              ├─ Vtot_CP = Vcst + Vvar (:449-455)
              ├─ CirculationFromPolarData (:458)                               AFI_ComputeAirfoilCoefs per CP (L4: host)
              └─ relaxation, MeanGamma, convergence test (:462-467)
```

**Sizes.**
- Targets nCP_tot = Σ nSpan, which is B·S: 87 for BAR, 150 for the B = 3, S = 50 example of olaf-placement §1.
- Rings nCP_tot × iNWStart (= 2, `FVW.f90:103-104`), 4 segments each.
- Pairs per iteration: 4·iNWStart·nCP_tot², i.e. 60 552 (BAR), 180 000 (S = 50), 720 000 (nCP = 300).
- **Independent of the wake length**: the cost is identical in the nNW = 120 and 300 cases.
- It is ×iterations (≤ 30) and ×(1 + `NumCrctn`) in the glue.

**Placement.**
- **Host execution space**, in every build:
  - L1 fails: the work scales with nSpan², not with the wake;
  - L4 fails for the iteration as a whole: `CirculationFromPolarData` calls AirfoilInfo;
  - the data are tiny: 87 targets.
- Kokkos rather than a Fortran-only restructure, because:
  - the CUDA builds have Fortran `OPENMP=OFF` but Kokkos-OpenMP in their host space (`kokkos-5.3.0-cuda-*` defines `KOKKOS_ENABLE_OPENMP`), so only a Kokkos host kernel gets threads on the GPU nodes;
  - it is one more Fortran-free piece (roadmap).

**What crosses.** Nothing crosses to the device. Fortran → C++ per solve, by pointer, no copy (host space):
- `x%W(:)%r_NW(:, :, 1:iNWStart+1)`, `x%W(:)%Eps_NW(1, :, 1:iNWStart)`, `m%W(:)%CP`;
- per iteration, `GammaLastIter` (nCP_tot) in and `Vvar` (3·nCP_tot) out.

**Who reads the state outside the unit.**
- The unit keeps no state between calls: the cache is rebuilt from the arguments on every call.
- `FVW_CopyContState`/`AD_CopyContState` from the glue (`FAST_Funcs.f90:459`, `FAST_Solver.f90:1275`, `FAST_Subs.f90:5299`) and the rollback at the same `t` (`FVW.f90:597-615`) therefore cannot see stale data. **A7 holds by construction.**
- The API-visible outputs (`z%Gamma_LL`, `m%W%Vtot_CP/Vind_CP = -9999`, `alpha_LL`, `Vreln_LL`) stay written by the unchanged Fortran.

**Fusion: F10 split into two variants** (olaf-placement §5 listed only the R1 form). Measured with the pin's own routines (`scripts/unit4-mb`, job 1923596, one core):

| nCP_tot | Pin, per iteration | F10-B: cache build (per solve) | F10-B: replay per iteration | F10-R1: matrix build | F10-R1: mat-vec per iteration | F10-B vs pin | F10-R1 vs pin |
|---|---|---|---|---|---|---|---|
| 87 (Vatistas, the default) | 3.08 ms | 3.45 ms | 0.164 ms | 3.52 ms | 0.035 ms | **bitwise** (all RegFunction 0–4) | 5.4e-16 relative |
| 150 | 9.08 ms | 9.81 ms | 0.493 ms | 9.63 ms | 0.155 ms | bitwise | 8.2e-16 |
| 300 | 36.3 ms | 35.1 ms | 1.98 ms | 38.6 ms | 0.63 ms | bitwise | 8.3e-16 |

- **F10-B (geometry cache + Γ-tail replay), class B.**
  - Everything in `ui_seg_11` before `FVW_BiotSavart.f90:125` is Γ-independent: `norm_a`, `norm_b`, `denominator`, the guards, `crossprod`, `norm2_r0` and the regularization `Kv`. Cache `Kv`, `norm_a+norm_b`, `denominator`, `crossprod(3)` and a guard mask once per solve.
  - Per iteration, compute `((Γ·fourpi_inv)·Kv)·(na+nb)/den` and `·crossprod`, then accumulate per target in the pin's order. That is the same operations in the same order, so the result is bitwise.
  - BAR projection: 300 × 3.45 ms + 7720 × 0.164 ms = **2.3 s serial (10.3×)**.
  - Parallel over targets: the per-target order is unchanged, so the result is still B. The ideal on 8 threads is about 0.3 s; with per-iteration launch latency, **≈ 0.4 s** is the estimate to measure (P2/P4).
- **F10-R1 (influence matrix, mat-vec), class R1.**
  - Only 1.7× better than F10-B serially (1.33 s), because the build dominates either way.
  - Every rounding difference can flip `maxval(abs(DGamma))/MeanGamma < ConvCrit` (`FVW_Wings.f90:467`) and change the iteration count (card audit C4).
  - **Not recommended.**
- Do not fuse the iteration with `CirculationFromPolarData`. It needs AirfoilInfo (host, L4), and its output decides convergence on the host.

**Traps** (full list in `tests/circulation-solve/unit.md`, audit C1–C12):
- the `!OMP` typo makes the pin serial (`FVW_BiotSavart.f90:527-528, 548-549`);
- the left-associated Γ product (`:125`);
- the +0 adds of skipped segments, which are exact to skip because `Vvar` starts at +0;
- convergence-flip amplification;
- ring corner order;
- `Eps_NW` row 1 for all four segments: not unit 2's segment pack;
- FMA must stay off in both modes.

### 2.2 Wind and tower at the OLAF wind points (4b)

**Call graph** (missing from olaf-placement §1, which starts at `FVW_UpdateStates`):
```
AD_UpdateStates (AeroDyn.f90:1942)
 ├─ do i=1,size(u): AD_CalcWind(utimes(i), u(i), …, m%Inflow(i))   (:1976-1979)   ×2: all r_wind points → InflowWakeVel
 ├─ do i=2,1,-1: AD_CalcWind(utimes(i), uInterp, …, m%Inflow(1))   (:1993-2003)   ×2: wake part recomputed, result unchanged (H2)
 ├─ do i=1,size(u): SetInputsForFVW(p, u(i), i, m)                 (:2026-2029)   ×2: V_wind = InflowWakeVel; TwrInflArray over all points (:4226-4235)
 └─ FVW_UpdateStates
AD_CalcOutput (AeroDyn.f90:2329)
 ├─ AD_CalcWind (:2365)                                                            ×1
 └─ SetInputsForFVW (:2378)                                                        ×1 (TwrInflArray over all points)
```

**Sizes.**
- `r_wind` = LL CPs + NW at maximum size, B(S+1)(nNWMax+1) + FW, 2B(nFWMax+1) + all grid points + source panels (`FVW_Subs.f90:568-638`).
- BAR r8: 87 + 10 890 + 66 + 9600 = 20.6k points. BAR perf: 36.9k points.
- Points with z ≤ 0 are skipped by the uniform field (`IfW_FlowField.f90:182`). `TwrInflArray` visits every point.

**The IEEE trap (H1).**
- `UniformField_GetVel` is `pure` and contains `use, intrinsic :: ieee_arithmetic` (`IfW_FlowField.f90:386`), but uses nothing from it.
- gfortran then wraps every call in `_gfortran_ieee_procedure_entry/exit`, which saves and restores the FP environment.
- In the profile, 7.95 s of `UniformField_GetVel`'s 9.38 s is that wrapper; `pow` is 0.98 s.
- The microbenchmark gives 14.4× per point and bitwise-identical output (`ieeemb`).
- The same wrapper shows up in the pin's `ui_leaves_nograd` (0.72 s on `bar-r8-tree.pin8`), which the ported treecode already avoids.

**Dead work (H2–H4), class B.**
- Which slices are read: only `u(1)`'s wake slice (`FVW.f90:645`), the LL slice through `uInterp` (`:1378`), the grid slice of `CalcOutput`'s `u` (`:1618`, used only when `bDoGrid`, `:1621-1627`) and the panel slice (`FVW_Subs.f90:2105`).
- H2: the `:1996` loop writes `m%Inflow(1)` last with `utimes(1)`, so its wake slice is bitwise the `:1977` value.
- H3: the tower influence on unread slices.
- H4: grid points outside the grid window. On BAR the grid is written only at t ∈ [1.2, 1.4] s, yet its 9600 points are evaluated 5× per step for 30 s.

**Placement.** Host:
- InflowWind is Lane D (L4);
- the tower influence reads AeroDyn's tower mesh.

**State.** None, so there is no A7 surface.

**Gain.**
- H1 takes the wind from 9.7 to about 1.0 s (r8) and from 12.2 to about 1.3 s (perf).
- H2 + H4 then cut the evaluated points per step from 5 × (wake + grid) to 3 × wake: ×0.29 (r8) to ×0.35 (perf). That gives **≈ 0.3–0.5 s**.
- H3 visits 0.51× (r8) to 0.58× (perf) of today's tower points. The tower goes from 1.45 to about 0.75 s with 8 threads, and from 4.8 to about 2.5 s serial (GPU nodes).

### 2.3 Wake time-stepping and phase-2 residency (deferred again)
- CCSD means/dxdt/eps + Euler + `PropagateWake` take 0.05 s (r8) and 0.08–0.12 s (perf) per run: **≤ 0.25%** on the CPU.
- On the GPU after unit 4, it would be about 1% of `bar-r8-direct`'s projected 8.4 s.
- F3, F4 and F9 (olaf-placement §5) remain legal as written; they just buy nothing yet.
- The pack is 0.3–0.75 s per run. F2 and F5 would remove it, but only with device-owned `x` (architecture §5: callers can read `Sgmt`/`Part`).
- **When residency does matter:** at wake sizes of 10⁵–10⁶ markers and once InflowWind evaluates the wake points on the device. Until then, the wind (O(N), 3–5 evaluations per step, with pow) is a larger O(N) host term than the time-stepping (O(N), about 10 passes, no transcendental).
- Its hardest gate stays A7: the glue copies AD's whole state CURR→PRED before every `AD_UpdateStates` (`FAST_Funcs.f90:459`) and once per correction (`FAST_Solver.f90:1275`), and PRED→CURR after the step (`FAST_Subs.f90:5299`), with generated Fortran deep copies.
- **Recommended design for that later unit:**
  - write the lattice back to `x` at every API exit (one D2H per `UpdateStates`, about 56 B × N_P, ≈ 2–3 ms at 10⁶);
  - validate it by content at the next entry, unit 2's slab `memcmp`.
  - This makes every glue copy and rollback a cache miss, never a wrong answer. A dirty-flag lazy sync would need a hook in generated code the overlay cannot touch.

### 2.4 What else the profile shows (not in unit 4)
- After unit 4, the GPU runs are kernel-bound: 4.3 of about 8.4 s for `bar-r8-direct`. The next items there:
  - the LL kernel's 2 calls per step with the same far wake (olaf-placement F5-ii, class R);
  - the direct wake kernel's shape.
- On the CPU, `bar-r8-direct` and `bar-perf-part` stay 85–90% kernel after unit 4 (the O(N²) direct paths). That is expected: those are the methods the treecode replaces.
- `TwrInflArray`'s remaining serial cost on GPU nodes is about 2.5 s after H3. With Fortran OpenMP off it is the largest host term left on a GPU node. Porting it is AeroDyn work (a Kokkos host loop or the `!$OMP` region in a C++ host kernel) for a later unit.
- `AD_CalcOutput`'s wake-point wind and tower evaluation (1 of the 5 and 1 of the 3) is read only on grid-output steps. Skipping it needs an FVW-side query ("is a grid due?") that AeroDyn can call: listed, not in 4b.

## 3. Recommended unit-4 boundary

### 3.1 Split
- **4a `circulation-solve`** (OLAF, `FVW_Wings.f90` overlaid; `ofk_kernels` gains `circulation_{strict,fast}.cpp`, both without FMA).
  - Shim only the `Vvar` block of `Wings_ComputeCirculationPolarData` (`FVW_Wings.f90:422-447`): `ofk_cs_build` once before the `do while`, `ofk_cs_vvar` per iteration. Everything else in the routine stays Fortran.
  - `ui_quad_n1` stays as it is.
  - Gates: claim 1 (bitwise, including per-solve iteration counts), R6, R7, R8 (bitwise, so drift 0), P2, P4 (serial fallback), P6, A1/A2/A5, A7 (stateless: the test still runs).
  - Card: `tests/circulation-solve/unit.md`.
- **4b `wake-wind`** (AeroDyn + one FVW pack routine; host Fortran only, no Kokkos): H2–H4.
  - **InflowWind is an external library call (user decision, 2026-11-17).** 4b changes only *how often and where* AeroDyn calls `IfW_FlowField_GetVelAcc`. It never edits, overlays or ports InflowWind.
  - So **H1 (the IEEE fix in `UniformField_GetVel`) is out of scope.** It is kept as an upstream suggestion for InflowWind's owners (audit #15). Without H1, each remaining wind point still costs the ~14× IEEE overhead. The wind projection in §3.2 assumed H1, so it is too optimistic: re-measure after H2–H4.
  - Separate gates: R6 (including `aerodyn_utest` and glue OLAF cases with `NumCrctn` ≥ 1); outputs and VTK bitwise from line 8; P2.
  - Its audit rows are W2–W6 in the same card. W1, the H1 row, is upstream only.
- Order: 4a first. It is the larger win, it is self-contained in OLAF, and it needs no other lane.

### 3.2 Expected gains (projection from the measured buckets)

| Case | Now | After 4a | After 4a + 4b |
|---|---|---|---|
| CPU `bar-r8-tree` (ov-fomp, 8 threads; pin+OMP 44.2 s) | 40.7 s | ≈ 17.4 s | **≈ 7.3 s** |
| CPU `bar-perf-tree` (pin+OMP 54.7 s) | 48.4 s | ≈ 25.2 s | **≈ 12.4 s** |
| CPU `bar-r8-direct` | 74.8 s | ≈ 51.4 s | ≈ 41.3 s (kernel 38.4 s) |
| CPU `bar-perf-part` | 110.7 s | ≈ 87.2 s | ≈ 74.4 s (kernel 69.9 s) |
| CPU `bar-perf-direct` | 146.2 s | ≈ 122.9 s | ≈ 110.1 s (kernel 105.9 s) |
| GPU `bar-r8-direct` (H100) | 37.8 s | ≈ 18.2 s | **≈ 8.4 s** |
| GPU `bar-r8-tree` (A100 MIG) | 34.1 s | ≈ 15.0 s | **≈ 5.4 s** |
| GPU `bar-perf-direct` (H100) | 46.2 s | ≈ 26.4 s | **≈ 13.4 s** |

Assumptions:
- 4a costs 0.4 s threaded (the GPU nodes get Kokkos-OpenMP threads too).
- 4b's factors are from §2.2.
- Kernels and everything else are unchanged.
- The GPU host scales at 0.84× (H100) or 0.82× (MIG) of `ov-omp` (§1.3).
- These are projections, not measurements. The port must measure P2 and P3.

### 3.3 Hardest gates
1. **4a claim 1 with iteration counts (C4).** Bitwise for the whole solve, in *both* modes. FMA anywhere in the cache changes Γ by ulps, which can change the iteration count and then Γ by up to 1e-3 relative. That is a declared exception: fast = strict here. The cost is nil, because the kernel is tiny.
2. **4a P4 / `N_min`.** About 26 launches per solve at 87 targets. The threaded launch latency must beat 0.164 ms of serial replay, or the solve falls back to serial. Measure on the r-test OLAF cases, the smallest of which have a few spanwise stations.
3. **4b A2/A4 and the deadness proofs.** Each of H2–H4 is legal only because no reader exists. The tests must enumerate the readers (card W2–W4) and include the glue with corrections, grid outputs on and off, and source panels.
4. A7 is trivially satisfied for unit 4 (no persistent state). **It is the gate that will decide the residency unit** (§2.3).

## 4. Corrections made elsewhere

- **olaf-placement.md:**
  - The §0 row "Lifting line and circulation solve" and §5 "What we do not fuse" called the circulation solve "small". It is small in data but the **largest single cost** on the tree cases (49–58% of the CPU run; 16–58% over all five cases).
  - §0 "Wind exchange" and §2 omitted AeroDyn's host evaluation at all wind points (5 wind + 3 tower per step).
  - F10's legality was R1 only; F10-B is class B.
  - §1's anatomy now points to §2.2 here for the AeroDyn wrapper.
  - Edited in place, marked "(corrected 2026-11-17, olaf-unit4.md)".
- **architecture.md §5** named the host pack "the GPU path's main Amdahl term". It is 0.3–0.75 s of 34–46 s. Edited in place with the measured split.

## 5. Audit additions (beyond olaf-placement §6 #1–14)

| # | Trap | Where | Action |
|---|---|---|---|
| 15 | `use, intrinsic :: ieee_arithmetic` in a hot per-point routine costs 14× (gfortran saves and restores the FP environment on each call); module-level `USE IEEE_ARITHMETIC` in `NWTC_IO` (`NWTC_IO.f90:25`) can do the same to procedures that see it | `IfW_FlowField.f90:386`; observed in the pin's `ui_leaves_nograd` | Grep every hot routine's scope for IEEE modules before timing it; C++ ports avoid it by construction |
| 16 | Static scans miss iteration-multiplied host loops: `kernel-inventory.py` rates `Wings_ComputeCirculationPolarData` "GLUE", wake 0, and `ui_quad_n1` HOST, yet together they are 58% of the run | `FVW_Wings.f90:314`, `FVW_BiotSavart.f90:513` | Profile before choosing a unit (olaf-placement §7 already says call frequency wins; this adds *convergence-loop* multiplicity) |
| 17 | Discrete decisions downstream of a reordered sum (convergence tests, iteration caps) turn ulp changes into ConvCrit-sized changes | `FVW_Wings.f90:467-471` | Kernels feeding such tests are class B in both modes |
| 18 | `r_wind` is packed at **maximum** size and always contains every grid point; consumers read slices | `FVW_Subs.f90:568-638` | Any change to wind evaluation must be slice-aware (card W3) |

## 6. Reproduce

| What | Command | Output |
|---|---|---|
| Sampler + buckets | `gcc -O2 -fPIC -shared -o libwsamp.so scripts/wsamp/wsamp.c -lrt`; `sbatch -J ofk-unit4-prof scripts/wsamp/prof.sbatch` (job 1923289); `scripts/wsamp/wsamp_report.py <run>/prof 500 --buckets [--csv label]` | `opt/build/u4-prof/run-1923289/<case>.<kf8\|ko8\|pin8>/` |
| Call counts | `tools/openfast-build.sh fomp opt/build/u4-fomp-pg -DCMAKE_Fortran_FLAGS=-pg -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg` (`OFK_TARGETS=aerodyn_driver`); `sbatch -J ofk-unit4-gprof scripts/wsamp/gprof.sbatch` (job 1923304) | `opt/build/u4-prof/gprof-1923304/<case>/{flat,gprof}.txt` |
| F10 and IEEE microbenchmarks | `sbatch -J ofk-unit4-mb scripts/unit4-mb/mb.sbatch` (job 1923596; builds with `scripts/unit4-mb/build.sh`) | `opt/build/u4-prof/ofk-unit4-mb.1923596.log` |
| GPU numbers used | jobs 1860779 (H100: `iv-gpu.1860779.log`), 1887850 (A100 MIG: `u3-gpu.ofk-treecode-mig.1887850.log`) | — |
