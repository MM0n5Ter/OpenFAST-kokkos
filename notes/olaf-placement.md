# OLAF: what runs on the CPU, what runs on the GPU, and what to fuse

Drafted 2026-11-14 from a read of `deps/openfast` @ `12fa3fc` (file:line references are to that pin) and from `scripts/kernel-inventory.py`.
This note **decides** placement; the script only **hints** (see §7).

## 0. Summary

| Group | Routines | Decision | Why |
|---|---|---|---|
| **Wake-sized kernels** | `ui_seg`, `ui_part_nograd`, `ui_tree_part`(+`ui_leaves_nograd`), `ui_tree_segment`, `ui_quad_src_nn` | **GPU, month 1** | O(N·M) or O(M log N) work over 10⁴–10⁶ wake elements. Over 90% of OLAF time at scale. |
| **Packing into kernel inputs** | `PackPanelsToSegments` → `LatticeToSegments` (+ mirror), `SegmentsToPart(Wrap)`, `PackConvectingPoints`, `UnPackInducedVelocity` | **GPU, month 1, fused** (F1, F2) | O(N) data movement executed 2–5× per step. On the host it forces a full host→device upload of the packed arrays on every call. |
| **Treecode build** | `grow_tree_part(_parallel/_substep)`, `compact_tree_part`, `cut_tree` | **Host first, GPU in week 3–4** (Lane B) | Pointer-based recursive octree. Phase 1 builds it on the host, flattens it, uploads it, and runs the traversal on the GPU. Then a Morton/LBVH build on the device. |
| **Wake time-stepping** | `FVW_CalcContStateDeriv` (means + dxdt), `FVW_Euler1`/`RK4` update, `PropagateWake`, `FakeGroundEffect`, `FVW_ContStates_Interp` | **Host in month 1; GPU in phase 2** (F3, F4) | Each one is O(N), cheap and memory-bound. They matter only once the kernels are on the GPU (Amdahl), and only if the wake state becomes device-resident. That changes who owns `x` (§4). |
| **Wind exchange** | `SetRequestedWindPoints`, `DistributeRequestedWind_*` | **Host, until InflowWind is ported** (Lane D) | The wind comes from InflowWind, which runs on the CPU. This is the one unavoidable D2H/H2D per step: 2 × 24 B × N. |
| **Lifting line and circulation solve** | `Wings_Panelling`, `Wings_ComputeCirculation(PolarData)`, `CirculationFromPolarData`, `ui_quad_n1`, `Map_LL_NW`, `Map_NW_FW`, `PlaceTipRoot`, UA wrappers | **Host, permanently for this phase** | O(nSpan) to O(nSpan²) with nSpan of about 10–100 per blade. Calls AirfoilInfo/UnsteadyAero (Fortran, CPU), uses ModMesh and has iterative convergence logic. Too small for a GPU. |
| **Setup, I/O, outputs** | `FVW_Init*`, `FVW_IO.f90`, `WriteVTKOutputs`, `FVW_AeroOuts`, `srcPnl_init/solve` (LAPACK), grid-output packing | **Host** | Run once, run rarely, or do I/O. |

**The finding that matters most for scheduling:** lifting-line induction is **not** cheap. The roadmap has it in Lane A week 3 ("block per CP"). This note sizes it as a peer of the wake kernel.
- `VelocityMethod(2)` is hard-wired to `idVelocityBasic` (`FVW_IO.f90:79`), so every lifting-line call is a full-wake pack plus a **direct** `ui_seg` over all N_seg segments × N_LL control points.
- There are 2–4 such calls per aero step (§1).
- With 10⁶ segments and 150 LL points, that is 3–6 × 10⁸ segment interactions per step, the same order as the wake treecode.
- The kernel has an awkward shape: few targets, many sources (§3).
- So the LL path must reach the GPU in the same milestone as `WakeInducedVelocities`, not after it.
- A block per CP alone (150 blocks) under-fills an H100. Split the sources across blocks too (§3).

## 1. Anatomy of one OLAF step (who calls the wake-sized work, and how often)

AeroDyn calls `FVW_UpdateStates` once per aero step (DTaero), then `FVW_CalcOutput`. The wake is integrated every DTfvw (≥ DTaero, `FVW.f90:598`).

```
FVW_UpdateStates (FVW.f90:567)
 ├─ Wings_Panelling, Map_LL_NW, Map_NW_FW                         host, O(nSpan)
 ├─ [UA] CalculateInputsAndOtherStatesForUA (:637) ── LL #a        (only if UA on)
 ├─ if ComputeWakeInduced (every DTfvw):
 │   ├─ DistributeRequestedWind_NWFW (:645)                        O(N) unpack of IfW wind
 │   ├─ FVW_Euler1 (:653)  or FVW_RK4 (:655, 4 stages)
 │   │    └─ FVW_CalcContStateDeriv (:801)
 │   │         ├─ WakeInducedVelocities (:855) ─────────── WAKE  ×1 (Euler) / ×4 (RK4)
 │   │         ├─ 4 mean-velocity reductions (:858-923)       O(N)
 │   │         ├─ dxdt = Vwnd + Vind (:951-954)               O(N)
 │   │         └─ Eps rates (:983-1007)                       O(N)
 │   │    └─ x += dt*dxdt (:1085-1092)                        O(N)
 │   └─ PropagateWake (:668) [and again on m%x1 if overcycling]   O(N) shift of r, Gamma, Eps, dxdt
 ├─ [overcycling] FVW_ContStates_Interp (:706)                     O(N)
 ├─ Wings_Panelling, Map_LL_NW, Map_NW_FW
 ├─ FVW_CalcConstrStateResidual (:727)
 │    └─ Wings_ComputeCirculationPolarData (FVW_Wings.f90:314)
 │         ├─ LiftingLineInducedVelocities(iDepthStart=iNWStart+1) (:391) ── LL #b
 │         └─ do while: ui_quad_n1 over nSpan² × iNWStart + polar lookup     host, small
 ├─ [UA] CalculateInputsAndOtherStatesForUA (:735) ── LL #c, UA_UpdateState, UA_SetGammaDyn
 ├─ FakeGroundEffect (:750)                                        O(N)
 └─ SetRequestedWindPoints (:753)                                  O(N) pack for IfW
FVW_CalcOutput (FVW.f90:1472)
 └─ CalcOutputForAD (:1426) → Wings_Panelling, LiftingLineInducedVelocities(iDepthStart=1) (:1457) ── LL #d
    └─ WriteVTKOutputs (:1551) → [if grid outputs] InducedVelocitiesAll_Init/OnGrid_Calc/End  (every DTvtk only)
```

| Path | Per aero step | Per call: build | Per call: evaluate |
|---|---|---|---|
| WAKE (`WakeInducedVelocities`, `FVW_Subs.f90:1379`) | 1 per DTfvw (×4 with RK4) | pack all segments → particles → tree (default `VelocityMethod(1)=2`, `FVW_IO.f90:74`) | tree traversal at M = free wake markers |
| LL (`LiftingLineInducedVelocities`, `FVW_Subs.f90:1488`) | **2** (b, d), **4** with UA (a, b, c, d) | pack all segments from `iDepthStart` (no particles, no tree) | direct `ui_seg`, M = nSpan·nWings (~150) |
| GRID (`InducedVelocitiesAll_OnGrid_Calc`, `:1133`) | 1 per DTvtk, only with grid outputs | same as WAKE | M = grid points |
| SRC panels (`srcPnl_ExtVelocities_OnPanels`, `:2081`) | every `nSrcPnlUpdate` steps | same as WAKE | M = number of panels |

Every path rebuilds its inputs from the lattice state `x%W(iW)%{r,Gamma,Eps}_{NW,FW}`, and nothing is reused between paths. That is the root of most fusion opportunities in §5.

**Sizes** (B blades/wings, S spanwise panels per blade, Nn/Nf near/far wake panels, `FWnSpan = 1`, `FVW_Subs.f90:59`):
- markers N_P ≈ B(S+1)(Nn+1) + 2B(Nf+1);
- segments N_S ≈ 2·B(S+1)(Nn+1) + 3B(Nf+1), doubled with `ShearModel = mirror`;
- particles = N_S · `PartPerSegment` (default 1);
- LL targets M_LL = B·S.

Example: B=3, S=50, Nn=1000 gives N_P ≈ 1.5×10⁵, N_S ≈ 3×10⁵ and M_LL = 150. One LL call is then 4.6×10⁷ segment interactions.

## 2. The CPU/GPU boundary for month 1, and what crosses it

Month 1 keeps the Fortran ownership of `x` (registry state). The glue copies it (`FVW_CopyContState`, predictor-corrector in the OpenFAST glue, `m%x1/x2` for overcycling) and writes it to VTK. The GPU sees it through the boundary below.

| Boundary call | H→D per call | D→H per call | Notes |
|---|---|---|---|
| WAKE | lattice `r` (3·N_P), `Gamma` and `Eps` (~4·N_P): **~56 B × N_P**; packing happens on the device | `Vind` at free markers: 24 B × N_P | The tree is built on the host in phase 1, so the flattened tree also goes H→D (~(10·3+8)·8 B per node plus the leaf buffer, 56 B × N_part). Phase 2 removes this. |
| LL | same lattice, ~56 B × N_P | 24 B × M_LL (tiny) | 2–4 calls per step with **the same far wake**. See F5 for upload-once. |
| Wind (IfW) | `V_wind` 24 B × N_P (host→device, phase 2 only) | positions 24 B × N_P | Unavoidable until IfW is on the device |

Order of magnitude: N_P = 10⁶ gives ~56 MB per upload, about 2–3 ms over PCIe Gen4/5. Treat this as the budget the kernels must amortize. The performance gate (criteria.md P3) reports transfer time separately from kernel time.

## 3. Kernel shapes → parallel pattern (one source, two backends)

| Shape | Where | Kokkos pattern | Order of the per-CP sum |
|---|---|---|---|
| **Square** (M ≈ N, both ≥ 10⁴) | WAKE direct (`VelocityMethod` 1/3), GRID, SRC | `RangePolicy` over targets; inner serial loop over sources (as in the Fortran OMP code); source tiles in team scratch on the GPU | **Same as Fortran** → bitwise on Serial/OpenMP; ≤ 1e-12 relative on CUDA (FMA, `sqrt`, `exp`) |
| **Skinny** (M ≈ 10², N ≥ 10⁵) | **LL** (`ui_seg`, 2–4× per step) | `TeamPolicy`: league over targets × chunks of sources, `TeamThreadRange` + `parallel_reduce`, then a deterministic second-stage sum over chunks. A target-only `RangePolicy` gives 150 threads on a GPU. | **Reordered** (tree sum per target). Deterministic run to run, but not Fortran order → tolerance class R1 (criteria.md) |
| **Tree traversal** (M targets, log-N work each, divergent) | WAKE default (`VelocityMethod` 2), LL if ever switched to tree | Flat node arrays with explicit stack and no recursion. Targets sorted by Morton key so a warp walks similar paths. Leaves use the packed `LeafBuf` layout that upstream already uses (`FVW_VortexTools.f90:644`) | Not compared with the Fortran tree bitwise (harness.md §1). The tree error must stay within the Fortran tree's error vs direct. |
| **Elementwise / shift / reduction** | pack/unpack, dxdt, Euler, Propagate, means | `RangePolicy` / `MDRangePolicy`; the means are `parallel_reduce` | Elementwise: bitwise. Means: reordered, R1. |

Rules for the "square" kernels, so that L0 can be bitwise on CPU:
1. One accumulator triple per target, initialized to 0. Sources are visited in packed order. The result is added to `Uind` once ("Side effects!!!": accumulate, never overwrite).
2. Hoist `select case (RegFunction)` out of the loop (one kernel instantiation per regularization, as the Fortran does by copy-paste). Template on the regularization functor.
3. Keep every guard exactly: `denominator > PRECISION_UI`, `norm2_orth > PRECISION_UI`, `norm2_r0 > PRECISION_UI`, `exp_value < MIN_EXP_VALUE`, the particle `r3 > PART_REG_CUT3*rc3` skip, and the compact-kernel `merge` + clamps. See the audit list in §6.

## 4. Who owns the wake state (the residency decision)

- **Phase 1 (month 1): Fortran owns `x`; the device holds a persistent mirror.**
  - Upload at the boundary (§2), with buffers allocated once and grown by doubling (`nNWMax`/`nFWMax` give the final size at init, so allocate once in `FVW_Init`).
  - API risk: none. Every Fortran consumer of `x` keeps working.
- **Phase 2: the device owns the wake lattice between steps, and Fortran gets a lazily synced mirror.**
  - The whole `ComputeWakeInduced` branch (CCSD + Euler/RK4 + Propagate + FakeGroundEffect) runs on the device. Only these cross:
    - the first `iNWStart+1` panels written by `Map_LL_NW`/`Map_NW_FW` (host → device, O(nSpan));
    - the wind (§2);
    - the full lattice, only when someone reads it on the host: VTK output, `FVW_CopyContState` from the glue, `FVW_End`, linearization.
  - It needs a "host dirty / device dirty" flag in the unit's internal handle and a sync call at every Fortran entry point that reads `x`.
  - That list of entry points comes from `scripts/callers.sh` plus a grep for `x%W` outside the unit. It is the API-gate item for phase 2.

## 5. Fusion plan

Each candidate is classified by its **legality**:
- **B**: bitwise-preserving on CPU backends (same arithmetic, same order);
- **R**: reorders a floating-point sum (tolerance class R1);
- **A**: algorithmic change (different but valid approximation, tolerance class R2, opt-in only).

| ID | Fuse what | Where today | Gain | Legality | Phase |
|---|---|---|---|---|---|
| **F1** | `LatticeToSegments` (NW + FW, all wings) + ground mirror + `SegmentsToPart` → **one pack kernel** writing segments and particles directly | `PackPanelsToSegments` `FVW_Subs.f90:913`, `LatticeToSegments` `FVW_VortexTools.f90:283`, `SegmentsToPart` `:396` | Removes 3 passes and 2 intermediate arrays per call; keeps packing on the device so only the lattice is uploaded. Blocker: the `iHead*` running counters are serial. The offsets have a **closed form** (per wing/depth/span; `CountSegments` already computes them), so each thread computes its own output index. | B (the same arithmetic per element; output order identical to the Fortran pack) | 1 |
| **F2** | `PackConvectingPoints` + kernel + `UnPackInducedVelocity` → the kernel reads targets **from the lattice** and writes `Vind_NW/FW` in lattice layout | `FVW_Subs.f90:1416`, `:1448` | Removes 2 O(N) copies. The free-wake targets *are* the lattice points `r_NW(:, :, 1:nNWEff+1)`, `r_FW(:, :, 1:nFWEff+1)`, so target i maps to (wing, span, depth) in closed form. | B | 1 |
| **F3** | `CalcContStateDeriv` means + dxdt + `Euler1` update + `Eps` update → 1 reduction kernel (the 4 means) + 1 elementwise kernel | `FVW.f90:858-1007`, `:1085-1092` | 8+ full-array passes become 2. Needs phase-2 residency to matter. | R for the 4 means (`VmeanNWFree`, `VmeanNW`, `VmeanNWFixed`, `VmeanFWFree`); B elsewhere. `VmeanNWFixed` is debug-only (`:892` "TODO REMOVE FOR DEBUG ONLY") but it is computed, so compute it too. | 2 |
| **F4** | `PropagateWake`: 8 array shifts (r/Gamma/Eps × NW/FW, plus `dxdt`) → 1 launch with double buffers | `FVW_Subs.f90:325` | 8 launches become 1. The in-place descending shift has a race on the GPU, so write to the second buffer and swap. A ring buffer would remove the copy entirely, but it changes the layout of `x` that callers see (API), so it waits for full device ownership. | B | 2 |
| **F5** | **LL upload-once and far-wake reuse.** Within one step, LL #b/#c/#d see the same far wake. Only the first `iNWStart` panels change (`Map_LL_NW`), plus `FakeGroundEffect`. | `FVW_Wings.f90:391`, `FVW.f90:1457`, `:1749` | (i) Upload the lattice once per step, then patch the first panels: B, needs a version counter bumped by every writer of `x`. (ii) Cache `Uind_LL(far wake)` and add the near-panel part per call: saves 1–3 of the 2–4 full LL sums per step. | (i) B. (ii) R: splits the sum. | (i) 1, (ii) 2 |
| **F6** | `SetRequestedWindPoints` pack + D2H, and `DistributeRequestedWind_NWFW` unpack + H2D → one strided copy each | `FVW_Subs.f90:568`, `:668` | Halves the boundary traffic | B | 2 (with Lane D) |
| **F7** | WAKE + GRID + SRC evaluation in one launch when they share a step (same sources, concatenated targets) | `FVW.f90:1628-1640` | One tree and segment build instead of 2–3 on output steps | B for direct; tree as in §3 | 2 |
| **F8** | Fuse `ui_quad_src_nn` into the same target loop as the vortex kernel (one read of the target, two source sets) | `InducedVelocitiesAll_Calc` `FVW_Subs.f90:1323` | One pass over targets instead of two | R (two partial sums become one) unless kept as two accumulators added in the Fortran order: then B | 1 (cheap) |
| **F9** | Euler update **and** `PropagateWake` shift as one out-of-place pass: `x'(:,:,a) = x(:,:,a-1) + dt·dxdt(:,:,a-1)` for the rows Euler touched (a−1 ≤ nNW+1; `Eps` to nNW), and a plain copy for the rest. Same for `dxdt`. | `FVW.f90:1085-1092` then `:665-668` → `FVW_Subs.f90:325`. Only `PrepareNextTimeStep` (counter bumps) runs in between. | Two full read+write passes over r/Eps become one, and the in-place descending shift (a race on a GPU) goes away. Euler1 only; RK4 keeps its stages. | B (identical single operation per element) | 2 (with F3/F4) |
| **F10** | **Hoist the geometry out of the circulation iteration** (host). `Vvar` is linear in Γ (`ui_seg_11`: `Kv = SegGamma·…`) and the ring geometry and `Eps` are fixed inside the `do while`. Build the influence matrix A (3·nCP × nCP·iNWStart) once per solve, then `Vvar = A·GammaLastIter` per iteration. | `FVW_Wings.f90:421-469` (`ui_quad_n1` per span × depth × iteration, `FVW_BiotSavart.f90:513`) | Up to `CircSolvMaxIter` = 30 (`FVW_IO.f90:47`) geometry evaluations per solve become 1, plus 30 small mat-vecs. Note: `ui_quad_n1`'s directives are written `!OMP`, not `!$OMP` (`FVW_BiotSavart.f90:527-528`), so today this loop is **serial**. | R1 (Γ is multiplied after the kernel sum instead of inside it) | 1 (host-only change; can go in before any GPU work) |

What we do **not** fuse:
- Anything with the circulation solve loop. It is host, small, and couples to AirfoilInfo/UA.
- Tree build with evaluation: the tree is built once and evaluated at M targets.
- Across `FVW_UpdateStates` and `FVW_CalcOutput`: the glue may run other modules in between and may call `CalcOutput` more than once (predictor-corrector, Jacobians in `AD_JacobianPInput`). State must be re-validated at every API entry.

## 6. Audit list: traps found while reading (each needs an L0 case or a note in the unit card)

| # | Trap | Where | What the port must do |
|---|---|---|---|
| 1 | `MINNORM = 1e-4` is a default-kind literal, so **its value depends on the compile flags**. OpenFAST's double build uses `-fdefault-real-8 -fdefault-double-8` (`cmake/OpenfastFortranOptions.cmake:138-140`; Intel and LLVM Flang have equivalents), so there it is exactly `1e-4`. In the single build it is `1e-4` as a single. *Corrected 2026-11-14:* this row first said `9.9999997e-05`, true only without that flag. Hard-coding `double(1e-4f)` fails L0 at relative error 1.0 (found by the `biotsavart-direct` port). | `FVW_BiotSavart.f90:20`, `:377` | **Pass flag-dependent constants from the Fortran shim**, which is correct under any flags. Build the oracle with the pin's exact flags. L0 straddles the threshold. |
| 2 | `PRECISION_UI = epsilon(1.0_ReKi)/100` depends on precision | `:12` | `std::numeric_limits<ReKi>::epsilon()/100`, picked up from `OPENFAST_DOUBLE_PRECISION` |
| 3 | `MINDENOM = 0` is added to the denominator (a no-op, but keep the add for bitwise CPU) | `:18`, `:192` | Keep it as is |
| 4 | `SegConnct` holds **1-based point indices**; rows 3/4 are depth/span; the mirror adds `nSegP` | `FVW_Subs.f90:980` | Subtract 1 at use, never in storage. Rows 3/4 are used by the segment tree. |
| 5 | `Uind` is accumulated ("Side effects!!!"); the callers prefill. `m%Uind(:,nCPs+1:)=1000` is a debug sentinel. | `FVW_Subs.f90:1403-1404` | Touch only `[0,nCPs)`. L0 checks sentinel tails. |
| 6 | `ui_seg` regularization branches differ in **guard style**: Vatistas uses `cycle` (`<=`), the others nested `if` (`>`). They are equivalent for finite inputs but **not for NaN**: Vatistas lets a NaN through to `Uind`, while the nested `>` style skips it. Keep each one literally. | `FVW_BiotSavart.f90:299-306` | Template per branch; don't "unify" guards |
| 7 | Particle Exp kernel skips `exp` when `r3 > PART_REG_CUT3*rc3`; Compact uses `merge` with both arms evaluated plus clamps | `:415`, `:437-444` | Same branch structure. The clamps exist to avoid 0/0 in the discarded arm, which matters on GPUs too (NaN × 0 = NaN). |
| 8 | The order-2 multipole term keeps a loop form: "BAR_OLAF test fails (very slightly) with unrolled version" | `FVW_VortexTools.f90:1612` | Evidence that r-test is sensitive at this level. Keep the term order, and expect the tree path to need tolerance-level comparison. |
| 9 | `LatticeToSegments` **includes or omits** shed segments depending on `bShedVorticity`, `bShedLastVorticity` and `LastNWShed`; the FW uses the Eps hack (`bHackEpsilon`, tip/root epsilon in rows 2/3) | `FVW_VortexTools.f90:283-394`, `FVW_Subs.f90:872` | The closed-form offsets in F1 must reproduce `CountSegments` for all 4 flag combinations. L0 enumerates them. |
| 10 | `PropagateWake` nullifies with `-999.9` and the code checks `< -999` in places | `FVW_Subs.f90:325-400` | Preserve sentinels in any buffer the host can read |
| 11 | OLAF+OpenMP is non-deterministic (#3401/#3402); baselines use `OMP_NUM_THREADS=1`. Also, with OpenMP on, `!$OMP SIMD reduction` (`FVW_BiotSavart.f90:385,430`) **reorders the Fortran's own particle sums**: about 1.5e-15 against the OpenMP-off build. | testing.md | The class-B oracle is the pin built with **OPENMP=OFF** and its exact flags. Kokkos-Serial/OpenMP L0 compares against that. |
| 12 | `T_Tree` holds **pointers into `Part`** (`Tree%Part%P => PartP(:,1:nPart)`) | `FVW_VortexTools.f90:577-579` | The device tree owns its own copies; never alias Fortran memory from C++ beyond one call |
| 13 | More default-kind literals (`scripts/literal-kinds.py`, 12 hits on OLAF; like #1, each is **double in the double build** because of `-fdefault-real-8`):<br>- `CoreSpreadAlpha = 1.25643` sets the core-growth rate `visc_fact` of every wake panel, every step;<br>- the degenerate-domain threshold `1.0e-5` of both tree builds;<br>- the `PlaceTipRoot` constants `0.15 / 0.85 / 3.14 / 1.3 / 1.7 / 1e-6`, which are on the host. | `FVW_Subs.f90:56` → `FVW.f90:992`; `FVW_VortexTools.f90:814, :1179`; `:152-178` | Take them from the Fortran (shim argument or a `ReKi` value), not as `double(xf)`. Re-run the script on every new file. |
| 14 | `UpdateStates` can be called **twice at the same `t`** (`bReevaluation` → `RollBackPreviousTimeStep`, `FVW.f90:597-610`), and the glue copies `x` (`FVW_CopyContState`, `:649`, `:680`) | `FVW.f90` | Phase-2 residency must handle rollback and copies (criteria.md A7). Phase 1 is unaffected because Fortran owns `x`. |

## 7. How `scripts/kernel-inventory.py` relates to this note

`scripts/kernel-inventory.py [--md] [--json out] [files]` scans every routine and prints:
- size and loop nest;
- **wake-scale loops**: loops bounded by nNW/nFW/nSeg/nCPs/nPart…, as opposed to nSpan/nWings;
- OpenMP directives and allocations;
- `iHead` running counters (serial pack pattern), recursion and pointer trees.

It then gives a hint: `DEVICE`, `DEVICE*` (device after removing a blocker), `GLUE` or `HOST`. On OLAF: 139 routines; 14 DEVICE, 12 DEVICE*, 11 GLUE, 102 HOST.

Where the hint and this note disagree, this note wins, and the reason is call frequency, which a static scan cannot see:
- `curl_regular_grid`, `EllipsoidPanels` and `ReadAndInterpGamma` are DEVICE by shape but run at init or for output only → host.
- `ui_quad_n1` has a wake-looking loop (`nCPs`), but its only caller passes the nSpan lifting-line points inside the circulation iteration (`FVW_Wings.f90:441`) → host.

Use the script on the next module first, then write that module's version of this note.
