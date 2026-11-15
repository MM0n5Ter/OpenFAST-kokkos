# treecode

The OLAF particle treecode (`VelocityMethod(1) = 2`, OLAF's default wake path): `grow_tree_part` and `ui_tree_part` (`FVW_VortexTools.f90`) ported to C++/Kokkos behind the frozen Fortran interfaces of `InducedVelocitiesAll_Init/Calc` and `LiftingLineInducedVelocities`. It lives in unit 2's per-instance handle. Unit 3 of 3 (architecture.md §7).

- **Owner**: agent session `199ad883` (subagent, unit 3). Date: 2026-11-14.
- **Toolchain**: Kokkos 5.3.0 (C++20), gcc/gfortran 12.3.0, cmake 3.31.12, CUDA 12.3, OpenBLAS 0.3.28.

## Status / resume (checkpoint 2026-11-14 ~21:25, requested by the parent: repo layout migration)

### Done, with evidence
All CPU evidence comes from the login node (Xeon 8480+) unless a job is named.
- **Code.** The strict and fast tree is in the fork; see the file list below. It compiles on OpenMP (`opt/build/u3-omp`, `u3-iv-omp`, `u3-bsd-omp`) and on CUDA sm_80 and sm_90 (`u3-cuda-a100`, `u3-iv-cuda-a100`, `u3-cuda-h100`).
- **L0 `l0_tc` (`opt/build/u3-omp`): 951/951 PASS, seed 1** (40 random cases plus the fixed set, including 20000-particle clouds on the parallel build path); 959/959 on seed 7 (80 random cases).
  - Claim 1, strict vs the Fortran tree (pin flags, OPENMP=OFF): **bitwise in every case**, 0 cases with bit differences. This includes -0, NaN tails, ld=4 CPs, and every RegFunctionPart / BranchFactor / DistanceDirect / cloud kind.
  - R2, fast vs direct, ≤ max(1.1 × Fortran-tree vs direct, 1e-12): **PASS**. Worst ratio 0.909, i.e. E(fast) = E(Fortran) in every case, because the tree is the same tree.
  - 5 Fortran-singular cases are excluded and reproduced (trap T4 below).
  - Claim 3 (information only), fast vs strict: worst 5.4e-14 of Σ|c|; acceptance flips 0.
  - P6 at L0: 0 allocations in 20 steady grow+eval steps.
- **R7 `r7_threads.sh`: PASS.** The strict and fast dumps are md5-identical across processes and across 1/4/8 threads.
- **Unit 2 `l0_iv` with the tree ported (`opt/build/u3-iv-omp`): 13/13 PASS.** It goes through the real `InducedVelocitiesAll_Init/Calc` and `LiftingLineInducedVelocities` shims.
  - TreePart in strict: 380 comparisons per run, bitwise against the oracle, packed and lattice runs both.
  - Strict vs fast: 1.9e-15.
  - The A7 sequence now also runs with wake VelocityMethod 2: 544 comparisons, bitwise against fresh instances and the oracle.
  - P6 with VelocityMethod 2: 0 allocations.
- **First CPU sweep (job 1862730, med, 8 threads).** It ran the *serial*-build binary, now superseded. Fortran+OpenMP vs fast (build / eval / total, s):

  | N | Fortran+OpenMP | fast |
  |---|---|---|
  | 1e5 | 0.0109 / 0.0256 / 0.0365 | 0.0157 / 0.0225 / 0.0382 |
  | 1e6 | 0.156 / 0.253 / 0.410 | 0.221 / 0.221 / 0.443 |

  The fast evaluation is 1.14× the Fortran, but the build was serial. The parallel build is in the current binary, re-measured by job 1862863.

- **Results that arrived at the checkpoint (21:25).**
  - **Job 1862863: CPU sweep, med (Xeon 8480+), 8 threads, current binary with the parallel host build.** Total s per call (build + eval), best of 3:

    | N | Fortran+OpenMP | Kokkos fast | Kokkos strict | fast speedup |
    |---|---|---|---|---|
    | 1e3 | 4.25e-4 | 2.69e-4 | 3.79e-4 | 1.58× |
    | 3e3 | 1.67e-3 | 1.02e-3 | 1.41e-3 | 1.64× |
    | 1e4 | 4.04e-3 | 3.05e-3 | 3.90e-3 | 1.33× |
    | 3e4 | 1.17e-2 | 9.82e-3 | 1.24e-2 | 1.19× |
    | 1e5 | 3.91e-2 | 3.31e-2 | 4.23e-2 | 1.18× |
    | 3e5 | 1.21e-1 | 1.02e-1 | 1.36e-1 | 1.18× |
    | 1e6 | 4.55e-1 | 3.65e-1 | 4.30e-1 | 1.25× |

    At 1e6 the split is: Fortran build 0.195 / eval 0.260; fast build 0.131 / eval 0.234. Strict is at 0.94–1.06× the Fortran. P1/P2 at kernel level: fast ≥ 1.18× at every N.
  - **Job 1862865: A100 MIG, `u3-cuda-a100`.**
    - l0_tc 3/3. **Claim 2:** strict CUDA vs Fortran worst 1.8e-16 of Σ|c| (bound 1e-13).
    - R2 PASS on the GPU (worst ratio 0.909, 5 singular cases reproduced).
    - Claim 3 on the GPU: 4.6e-14.
    - P6 on the GPU: 0 allocations in steady steps.
    - R7 on the GPU: md5-identical across processes.
    - Unit 2 `l0_iv` on the GPU with the tree ported: 12/12; TreePart strict vs oracle 4.4e-17, fast vs strict 1.9e-15; A7 (544) bitwise.
- **`opt/build/fork3-fomp` built completely** (aerodyn_driver, openfast, regression_tests).

### Half-done
- **The fast-mode device build.** Declarations only: `DevBuild<T>` and `fast::build_device` in `treecode_common.hpp`. **Nothing is implemented and nothing calls them.** Today fast mode on CUDA builds the tree on the host (in parallel) and uploads the flat tree. The design is below ("Next steps" 3).
- **Fork builds.**
  - `opt/build/fork3-fomp` (Kokkos-OpenMP + Fortran OPENMP=ON) is fully built from the pre-migration fork sources.
  - `opt/build/fork3-off` was **not started**: the build chain broke when I edited `build_u3.sh` while it ran.
  - `fork3-cuda-a100` and `fork3-cuda-h100` were not started.
- **Not done:** system level (R6, R8, claim 1 on the BAR case, P2, P6 on a full run), A1/A2/A5 on the fork build, GPU evidence, mutation tests, the segment tree.

### Slurm jobs (logs in `opt/build/`)
| Job | What | Log / output | Look for |
|---|---|---|---|
| 1862730 `ofk-treecode-sweep` (done) | CPU sweep, serial-build binary | `u3-cpu.ofk-treecode-sweep.1862730.log`, `u3-sweep-1862730/` | superseded by 1862863 |
| 1862863 `ofk-treecode-sweep2` (med, **done**, results above) | P1/P5 CPU sweep N=1e3..1e6, Fortran+OpenMP 8 threads vs Kokkos-OpenMP fast/strict, with the per-phase profile | `u3-cpu.ofk-treecode-sweep2.1862863.log`, `u3-sweep-1862863/sweep.csv` | `csv,` rows; fast/Fortran at 1e5, 1e6 |
| 1862864 `ofk-treecode-l0` (med) | CPU L0 of units 3, 2 (tree ported) and 1 (`u3-bsd-omp`, unit 1 must stay 13/13) | `u3-cpu.ofk-treecode-l0.1862864.log`, `u3-l0-1862864/*.ctest.log` | `tests passed` lines |
| 1862865 `ofk-treecode-mig` (interactive_gpu, A100 MIG, **done**, results above) | claim 2 (strict CUDA vs Fortran, `OFK_L0_RTOL=1e-13`), R2/claim 3/P6 on the GPU, R7 GPU, unit 2 `l0_iv` on the GPU (`u3-iv-cuda-a100`). The full-case GPU block is skipped: `fork3-cuda-a100` does not exist yet. | `u3-gpu.ofk-treecode-mig.1862865.log`, `u3-mig-a100-1862865/` | `claim 1/2 ... worst err`, `tests passed` |
| 1862690 `ofk-treecode-perf` (a100 full GPU, pending) | GPU sweep, port only, fast and strict, profile split | `u3-gpu.ofk-treecode-perf.1862690.log`, `u3-perf-a100-1862690/sweep.csv` | `a100,csv,` rows; `report` lines (build / treeH2D / cpH2D / kernel / D2H) |
| 1862732 `ofk-treecode-perf` (h100,h200, pending; queue estimate 2026-11-16) | same on sm_90 | `u3-gpu.ofk-treecode-perf.1862732.log` | same |

The GPU perf jobs copy `perf_tc` from `opt/build/u3-cuda-{a100,h100}` **when they start**. After the migration, rebuild those trees from the new layout before the jobs start, or cancel and resubmit. GPU scripts already prepend the portable OpenBLAS (`opt/openblas-0.3.28-portable/lib`); the module's OpenBLAS dies with SIGILL on the EPYC GPU nodes.

### Files created or modified in `fork/` (all uncommitted)
| File | Change |
|---|---|
| `fork/modules/kokkos/src/treecode_common.hpp` | **new**: flat-tree types (preorder plus skip), build fragments/tasks, `TreeConsts`, `EvalArgs`; `DevBuild` and `fast::build_device` (declared, not implemented) |
| `fork/modules/kokkos/src/treecode_impl.hpp` | **new**: the host build (Fortran-faithful, parallel subtrees spliced in preorder), the multipole (statement by statement), the stackless traversal; it reuses unit 1's `part_term` for leaves |
| `fork/modules/kokkos/src/treecode_strict.cpp`, `treecode_fast.cpp` | **new**: the two mode TUs (6 lines each) |
| `fork/modules/kokkos/src/induced_velocities.cpp` | **modified** (unit 2's handle): `kMethodTree = 2`; tree state in `Handle<T>` (`HostTree`, device `dnr/dmom/dni/dleaf`, `treeValid/treeMode/treeThr`); `sources_packed` accepts method 2 (real RegParam, host staging only on device backends); `ensureR`, `tree_build`, `flat_tree`; `eval` gains bf/dd and the tree branch; C API `ofk_iv_tree_build`, `ofk_iv_eval_tree`, `ofk_tc_enabled`; config `OFK_TC`; tree fields in the profile line; `sources_valid` / `release_sources` know the tree |
| `fork/modules/kokkos/src/induced_velocities.hpp` | **modified**: declarations of the three new C functions |
| `fork/modules/kokkos/CMakeLists.txt` | **modified**: `treecode_strict.cpp` added to `OFK_STRICT_SOURCES`, `treecode_fast.cpp` to `OFK_FAST_SOURCES`, header comment |
| `fork/modules/aerodyn/src/FVW_Subs.f90` | **modified** (+80 lines, all under `#ifdef OPENFAST_KOKKOS`): bind(C) interfaces `ofk_iv_tree_build`, `ofk_iv_eval_tree`, `ofk_tc_enabled`; private helpers `ofk_tc_grow` (replaces `grow_tree_part`) and `ofk_tc_ui` (replaces `ui_tree_part`); shims in `InducedVelocitiesAll_Init` (TreePart), `InducedVelocitiesAll_Calc` (TreePart) and `LiftingLineInducedVelocities` (TreePart, with DistanceDirect). Each one falls back to the original Fortran when `OFK_TC=fortran`. |

**Outside `fork/`:**
- `tests/treecode/*` is new: `CMakeLists.txt`, `tc_api.f90`, `l0_tc.f90`, `perf_tc.f90`, `r7_threads.sh`, `make_cases.sh`, `build_u3.sh`, `cpu.sbatch`, `gpu.sbatch`, this card.
- `tests/induced-velocities/l0_iv.f90` is **adapted**:
  - `scale_sum` gives method 2 the particle Σ|c| scale. Before, method 2 was "exact", because both sides ran the same Fortran.
  - `a7_sequence(velW)` is called for wake methods 1 and 2.
- `fork/FVW_VortexTools.f90` is **not** modified: `grow_tree_part` and `ui_tree_part` stay public Fortran for direct callers such as `FVW_Tests`.

### Exact next steps (on the new layout)
1. Rebuild the test trees from the new layout: `u3-omp`, `u3-iv-omp`, `u3-bsd-omp`, `u3-cuda-a100`, `u3-iv-cuda-a100`, `u3-cuda-h100` (see `build_u3.sh`; copy the script before running it). Rerun `ctest` in `u3-omp` and `u3-iv-omp`; expect the results above. Read jobs 1862863, 1862864 and 1862865.
2. Build the fork trees `fork3-off` (OPENFAST_KOKKOS=OFF), `fork3-fomp` (Kokkos-OpenMP + OPENMP=ON) and `fork3-cuda-{a100,h100}`. Then:
   - `sbatch -J ofk-treecode-sys tests/treecode/cpu.sbatch sys` runs:
     - R6;
     - claim 1 at system level (strict text output = fork3-off text output on `iv-cases/bar-r8-tree`);
     - R8 against base₁ `iv-r8-1860700/bar-r8-tree.fomp1`, with the FMA envelope `iv-r8ffma-1860993/bar-r8-tree.ffma`;
     - P2 on `bar-r8-tree` and `u3-cases/bar-perf-tree` (fork3-fomp vs pin+OpenMP, 8 threads, plus an `OFK_TC=fortran` A/B);
     - P6 full run.
   - Rerun the `mig` job once `fork3-cuda-a100` exists: GPU full case, R8 on the GPU, R7 on two GPU processes.
   - **Note for R8:** base₁ is the pin with OPENMP=ON, whose `ui_leaves_nograd` `!$OMP SIMD` reductions reorder the leaf sums. So strict (= OPENMP=OFF Fortran) will **not** be identical to base₁. Claim 1 at system level is strict vs fork3-off. The OFF-vs-OMP drift is itself a natural envelope; report it.
3. **The fast-mode device build** (design decided, not coded): level-synchronous on the device, one TeamPolicy league per level over the active nodes.
   - K1, a team per node: a custom team reduction of {Σw·p, Σw, maxReg, min/max xyz, cnt[8]}; then the dead test; then a team Sum over the 30 moments; then status, written as a BFS record. K1 also writes, per particle position, nodeOf, level and octant.
   - K2: a scan of nChild over the level's nodes gives the child slots and nNext. This is one scalar D2H per level.
   - K3: a global `parallel_scan` with a struct of 8 ints (per-octant exclusive rank) over particle positions in split nodes, plus the node base. Integers, so exact.
   - K4: a stable scatter `perm → perm2`.
   - K5: child records, G_child = G + (0.5·R/2)·(±1) and R/2. This is exactly the Fortran's float recursion, so the octant decisions and the topology equal the Fortran's (min/max/counts are exact). Only centers and moments round differently.
   - After the levels: a bottom-up pass (subtree size, subtree leaves) and a top-down pass (preorder index, leaf start), then a write into the same flat preorder/skip arrays plus the leaf gather.
   - No thrust sort: Kokkos 5.3's `sort_by_key` uses thrust temporaries, which are per-call device allocations (P6), and keys of 3 bits × 21 levels would limit the depth.
   - n ≤ 1 stays on the host. Use it only on device backends; host backends keep the parallel host build, which is a loop-shape choice per backend.
   - Gate it with l0_tc's R2 and with claim 3 on the GPU. Then resubmit the GPU perf jobs.
4. Then mutation checks:
   - octant `>` → `>=`;
   - powi chain r^5 = r^4·r;
   - `pow` → `sqrt` in the strict build;
   - leaves before branches swapped;
   - degenerate threshold hard-coded as `1e-5f`.
   Expect l0_tc claim 1 to fail for each. Then A1/A2/A5/A7 on fork3-*, then E1–E6.
5. The segment tree (`grow_tree_segment`, `ui_tree_segment`) was not started; see "Segment tree: what it needs" below.

## Contract (frozen, pin `12fa3fc`)

```fortran
subroutine grow_tree_part(Tree, nPart, PartP, PartAlpha, PartRegFunction, PartRegParam, iStep)    ! FVW_VortexTools.f90:557
   type(T_Tree), intent(inout), target :: Tree     ! Tree%Part%P => PartP(:,1:nPart): POINTERS into Part (audit 12)
subroutine ui_tree_part(Tree, icp_end, CPs, BranchFactor, DistanceDirect, Uind, ErrStat, ErrMsg)  ! :1512
   real(ReKi), dimension(:,:), intent(inout) :: Uind  ! ACCUMULATED: Uind(:,icp) = Uind(:,icp) + Uind_tmp
```

Every in-tree caller holds `T_Tree` as a **local** variable and passes it only to `InducedVelocitiesAll_Init/Calc/End` or to the LL tree case:
- `FVW_Subs.f90:1389`, `:1502`, `:2090` (pin);
- `FVW.f90:1572`.

So the port can leave the Fortran tree un-grown in Kokkos builds. `cut_tree` on an un-grown tree is a no-op. `FVW_Tests` calls `grow_tree_part`/`ui_tree_part` directly; those stay Fortran.

Entry points and the handle (unit 2 lesson 1): the key is `m%Sgmt`, which `Init`, `Calc` and the LL path all receive.

| Entry point | What the port does |
|---|---|
| `InducedVelocitiesAll_Init`, TreePart | The Fortran `SegmentsToPartWrap` still runs: `Part` is an API output. Then `ofk_tc_grow` builds the particle source set from `Part` plus the tree. |
| `InducedVelocitiesAll_Calc`, TreePart | `ofk_tc_ui(…, TreeBranchFactor(1), DistanceDirect = 0)`. If the handle has no tree for `(method 2, Part%nAct)` it grows one first (Calc without Init). |
| `LiftingLineInducedVelocities`, TreePart (`VELOCITYMETHODLL=2`, advanced option) | grow, then `ofk_tc_ui(…, TreeBranchFactor(2), DistanceDirect = 2.2·MaxWingLength)` |
| `InducedVelocitiesAll_End` | unchanged (`cut_tree` on the un-grown local tree) |

## Design

### Strict: the Fortran tree, bit for bit, in GPU-legal form
- **Build (host).** A depth-first recursion over index ranges, doing per node exactly what `grow_tree_part_substep` does, in `node%iPart` order (increasing particle index, kept by a stable counting partition):
  - the weight norm `pow(s, 0.5)` (gfortran emits `call pow` for `**0.5_ReKi`);
  - the barycenter;
  - `EqualRealNos` dead test (tolerance 50·ε);
  - moments in statement order;
  - the leaf rule `R ≤ maxReg ∨ n ≤ 32`;
  - octants from the geometric center with `>`;
  - the degenerate-extent rule (the threshold `1.0e-5` is passed from the shim as `real(1.0e-5, ReKi)`);
  - singletons as the node's leaves, in octant order;
  - children in octant order, with `G + (0.5·halfSize)·(±1)` and `R/2`;
  - root R = max extent × `1.001_ReKi`.
- **Parallel subtrees.** Like the Fortran's `grow_tree_part_parallel`, subtrees below depth 2 are built in parallel (n ≥ 8192) and spliced in preorder. Every node's arithmetic depends only on its own particles, so the result is the serial result. That is checked bitwise by L0 (20000-particle clouds) and by R7 (1/4/8 threads).
- **Does the Fortran build depend on thread count or OpenMP task order? No.** Each subtree is independent and per-node order is fixed. Strict is defined against the OPENMP=OFF build only because, with `-fopenmp`, the leaf kernel's `!$OMP SIMD reduction` reorders the Fortran's own leaf sums (audit 11).
- **Flat tree.**
  - Nodes are in depth-first **preorder**, with `skip(i)` = the first node after i's subtree. Dead nodes are dropped: they do no arithmetic in the Fortran.
  - Per node, node-major: center, radius, maxRegParam (5); Moments (30); leafStart, leafN, skip (3).
  - Leaves are packed in preorder, SoA with 7 columns, as `compact_tree_part` does.
- **Evaluation.**
  - One thread per target. A stackless loop: `i = open ? i+1 : skip(i)`. It visits nodes in exactly `ui_tree_part_11`'s recursion order.
  - MAC: `max(bf·R + f·maxReg, dd)`, `r = sqrt((dx²+dy²)+dz²)`, open if `r < dist`.
  - Open: the node's leaf bucket is summed in its own accumulator (unit 1's `part_term`, the same body as `ui_leaves_nograd`), then added. Far: the multipole, statement by statement:
    - the order-2 loop form is kept (audit 8);
    - `r**3/5/7` use gfortran's powi chains `r·r²`, `r²·r³`, `r³·r⁴` (from the assembly);
    - `3·fourpi_inv` and `7.5·fourpi_inv` are constant-folded, as gfortran does.
  - The target adds `Uind + Uind_tmp` once, also when the tree is empty (`-0 + +0`).
- On CUDA, strict uses the same host-built tree. Only `exp` and the device `pow` (build) differ.

### Fast (class R2): the same tree, FMA, vectorized leaves
- The **same topology rules**. Today the fast mode reuses the host build compiled with FMA and `sqrt` for the norm. It is the same tree up to the `EqualRealNos` edge.
- Fast-mode evaluation on host backends: each target's leaf sums go to 8 fixed lanes across all of its buckets, with `omp simd` and a fixed pairwise tree at the end. On device backends: 1 lane.
- **Why this design.** R2 asks the fast tree's error vs direct to be ≤ 1.1 × the Fortran tree's. With the same tree, E(fast) = E(Fortran) up to rounding: L0 worst ratio 0.909 = 1/1.1. So R2 holds by construction rather than by luck on each input. A different tree (Morton/LBVH, another leaf size or MAC) would have needed a tighter MAC to be safe on every input, at a cost in speed.
- The planned device build (next step 3) keeps this property: octant decisions use the Fortran's float recursion for the geometric centers, so the topology is the Fortran's.
- The GPU speed comes from the evaluation (thread per target, no stack) and, once built, from the device build. The host build is the GPU path's Amdahl term (measurement pending, jobs 1862690/1862732).

### Device state (P6)
- In unit 2's handle (key `m%Sgmt`):
  - host build buffers: `std::vector`s, capacity reused;
  - device tree arrays, sized on the first call from the capacity hint `capPart = size(Part%P,2)`: nodes = capPart/4, leaves = capPart;
  - growth by doubling as the fallback.
- On host backends the kernels read the host buffers through unmanaged views (no copy). On device backends the flat tree is copied in.
- L0 P6: 0 allocations after warm-up (l0_tc and l0_iv with VelocityMethod 2).

## Audit list → L0 case
| # | Trap | Port | L0 |
|---|---|---|---|
| 8 (placement) | order-2 multipole loop form ("BAR_OLAF fails with unrolled") | kept literally | claim 1 bitwise, every case |
| 11 | `!$OMP SIMD` in `ui_leaves_nograd` reorders with OpenMP on | oracle = OPENMP=OFF | all |
| 12 | `T_Tree` points into `Part` | the device tree owns copies; Fortran tree not grown | — |
| 13 | default-kind literal `1.0e-5` (degenerate extent) | passed from the shim as `real(1.0e-5, ReKi)` | kinds 6/7: clumps of extent 0.9e-5 / 1.1e-5 |
| T1 (new) | `**0.5_ReKi` → `pow` (not `sqrt`) | strict uses `pow` | claim 1 |
| T2 (new) | `r**n` powi chains | reproduced | claim 1 (mutation pending) |
| T3 (new) | dead nodes via `EqualRealNos` (absolute 1.1e-14 on Σ\|α\|) | same test | kinds 8 (zero α), 9 (tiny α) |
| T4 (new) | **Fortran singularity**: a zero-radius node (all particles coincident, or n = 1) with RegFunctionPart = 0 and DistanceDirect = 0. A target at its center fails `r < max(0,0)`, so it gets the multipole at r ≈ 1e-15 (≈ 1e31). | reproduced (strict bitwise); excluded from R2 and counted | kind 5 (identical clump), targets on particles |
| T5 (new) | RegParam is used by the build (maxRegParam) even for RegNone | the source set copies the real RegParam for method 2 | reg 0 cases |
| T6 (new) | n = 0: `Uind = Uind + 0` still runs (turns -0 into +0) | same | n = 0 cases with -0 in U |
| T7 (new) | root n = 1: radius 0, monopole only, maxReg = RegParam(1) | same | n = 1 |
| T8 (new) | NCRIT = 32 boundary | — | n = 31, 32, 33 |
| T9 (new) | grid points on cell boundaries (`>` vs `>=`) | `>` | kind 3 |

## Gate results (partial, checkpoint)
| ID | Result | Evidence |
|---|---|---|
| L0 claim 1 (strict vs Fortran OPENMP=OFF) | **PASS** (CPU) | l0_tc 951/951 (seed 1), 959/959 (seed 7), 0 bit-diff cases; l0_iv TreePart bitwise through the shims |
| L0 claim 2 (strict CUDA) | **PASS** (A100) | job 1862865: 1.8e-16 of Σ\|c\| (bound 1e-13); R2, R7, P6 and unit 2's l0_iv also pass on the GPU |
| R2 (fast vs direct ≤ 1.1 × Fortran tree vs direct) | **PASS** (CPU) | worst ratio 0.909; all regs × BF {1, 1.5, 2, 3} × DD {0, 0.3L, 2.2L}; edge cases listed above; 5 singular cases reproduced |
| claim 3 (info) | 5.4e-14 of Σ\|c\| | no acceptance flips seen |
| R7 | **PASS** (CPU) | r7_threads.sh: md5-identical, 1/4/8 threads, 2 processes |
| P6 (L0) | **PASS** (CPU) | 0 allocations in steady steps |
| A7 (L0) | **PASS** (CPU) | l0_iv `a7_sequence(2)`: re-evaluation, rollback, misc copies, orphaned handle |
| R6, R8, P2, P6 (full run), A1/A2/A5 | NOT RUN | next steps 2 |
| P1/P5 | partial | CPU sweep (job 1862863): fast 1.18–1.64× the Fortran+OpenMP at every N from 1e3 to 1e6; GPU sweep pending (jobs 1862690/1862732) |

## Segment tree: what it needs (not started)
- `grow_tree_segment` (`:947`) and `ui_tree_segment` (`:1777`) have the same tree skeleton with segments as sources. The per-node center and moments are built from segment midpoints and their Γ·(P2−P1).
- The leaves call the segment kernel. That is unit 1's `seg_term`, with `RegFunction` 0–4 and the `SegConnct` indirection: 1-based, rows 3/4 unused by the tree.
- The device source set would be unit 2's segment SoA. The flat preorder/skip layout, the traversal and `multipole` are reusable as is.
- The build needs its own substep translation (`:1066`), including its own literal `1.0e-5` (`:1179`). Its leaf kernel order must follow `ui_tree_segment_11` (`:1807`).
- It is VelocityMethod 4: not the default, and used for LL only on request.

## Lessons so far (for the final report)
1. **R2 is easiest to satisfy by not changing the approximation.** Keeping the Fortran's topology rules (which are cheap to reproduce exactly: min/max/counts and a float recursion for cell centers) gives E(fast) = E(Fortran) and turns R2 into a rounding-level check. A different, GPU-friendlier tree would trade a guaranteed R2 for speed.
2. **"Strict GPU-legal" does not need a stack:** a preorder layout with skip pointers gives the recursion's exact visiting order as a flat loop.
3. **The Fortran tree has a singular corner** (T4) that makes R2's ratio meaningless there. The R2 test must classify such cases instead of dividing by them.
4. **Test-side traps cost time:** `NWTC_Library%Pi` is 0 until `NWTC_Init`, and editing a bash script while it runs breaks the rest of the script.
5. Fortran's "parallel tree build" is still serial at the top two levels; so is any Fortran-faithful host build. That serial part is the GPU path's Amdahl term, which motivates the device build.
