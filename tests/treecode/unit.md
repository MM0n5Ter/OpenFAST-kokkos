# treecode

The OLAF particle treecode (`VelocityMethod(1) = 2`, OLAF's default wake path): `grow_tree_part` and `ui_tree_part` (`FVW_VortexTools.f90`) ported to C++/Kokkos behind the frozen Fortran interfaces of `InducedVelocitiesAll_Init/Calc` and `LiftingLineInducedVelocities`. It lives in unit 2's per-instance handle. Unit 3 of 3 (architecture.md §7).


- **Owner**: agent session `199ad883` (subagents, unit 3). Dates: 2026-11-14 (port, checkpoint), 2026-11-15 (overlay, system level, device build, mutations).
- **Toolchain**: Kokkos 5.3.0 (C++20), gcc/gfortran 12.3.0, cmake 3.31.12, CUDA 12.3, OpenBLAS 0.3.28.

## Status (final, 2026-11-15)

**Done.** Every gate in scope passes on CPU and on the A100 (MIG); the GPU timing sweeps with the device build are queued (P1/P3/P5 on full GPUs). The segment tree (VelocityMethod 4) is deferred (below).

- **Code** (overlay layout, `src/`): the strict and fast tree from the 2026-11-14 checkpoint, plus the **fast-mode device build** (`src/modules/kokkos/src/treecode/treecode_dbuild.hpp`, level-synchronous K1-K5, no sort library, no atomics). Selection (`OFK_TC_DEVBUILD`): `auto` (default) = device backends, fast mode, n >= `OFK_TC_DEVBUILD_MIN` (4096); `1` = always, also on host backends (to test the algorithm on the CPU); `0` = never (host build, uploaded). Strict always builds on the host. `OFK_TC_DEVBUILD_CHECK=1` runs the host fast build as well and requires the same topology bit for bit.
- **Host backends are unchanged** by the device build (same code path in `tree_build`/`flat_tree` unless `OFK_TC_DEVBUILD=1`). So the CPU system results measured with `ov-fomp` (built before the device build) hold for the final source; `ov3-fomp` (final source) gives bit-identical bar-r8-tree outputs in fast and strict (job 1887985).
- **Trees** (`opt/build/`): test trees `u3d-omp`, `u3d-iv-omp`, `u3d-cuda-a100`, `u3d-iv-cuda-a100`, `u3d-cuda-h100`; OpenFAST `ov3-cuda-a100`, `ov3-fomp`, `ov3-sp-omp` (single precision, aerodyn_driver only). Frozen timing binaries: `u3d-frozen-cuda-{a100,h100}/perf_tc`. The old `u3-*` and `fork3-*` trees were built from `fork/` and are superseded (u3-cuda-a100/h100 are still held by queued jobs 1875922/1862732).
- **Open** (none blocks the gates):
  1. GPU sweeps 1887932 (A100) and 1887933 (H100) are pending: they give the full-GPU device-build timing and the evidence for `OFK_TC_DEVBUILD_MIN` (the sweep has a host-build/device-build A/B at every N). On the MIG slice the device build already wins at 1e4 (6.7 ms vs 16.6 ms per build, noisy CPU) and at 1e6 (22 ms vs 86 ms).
  2. Job 1862732 (H100, old host-build binary `u3-cuda-h100`, from the checkpoint) is superseded by 1887933, which also measures the host build (`devbuild0` rows). The parent may cancel it.
  3. The device build costs about 1.4 ms per build at 21k particles on the MIG slice (BAR: 0.42 s over 303 builds): roughly 10 kernel launches and 2 scalar D2H per level. Launching K1/K3 with an upper-bound league size would remove one D2H per level, if the sweep shows that small-N latency matters.
  4. E1-E6: E1 below; E2/E5/E6 are the parent's ledger rows (this session wrote none, as instructed).

### Slurm jobs (logs in `opt/build/`)
| Job | What | Log / output | Result |
|---|---|---|---|
| 1862730, 1862863, 1862864, 1862865 (checkpoint) | CPU sweeps, CPU L0, MIG L0 of the fork-built code | `u3-cpu.*`, `u3-gpu.ofk-treecode-mig.1862865.log` | see the gate table |
| 1862690 `ofk-treecode-perf` (a100, done) | GPU sweep, **host build** uploaded (checkpoint binary) | `u3-gpu.ofk-treecode-perf.1862690.log`, `u3-perf-a100-1862690/` | 1e6: build 61.8 ms, eval 5.5 ms per call |
| 1862732 `ofk-treecode-perf` (h100,h200, **pending**) | same, sm_90, old binary | `u3-gpu.ofk-treecode-perf.1862732.log` | superseded by 1887933 |
| 1863342 `ofk-infra-sys` (med, done; the overlay check) | ctest `aerodyn_utest\|OLAF` 17/17; BAR text outputs | `ov-cpu.ofk-infra-sys.1863342.log`, `ov-sys-1863342/` | R6, claim 1 at system level |
| 1887061 `ofk-treecode-sys` (med, done) | R8/R7 on `ov-fomp`, P2 (2 reps, `OFK_TC=fortran` A/B), P6 full run | `u3-cpu.ofk-treecode-sys.1887061.log`, `u3-sys-1887061/` | all PASS |
| 1887678 `ofk-treecode-migdb` (MIG, done) | first CUDA run of the device build: l0_tc, l0_tc_devbuild, perf A/B | `u3d-gpu.ofk-treecode-migdb.1887678.log` | PASS |
| 1887850 `ofk-treecode-mig` (MIG, done) | GPU gates on `u3d-*`/`ov3-cuda-a100`: l0_tc 6/6, l0_iv 12/12 (x2), BAR on GPU: R7, R8, P6 | `u3-gpu.ofk-treecode-mig.1887850.log`, `u3-mig-a100-1887850/` | all PASS |
| 1887932 `ofk-treecode-perf` (a100, **pending**) | P1/P3/P5 sweep 1e3..1e6, fast/strict + host/device-build A/B, `U3=u3d-frozen` | `u3-gpu.ofk-treecode-perf.1887932.log`, `u3-perf-a100-1887932/sweep.csv` | `a100,csv,` and `a100,devbuild{0,1},csv,` rows; `report` lines |
| 1887933 `ofk-treecode-perf` (h100,h200, **pending**) | same on sm_90 | `u3-gpu.ofk-treecode-perf.1887933.log`, `u3-perf-h100-1887933/` | same |
| 1887985 `ofk-treecode-ov3check` (med, done) | `ov3-fomp` (final source) vs `ov-fomp` on bar-r8-tree; device build forced on the CPU inside OpenFAST; single-precision run | `u3d-cpu.ofk-treecode-ov3check.1887985.log`, `u3d-ov3check-1887985/` | fast and strict IDENTICAL to job 1887061's; forced device build R8 2.2e-13, 0/542; single precision 301 finite data lines |

### Files (this unit; uncommitted)
| File | What |
|---|---|
| `src/modules/kokkos/src/treecode/treecode_common.hpp` | flat-tree types, build fragments, `TreeConsts`, `EvalArgs`, `DevBuild` (device build buffers) |
| `src/modules/kokkos/src/treecode/treecode_impl.hpp` | host build (Fortran-faithful, parallel subtrees), multipole, stackless traversal; includes the device build in the fast TU |
| `src/modules/kokkos/src/treecode/treecode_dbuild.hpp` | **new 2026-11-15**: the fast-mode device build (`build_device`, `write_device`) |
| `src/modules/kokkos/src/treecode/treecode_{strict,fast}.cpp` | the two mode TUs |
| `src/modules/kokkos/src/induced-velocities/induced_velocities.{cpp,hpp}` | unit 2's handle: tree state, `tree_build` (host or device build), `flat_tree`, C API `ofk_iv_tree_build`, `ofk_iv_eval_tree`, `ofk_tc_enabled`; `OFK_TC`, `OFK_TC_DEVBUILD[_MIN,_CHECK]` |
| `src/modules/kokkos/{CMakeLists.txt,ofk_kernels.cmake}` | `treecode_strict.cpp` strict, `treecode_fast.cpp` fast |
| `src/modules/aerodyn/src/FVW_Subs.f90` | shims (+80 lines under `#ifdef OPENFAST_KOKKOS`), see Contract |
| `tests/treecode/*` | `CMakeLists.txt` (+ `l0_tc_devbuild[_seed7]`, `l0_tc_r7_devbuild`), `tc_api.f90`, `l0_tc.f90` (kinds 6/7 tightened, see audit 13), `perf_tc.f90`, `r7_threads.sh`, `make_cases.sh`, `build_u3.sh` (`U3=` prefix), `cpu.sbatch` (overlay paths), `gpu.sbatch` (`U3=`/`OV=` trees, device-build runs, P6 bytes > 0), `mutate.py` (new), this card |
| `tests/induced-velocities/l0_iv.f90` | method-2 scale and A7 sequence for wake method 2 (checkpoint) |

E1 (size): C++ of the treecode 1244 lines (`treecode_*`; the device build 557 of them), handle changes in `induced_velocities.cpp` about 250 lines, Fortran shims 80 lines; Fortran in scope (`grow_tree_part*`, `compact_tree_part`, `ui_tree_part*`, `ui_leaves_nograd`) about 550 lines of `FVW_VortexTools.f90`.

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
- The **same topology rules**. On host backends the fast mode uses the host build compiled with FMA and `sqrt` for the norm; on device backends the device build below. Both give the Fortran's tree up to the `EqualRealNos` edge.
- Fast-mode evaluation on host backends: each target's leaf sums go to 8 fixed lanes across all of its buckets, with `omp simd` and a fixed pairwise tree at the end. On device backends: 1 lane.
- **Why this design.** R2 asks the fast tree's error vs direct to be ≤ 1.1 × the Fortran tree's. With the same tree, E(fast) = E(Fortran) up to rounding: L0 worst ratio 0.909 = 1/1.1. So R2 holds by construction rather than by luck on each input. A different tree (Morton/LBVH, another leaf size or MAC) would have needed a tighter MAC to be safe on every input, at a cost in speed.
- The device build keeps this property: octant decisions use the Fortran's float recursion for the geometric centers, so the topology is the Fortran's (checked bit for bit by `OFK_TC_DEVBUILD_CHECK`).
- The GPU speed comes from the evaluation (thread per target, no stack) and from the device build. With the host build, the build was 92% of the A100's time per call at 1e6 (job 1862690: 61.8 of 67.3 ms).

### Fast, device backends: the device build (`treecode_dbuild.hpp`)
- **Level-synchronous**, one BFS record per node; a level's records are contiguous, a record's children contiguous and in octant order. Per level:
  - A: chunk layout (records cut into chunks of 512 positions; exclusive scan, one scalar D2H);
  - K1, a team per chunk: Σw·p, Σw, maxReg, min/max xyz, octant per particle (`>` against the record's geometric center), octant counts;
  - K2, a team per record: combine the chunks; the Fortran's decisions (dead, leaf rule, degenerate extent), barycenter, each chunk's octant base (integer team scans);
  - K3/K4: moments about the barycenter, per chunk then per record;
  - K5: children (scan of child counts, one scalar D2H), a stable partition by octant (byte-packed `uint64` team scan per 128-position tile, integer and exact), child records with `G + (0.5·(R/2))·(±1)`, `R/2`.
  - Then bottom-up (subtree sizes and leaves) and top-down (preorder index, leaf start) passes, and `write_device`: the same flat arrays as the host build, leaves gathered in the host build's order.
- **Determinism (R7):** float sums are per-thread strided loops plus a fixed pairwise tree in team scratch over a fixed team size (128 on CUDA, 1 on host backends); every partition is integer arithmetic. No `sort_by_key` (thrust temporaries are per-call allocations, P6), no atomics.
- **Exactness:** min/max, counts and the geometric centers (`±q5` with q5 = 0.5·R/2: FMA cannot change `q5·(±1)`) equal the host build's, so the topology is identical; only barycenters and moments round differently. n <= 1 stays on the host.
- **Memory (P6):** `DevBuild` in the handle, sized on the first call from `capPart` (records capPart/4, chunks capPart/512 + records); records grow by doubling (contents kept) as the fallback.

### Device state (P6)
- In unit 2's handle (key `m%Sgmt`):
  - host build buffers: `std::vector`s, capacity reused;
  - device tree arrays, sized on the first call from the capacity hint `capPart = size(Part%P,2)`: nodes = capPart/4, leaves = capPart;
  - growth by doubling as the fallback.
- On host backends the kernels read the host buffers through unmanaged views (no copy). On device backends the flat tree is written by the device build (fast) or copied in (strict, or n below `OFK_TC_DEVBUILD_MIN`); the device build uploads the particles (7 columns) once per build instead.
- P6: 0 allocations after warm-up in L0 (l0_tc and l0_iv with VelocityMethod 2, both builds) and in full BAR runs (CPU and GPU).

## Audit list → L0 case
| # | Trap | Port | L0 |
|---|---|---|---|
| 8 (placement) | order-2 multipole loop form ("BAR_OLAF fails with unrolled") | kept literally | claim 1 bitwise, every case |
| 11 | `!$OMP SIMD` in `ui_leaves_nograd` reorders with OpenMP on | oracle = OPENMP=OFF | all |
| 12 | `T_Tree` points into `Part` | the device tree owns copies; Fortran tree not grown | — |
| 13 | default-kind literal `1.0e-5` (degenerate extent) | passed from the shim as `real(1.0e-5, ReKi)` | kinds 6/7: clumps of extent 1e-5·(1 ∓ 1e-9). The earlier 0.9e-5/1.1e-5 clumps let a hard-coded `1.0e-5f` survive (mutation check); now caught |
| T1 (new) | `**0.5_ReKi` → `pow` (not `sqrt`) | strict uses `pow` | claim 1 |
| T2 (new) | `r**n` powi chains | reproduced | claim 1 (mutation: caught, 134 cases) |
| T3 (new) | dead nodes via `EqualRealNos` (absolute 1.1e-14 on Σ\|α\|) | same test | kinds 8 (zero α), 9 (tiny α) |
| T4 (new) | **Fortran singularity**: a zero-radius node (all particles coincident, or n = 1) with RegFunctionPart = 0 and DistanceDirect = 0. A target at its center fails `r < max(0,0)`, so it gets the multipole at r ≈ 1e-15 (≈ 1e31). | reproduced (strict bitwise); excluded from R2 and counted | kind 5 (identical clump), targets on particles |
| T5 (new) | RegParam is used by the build (maxRegParam) even for RegNone | the source set copies the real RegParam for method 2 | reg 0 cases |
| T6 (new) | n = 0: `Uind = Uind + 0` still runs (turns -0 into +0) | same | n = 0 cases with -0 in U |
| T7 (new) | root n = 1: radius 0, monopole only, maxReg = RegParam(1) | same | n = 1 |
| T8 (new) | NCRIT = 32 boundary | — | n = 31, 32, 33 (mutation 32→31 caught) |
| T9 (new) | grid points on cell boundaries (`>` vs `>=`) | `>` | kind 3 |

## Gate results (final, 2026-11-15)
CPU = Xeon 8480+ (login node or `med`); GPU = A100 MIG 2g.20gb on `interactive_gpu` (EPYC 7443 host). L0 logs: `opt/build/u3d-omp/ctest.login2.log`, `u3d-iv-omp/ctest.login{,.db}.log`, job logs as named.

| ID | Result | Evidence |
|---|---|---|
| Claim 1, L0 (strict vs Fortran OPENMP=OFF) | **PASS** (CPU) | `u3d-omp` l0_tc 951/951 (seed 1), 959/959 (seed 7, 80 random): 0 cases with bit differences. l0_iv TreePart bitwise through the shims (`u3d-iv-omp` 15/15) |
| Claim 1, system | **PASS** | bar-r8-tree text output from line 8: ov-omp strict (Kokkos tree) = ov-off, 301 lines (job 1863342). Also ov-fomp strict = base₁ (pin OPENMP=ON, 1 thread) (job 1887061) |
| Claim 2 (strict CUDA vs Fortran) | **PASS** (A100) | job 1887850: 1.21e-16 (seed 1), 1.80e-16 (seed 7) of Σ\|c\| (bound 1e-13); l0_iv TreePart 4.35e-17 |
| R2 (fast vs direct ≤ 1.1 × Fortran tree vs direct) | **PASS** CPU and GPU, host and device build | worst ratio 0.9091 in every run: `u3d-omp` l0_tc, l0_tc_devbuild (seed 1, 7); job 1887850 the same four. 5 (7 with the 20000-particle clouds) Fortran-singular cases reproduced and excluded (T4) |
| Claim 3 (fast vs strict, info) | 3.4e-14 / 3.7e-14 CPU (host / device build); 3.8e-14 / 4.6e-14 GPU | 0 acceptance flips (> 1e-12) anywhere. Device-build topology = host build's, bit for bit, in every L0 case (`OFK_TC_DEVBUILD_CHECK=1`, CPU and GPU) |
| R6 | **PASS** | ctest `aerodyn_utest\|OLAF` 17/17 on ov-off, ov-omp strict, fast, fast `OFK_TC=fortran` (job 1863342 + the py_ad_B1n2_OLAF fix, notes/overlay.md §4.3); GPU 11/11 (job 1863556) |
| R7 | **PASS** | L0: md5-identical dumps across processes and 1/4/8 threads, host build and device build, CPU (`u3d-omp`) and GPU (job 1887850). System: bar-r8-tree fast 1 vs 8 threads identical (job 1887061); GPU fast two processes identical, auto and device build (job 1887850) |
| R8 (bar-r8-tree, 300 steps × 542 channels; base₁ = pin OPENMP=ON 1 thread, envelope = pin with FMA, max 6.39e-13) | **PASS** | CPU (job 1887061, ov-fomp): strict 0; fast 1 and 8 threads 2.08e-13 (worst ratio 0.37); 0/542 over 2×. GPU (job 1887850): strict 5.3e-14, fast (auto) 4.29e-13, fast device build 6.7e-14; 0/542. No channel ever > 1e-6 |
| OFF-vs-OMP drift (info) | **0** | ov-off (OPENMP=OFF) = base₁ (OPENMP=ON) on every channel; the leaf `!$OMP SIMD` reordering does not show at system level in this case |
| P1/P5 (kernel, grow + eval per call) | CPU **PASS**; GPU partial | CPU, 8 threads (job 1862863; host code unchanged since): fast 1.18–1.64× Fortran+OpenMP at N = 1e3..1e6 (1e6: 0.365 vs 0.455 s). A100 host build (job 1862690): 1e6 67.3 ms (eval 5.5 ms vs Fortran 260 ms on 8 Xeon cores, 47×; total 6.8×). Device build on the MIG slice (job 1887678): 1e6 build 22.2 ms vs 86.0 ms host build, eval 11.2 ms. Full A100/H100 with the device build: jobs 1887932/1887933 (pending) |
| P2 (CPU users lose nothing, 8 threads, 2 reps) | **PASS** | job 1887061, ov-fomp vs openfast-pin-fomp: bar-r8-tree 41.0 vs 44.1 s (**1.08×**; `OFK_TC=fortran` 1.01×, strict 1.04×); bar-perf-tree 47.9 vs 53.9 s (**1.13×**; `OFK_TC=fortran` 1.03×, strict 1.06×). The tree is 11–17% of these runs (bar-r8-tree: eval 4.0 s + build 0.30 s of 40.3 s; bar-perf-tree: 7.5 + 0.49 s of 48.0 s) |
| P3 (GPU split) | partial | BAR on MIG (job 1887850): tree build 0.42 s / 303 builds, particle H2D 0.044 s, of 34.1 s wall; perf sweep `report` lines pending (1887932/3) |
| P4 | not rerun | the r-test cases are below `OFK_TC_DEVBUILD_MIN` on most steps; unit 2's P4 stands |
| P6 | **PASS** | L0: 0 allocations in 20 steady grow+eval steps (CPU and GPU, host and device build). Full run, CPU (job 1887061): 4 Kokkos allocations, last after evaluation 2 of 903. Full run, GPU with the device build (job 1887850): **0 allocations with bytes > 0 after evaluation 3** (12372 zero-byte `TeamScratchMemory` events ignored, architecture.md §5) |
| A1 | **PASS** | `nm` of `libaerodynlib.a`: no `_MOD_` symbol of ov-pin missing in ov-off, ov-omp, ov-fomp, ov-cuda-a100; the one added is unit 2's `fvw_subs::fvw_kokkos_release` |
| A2 | **PASS** | the three overlaid Fortran files preprocess (`-DOPENFAST_DOUBLE_PRECISION`, no `OPENFAST_KOKKOS`) to exactly the pin's; every other change is in `src/modules/kokkos/` and `tests/` |
| A5 | **PASS** | OFF build = pin: `tests/overlay/a5.sh ov-pin ov-off` (`ov-a5.log`): 76/76 `.mod` identical, 17/19 objects identical, the rest differ only in line-number immediates (the Fortran sources are unchanged since). `DOUBLE_PRECISION=OFF` with Kokkos-OpenMP: `ov3-sp-omp` (`tools/openfast-build.sh omp … -DDOUBLE_PRECISION=OFF`, aerodyn_driver) compiles and links (log `ov3-sp-omp.build.log`); it runs bar-r8-tree to the end with finite outputs (job 1887985) |
| A6 | **PASS** (code) | no non-const `static` in `treecode_*`; device build buffers live in the per-instance handle |
| A7 | **PASS** CPU and GPU | l0_iv `a7_sequence` with wake method 2 (re-evaluation, rollback, misc copies, orphaned handle): 544 comparisons bitwise in strict, CPU (`u3d-iv-omp`) and GPU (job 1887850), with the default and the forced device build |
| Mutation checks (L0 sensitivity) | **12/12 caught** | `tests/treecode/mutate.py` on a copy of `src/`, log `opt/build/u3d-mutations2.log`: octant `>=`, powi r⁵ = r⁴·r, pow→sqrt, singleton order, `1.0e-5f`, dead-node tolerance 0, NCRIT 31 (claim 1 fails); fast-only MAC halved, device-build radius ×0.25 (R2 fails); device-build octant `>=`, singleton order, child-center sign (topology check fails). First run 11/12: `1.0e-5f` survived until kinds 6/7 were tightened |
| E1–E6 | E1 above | E2/E5/E6 in the parent's ledger |

## Segment tree: deferred (out of scope for unit 3; not started)
- `grow_tree_segment` (`:947`) and `ui_tree_segment` (`:1777`) have the same tree skeleton with segments as sources. The per-node center and moments are built from segment midpoints and their Γ·(P2−P1).
- The leaves call the segment kernel. That is unit 1's `seg_term`, with `RegFunction` 0–4 and the `SegConnct` indirection: 1-based, rows 3/4 unused by the tree.
- The device source set would be unit 2's segment SoA. The flat preorder/skip layout, the traversal and `multipole` are reusable as is.
- The build needs its own substep translation (`:1066`), including its own literal `1.0e-5` (`:1179`). Its leaf kernel order must follow `ui_tree_segment_11` (`:1807`).
- It is VelocityMethod 4: not the default, and used for LL only on request.

## Lessons
1. **R2 is easiest to satisfy by not changing the approximation.** Keeping the Fortran's topology rules (which are cheap to reproduce exactly: min/max/counts and a float recursion for cell centers) gives E(fast) = E(Fortran) and turns R2 into a rounding-level check. A different, GPU-friendlier tree would trade a guaranteed R2 for speed.
2. **"Strict GPU-legal" does not need a stack:** a preorder layout with skip pointers gives the recursion's exact visiting order as a flat loop.
3. **The Fortran tree has a singular corner** (T4) that makes R2's ratio meaningless there. The R2 test must classify such cases instead of dividing by them.
4. **Test-side traps cost time:** `NWTC_Library%Pi` is 0 until `NWTC_Init`, and editing a bash script while it runs breaks the rest of the script.
5. Fortran's "parallel tree build" is still serial at the top two levels; so is any Fortran-faithful host build. That serial part is the GPU path's Amdahl term, which motivates the device build.
6. **A device rewrite of a host algorithm should be testable on the CPU and checkable against the host.** The device build is written over the execution space (team size 1 on host backends), so the whole CPU L0 exercised it in seconds on the login node; the first CUDA run then passed. `OFK_TC_DEVBUILD_CHECK` compares the exact parts (topology, radii, leaf copies) with the host build. It is the only check that catches device-only structural mutations: R2 passes for any valid tree.
7. **Mutation checks find test gaps that "covers every branch" misses.** The degenerate-extent threshold had straddling cases (0.9e-5 / 1.1e-5), and a `1.0e-5f` literal still survived. A straddle must be tighter than the smallest plausible wrong constant: for a kind trap that is the float-vs-double gap (relative 2.5e-9).
8. **Amdahl at system level.** The kernel is 1.2–1.6× faster on CPU and the tree build 4× faster on the GPU, but the BAR runs gain 1.08–1.13× on CPU, because the tree is 11–17% of the step. The other host Fortran dominates there (phase-2 residency).
9. **Freeze the binaries a queued timing job will measure** (`u3d-frozen-*`). Otherwise a rebuild between submission and start silently changes what gets timed.

