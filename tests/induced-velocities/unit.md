# induced-velocities

The direct-method paths of OLAF's induced-velocity driver, ported to C++/Kokkos behind the frozen Fortran interfaces, with a persistent per-instance device handle. Unit 2 of 3 (architecture.md §7). The treecode paths (`VelocityMethod` 2 and 4) stay Fortran (unit 3).

- **Status (v1, 2026-11-14)**: in the fork, built in all three configurations (OFF, Kokkos-OpenMP, Kokkos-CUDA). See the gate table.
- **Owner**: agent session `199ad883` (subagent, unit 2). Date: 2026-11-14.
- **Toolchain**: Kokkos 5.3.0 (C++20), gcc/gfortran 12.3.0, cmake 3.31.12, CUDA 12.3, OpenBLAS 0.3.28.
- Builds on unit 1 (`tests/biotsavart-direct/unit.md`): the kernels (`seg_term`/`part_term`, square and skinny shapes, strict/fast TUs) are reused unchanged; unit 1's loop shapes gained entry points that take a pre-packed SoA block and a persistent scratch buffer.

## Status / resume (2026-11-14 checkpoint)

### Results that arrived after the checkpoint (read by the parent session, 2026-11-14 20:30)
- **R8 (job 1860993): PASS.** The envelope is the pin built with FMA (`-march=x86-64-v3`) at 1 thread.
  - Tree-wake case: envelope max 6.39e-13; fast drift 1.86e-13; 0/542 channels over 2× envelope; no channel ever above 1e-6.
  - Direct-wake case: envelope max 2.60e-13; fast drift 9.70e-14; 0/542 over 2× envelope.
  - Strict: 0 drift in both cases.
  - The Fortran 1-vs-8-thread envelope is 0 in both cases (job 1860700), so the FMA envelope is the meaningful one.
- **P2, tree wake (job 1860992), with Fortran OpenMP kept** (fork built OPENMP=ON + Kokkos-OpenMP): **PASS**.
  - Fortran 8 threads: 44.92 s; fork fast: 44.07 s (1.02×); fork strict: 44.11 s.
  - The earlier 2× slowdown came from the OPENMP=OFF build configuration.
- **P2, direct wake, same configuration:** fork fast 145.7 s vs Fortran+OpenMP 299.7 s (job 1860602), i.e. **2.06×**.
- **P6, full run on CPU (job 1860992): PASS.**
  - Kokkos memory tool: 4 allocations (packed sources) or 5 (lattice sources), all made by evaluation 3; none in the remaining ~900 evaluations.
  - Logs: `opt/build/iv-perf2-1860992/memtool.*`.
- **Particle-wake BAR case (job 1860602):** Fortran+OMP8 2.04 min, fast 2.00 min, strict 2.12 min.
- **GPU L0 on an A100 MIG slice (job 1862373, sm_80 build, `opt/build/a100mig.1862373.log`):**
  - 12/12 ctest pass.
  - **Claim 2:** strict GPU vs Fortran worst 1.63e-16 (Basic) / 4.35e-17 (Part); pack and A7 bitwise.
  - **Claim 3:** fast vs strict 8.25e-14.
  - **P6 (L0):** 0 device allocations in 20 steady steps.
  - **A7:** bitwise.
- **GPU full cases (job 1862603, A100 MIG, `fork-cuda-a100` with the portable OpenBLAS; log `opt/build/a100cases.1862603.log`):**
  - The first attempt (1862373) died with "Illegal instruction" in `dgemm_`, because the cluster openblas is a zen4/AVX-512 build and the A100 nodes are Zen3. Fixed with `opt/openblas-0.3.28-portable` on `LD_LIBRARY_PATH`.
  - **R8 on GPU: PASS.** Against the FMA envelope:
    - direct wake: strict 0, fast 9.45e-14 (0/542 channels over 2×);
    - tree wake: strict 0, fast 1.38e-13 (0/542).
    - No channel ever exceeds 1e-6.
  - **R7 on GPU: PASS.** Two processes give identical data lines in both cases; only the header timestamp differs. The job script's whole-file `cmp` said "DIFFER", so the script now compares data lines.
  - **P6 on GPU: PASS in substance.** Our own allocations: 15 events, the last one at evaluation 3.
    - The memory tool also logged 555 `Kokkos::CudaSpace::TeamScratchMemory` "allocations" of **0 bytes** up to the last evaluation.
    - That is Kokkos 5.3.0 behaviour: `CudaInternal::resize_team_scratch_space` "re-allocates" 0 bytes on every team launch that asks for no level-1 scratch (`core/src/Cuda/Kokkos_Cuda_Instance.cpp:434-438`). No device memory is allocated.
    - `memtool` should ignore 0-byte events.
  - **Timing:** not meaningful on a 2g.20gb MIG slice (direct-wake BAR: fast 41 s, strict 50 s). The full-GPU timings come from the h100 jobs.
- **H100 NVL (job 1860779, fork-cuda build, log `opt/build/iv-gpu.1860779.log`; parent session 2026-11-14 23:45):**
  - L0: 12/12 pass, with the same numbers as on the A100 (claim 2: 1.63e-16; claim 3: 8.25e-14).
  - **P1/P3, direct-wake BAR perf case:** whole run 46.4 s on GPU, vs Fortran+OpenMP (8 threads) 299.7 s, i.e. **6.5×**.
  - **Profile, BAR r8-direct, 37.8 s total:**
    - GPU kernel 4.34 s;
    - H↔D copies and packing under 0.1 s (400.6 MB H2D);
    - the remaining ~33 s is host Fortran (wake stepping, lifting line, rest of AeroDyn).
    - So on the GPU this unit is now **host-bound (Amdahl)**. The next gains are residency phase 2 and porting the remaining host work, not faster kernels.
  - Tree-wake case: 73.7 s, slower than Fortran+OpenMP (44.9 s). This build predates the ported treecode and has `OPENMP=OFF`, so the Fortran tree runs on one host thread.
  - The H100 nodes ran fine with the cluster openblas (job 1860007 too). They are not Zen3.
  - Unit 1's job 1860007 on the H100: the fork-cuda OLAF ctests pass 13/13 in both modes, and R7 holds across processes.
- **Still queued:** job 1860223 (full GPU plan; estimated start 2026-11-16 06:09, about 950 jobs ahead on h100/h200).

Session stopped by the user's disconnect at about 18:55. Code is complete and CPU-verified; GPU evidence and some CPU timing/R8 runs are pending in Slurm. Nothing is committed (rule: no commits).

**Done, with evidence**
- Code: handle + lattice mirror + device pack + shims (Files below), built in fork-off / fork-omp / fork-cuda, and in the test trees `opt/build/iv-omp`, `opt/build/iv-cuda`.
- L0 on CPU (`opt/build/iv-omp`, job 1860600 and login node): 15/15 ctest pass. Claim 1 (Fortran OPENMP=OFF vs strict): bitwise in both source paths (packed, lattice), 4816 comparisons, every VelocityMethod 1-4 x pack flags x RegFunction, sentinels. Claim 3 (strict vs fast): worst 8.2e-14 of sum|c| (Basic), 3.1e-15 (Part), PASS at 1e-12. Fortran's own FMA spread on the same inputs (`l0_iv_ffma`): 8.2e-14 (Basic), 6.0e-14 (Part). R7: bitwise in process and across processes / 1, 4, 8 threads, both modes. A7 sequence: bitwise vs fresh instances and vs the oracle. P6 at L0: 0 allocations in 20 steady steps. Mutations M1, M2, M3, M5 caught.
- R6 (job 1860600, log `opt/build/iv-cpu.ofk-induced-velocities-r6.1860600.log`): 16/17 on fork-off, fork-omp fast and strict alike (`py_ad_B1n2_OLAF` fails everywhere, upstream path assumption). `.outb` strict vs OFF identical in all 17; fast identical in all aerodyn_driver outputs (output resolution only, see lesson 5).
- R8, tree-wake BAR case (job 1860700, full-precision text output, 300 steps): Fortran OMP 1 vs 8 identical (envelope 0); fork-off and kk-strict identical to Fortran; kk-fast drift max 1.9e-13 of channel max (AB1N003Vindy), never > 1e-6; kk-fast at 2 and 8 threads identical. By the letter of R8 (<= 2 x envelope = 0): FAIL; see lesson 6.
- P2 direct-wake BAR (`bar-perf-direct`, job 1860602): Fortran+OpenMP 8 threads 299.7 s; fork-omp fast 153.8 s (1.95x); strict 294.2 s (1.02x). Profile (fast, packed): kernel 110.1 s, pack 0.10 s, rest (host Fortran) ~46.5 s of 156.7 s.
- A1 (no `_MOD_` symbol removed; one added: `fvw_subs_MOD_fvw_kokkos_release`), A5 (fork-off: `fvw_subs.mod` md5 = pin, `.text` differs only in line-number immediates, +59/+187).

**Slurm jobs (logs under `opt/build/`)**

| Job | What | Log / output | Look for |
|---|---|---|---|
| 1860223 `ofk-induced-velocities-gpu1` (h100/h200, 2 h, est. start 2026-11-16 06:09) | `gpu_steps.sh`: L0 on CUDA, BAR cases on GPU (fast/strict, profile, memtool), fork-cuda r-tests fast/strict | `iv-gpu.1860223.log`, `iv-gpu-1860223/` | claim 2 (strict CUDA vs Fortran <= 1e-14), claim 3, R7, A7, P6 `last_alloc_after_eval` in `memtool.*.txt`, P3 split in `cases/*.prof/stdout.log`, `.outb` of `cuda-*` vs `iv-r6-1860600/off` |
| 1860779 `ofk-induced-velocities-gpucore` (30 min, backfill) | `gpu_core.sh`: the core subset of the above | `iv-gpu.1860779.log`, `iv-gpu-1860779/` | same |
| 1860602 `ofk-induced-velocities-perf` (running) | P1/P2/P3 CPU: `bar-perf-direct`, `bar-perf-part` | `iv-cpu.ofk-induced-velocities-perf.1860602.log`, `iv-perf-1860602/` | `Total Real Time` in each `stdout.log` (the first `time.txt` is empty: `bc` missing, fixed) |
| 1860700 `ofk-induced-velocities-r8` (running) | R8 on `bar-r8-tree` (done) and `bar-r8-direct` | `iv-cpu.ofk-induced-velocities-r8.1860700.log`, `iv-r8-1860700/` | the `=== R8 bar-r8-direct` block |
| 1860993 `ofk-induced-velocities-r8ffma` (after 1860700) | R8 with the pin built with FMA (`openfast-pin-ffma`) as envelope | `iv-cpu.ofk-induced-velocities-r8ffma.1860993.log` | drift of kk-fast vs the FMA envelope per channel |
| 1860991 `ofk-induced-velocities-build-fomp` | builds `fork-omp-fomp` (Kokkos-OpenMP + Fortran OpenMP), aerodyn_driver | `iv-build.ofk-induced-velocities-build-fomp.1860991.log` | `build fork-omp-fomp OK` |
| 1860992 `ofk-induced-velocities-perf2` (after 1860991) | P2 with Fortran OpenMP kept (tree and direct), P6 memtool on CPU (new capacity hints) | `iv-cpu.ofk-induced-velocities-perf2.1860992.log`, `iv-perf2-1860992/` | tree case kf-fast8 vs fomp8 time; `last_alloc_after_eval` (expect <= a few evals) |
| 1860994 `ofk-induced-velocities-build-omp` (after 1860602, 1860700, 1860993) | rebuilds `fork-omp` with the final code | `iv-build.ofk-induced-velocities-build-omp.1860994.log` | `build fork-omp OK` |
| 1860007 (unit 1) | unit 1's gpu2, still pending | — | — |

**Half-done / stale**
- `opt/build/fork-omp` is **stale** (built before `ofk_iv_reserve`/`capCP`): arithmetic identical, P6 behaviour differs (doubling only). Rebuilt by job 1860994. `opt/build/fork-off`, `fork-cuda`, `iv-omp`, `iv-cuda` are current. `fork-omp-fomp` is half-built (job 1860991 finishes it).
- `openfast` and the C-binding `.so` in fork-omp also predate the last change (job 1860994 rebuilds everything).
- The gate table below is filled from CPU evidence only; GPU rows are PENDING.

**Next steps for whoever resumes, in order**
1. When 1860700 / 1860993 finish: read the R8 blocks; record per-channel max drift and first step > 1e-6 for `bar-r8-direct` (and confirm whether Fortran OMP 1 vs 8 is again identical there; lesson 6 currently claims it only for the tree case).
2. When 1860992 finishes: P2 for the tree case with Fortran OpenMP kept; P6 on CPU from `memtool.bar-r8-direct.{packed,lattice}.txt` (expect the last allocation within the first few evaluations).
3. When a GPU job runs: fill claims 2/3, R7, A7, P6, P3 (kernel / H2D+D2H / host from the `[ofk_iv]` line and `Total Real Time`), fork-cuda R6 and `.outb` vs OFF; compare `cases/bar-r8-*.{fast,strict}/ad_driver.out` with `iv-r8-1860700/*.fomp1/ad_driver.out` using `r8.py` (GPU R8).
4. Then re-run R6 on the rebuilt fork-omp (optional; arithmetic unchanged) and write E1-E6 (compute hours: `tools/cost/compute_hours.sh 2026-11-14 ofk-induced-velocities`).
5. Decide (human): R8 envelope definition when OMP_1 = OMP_N (lesson 6), and whether `fork-omp` should be built with `OPENMP=ON` (lesson 7).

## Files

| Where | What |
|---|---|
| `fork/modules/aerodyn/src/FVW_Subs.f90` | Edited in place (+187 lines, all under `#ifdef OPENFAST_KOKKOS`): module-level `bind(C)` interfaces (private), three private helpers (`ofk_iv_set_sources`, `ofk_iv_sources_from_packed`, `ofk_iv_direct`), one public `FVW_Kokkos_Release(m)`; shims in `InducedVelocitiesAll_Init` (source set after the Fortran pack), `InducedVelocitiesAll_Calc` (`idVelocityBasic`/`idVelocityPart` cases), `LiftingLineInducedVelocities` (same two cases). Every original line is kept in `#else` or unchanged. `WakeInducedVelocities`, `InducedVelocitiesAll_OnGrid_Calc` and `srcPnl_ExtVelocities_OnPanels` reach the port through Init/Calc without edits. |
| `fork/modules/aerodyn/src/FVW.f90` | +4 lines in `FVW_End`: `FVW_Kokkos_Release(m)` (frees the handle). |
| `src/modules/kokkos/src/induced-velocities/induced_velocities.{hpp,cpp}` | C API; per-instance handles in a guarded process-wide map; lattice mirror with slab-wise validation; source sets (lattice or packed); evaluation; counters and profile; finalize hook. |
| `src/modules/kokkos/src/induced-velocities/induced_velocities_common.hpp`, `induced_velocities_pack.cpp` | Device packing (F1): lattice → SoA segments (closed-form offsets, mirror) → particles. Compiled without contraction in every mode (class B). |
| `src/modules/kokkos/src/biotsavart-direct/biotsavart_direct_{impl,common}.hpp`, `biotsavart_direct.cpp`, new `biotsavart_direct_consts.hpp` | Unit 1, adapted: `run_seg_soa`/`run_part_soa` (pre-packed SoA + persistent scratch), `pack_*_into`, the skinny shape's partial sums in a caller-provided buffer, `make_consts` moved to a header, `ofk_bsd::tuning()` exported. Unit 1's entry points and results are unchanged (unit 1 L0 re-run: 13/13 as before). |
| `fork/modules/kokkos/CMakeLists.txt` | Two new sources; the pack TU is in `OFK_STRICT_SOURCES` (no FMA, `--fmad=false`). |
| `tests/induced-velocities/` | `CMakeLists.txt` (oracle from the pin with three modules renamed, port from the fork, L0 tests, memtool), `l0_iv.f90`, `l0_iv_util.f90`, `memtool.cpp` (Kokkos Tools memory events), `make_cases.sh` (scratch BAR cases), `run_case.sh`, `r8.py`, `cpu.sbatch`, `gpu.sbatch` + `gpu_steps.sh`. |

Unit 1's test directory is **not** modified.

## Contract (frozen, pin `12fa3fc`)

```fortran
subroutine InducedVelocitiesAll_Init(p, x, m, Sgmt, Part, Tree, Panl, ErrStat, ErrMsg, allocPart)   ! FVW_Subs.f90:1237
   type(FVW_ParameterType), intent(in), target :: p;  type(FVW_ContinuousStateType), intent(in) :: x
   type(FVW_MiscVarType),   intent(in), target :: m   ! NOTE intent(in): no handle can be stored in m here
   type(T_Sgmt), intent(inout) :: Sgmt                ! API OUTPUT: the packed segments (read by Calc, OnGrid_Calc, trees)
   type(T_Part), intent(inout) :: Part                ! API OUTPUT for VelocityMethod 2/3
   type(T_Tree), intent(out) :: Tree; type(T_Panl), intent(out) :: Panl; logical, intent(in) :: allocPart
subroutine InducedVelocitiesAll_Calc(CPs, nCPs, Uind, p, Sgmt, Part, Tree, Panl, ErrStat, ErrMsg)    ! :1285
   real(ReKi), dimension(:,:), intent(in) :: CPs      ! 3 x nCPs++ (sections allowed)
   real(ReKi), dimension(:,:), intent(inout) :: Uind  ! 3 x nCPs++, ACCUMULATED ("Side effects!!!")
   type(T_Sgmt), intent(in) :: Sgmt; type(T_Part), intent(in) :: Part   ! NOTE: no m, no x
subroutine InducedVelocitiesAll_End(p, Tree, Part, Panl, ErrStat, ErrMsg, deallocPart)               ! :1331  NOTE: no Sgmt
subroutine LiftingLineInducedVelocities(p, x, InductionAtCP, iDepthStart, m, ErrStat, ErrMsg)        ! :1488
   ! writes m%Sgmt (pack from iDepthStart), m%Part (VelocityMethod(2)=3), m%W(:)%Vind_CP, Vind_LL
subroutine WakeInducedVelocities(p, x, m, ErrStat, ErrMsg)                                           ! :1379
   ! writes m%CPs, m%Uind (tail m%Uind(:,nCPs+1:)=1000), m%Sgmt, m%Part, m%W(:)%Vind_NW/FW (-9999 / 2222 sentinels)
```

Callers (pin): `FVW.f90:855` (CCSD → Wake), `:1457` (CalcOutputForAD → LL, iDepthStart=1), `:1628-1640` (grid outputs: Init, OnGrid_Calc per grid, End), `:1749` (UA → LL); `FVW_Wings.f90:391` (circulation solve → LL, iDepthStart=iNWStart+1); `FVW_Subs.f90:2098-2102` (source panels); `FVW_Tests.f90:822,831`. All compile and run unchanged.

## Design

### The handle (goal 1, P6): option (b), a module-level map keyed per instance

Every entry point of the unit receives the instance's `T_Sgmt` (`m%Sgmt` in every in-tree caller): Init and Calc as `Sgmt`, LL and Wake through `m`. Its address is the key. The C++ side holds `std::unordered_map<const void*, Entry>` (one `Handle<double>`/`Handle<float>` per key), guarded by a mutex, never destroyed (no static-destruction race with `Kokkos::finalize`), emptied by a `Kokkos::push_finalize_hook`. The key is passed as a Fortran assumed-type scalar (`type(*)`), i.e. the address, with no `target` attribute needed and no interface change.

A handle owns: the lattice mirror (device) + its host shadow (pinned), the block table, the SoA source blocks (segments, particles) and their pinned staging (packed mode), the target/output buffers and pinned staging, the skinny scratch. All grow by doubling (`ensure1`/`ensure2`), never per call. Nothing in a handle is trusted without validation: the lattice is compared slab by slab, bitwise (`memcmp`), with the host shadow on every call, and only changed slabs are uploaded.

Why not (a), a registry field in `FVW_MiscVarType`:
- `InducedVelocitiesAll_Init` has `m` as `intent(in)` and `Calc`/`End` do not receive `m` at all: a field in `m` could only be reached through `Sgmt` anyway, so a `T_Sgmt` field would be needed;
- a registry field is **copied** by the generated `FVW_CopyMisc` (two instances share one handle: a double free at `FVW_End` or a cross-instance race) and **packed/unpacked** for checkpoints (a handle id from another process would alias a live handle);
- it changes `FVW_Types.f90`, so the L0 oracle (pin) and the port could no longer share the registry types in one process, and A5 (OFF build identical to the pin) gets a large generated diff.

| Hazard | What happens with the address-keyed map |
|---|---|
| Glue copies misc (`FVW_CopyMisc`) | The copy's `Sgmt` has a new address → a new, independent handle, created lazily. The original is untouched. L0 A7: copy and original used alternately, both bitwise equal to a fresh instance. |
| Glue copies state (`FVW_CopyContState`, predictor-corrector, `m%x1/x2`) | Irrelevant to identity: the lattice is validated by content on every call, whichever `x` is passed. |
| Re-evaluation / rollback at the same `t` (`FVW.f90:597-610`) | Same: content validation sees the restored lattice (L0 A7: x1, x2, x3, x1 again). No version counter, no assumption on `t`. |
| Several rotors / instances in one process | One handle per `Sgmt`; AeroDyn has one FVW instance for all rotors (`ad_QuadRotor_OLAF`); FAST.Farm/Python get one per OpenFAST instance. The map is mutex-guarded; a handle is used by its instance's thread only. (Concurrent Kokkos-OpenMP dispatch from several host threads, as FAST.Farm's OpenMP over turbines would do, is a Kokkos question, not a handle one; noted for unit 3.) |
| C binding loaded from Python | Kokkos is initialized lazily (`ofk_runtime`), the handle lazily; R6 runs `HelicalWake_OLAF_py`/`EllipticalWing_OLAF_py` through `libopenfastlib`. |
| `*_End` never called / copies destroyed by generated code | Handle orphaned: freed by the finalize hook at exit; the number of live handles is capped (`OFK_IV_MAX_HANDLES`, default 64, LRU eviction; an evicted live instance just rebuilds lazily). A reused address finds the stale handle, whose content is re-validated, so it is only a cache, never a source of wrong results (L0 with `OFK_IV_MAX_HANDLES=1`: all pass). |
| `InducedVelocitiesAll_End` has no `Sgmt` | It cannot release the source set; Calc instead checks `(method, nSeg | nPart)` against the handle's set and rebuilds from the host `Sgmt`/`Part` if they differ (`ofk_iv_sources_valid`). Init and LL always rebuild. |

### Source sets: one source, two paths (goal 2, F5, F1)

The Fortran pack (`PackPanelsToSegments`, `SegmentsToPartWrap`) still runs on the host in every mode: `Sgmt` and `Part` are API outputs of Init/LL (read by Calc, OnGrid_Calc, the trees, and copied by `FVW_CopyMisc`), so they stay exact. The device source set is built from one of:
- **lattice** (default on device backends): per wing, the six lattice arrays `r/Gamma/Eps × NW/FW` are synchronized slab by slab (depth slabs; only the active depths) into a flat mirror; then **one device kernel packs lattice → SoA segments** directly in the kernels' layout, with the mirror image, in the Fortran's segment order (closed-form `iHeadC`, checked against the Fortran's `nSeg` in every call; a mismatch is fatal), and a second kernel converts segments → particles (`k = j·nPartPerSeg + iSub`). This is F1 on the device, and F5(i): the source set is uploaded once and patched, not repacked and re-uploaded per call.
- **packed** (default on host backends, where it costs one O(N) copy and no transfer): the Fortran-packed `Sgmt`/`Part` are copied into the SoA block (host-parallel), and uploaded (device builds) column by column.

`OFK_IV_SOURCES=packed|lattice` switches (the unfused path stays available, criteria F5). Both paths are L0-tested on both backends.

| Fusion | Class | Evidence |
|---|---|---|
| F1 lattice → segments → particles (device) | **B** (elementwise, Fortran operation order, no contraction in any mode) | L0 strict lattice = Fortran, bitwise, all flag combinations; mutations M1, M2, M5 caught |
| F5(i) upload once + patch (slab validation) | **B** (a cache; the arithmetic is unchanged) | L0 A7; counters `lattice compared/changed/uploaded` |
| F5(ii) cached far-wake sum | not done (class R) | — |

### Evaluation

`ofk_iv_eval` runs unit 1's `run_seg_soa`/`run_part_soa` on the handle's SoA, in the mode `OFK_MODE` selects. Host backends: the kernel reads `CPs` and accumulates into `Uind` directly (unmanaged views on the Fortran memory, any leading dimension): no copy at all. Device backends: CPs are staged (pinned) and uploaded; `U` on the device starts at +0, the kernel computes `U = U + s = s` exactly, the result comes back to pinned staging and the host does `Uind = Uind + s`: the same single addition as the Fortran, without uploading `Uind`. The source-panel term (`ui_quad_src_nn`) stays Fortran and is added after, as in the Fortran.

## Audit list → L0 case

| # (olaf-placement §6) | Trap | What the port does | L0 case |
|---|---|---|---|
| 4 | `SegConnct` 1-based, rows 3/4 depth/span, mirror adds `nSegP` | Packed path offsets at use and validates every index; lattice path never builds connectivity (closed form) | every case, both paths; `Sgmt` compared bitwise |
| 5 | `Uind` accumulated; `m%Uind(:,nCPs+1:)=1000`; Vind sentinels -9999/2222 | Only columns `1:nCPs` touched; sentinels written by the unchanged Fortran | `wake Uind` (tail bitwise), `calc #1` (nonzero `Uind`, `-0`, NaN tail), `calc #2 strided` (skipped columns NaN), Vind_NW/FW full arrays |
| 9 | `LatticeToSegments` includes/omits shed segments by `bShedVorticity`, `bShedLastVorticity`, `LastNWShed`; FW eps hack | Closed-form offsets per block = `CountSegments`, checked every call; `hackEps` per block | `flags` cases: VelocityMethod 1..4 × FWShed T/F × mirror T/F × nNW<nNWMax or = (LastNWShed); `iDepthStart` 1/2/3, `nNW < iDepthStart`, `nFW=0`, empty wake, one-panel span; random cases |
| 11 | OLAF+OpenMP non-deterministic; Fortran `!$OMP SIMD` reorders | Oracle = pin, OPENMP=OFF flags | all `l0_iv_*` |
| 14 | `UpdateStates` twice at the same `t`; glue copies `x` | Content-validated mirror | A7 sequence |
| (new) | Calc without a preceding Init on this `Sgmt` | `ofk_iv_sources_valid` → rebuild from host arrays | (by construction; End has no `Sgmt`) |
| (new) | `type(*)` key / copies / orphaned handles | address-keyed, LRU, finalize hook | A7 copy cases; cap-1 run |

Mutation check (L0 must fail): M1 FW eps hack ignored → 494 failures; M2 `shedLast` ignored in the offsets → 413; M3 lattice validation skipped (stale mirror) → 617; M5 `DP*(Γ/n)` instead of `(DP·Γ)/n` in the particle pack → 41 (ulp-level). M4 (`(0.5+i)·L` → `0.5L + iL`) survived: an equivalent mutant (`0.5L` is exact, so both round once).

## Gate results

| ID | Result | Evidence |
|---|---|---|
| L0 claim 1 (Fortran vs strict, CPU) | **PASS** | bitwise, packed and lattice sources, all VelocityMethod x flags x RegFunction; `iv-omp` ctest |
| L0 claim 2 (strict CUDA) | PENDING | GPU jobs 1860223 / 1860779 |
| L0 claim 3 (strict vs fast) | **PASS** (CPU) | 8.2e-14 Basic, 3.1e-15 Part (bound 1e-12); Fortran's own FMA spread 8.2e-14 / 6.0e-14. CUDA PENDING |
| R6 | **PASS** (CPU) | job 1860600: 16/17 on all three, same failure; strict `.outb` = OFF in 17/17. CUDA PENDING |
| R7 | **PASS** (CPU) | in process and 1/4/8 threads (md5), both modes; system level: kk-fast 2 vs 8 threads identical (R8 job). CUDA PENDING |
| R8 | **FAIL by the letter / PASS in substance (CPU, tree wake)** | envelope OMP1 vs OMP8 = 0; kk-fast drift 1.9e-13, never > 1e-6 in 300 steps; direct case, FMA envelope and GPU PENDING |
| P6 | **PASS at L0 (CPU)**; full run PENDING | L0: 0 allocations in 20 steady steps. Full run with the old binary (doubling only): last allocation after eval 72/93 (lesson 3, fixed by capacity hints; re-measured by job 1860992 and the GPU jobs) |
| P1/P3 | partial (CPU) | `bar-perf-direct`, 8 threads: whole run 299.7 s (Fortran+OpenMP) vs 153.8 s (fast); fast profile: kernel 110.1 s, pack 0.1 s, host 46.5 s, no transfers on CPU. GPU split PENDING |
| P2 | **PASS** direct wake (1.95x); **FAIL** tree wake with the `OPENMP=OFF` fork build (2x slower: Fortran treecode loses OpenMP); re-measured with `fork-omp-fomp` by job 1860992 |
| P4 | marginal | r-test times (ctest -j 8, noisy): within 1.10x except `ad_EllipticalWingInf_OLAF` 0.39 → 0.49 s (Kokkos init) |
| A1 | **PASS** | no `_MOD_` removed; added `__fvw_subs_MOD_fvw_kokkos_release` (Kokkos builds only); C symbols `ofk_iv_*` |
| A2 | **PASS** | fork diff: `FVW_Subs.f90`, `FVW.f90` (FVW_End, 4 lines), unit-1 files, `modules/kokkos/`; every caller compiles unchanged; 3 configurations build |
| A5 | **PASS** | fork-off: `fvw_subs.mod`/`fvw.mod` md5 = pin; `.text` differs only in line-number immediates; `l0_iv_kkoff` bitwise |
| A7 | **PASS** (CPU) | L0 A7 sequence: re-evaluation, changed x, rollback, misc copy alternating with original, orphaned copy, LRU eviction (cap 1) |

## Lessons

Short version (the report of session `199ad883` has the long one).

1. **The Fortran signatures decide where per-instance state can live.** `Init` has `m` as `intent(in)`, `Calc` has no `m`, `End` has no `Sgmt`. The only object every entry point sees is `T_Sgmt`, so its address is the key. Rule for unit 3 and later: list, per entry point, which instance objects arrive and with which intent, before choosing the handle's home.
2. **Validate content, never identity.** A bitwise slab compare against a host shadow makes every A7 hazard (glue copies of `x`, rollback, re-evaluation, key reuse after a free, orphaned handles) a cache question instead of a correctness question, with no version counter to thread through callers. The cost is O(lattice) host reads per call, about the size of the upload it saves, so it pays only when most slabs are unchanged (LL #b → #d → next wake).
3. **"Grow by doubling" is not enough for P6 in a time loop whose problem grows.** OLAF's wake grows by one panel per step until `nNWMax`, so a doubling buffer reallocates for the first `nNWMax` steps (memtool: last allocation after evaluation 72 of 93 in a 30-step run). The Fortran already allocates its misc arrays at the maximum (`FVW_InitMiscVarsPostParam`): pass those extents as capacity hints on the first call (`ofk_iv_reserve`, `capCP`).
4. **API-visible intermediates block the host-side half of a fusion.** `Sgmt`/`Part` are `intent(inout)` outputs of `Init` and misc state of `LL`, so the Fortran pack must still run on the host even though the device now packs from the lattice (F1 is done on the device, but the host copy stays as the "side output" of criteria F1). That host pack is the GPU path's Amdahl term.
5. **`.outb` cannot carry R8 or claim 1.** It stores channels as 16-bit scaled integers: `.outb` identity means agreement to about 1.5e-5 of the channel range. R8 needs text output at full precision (`OutFmt "ES24.16E3"`, `OutFileFmt 1`); with it, strict = Fortran exactly, over 300 steps, on both BAR cases.
6. **R8's envelope can be zero.** OLAF's Fortran OpenMP loops are per target (fixed order), so `OMP_NUM_THREADS=1` vs `8` gives identical output on both BAR cases (tree and direct). Any reordered or FMA kernel then "fails" R8 at 2 × 0. A non-zero envelope of the same kind as the fast mode's difference is the pin built with FMA (`-march=x86-64-v3`).
7. **Build configuration is part of P2.** The fork's `OPENFAST_KOKKOS=ON` builds of unit 1 used `OPENMP=OFF` for the Fortran, so every non-ported OpenMP loop (the treecode) runs single-threaded: the tree-wake BAR case is 2× slower than the pin with OpenMP. Kokkos-OpenMP and Fortran OpenMP together (`fork-omp-fomp`) is the configuration to compare.
8. Smaller ones: concurrent runs with `OMP_PROC_BIND=spread` pin to the same cores (3× slowdown); Fortran `get_environment_variable` into a 64-character buffer silently truncates paths (unit 1's `l0_util.f90` has this; a stray file `opt/build/biots` is its trace); `type(*)` assumed-type dummies give one C entry per routine with a run-time real kind instead of `_r8`/`_r4` pairs; nvcc contracts by default, so every new TU with device arithmetic must be put in the strict or fast list explicitly.
