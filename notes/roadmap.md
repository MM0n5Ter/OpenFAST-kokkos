# Roadmap: one month on OLAF, then the whole repository

The target is the whole repository Fortran-free ([decisions.md](decisions.md) #2). The first month ports OpenFAST's hottest algorithm, OLAF wake induction. Four people work in four lanes.

## Rules for parallel work
- Each unit keeps its interface, so a lane porting a caller never waits for the callee's port, and the reverse.
- Lanes conflict only on shared files: CMake files, registry `.txt` files and module-level variables. Each shared file has one owner. Lane A owns the build.
- One unit = one branch = one PR. The main branch always builds and passes L1. A unit merges only when all three gates pass.
- Every unit has a **unit card** (`tests/<unit>/unit.md`, template in [harness.md](harness.md)):
  - its files;
  - its contract;
  - its callers and the downstream consumers that touch it;
  - its L0–L3 tests;
  - owner and status;
  - agent cost (tokens, $, E_k) and human interventions.

## The month (4 people × 4 weeks)

```
            Week 1: foundations      Week 2: direct kernels     Week 3: treecode + gates    Week 4: results
          ┌──────────────────────┬──────────────────────────┬──────────────────────────┬─────────────────────────┐
 Lane A   │ Kokkos installs (OMP,│ ui_seg, ui_part_nograd   │ GPU tuning of direct     │ CPU/GPU scaling runs    │
 build/   │ CUDA); toy example   │ in Kokkos behind the     │ kernels; LL induction    │ vs Fortran+OpenMP;      │
 kernels  │ passes on both;      │ Fortran shims;           │ (block per CP);          │ cost/agent metrics      │
          │ OPENFAST_KOKKOS opt. │ aerodyn_utest passes     │ persistent device bufs   │                         │
          ├──────────────────────┼──────────────────────────┼──────────────────────────┼─────────────────────────┤
 Lane B   │ contract for         │ treecode design: linear  │ treecode in Kokkos behind│ treecode tuning;        │
 treecode │ InducedVelocitiesAll │ (Morton) octree, stack-  │ InducedVelocitiesAll_*;  │ tree-vs-direct accuracy │
          │ _Init/Calc/End       │ less traversal           │ tree-vs-direct tests     │ study                   │
          ├──────────────────────┼──────────────────────────┼──────────────────────────┼─────────────────────────┤
 Lane C   │ baseline r-test at   │ L0 harness: oracle build │ gate runner; L2 on all   │ regression report:      │
 tests    │ OMP=1; scaled OLAF   │ + capture/replay of      │ OLAF cases; L3 Python    │ pass rate, max rel err  │
          │ case; profiling      │ InducedVelocitiesAll_Calc│ C-binding cases          │                         │
          ├──────────────────────┼──────────────────────────┼──────────────────────────┼─────────────────────────┤
 Lane D   │ API surface + symbol │ InflowWind FlowField     │ IfW unit tests + ifw     │ animation: Fortran vs   │
 IfW +    │ baselines; VTK→PNG→  │ GetVelAcc in Kokkos (at  │ r-tests; ADI C binding / │ Kokkos wake side by     │
 anim.    │ MP4 pipeline         │ wake points)             │ pyOpenFAST check         │ side, with timers       │
          └──────────────────────┴──────────────────────────┴──────────────────────────┴─────────────────────────┘
 Milestones: M1 (end W1) baseline numbers + toy passes on both Kokkos builds │ M2 (W2) direct kernels pass L0 + L1
             M3 (W3) OLAF unit passes API + regression gates                 │ M4 (W4) performance + animation + write-up
```

**Hand-written lines touched in the month:**
- OLAF induction core: about 3.6k (`FVW_BiotSavart` 724, `FVW_VortexTools` 2,424, induction drivers in `FVW_Subs` about 450).
- InflowWind `IfW_FlowField`: about 2.2k.

**Definition of done**
1. **API:** the OLAF r-test cases, `aerodyn_utest` and the ADI C binding run with zero source edits outside the unit's files. Exported symbols and output headers are unchanged.
2. **Regression:** L0 passes. `aerodyn_utest` passes. All 16 OLAF r-test cases pass at r-test tolerances.
3. **Performance:** time per step against the number of wake markers, comparing Fortran+OpenMP, Kokkos-OpenMP and Kokkos-CUDA, on the scaled case.
4. **Animation (add-on):** one wake case rendered from both builds.

**Fallbacks and risks**
- **Treecode fallback.** If the Kokkos treecode is not ready by M3, ship the direct N² sums on GPU (`VelocityMethod` 1 or 3, both valid OLAF options). Report where they cross over the CPU treecode, and finish the tree in week 4.
- **Non-determinism.** OLAF with OpenMP is not deterministic, so baselines use `OMP_NUM_THREADS=1`.
- **Tiny r-test cases.** They cannot show speed, so timing needs our scaled case.
- **Build risk.** The Kokkos build is unverified until M1. The fallbacks are in [build-kokkos-fortran.md](build-kokkos-fortran.md) §5.

## After the month: the whole repository

| Phase | Units | Hand-written lines | Notes |
|---|---|---|---|
| 1 | Rest of OLAF/FVW; all of InflowWind; **registry retarget** (57 registry files → C++ types + interop) | ~20k + tool | The registry retarget unblocks every later module |
| 2 | FAST.Farm AWAE + WakeDynamics; SeaState; NWTC hot code (ModMesh, ModMesh_Mapping, NWTC_Num, interpolation) | ~45k | FAST.Farm wake planes give a farm-scale animation |
| 3 | Full turbine: AeroDyn rest (BEMT, UA, DBEMT, AirfoilInfo, AeroAcoustics), ElastoDyn, BeamDyn, ServoDyn (keep the controller DLL API), HydroDyn, MoorDyn, SubDyn; glue `openfast-library` (keep the `FAST_Library` C API) | ~215k | Gate: the full r-test suite |
| 4 | Small modules, TurbSim, FAST.Farm glue, the rest of NWTC (IO, system files); delete the shims | ~100k | `map` is already MAP++ (C), so only its Fortran wrapper is ported |

For phases 2–4 we report the **rate** measured in the month and in phase 1 (kLOC per week per lane, $ per kLOC, interventions per unit), and project the cost of the rest.
