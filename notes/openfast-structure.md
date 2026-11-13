# OpenFAST structure: what the units of work are

Measured on `deps/openfast` @ `12fa3fc` with `scripts/`.

## Size

`scripts/loc-by-module.sh` (lines; "types" = generated `*_Types.f90`; "reg" = registry files):

| Directory | Hand-written Fortran | Types | C/C++ | reg |
|---|---|---|---|---|
| nwtc-library | 64,702 | 1,973 | 4,523 | 3 |
| aerodyn | 49,586 | 29,427 | 0 | 10 |
| hydrodyn | 48,370 | 16,837 | 0 | 8 |
| subdyn | 37,115 | 5,486 | 0 | 1 |
| openfast-library (glue) | 20,602 | 12,001 | 85 | 2 |
| glue-codes/fast-farm | 20,306 | 2,692 | 0 | 2 |
| moordyn | 17,174 | 5,915 | 0 | 1 |
| elastodyn | 16,870 | 8,212 | 0 | 1 |
| beamdyn | 15,042 | 3,897 | 0 | 1 |
| servodyn | 14,828 | 10,213 | 0 | 2 |
| turbsim | 14,214 | 299 | 0 | 0 |
| inflowwind | 13,581 | 5,392 | 0 | 1 |
| seastate | 12,391 | 3,147 | 0 | 0 |
| small modules and glue (awae, icefloe, icedyn, simple-elastodyn, soildyn, aerodisk, feamooring, extptfm, wakedynamics, map wrapper, externalinflow, extloads, lindyn, orcaflex-interface, version, labview, openfast) | ~44k | ~35k | map: 39k C (MAP++) | |
| **Total** | **388,695** | **140,499** | **53,094** | **57** |

C/C++ already in the tree: `openfast-registry` (the type generator, 4.2k), `glue-codes/openfast-cpp` (3.4k), MAP++ (C, 39k), Simulink glue, `awae` C helpers.

## Every module uses one framework API
- Each module exposes `<Mod>_Init / _UpdateStates / _CalcOutput / _End`, plus `_CalcContStateDeriv` and the Jacobian routines.
- Their arguments are registry-typed: `InitInput, u, p, x, xd, z, OtherState, m, y`.
- `scripts/api-surface.sh` finds the full set in 36 modules and submodules: AD, ADI, BEMT, DBEMT, UA, AA, FVW, ED, BD, SrvD, StC, HydroDyn, SeaSt, MD, SD, IfW, AWAE, WD, ….
- So one contract template covers every module, and the C++ side of each contract is generated once the registry is retargeted.

## Dependency graph (`scripts/module-deps.py`)

```
 nwtc-library   fan-in 28: everyone
 version        fan-in 14
 inflowwind     fan-in 7: aerodyn, aerodisk, awae, seastate, externalinflow, openfast-library, fast-farm
 seastate       fan-in 5: aerodyn, hydrodyn, moordyn, openfast-library, labview
 aerodyn        fan-in 3: openfast-library, fast-farm, labview
                fan-out: inflowwind, lindyn, seastate, nwtc-library, version
 every other physics module: nwtc-library (+ version) only
 openfast-library (glue): fan-out 22, everything
```

Physics modules are independent of each other, so they can be ported in any order and in parallel.

## Inside AeroDyn *(v5.0.0, re-measure on the pin)*

| Part | Lines | Depends on |
|---|---|---|
| AeroDyn core (`AeroDyn.f90`, IO, rotor/tower/nacelle) | ~14.1k | everything below |
| BEMT | ~4.1k | UA, DBEMT, AFI |
| DBEMT | 0.9k | — |
| UnsteadyAero (UA) | 4.3k | AFI |
| AirfoilInfo (AFI) | 2.2k | stateless library |
| **FVW (OLAF)** | ~8.3k | induction core (~3.4k) is a leaf; `FVW_Wings` uses UA and AFI |
| AeroAcoustics | ~3.4k | — |
| ADI wrapper + C binding | ~3.4k | AeroDyn + InflowWind |

AeroDyn picks the wake model with `Wake_Mod`: 0 none, 1 BEMT, 3 FVW/OLAF.

## OLAF: from the entry points to the hot kernels

OLAF files on the pin:
- `FVW.f90` (1,873 lines);
- `FVW_Subs.f90` (2,224);
- `FVW_VortexTools.f90` (2,424; the treecode, 38 OpenMP directives);
- `FVW_BiotSavart.f90` (724; the direct kernels, 37 OpenMP directives);
- `FVW_Tests.f90` (961);
- `FVW_IO.f90`, `FVW_Wings.f90`, `FVW_Types.f90` (generated).

**Layer 1: who calls OLAF.** Only AeroDyn, and only when `Wake_Mod = 3`:

| AeroDyn (`AeroDyn.f90`) | calls |
|---|---|
| `AD_Init` → `Init_OLAF` (:530) | `FVW_Init` (:5722) |
| `AD_UpdateStates` | `SetInputsForFVW` (:2027), `FVW_UpdateStates` (:2031) |
| `AD_CalcOutput` path | `SetInputsForFVW` (:2378), `FVW_CalcOutput` (:2381), `SetOutputsFromFVW` (:2384) |
| `AD_End` | `FVW_End` (:1891) |
| `AD_JacobianPInput` | `FVW_CalcOutput` (:7636, :7655) |

AeroDyn itself is reached from the OpenFAST glue (`openfast-library`), the AeroDyn driver, the ADI C binding (`AeroDyn_Inflow_C_Binding.f90`, used by pyOpenFAST), FAST.Farm and the LabVIEW wavetank.

**Layer 2: OLAF's framework API** (`FVW.f90`):
- `FVW_Init` :40, `FVW_End` :522;
- `FVW_UpdateStates` :567, `FVW_CalcContStateDeriv` :801, `FVW_Euler1` :1060;
- `FVW_CalcConstrStateResidual` :1353, `FVW_CalcOutput` :1472.

**Layer 3: the path to the hot core**

```
FVW_UpdateStates (FVW.f90:567)
 ├─ time integration: FVW_Euler1 (:1060, default) or RK4
 │   └─ FVW_CalcContStateDeriv (:801)
 │       └─ WakeInducedVelocities (called :855; FVW_Subs.f90:1379)        "expensive N^2 call"
 │           └─ InducedVelocitiesAll_Init / _Calc / _End (FVW_Subs.f90:1237 / :1285 / :1334)
 │               ├─ VelocityMethod 1: ui_seg            (FVW_BiotSavart.f90:135)   direct segments
 │               ├─ VelocityMethod 2: grow_tree_part (FVW_VortexTools.f90:557)
 │               │                    + ui_tree_part (:1512)                         particle treecode
 │               ├─ VelocityMethod 3: ui_part_nograd    (FVW_BiotSavart.f90:360)   direct particles
 │               ├─ VelocityMethod 4: grow_tree_segment (:947) + ui_tree_segment     segment treecode
 │               └─ source panels:    ui_quad_src_nn    (FVW_BiotSavart.f90:687)
 ├─ wake propagation, near-wake → far-wake bookkeeping
 └─ circulation solve: FVW_CalcConstrStateResidual (:1353) → FVW_Wings (uses UA, AFI)
     └─ LiftingLineInducedVelocities (called :1457, :1749; FVW_Subs.f90:1488) → ui_seg
FVW_CalcOutput (:1472)
 └─ InducedVelocitiesAll_OnGrid_Calc (called :1633; FVW_Subs.f90:1133)      grid outputs
```

**Porting boundary for the month:** `InducedVelocitiesAll_Init/Calc/End` and the direct kernels.
- Every route into the expensive math goes through them: wake self-induction, lifting-line induction and grid outputs.
- Their inputs are flat arrays (`Sgmt%Points`, `Sgmt%Connct`, `Sgmt%Gamma`, `Part%P`, `Part%Alpha`, …) plus scalars, and their output is `Uind(3, nCPs)`.
- The tree types `T_Tree`/`T_Node` (`FVW_VortexTools.f90:67-88`) are pointer-based and recursive. They are internal to the unit and are replaced by a flat (Morton-ordered) octree on the C++ side.
- **Things a port must preserve** (from the source comments and code):
  - `Uind` is *accumulated* into ("Side effects!!!");
  - arrays are "3 x nCPs++": only the first `nCPs` columns may be touched;
  - `SegConnct` holds **1-based** point indices.

## External API surface (`scripts/api-surface.sh`)

`bind(C)` symbols per file:

| File | C symbols | Used by |
|---|---|---|
| `openfast-library/src/FAST_Library.f90` | 22 | openfast-cpp (AMR-Wind/ExaWind), Simulink, LabVIEW, pyOpenFAST (`fastlib`) |
| `map/src/map.f90` | 17 | MAP++ (C) |
| `seastate/src/SeaState_C_Binding.f90` | 14 | pyOpenFAST |
| `aerodyn/src/AeroDyn_Inflow_C_Binding.f90` | 12 | pyOpenFAST, LabVIEW wavetank; **the path into OLAF from Python** |
| `inflowwind/src/IfW_C_Binding.f90` | 6 | pyOpenFAST |
| `moordyn`, `hydrodyn` C bindings | 5 each | pyOpenFAST |
| `glue-codes/labview/src/WaveTank.f90` | 4 | LabVIEW |
| `servodyn/src/BladedInterface.f90` | 1, plus the DLL call | ROSCO and other controllers |
