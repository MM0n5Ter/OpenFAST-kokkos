# `src/`: the port, as an overlay on the pinned OpenFAST

`deps/openfast` (pin `12fa3fc`) is never edited. The port is the files below, applied on top of it at configure time (`notes/overlay.md`; build with `tools/openfast-build.sh <config> <builddir>`). There are three kinds of file:

```
src/
├── cmake/                       build integration (no OpenFAST code)
│   ├── ofk_overlay.cmake        injected with -DCMAKE_PROJECT_INCLUDE: swaps the overlay files into
│   │                            OpenFAST's module targets and links ofk_kernels
│   ├── ofk_options.cmake        every OFK_* build option, in one place
│   └── ofk_sources.cmake        for test trees: "the port's copy of <file>" (overlay if present, else the pin)
│
├── modules/<mod>/src/*.f90      Fortran overlay: full copies of OpenFAST files at OpenFAST's own paths
│                                (the path must mirror deps/openfast; that is how the overlay finds them)
│
└── modules/kokkos/              new code: the ofk_kernels library (C++20, Kokkos 5.3)
    ├── CMakeLists.txt           stand-alone / test-tree entry; layout comment
    ├── ofk_kernels.cmake        the target, and the strict/fast flag set of every translation unit
    └── src/
        ├── common/              shared by every unit: Kokkos lifetime, OFK_MODE (ofk_runtime), vector libm
        └── <unit>/              one directory per ported unit, same name as tests/<unit>/
```

## Fortran overlay files

Every change sits under `#ifdef OPENFAST_KOKKOS`. Preprocessed without it, each file is the pin's (gate A2).

| File (`modules/aerodyn/src/`) | Changed by unit(s) | What the shim does |
|---|---|---|
| `FVW_BiotSavart.f90` | biotsavart-direct | `ui_seg`, `ui_part_nograd` bodies call `ofk_biotsavart_*` |
| `FVW_Subs.f90` | induced-velocities, treecode, circulation-solve | `InducedVelocitiesAll_*`, `LiftingLineInducedVelocities`, the particle tree path, and the release of the per-instance handles (`FVW_Kokkos_Release`) |
| `FVW.f90` | induced-velocities | `FVW_End` releases the instance's handles |
| `FVW_Wings.f90` | circulation-solve | the `Vvar` block of `Wings_ComputeCirculationPolarData` calls `ofk_cs_*` |

## C++ units (`modules/kokkos/src/<unit>/`)

| Directory | Ports | Main files | Card |
|---|---|---|---|
| `common/` | — | `ofk_runtime.{hpp,cpp}`, `ofk_vector_libm.hpp` | — |
| `biotsavart-direct/` | `FVW_BiotSavart::ui_seg`, `ui_part_nograd` (direct Biot–Savart kernels, also used by the next two units) | `biotsavart_direct.cpp` (C API), `_body.hpp` (physics), `_strict.cpp` / `_fast.cpp` | `tests/biotsavart-direct/unit.md` |
| `induced-velocities/` | `InducedVelocitiesAll_*`, `LiftingLineInducedVelocities` | `induced_velocities.cpp` (handles, lattice mirror, sources), `_pack.cpp` (device packing) | `tests/induced-velocities/unit.md` |
| `treecode/` | `grow_tree_part`, `ui_tree_part` | `treecode_impl.hpp`, `treecode_dbuild.hpp` (device build), `_strict.cpp` / `_fast.cpp` | `tests/treecode/unit.md` |
| `circulation-solve/` | the `Vvar` block of the circulation solve | `circulation.cpp` (handles, C API), `_body.hpp`, `_strict.cpp` / `_fast.cpp` / `_fastfma.cpp` | `tests/circulation-solve/unit.md` |

Conventions (`notes/architecture.md`):
- **One body, several translation units.** `<unit>_body.hpp` holds the physics, and each `<unit>_<mode>.cpp` compiles it with one flag set. Every TU must be listed as strict or fast in `ofk_kernels.cmake`, because nvcc contracts to FMA by default.
- **Includes.** The include root is `modules/kokkos/src/`. A cross-unit include names its directory (`#include "biotsavart-direct/biotsavart_direct_body.hpp"`); same-directory includes stay bare.
- **C API.** Each unit's entry points are `extern "C" ofk_<prefix>_*`, declared in its `<unit>.hpp` and called by the Fortran shims through `bind(C)` interfaces.
- **A new unit** gets `src/modules/kokkos/src/<unit>/` and `tests/<unit>/`, its TUs in `ofk_kernels.cmake`, and a row in both tables above.
