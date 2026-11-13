# OpenFAST-kokkos

An agentic, **API-preserving**, repository-level migration of [OpenFAST](https://github.com/OpenFAST/openfast) from Fortran to C++/Kokkos.
- The agent replaces OpenFAST one API unit at a time, behind frozen interfaces. It starts with the hottest algorithm: OLAF free-vortex-wake induction (Biot–Savart direct sums and the particle treecode).
- After every step the repository still builds, and its callers and downstream codes keep working with zero source edits: the OpenFAST glue, the module drivers, the C bindings used by pyOpenFAST, the `FAST_Library` API used by AMR-Wind/ExaWind, Simulink and LabVIEW, the controller DLL interface used by ROSCO, and the input and output files read by WEIS.
- Every step must pass three gates: **API**, **Regression** and **Performance**.
- The one-month target is the OLAF unit on CPU and GPU from one Kokkos source. The end state is a Fortran-free OpenFAST.

| Path | Contents |
|---|---|
| `notes/` | Working notes: decisions, landscape, OpenFAST structure, downstream consumers, testing, roadmap, **how to build Kokkos with Fortran**, **how to build the test harness** ([index](notes/README.md)) |
| `deps/openfast/` | Submodule: upstream OpenFAST, pinned to `dev` `12fa3fc` (2026-11-06). It is the frozen baseline and the oracle for every comparison. **Never edit it.** |
| `scripts/` | Analysis scripts that produce the numbers in `notes/`. Each one defaults to `deps/openfast`. |
| `examples/fortran-kokkos/` | Minimal template of the pattern every unit follows: the frozen Fortran, the port behind the same interface, the C++/Kokkos kernel, an L0 differential test and a shared-library (C binding) check |

Planned top-level directories, added when their work starts (see [roadmap](notes/roadmap.md)):
- `tools/`: gate runner (`api`, `reg`, `perf`), capture/replay probes, interface extractor.
- `tests/`: L0 tests per unit, the scaled OLAF performance case, the downstream (L3) workflows.
- The OpenFAST fork where the port happens. Where it lives is an open decision (see [decisions](notes/decisions.md)).

## Getting started

```bash
git clone --recurse-submodules https://github.com/MM0n5Ter/OpenFAST-kokkos.git
# or, in an existing clone:
git submodule update --init --depth 1 deps/openfast
# regression-test data (large; needed only to run r-test). Uses the commit OpenFAST pins (dfb039d):
git -C deps/openfast submodule update --init --depth 1 reg_tests/r-test
```

The submodule is shallow (`shallow = true` in `.gitmodules`), which is enough for reading and building the pinned tree.

Building: [notes/build-kokkos-fortran.md](notes/build-kokkos-fortran.md). Testing: [notes/harness.md](notes/harness.md).

## Scripts

| Script | What it reports |
|---|---|
| `scripts/loc-by-module.sh` | Hand-written Fortran, registry-generated `*_Types.f90`, C/C++ lines and registry files per module |
| `scripts/module-deps.py [--edges] [--json out]` | Module dependency graph from Fortran `use` statements: fan-in, fan-out and dependencies per module |
| `scripts/callers.sh [-o OWNER] ROUTINE…` | For each routine: files and module directories outside its owner (default `aerodyn`) that call it |
| `scripts/rtest-inventory.sh [PATTERN]` | Regression cases per registration function, and the cases whose name matches a pattern (default `olaf`) |
| `scripts/api-surface.sh [-v]` | External API in the source: `bind(C)` symbols per file, framework entry points per module, the controller DLL interface |

## Status

Week 1 (foundations) has not started. See [roadmap](notes/roadmap.md).
- The Fortran side of `examples/fortran-kokkos` was compiled and run against a plain C++ stand-in for the kernel.
- The Kokkos side has **not been built yet**. The build commands in the notes are unverified until week 1.
