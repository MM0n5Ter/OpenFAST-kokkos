# OpenFAST-kokkos

An agentic, **API-preserving**, repository-level migration of [OpenFAST](https://github.com/OpenFAST/openfast) from Fortran to C++/Kokkos.
- The agent replaces OpenFAST one API unit at a time, behind frozen interfaces. It starts with the hottest algorithm: OLAF free-vortex-wake induction (Biot–Savart direct sums and the particle treecode).
- After every step the repository still builds, and its callers and downstream codes keep working with zero source edits: the OpenFAST glue, the module drivers, the C bindings used by pyOpenFAST, the `FAST_Library` API used by AMR-Wind/ExaWind, Simulink and LabVIEW, the controller DLL interface used by ROSCO, and the input and output files read by WEIS.
- Every step must pass three gates: **API**, **Regression** and **Performance**.
- The one-month target is the OLAF unit on CPU and GPU from one Kokkos source. The end state is a Fortran-free OpenFAST.

| Path | Contents |
|---|---|
| `notes/` | Working notes: decisions, landscape, OpenFAST structure, downstream consumers, testing, roadmap, **how to build Kokkos with Fortran**, **how to build the test harness** ([index](notes/README.md)) |
| `fork/` | Submodule: the fork [`MM0n5Ter/openfast`](https://github.com/MM0n5Ter/openfast), branch `kokkos`, started at the same pin. **The port happens here, in place.** Architecture: [notes/architecture.md](notes/architecture.md). |
| `deps/openfast/` | Submodule: upstream OpenFAST, pinned to `dev` `12fa3fc` (2026-11-06). It is the frozen baseline and the oracle for every comparison. **Never edit it.** |
| `scripts/` | Analysis scripts that produce the numbers in `notes/`. Each one defaults to `deps/openfast`. |
| `examples/fortran-kokkos/` | Minimal template of the pattern every unit follows: the frozen Fortran, the port behind the same interface, the C++/Kokkos kernel, an L0 differential test and a shared-library (C binding) check |

Planned top-level directories, added when their work starts (see [roadmap](notes/roadmap.md)):
- `tools/`: gate runner (`api`, `reg`, `perf`), capture/replay probes, interface extractor.
- `tests/`: L0 tests per unit, the scaled OLAF performance case, the downstream (L3) workflows.
- (the fork is now `fork/`, see the table above)

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
| `scripts/kernel-inventory.py [--md] [--json out] [FILES]` | Per routine: size, loop nest, wake-scale loops, OpenMP, allocations, serial `iHead` counters, recursion/pointer trees, and a CPU/GPU placement hint (default: the OLAF files) |
| `scripts/literal-kinds.py [FILES]` | Default-kind real literals, whose value depends on compile flags. OpenFAST's double build uses `-fdefault-real-8`, so `MINNORM=1e-4` is exactly 1e-4 there. A port passes each one from the Fortran instead of hard-coding it (criteria R3). |
| `tools/cost/session_cost.py`, `tools/cost/compute_hours.sh` | **Cost ledger.** Money (tokens → USD API-equivalent, GPU/CPU hours) and time (agent wall/active hours, human prompts) per session, appended to `tools/cost/ledger.csv` ([method](tools/cost/README.md)) |

## Status

See [roadmap](notes/roadmap.md).
- **2026-11-14:** OLAF analysed routine by routine: what runs on the CPU, what runs on the GPU, and what to fuse ([olaf-placement.md](notes/olaf-placement.md)). The gates are refined into numbered checks ([criteria.md](notes/criteria.md)). Agent skills are in `.claude/skills/` (`ofk-analyze-unit`, `ofk-port-unit`).
- **Kokkos 5.3.0** (decision #6) was built from the release tarball, as two shared C++20 installs under `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/`: `kokkos-5.3.0-omp` and `kokkos-5.3.0-cuda-h100` (sm_90).
- `examples/fortran-kokkos` on 5.3.0:
  - **passes on Kokkos-OpenMP**: L0 206/206 at relative error 0 (bitwise), and ctest 2/2;
  - **passes on the H100** (job 1858536): worst relative error 2.3e-18, `capi_smoke` OK.
  - Earlier, on 4.5.1, it passed on the H100: worst relative error 2.3e-18.
- **Unit 1 of 3 done (v1)** ([unit card](tests/biotsavart-direct/unit.md)): the direct Biot–Savart kernels `ui_seg` and `ui_part_nograd`. They are in the fork behind `-DOPENFAST_KOKKOS=ON` ([architecture](notes/architecture.md)).
  - **strict** mode is bitwise to the Fortran on CPU. On the H100 it differs only through the device `exp` (≤ 9e-16).
  - **fast** mode runs 1.07–3.95× Fortran+OpenMP on 8 CPU cores, and up to ~120× on one H100.
  - In strict mode, the r-test OLAF outputs are identical to the OFF build.
  - Open:
    - the FMA policy for fast mode ([decisions](notes/decisions.md));
    - per-call device allocation, which unit 2 fixes.
- **Unit 2 in progress:** `InducedVelocitiesAll_*` and lifting-line induction, with persistent device buffers and the first end-to-end and R8 runs.
