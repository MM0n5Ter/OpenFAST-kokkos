# OpenFAST-kokkos: instructions for agents

API-preserving migration of OpenFAST (Fortran) to C++/Kokkos, one API unit at a time. Start with `README.md` and `notes/README.md`.

- `deps/openfast` is the frozen baseline and oracle (pin `12fa3fc`). **Never edit it.** Cite code as `file:line` on the pin.
- The port lives in this repo as an **overlay** on the pin (`notes/overlay.md`). Changed OpenFAST files are full copies at their own paths under `src/modules/<mod>/src/`, with the Fortran shims in place. The C++/Kokkos code lives in `src/modules/kokkos/`. Build with `tools/openfast-build.sh <config> <builddir>`. `fork/` is the retired layout: do not use it. Unit cards and L0 tests go in `tests/<unit>/`. Read `notes/architecture.md` before porting: one source, strict and fast modes, the equivalence chain, per-instance handles.
- Every unit passes the gates in `notes/criteria.md` (API A*, Regression R*, Performance P*, Placement L*, Fusion F*). Record results in `tests/<unit>/unit.md`.
- Skills:
  - `ofk-analyze-unit`: placement and fusion analysis before porting;
  - `ofk-port-unit`: kernel, shim, L0 test and gates.
- Worked example of an analysis: `notes/olaf-placement.md`. Template for a port: `examples/fortran-kokkos/`.
- Toolchain: `module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0`. Kokkos **5.3.0** installs are in `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/` (`kokkos-5.3.0-omp`, `kokkos-5.3.0-cuda-h100`, and `kokkos-5.3.0-cuda-a100` for quick correctness runs on `interactive_gpu`). Code is C++20 and uses `find_package(Kokkos 5.3...<6)`. Compile CUDA on the login node; only GPU *runs* need `sbatch`.
- Build trees go outside the repo or in `build*/` (git-ignored).
- Numbers in `notes/` come from `scripts/`. Re-run the script when you change a number.
- **Track cost in money and time.** End every session with a row in `tools/cost/ledger.csv` (`tools/cost/session_cost.py … --ledger tools/cost/ledger.csv --unit <unit> --phase <phase>`; see `tools/cost/README.md`). Name Slurm jobs `ofk-<unit>-<what>` so GPU-hours can be counted. Report `$` as "USD (API-equivalent)".
