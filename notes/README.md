# Notes

These are working notes. Numbers were measured on `deps/openfast` at `12fa3fc` with the scripts in `../scripts/`, unless a line is marked *(v5.0.0)* or *(verify)*.

- [decisions.md](decisions.md): the story, design decisions and open TODOs
- [landscape.md](landscape.md): related work, including GPU work on OpenFAST/OLAF, and where we fit
- [openfast-structure.md](openfast-structure.md): module sizes, dependency graph, the framework API, inside AeroDyn, and OLAF's entry points down to the hot kernels
- [downstream.md](downstream.md): projects built on OpenFAST, how each one couples to it, and the API surface each one freezes
- [testing.md](testing.md): OpenFAST's existing tests, our test levels L0–L3, and why the tests are not ported first
- [roadmap.md](roadmap.md): the one-month plan (4 people × 4 weeks) and the long-term phases to a Fortran-free repo
- [build-kokkos-fortran.md](build-kokkos-fortran.md): **how to build Kokkos and link it with OpenFAST's Fortran**, and the rules at the Fortran/C++ boundary
- [criteria.md](criteria.md): **the gates as numbered checks**: API (A1–A8), Regression with equivalence classes B/R1/R2/S (R1–R7), Performance (P1–P6), Placement (L1–L4), Fusion legality (F1–F5), effort metrics (E1–E4)
- [architecture.md](architecture.md): **the shape of every ported unit**: code layout (fork, `ofk_kernels`), one source with strict and fast modes, the equivalence chain, per-instance device handles, and the three units that establish it
- [olaf-placement.md](olaf-placement.md): **OLAF, CPU vs GPU per routine**: the step anatomy with calls per step, the host↔device boundary, kernel shapes, the fusion plan (F1–F8) and the porting audit list
- [harness.md](harness.md): **how to build the test harness for each translated file or library**: L0 oracle tests, regression, API checks and performance runs
