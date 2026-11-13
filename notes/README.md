# Notes

These are working notes. Numbers were measured on `deps/openfast` at `12fa3fc` with the scripts in `../scripts/`, unless a line is marked *(v5.0.0)* or *(verify)*.

- [decisions.md](decisions.md): the story, design decisions and open TODOs
- [landscape.md](landscape.md): related work, including GPU work on OpenFAST/OLAF, and where we fit
- [openfast-structure.md](openfast-structure.md): module sizes, dependency graph, the framework API, inside AeroDyn, and OLAF's entry points down to the hot kernels
- [downstream.md](downstream.md): projects built on OpenFAST, how each one couples to it, and the API surface each one freezes
- [testing.md](testing.md): OpenFAST's existing tests, our test levels L0–L3, and why the tests are not ported first
- [roadmap.md](roadmap.md): the one-month plan (4 people × 4 weeks) and the long-term phases to a Fortran-free repo
- [build-kokkos-fortran.md](build-kokkos-fortran.md): **how to build Kokkos and link it with OpenFAST's Fortran**, and the rules at the Fortran/C++ boundary
- [harness.md](harness.md): **how to build the test harness for each translated file or library**: L0 oracle tests, regression, API checks and performance runs
