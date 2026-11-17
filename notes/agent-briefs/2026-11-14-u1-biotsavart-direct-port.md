<!-- Brief given to subagent af4ce10c (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are doing **unit 1 of 3** in the OpenFAST → C++/Kokkos port. Repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. The three units establish the port's architecture. The user wants us to **learn something from each unit**, so your lessons section matters as much as the code.

## Read first
1. `CLAUDE.md`, and `notes/architecture.md`, which is the design you are implementing; read all of it.
2. `notes/criteria.md`: the "Modes and the equivalence chain" paragraph and the class table.
3. `.claude/skills/ofk-port-unit/SKILL.md`.
4. `tests/biotsavart-direct/unit.md`. A previous agent already ported `ui_seg` and `ui_part_nograd` as v0, staged outside the fork. Its findings:
   - On Kokkos-OpenMP, L0 is bitwise in every regularization branch.
   - A1/A2/A5 pass; the pin's `FVW_RunTests` passes through the port.
   - **P2 fails** for `ui_part_nograd`: 0.48–0.91× of Fortran+OpenMP, because the Fortran vectorizes with `!$OMP SIMD reduction` (`FVW_BiotSavart.f90:385,430`).
   - **P6 fails**: device allocation on every call. That is out of scope here; unit 2 adds the handle.
   - `MINNORM` must be passed from the Fortran, because of `-fdefault-real-8`.
   - The class-B oracle is the pin built with OPENMP=OFF.
   - Linking NWTC needs `module load openblas/0.3.28` and `-DBLA_VENDOR=OpenBLAS`.

   v0 code is in `modules/kokkos/` and `port/modules/aerodyn/src/FVW_BiotSavart.f90` in the parent repo. v0 tests and build recipe are in `tests/biotsavart-direct/`. Build trees are in `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/` (`openfast-pin`, `biotsavart-direct-*`).

## What to do
**A. Move the port into the fork, in place.**
- `fork/` is a submodule: `MM0n5Ter/openfast`, local branch `kokkos`, at the same pin `12fa3fc` as `deps/openfast`.
- Put `ofk_kernels` (runtime plus this unit) in `fork/modules/kokkos/`.
- Put the shim **in place** in `fork/modules/aerodyn/src/FVW_BiotSavart.f90`: same module and public interface, `#ifdef OPENFAST_KOKKOS` → `ofk_` call, `#else` the original body.
- Once `tests/biotsavart-direct/` builds against the fork, delete the parent repo's `modules/` and `port/` directories. They are untracked staging.

**B. Integrate into the fork's own CMake** (decision #5: extend OpenFAST's CMake, never generate a new build system).
- Add `option(OPENFAST_KOKKOS ...)` (default OFF), `find_package(Kokkos 5.3...<6)`, `add_subdirectory(modules/kokkos)`.
- Define `OPENFAST_KOKKOS` for the AeroDyn Fortran sources, and link `aerodynlib` (and whatever needs it) to `ofk_kernels`.
- With the option **OFF**, the fork must build exactly like the pin.
- Build the full fork in three configurations, outside the repo under `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/fork-{off,omp,cuda}`:
  - OFF;
  - ON with `kokkos-5.3.0-omp`;
  - ON with `kokkos-5.3.0-cuda-h100`, compiled on the login node with `nvcc_wrapper`, `NVCC_WRAPPER_DEFAULT_COMPILER=g++`.
  - At least `aerodyn_driver`, `openfast`, and the AeroDyn C binding `.so` must build and link. Record any linker-language or Fortran↔C++ link issues; they will recur for every module.

**C. Retrofit to the architecture** (`notes/architecture.md` §2–3).
- One interaction body per kernel in `biotsavart_direct_body.hpp`, templated on the regularization branch. Square, skinny (team) and any other shape call that same body; there is no duplicated physics.
- Two translation units, `_strict.cpp` (FMA off: `-ffp-contract=off`, nvcc `--fmad=false`) and `_fast.cpp` (FMA on), selected at runtime by `OFK_MODE=strict|fast`. The default is fast, and the mode is printed once.
- **strict**: the Fortran's per-target serial order. It must be bitwise against the OPENMP=OFF Fortran on Serial/OpenMP.
- **fast**: a fixed blocked order, i.e. a constant number of partial sums or lanes, independent of thread count, so results are deterministic across runs **and across thread counts**. It must vectorize on CPU (check with `-fopt-info-vec` or by timing) and use FMA. You may change data layout, e.g. copy to SoA inside the call.
- Keep each regularization guard literally, including the NaN behaviour difference between Vatistas and the others.

**D. Gates** (record PASS/FAIL with evidence in `tests/biotsavart-direct/unit.md`; rewrite that card for v1 and keep a short v0 history).
- **Chain claim 1:** strict vs Fortran (OPENMP=OFF, pin flags) on Serial/OpenMP = 0, every branch, seeds 1 and 7.
- **Chain claim 3:** fast vs strict on OpenMP: worst error relative to Σ|contributions|. Declare the class (R1, bound 1e-12) and report the measured worst value.
- **R7:** both modes bitwise run to run; fast also bitwise across 1, 4 and 8 threads.
- **Chain claims 2 and 3 on CUDA, and R7 on GPU:** submit one `sbatch` with `-p h100,h200 --gres=gpu:1 -c 8 -t 00:20:00 -J ofk-biotsavart-direct-gpu`. Both partitions are sm_90; the queue can take hours.
  - Submit as soon as the CUDA build works, then continue with other work.
  - If it hasn't run when you are done, report the job id and log path; don't wait idle.
  - Job 1858779 (v0 GPU L0) is also pending. If it has finished by then, record its result as v0 history.
- **P2:** fast mode on Kokkos-OpenMP vs the Fortran built with `-fopenmp`, 8 threads, at 10⁴ and 10⁵ targets × sources, all regularization branches. Run it on the `med` partition with `-c 8`, because `interactive_cpu` caps users at 4 CPUs. Also report the strict-mode cost.
- **A1/A2/A5** on the fork build. **R6:** run `aerodyn_utest` and `ctest -R OLAF` on fork-OFF vs fork-OMP. This needs r-test data:
  - Initialize it once with `git -C deps/openfast submodule update --init --depth 1 reg_tests/r-test`, as `README.md` describes.
  - Reuse that copy for the fork, e.g. a symlink or a `--reference` clone. Do not commit it into the fork.
  - If r-test doesn't fit the time, say so and leave R6 for unit 2, which needs it anyway.

## Rules
- **Never edit `deps/openfast`.** Initializing its r-test submodule is allowed.
- **No git commits or pushes**, in the fork or the parent. The user hasn't asked for them yet.
- Toolchain: `source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28`. Kokkos 5.3.0 installs are under `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/`. Code is C++20. Compile everything on the login node.
- The PI-wide GPU limit is 16. Name every Slurm job `ofk-biotsavart-direct-<what>`.
- Write only in `fork/` and `tests/biotsavart-direct/`, plus the deletion of `modules/` and `port/` in the parent. If a shared doc (`notes/*`, skills, `CLAUDE.md`, `README.md`, `examples/`) should change, put the exact suggestion in your report instead of editing it.

## Report back (concise)
1. **Files changed or created**, and the fork diff summary (`git -C fork status --short`, `git -C fork diff --stat`).
2. **Gate table** (ID, PASS/FAIL, evidence), the chain-claim numbers, and P2 timings for strict and fast.
3. **Slurm job ids and log paths.**
4. **Lessons for units 2 and 3:** what in the architecture worked, what didn't, what you would change in `notes/architecture.md` or the skill, and what cost the most time. Be specific, with file:line where useful.
5. **Suggested exact edits to shared docs.**
6. **Rough wall-clock time spent.**
