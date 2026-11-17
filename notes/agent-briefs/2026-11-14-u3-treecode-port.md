<!-- Brief given to subagent aeb5f78a (session 199ad883), verbatim. Cost: tools/cost/ledger.csv. -->

You are doing **unit 3 of 3** in the OpenFAST → C++/Kokkos port. Repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. The three units establish the port's architecture, and the user wants us to **learn something from each unit**, so your lessons section matters as much as the code.

## Read first
1. `CLAUDE.md`.
2. **`notes/architecture.md`**, all of it. It includes the lessons of units 1 and 2:
   - ODR namespaces;
   - the FMA finding;
   - the CPU vectorization recipe;
   - the handle keyed by `m%Sgmt`, with content re-checks and first-call sizing;
   - the like-for-like P2 build;
   - the R8 envelope.
3. `notes/criteria.md`, especially class **R2** (approximation) and R8.
4. `.claude/skills/ofk-port-unit/SKILL.md`.
5. `notes/olaf-placement.md`: §1 (wake induction ×4 with RK4, default `VelocityMethod(1)` = particle treecode), §3 (tree shape), §6 audit items 8, 12 and 13 (tree threshold `1.0e-5`, multipole order-2 loop form "BAR_OLAF fails with unrolled version", `T_Tree` pointers into `Part`).
6. Unit 2's card `tests/induced-velocities/unit.md` and its code `fork/modules/kokkos/src/induced_velocities.*`. You extend its handle. Unit 1's card and code are `tests/biotsavart-direct/`, `fork/modules/kokkos/src/biotsavart_direct_*`.

## State
- `fork/` is the submodule `MM0n5Ter/openfast`, local branch `kokkos`, at pin `12fa3fc`.
  - It holds **uncommitted** unit 1 and 2 changes. Keep them working: unit 1 and unit 2 L0 and their r-tests must still pass after your changes.
- `deps/openfast` is the frozen pin. **Never edit it.**
- Unit 2 status: all CPU gates pass.
  - Strict is bitwise.
  - Fast is ≤ 8e-14 against strict.
  - R8 passes against the FMA-built pin's envelope.
  - P2: 2.06× on the direct wake, and level (1.02×) on the tree wake, because the tree is still Fortran.
  - P6 passes.
  - GPU correctness is running now on an A100 MIG slice (job 1862373). Its results may reach you by message.
- **Open user decisions; treat as defaults, don't spend time on them:**
  - FMA stays on in fast mode, with the bound = 1.1 × the Fortran's own FMA spread where the Fortran itself moves.
  - The R8 envelope is the pin built with FMA when the OMP_1-vs-OMP_N envelope is 0.

## Unit 3: `treecode`
The OLAF particle treecode in `fork/modules/aerodyn/src/FVW_VortexTools.f90` (pin lines):
- `grow_tree_part` (:557) with its `_rec`, `_substep` and `_parallel` helpers;
- `ui_tree_part` (:1512) / `ui_tree_part_11`;
- whatever `InducedVelocitiesAll_Init/Calc` (`FVW_Subs.f90`) need to dispatch to it;
- the multipole moments.

This is OLAF's **default wake path** (`VelocityMethod(1)`), so it is where the real speedup is. **Do the particle tree first.** Do the segment tree (`grow_tree_segment`, `ui_tree_segment`) only if time remains; otherwise write down precisely what it needs.

**Modes:**
- **strict:** the Fortran tree, bit for bit on the CPU: same tree structure, traversal order, multipole sums and acceptance criterion. Evaluated in a GPU-legal form (flat node arrays, explicit stack instead of recursion and pointers). It must be bitwise against the OPENMP=OFF Fortran.
  - Check first whether the Fortran's tree build depends on thread count or OpenMP task order. If it does, say so and define strict against the OPENMP=OFF build.
- **fast:** free to use a GPU-friendly tree, e.g. Morton/space-filling-curve sorted particles, a built-on-device or bottom-up tree, or a different leaf size. It is class **R2**:
  - its error against the **direct sum** must be ≤ 1.1 × the Fortran tree's error against the direct sum on the same inputs (criteria R2);
  - it must stay deterministic run to run and across thread counts (R7).
  - Pick the design so the GPU wins at 10⁵–10⁶ particles. Measure, don't guess.
- Reuse unit 1's interaction bodies for near-field direct sums. Don't duplicate the physics.
- Device buffers (tree arrays, sorted particles) live in unit 2's handle and follow its sizing rule. No per-call allocation (P6).

## Gates (record in `tests/treecode/unit.md` with PASS/FAIL and evidence)
- **L0:**
  - claim 1 (strict vs Fortran OPENMP=OFF, bitwise);
  - R2 (fast-tree error vs direct ≤ 1.1 × Fortran-tree error vs direct, every regularization, branch factor and `DistanceDirect` value used, and edge cases such as coincident particles, the degenerate-domain threshold `1.0e-5`, one particle, an empty tree);
  - R7.
- **System level:**
  - r-test OLAF cases on the fork (strict `.outb` identical to OFF);
  - **R8** on the BAR tree-wake case `opt/build/iv-cases/bar-r8-tree`, made by unit 2, against base₁ from `opt/build/iv-r8-1860700/bar-r8-tree.fomp1` and the FMA envelope from `opt/build/iv-r8ffma-1860993/bar-r8-tree.ffma`. Use full-precision text output, as unit 2 did.
- **P2:** like for like, i.e. a fork built `OPENMP=ON` plus Kokkos-OpenMP vs the pin with OpenMP, 8 threads, on the tree-wake BAR case.
- **P1/P5:** a scaling sweep of particles N = 10³…10⁶ for tree build plus evaluation:
  - Fortran+OpenMP (8 threads);
  - Kokkos-OpenMP;
  - Kokkos-CUDA, with the crossover point and the per-call split (build, evaluation, H↔D).
- **A1/A2/A5/A7** on the fork build. **P6.**

## Building and running
- Toolchain: `source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28`. Compile everything on the login node.
- Kokkos 5.3.0 installs are in `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/`: `-omp`, `-cuda-h100` (sm_90) and `-cuda-a100` (sm_80).
- **Use your own new build trees,** e.g. `opt/build/u3-*` and `opt/build/fork3-{off,omp,fomp,cuda-h100,cuda-a100}`. Do **not** rebuild or modify existing trees (`fork-*`, `iv-*`, `bsd-*`, `biotsavart-direct-*`): queued GPU jobs (1860007, 1860223, 1860779) and the running A100 job 1862373 use them.
- **GPU correctness:** `-p interactive_gpu --gres=gpu:1 -c 4`, at most 4 h. These are A100 MIG slices with no queue; use the `-cuda-a100` builds. A template is `tests/induced-velocities/gpu_a100_mig.sbatch`.
- **GPU timing** needs a full GPU:
  - Submit early to `-p h100,h200` (sm_90 build) **and** separately to `-p a100` (sm_80 build), both `--gres=gpu:1 -c 8`. Report whichever runs.
  - The h100 queue currently estimates 2026-11-16, so don't block on it. If timing jobs are still pending when you are done, report their ids and log paths.
- **CPU timing:** the `med` partition, `-c 8`, never two concurrent runs with `OMP_PROC_BIND=spread`.
- Job names: `ofk-treecode-<what>`. The PI-wide limit is 16 GPUs.
- Count `cpu_core_h` and `gpu_h` from `sacct` (`tools/cost/compute_hours.sh 2026-11-14 ofk-treecode`); don't estimate.

## Rules
- Don't edit `deps/openfast`.
- **No git commits or pushes.**
- Write only in `fork/` and `tests/treecode/`, plus `tests/induced-velocities/` or `tests/biotsavart-direct/` only where their tests need adapting; say so if you do.
- Put suggested edits to shared docs in your report instead of making them.
- **If the session is about to end** (you'll get a message), checkpoint as unit 2 did:
  - move remaining work into `sbatch` jobs with dependencies;
  - write a "Status / resume" section at the top of `tests/treecode/unit.md`;
  - then report.

## Report back (concise)
1. **Files and fork diff summary.**
2. **Gate table with evidence**, the tree design decision for fast mode (and why), and the R2 and R8 numbers.
3. **Scaling sweep table and crossover.**
4. **Slurm jobs and log paths.**
5. **Lessons for the architecture, now that all three units exist:** what generalizes to other OpenFAST modules, what was OLAF-specific, what you'd change in `architecture.md` or the skill, and what cost the most time.
6. **Suggested exact edits to shared docs.**
7. **Wall-clock time spent.**
