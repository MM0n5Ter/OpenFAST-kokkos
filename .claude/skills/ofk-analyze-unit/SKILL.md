---
name: ofk-analyze-unit
description: Analyze an OpenFAST Fortran module or unit before porting it to Kokkos - map the per-step call graph, decide CPU vs GPU placement per routine, list fusion candidates with their legality class, and record porting traps. Use when starting any new unit/module (e.g. "analyze InflowWind FlowField", "which parts of BEMT should go on the GPU", "where can we fuse kernels in X").
---

# Analyze a unit for porting (placement + fusion)

Output: `notes/<module>-placement.md`, in the shape of `notes/olaf-placement.md` (the worked example). Read that note first and copy its section structure.
Criteria the analysis must apply: `notes/criteria.md` gates **L** (placement) and **F** (fusion legality), and equivalence classes **B/R1/R2/S**.

Ground rules:
- Read the pinned source in `deps/openfast` only. **Never edit it.**
- Cite every claim as `file:line` on the pin.
- Static shape is not enough: placement depends on **how often** a routine runs, and only the call graph shows that.

## Steps

1. **Inventory.** Run `scripts/kernel-inventory.py <files>`, and `--md` for the table.
   - Note routines with `wake`-scale loops (the column counts loops whose bound scales with the problem: markers, segments, points, grid), `iHead` counters, recursion or pointer trees, and OpenMP directives.
   - For a non-OLAF module, first extend `RE_WAKE` in the script with that module's large-dimension names, e.g. the number of wind points or grid points, and say so in the note.
2. **Step anatomy.**
   - Start from the module's framework entry points (`scripts/api-surface.sh`; `<Mod>_UpdateStates`, `_CalcOutput`, `_CalcConstrStateResidual`).
   - Draw the call tree down to the large loops.
   - For every large-loop path, record **calls per time step** (and how integrator options such as RK4 or UA multiply them) and what each call rebuilds from state.
   - Read the input-file defaults (`*_IO.f90` `ReadVarWDefault`) and any hard-wired options. In OLAF, `VelocityMethod(2)` is hard-wired, which made the lifting-line path expensive.
3. **Sizes.** Write formulas for each large dimension in terms of input parameters, plus one realistic example (not an r-test size; those are tiny).
4. **Placement table.** For each routine, decide DEVICE (month 1), DEVICE (phase 2, needs residency), or HOST, applying L1–L4. Where your decision overrides the script's hint, say why (call frequency, callee must stay on host, …).
5. **Boundary.** For month 1, list what crosses host↔device per call and per step, in bytes as a function of N. Name the state owner, and list the Fortran entry points that read that state (`scripts/callers.sh`, grep for the state's components outside the unit).
6. **Kernel shapes.** Classify each device kernel as square, skinny, tree or elementwise, then give the Kokkos pattern and whether the per-target sum order matches the Fortran (class B or R1).
7. **Fusion candidates.** Table: ID, what is fused, where (file:line), gain, legality (B/R/A per criteria F1–F5), phase. Include the "do not fuse" list with reasons, especially across framework API boundaries.
8. **Audit list.** Every trap a translator could get wrong, each with the action and the L0 case that checks it. Look for:
   - default-kind literals in `ReKi` parameters. Their value depends on compile flags: OpenFAST's double build uses `-fdefault-real-8`, so `1e-4` is a double there. Run `scripts/literal-kinds.py <files>`, and mark each hit "pass from Fortran";
   - `!$OMP SIMD reduction` and other directives that reorder the Fortran's own sums when OpenMP is on. They decide which oracle build is class B;
   - OpenMP directives that are commented out by a typo (`!OMP` instead of `!$OMP`), so a routine assumed parallel is serial: `grep -n '^\s*!OMP' <files>`;
   - `epsilon`/`tiny`-based constants;
   - 1-based indices stored in integer arrays;
   - accumulate-vs-overwrite ("side effects");
   - "3 × n++" over-sized arrays and sentinel values (`-999`, `-9999`, `1000`);
   - guard style per branch;
   - `merge` (both arms evaluated);
   - Fortran pointers aliasing caller arrays;
   - comments that say an r-test is sensitive.
9. **Hand off.**
   - Add the note to `notes/README.md`.
   - Open or append the unit card `tests/<unit>/unit.md` (template in `notes/harness.md` §0), with the placement decisions, the declared equivalence class per kernel and the audit list.
10. **Record cost (money and time).** Before ending the session, append the ledger row (criteria.md E2, E5, E6; method in `tools/cost/README.md`):
    ```bash
    tools/cost/session_cost.py --session <id-prefix> --ledger tools/cost/ledger.csv --unit <unit> --phase analysis --note "<one line>"
    ```
    - Ask the user for `--human-min` if they haven't given it; leave it empty rather than guess.
    - Copy the session's `$` (API-equivalent), `active_h` and `human_prompts` into the unit card.

## Done when
- Every wake-scale routine has a decision and a reason.
- Every fusion has a legality class.
- Every audit item has an L0 case name, even if the test isn't written yet.
- The note states the one or two findings that change the plan, at the top.
