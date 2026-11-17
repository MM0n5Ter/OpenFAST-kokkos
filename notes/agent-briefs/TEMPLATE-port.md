# Brief template: one phase of one unit

Copy, fill every `<...>`, and delete what does not apply. Point to files by path, and do not paste their contents. Keep the brief under about 5,000 characters. Process: `notes/process.md`.

---

You are doing **<phase: analysis | port + L0 | system gates | record>** for unit **<n> `<unit>`** of the OpenFAST-kokkos port. Repo: `/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos`. The repo's `CLAUDE.md` and the `ofk-port-unit` skill apply.

**Goal.** <one or two sentences: the routine(s), the change, and what "done" means for this phase>

**Read first** (in this order; grep and read ranges rather than whole files):
- `tests/<unit>/unit.md`: the contract, classes, audit list, and the "Resume" section if a previous phase left one;
- `<notes/... the analysis section that chose this unit>`;
- `notes/process.md` §<3|4>, plus the skill's lessons for <relevant earlier unit>.

**Already decided** (do not reopen): <decisions.md #numbers and one line each, e.g. #12 no CPU↔GPU bitwise; #14 unit boundary; #15 FMA default>.
**Ruled out:** <approaches tried or rejected, with the reason>.

**Scope.**
- May change: <files and dirs>. Never: `deps/openfast`, other units' build trees, <trees held by queued jobs>.
- Build trees: prefix `<ovN>-*` / `<uN>-*`. At most `-j 12` in total, and chain builds with `|| exit 1`.
- Jobs: name them `ofk-<unit>-<what>`. Use `med` for timing and `interactive_gpu` for GPU correctness. Do not submit full-GPU jobs without asking <or: allowed: ...>.
- No commits and no ledger rows. The orchestrator does both.

**Gates for this phase:** <IDs from notes/criteria.md, and where each runs (process.md §4 order: L0 → med → R6 → MIG)>.

**Experiment (if any).** Variants <...>. Measure <...>. **Decision rule, fixed now:** <...>. The rule must include every limit already in force (claim-3 per case, R6, R8), and must judge speed end to end. Budget for the experiment: <N requests / $>.

**Output discipline.** Scripts print one verdict line per gate and write details to log files. Read logs with `grep`, never `cat` a whole log. Reuse the existing gate scripts in `tests/<earlier unit>/` before writing new ones.

**Budget.** About <N> requests / <$> USD (API-equivalent). If you reach it, write the state into the card's "Resume" section and hand back.

**Hand back:**
1. a gate table (ID, PASS/FAIL, evidence path or job ID);
2. the bitwise claims, split into A (equality with the Fortran: which mode, which outputs) and B (reproducibility);
3. open issues, phrased as questions for the user or orchestrator;
4. changed files, job IDs, and build trees that are safe or unsafe to reuse.
