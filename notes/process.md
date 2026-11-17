# The porting process: plan, translate, verify, record

How one unit goes from "this is slow" to a pushed, gated port. Written 2026-11-17 from units 1–4a (`notes/runs.md` has their numbers). Rules that already live elsewhere are linked, not repeated: gates in [criteria.md](criteria.md), design in [architecture.md](architecture.md), decisions in [decisions.md](decisions.md), test levels in [testing.md](testing.md) and [harness.md](harness.md).

## 0. Roles

| Who | Does | Does not |
|---|---|---|
| **User** | Picks scope and policy (e.g. FMA default, unit boundaries, what to push). Answers decision points. | — |
| **Orchestrator** (the main session) | Plans the unit, writes the brief, launches one agent per phase, verifies each hand-back against the logs, writes the ledger rows, commits, pushes when asked, reports. | Port code itself; decide policy the user owns. |
| **Phase agent** (fresh context, `Agent` tool) | One phase of one unit (§2–§4), within the brief's scope and budget. Writes its results into the repo (card, notes, tests). | Commit, push, write ledger rows, change decided policy, edit `deps/openfast`. |

One unit is ported at a time. Each agent learns from the previous one only through the repo: the unit cards, the `ofk-port-unit` skill's lessons, and the notes.

## 1. Plan

1. **Measure before choosing.** Profile the current code on the real cases (BAR r8/perf; `scripts/wsamp/`, gprof). The unit is whatever dominates wall time, not whatever looks parallel. Unit 4's analysis overturned the plan: the circulation solve was 58% of the run, and wake residency, the expected next unit, was 0.25%.
2. **Analysis phase** (`ofk-analyze-unit` skill), output `notes/<module>-placement.md` or a unit note:
   - the call graph with calls per step;
   - placement per routine (L1–L4);
   - fusion candidates with their legality (F1–F5);
   - an audit list where every trap gets an L0 case name.
3. **Unit card draft** `tests/<unit>/unit.md`:
   - the frozen contract (signatures from the pin);
   - the equivalence class per output (B/R1/R2);
   - the audit list;
   - the gates in scope.
4. **Decision points go to the user before the port**, each with a recommendation. Examples: where the unit boundary stops (decision #14: modules we do not own are library calls), and which variant ships if an experiment is involved.
5. **Experiments get a decision rule fixed in advance, and the rule must include every limit already in force.** 4a's rule left out the per-case claim-3 limit, so the rule's winner (F1) was overturned (decision #15). Judge a speed gain end to end: a kernel-only gain does not justify losing exactness.

## 2. Brief (orchestrator → phase agent)

Template: [agent-briefs/TEMPLATE-port.md](agent-briefs/TEMPLATE-port.md). Every past brief is kept verbatim in [agent-briefs/](agent-briefs/). A brief states:
- **goal and phase**: one of analyze / port + L0 / system gates / record;
- **files to read first**, by path. Do not paste their contents;
- **what is already decided** (decisions.md numbers) and what is ruled out;
- **scope limits**: files it may change, build trees it must not touch, jobs it may submit (names `ofk-<unit>-<what>`), the GPU limit;
- **the decision rule** for any experiment, complete (§1.5);
- **budget**: a request and USD cap. At the cap the agent writes its state into the card's "Resume" section and hands back;
- **hand-back format**: gate table with evidence paths, open issues as explicit questions, the list of changed files, the job IDs. No commits and no ledger rows.

## 3. Translate (port phase)

Follow `ofk-port-unit` (`.claude/skills/ofk-port-unit/SKILL.md`). The order matters:
1. **L0 oracle first.** Build the pin's own files renamed (`X` → `X_ref`) into `oracle_<unit>`. Run ref vs ref **per case**: the pin itself reads uninitialised memory in some corners (4a). Every audit item gets a fixed case, then add randomized cases with a logged seed.
2. **Kernel.** One body in `<unit>_body.hpp`, instantiated per mode (`_strict.cpp`, `_fast.cpp`, plus any experiment variant). Strict reproduces the Fortran's order with no FMA. Pass flag-dependent constants from the Fortran.
3. **Shim.** An overlay copy of the Fortran file with `#ifdef OPENFAST_KOKKOS`. The caller does not change (A2). State goes in a per-instance handle released by the module's `End`.
4. **Mutation check.** Break each guard, order and index on a **copy** of `src/`. L0 or A7 must catch every mutation. 4a: 8/8; unit 3: 12/12.
5. **Accept empty extents.** The r-test set contains degenerate cases (e.g. `ad_Sphere_OLAF` has no lifting line).

## 4. Verify (gates)

Run cheapest first, and stop at the first failure:

| Step | Where | Gates | Typical time |
|---|---|---|---|
| 1. L0 | login node | claim 1 (strict = pin, bitwise, iteration counts included), A7, P6 at L0, R7 at L0 (threads, processes) | seconds to minutes |
| 2. Whole runs, CPU | `med` (`sbatch`, 8 cores) | system claim 1 (strict = base₁ from line 8), R8 (drift ≤ 2× the FMA envelope), R7, P2 (vs `openfast-pin-fomp`, same node), P4, P6 (memtool) | 10–60 min |
| 3. Regression | `med` | R6: `ctest -R 'aerodyn_utest\|OLAF\|5MW_MRSemi'`, in strict and in default mode | about 30 min |
| 4. GPU correctness | `interactive_gpu` (A100 MIG, no queue) | claim 2 (strict GPU vs strict CPU), CUDA build and link, R7 on the GPU | minutes |
| 5. GPU performance | `a100` / `h100` (queues of hours to days) | P1/P3 on full GPUs | submit early, never block on it |
| 6. API | login node | A1 (`nm` symbols), A2 (preprocess without `OPENFAST_KOKKOS` = pin), A5 (`.mod`/objects identical), single-precision build | minutes |

The equivalence chain these gates check is in architecture.md §4. **State bitwise claims in two separate parts** (the template is the "Bitwise and reproducibility claims" section of `tests/circulation-solve/unit.md`):
- **A. Equality with the Fortran:** for each mode, say which outputs are identical to which Fortran build. Say where the default mode is *not* identical, and why (e.g. other units' fast mode).
- **B. Reproducibility:** the same build gives the same bits across threads and processes. B says nothing about A.

Any case that exceeds a limit is **kept as a fixed, reproducible L0 case** (seed or generator state) and recorded in the card. Do not loosen the limit, and call an explanation plausible until it is shown.

## 5. Record

The phase agent writes into the repo:
- **unit card**: the status, the bitwise claims (A/B), a gate table with evidence (log path, job ID), Slurm jobs, lessons, and "Resume";
- **`ofk-port-unit` skill**: lessons that apply to later units (keep them short; the card keeps the detail);
- **notes**: corrections of earlier notes, marked "(corrected <date>)".

The orchestrator:
- **verifies the hand-back**: reads the key logs, re-runs L0 on the final source itself, and checks that the claims match the evidence. 4a's L0 logs were 2 minutes older than its last source edit, so they were re-run;
- writes the **ledger row** (`tools/cost/session_cost.py --agents sub --agent <id>`, CPU/GPU hours from `tools/cost/compute_hours.sh`), plus a row for its own orchestration;
- cancels jobs that nothing depends on any more (4a's optional profiling job, 16 of its 35.7 core-hours);
- **commits** with a message that names the gates, and pushes only when the user asks;
- reports to the user, with the decisions that are theirs to make.

## 6. Cost control

Measured in `notes/runs.md` §2: about half of all spend is **cache reads**, which cost (requests) × (context size). A long session therefore costs roughly quadratically in its length: 4a's port agent averaged 423k tokens of context per request over 254 requests. The rules, with the expected effect:
1. **One phase per agent** (§2), so no agent carries the whole unit. Expected: average context about 150k instead of 300–420k, so cache-read cost about halves.
2. **A shared gate runner** instead of per-unit job scripts. Each unit wrote 5–7 sbatch files and helpers (4a: about 1,500 lines of tests and scripts). Use parameterized templates in `tests/common/` (to build) that print one verdict line per gate and write details to log files.
3. **Quiet tools.** Bash output was 150–515k characters per agent. Scripts print verdicts; agents `grep` logs instead of `cat`-ing them, and never read whole Fortran files (cite `file:line`, read ranges).
4. **Experiments are scoped up front** (§1.5), with their own budget line in the brief.
5. **A budget in every brief**, with a checkpoint hand-back at the cap.
6. **Model by phase (trial pending):** Opus for analysis and translation; try Sonnet (half the price) for the system-gate and record phases. Compare its verification quality on one unit before adopting it.
7. **Orchestrator discipline:** verify with targeted `grep` on logs, use one background wait per job group instead of polling, and do not re-derive what the card already says.

Target: **about $10–15 (API-equivalent) per mid-size port unit**, against $16–40 so far.

## 7. Where things live

| What | Where |
|---|---|
| Gates and classes | `notes/criteria.md` |
| Design rules (modes, FMA, handles, equivalence chain) | `notes/architecture.md` |
| Decisions with dates and owners | `notes/decisions.md` |
| Skills (agent procedures) | `.claude/skills/ofk-analyze-unit/`, `.claude/skills/ofk-port-unit/` |
| Briefs given to agents (verbatim) and the template | `notes/agent-briefs/` |
| Per-unit evidence | `tests/<unit>/unit.md` |
| Run statistics (cost, time, compute, outcomes) | `notes/runs.md`, generated parts from `tools/cost/summary.py` |
| Raw cost ledger and method | `tools/cost/ledger.csv`, `tools/cost/README.md` |
| Build configs | `tools/openfast-build.sh`, `src/cmake/` |
