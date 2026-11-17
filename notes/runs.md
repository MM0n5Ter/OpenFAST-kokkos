# Run statistics: units 1–4a (2026-11-14 to 2026-11-17)

What each unit delivered, what it cost, and where the cost went. The cost tables come from `tools/cost/summary.py` over `tools/cost/ledger.csv` (re-run it to refresh them). Money is **USD (API-equivalent)**: tokens at list prices (`tools/cost/prices.json`), not what a subscription bills. Process: [process.md](process.md). Briefs: [agent-briefs/](agent-briefs/).

## 1. Outcomes per unit

| Unit | What was ported | Agents (fresh contexts) | USD | Agent active h | Compute (CPU core-h / GPU-h) | Result | Evidence |
|---|---|---|---|---|---|---|---|
| OLAF analysis | placement and fusion plan for units 1–3 | main session | 8.47 | 0.55 | — | `notes/olaf-placement.md` | |
| 1 `biotsavart-direct` | `ui_seg`, `ui_part_nograd` (direct Biot–Savart) | 2 (v0 standalone, v1 in the build) | 24.32 | 1.62 | 3.7 / 0.04 | claim 1 bitwise; fast 1.07–3.95× vs Fortran+OpenMP at 10⁴–10⁵; raised the FMA decision (→ keep FMA, bound = 1.1× the Fortran's own FMA spread) | `tests/biotsavart-direct/unit.md` |
| 2 `induced-velocities` | `InducedVelocitiesAll_*`, LL induction, device-resident sources | 1 (checkpoint at user disconnect, finished in the main session) | 16.21 | 0.89 | 7.4 / 0.00 | R8 0/542; P2 1.02× (tree) / 2.06× (direct) on CPU; GPU L0 and full cases pass; P6 0 steady allocations | `tests/induced-velocities/unit.md` |
| 3 `treecode` | particle treecode, host and device tree build | 2 (port, finish) | 20.13 | 1.39 | 3.4 / 0.13 | claim 1 bitwise; claim 2 1.8e-16; R8 0/542; mutations 12/12; GPU eval about 47× at 10⁶; GPU sweeps queued | `tests/treecode/unit.md` |
| 4 analysis | wall-clock profile after units 1–3 | 1 | 10.05 | 0.52 | 3.9 / 0 | the circulation solve is 58% of the run, the wind 24%; residency deferred (0.25%) | `notes/olaf-unit4.md` |
| 4a `circulation-solve` | LL circulation iteration (F10-B), plus an FMA experiment | 1 | 30.38 | 2.47 | 35.7 / 0.06 | strict and F0 bitwise incl. iteration counts; R6 18/18; R8 0/542; mutations 8/8; whole run 2.52–2.59× vs Fortran+OpenMP; default F0 (decision #15) | `tests/circulation-solve/unit.md` |
| infra | Kokkos builds, cost tooling, fork → overlay migration, orchestration | 1 agent + main session | 47.78 | 5.39 | — | overlay layout (`notes/overlay.md`); ledger; this review | |
| **Total** | | 9 agents + main | **157.33** | 12.84 | 54.1 / 0.23 | | |

Cumulative effect on the BAR cases (CPU, 8 threads, against the Fortran with OpenMP on the same node): `bar-r8-tree` 44.0 → 17.0 s, `bar-r8-direct` 128.2 → 50.8 s. What remains is wind and tower work in AeroDyn (unit 4b).

**Restarts and rework** (counted in the cost above):
- unit 1 was done twice: a standalone v0, then v1 moved into the build;
- unit 2's agent was cut off by a user disconnect, and the main session finished it;
- unit 3 was split into a port and a finish session, around the layout migration;
- the overlay migration (13.70) replaced the fork layout after one day;
- 4a's first R6 run found an abort (`ad_Sphere_OLAF`, no lifting line) and was cancelled; its FMA experiment needed pin builds with FMA and a ptrace counter.

## 2. Cost

### 2.1 Tables (`tools/cost/summary.py`, 2026-11-17)

| Unit | sessions | USD | of which cache reads | active h | requests | avg context per request | CPU core-h | GPU-h |
|---|---|---|---|---|---|---|---|---|
| olaf (analysis) | 1 | 8.47 | 3.35 (40%) | 0.55 | 83 | 202k | 0.0 | 0.00 |
| biotsavart-direct | 2 | 24.32 | 13.43 (55%) | 1.62 | 266 | 252k | 3.7 | 0.04 |
| induced-velocities | 1 | 16.21 | 9.12 (56%) | 0.89 | 138 | 330k | 7.4 | 0.00 |
| treecode | 2 | 20.13 | 9.35 (46%) | 1.39 | 181 | 258k | 3.4 | 0.13 |
| circulation-solve | 2 | 40.43 | 27.01 (67%) | 2.99 | 367 | 368k | 39.6 | 0.06 |
| infra | 4 | 47.78 | 18.73 (39%) | 5.39 | 473 | 198k | 0.0 | 0.00 |
| **total** | 12 | **157.33** | 80.98 (51%) | 12.84 | 1508 | 269k | 54.1 | 0.23 |

Per session, the most expensive was the 4a port: 254 requests at 423k tokens of context each, so 21.51 of its 30.38 went on cache reads.

Two rows were added late, on 2026-11-17:
- the overlay migration agent (13.70);
- the main session's orchestration not covered by earlier rows (28.34). The main session costs 42.54 in total, of which 14.20 was already in earlier rows.

Before these rows, the ledger understated the total by 42.04.

### 2.2 Where the money goes
- **Cache reads: 51% overall, up to 71% per session.** Each request re-reads the whole context at 0.20 USD per million tokens. The cost is (requests) × (context size), and context grows through the session, so a session's cost grows roughly with the square of its length. 4a: twice the requests of unit 3 and 1.5× the context, so about 3× the cache-read cost.
- **What fills the context:**
  - tool output, mostly Bash: 150–515k characters per agent (4a: 515k over 243 calls);
  - files the agent wrote: 26–106k characters;
  - the agent's own reasoning and output: 110–290k tokens.
  - Whole-file reads were rare after unit 1. The bulk is command and job-log output.
- **Output tokens:** 2–6 USD per session. Not the main lever.
- **Per-unit job scaffolding:** each port wrote 5–7 sbatch files and helper scripts. 4a wrote about 1,500 lines of tests and scripts, beyond its 826 lines of C++.
- **Compute** is cheap next to tokens: 54 CPU core-hours and 0.23 GPU-hours in total. The real cost of GPUs is queue time (A100 about 7 h, H100 days), which `interactive_gpu` (MIG, no queue) avoids for correctness checks.

### 2.3 Is about 10–40 USD per unit high?
- **Per port unit:** 16–40 USD API-equivalent (mean about 25), plus about 10 for each analysis.
- **Overhead:** orchestration and infra are another 30% of the total.
- **Avoidable spend, judged from the drivers above:**
  - long single-agent sessions (4a);
  - repeated job scaffolding;
  - verbose command output;
  - rework from layout changes (one-off).
- **Levers** (process.md §6):
  - one phase per agent;
  - a shared gate runner;
  - quiet scripts;
  - experiments scoped up front;
  - a budget per brief;
  - a Sonnet trial for the mechanical phases.
- **Target:** about 10–15 USD per mid-size port unit. Unit 4b is the first to run under these rules, and its row will show whether they work.

## 3. Compute and queues

| Partition | Used for | Experience |
|---|---|---|
| login node | builds (`-j 12` total), L0 | load 20+ at times: no timing there |
| `med` (Xeon 8480+) | system gates, P2 timing, R6 | starts within minutes |
| `interactive_gpu` (A100 MIG 2g.20gb) | GPU correctness | no queue, 4 h limit; not for timing |
| `a100`, `h100`/`h200` | GPU performance | hours to days of queue (CVPR rush, fair-share 0.035); 4 GPU jobs pending since 2026-11-14–16 |

## 4. Regenerate
```bash
tools/cost/summary.py                       # §2.1 tables (by unit, by phase, by session)
tools/cost/compute_hours.sh 2026-11-01      # Slurm CPU/GPU hours of ofk-* jobs
tools/cost/session_cost.py --session <id>   # one session; --agents sub --agent <prefix> for one agent
```
