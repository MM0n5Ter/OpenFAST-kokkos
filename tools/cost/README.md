# Cost ledger: what each unit costs in money and time

Every agent session that works on this repository ends with one ledger row (criteria.md **E2, E5, E6**). The rows give the per-unit cost. Summed over the month, they give the rates ($ per kLOC, hours per unit) that roadmap.md uses to project the cost of porting the rest of OpenFAST.

| What | Measured by | Column(s) |
|---|---|---|
| **Money: model tokens** | `session_cost.py`: usage of every model request in the Claude Code transcript (main session plus subagents), deduplicated by `requestId`, priced with `prices.json` | `input_tok`, `output_tok`, `cache_write_tok`, `cache_read_tok`, `usd_api_equiv` |
| **Money: compute** | `compute_hours.sh`: Slurm accounting for jobs named `ofk-*` | `cpu_core_h`, `gpu_h` |
| **Time: agent** | `session_cost.py`. `wall_h` runs from the first to the last event. `active_h` sums the gaps between events that are ≤ 5 min (`--idle`), so overnight idle doesn't count. | `wall_h`, `active_h` |
| **Time: people** | Entered by hand. Only the person knows how long they spent reading, reviewing or unblocking. | `human_minutes` |
| **Interventions** | `session_cost.py` counts typed and queued human prompts. The unit card adds one line per intervention saying *why* (E3). | `human_prompts` |
| Size of the request | `requests`, `tool_calls` | for $ per request and for spotting runaway loops |

## Use

```bash
# at the end of a session (run from the directory Claude Code was started in, or pass --dir)
tools/cost/session_cost.py --session <id-prefix> --ledger tools/cost/ledger.csv \
    --unit fvw-induction --phase port --note "ui_seg kernel + L0" --human-min 20 \
    --gpu-h "$(tools/cost/compute_hours.sh 2026-11-14 | awk -F'gpu_h=' '/TOTAL/{print $2}')"

tools/cost/session_cost.py --since 2026-11-01            # overview of all sessions
tools/cost/compute_hours.sh 2026-11-01                   # GPU and CPU hours of ofk-* jobs
```

- The session id is the transcript file name in `~/.claude/projects/<folder>/`. The folder is named after the directory Claude Code was started in. Inside a session, it is the directory under the session's scratchpad path.
- Name every Slurm job `ofk-<unit>-<what>` (`sbatch -J ofk-fvw-l0-h100 …`) so `compute_hours.sh` can find it.

## Rules
- **Report `$` as "USD (API-equivalent)".** The prices in `prices.json` are Anthropic list prices. A subscription plan bills differently. The equivalent is what makes runs and units comparable, and what an outside reader can reproduce. Check the prices against the official pricing page before a report, and bump `as_of`.
- Cache reads dominate the token count of long agent sessions (often over 95% of tokens) but cost 0.05× input on Opus 5.5. Keep the four token kinds separate: a cost change usually comes from cache writes or output.
- **One unit = one row per session phase** (`analysis`, `port`, `debug`, `review`, `infra`). Infra work that serves every unit, like this tooling, uses `--unit infra`.
- A session that serves two units gets two rows with `--note` saying the split. Don't split tokens by guesswork.
- The ledger is append-only and lives in git: `ledger.csv`.

## Limits (what this does not see)
- Sessions run from another config directory, another machine or the cloud: run the script where their transcripts are and append the row there, or pass `--dir`.
- Human time is self-reported.
- GPU-hours count allocation, not use: an idle allocated GPU is still billed to the group's limit.

## Splitting one session into several rows
A session often covers several units, and subagents run inside it. Split it so that no token is counted twice:
- `--agents sub --agent <id-prefix>` takes one subagent's transcript only (`<session>/subagents/agent-<id>.jsonl`). The row's `session` column is `<session8>/<agent>`. Use this for a unit delegated to a subagent.
- `--agents main` takes the parent transcript only.
- `--after T` and `--before T` (ISO, e.g. `2026-11-14T15:23:05Z`) keep only events in that window. Split at the timestamp of the human prompt that started the new phase.

Check: the parts must add up to the plain `--session` total (2026-11-14: $8.66 + $2.74 + $8.23 = $19.63).
