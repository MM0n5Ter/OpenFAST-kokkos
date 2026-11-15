#!/usr/bin/env python3
"""Money and time spent by agent sessions, from Claude Code transcripts.

Reads ~/.claude/projects/<project>/<session>.jsonl (plus <session>/subagents/*.jsonl when present)
and reports, per session:
  money  tokens by kind (input, output, cache write 5m/1h, cache read) per model, and USD at the
         API list prices in prices.json ("USD (API-equivalent)"; a subscription bills differently)
  time   wall clock (first->last event), active agent time (sum of gaps between consecutive events
         that are <= --idle seconds), number of model requests and tool calls
  people human prompts (typed or queued messages; task notifications and tool results excluded)

Usage:
  tools/cost/session_cost.py                              # all sessions of this project
  tools/cost/session_cost.py --session 199ad883            # one session (prefix match)
  tools/cost/session_cost.py --since 2026-11-14 --json out.json
  tools/cost/session_cost.py --session 199ad883 --ledger tools/cost/ledger.csv \
      --unit olaf --phase analysis --note "placement + cost tooling"

The ledger row leaves human_minutes and gpu_hours for you to fill in (or pass --human-min/--gpu-h;
tools/cost/compute_hours.sh gives the GPU and CPU hours from Slurm).
"""
import argparse, csv, datetime as dt, glob, json, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PRICES = os.path.join(HERE, "prices.json")
LEDGER_FIELDS = ["date", "unit", "phase", "session", "models", "requests", "tool_calls", "human_prompts",
                 "input_tok", "output_tok", "cache_write_tok", "cache_read_tok", "usd_api_equiv",
                 "wall_h", "active_h", "human_minutes", "cpu_core_h", "gpu_h", "note"]


def project_dir():
    slug = os.getcwd().replace("/", "-")
    base = os.path.expanduser("~/.claude/projects")
    # The transcript folder is named after the directory Claude Code was started in, which can be a parent.
    cands = [d for d in glob.glob(os.path.join(base, "*")) if slug.startswith(os.path.basename(d))]
    return max(cands, key=len) if cands else None


def ts(s):
    return dt.datetime.fromisoformat(s.replace("Z", "+00:00"))


def load(path):
    with open(path, errors="replace") as f:
        for line in f:
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def is_human_text(text):
    t = text.lstrip()
    return bool(t) and not t.startswith(("<task-notification", "<system-reminder", "Base directory for this skill"))


def analyze(files, idle, prices, after=None, before=None):
    reqs = {}                       # requestId -> (model, usage, speed); last streamed line wins
    times, humans, tool_calls = [], set(), 0
    for p in files:
        for o in load(p):
            if "timestamp" in o:
                try:
                    t0 = ts(o["timestamp"])
                except ValueError:
                    t0 = None
                if t0 and ((after and t0 < after) or (before and t0 >= before)):
                    continue                # outside the --after/--before window
                if t0:
                    times.append(t0)
            elif after or before:
                continue
            typ = o.get("type")
            if typ == "assistant":
                m = o.get("message", {})
                key = o.get("requestId") or m.get("id")
                if key and m.get("usage"):
                    reqs[key] = (m.get("model", "?"), m["usage"])
                tool_calls += sum(1 for b in m.get("content", []) or [] if isinstance(b, dict) and b.get("type") == "tool_use")
            elif typ == "user" and not o.get("isSidechain"):
                c = o.get("message", {}).get("content")
                if isinstance(c, str) and is_human_text(c):
                    humans.add(c.strip()[:200])
            elif typ == "queue-operation" and o.get("operation") == "enqueue" and is_human_text(o.get("content", "")):
                humans.add(o["content"].strip()[:200])
    tok = defaultdict(lambda: defaultdict(int))
    usd, unknown = 0.0, set()
    for model, u in reqs.values():
        if model == "<synthetic>":   # locally generated placeholder messages, no billed tokens
            continue
        cc = u.get("cache_creation") or {}
        w1h = cc.get("ephemeral_1h_input_tokens", 0)
        w5m = cc.get("ephemeral_5m_input_tokens", u.get("cache_creation_input_tokens", 0) - w1h)
        t = tok[model]
        t["input"] += u.get("input_tokens", 0)
        t["output"] += u.get("output_tokens", 0)
        t["cache_write_5m"] += w5m
        t["cache_write_1h"] += w1h
        t["cache_read"] += u.get("cache_read_input_tokens", 0)
        pr = prices["models"].get(model)
        if not pr:
            unknown.add(model)
            continue
        mult = prices.get("fast_mode_multiplier", {}).get(model, 1.0) if u.get("speed") == "fast" else 1.0
        usd += mult * sum(pr[k] * v for k, v in (("input", u.get("input_tokens", 0)), ("output", u.get("output_tokens", 0)),
                                                  ("cache_write_5m", w5m), ("cache_write_1h", w1h),
                                                  ("cache_read", u.get("cache_read_input_tokens", 0))) if v) / 1e6
    times.sort()
    wall = (times[-1] - times[0]).total_seconds() if len(times) > 1 else 0.0
    active = sum(g for g in ((b - a).total_seconds() for a, b in zip(times, times[1:])) if g <= idle)
    return dict(start=times[0].isoformat() if times else None, requests=len(reqs), tool_calls=tool_calls,
                human_prompts=len(humans), tokens={m: dict(v) for m, v in tok.items()}, usd_api_equiv=round(usd, 4),
                wall_h=round(wall / 3600, 3), active_h=round(active / 3600, 3), unpriced_models=sorted(unknown))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", help="transcript folder (default: derived from the current directory)")
    ap.add_argument("--session", action="append", help="session id or prefix (repeatable)")
    ap.add_argument("--since", help="only sessions with events on/after this date (YYYY-MM-DD)")
    ap.add_argument("--after", help="count only events at/after this ISO time (e.g. 2026-11-14T15:20:00Z); "
                    "use it to split one session into several ledger rows without double counting")
    ap.add_argument("--before", help="count only events before this ISO time")
    ap.add_argument("--agents", choices=["all", "main", "sub"], default="all",
                    help="all = main transcript + subagents (default); main = main only; sub = subagents only")
    ap.add_argument("--agent", help="with --agents sub: only subagent transcripts whose id starts with this prefix")
    ap.add_argument("--idle", type=float, default=300, help="gaps longer than this (s) are not active time")
    ap.add_argument("--json", help="write the per-session results here")
    ap.add_argument("--ledger", help="append one row per session to this CSV")
    ap.add_argument("--unit", default="")
    ap.add_argument("--phase", default="", help="analysis | port | debug | review | infra")
    ap.add_argument("--note", default="")
    ap.add_argument("--human-min", default="", help="human attention in minutes (only a person knows this)")
    ap.add_argument("--gpu-h", default="")
    ap.add_argument("--cpu-h", default="")
    a = ap.parse_args()
    d = a.dir or project_dir()
    if not d or not os.path.isdir(d):
        sys.exit("transcript folder not found; pass --dir ~/.claude/projects/<project>")
    prices = json.load(open(PRICES))
    results = {}
    for main_file in sorted(glob.glob(os.path.join(d, "*.jsonl"))):
        sid = os.path.basename(main_file)[:-6]
        if a.session and not any(sid.startswith(s) for s in a.session):
            continue
        subs = sorted(glob.glob(os.path.join(d, sid, "subagents", "*.jsonl")))
        if a.agent:
            subs = [s for s in subs if os.path.basename(s)[len("agent-"):].startswith(a.agent)]
        files = {"all": [main_file] + subs, "main": [main_file], "sub": subs}[a.agents]
        if not files:
            continue
        after = ts(a.after) if a.after else None
        before = ts(a.before) if a.before else None
        r = analyze(files, a.idle, prices, after, before)
        if not r["requests"]:
            continue
        if a.since and (not r["start"] or r["start"][:10] < a.since):
            continue
        r["subagent_files"] = len(files) - (a.agents != "sub")
        key = sid if a.agents != "sub" else f"{sid[:8]}/{a.agent or 'sub'}"
        results[key] = r
    for sid, r in results.items():
        tk = defaultdict(int)
        for v in r["tokens"].values():
            for k, n in v.items():
                tk[k] += n
        print(f"{sid[:8]}  start {r['start'][:16] if r['start'] else '?'}  ${r['usd_api_equiv']:.2f} (API-equiv)  "
              f"wall {r['wall_h']:.2f} h  active {r['active_h']:.2f} h  requests {r['requests']}  "
              f"tools {r['tool_calls']}  human prompts {r['human_prompts']}")
        print(f"          tokens: in {tk['input']:,}  out {tk['output']:,}  cache-write {tk['cache_write_5m'] + tk['cache_write_1h']:,}  "
              f"cache-read {tk['cache_read']:,}  models {', '.join(r['tokens'])}"
              + (f"  UNPRICED: {r['unpriced_models']}" if r["unpriced_models"] else ""))
    if a.json:
        json.dump(results, open(a.json, "w"), indent=1)
    if a.ledger:
        new = not os.path.exists(a.ledger)
        with open(a.ledger, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=LEDGER_FIELDS)
            if new:
                w.writeheader()
            for sid, r in results.items():
                tk = defaultdict(int)
                for v in r["tokens"].values():
                    for k, n in v.items():
                        tk[k] += n
                w.writerow(dict(date=(r["start"] or "")[:10], unit=a.unit, phase=a.phase, session=sid if "/" in sid else sid[:8],
                                models="+".join(r["tokens"]), requests=r["requests"], tool_calls=r["tool_calls"],
                                human_prompts=r["human_prompts"], input_tok=tk["input"], output_tok=tk["output"],
                                cache_write_tok=tk["cache_write_5m"] + tk["cache_write_1h"], cache_read_tok=tk["cache_read"],
                                usd_api_equiv=r["usd_api_equiv"], wall_h=r["wall_h"], active_h=r["active_h"],
                                human_minutes=a.human_min, cpu_core_h=a.cpu_h, gpu_h=a.gpu_h, note=a.note))
        print(f"appended {len(results)} row(s) to {a.ledger}")


if __name__ == "__main__":
    main()
