#!/usr/bin/env python3
"""Per-unit and per-phase totals of tools/cost/ledger.csv, as Markdown tables.

  tools/cost/summary.py [--ledger tools/cost/ledger.csv]

Money is USD at API list prices (API-equivalent, see README.md). Cache reads are split out because they
dominate: their cost is (requests) x (context size), so it grows with session length.
"""
import argparse, collections, csv, json, os

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--ledger", default=os.path.join(here, "ledger.csv"))
ap.add_argument("--prices", default=os.path.join(here, "prices.json"))
a = ap.parse_args()
rows = list(csv.DictReader(open(a.ledger)))
price = json.load(open(a.prices))["models"]

def num(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

def cache_read_usd(r):
    m = (r.get("models") or "claude-opus-5-5").split(";")[0]
    return num(r["cache_read_tok"]) * price.get(m, price["claude-opus-5-5"])["cache_read"] / 1e6

keys = ["usd_api_equiv", "active_h", "requests", "output_tok", "cache_read_tok", "cpu_core_h", "gpu_h"]
def agg(group):
    t = collections.OrderedDict()
    for r in rows:
        g = group(r); s = t.setdefault(g, collections.Counter())
        for k in keys: s[k] += num(r[k])
        s["cr_usd"] += cache_read_usd(r); s["sessions"] += 1
    return t

def table(title, t):
    print(f"\n### {title}\n")
    print("| | sessions | USD (API-eq.) | of which cache reads | active h | requests | avg context read per request | CPU core-h | GPU-h |")
    print("|---|---|---|---|---|---|---|---|---|")
    tot = collections.Counter()
    for g, s in t.items():
        tot.update(s)
        print(f"| {g} | {s['sessions']:.0f} | {s['usd_api_equiv']:.2f} | {s['cr_usd']:.2f} ({100*s['cr_usd']/max(s['usd_api_equiv'],1e-9):.0f}%) | "
              f"{s['active_h']:.2f} | {s['requests']:.0f} | {s['cache_read_tok']/max(s['requests'],1)/1e3:.0f}k | {s['cpu_core_h']:.1f} | {s['gpu_h']:.2f} |")
    print(f"| **total** | {tot['sessions']:.0f} | **{tot['usd_api_equiv']:.2f}** | {tot['cr_usd']:.2f} ({100*tot['cr_usd']/max(tot['usd_api_equiv'],1e-9):.0f}%) | "
          f"{tot['active_h']:.2f} | {tot['requests']:.0f} | {tot['cache_read_tok']/max(tot['requests'],1)/1e3:.0f}k | {tot['cpu_core_h']:.1f} | {tot['gpu_h']:.2f} |")

table("By unit", agg(lambda r: r["unit"]))
table("By unit and phase", agg(lambda r: f"{r['unit']} / {r['phase']}"))
table("By session (one row per agent)", agg(lambda r: f"{r['date']} {r['unit']} / {r['phase']} ({r['session']})"))
