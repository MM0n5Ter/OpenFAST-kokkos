#!/usr/bin/env python3
"""Directory-level dependency graph of OpenFAST from Fortran `use` statements.

Each Fortran `module X` is assigned to the directory that defines it
(modules/<name> or glue-codes/<name>). An edge A -> B means a file in A
`use`s a module defined in B. Module names defined in more than one
directory (e.g. platform-specific Sys* files) are skipped.

Usage: scripts/module-deps.py [--root deps/openfast] [--edges] [--json out.json]
"""
import argparse, json, os, re, sys
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
here = os.path.dirname(os.path.abspath(__file__))
ap.add_argument("--root", default=os.path.join(here, "..", "deps", "openfast"))
ap.add_argument("--edges", action="store_true", help="print every A -> B edge")
ap.add_argument("--json", help="write {dir: {uses: [...], used_by: [...]}} here")
a = ap.parse_args()

MOD = re.compile(r"^\s*module\s+(?!procedure\b)(\w+)", re.I | re.M)
USE = re.compile(r"^\s*use\s*(?:,\s*\w+\s*::\s*|::\s*|\s)\s*(\w+)", re.I | re.M)

files = defaultdict(list)                       # dir -> [paths]
for top in ("modules", "glue-codes"):
    base = os.path.join(a.root, top)
    for name in sorted(os.listdir(base)):
        for dp, _, fs in os.walk(os.path.join(base, name)):
            files[f"{top}/{name}"] += [os.path.join(dp, f) for f in fs
                                       if f.lower().endswith((".f90", ".f"))]

def read(p):
    with open(p, errors="replace") as f:
        return f.read()

owner = defaultdict(set)                        # module name -> {dirs}
for d, ps in files.items():
    for p in ps:
        for m in MOD.findall(read(p)):
            owner[m.lower()].add(d)
unique = {m: next(iter(ds)) for m, ds in owner.items() if len(ds) == 1}

uses = defaultdict(set)
for d, ps in files.items():
    for p in ps:
        for m in USE.findall(read(p)):
            o = unique.get(m.lower())
            if o and o != d:
                uses[d].add(o)
used_by = defaultdict(set)
for d, os_ in uses.items():
    for o in os_:
        used_by[o].add(d)

dirs = sorted(set(files) | set(used_by))
print(f"{'fan_in':>6} {'fan_out':>7}  dir: uses")
for d in sorted(dirs, key=lambda d: (-len(used_by[d]), d)):
    if not files[d]:
        continue
    deps = " ".join(x.split("/")[1] for x in sorted(uses[d]))
    print(f"{len(used_by[d]):>6} {len(uses[d]):>7}  {d}: {deps}")
if a.edges:
    for d in dirs:
        for o in sorted(uses[d]):
            print(f"{d} -> {o}")
if a.json:
    with open(a.json, "w") as f:
        json.dump({d: {"uses": sorted(uses[d]), "used_by": sorted(used_by[d])} for d in dirs}, f, indent=1)
