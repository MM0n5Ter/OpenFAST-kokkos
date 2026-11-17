#!/usr/bin/env python3
"""FMA experiment of unit circulation-solve: compare per-solve iteration counts against a reference.
  iters_cmp.py <ref> <label>=<file> [...]
Files: "<solve> <iterations> [cap]" per line ('#' comments): itertrace (pin binaries) or OFK_CS_ITERS (port).
Prints solves, mean iterations, solves whose count differs, max |delta|, first differing solve, delta histogram."""
import collections, sys
def load(f):
    return [int(l.split()[1]) for l in open(f) if l.strip() and not l.startswith('#')]
ref = load(sys.argv[1])
print(f"ref {sys.argv[1]}: {len(ref)} solves, {sum(ref)} iterations, mean {sum(ref)/max(1,len(ref)):.2f}, at cap 30: {sum(1 for x in ref if x >= 30)}")
for a in sys.argv[2:]:
    lab, f = a.split('=', 1)
    try:
        v = load(f)
    except OSError as e:
        print(f"  {lab}: missing ({e})"); continue
    n = min(len(v), len(ref))
    d = [v[i] - ref[i] for i in range(n)]
    nd = sum(1 for x in d if x)
    first = next((i + 1 for i, x in enumerate(d) if x), None)
    h = collections.Counter(x for x in d if x)
    print(f"  {lab}: {len(v)} solves{' (LENGTH DIFFERS)' if len(v) != len(ref) else ''}, {sum(v)} iterations, "
          f"differing solves {nd}/{n}, max |delta| {max((abs(x) for x in d), default=0)}, first at solve {first}, "
          f"deltas {dict(sorted(h.items()))}")
