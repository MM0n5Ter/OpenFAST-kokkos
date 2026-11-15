"""Compare the channel data of OpenFAST .outb files between runs (the files themselves always
differ: the header carries a timestamp). Channels are stored quantized, so "bitwise" here means
identical at output resolution.
  python cmp_outb.py deps/openfast/reg_tests/lib <ref dir> <dir> [<dir> ...]
Each dir holds the .outb files at the same relative paths (r6.sbatch / gpu2.sbatch copy them).
Needs numpy (e.g. $OPT/venv-rtest/bin/python)."""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, sys.argv[1])
from fast_io import load_output  # noqa: E402  (OpenFAST reg_tests/lib)

ref, others = sys.argv[2], sys.argv[3:]
for f in sorted(glob.glob(os.path.join(ref, "**", "*.outb"), recursive=True)):
    rel = os.path.relpath(f, ref)
    a = load_output(f)[0]
    row = [rel]
    for o in others:
        g = os.path.join(o, rel)
        if not os.path.exists(g):
            row.append("missing")
            continue
        b = load_output(g)[0]
        if a.shape != b.shape:
            row.append("shape!")
        elif np.array_equal(a, b, equal_nan=True):
            row.append("identical")
        else:
            row.append(f"maxrel={np.nanmax(np.abs(a - b)) / (np.nanmax(np.abs(a)) or 1.0):.1e}")
    print("  ".join(row))
