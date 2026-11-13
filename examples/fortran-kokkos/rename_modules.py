#!/usr/bin/env python3
"""Rename Fortran modules in a frozen source file, to build an in-process oracle.

    module X / end module X / use X [, only: ...]  ->  ... X_ref ...

for each X given (case-insensitive, as Fortran is). Everything else is copied
unchanged, so the oracle is the frozen code. Module procedures then get mangled
as __x_ref_MOD_* and do not clash with the port's __x_MOD_*.

Limits: a file that also defines bind(C) names or non-module procedures would
still clash; check with `nm` (see notes/harness.md).

Usage: rename_modules.py IN OUT MODULE [MODULE ...]
"""
import re, sys

src, dst, *mods = sys.argv[1:]
if not mods:
    sys.exit(__doc__)
names = "|".join(re.escape(m) for m in mods)
pat = re.compile(
    rf"^([ \t]*(?:end[ \t]*)?module[ \t]+|[ \t]*use[ \t]*(?:::[ \t]*)?)({names})\b(?![ \t]*procedure)",
    re.IGNORECASE | re.MULTILINE)
with open(src) as f:
    text = f.read()
out, n = pat.subn(lambda m: f"{m.group(1)}{m.group(2)}_ref", text)
if n == 0:
    sys.exit(f"rename_modules: no module statement for {mods} in {src}")
with open(dst, "w") as f:
    f.write(out)
