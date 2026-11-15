#!/usr/bin/env python3
"""Find real literals without a kind suffix (default kind = single precision in gfortran).

Why it matters for a port: the value of a default-kind literal depends on the compile flags.
OpenFAST's DOUBLE_PRECISION build adds `-fdefault-real-8 -fdefault-double-8`
(cmake/OpenfastFortranOptions.cmake:138-140), so there `MINNORM=1e-4` is exactly the double 1e-4.
Without that flag it would be float(1e-4), which is what the "value-as-stored" column shows.
Each hit goes into the unit card's audit list. The C++ side takes the value from the Fortran shim
(or as ReKi) instead of hard-coding either form. Corrected 2026-11-14 by the biotsavart-direct
port: hard-coding double(1e-4f) failed L0.

Usage: scripts/literal-kinds.py FILE...   (default: OLAF files)
Output: file:line  literal  value-if-compiled-WITHOUT-fdefault-real-8  context
"""
import os, re, struct, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "deps", "openfast", "modules", "aerodyn", "src")
DEFAULT = [os.path.join(SRC, f) for f in ("FVW.f90", "FVW_Subs.f90", "FVW_BiotSavart.f90",
                                          "FVW_VortexTools.f90", "FVW_Wings.f90")]

# A real literal: digits with a '.' and/or an exponent, NOT followed by _kind and not d-exponent.
RE_LIT = re.compile(r"(?<![\w.])((?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?)(?![\w.]|_)")


def as_float32(s):
    v = float(s)
    return struct.unpack("f", struct.pack("f", v))[0]


def code_part(line):
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "!":
            return line[:i]
    return line


def strip_strings(s):
    return re.sub(r"'[^']*'|\"[^\"]*\"", "''", s)


def main(paths):
    hits = 0
    for p in paths or DEFAULT:
        with open(p, errors="replace") as f:
            for n, raw in enumerate(f, 1):
                code = strip_strings(code_part(raw))
                if re.match(r"\s*(print|write|format)\b", code, re.I):
                    continue
                for m in RE_LIT.finditer(code):
                    lit = m.group(1)
                    if not re.search(r"[.eE]", lit):
                        continue  # integer literal
                    v32 = as_float32(lit)
                    if v32 == float(lit):
                        continue  # exactly representable (0.5, 2.0, 1.5 ...): harmless
                    hits += 1
                    print(f"{os.path.basename(p)}:{n}  {lit:<10} -> {v32:.17g}  | {raw.strip()[:90]}")
    print(f"# {hits} default-kind real literals that are not exact in single precision", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
