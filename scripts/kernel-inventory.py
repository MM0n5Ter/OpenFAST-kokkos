#!/usr/bin/env python3
"""Per-routine inventory of Fortran source, with a first-pass CPU/GPU placement hint.

Usage:
  scripts/kernel-inventory.py [FILES...] [--md] [--json OUT]
  (default FILES: the OLAF sources in deps/openfast/modules/aerodyn/src)

For every subroutine/function (including contained ones) it reports:
  lines, max loop-nest depth, OpenMP directives, allocate statements,
  "iHead" running counters (serial packing pattern), select-case dispatch,
  recursion, pointer use, STOP/print (non-device code), registry-typed args
  (p/x/m/z/u/y: framework types that cannot cross to the device as is),
  and the routines it calls.

The placement hint is a heuristic, not a decision. It flags what a human (or the
porting agent) must look at; notes/olaf-placement.md records the decisions.
  DEVICE      has a loop over a wake-sized range (nNW, nFW, nSeg, nCPs, nPart...), no I/O -> Kokkos kernel
  DEVICE*     like DEVICE but has a blocker to rewrite first (iHead counter, recursion, pointer tree)
  GLUE        orchestration over registry types (calls kernels, little arithmetic) -> stays host, drives device
  HOST        I/O, error handling, small serial logic, solvers on O(nSpan) data
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SRC = os.path.join(HERE, '..', 'deps', 'openfast', 'modules', 'aerodyn', 'src')
DEFAULT_FILES = ['FVW.f90', 'FVW_Subs.f90', 'FVW_BiotSavart.f90', 'FVW_VortexTools.f90', 'FVW_Wings.f90', 'FVW_IO.f90']

RE_START = re.compile(r'^\s*(?:(?:pure|elemental|recursive|module|impure)\s+)*'
                      r'(?:(?:real|integer|logical|character|type)\s*(?:\([^)]*\))?\s*)?'
                      r'(subroutine|function)\s+(\w+)', re.I)
RE_END = re.compile(r'^\s*end\s*(subroutine|function)\b', re.I)
RE_DO = re.compile(r'^\s*(?:\w+\s*:\s*)?do\b(?!\s*while)', re.I)
RE_DOWHILE = re.compile(r'^\s*(?:\w+\s*:\s*)?do\s+while\b', re.I)
RE_ENDDO = re.compile(r'^\s*end\s*do\b', re.I)
RE_CALL = re.compile(r'\bcall\s+(\w+)', re.I)
RE_OMP = re.compile(r'^\s*!\$omp', re.I)
RE_ALLOC = re.compile(r'\b(?<!de)allocate\s*\(|\bAllocAry\s*\(', re.I)
RE_IHEAD = re.compile(r'\b(iHead\w*)\s*=\s*\1\s*\+', re.I)
RE_SELECT = re.compile(r'^\s*select\s+case', re.I)
RE_STOP = re.compile(r'^\s*(stop\b|print\b|write\s*\(|call\s+WrScr)', re.I)
RE_PTR = re.compile(r'\bpointer\b|=>', re.I)
RE_REGTYPE = re.compile(r'type\s*\(\s*(FVW_\w+Type|T_Tree|T_Node|T_Sgmt|T_Part|T_Panl|GridOutType)\s*\)', re.I)
RE_DEVVER = re.compile(r'^\s*if\s*\(\s*DEV_VERSION', re.I)
# Loop bounds that scale with the wake (N markers / segments / particles / CPs), as opposed to
# the blade span (nSpan ~ 10-100) or the number of wings. Only wake-scale loops are worth a GPU.
RE_DOBOUND = re.compile(r'^\s*(?:\w+\s*:\s*)?do\s+\w+\s*=\s*(.*)$', re.I)
RE_WAKE = re.compile(r'\b(nNW\w*|nFW\w*|nSeg\w*|nCPs|icp_?end|iCPEnd|iSegEnd|nPart|Part%n|nDepth|nP|nL|iEnd|'
                     r'nx|ny|nz|size\(\s*\w+\s*,\s*[23]\s*\))', re.I)
NOISE_CALLS = {'seterrstat', 'allocary', 'tic', 'toc', 'wrscr', 'cleanup', 'find_nan_2d', 'find_nan_1d', 'find_nan_3d'}


def strip_comment(line):
    if RE_OMP.match(line):
        return line
    out, q = [], None
    for ch in line:
        if q:
            if ch == q:
                q = None
        elif ch in '\'"':
            q = ch
        elif ch == '!':
            break
        out.append(ch)
    return ''.join(out)


def scan(path):
    routines, stack = [], []
    with open(path, errors='replace') as f:
        lines = f.readlines()
    for n, raw in enumerate(lines, 1):
        line = strip_comment(raw)
        if not line.strip():
            continue
        m = RE_START.match(line)
        if m and not re.match(r'^\s*end', line, re.I):
            r = dict(file=os.path.basename(path), name=m.group(2), start=n, kind=m.group(1).lower(),
                     recursive=bool(re.search(r'\brecursive\b', line, re.I)), parent=stack[-1]['name'] if stack else None,
                     omp=0, alloc=0, ihead=0, select=0, io=0, ptr=0, regtypes=set(), calls=set(),
                     depth=0, maxdepth=0, dowhile=0, devversion=0, wake=0)
            stack.append(r)
            continue
        if RE_END.match(line) and stack:
            r = stack.pop()
            r['end'] = n
            r['lines'] = n - r['start'] + 1
            routines.append(r)
            continue
        if not stack:
            continue
        r = stack[-1]
        if RE_OMP.match(line):
            r['omp'] += 1
            continue
        if RE_DOWHILE.match(line):
            r['dowhile'] += 1
            r['depth'] += 1
        elif RE_DO.match(line):
            b = RE_DOBOUND.match(line)
            if b and RE_WAKE.search(b.group(1)):
                r['wake'] += 1
            r['depth'] += 1
            r['maxdepth'] = max(r['maxdepth'], r['depth'])
        elif RE_ENDDO.match(line) or re.match(r'^\s*enddo\b', line, re.I):
            r['depth'] = max(0, r['depth'] - 1)
        # one-line "do ...; ...; enddo;" nests
        inline = len(re.findall(r';\s*do\b', line, re.I))
        if inline:
            r['maxdepth'] = max(r['maxdepth'], r['depth'] + inline)
        r['alloc'] += len(RE_ALLOC.findall(line))
        r['ihead'] += len(RE_IHEAD.findall(line))
        r['select'] += bool(RE_SELECT.match(line))
        r['io'] += bool(RE_STOP.match(line))
        r['ptr'] += bool(RE_PTR.search(line))
        r['devversion'] += bool(RE_DEVVER.match(line))
        for t in RE_REGTYPE.findall(line):
            r['regtypes'].add(t)
        for c in RE_CALL.findall(line):
            if c.lower() not in NOISE_CALLS:
                r['calls'].add(c)
    return routines


def classify(r):
    blockers = []
    if r['ihead']:
        blockers.append('iHead counter')
    if r['recursive']:
        blockers.append('recursion')
    if r['ptr'] and any(t.lower() in ('t_tree', 't_node') for t in r['regtypes']):
        blockers.append('pointer tree')
    io_heavy = r['io'] - r['devversion'] > 2
    if r['file'] == 'FVW_IO.f90' or io_heavy and r['maxdepth'] < 2:
        return 'HOST', blockers
    if r['name'].lower().startswith(('print', 'find_', 'have_nan')):
        return 'HOST', blockers                        # debug helpers
    if r['wake'] and not io_heavy:
        return ('DEVICE*' if blockers else 'DEVICE'), blockers
    if len(r['calls']) >= 3 and r['regtypes']:
        return 'GLUE', blockers
    if r['maxdepth'] >= 1 and not r['regtypes'] and r['io'] == 0 and r['recursive']:
        return 'DEVICE*', blockers                     # tree traversal: per-CP work, wake-scale by construction
    return 'HOST', blockers


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='*')
    ap.add_argument('--md', action='store_true', help='markdown table')
    ap.add_argument('--json', help='write JSON here')
    ap.add_argument('--min-lines', type=int, default=8, help='skip routines shorter than this (Failed/CleanUp helpers)')
    a = ap.parse_args()
    files = a.files or [os.path.join(DEFAULT_SRC, f) for f in DEFAULT_FILES]
    rows = []
    for f in files:
        for r in scan(f):
            if r['lines'] < a.min_lines:
                continue
            r['hint'], r['blockers'] = classify(r)
            rows.append(r)
    rows.sort(key=lambda r: (r['file'], r['start']))
    if a.json:
        with open(a.json, 'w') as fo:
            json.dump([{**r, 'regtypes': sorted(r['regtypes']), 'calls': sorted(r['calls'])} for r in rows], fo, indent=1)
    if a.md:
        print('| File:line | Routine | Lines | Nest | Wake loops | OMP | Alloc | Hint | Blockers |')
        print('|---|---|---|---|---|---|---|---|---|')
        for r in rows:
            print(f"| {r['file']}:{r['start']} | `{r['name']}` | {r['lines']} | {r['maxdepth']} | {r['wake']} | {r['omp']} | {r['alloc']} "
                  f"| {r['hint']} | {', '.join(r['blockers'])} |")
    else:
        print(f"{'file:line':28s} {'routine':34s} {'lines':>5s} {'nest':>4s} {'wake':>4s} {'omp':>3s} {'alloc':>5s} {'ihead':>5s} hint      blockers")
        for r in rows:
            print(f"{r['file']+':'+str(r['start']):28s} {r['name'][:34]:34s} {r['lines']:5d} {r['maxdepth']:4d} {r['wake']:4d} {r['omp']:3d} "
                  f"{r['alloc']:5d} {r['ihead']:5d} {r['hint']:9s} {', '.join(r['blockers'])}")
        from collections import Counter
        c = Counter(r['hint'] for r in rows)
        print(f"\n{len(rows)} routines: " + ', '.join(f'{k} {v}' for k, v in sorted(c.items())), file=sys.stderr)


if __name__ == '__main__':
    main()
