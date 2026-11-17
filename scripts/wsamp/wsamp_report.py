#!/usr/bin/env python3
"""Resolve wsamp stacks and print self/inclusive time per routine, and bucket totals.
usage: wsamp_report.py <prefix> <hz> [--top N] [--buckets] [--callers SYM] [--children SYM]"""
import sys, struct, subprocess, bisect, collections, re, os
pre, hz = sys.argv[1], float(sys.argv[2]); args = sys.argv[3:]
top = int(args[args.index('--top')+1]) if '--top' in args else 60
# maps
regions = []  # (start, end, off, path)
for l in open(pre + '.maps'):
    p = l.split()
    if len(p) < 6 or 'x' not in p[1]: continue
    s, e = (int(x, 16) for x in p[0].split('-')); regions.append((s, e, int(p[2], 16), p[5]))
regions.sort()
symcache = {}
def syms(path):
    if path in symcache: return symcache[path]
    out = []
    for flags in (['-n', '--defined-only'], ['-n', '--defined-only', '-D']):
        try: r = subprocess.run(['nm', *flags, path], capture_output=True, text=True).stdout
        except Exception: r = ''
        for l in r.splitlines():
            q = l.split(None, 2)
            if len(q) == 3 and q[1] in 'tTwWiI': out.append((int(q[0], 16), q[2]))
        if out: break
    out.sort(); symcache[path] = ([a for a, _ in out], [n for _, n in out]); return symcache[path]
def is_exec(path):
    with open(path, 'rb') as f: h = f.read(18)
    return struct.unpack('<H', h[16:18])[0] == 2
execcache = {}
rescache = {}
def resolve(pc):
    if pc in rescache: return rescache[pc]
    i = bisect.bisect_right(regions, (pc, 1 << 64)) - 1
    name = '??'
    if i >= 0 and regions[i][0] <= pc < regions[i][1]:
        s, e, off, path = regions[i]
        if os.path.exists(path):
            if path not in execcache: execcache[path] = is_exec(path)
            a = pc if execcache[path] else pc - s + off
            ad, nm = syms(path); j = bisect.bisect_right(ad, a) - 1
            name = (nm[j] if j >= 0 else '??') + ('' if execcache[path] else '@' + os.path.basename(path).split('.so')[0])
    rescache[pc] = name; return name
def clean(n):
    m = re.match(r'__(\w+?)_MOD_(\w+)', n)
    if m: n = m.group(2) + ' [' + m.group(1) + ']'
    n = re.sub(r'\._omp_fn\.\d+', ' (omp)', n)
    return n
data = open(pre + '.stacks', 'rb').read(); pos = 0; stacks = []
while pos < len(data):
    d, = struct.unpack_from('<i', data, pos); pos += 4
    pcs = struct.unpack_from('<%dQ' % d, data, pos); pos += 8 * d
    st = [clean(resolve(pc if k == 0 else pc - 1)) for k, pc in enumerate(pcs)]
    k = next((i for i, x in enumerate(st) if x.endswith('@libc') and ('sigaction' in x or 'restore_rt' in x or 'killpg' in x)), -1)
    if k >= 0 and k + 1 < len(st):  # frame after the trampoline is the interrupted PC (no -1 adjustment)
        st = [clean(resolve(pcs[k + 1]))] + st[k + 2:]
    stacks.append(st)
N = len(stacks); dt = 1.0 / hz
print(f'samples {N} = {N*dt:.2f} s wall (main thread)')
selfc = collections.Counter(s[0] for s in stacks if s)
incl = collections.Counter()
for s in stacks:
    for f in set(s): incl[f] += 1
if '--buckets' not in args or '--all' in args:
    print('\n== inclusive (top %d)' % top)
    for f, c in incl.most_common(top): print(f'{c*dt:9.2f} s {100*c/N:6.1f}%  self {selfc[f]*dt:8.2f}  {f}')
    print('\n== self (top %d)' % top)
    for f, c in selfc.most_common(top): print(f'{c*dt:9.2f} s {100*c/N:6.1f}%  {f}')
for key in ('--callers', '--children'):
    if key in args:
        tgt = args[args.index(key)+1]; cc = collections.Counter()
        for s in stacks:
            for k, f in enumerate(s):
                if f.startswith(tgt):
                    o = s[k+1] if key == '--callers' else (s[k-1] if k > 0 else '<self>')
                    if key == '--callers' and k+1 >= len(s): o = '<top>'
                    cc[o] += 1; break
        print(f'\n== {key[2:]} of {tgt}')
        for f, c in cc.most_common(30): print(f'{c*dt:9.2f} s {100*c/N:6.1f}%  {f}')
# buckets: first rule (in order) whose pattern matches any frame of the stack (frames lower-cased)
KER = r'^ofk_iv_eval|^ofk_tc_ui|^ofk_iv_tree_build|ofk_bsd::|ofk_tc::|^ui_seg \[|^ui_part_nograd|^ui_tree_part|^grow_tree_part|^ui_leaves|^ui_tree_segment'
B = [
 ('out', 'output/VTK/grid', [r'^writevtkoutputs']),
 ('kW', 'kernel: wake (tree build+eval / direct)', [r'^wakeinducedvelocities', KER]),
 ('kL', 'kernel: LL induction', [r'^liftinglineinducedvelocities', KER]),
 ('circ', 'circ. solve: ui_quad_n1 iterations', [r'^ui_quad_n1']),
 ('pkW', 'pack/unpack: wake (Fortran pack + ofk sources)', [r'^wakeinducedvelocities']),
 ('pkL', 'pack/unpack: LL (Fortran pack + ofk sources)', [r'^liftinglineinducedvelocities']),
 ('circR', 'circ. solve: rest (polar lookup, loop)', [r'^wings_computecirculation']),
 ('ccsd', 'CCSD means/dxdt/eps', [r'^fvw_calccontstatederiv']),
 ('upd', 'Euler/RK4 update', [r'^fvw_euler1|^fvw_rk4']),
 ('prop', 'PropagateWake', [r'^propagatewake']),
 ('fvw', 'FVW other (panelling, maps, interp, copies)', [r'^fvw_updatestates|^fvw_calcoutput|^fvw_end|^fvw_init']),
 ('wind', 'AD_CalcWind (IfW at wake + LL pts)', [r'^ad_calcwind']),
 ('twr', 'SetInputsForFVW (TwrInflArray on wake pts)', [r'^setinputsforfvw']),
 ('ad', 'AeroDyn other', [r'\[aerodyn\]|\[aerodyn_inflow\]']),
 ('init', 'init / driver / I/O / other', [r'.']),
]
cb = collections.Counter(); ex = collections.defaultdict(collections.Counter)
for s_ in stacks:
    lw = [f.lower() for f in s_]
    for key, nm, pats in B:
        if all(any(re.search(p, f) for f in lw) for p in pats): cb[key] += 1; ex[key][s_[0] if s_ else '?'] += 1; break
if '--csv' in args:
    lab = args[args.index('--csv')+1]
    print('csv,' + lab + ',' + f'{N*dt:.2f},' + ','.join(f'{cb[k]*dt:.2f}' for k, _, _ in B))
    print('csvhdr,label,total,' + ','.join(k for k, _, _ in B))
else:
    print('\n== buckets (first matching rule, top-down priority)')
    for key, nm, _ in B:
        c = cb[key]; print(f'{c*dt:9.2f} s {100*c/N:6.1f}%  {nm}   [{", ".join(f"{k[:60]}:{v*dt:.2f}" for k, v in ex[key].most_common(3))}]')
