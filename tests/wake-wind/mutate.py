#!/usr/bin/env python3
"""Mutation checks of the unit-4b L0 (each predicate broken once). Works on a COPY of the repo (src/ and tools/ copied,
deps/ linked), never on the repo. Each mutation: patch the copy, rebuild aerodyn_driver, run l0_ww.sh, restore.
  mutate.py <copy root> <build dir> <work dir> [mutation indices, default all]      verdict lines: MUT <name>: CAUGHT (<failing cases>) | MISSED"""
import os, shutil, subprocess, sys
root, build, work = sys.argv[1:4]
repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AD = 'src/modules/aerodyn/src/AeroDyn.f90'; FS = 'src/modules/aerodyn/src/FVW_Subs.f90'
M = [  # name, file, old, new, cases
 ('H4/H3 window edge: [tStart, tEnd] instead of [tStart-2DT, tEnd+DT]', FS,
  "- 2.0_DbKi*p%DTaero .and. &\n                     tq <= real(m%GridOutputs(iGrid)%tEnd,DbKi)   + p%DTaero",
  "- 0.0_DbKi*p%DTaero .and. &\n                     tq <= real(m%GridOutputs(iGrid)%tEnd,DbKi)", 'gedge gend base'),
 ('H3 slice boundary: LL slice one point short', AD, "      call slice(1, iLL_end)\n", "      call slice(1, iLL_end-1)\n", 'base'),
 ('H3 slice boundary: panel slice starts one point late', AD,
  "   call slice(iGrid_start(p%FVW%nGridOut+1), size(m%FVW%r_wind,2))", "   call slice(iGrid_start(p%FVW%nGridOut+1)+1, size(m%FVW%r_wind,2))", 'srcpnl'),
 ('W3 count: far wake left out of the packing', FS, "   if (p%nFWMax>0) then\n      do iW=1,p%nWings\n         iWake_end", "   if (.false.) then\n      do iW=1,p%nWings\n         iWake_end", 'fwlive'),
 ('H3 tIndx: wake slices kept for tIndx=2 instead of 1', AD, "   if (tIndx == 1) then\n      call slice(1, iWake_end)", "   if (tIndx == 2) then\n      call slice(1, iWake_end)", 'base'),
 ('H2 pass: wake skipped on the last pass (i=1) instead of i=2', AD, "ofkSkipWake=ofkWW)", "ofkSkipWake=ofkWW .and. i==1)", 'grid3d'),
]
os.makedirs(work, exist_ok=True)
for t in ('src', 'tools'):
    shutil.rmtree(os.path.join(root, t), ignore_errors=True); shutil.copytree(os.path.join(repo, t), os.path.join(root, t))
if not os.path.exists(os.path.join(root, 'deps')): os.symlink(os.path.join(repo, 'deps'), os.path.join(root, 'deps'))
os.makedirs(os.path.join(root, 'tests'), exist_ok=True)
sel = [int(x) for x in sys.argv[4:]] or range(len(M))
for k, (name, f, old, new, cases) in enumerate(M):
    if k not in sel: continue
    p = os.path.join(root, f); orig = open(os.path.join(repo, f)).read()
    assert orig.count(old) == 1, name
    open(p, 'w').write(orig.replace(old, new))
    b = subprocess.run(['bash', os.path.join(root, 'tools/openfast-build.sh'), 'omp', build],
                       env=dict(os.environ, OFK_TARGETS='aerodyn_driver'), capture_output=True, text=True)
    open(p, 'w').write(orig)
    if b.returncode:
        print(f'MUT {name}: BUILD FAILED (see {build}.build.log)'); continue
    r = subprocess.run(['bash', os.path.join(repo, 'tests/wake-wind/l0_ww.sh'), os.path.join(build, 'modules/aerodyn/aerodyn_driver'),
                        os.path.join(work, f'm{k}')] + cases.split(), capture_output=True, text=True)
    fails = [l.split()[2].rstrip(':') for l in r.stdout.splitlines() if ': FAIL' in l]
    print(f"MUT {name}: {'CAUGHT (' + ' '.join(fails) + ')' if fails else 'MISSED'}", flush=True)
