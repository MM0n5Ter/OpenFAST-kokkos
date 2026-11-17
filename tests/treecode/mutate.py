#!/usr/bin/env python3
"""Mutation checks of unit treecode's L0 (unit.md 'Exact next steps' 4). Works on a COPY of src/
(opt/build/u3d-mut-src), never on the repo. Results of unit 3: opt/build/u3d-mutations2.log. Each mutation: patch, rebuild l0_tc, run, restore.
  mutate.py <copy root> <build dir>"""
import os, subprocess, sys, re
root, build = sys.argv[1], sys.argv[2]
K = os.path.join(root, 'src/modules/kokkos/src')
IMPL, IV, DB = (os.path.join(K, f) for f in ('treecode/treecode_impl.hpp', 'induced-velocities/induced_velocities.cpp', 'treecode/treecode_dbuild.hpp'))
M = [  # name, file, old, new, env, expectation
 ('octant > -> >= (host build)', IMPL, "(px(i) > G[0] ? 1 : 0)", "(px(i) >= G[0] ? 1 : 0)", {}, 'claim 1 fails'),
 ('powi r^5 = r^4*r', IMPL, "const T r5 = r2 * r3;", "const T r5 = r4 * r;", {}, 'claim 1 fails'),
 ('pow -> sqrt in strict build', IMPL, "wLoc = std::pow(s, T(0.5));", "wLoc = std::sqrt(s);", {}, 'claim 1 fails'),
 ('singleton leaves in reverse octant order', IMPL, "for (int o = 0; o < 8; ++o)\n        if (cnt[o] == 1) push_leaf", "for (int o = 7; o >= 0; --o)\n        if (cnt[o] == 1) push_leaf", {}, 'claim 1 fails'),
 ('degenerate threshold hard-coded 1.0e-5f', IV, "tc.degenerate = thr;", "tc.degenerate = T(1.0e-5f);", {}, 'claim 1 fails'),
 ('dead-node tolerance 0', IV, "tc.eqTol = T(100) * std::numeric_limits<T>::epsilon() / T(2);", "tc.eqTol = T(0);", {}, 'claim 1 fails'),
 ('NCRIT 32 -> 31', IV, "tc.ncrit = 32;", "tc.ncrit = 31;", {}, 'claim 1 fails'),
 ('fast only: MAC radius halved', IMPL, "const T dist = fmax_f(bf * nrv(i, 3) + ff * nrv(i, 4), dd0);", "const T dist = fmax_f((kStrict ? bf : T(0.5) * bf) * nrv(i, 3) + ff * nrv(i, 4), dd0);", {}, 'R2 fails'),
 ('device build: octant > -> >=', DB, "const int o = (x > g0 ? 1 : 0)", "const int o = (x >= g0 ? 1 : 0)", {'OFK_TC_DEVBUILD': '1', 'OFK_TC_DEVBUILD_CHECK': '1'}, 'topology check fails'),
 ('device build: singletons reverse order', DB, "for (int o = 0; o < 8; ++o) {\n            const int cnt = ri(k, I_CNT + o);\n            if (cnt == 1) {", "for (int o = 7; o >= 0; --o) {\n            const int cnt = ri(k, I_CNT + o);\n            if (cnt == 1) {", {'OFK_TC_DEVBUILD': '1', 'OFK_TC_DEVBUILD_CHECK': '1'}, 'topology check fails'),
 ('device build: child center sign flipped (x)', DB, "rg2(j, 0) = g0 + q5 * ((o & 1) ? T(1) : T(-1));", "rg2(j, 0) = g0 + q5 * ((o & 1) ? T(-1) : T(1));", {'OFK_TC_DEVBUILD': '1', 'OFK_TC_DEVBUILD_CHECK': '1'}, 'topology check fails'),
 ('device build, no check: stored node radius x0.25 (MAC too loose)', DB, "rnr(k, 3) = R;", "rnr(k, 3) = T(0.25) * R;", {'OFK_TC_DEVBUILD': '1'}, 'R2 fails'),
]
env0 = dict(os.environ, OMP_NUM_THREADS='4', OFK_L0_BIG='0', OFK_L0_INFO='1', OFK_L0_RTOL='0')
for name, f, old, new, env, expect in M:
    src = open(f).read()
    assert src.count(old) == 1, (name, src.count(old))
    open(f, 'w').write(src.replace(old, new))
    try:
        b = subprocess.run(['cmake', '--build', build, '-j', '12', '--target', 'l0_tc'], capture_output=True, text=True)
        if b.returncode:
            print(f'MUTATION {name}: BUILD FAILED\n{b.stdout[-2000:]}'); continue
        r = subprocess.run([os.path.join(build, 'l0_tc')], capture_output=True, text=True, env=dict(env0, **env))
        out = r.stdout + r.stderr
        c1 = re.search(r'claim 1/2 .*', out); r2 = re.search(r'R2: .*', out); s = re.search(r'l0_tc: \d+ passed, \d+ failed', out)
        chk = 'OFK_TC_DEVBUILD_CHECK' in out and 'differs' in out
        print(f'MUTATION {name} (expect: {expect}) rc={r.returncode}')
        print('   ', c1.group(0) if c1 else '-', '|', (r2.group(0)[:150] if r2 else '-'), '|', s.group(0) if s else '-',
              '| topology check:', 'FAILED (' + re.search(r'differs[^\n]*', out).group(0)[:90] + ')' if chk else 'ok/not run')
    finally:
        open(f, 'w').write(src)
b = subprocess.run(['cmake', '--build', build, '-j', '12', '--target', 'l0_tc'], capture_output=True, text=True)
r = subprocess.run([os.path.join(build, 'l0_tc')], capture_output=True, text=True, env=env0)
print('UNMUTATED:', re.search(r'claim 1/2 .*', r.stdout).group(0), '|', re.search(r'l0_tc: \d+ passed, \d+ failed', r.stdout).group(0))
