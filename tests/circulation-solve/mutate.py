#!/usr/bin/env python3
"""Mutation checks of unit circulation-solve's L0 (criteria: L0 sensitivity). Works on a COPY of src/,
never on the repo (another build may read src/). Each mutation: patch the copy, rebuild l0_cs, run it, restore.
  mutate.py <copy root> <build dir of tests/circulation-solve configured with -DOFK_SRC_OVERRIDE... (see unit.md)>
The copy's tree must be what the build compiles (CMake ofk_sources.cmake resolves src/ relative to itself)."""
import os, re, subprocess, sys
root, build = sys.argv[1], sys.argv[2]
K = os.path.join(root, 'src/modules/kokkos/src', 'circulation-solve')
BODY, ENTRY, COMMON = (os.path.join(K, f) for f in ('circulation_body.hpp', 'circulation.cpp', 'circulation_common.hpp'))
M = [  # name, file, old, new, expectation
 ('wrong regularization branch (LambOseen runs Vatistas)', BODY,
  "case kRegLambOseen: build_reg<kRegLambOseen>(g, par); break;", "case kRegLambOseen: build_reg<kRegVatistas>(g, par); break;", 'claim 1 fails'),
 ('Gamma order: fourpi_inv*Kv first, then Gamma (F10-R1 style, C2)', BODY,
  "        T K = G * f;\n        K = K * kv[l];", "        T K = f * kv[l];\n        K = K * G;", 'claim 1 fails'),
 ('cache not invalidated: reused after the first build (Eps / geometry changes, C12)', ENTRY,
  "  const auto t0 = clk::now();\n  if (h.mode == static_cast<int>(ofk::Mode::strict)) strict::build(g, h.parBuild);",
  "  if (h.st.nBuild > 0) { h.built = true; return; }\n  const auto t0 = clk::now();\n  if (h.mode == static_cast<int>(ofk::Mode::strict)) strict::build(g, h.parBuild);", 'A7 fails'),
 ('off-by-one: last ring (kCP = nCP_tot) not replayed', BODY,
  "for (int kc = 0, r = 0; kc < g.nT; ++kc)", "for (int kc = 0, r = 0; kc + 1 < g.nT; ++kc)", 'claim 1 fails'),
 ('segment 3-4 taken as 4-3 (C5)', BODY,
  "constexpr int kSegA[4] = {0, 2, 1, 3};\nconstexpr int kSegB[4] = {1, 3, 2, 0};",
  "constexpr int kSegA[4] = {0, 3, 1, 3};\nconstexpr int kSegB[4] = {1, 2, 2, 0};", 'claim 1 fails'),
 ('Eps_NW row 2 instead of row 1 (C7)', ENTRY,
  "P[12] = eNW[std::size_t(ldE1) * (", "P[12] = eNW[1 + std::size_t(ldE1) * (", 'claim 1 fails'),
 ('guards dropped: early-return segments added (C3)', BODY,
  "so[t] = skip ? T(-1) : norm_a + norm_b;", "so[t] = norm_a + norm_b;", 'claim 1 fails'),
 ('Gamma index off by one depth (ring r uses GammaLastIter(r/iNWStart) of the next station at depth 2)', BODY,
  "      const int j = 4 * r + sg;\n      const SegEnds<T> e = seg_ends(g, j);\n      const std::size_t o = std::size_t(j) * np",
  "      const int j = 4 * r + sg;\n      const T Gx = (d == 1 && kc + 1 < g.nT) ? gam[kc + 1] : G; (void)Gx;\n      const SegEnds<T> e = seg_ends(g, j);\n      const std::size_t o = std::size_t(j) * np", 'claim 1 fails (if Gx used)'),
]
# the last one only compiles a dead variable unless G is replaced; make it effective:
M[-1] = (M[-1][0], BODY, "        T K = G * f;\n        K = K * kv[l];", "        T K = ((d == 1 && kc + 1 < g.nT) ? gam[kc + 1] : G) * f;\n        K = K * kv[l];", 'claim 1 fails')
env0 = dict(os.environ, OMP_NUM_THREADS='4', OFK_L0_CASES='60', OFK_L0_INFO='1')
def run():
    b = subprocess.run(['cmake', '--build', build, '-j', '12', '--target', 'l0_cs'], capture_output=True, text=True)
    if b.returncode:
        return 'BUILD FAILED: ' + b.stdout[-800:]
    r = subprocess.run([os.path.join(build, 'l0_cs')], capture_output=True, text=True, env=env0)
    out = r.stdout + r.stderr
    v = re.findall(r'variant (strict|F0-parfor)\s*: cases (\d+), with any differing value (\d+), with a different iteration count (\d+)', out)
    s = re.search(r'l0_cs: (PASS|FAIL)[^\n]*', out)
    a7 = len(re.findall(r'FAIL A7', out))
    return f"{s.group(0) if s else '?'} | " + '; '.join(f'{n}: {d}/{c} cases differ, {i} with other iteration counts' for n, c, d, i in v) + f' | A7 failures: {a7}'
for name, f, old, new, expect in M:
    src = open(f).read()
    assert src.count(old) == 1, (name, src.count(old))
    open(f, 'w').write(src.replace(old, new))
    try:
        print(f'MUTATION {name} (expect: {expect})\n    {run()}', flush=True)
    finally:
        open(f, 'w').write(src)
print('UNMUTATED:', run())
