"""R8 (criteria.md): trajectory drift per output channel against base1 (Fortran, OMP_NUM_THREADS=1).
  python r8.py <fork>/reg_tests/lib <base1.outb> <baseN.outb> <label>=<run.outb> [...]
For every channel: max |run - base1| relative to max|base1| of the channel, the envelope
max |baseN - base1| (same normalization), the ratio, and the first output step where
|run - base1| > 1e-6 * max|base1|. Prints the worst channels and the PASS rule
max|port - base1| <= 2 * max|baseN - base1| per channel (channels with a zero envelope are listed)."""
import sys
import numpy as np
sys.path.insert(0, sys.argv[1])
from fast_io import load_output  # noqa: E402

b1, info = load_output(sys.argv[2])[:2]
bN = load_output(sys.argv[3])[0]
names = info['attribute_names']
scale = np.maximum(np.max(np.abs(b1), axis=0), 1e-300)
env = np.max(np.abs(bN - b1), axis=0) / scale
print(f"base1 {sys.argv[2]}\nbaseN {sys.argv[3]}: channels={b1.shape[1]} steps={b1.shape[0]} "
      f"envelope max={env[1:].max():.3e} zero-envelope channels={(env[1:] == 0).sum()}")
for arg in sys.argv[4:]:
    lab, f = arg.split('=', 1)
    r = load_output(f)[0]
    if r.shape != b1.shape:
        print(f"{lab}: shape {r.shape} != {b1.shape}"); continue
    d = np.abs(r - b1) / scale
    dm = d.max(axis=0)
    first = np.array([np.argmax(d[:, k] > 1e-6) if (d[:, k] > 1e-6).any() else -1 for k in range(d.shape[1])])
    ok = dm <= 2 * env
    ch = range(1, d.shape[1])  # skip Time
    nfail = sum(1 for k in ch if not ok[k])
    nexceed = sum(1 for k in ch if first[k] >= 0)
    print(f"{lab}: max drift {dm[1:].max():.3e} (channel {names[1 + int(np.argmax(dm[1:]))]}); "
          f"channels > 2x envelope: {nfail}/{d.shape[1]-1}; channels ever > 1e-6: {nexceed}; "
          f"earliest step > 1e-6: {min([first[k] for k in ch if first[k] >= 0], default=-1)}")
    worst = sorted(ch, key=lambda k: -dm[k])[:8]
    for k in worst:
        print(f"    {names[k]:>22s}  drift {dm[k]:.3e}  envelope {env[k]:.3e}  ratio "
              f"{(dm[k] / env[k]) if env[k] > 0 else float('inf'):.2f}  first step >1e-6: {first[k]}")
