#!/usr/bin/env python3
"""Write a small periodic (FileID 8) TurbSim .bts (Grid3D) wind file for the 4b L0 (BAR rotor fits: z 5..265 m, y +-130 m).
   mkbts.py <out.bts>   Deterministic (no RNG): shear + travelling sine gusts, so every grid cell differs."""
import struct, sys, math
nz, ny, nt, dz, dy, dt, uhub, zhub, zb = 27, 27, 64, 10.0, 10.0, 0.1, 8.0, 140.0, 5.0
slope, off = [1000.0]*3, [0.0]*3
desc = b'ofk 4b L0 synthetic Grid3D'
with open(sys.argv[1], 'wb') as f:
    f.write(struct.pack('<h4i', 8, nz, ny, 0, nt))
    f.write(struct.pack('<6f', dz, dy, dt, uhub, zhub, zb))
    f.write(struct.pack('<6f', slope[0], off[0], slope[1], off[1], slope[2], off[2]))
    f.write(struct.pack('<i', len(desc)) + desc)
    for it in range(nt):
        for iz in range(nz):
            z = zb + iz*dz
            for iy in range(ny):
                y = -0.5*(ny-1)*dy + iy*dy
                ph = 2*math.pi*it/nt
                u = uhub*(z/zhub)**0.2 + 0.6*math.sin(ph + 0.05*y + 0.03*z)
                v = 0.4*math.sin(2*ph + 0.02*z)
                w = 0.3*math.cos(ph + 0.04*y)
                f.write(struct.pack('<3h', *(int(round(c*s)) for c, s in zip((u, v, w), slope))))
