#!/usr/bin/env python3
"""
What one ulp costs: compare job_restart_cadence's leg P (checkpointing every
EVERY steps, canonicalising) against leg R (the same cadence with
FESOM_RESTART_CANONICALIZE=0) at each checkpoint, so the growth of the
perturbation can be read rather than assumed.

Leg R is bit-identical to leg Q (one checkpoint at the end), so R is the
un-perturbed trajectory and every P-R difference is the canonicalisation.

Run after job_restart_cadence; takes a few minutes, one comparison per
checkpoint. Point counts matter as much as the maxima: the whole wet domain
differs ten steps in, which is the global SSH solve landing inside its own 1e-5
tolerance, not chaos spreading from the partition boundaries.
"""
import netCDF4 as nc, numpy as np, glob, os
base = "/work/ab0995/a270088/port_paper_v2/runs/restart_cadence"
print(f"{'step':>5} {'temp max|d|':>12} {'temp npts':>10} {'u max|d|':>12} {'u npts':>10}")
for i, f in enumerate(sorted(glob.glob(base + "/P/*.restart.nc"))):
    b = os.path.basename(f)
    A = nc.Dataset(f); B = nc.Dataset(base + "/R/" + b)
    out = []
    for v in ("temp", "u"):
        a = np.asarray(A[v][:], dtype=np.float64)
        c = np.asarray(B[v][:], dtype=np.float64)
        d = np.abs(a - c)
        out += [d.max(), int((d != 0).sum())]
    A.close(); B.close()
    print(f"{(i+1)*10:>5} {out[0]:12.3e} {out[1]:>10} {out[2]:12.3e} {out[3]:>10}", flush=True)
