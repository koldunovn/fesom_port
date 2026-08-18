#!/usr/bin/env python3
"""Compare two restart files field by field. Exit 0 only if every field agrees
exactly.

    compare_restart.py A.restart.nc B.restart.nc

The gate for the restart port is a maximum absolute difference of 0: a run
interrupted and resumed must produce the same persistent state as the run that
was never interrupted. Anything else means a piece of state is missing from the
file, and the run is a different run from the one the reference did.
"""
import sys

import netCDF4
import numpy as np


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    a, b = (netCDF4.Dataset(p) for p in argv[1:])

    va, vb = set(a.variables), set(b.variables)
    bad = 0
    if va != vb:
        print(f"variable sets differ: only in A {sorted(va - vb)}, "
              f"only in B {sorted(vb - va)}")
        bad += 1

    for att in ("restart_format_version", "sum_nlevels_nod2D"):
        if a.getncattr(att) != b.getncattr(att):
            print(f"{att}: {a.getncattr(att)} vs {b.getncattr(att)}")
            bad += 1

    width = max((len(v) for v in va & vb), default=1)
    for name in sorted(va & vb):
        xa = np.asarray(a.variables[name][:], dtype=np.float64)
        xb = np.asarray(b.variables[name][:], dtype=np.float64)
        if xa.shape != xb.shape:
            print(f"{name:<{width}}  SHAPE {xa.shape} vs {xb.shape}")
            bad += 1
            continue
        d = np.abs(xa - xb)
        mx = float(d.max()) if d.size else 0.0
        nd = int((d != 0).sum())
        flag = "" if mx == 0.0 else "   <-- DIFFERS"
        print(f"{name:<{width}}  max|A-B| = {mx:.17g}  ({nd} of {d.size} entries){flag}")
        if mx != 0.0:
            bad += 1

    print()
    if bad:
        print(f"FAIL: {bad} field(s)/attribute(s) differ")
        return 1
    print("PASS: the restarted run reproduces the uninterrupted run exactly")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
