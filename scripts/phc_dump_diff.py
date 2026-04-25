#!/usr/bin/env python3
"""Compare post-PHC-load surface T,S between two runs (1-rank ref vs N-rank test).

Usage: phc_dump_diff.py <ref_dir> <tst_dir>
"""
import glob
import sys

import numpy as np


def load(d):
    rows = []
    for f in sorted(glob.glob(f"{d}/phc_dump_postload_rank*.txt")):
        rows.append(np.loadtxt(f, comments="#", ndmin=2))
    a = np.concatenate(rows)
    a = a[np.argsort(a[:, 0].astype(np.int64))]
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx]


ref = load(sys.argv[1])
tst = load(sys.argv[2])
print(f"ref={len(ref)} tst={len(tst)}")
common = np.intersect1d(ref[:, 0], tst[:, 0])
ref = ref[np.isin(ref[:, 0], common)]
tst = tst[np.isin(tst[:, 0], common)]
ref = ref[np.argsort(ref[:, 0])]
tst = tst[np.argsort(tst[:, 0])]
assert np.array_equal(ref[:, 0], tst[:, 0])
for k, name in enumerate(["T_surface", "S_surface"]):
    a = ref[:, 1 + k]; b = tst[:, 1 + k]
    d = np.abs(b - a)
    bins = [0, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1e-1, 1.0, 10.0, 1e9]
    h, _ = np.histogram(d, bins=bins)
    print(f"\n{name} |diff| histogram")
    for lo, hi, c in zip(bins[:-1], bins[1:], h):
        if c: print(f"  [{lo:9.1e},{hi:9.1e})  : {c:8d}")
    if d.max() > 1e-9:
        idx = np.argsort(-d)[:5]
        print(f"  top-5 by |diff|:")
        for i in idx:
            print(f"    gid={int(ref[i,0]):7d}  ref={a[i]:.6f}  tst={b[i]:.6f}  diff={b[i]-a[i]:+.4f}")
