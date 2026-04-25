#!/usr/bin/env python3
"""Per-column histogram of |diff| at Q point — find true SST disagreement."""
import glob
import sys

import numpy as np

ref_dir = sys.argv[1]
tst_dir = sys.argv[2]


def load(directory):
    rows = []
    for f in sorted(glob.glob(f"{directory}/evp_dump_step1_Q_node_rank*.txt")):
        rows.append(np.loadtxt(f, comments="#", ndmin=2))
    a = np.concatenate(rows)
    a = a[np.argsort(a[:, 0].astype(np.int64))]
    # Dedup by gid (overlap-partition)
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx]


ref = load(ref_dir)
tst = load(tst_dir)
print(f"ref={len(ref)} tst={len(tst)}")
gid_r = ref[:, 0].astype(np.int64)
gid_t = tst[:, 0].astype(np.int64)
common = np.intersect1d(gid_r, gid_t)
ref = ref[np.isin(gid_r, common)]
tst = tst[np.isin(gid_t, common)]
ref = ref[np.argsort(ref[:, 0])]
tst = tst[np.argsort(tst[:, 0])]
assert np.array_equal(ref[:, 0], tst[:, 0])

cols = ["a_ice", "m_ice", "m_snow", "elevation", "srfoce_temp"]
for k, name in enumerate(cols):
    a = ref[:, 1 + k]
    b = tst[:, 1 + k]
    d = np.abs(b - a)
    bins = [0, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1e-1, 1.0, 10.0, 1e9]
    h, _ = np.histogram(d, bins=bins)
    print(f"\n{name}: |diff| histogram")
    for lo, hi, c in zip(bins[:-1], bins[1:], h):
        if c > 0:
            print(f"  [{lo:9.1e}, {hi:9.1e})  : {c:8d}")
    if d.max() > 1e-9:
        idx = np.argsort(-d)[:5]
        print(f"  top-5 by |diff|:")
        for i in idx:
            g = int(ref[i, 0])
            print(f"    gid={g:7d}  ref={a[i]:.6e}  tst={b[i]:.6e}  diff={b[i]-a[i]:+.3e}")
