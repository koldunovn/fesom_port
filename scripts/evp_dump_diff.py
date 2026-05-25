#!/usr/bin/env python3
"""Diff EVP dump files between a 1-rank reference run and an N-rank test run.

Usage: evp_dump_diff.py <ref_dir> <tst_dir>

Per-step (s1, s2) and per-point (Q, U0, P, 1, 2, 3, 4, END) histograms of
|diff|. Report top-5 diverging gids per (step, point, column).
"""
import glob
import sys

import numpy as np

POINTS = [
    ("Q",   "node", ["a_ice", "m_ice", "m_snow", "elevation", "srfoce_temp"]),
    ("U0",  "node", ["u_ice", "v_ice"]),
    ("P",   "node", ["inv_mass", "inv_areamass", "rhs_a", "rhs_m"]),
    ("P",   "elem", ["ice_strength"]),
    ("F",   "node", ["stress_atmice_x", "stress_atmice_y", "u_w", "v_w"]),
    ("1",   "elem", ["sigma11", "sigma12", "sigma22"]),
    ("2",   "node", ["u_rhs", "v_rhs"]),
    ("3",   "node", ["u_ice", "v_ice"]),
    ("4",   "node", ["u_ice", "v_ice"]),
    ("END", "node", ["u_ice", "v_ice"]),
]
STEPS = [1, 2]
THRESHOLD = 1e-9    # below this, treat as FP-summation noise


def load(directory, step, point, cls):
    files = sorted(glob.glob(
        f"{directory}/evp_dump_s{step}_{point}_{cls}_rank*.txt"))
    if not files:
        return None
    rows = []
    for f in files:
        a = np.loadtxt(f, comments="#", ndmin=2)
        if a.size:
            rows.append(a)
    if not rows:
        return np.empty((0, 0))
    a = np.concatenate(rows)
    a = a[np.argsort(a[:, 0].astype(np.int64))]
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx]


def main(ref_dir, tst_dir):
    print(f"Ref:  {ref_dir}")
    print(f"Test: {tst_dir}")
    print()
    print(f"{'step':>4}  {'point':>4}  {'array':>5}  {'col':>14}  "
          f"{'max|diff|':>12}  {'max|val|':>12}  {'rel':>10}  {'#sig':>8}")
    print("-" * 92)
    for step in STEPS:
        any_dump_for_step = False
        for point, cls, hdr in POINTS:
            ref = load(ref_dir, step, point, cls)
            tst = load(tst_dir, step, point, cls)
            if ref is None or tst is None:
                continue
            any_dump_for_step = True
            common_g = np.intersect1d(ref[:, 0], tst[:, 0])
            ref_c = ref[np.isin(ref[:, 0], common_g)]
            tst_c = tst[np.isin(tst[:, 0], common_g)]
            ref_c = ref_c[np.argsort(ref_c[:, 0])]
            tst_c = tst_c[np.argsort(tst_c[:, 0])]
            for k, h in enumerate(hdr):
                a = ref_c[:, 1 + k]
                b = tst_c[:, 1 + k]
                d = np.abs(b - a)
                mad = d.max() if len(d) else 0.0
                mav = max(np.abs(a).max(), np.abs(b).max()) if len(a) else 0.0
                rel = mad / mav if mav > 0 else 0.0
                n_sig = int((d > THRESHOLD).sum())
                flag = "  <-- SIGNAL" if mad > THRESHOLD and (rel > 1e-10 or mav > 1e-12) else ""
                print(f"  s{step}  {point:>4}  {cls:>5}  {h:>14}  "
                      f"{mad:12.3e}  {mav:12.3e}  {rel:10.3e}  {n_sig:>8d}{flag}")
        if not any_dump_for_step:
            print(f"  s{step}: no dumps found")
        print()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
