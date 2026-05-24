#!/usr/bin/env python3
"""
KPP Phase K0 byte-identity check.

Compares the snapshot .nc files of two C-port runs and asserts every data
variable is bit-for-bit identical (max abs diff == 0.0). Used to prove that
adding the KPP scaffolding (struct alloc, dispatch, dbsfc array, dump harness)
does NOT change the default PP path — i.e. KPP-off == HEAD.

Usage:
    python3 kpp_byteident_check.py <run_dir_A> <run_dir_B>

Compares every snap_NNNNNN.nc present in BOTH dirs. Exit 0 iff identical.
"""
import sys, os, glob
import numpy as np
import netCDF4

# Coordinate / topology vars are static (mesh) — still compared (they must match),
# but listed here so a mismatch there is reported distinctly from physics drift.
COORD_VARS = {"lon", "lat", "zbar", "Z", "elem_nodes",
              "nlevels_nod2D", "nlevels", "time"}


def compare_file(pa, pb):
    da = netCDF4.Dataset(pa)
    db = netCDF4.Dataset(pb)
    names = sorted(set(da.variables) & set(db.variables))
    worst = 0.0
    worst_var = None
    rows = []
    for v in names:
        a = np.ma.filled(da.variables[v][:], np.nan).astype(np.float64)
        b = np.ma.filled(db.variables[v][:], np.nan).astype(np.float64)
        if a.shape != b.shape:
            rows.append((v, "SHAPE", a.shape, b.shape))
            worst = np.inf
            worst_var = v
            continue
        # treat NaN==NaN as equal (masked land); flag NaN-vs-number
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.abs(np.where(both_nan, 0.0, a - b))
        mismatched_nan = np.isnan(diff).sum()
        md = float(np.nanmax(diff)) if diff.size else 0.0
        if mismatched_nan:
            md = np.inf
        kind = "COORD" if v in COORD_VARS else "data"
        rows.append((v, kind, md, mismatched_nan))
        if md > worst:
            worst = md
            worst_var = v
    da.close()
    db.close()
    return worst, worst_var, rows


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    A, B = sys.argv[1], sys.argv[2]
    snaps_a = {os.path.basename(p) for p in glob.glob(os.path.join(A, "snap_*.nc"))}
    snaps_b = {os.path.basename(p) for p in glob.glob(os.path.join(B, "snap_*.nc"))}
    common = sorted(snaps_a & snaps_b)
    if not common:
        print(f"FAIL: no common snap_*.nc in\n  {A}\n  {B}")
        print(f"  A has {sorted(snaps_a)}\n  B has {sorted(snaps_b)}")
        sys.exit(1)

    overall = 0.0
    for snap in common:
        worst, worst_var, rows = compare_file(os.path.join(A, snap),
                                              os.path.join(B, snap))
        status = "IDENTICAL" if worst == 0.0 else "DIFFERS"
        print(f"\n=== {snap}: {status} (worst |Δ|={worst:.3e} in {worst_var}) ===")
        for v, kind, md, *extra in rows:
            if kind == "SHAPE":
                print(f"  {v:18s} SHAPE MISMATCH {md} vs {extra[0]}")
            elif md != 0.0:
                tag = " <-- NONZERO" if kind == "data" else " (coord)"
                print(f"  {v:18s} [{kind}] max|Δ|={md:.6e}{tag}")
        overall = max(overall, worst)

    print("\n" + "=" * 60)
    if overall == 0.0:
        print("RESULT: BYTE-IDENTICAL across all snapshots. PASS.")
        sys.exit(0)
    else:
        print(f"RESULT: NOT identical (worst |Δ|={overall:.6e}). FAIL.")
        sys.exit(1)


if __name__ == "__main__":
    main()
