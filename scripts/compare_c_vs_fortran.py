#!/usr/bin/env python3
"""Side-by-side comparison: C-port monthly mean output vs FESOM-Fortran
monthly mean output for the same model year.

Both runs at CORE2 mesh, JRA55 1958 forcing, but possibly different dt
(C-port=500s, Fortran=1800s). We compare against monthly-mean fields:
agreement to 0.5°C / 0.5 PSU / 5% is "qualitatively correct".

Usage:
    python3 compare_c_vs_fortran.py <c_port_dir> <fortran_dir> [--year 1958]
"""
import argparse
import pathlib
import sys

import netCDF4 as nc
import numpy as np


FILL_THRESH = 1e30  # any |val| > this is a netCDF fill


def stats(arr):
    """Return (min, max, mean, std, count) over non-zero non-fill cells."""
    a = np.asarray(arr)
    mask = (a != 0) & (np.abs(a) < FILL_THRESH)
    if mask.sum() == 0:
        return (np.nan,) * 5
    sel = a[mask]
    return float(sel.min()), float(sel.max()), float(sel.mean()), float(sel.std()), int(mask.sum())


def diff_stats(f, c):
    """Compute per-cell diff stats where BOTH arrays have valid data."""
    a = np.asarray(f); b = np.asarray(c)
    # Handle Fortran-at-interfaces vs C-at-mid layer-count mismatch
    # by trimming to the common leading layers. e.g. Kv (Fortran=48, C=47).
    if a.ndim == 2 and b.ndim == 2 and a.shape[0] != b.shape[0]:
        n = min(a.shape[0], b.shape[0])
        a = a[:n]; b = b[:n]
    mask = (a != 0) & (b != 0) & (np.abs(a) < FILL_THRESH) & (np.abs(b) < FILL_THRESH)
    if mask.sum() == 0:
        return None
    d = a[mask] - b[mask]
    # spatial correlation across the valid mask
    af = a[mask].astype(np.float64); bf = b[mask].astype(np.float64)
    af -= af.mean(); bf -= bf.mean()
    denom = np.sqrt((af*af).sum() * (bf*bf).sum())
    corr = float((af*bf).sum() / denom) if denom > 0 else float('nan')
    return {
        "n":          int(mask.sum()),
        "mean":       float(d.mean()),
        "absmean":    float(np.abs(d).mean()),
        "max_abs":    float(np.abs(d).max()),
        "std":        float(d.std()),
        "corr":       corr,
    }


# (c_port_name, fortran_name); both expected to be (time, ...) shaped.
VARS_3D = [
    ("temp",    "temp"),
    ("salt",    "salt"),
    ("u",       "u"),
    ("v",       "v"),
    ("w",       "w"),
    ("Kv",      "Kv"),
    ("Av",      "Av"),
]
VARS_2D = [
    ("ssh",     "ssh"),
    ("sst",     "sst"),
    ("sss",     "sss"),
    ("a_ice",   "a_ice"),
    ("m_ice",   "m_ice"),
    ("m_snow",  "m_snow"),
    ("uice",   "uice"),
    ("vice",   "vice"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("c_port_dir", type=pathlib.Path)
    ap.add_argument("fortran_dir", type=pathlib.Path)
    ap.add_argument("--year", type=int, default=1958)
    args = ap.parse_args()

    print(f"C-port  dir : {args.c_port_dir}")
    print(f"Fortran dir : {args.fortran_dir}")
    print(f"year        : {args.year}")
    print()

    # Read each Fortran monthly file once (whole year). Slice per-month.
    # Read each C-port .monthly.nc once. Slice per-month.
    print(f"{'Var':6s} {'mon':>4s} | {'F range':>20s} {'F mean':>9s} | {'C range':>20s} {'C mean':>9s} | "
          f"{'Δmean':>9s} {'|Δ|mean':>9s} {'maxabs':>10s} {'corr':>7s}")
    print("-" * 130)

    for c_name, f_name in VARS_3D + VARS_2D:
        f_path = args.fortran_dir / f"{f_name}.fesom.{args.year}.nc"
        c_path = args.c_port_dir / f"{c_name}.fesom.{args.year}.monthly.nc"
        if not f_path.exists():
            print(f"{c_name:6s} {'-':>4s} | <missing Fortran file: {f_path.name}>")
            continue
        if not c_path.exists():
            print(f"{c_name:6s} {'-':>4s} | <missing C port file: {c_path.name}>")
            continue
        f_ds = nc.Dataset(f_path)
        c_ds = nc.Dataset(c_path)
        n_months = min(len(f_ds["time"]), len(c_ds["time"]))
        for m in range(n_months):
            f_arr = np.asarray(f_ds[f_name][m, ...])
            c_arr = np.asarray(c_ds[c_name][m, ...])
            # Fortran 2D vars may be (time, nod2); C may be (time, nod2). Both should align.
            f_min, f_max, f_mean, _, fN = stats(f_arr)
            c_min, c_max, c_mean, _, cN = stats(c_arr)
            d = diff_stats(f_arr, c_arr)
            if d is None:
                print(f"{c_name:6s} {m+1:>4d} | <no overlapping valid cells>")
                continue
            print(f"{c_name:6s} {m+1:>4d} | "
                  f"[{f_min:+8.2f},{f_max:+8.2f}] {f_mean:+9.4f} | "
                  f"[{c_min:+8.2f},{c_max:+8.2f}] {c_mean:+9.4f} | "
                  f"{d['mean']:+9.4f} {d['absmean']:9.4f} {d['max_abs']:10.3e} {d['corr']:7.4f}")
        print()


if __name__ == "__main__":
    main()
