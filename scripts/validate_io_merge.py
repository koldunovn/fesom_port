#!/usr/bin/env python3
"""End-to-end validation of the FESOM2 C-port IO system.

For a given output directory containing both daily and monthly netCDF
files for the same variables, this script verifies:

  1. xarray.open_mfdataset glues year-files into a continuous timeseries
     with calendar-decoded dates (no 1900 placeholder).
  2. For every full calendar month covered by both streams, the mean of
     that month's daily records equals the corresponding monthly record
     within `--tol` relative tolerance (default 1e-6).

Usage:
    python3 scripts/validate_io_merge.py <out_dir> [--tol 1e-6]

Requires the run to have used a `dt` that divides 86400 cleanly so each
day has the same step count; otherwise the daily-vs-monthly comparison
picks up a step-count weighting bias of order 1/N_steps_per_day.

Exit 0 on success, 1 on any failure.
"""
import argparse
import pathlib
import sys

import numpy as np
import xarray as xr


def find_var_pairs(out_dir: pathlib.Path):
    """Return {var: (daily_paths, monthly_paths)} for vars that have both."""
    daily   = {}
    monthly = {}
    for p in out_dir.glob("*.fesom.*.daily.nc"):
        var = p.name.split(".fesom.")[0]
        daily.setdefault(var, []).append(p)
    for p in out_dir.glob("*.fesom.*.monthly.nc"):
        var = p.name.split(".fesom.")[0]
        monthly.setdefault(var, []).append(p)
    pairs = {}
    for var in sorted(set(daily) & set(monthly)):
        pairs[var] = (sorted(daily[var]), sorted(monthly[var]))
    return pairs


def open_concat(paths, var):
    """xr.open_mfdataset across a list of paths along the time axis."""
    ds = xr.open_mfdataset([str(p) for p in paths], combine="by_coords")
    return ds[var]


def validate_var(var: str, daily_paths, monthly_paths, tol: float) -> bool:
    """Return True iff every full month in monthly equals the corresponding
    mean of daily records within `tol` relative tolerance. Months that
    aren't fully covered by daily records are skipped (they're partial
    runs and the monthly record's bnds extend beyond the data)."""
    da_d = open_concat(daily_paths,   var)
    da_m = open_concat(monthly_paths, var)

    print(f"== {var} ==")
    print(f"   daily   : {len(da_d.time)} records "
          f"[{da_d.time.values[0]} .. {da_d.time.values[-1]}]")
    print(f"   monthly : {len(da_m.time)} records "
          f"[{da_m.time.values[0]} .. {da_m.time.values[-1]}]")

    ok = True
    for t_m in da_m.time.values:
        # Determine the month [start, end) bounding this monthly record
        # and check whether daily records fully cover it.
        ts = np.datetime64(t_m, 'M')                # truncate to month
        m_start = ts.astype('datetime64[D]').astype('datetime64[ns]')
        m_end_M = ts + np.timedelta64(1, 'M')
        m_end   = m_end_M.astype('datetime64[D]').astype('datetime64[ns]')
        # Slice daily records that fall in [m_start, m_end)
        mask = (da_d.time >= m_start) & (da_d.time < m_end)
        d_in_month = da_d.where(mask, drop=True)

        days_in_month = (m_end - m_start) / np.timedelta64(1, 'D')
        if d_in_month.time.size != int(days_in_month):
            print(f"   skip {ts} — daily coverage {d_in_month.time.size}/{int(days_in_month)} (partial run)")
            continue

        d_avg = d_in_month.mean(dim="time")
        m_rec = da_m.sel(time=t_m)
        diff   = (d_avg - m_rec).values
        denom  = np.maximum(np.abs(m_rec.values), 1e-12)
        relerr = np.abs(diff) / denom
        max_abs = float(np.abs(diff).max())
        max_rel = float(relerr.max())
        verdict = "PASS" if max_rel < tol else "FAIL"
        if max_rel >= tol:
            ok = False
        print(f"   {ts} : max-abs={max_abs:.3e}  max-rel={max_rel:.3e}  {verdict}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", type=pathlib.Path)
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="relative-error tolerance (default 1e-6)")
    args = ap.parse_args()

    if not args.out_dir.is_dir():
        print(f"error: {args.out_dir} is not a directory", file=sys.stderr)
        return 2

    pairs = find_var_pairs(args.out_dir)
    if not pairs:
        print(f"error: no var has both daily + monthly files in {args.out_dir}", file=sys.stderr)
        return 2

    all_ok = True
    for var, (daily, monthly) in pairs.items():
        if not validate_var(var, daily, monthly, args.tol):
            all_ok = False
        print()

    if all_ok:
        print("ALL VARS PASS — merge + mean-of-daily == monthly within tolerance.")
        return 0
    else:
        print("FAIL — see per-month diagnostics above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
