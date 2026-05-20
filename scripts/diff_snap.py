#!/usr/bin/env python3
"""Bit-identity gate for Tasks 3/4 of the IO-system overhaul.

Compares the snap_*.nc files in two directories field-by-field. Fields
are deemed bit-identical if their numpy arrays are equal under
`np.array_equal` (no tolerance — numerical drift indicates a code change
that altered physics).

Usage:
    python3 scripts/diff_snap.py <pre_dir> <post_dir>

Exit 0 on full match, 1 otherwise. Mismatches print field name + max
absolute difference + lat/lon of the worst node so a follow-up debug
run can probe the divergent column.
"""
import sys
import pathlib

import netCDF4 as nc
import numpy as np


def diff_file(pre_path: pathlib.Path, post_path: pathlib.Path) -> bool:
    if not post_path.exists():
        print(f"MISSING: {post_path}")
        return False
    with nc.Dataset(pre_path) as pre, nc.Dataset(post_path) as post:
        ok = True
        common = sorted(set(pre.variables) & set(post.variables))
        for name in common:
            a = np.asarray(pre.variables[name][...])
            b = np.asarray(post.variables[name][...])
            if a.shape != b.shape:
                print(f"  {name}: SHAPE {a.shape} vs {b.shape}")
                ok = False
                continue
            if np.array_equal(a, b):
                continue
            diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
            max_diff = float(diff.max())
            argmax_flat = int(diff.argmax())
            print(f"  {name}: max-abs-diff={max_diff:.3e}  argmax_flat={argmax_flat}")
            ok = False
        only_pre  = sorted(set(pre.variables)  - set(post.variables))
        only_post = sorted(set(post.variables) - set(pre.variables))
        if only_pre:  print(f"  vars only in pre:  {only_pre}");  ok = False
        if only_post: print(f"  vars only in post: {only_post}"); ok = False
    return ok


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    pre_dir  = pathlib.Path(sys.argv[1])
    post_dir = pathlib.Path(sys.argv[2])

    pre_files = sorted(pre_dir.glob("snap_*.nc"))
    if not pre_files:
        print(f"No snap_*.nc in {pre_dir}", file=sys.stderr)
        return 2

    print(f"Comparing {len(pre_files)} file(s) in {pre_dir} vs {post_dir}")
    print(f"Pre  SHA: {(pre_dir/'SHA.txt').read_text().strip()}"
          if (pre_dir/'SHA.txt').exists() else "Pre  SHA: (none)")
    print(f"Post SHA: {(post_dir/'SHA.txt').read_text().strip()}"
          if (post_dir/'SHA.txt').exists() else "Post SHA: (none)")
    print()

    all_ok = True
    for pre in pre_files:
        post = post_dir / pre.name
        print(f"== {pre.name} ==")
        if not diff_file(pre, post):
            all_ok = False
        print()
    if all_ok:
        print("ALL FIELDS BIT-IDENTICAL — Tasks 3/4 gate passes.")
        return 0
    else:
        print("DIVERGENCE — block on Task 5.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
