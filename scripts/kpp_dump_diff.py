#!/usr/bin/env python3
"""Diff KPP dump files between the C port and the Fortran reference.

Usage: kpp_dump_diff.py <ref_dir> <tst_dir> [tag_substr]

Mirrors scripts/evp_dump_diff.py for the KPP dual-instrumentation
(reference_evp_dump_diagnostic). Both codes write, gated by FESOM_KPP_DUMP_DIR,
per-node files:

    kpp_dump_s<step>_<tag>_rank<R>.txt
    # step=.. tag=.. rank=R N=.. ncomp=C   (header)
    <gid> <v0> <v1> ... <v(C-1)>           (one line per owned node, 1-based gid)

The script auto-discovers every <tag> present in BOTH dirs, merges all ranks
(gid-keyed, deduped), intersects gids, and reports per-column max|diff|,
max|val|, relative diff, and the count above THRESHOLD. A `<-- SIGNAL` flag
marks real divergence (above FP-summation noise). Column names come from
COLS[tag] when known, else c0,c1,...

Bit-identical where the algebra is exact; <=1e-10 rel where lookup-table
interpolation or summation order legitimately differ (per the plan).
"""
import glob
import os
import re
import sys

import numpy as np

THRESHOLD = 1e-9   # below this, FP-summation noise

# Optional pretty column names for the small fixed-layout tags. Full-column
# tags (dVsq, dbsfc, viscA, diffKt, diffKs, ghats — one column per layer
# interface) fall back to c0,c1,... where cK == layer interface K (0-based).
# Keep in sync with the C dump call sites + the Fortran kpp_dump_* helpers.
COLS = {
    "prestep": ["ustar", "Bo"],
    "bldepth": ["hbl", "kbl", "bfsfc", "stable", "caseA"],
}

FNAME_RE = re.compile(r"kpp_dump_s(\d+)_(.+)_rank(\d+)\.txt$")


def discover(directory):
    """Return {(step, tag): [files]} for a dump dir."""
    out = {}
    for f in glob.glob(os.path.join(directory, "kpp_dump_s*_rank*.txt")):
        m = FNAME_RE.search(os.path.basename(f))
        if not m:
            continue
        step, tag = int(m.group(1)), m.group(2)
        out.setdefault((step, tag), []).append(f)
    return out


def load(files):
    rows = []
    for f in files:
        a = np.loadtxt(f, comments="#", ndmin=2)
        if a.size:
            rows.append(a)
    if not rows:
        return np.empty((0, 0))
    a = np.concatenate(rows)
    a = a[np.argsort(a[:, 0].astype(np.int64))]
    _, idx = np.unique(a[:, 0], return_index=True)   # dedup halo-shared gids
    return a[idx]


def main(ref_dir, tst_dir, tag_filter=None):
    print(f"Ref (Fortran?): {ref_dir}")
    print(f"Test (C?):      {tst_dir}")
    ref_index = discover(ref_dir)
    tst_index = discover(tst_dir)
    common = sorted(set(ref_index) & set(tst_index))
    if tag_filter:
        common = [(s, t) for (s, t) in common if tag_filter in t]
    if not common:
        print("No common (step, tag) dumps found.")
        print(f"  ref tags: {sorted({t for _, t in ref_index})}")
        print(f"  tst tags: {sorted({t for _, t in tst_index})}")
        sys.exit(1)

    print()
    print(f"{'step':>4}  {'tag':>10}  {'col':>12}  {'max|diff|':>12}  "
          f"{'max|val|':>12}  {'rel':>10}  {'#sig':>7}")
    print("-" * 84)
    worst = 0.0
    for step, tag in common:
        ref = load(ref_index[(step, tag)])
        tst = load(tst_index[(step, tag)])
        if ref.size == 0 or tst.size == 0:
            continue
        g = np.intersect1d(ref[:, 0], tst[:, 0])
        r = ref[np.isin(ref[:, 0], g)]
        t = tst[np.isin(tst[:, 0], g)]
        r = r[np.argsort(r[:, 0])]
        t = t[np.argsort(t[:, 0])]
        ncol = min(r.shape[1], t.shape[1]) - 1
        names = COLS.get(tag, [f"c{i}" for i in range(ncol)])
        for k in range(ncol):
            a = r[:, 1 + k]
            b = t[:, 1 + k]
            d = np.abs(b - a)
            mad = float(d.max()) if d.size else 0.0
            mav = float(max(np.abs(a).max(), np.abs(b).max())) if a.size else 0.0
            rel = mad / mav if mav > 0 else 0.0
            n_sig = int((d > THRESHOLD).sum())
            nm = names[k] if k < len(names) else f"c{k}"
            flag = "  <-- SIGNAL" if mad > THRESHOLD and (rel > 1e-10 or mav > 1e-12) else ""
            print(f"  s{step}  {tag:>10}  {nm:>12}  {mad:12.3e}  "
                  f"{mav:12.3e}  {rel:10.3e}  {n_sig:>7d}{flag}")
            worst = max(worst, mad)
        print()
    print("=" * 84)
    print(f"worst max|diff| across all dumps: {worst:.3e}")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None)
