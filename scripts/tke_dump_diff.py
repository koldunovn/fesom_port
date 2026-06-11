#!/usr/bin/env python3
"""Diff TKE dump files between the C port and the Fortran reference.

Usage: tke_dump_diff.py <ref_dir> <tst_dir> [tag_substr] [--cols]

Clone of scripts/ale_dump_diff.py for the TKE-port dual instrumentation.
Both codes write, gated by FESOM_TKE_DUMP_DIR (TKE calls 1..FESOM_TKE_DUMP_STEPS,
default 3):

    tke_dump_s<step>_<tag>_rank<R>.txt
    # step=.. tag=.. rank=R N=.. ncomp=C   (header)
    <gid> <v0> <v1> ... <v(C-1)>           (one line per owned entity, 1-based gid)

Node tags (ncomp=nl unless noted): normstress(1) vshear2 bvfreq2 dztrr tkeold
(inputs at the cvmix call), tke tkeav tkekv (post-loop outputs), tbpr tspr tdif
tdis twin tiwf tbck ttot lmix pr (budget diagnostics), kv (final wired Kv).
Element tag: av (final wired Av, nl cols).

Default output: ONE summary line per (step, tag) — max|diff| over all columns,
the worst column index, max|val|, relative diff, #entries above THRESHOLD.
--cols prints the per-column detail for the selected tags.
"""
import glob
import os
import re
import sys

import numpy as np

THRESHOLD = 1e-9   # below this, FP-summation noise

FNAME_RE = re.compile(r"tke_dump_s(\d+)_(.+)_rank(\d+)\.txt$")


def discover(directory):
    out = {}
    for f in glob.glob(os.path.join(directory, "tke_dump_s*_rank*.txt")):
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
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx]


def main(ref_dir, tst_dir, tag_filter=None, per_col=False):
    print(f"Ref (Fortran): {ref_dir}")
    print(f"Test (C):      {tst_dir}")
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
    print(f"{'step':>4}  {'tag':>11}  {'worst-col':>9}  {'max|diff|':>12}  "
          f"{'max|val|':>12}  {'rel':>10}  {'#sig':>8}")
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
        a = r[:, 1:1 + ncol]
        b = t[:, 1:1 + ncol]
        d = np.abs(b - a)
        col_mad = d.max(axis=0) if d.size else np.zeros(ncol)
        kworst = int(col_mad.argmax()) if ncol else 0
        mad = float(col_mad.max()) if d.size else 0.0
        mav = float(max(np.abs(a).max(), np.abs(b).max())) if a.size else 0.0
        rel = mad / mav if mav > 0 else 0.0
        n_sig = int((d > THRESHOLD).sum())
        flag = "  <-- SIGNAL" if mad > THRESHOLD and (rel > 1e-10 or mav > 1e-12) else ""
        print(f"  s{step}  {tag:>11}  c{kworst:<8d}  {mad:12.3e}  "
              f"{mav:12.3e}  {rel:10.3e}  {n_sig:>8d}{flag}")
        if per_col:
            for k in range(ncol):
                ka = a[:, k]; kb = b[:, k]; kd = d[:, k]
                kmad = float(kd.max()) if kd.size else 0.0
                kmav = float(max(np.abs(ka).max(), np.abs(kb).max())) if ka.size else 0.0
                krel = kmad / kmav if kmav > 0 else 0.0
                ksig = int((kd > THRESHOLD).sum())
                if kmad > 0:
                    print(f"      {'':>11}  c{k:<8d}  {kmad:12.3e}  "
                          f"{kmav:12.3e}  {krel:10.3e}  {ksig:>8d}")
        worst = max(worst, mad)
    print("=" * 84)
    print(f"worst max|diff| across all dumps: {worst:.3e}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--cols"]
    per_col = "--cols" in sys.argv[1:]
    if len(args) not in (2, 3):
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(args[0], args[1], args[2] if len(args) == 3 else None, per_col)
