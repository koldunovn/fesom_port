#!/usr/bin/env python3
"""K6 blmix algebra isolation: where ALL blmix inputs that differ C-vs-Fortran
(the bldepth outputs hbl/bfsfc/stable/caseA/kbl + ustar) match, blmc/dkm1/ghats
MUST match — every other blmix input (interior viscA/diffK, mesh hnode/Z/zbar, the
wscale table, cg) is already bit-identical. So partition nodes by input agreement
and report blmc/dkm1/ghats agreement in the matching subset. A clean match there
proves the algebra; the live-run diffs elsewhere are forcing-driven (bfsfc differs
at ~98% of nodes, project_forcing_step1_diff).

Usage: kpp_blmix_xcheck.py <ref_dir(Fortran)> <tst_dir(C)>
"""
import glob, os, re, sys
import numpy as np

FN = re.compile(r"kpp_dump_s(\d+)_(.+)_rank(\d+)\.txt$")

def load(directory, tag):
    out = {}
    for f in glob.glob(os.path.join(directory, f"kpp_dump_s*_{tag}_rank*.txt")):
        if not FN.search(os.path.basename(f)):
            continue
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                out[int(p[0])] = np.array([float(x) for x in p[1:]])
    return out

def col(d, gids, c):
    return np.array([d[g][c] for g in gids])

def colmax(d, gids):
    """max over all columns per gid -> (n,) of per-node max|.| is done by caller."""
    return np.array([d[g] for g in gids])

def main():
    ref, tst = sys.argv[1], sys.argv[2]
    bd_r, bd_t = load(ref, "bldepth"), load(tst, "bldepth")
    pr_r, pr_t = load(ref, "prestep"), load(tst, "prestep")
    tags = ["blmc_m", "blmc_t", "blmc_s", "dkm1", "bl_ghats"]
    out_r = {t: load(ref, t) for t in tags}
    out_t = {t: load(tst, t) for t in tags}

    gids = set(bd_r) & set(bd_t) & set(pr_r) & set(pr_t)
    for t in tags:
        gids &= set(out_r[t]) & set(out_t[t])
    gids = sorted(gids)
    n = len(gids)
    print(f"common gids: {n}")

    # bldepth cols: 0 hbl, 1 kbl, 2 bfsfc, 3 stable, 4 caseA ; prestep: 0 ustar
    dhbl = np.abs(col(bd_r, gids, 0) - col(bd_t, gids, 0))
    dkbl = np.abs(col(bd_r, gids, 1) - col(bd_t, gids, 1))
    dbfs = np.abs(col(bd_r, gids, 2) - col(bd_t, gids, 2))
    dstb = np.abs(col(bd_r, gids, 3) - col(bd_t, gids, 3))
    dcsA = np.abs(col(bd_r, gids, 4) - col(bd_t, gids, 4))
    dust = np.abs(col(pr_r, gids, 0) - col(pr_t, gids, 0))

    match = (dhbl < 1e-9) & (dkbl < 0.5) & (dbfs < 1e-12) & \
            (dstb < 0.5) & (dcsA < 0.5) & (dust < 1e-12)
    nm = int(match.sum())
    print(f"input-matching nodes (hbl/bfsfc/stable/caseA/kbl/ustar all agree): "
          f"{nm}  ({100*nm/n:.2f}%)")
    if nm == 0:
        print("  (no fully-matching nodes — loosen, or rely on the logical argument)")

    midx = np.where(match)[0]
    for t in tags:
        Ar = np.array([out_r[t][gids[i]] for i in midx]) if nm else np.zeros((0,1))
        At = np.array([out_t[t][gids[i]] for i in midx]) if nm else np.zeros((0,1))
        if nm:
            d = np.abs(Ar - At)
            denom = np.maximum(np.abs(Ar), np.abs(At))
            rel = np.where(denom > 0, d / np.maximum(denom, 1e-30), 0.0)
            print(f"  [{t:9s}] matched-subset  max|d|={d.max():.3e}  max_rel={rel.max():.3e}  "
                  f"#(|d|>1e-9)={int((d>1e-9).sum())}")

    # sanity: how widespread are the input diffs (why the live-run blmc diffs are wide)
    print(f"\ninput-diff prevalence: dbfsfc>1e-12 at {int((dbfs>1e-12).sum())}/{n} "
          f"({100*(dbfs>1e-12).sum()/n:.0f}%); dustar>1e-12 at {int((dust>1e-12).sum())}/{n} "
          f"({100*(dust>1e-12).sum()/n:.0f}%); dstable!=0 at {int((dstb>0.5).sum())}")

if __name__ == "__main__":
    main()
