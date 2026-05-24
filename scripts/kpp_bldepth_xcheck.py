#!/usr/bin/env python3
"""K5 bldepth algebra isolation: confirm hbl/kbl/stable/caseA disagreements track
the ustar/Bo (forcing) INPUT disagreements, not a bldepth algebra bug.

dVsq/dbsfc/bvfreq are already bit-identical C-vs-Fortran (see kpp_dump_diff.py), so
the only differing bldepth inputs are ustar and Bo (forcing-derived). This script
merges the prestep (ustar,Bo) + bldepth (hbl,kbl,bfsfc,stable,caseA) dumps over all
ranks, partitions nodes by whether the inputs match, and reports output agreement in
the input-matching subset. If bldepth's algebra is faithful, the input-matching
subset must agree to FP precision.

Usage: kpp_bldepth_xcheck.py <ref_dir(Fortran)> <tst_dir(C)>
"""
import glob, os, re, sys
import numpy as np

FNAME_RE = re.compile(r"kpp_dump_s(\d+)_(.+)_rank(\d+)\.txt$")

def load(directory, tag):
    """gid -> np.array(ncomp) merged over ranks."""
    out = {}
    for f in glob.glob(os.path.join(directory, f"kpp_dump_s*_{tag}_rank*.txt")):
        if not FNAME_RE.search(os.path.basename(f)):
            continue
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                gid = int(p[0])
                out[gid] = np.array([float(x) for x in p[1:]])
    return out

def main():
    ref, tst = sys.argv[1], sys.argv[2]
    pre_r, pre_t = load(ref, "prestep"), load(tst, "prestep")
    bl_r,  bl_t  = load(ref, "bldepth"), load(tst, "bldepth")
    gids = sorted(set(pre_r) & set(pre_t) & set(bl_r) & set(bl_t))
    n = len(gids)
    print(f"common gids: {n}")

    ustar_r = np.array([pre_r[g][0] for g in gids]); ustar_t = np.array([pre_t[g][0] for g in gids])
    bo_r    = np.array([pre_r[g][1] for g in gids]); bo_t    = np.array([pre_t[g][1] for g in gids])
    hbl_r = np.array([bl_r[g][0] for g in gids]); hbl_t = np.array([bl_t[g][0] for g in gids])
    kbl_r = np.array([bl_r[g][1] for g in gids]); kbl_t = np.array([bl_t[g][1] for g in gids])
    stb_r = np.array([bl_r[g][3] for g in gids]); stb_t = np.array([bl_t[g][3] for g in gids])
    csA_r = np.array([bl_r[g][4] for g in gids]); csA_t = np.array([bl_t[g][4] for g in gids])

    du = np.abs(ustar_r - ustar_t)
    db = np.abs(bo_r - bo_t)
    # "input match" = ustar and Bo both agree to FP precision
    TOL_U, TOL_B = 1e-12, 1e-12
    match = (du < TOL_U) & (db < TOL_B)
    nm = int(match.sum())
    print(f"input-matching (|du|<{TOL_U} and |dBo|<{TOL_B}): {nm}  ({100*nm/n:.1f}%)")

    def report(mask, label):
        if mask.sum() == 0:
            print(f"  [{label}] (empty)"); return
        dh = np.abs(hbl_r[mask] - hbl_t[mask])
        dk = np.abs(kbl_r[mask] - kbl_t[mask])
        ds = np.abs(stb_r[mask] - stb_t[mask])
        dc = np.abs(csA_r[mask] - csA_t[mask])
        print(f"  [{label}] n={int(mask.sum())}")
        print(f"     hbl   max|d|={dh.max():.3e}  #(>1e-6)={int((dh>1e-6).sum())}")
        print(f"     kbl   max|d|={dk.max():.0f}    #(!=0)  ={int((dk>0).sum())}")
        print(f"     stable #(!=0)={int((ds>0).sum())}   caseA #(!=0)={int((dc>0).sum())}")

    report(match, "ustar&Bo match -> EXPECT hbl/kbl exact")
    report(~match, "ustar|Bo differ -> diffs expected")

    # correlation: does |dhbl| grow with input disagreement?
    dh_all = np.abs(hbl_r - hbl_t)
    big = (du > 1e-6) | (db > 1e-9)
    print(f"\nnodes with sizeable input diff (|du|>1e-6 or |dBo|>1e-9): {int(big.sum())}")
    print(f"  of those, #(|dhbl|>1m) = {int((dh_all[big]>1.0).sum())}")
    print(f"  nodes with SMALL input diff but |dhbl|>1m (unexplained?): "
          f"{int((dh_all[~big]>1.0).sum())}")

if __name__ == "__main__":
    main()
