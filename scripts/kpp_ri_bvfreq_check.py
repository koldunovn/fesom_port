#!/usr/bin/env python3
"""K3 ri_iwmix correctness check via bvfreq sign (step 1, uv=0 → shear=0).

At step 1 from rest, ri_iwmix reduces to a deterministic binary of sign(bvfreq):
  viscA = A_ver (=1e-4)            if bvfreq > 0   (frit=0, stable interior)
        = visc_sh_limit+A_ver (=5.1e-3) if bvfreq <= 0 (frit=1, static instability)
(diffK = K_ver / diff_sh_limit+K_ver analogously). This script confirms:
  (1) each code is self-consistent (viscA follows its own bvfreq sign), and
  (2) every C-vs-Fortran viscA disagreement is exactly a bvfreq sign flip
      (cross-code IC noise at weak stratification), NOT an algebra bug.

Usage: kpp_ri_bvfreq_check.py <fortran_dump_dir> <c_dump_dir>
"""
import glob
import sys
import numpy as np

A_VER, K_VER = 1.0e-4, 1.0e-5
VISC_SH, DIFF_SH = 5.0e-3, 5.0e-3
VISC_HI = VISC_SH + A_VER     # 5.1e-3
DIFF_HI = DIFF_SH + K_VER     # 5.01e-3


def load(directory, tag):
    files = sorted(glob.glob(f"{directory}/kpp_dump_s1_{tag}_rank*.txt"))
    rows = [np.loadtxt(f, comments="#", ndmin=2) for f in files]
    rows = [r for r in rows if r.size]
    a = np.concatenate(rows)
    a = a[np.argsort(a[:, 0].astype(np.int64))]
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx]


def interior_mask(viscA):
    """True for genuinely interior levels (exclude per-node surface/bottom edge
    copies + below-bottom zeros). ri_iwmix sets viscA(nzmin)=viscA(nzmin+1) and
    viscA(nzmax)=viscA(nzmax-1) as edge copies, so those follow the ADJACENT
    level's bvfreq, not the local one — they must not be sign-checked locally.
    The valid (nonzero) column block per node is [nzmin, nzmax]; interior is
    nzmin+1 .. nzmax-1."""
    v = viscA[:, 1:]
    nonzero = v != 0.0
    m = np.zeros_like(nonzero)
    for r in range(v.shape[0]):
        idx = np.nonzero(nonzero[r])[0]
        if idx.size >= 3:
            m[r, idx[1]:idx[-1]] = True   # drop first (nzmin) and last (nzmax)
    return m


def self_consistency(name, viscA, bvfreq, hi, lo, mask):
    v = viscA[:, 1:]
    b = bvfreq[:, 1:]
    expect = np.where(b > 0.0, lo, hi)
    bad = mask & ~np.isclose(v, expect, rtol=0, atol=1e-12)
    print(f"  {name}: interior self-consistency violations = {int(bad.sum())} "
          f"(of {int(mask.sum())} interior entries)")
    return int(bad.sum())


def main(fdir, cdir):
    fv, cv = load(fdir, "ri_viscA"), load(cdir, "ri_viscA")
    fb, cb = load(fdir, "ri_bvfreq"), load(cdir, "ri_bvfreq")
    g = np.intersect1d(fv[:, 0], cv[:, 0])
    def sel(a):
        a = a[np.isin(a[:, 0], g)]; return a[np.argsort(a[:, 0])]
    fv, cv, fb, cb = sel(fv), sel(cv), sel(fb), sel(cb)

    mf, mc = interior_mask(fv), interior_mask(cv)

    print("=== self-consistency (each code: interior viscA follows its bvfreq sign) ===")
    nbad = self_consistency("Fortran", fv, fb, VISC_HI, A_VER, mf)
    nbad += self_consistency("C      ", cv, cb, VISC_HI, A_VER, mc)

    print("\n=== cross-code: viscA disagreement vs bvfreq sign flip (interior) ===")
    v_f, v_c = fv[:, 1:], cv[:, 1:]
    b_f, b_c = fb[:, 1:], cb[:, 1:]
    valid = mf & mc   # interior in both
    visc_disagree = valid & ~np.isclose(v_f, v_c, rtol=0, atol=1e-12)
    sign_flip = valid & ((b_f > 0.0) != (b_c > 0.0))
    # disagreements NOT explained by a bvfreq sign flip = real algebra bugs
    unexplained = visc_disagree & ~sign_flip
    print(f"  interior entries compared : {int(valid.sum())}")
    print(f"  bvfreq sign flips (C vs F): {int(sign_flip.sum())}  "
          f"({100.0*sign_flip.sum()/max(valid.sum(),1):.2f}%)")
    print(f"  viscA disagreements       : {int(visc_disagree.sum())}")
    print(f"  UNEXPLAINED (disagree but bvfreq sign agrees): {int(unexplained.sum())}")
    ok = (nbad == 0) and (int(unexplained.sum()) == 0)
    print("\n" + ("PASS: ri_iwmix algebra correct; all disagreements are bvfreq sign flips."
                  if ok else "FAIL: unexplained disagreements or self-consistency violations."))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr); sys.exit(2)
    main(sys.argv[1], sys.argv[2])
