#!/usr/bin/env python3
"""Independent reference for the KPP init lookup tables + constants (Phase K1).

Recomputes Vtc, cg, deltaz, deltau and the full wm/ws tables directly from the
Fortran formulas (oce_ale_mixing_kpp.F90:115-219) — an implementation
independent of the C port — and compares against a dumped kpp_init_rank0.txt
(written by fesom_kpp_dump_init / the Fortran gate).

Usage:
    python3 kpp_init_reference.py <kpp_init_rank0.txt> [<kpp_init_rank0.txt> ...]

Reports max|diff| for the 4 constants and the wmt/wst tables vs this reference.
PASS if all <= 1e-12 relative (lookup interpolation/libm rounding tolerance).
"""
import sys
import numpy as np

# --- constants (verified vs Fortran + CORE2 namelists) ---
NNI, NNJ = 890, 480
EPSLN, EPSILON, VONK, CONC1 = 1.0e-40, 0.1, 0.4, 5.0
ZMIN, ZMAX, UMIN, UMAX = -4.0e-7, 0.0, 0.0, 0.04
RICR, CONCV = 0.3, 1.6
CSTAR, CONAM, CONCM = 10.0, 1.257, 8.380
CONC2, ZETAM = 16.0, -0.2
CONAS, CONCS, CONC3, ZETAS = -28.86, 98.96, 16.0, -1.0


def reference():
    Vtc = CONCV * np.sqrt(0.2 / CONCS / EPSILON) / VONK**2 / RICR
    cg = CSTAR * VONK * (CONCS * VONK * EPSILON) ** (1.0 / 3.0)
    deltaz = (ZMAX - ZMIN) / float(NNI + 1)
    deltau = (UMAX - UMIN) / float(NNJ + 1)
    wmt = np.zeros((NNI + 2, NNJ + 2))
    wst = np.zeros((NNI + 2, NNJ + 2))
    for i in range(NNI + 2):
        zehat = deltaz * i + ZMIN
        for j in range(NNJ + 2):
            usta = deltau * j + UMIN
            zeta = zehat / (usta**3 + EPSLN)
            if zehat >= 0.0:
                wm = VONK * usta / (1.0 + CONC1 * zeta)
                ws = wm
            else:
                if zeta > ZETAM:
                    wm = VONK * usta * (1.0 - CONC2 * zeta) ** 0.25
                else:
                    wm = VONK * (CONAM * usta**3 - CONCM * zehat) ** (1.0 / 3.0)
                if zeta > ZETAS:
                    ws = VONK * usta * (1.0 - CONC3 * zeta) ** 0.5
                else:
                    ws = VONK * (CONAS * usta**3 - CONCS * zehat) ** (1.0 / 3.0)
            wmt[i, j] = wm
            wst[i, j] = ws
    return {"Vtc": Vtc, "cg": cg, "deltaz": deltaz, "deltau": deltau,
            "wmt": wmt, "wst": wst}


def load_dump(path):
    consts = {}
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                p = line[1:].split()
                if len(p) == 2 and p[0] in ("Vtc", "cg", "deltaz", "deltau"):
                    consts[p[0]] = float(p[1])
                continue
            p = line.split()
            if len(p) == 4:
                rows.append((int(p[0]), int(p[1]), float(p[2]), float(p[3])))
    a = np.array(rows)
    wmt = np.zeros((NNI + 2, NNJ + 2))
    wst = np.zeros((NNI + 2, NNJ + 2))
    ii = a[:, 0].astype(int)
    jj = a[:, 1].astype(int)
    wmt[ii, jj] = a[:, 2]
    wst[ii, jj] = a[:, 3]
    return consts, wmt, wst


def relmax(a, b):
    d = np.abs(a - b)
    scale = np.maximum(np.abs(a), np.abs(b))
    scale[scale == 0] = 1.0
    return float(d.max()), float((d / scale).max())


def main(paths):
    ref = reference()
    print(f"reference constants: Vtc={ref['Vtc']:.17g} cg={ref['cg']:.17g} "
          f"deltaz={ref['deltaz']:.17g} deltau={ref['deltau']:.17g}")
    overall_ok = True
    for path in paths:
        print(f"\n=== {path} ===")
        consts, wmt, wst = load_dump(path)
        ok = True
        for nm in ("Vtc", "cg", "deltaz", "deltau"):
            if nm not in consts:
                print(f"  {nm}: MISSING in dump"); ok = False; continue
            ad = abs(consts[nm] - ref[nm])
            rd = ad / (abs(ref[nm]) or 1.0)
            flag = "" if rd <= 1e-12 else "  <-- EXCEEDS 1e-12"
            print(f"  {nm:7s} dump={consts[nm]:.17g} ref={ref[nm]:.17g} "
                  f"|Δ|={ad:.3e} rel={rd:.3e}{flag}")
            ok &= rd <= 1e-12
        for nm, tab in (("wmt", wmt), ("wst", wst)):
            mad, mrel = relmax(tab, ref[nm])
            flag = "" if mrel <= 1e-12 else "  <-- EXCEEDS 1e-12"
            print(f"  {nm}: max|Δ|={mad:.3e} max_rel={mrel:.3e}{flag}")
            ok &= mrel <= 1e-12
        print("  -> " + ("PASS" if ok else "FAIL"))
        overall_ok &= ok
    print("\n" + ("ALL PASS" if overall_ok else "SOME FAILED"))
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1:])
