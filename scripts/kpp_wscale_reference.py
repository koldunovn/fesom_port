#!/usr/bin/env python3
"""Independent reference for KPP wscale (Phase K2).

Reimplements wscale (oce_ale_mixing_kpp.F90:828-877) directly from the Fortran
— the 2D lookup-table bilinear interpolation (zehat<=zmax) and the stable
formula (zehat>zmax) — independent of the C port. Validates a wscale sweep
dump (`i j zehat ustar wm ws`) written by both codes' sweep harness.

Usage:
    python3 kpp_wscale_reference.py <kpp_wscale_rank0.txt> [...]

The lookup table is rebuilt from kpp_init_reference (already validated bit/4e-16
vs Fortran in K1), so any mismatch here isolates the wscale index/interp logic.
PASS if wm/ws max relative diff <= 1e-10.
"""
import sys
import numpy as np
import kpp_init_reference as kinit

NNI, NNJ = kinit.NNI, kinit.NNJ
ZMIN, ZMAX, UMIN = kinit.ZMIN, kinit.ZMAX, kinit.UMIN
VONK, CONC1, EPSLN = kinit.VONK, kinit.CONC1, kinit.EPSLN

_REF = kinit.reference()
WMT, WST = _REF["wmt"], _REF["wst"]
DELTAZ, DELTAU = _REF["deltaz"], _REF["deltau"]


def wscale(zehat, us):
    """Mirror of wscale(zehat, us, wm, ws)."""
    if zehat <= ZMAX:
        zdiff = zehat - ZMIN
        iz = int(zdiff / DELTAZ)        # INT (truncate toward 0)
        iz = min(iz, NNI)
        iz = max(iz, 0)
        izp1 = iz + 1
        udiff = us - UMIN
        ju = int(min(udiff / DELTAU, float(NNJ)))
        ju = max(ju, 0)
        jup1 = ju + 1
        zfrac = zdiff / DELTAZ - float(iz)
        ufrac = udiff / DELTAU - float(ju)
        fzfrac = 1.0 - zfrac
        wam = fzfrac * WMT[iz, jup1] + zfrac * WMT[izp1, jup1]
        wbm = fzfrac * WMT[iz, ju] + zfrac * WMT[izp1, ju]
        wm = (1.0 - ufrac) * wbm + ufrac * wam
        was = fzfrac * WST[iz, jup1] + zfrac * WST[izp1, jup1]
        wbs = fzfrac * WST[iz, ju] + zfrac * WST[izp1, ju]
        ws = (1.0 - ufrac) * wbs + ufrac * was
    else:
        u3 = us * us * us
        wm = VONK * us * u3 / (u3 + CONC1 * zehat + EPSLN)
        ws = wm
    return wm, ws


def main(paths):
    overall = True
    for path in paths:
        print(f"\n=== {path} ===")
        a = np.loadtxt(path, comments="#", ndmin=2)
        # columns: i j zehat ustar wm ws
        max_wm = max_ws = 0.0
        for row in a:
            zehat, us, wm_d, ws_d = row[2], row[3], row[4], row[5]
            wm_r, ws_r = wscale(zehat, us)
            for d, r, store in ((wm_d, wm_r, "wm"), (ws_d, ws_r, "ws")):
                rel = abs(d - r) / (max(abs(r), abs(d)) or 1.0)
                if store == "wm":
                    max_wm = max(max_wm, rel)
                else:
                    max_ws = max(max_ws, rel)
        print(f"  wm: max_rel={max_wm:.3e}   ws: max_rel={max_ws:.3e}")
        ok = max_wm <= 1e-10 and max_ws <= 1e-10
        print("  -> " + ("PASS" if ok else "FAIL  <-- EXCEEDS 1e-10"))
        overall &= ok
    print("\n" + ("ALL PASS" if overall else "SOME FAILED"))
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1:])
