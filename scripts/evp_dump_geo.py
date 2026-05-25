#!/usr/bin/env python3
"""Latitude-band breakdown of an EVP dump point's C-vs-Fortran |diff|.

Usage: evp_dump_geo.py <ref_dir> <tst_dir> <step> <point> <cls>
  e.g. evp_dump_geo.py .../evp_fdump/dump .../evp_cdump/dump 1 F node

Reports, per column of that dump point, the fraction of nodes/elems with
|diff|>1e-9 and max|diff| split by latitude band (gid->lat via the mesh diag).
Mirrors scripts/kpp_forcing_geo.py for the EVP dump format. For elem points the
lat is the mean of the 3 element-vertex lats.
"""
import glob
import os
import sys

import numpy as np
import netCDF4 as nc

COLS = {
    ("Q",   "node"): ["a_ice", "m_ice", "m_snow", "elevation", "srfoce_temp"],
    ("U0",  "node"): ["u_ice", "v_ice"],
    ("P",   "node"): ["inv_mass", "inv_areamass", "rhs_a", "rhs_m"],
    ("P",   "elem"): ["ice_strength"],
    ("F",   "node"): ["stress_atmice_x", "stress_atmice_y", "u_w", "v_w"],
    ("1",   "elem"): ["sigma11", "sigma12", "sigma22"],
    ("2",   "node"): ["u_rhs", "v_rhs"],
    ("3",   "node"): ["u_ice", "v_ice"],
    ("4",   "node"): ["u_ice", "v_ice"],
    ("END", "node"): ["u_ice", "v_ice"],
}
BANDS = [
    ("Antarctic <-60", lambda la: la < -60),
    ("S mid -60..-40", lambda la: (la >= -60) & (la < -40)),
    ("S sub -40..-23", lambda la: (la >= -40) & (la < -23)),
    ("Tropics |lat|<23", lambda la: np.abs(la) < 23),
    ("N sub 23..40", lambda la: (la >= 23) & (la < 40)),
    ("N mid 40..60", lambda la: (la >= 40) & (la < 60)),
    ("Arctic >60", lambda la: la > 60),
]


def load(directory, step, point, cls):
    out = {}
    for f in glob.glob(f"{directory}/evp_dump_s{step}_{point}_{cls}_rank*.txt"):
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            out[int(p[0])] = np.array([float(x) for x in p[1:]])
    return out


def mesh_latlon():
    for d in ("/scratch/a/a270088/fortran_2yr_dt1800/fesom.mesh.diag.nc",
              "/scratch/a/a270088/fortran_pp_2yr/fesom.mesh.diag.nc"):
        if os.path.exists(d):
            ds = nc.Dataset(d)
            lat = np.asarray(ds.variables["lat"][:])
            lon = np.asarray(ds.variables["lon"][:])
            elem = None
            for nm in ("elem", "elements", "face_nodes"):
                if nm in ds.variables:
                    elem = np.asarray(ds.variables[nm][:])
                    break
            if lat.max() < 3.2:
                lat = np.degrees(lat)
                lon = np.degrees(lon)
            return lat, lon, elem
    raise SystemExit("no mesh.diag.nc found")


def main(ref_dir, tst_dir, step, point, cls):
    cols = COLS[(point, cls)]
    R = load(ref_dir, step, point, cls)
    T = load(tst_dir, step, point, cls)
    g = sorted(set(R) & set(T))
    if not g:
        raise SystemExit(f"no common ids for s{step} {point} {cls}")
    lat, lon, elem = mesh_latlon()
    gi = np.array(g) - 1
    if cls == "elem":
        if elem is None:
            raise SystemExit("elem connectivity not in mesh diag; cannot map elem lat")
        ev = elem[:, gi] if elem.shape[0] == 3 else elem[gi, :].T  # (3, nelem_sel)
        la = np.nanmean(lat[(ev - 1).astype(int)], axis=0)
    else:
        la = lat[gi]
    print(f"ref={ref_dir}\ntst={tst_dir}\nstep {step} point {point} {cls}  n={len(g)}")
    for k, name in enumerate(cols):
        a = np.array([R[x][k] for x in g])
        b = np.array([T[x][k] for x in g])
        d = np.abs(a - b)
        nsig = int((d > 1e-9).sum())
        print(f"\n--- {name}: max|d|={d.max():.3e} max|val|={max(abs(a).max(),abs(b).max()):.3e} "
              f"#(|d|>1e-9)={nsig} ({100*nsig/len(g):.1f}%) ---")
        for nm, fn in BANDS:
            m = fn(la)
            if m.sum() > 0:
                frac = 100 * (d[m] > 1e-9).sum() / m.sum()
                print(f"  {nm:18s} n={int(m.sum()):6d}  frac(|d|>1e-9)={frac:5.1f}%  max|d|={d[m].max():.3e}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
