#!/usr/bin/env python3
"""Z9 end-to-end zstar climate validation: does C+zstar+KPP reproduce the NEW
Fortran zstar reference?

  Czs  = c_zstar_2yr         (C, FESOM_ALE=zstar, KPP, dt=1800, 1958-1959)
  F    = fortran_zstar_2yr   (Fortran which_ALE='zstar', KPP — the reference)
  Flin = fortran_linfs_2yr_b (Fortran linfs rerun, same frozen binary — the
                              COORDINATE CONTRAST: C-zstar vs F-linfs must show
                              a clearly LARGER signal than C-zstar vs F-zstar,
                              so the comparison resolves the coordinate; K9 pattern)

Reports area-weighted bias+RMS of annual-mean dSST/dSSS globally + by region.
Acceptance bar (plan Z9): KPP-class agreement — RMS ~0.005/0.002, biases ~0,
non-drifting 1958→1959. Clone of kpp_climate_compare.py.
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import nereus as nr

Czs  = "/work/ab0995/a270088/port/zstar/c_zstar_2yr"
F    = "/work/ab0995/a270088/port/zstar/fortran_zstar_2yr"
Flin = "/work/ab0995/a270088/port/zstar/fortran_linfs_2yr_b"
DIAG = f"{F}/fesom.mesh.diag.nc"
FILL = 1e29

mesh = nr.fesom.load_mesh(DIAG)
lon = np.asarray(mesh["lon"].values); lat = np.asarray(mesh["lat"].values)
if lon.max() > 180: lon = np.where(lon > 180, lon - 360, lon)
ea = np.asarray(nc.Dataset(DIAG).variables["elem_area"][:])
fn = np.asarray(nc.Dataset(DIAG).variables["face_nodes"][:]); fn = fn.T if fn.shape[0]==3 else fn
fn = fn-1 if fn.min()==1 else fn
na = np.zeros(len(lat))
for k in range(3): np.add.at(na, fn[:,k], ea/3.0)

def ann(path, var, suffix, year):
    a = np.asarray(nc.Dataset(f"{path}/{var}.fesom.{year}{suffix}.nc").variables[var][:])
    a = np.where(np.abs(a) < FILL, a, np.nan); return np.nanmean(a, 0)
def wmean(d, m): return np.sum(d[m]*na[m]) / np.sum(na[m])
def wrms(d, m):  return np.sqrt(np.sum(d[m]**2*na[m]) / np.sum(na[m]))

REGIONS = lambda wet: [
    ("GLOBAL",           wet),
    ("Tropics |lat|<23", wet & (np.abs(lat) < 23)),
    ("Midlat 40-60",     wet & (np.abs(lat) >= 40) & (np.abs(lat) < 60)),
    ("N high >60",       wet & (lat > 60)),
    ("S high <-60",      wet & (lat < -60)),
]

for year in (1958, 1959):
    print("="*84); print(f"YEAR {year}   dSST/dSSS (C-zstar − Fortran), area-weighted")
    print("   ZS = C-zstar vs F-zstar (headline) | LIN = C-zstar vs F-linfs (coordinate contrast)")
    print("="*84)
    for var, unit in [("sst","degC"), ("sss","PSU")]:
        fF  = ann(F,    var, "", year)
        fC  = ann(Czs,  var, ".monthly", year)
        try:    fL = ann(Flin, var, "", year); haveL = True
        except Exception: haveL = False
        wet = np.isfinite(fF) & np.isfinite(fC)
        dZ = np.where(wet, fC-fF, np.nan)
        dL = np.where(wet, fC-fL, np.nan) if haveL else None
        hdr = f"  d{var.upper()} [{unit}]  {'region':18s} {'ZS bias':>9s} {'ZS rms':>8s}"
        if haveL: hdr += f" | {'LIN bias':>9s} {'LIN rms':>8s}"
        print(hdr)
        for name, m in REGIONS(wet):
            row = f"               {name:18s} {wmean(dZ,m):+9.4f} {wrms(dZ,m):8.4f}"
            if haveL: row += f" | {wmean(dL,m):+9.4f} {wrms(dL,m):8.4f}"
            print(row)
        print()
