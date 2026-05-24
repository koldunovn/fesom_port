#!/usr/bin/env python3
"""K9 end-to-end KPP climate validation: does C+KPP reproduce the Fortran KPP
reference (the DEFAULT CORE2 config)?

  Ckpp = kpp_2yr_dt1800     (C + KPP, dt=1800 linfs, 1958-1959)
  F    = fortran_2yr_dt1800 (Fortran + KPP, dt=1800 linfs, 1958-1959) -- the reference
  Cpp  = pp_3yr_bvfreqfix    (C + PP, for the scheme contrast: PP vs the KPP reference)

Reports area-weighted bias+RMS of annual-mean dSST/dSSS globally + by region for
1958/1959, for Ckpp-F (the headline KPP validation) and Cpp-F (how far PP is from the
KPP reference -- KPP should be closer where the schemes differ, e.g. deep-convection
high latitudes). Mesh/area/region logic mirrors clim_validation_2yr.py.
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import nereus as nr

Ckpp = "/work/ab0995/a270088/port/kpp_2yr_dt1800"
F    = "/scratch/a/a270088/fortran_2yr_dt1800"
Cpp  = "/work/ab0995/a270088/port/pp_3yr_bvfreqfix"
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
    print("="*84); print(f"YEAR {year}   dSST/dSSS (C - Fortran-KPP), area-weighted")
    print("   Ckpp = C+KPP vs F+KPP (headline) | Cpp = C+PP vs F+KPP (scheme contrast)")
    print("="*84)
    for var, unit in [("sst","degC"), ("sss","PSU")]:
        fF = ann(F, var, "", year)
        fK = ann(Ckpp, var, ".monthly", year)
        try:    fP = ann(Cpp, var, ".monthly", year); haveP = True
        except Exception: haveP = False
        wet = np.isfinite(fF) & np.isfinite(fK)
        dK = np.where(wet, fK-fF, np.nan)
        dP = np.where(wet, fP-fF, np.nan) if haveP else None
        hdr = f"  d{var.upper()} [{unit}]  {'region':18s} {'KPP bias':>9s} {'KPP rms':>8s}"
        if haveP: hdr += f" | {'PP bias':>9s} {'PP rms':>8s}"
        print(hdr)
        for nm, m in REGIONS(wet):
            if m.sum() > 5:
                line = f"  {'':9s}{nm:18s} {wmean(dK,m):+9.4f} {wrms(dK,m):8.4f}"
                if haveP: line += f" | {wmean(dP,m):+9.4f} {wrms(dP,m):8.4f}"
                print(line)
