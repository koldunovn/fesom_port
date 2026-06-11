#!/usr/bin/env python3
"""T5 end-to-end TKE climate validation: does C+linfs+TKE reproduce the
Fortran linfs+cvmix_TKE reference?

  C    = c_tke_2yr           (C, FESOM_MIX_SCHEME=TKE, linfs, dt=1800, 1958-59)
  F    = fortran_linfs_tke   (Fortran mix_scheme='cvmix_TKE' — the reference,
                              one-knob copy of work_linfs_2yr_b, frozen binary)
  Fkpp = fortran_linfs_2yr_b (Fortran linfs+KPP, same frozen binary — the
                              SCHEME CONTRAST: C-TKE vs F-KPP must show a
                              clearly LARGER signal than C-TKE vs F-TKE, so
                              the comparison resolves the mixing scheme; the
                              Z9/K9 pattern)

Reports area-weighted bias+RMS of annual-mean dSST/dSSS globally + by region.
Acceptance bar (plan T5): KPP-class agreement — RMS ~0.005/0.002, biases ~0,
non-drifting 1958→1959. Clone of zstar_climate_compare.py.
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import nereus as nr

C    = "/work/ab0995/a270088/port/tke/c_tke_2yr"
F    = "/work/ab0995/a270088/port/tke/fortran_linfs_tke"
Fkpp = "/work/ab0995/a270088/port/zstar/fortran_linfs_2yr_b"
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
    print("="*84); print(f"YEAR {year}   dSST/dSSS (C-TKE − Fortran), area-weighted")
    print("   TKE = C-TKE vs F-TKE (headline) | KPP = C-TKE vs F-KPP (scheme contrast)")
    print("="*84)
    for var, unit in [("sst","degC"), ("sss","PSU")]:
        fF  = ann(F,    var, "", year)
        fC  = ann(C,    var, ".monthly", year)
        try:    fK = ann(Fkpp, var, "", year); haveK = True
        except Exception: haveK = False
        wet = np.isfinite(fF) & np.isfinite(fC)
        dT = np.where(wet, fC-fF, np.nan)
        dK = np.where(wet, fC-fK, np.nan) if haveK else None
        hdr = f"  d{var.upper()} [{unit}]  {'region':18s} {'TKE bias':>9s} {'TKE rms':>8s}"
        if haveK: hdr += f" | {'KPP bias':>9s} {'KPP rms':>8s}"
        print(hdr)
        for name, m in REGIONS(wet):
            row = f"               {name:18s} {wmean(dT,m):+9.4f} {wrms(dT,m):8.4f}"
            if haveK: row += f" | {wmean(dK,m):+9.4f} {wrms(dK,m):8.4f}"
            print(row)
        print()
