#!/usr/bin/env python3
"""PP re-baseline check after the K5 bvfreq N2 fix (commit 5acabea): does the
fixed C-PP still match Fortran-PP at least as well as the pre-fix C-PP did?

  Cnew = pp_3yr_bvfreqfix  (dt=1800 PP linfs, 1958-1960, WITH the N2smth_h + bottom-pad fix)
  Cold = clim_2yr_dt1800   (the validated pre-fix PP, 1958-1959; the 0.040/0.016 baseline)
  F    = fortran_pp_2yr    (dt=1800 PP linfs, 1958-1959; HAS N2smth_h, the default)

Reports area-weighted bias+RMS of annual-mean dSST/dSSS (C-F) globally + by region for
1958/1959, for BOTH Cnew and Cold, so the re-baseline shift is explicit. Year 1960 of Cnew
is a stability/drift check (no Fortran reference) — reports the global-mean SST/SSS drift
yr3-vs-yr2. Mirrors clim_validation_2yr.py's mesh/area/region logic.
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import nereus as nr

Cnew = "/work/ab0995/a270088/port/pp_3yr_bvfreqfix"
Cold = "/work/ab0995/a270088/port/clim_2yr_dt1800"
F    = "/scratch/a/a270088/fortran_pp_2yr"
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
    print("="*82); print(f"YEAR {year}   dSST/dSSS (C - Fortran), area-weighted   [Cnew=bvfreq-fix, Cold=pre-fix]")
    print("="*82)
    for var, unit in [("sst","degC"), ("sss","PSU")]:
        fF = ann(F, var, "", year)
        fN = ann(Cnew, var, ".monthly", year); fO = ann(Cold, var, ".monthly", year)
        wet = np.isfinite(fF) & np.isfinite(fN) & np.isfinite(fO)
        dN = np.where(wet, fN-fF, np.nan); dO = np.where(wet, fO-fF, np.nan)
        print(f"  d{var.upper()} [{unit}]  {'region':18s} {'new bias':>9s} {'new rms':>8s} | {'old bias':>9s} {'old rms':>8s}")
        for nm, m in REGIONS(wet):
            if m.sum() > 5:
                print(f"  {'':9s}{nm:18s} {wmean(dN,m):+9.4f} {wrms(dN,m):8.4f} | "
                      f"{wmean(dO,m):+9.4f} {wrms(dO,m):8.4f}")

# year-3 stability/drift (no Fortran ref): global-mean SST/SSS yr-over-yr
print("="*82); print("Cnew multi-year global-mean drift (stability check; no Fortran ref for 1960)")
print("="*82)
for var, unit in [("sst","degC"), ("sss","PSU")]:
    gm = []
    for y in (1958, 1959, 1960):
        a = ann(Cnew, var, ".monthly", y); m = np.isfinite(a)
        gm.append(wmean(np.where(m,a,np.nan), m))
    print(f"  {var.upper():4s}[{unit}] global mean: 1958={gm[0]:.4f}  1959={gm[1]:.4f}  1960={gm[2]:.4f}"
          f"   (yr3-yr2 drift {gm[2]-gm[1]:+.4f})")
