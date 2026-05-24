#!/usr/bin/env python3
"""DEFINITIVE 2-yr climate validation of the C port at HEAD 6d65bd3 (first full run with
the cyclic-seam IC fix) vs the Fortran+PP dt=1800 reference.

  C  = clim_2yr_dt1800   (dt=1800, PP, linfs, 1958-1959, all fixes incl. 6d65bd3)
  F  = fortran_pp_2yr    (dt=1800, PP, linfs, 1958-1959)
  B  = monthfix_1yr      ('before' the cyclic-seam fix — for the W Med 3-way)

Reports area-weighted MEAN (bias) AND RMS of annual-mean ΔSST/ΔSSS, globally + by region,
for 1958 and 1959 (drift check). Special focus: the western Mediterranean SSS, which the
cyclic-lon IC bug corrupted (flat ~37.35 vs PHC ~37.7) — confirm C-after now tracks Fortran.

Run: PYTHONPATH=/home/a/a270088/PYTHON \
  /work/ab0995/a270088/mambaforge/envs/nereus/bin/python scripts/clim_validation_2yr.py
"""
import pathlib, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import nereus as nr

C = "/work/ab0995/a270088/port/clim_2yr_dt1800"     # C after all fixes (HEAD 6d65bd3)
F = "/scratch/a/a270088/fortran_pp_2yr"             # Fortran+PP dt=1800 reference
B = "/work/ab0995/a270088/port/monthfix_1yr"        # C before cyclic-seam fix
DIAG = f"{F}/fesom.mesh.diag.nc"
OUT = pathlib.Path("/home/a/a270088/port2/fesom2_port/docs/validation_2yr_d1800")
OUT.mkdir(parents=True, exist_ok=True)
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
    a = np.asarray(nc.Dataset(f"{path}/{var}.fesom.{year}{suffix}.nc").variables[var][:])  # (12,n)
    a = np.where(np.abs(a) < FILL, a, np.nan); return np.nanmean(a, 0)

def wmean(d, m): return np.sum(d[m]*na[m]) / np.sum(na[m])
def wrms(d, m):  return np.sqrt(np.sum(d[m]**2*na[m]) / np.sum(na[m]))

REGIONS = lambda wet: [
    ("GLOBAL",            wet),
    ("Tropics |lat|<23",  wet & (np.abs(lat) < 23)),
    ("Subtropics 23-40",  wet & (np.abs(lat) >= 23) & (np.abs(lat) < 40)),
    ("Midlat 40-60",      wet & (np.abs(lat) >= 40) & (np.abs(lat) < 60)),
    ("N high >60",        wet & (lat > 60)),
    ("S high <-60",       wet & (lat < -60)),
    ("W Med",             wet & (lat>35.5)&(lat<44.5)&(lon>-6)&(lon<10)),
    ("seam-Med |lon|<2",  wet & (lat>35.5)&(lat<44.5)&(np.abs(lon)<2)),
]

print("="*78)
print("C (clim_2yr_dt1800, HEAD 6d65bd3)  −  Fortran (fortran_pp_2yr), dt=1800 PP linfs")
print("="*78)
for year in (1958, 1959):
    print(f"\n################################  YEAR {year}  ################################")
    for var, unit in [("sst","degC"), ("sss","PSU")]:
        cF = ann(F, var, "", year); cC = ann(C, var, ".monthly", year)
        wet = np.isfinite(cF) & np.isfinite(cC)
        d = np.where(wet, cC - cF, np.nan)
        print(f"\n  --- Δ{var.upper()} [{unit}] (C − F), area-weighted ---")
        print(f"  {'region':20s} {'mean(bias)':>11s} {'rms':>8s} {'|max|':>8s}   n")
        for nm, m in REGIONS(wet):
            if m.sum() > 5:
                print(f"  {nm:20s} {wmean(d,m):+11.4f} {wrms(d,m):8.4f} "
                      f"{np.nanmax(np.abs(d[m])):8.3f}  {int(m.sum()):>6d}")

# ---- western-Mediterranean 3-way: validates the cyclic-seam IC fix (6d65bd3) ----
print("\n" + "="*78)
print("WESTERN-MED SSS 3-WAY (validates cyclic-seam fix 6d65bd3), annual 1958")
print("  F=Fortran  B=C-before(monthfix_1yr)  A=C-after(clim_2yr, 6d65bd3)")
print("="*78)
sF = ann(F, "sss", "", 1958); sB = ann(B, "sss", ".monthly", 1958); sA = ann(C, "sss", ".monthly", 1958)
wet = np.isfinite(sF) & np.isfinite(sB) & np.isfinite(sA)
for nm, m in [("W Med", wet&(lat>35.5)&(lat<44.5)&(lon>-6)&(lon<10)),
              ("seam-Med |lon|<2", wet&(lat>35.5)&(lat<44.5)&(np.abs(lon)<2))]:
    if m.sum() > 3:
        print(f"\n  {nm} (n={int(m.sum())}):")
        print(f"    area-wtd SSS:   F={wmean(sF,m):.3f}   before={wmean(sB,m):.3f}   after={wmean(sA,m):.3f}")
        dB = np.where(wet, sB-sF, np.nan); dA = np.where(wet, sA-sF, np.nan)
        print(f"    ΔSSS vs F rms:  before={wrms(dB,m):.3f} -> after={wrms(dA,m):.3f}"
              f"   (mean bias {wmean(dB,m):+.3f} -> {wmean(dA,m):+.3f})")

# ---- maps: global ΔSST/ΔSSS + W Med zoom ----
fig = plt.figure(figsize=(16, 11))
for i,(var,vlim,unit) in enumerate([("sst",1.0,"degC"),("sss",0.5,"PSU")]):
    cF = ann(F,var,"",1958); cC = ann(C,var,".monthly",1958)
    wet = np.isfinite(cF)&np.isfinite(cC); d = np.where(wet, cC-cF, np.nan)
    ax = fig.add_subplot(2,2,i+1, projection=ccrs.Robinson()); ax.set_global()
    nr.plot(d, lon, lat, projection="rob", cmap="RdBu_r", vmin=-vlim, vmax=vlim, ax=ax,
            colorbar=True, colorbar_label=f"Δ{var} [{unit}]", title=f"Δ{var.upper()} = C − Fortran (1958)")
# W Med zoom for SSS: before vs after diff to Fortran
for j,(tag,run,suf) in enumerate([("before (monthfix_1yr)",B,".monthly"),("after (clim_2yr, 6d65bd3)",C,".monthly")]):
    s = ann(run,"sss",suf,1958); sf = ann(F,"sss","",1958)
    wet = np.isfinite(s)&np.isfinite(sf); d = np.where(wet, s-sf, np.nan)
    ax = fig.add_subplot(2,2,3+j, projection=ccrs.PlateCarree())
    ax.set_extent([-7,12,34,46], ccrs.PlateCarree())
    sc = ax.scatter(lon, lat, c=d, s=8, cmap="RdBu_r", vmin=-0.5, vmax=0.5, transform=ccrs.PlateCarree())
    ax.coastlines(resolution="50m"); ax.set_title(f"W Med ΔSSS: {tag} − Fortran")
    plt.colorbar(sc, ax=ax, label="ΔSSS [PSU]", shrink=0.7)
fig.suptitle("C port (HEAD 6d65bd3) vs Fortran+PP dt=1800 — annual 1958", fontsize=13)
fig.tight_layout(); fp = OUT/"clim_validation_2yr.png"; fig.savefig(fp, dpi=120, bbox_inches="tight")
print("\nwrote", fp)
