#!/usr/bin/env python3
"""Figures: C+KPP 5-yr (dt=1800) vs Fortran KPP — integral time series + spatial maps.

  C    = kpp_5yr_dt1800     (C port, KPP, dt=1800, 1958-1962)  [.monthly suffix]
  Fd12 = fortran_5yr_d1200  (Fortran, KPP, dt=1200, 1958-1962) -- the 5-yr reference (dt caveat)
  Fd18 = fortran_2yr_dt1800 (Fortran, KPP, dt=1800, 1958-1959) -- EXACT config, 2 yr

Fig 1: monthly global-mean SST, SSS, and NH/SH sea-ice area over 5 yr, all runs overlaid.
Fig 2/3: 5-yr-mean SST / SSS maps (C, Fd12, C-Fd12) + the exact-config 2-yr C-Fd18 diff.
"""
import pathlib, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import nereus as nr

C    = "/work/ab0995/a270088/port/kpp_5yr_dt1800"      # C port, KPP, dt1800, 5yr
Fk5  = "/scratch/a/a270088/fortran_kpp_5yr_d1800"       # Fortran KPP dt1800 5yr (EXACT match; may be running)
Fd18 = "/scratch/a/a270088/fortran_2yr_dt1800"          # Fortran KPP dt1800, 2yr (exact, available)
DIAG = f"{Fd18}/fesom.mesh.diag.nc"
OUT  = pathlib.Path("/home/a/a270088/port2/fesom2_port/docs/kpp_5yr_figures"); OUT.mkdir(parents=True, exist_ok=True)
FILL = 1e29
YEARS = range(1958, 1963)

mesh = nr.fesom.load_mesh(DIAG)
lon = np.asarray(mesh["lon"].values); lat = np.asarray(mesh["lat"].values)
if lon.max() > 180: lon = np.where(lon > 180, lon - 360, lon)
ea = np.asarray(nc.Dataset(DIAG).variables["elem_area"][:])
fn = np.asarray(nc.Dataset(DIAG).variables["face_nodes"][:]); fn = fn.T if fn.shape[0]==3 else fn
fn = fn-1 if fn.min()==1 else fn
na = np.zeros(len(lat))
for k in range(3): np.add.at(na, fn[:,k], ea/3.0)
wet_na_sum = None

def monthly(path, var, suf, year):
    f = f"{path}/{var}.fesom.{year}{suf}.nc"
    a = np.asarray(nc.Dataset(f).variables[var][:])              # (12, n)
    return np.where(np.abs(a) < FILL, a, np.nan)
def gmean(field):                                                # area-weighted, per month
    m = np.isfinite(field); return np.array([np.sum(np.where(m[i],field[i]*na,0))/np.sum(np.where(m[i],na,0)) for i in range(field.shape[0])])
def ice_area(a_ice, hemi):                                       # 1e12 m^2, per month
    sel = (lat>0) if hemi=="N" else (lat<0)
    return np.array([np.nansum(a_ice[i]*na*sel)/1e12 for i in range(a_ice.shape[0])])

# ---------- time series ----------
def series(path, suf, var, years):
    out=[]
    for y in years:
        try: out.append(monthly(path, var, suf, y))
        except Exception: break
    return np.concatenate(out,0) if out else None

runs = [("C KPP dt1800 (5yr)", C, ".monthly", "C0", "-"),
        ("Fortran KPP dt1800 (5yr)", Fk5, "", "C1", "--"),
        ("Fortran KPP dt1800 (2yr ref)", Fd18, "", "k", ":")]
fig, ax = plt.subplots(3,1, figsize=(11,11), sharex=True)
for label, path, suf, col, ls in runs:
    sst = series(path, suf, "sst", YEARS); sss = series(path, suf, "sss", YEARS)
    aice= series(path, suf, "a_ice", YEARS)
    if sst is not None:
        mo = np.arange(sst.shape[0])/12.0 + 1958
        ax[0].plot(mo, gmean(sst), col, ls=ls, label=label)
    if sss is not None:
        mo = np.arange(sss.shape[0])/12.0 + 1958
        ax[1].plot(mo, gmean(sss), col, ls=ls, label=label)
    if aice is not None:
        mo = np.arange(aice.shape[0])/12.0 + 1958
        ax[2].plot(mo, ice_area(aice,"N"), col, ls=ls, label=f"{label} NH")
        ax[2].plot(mo, ice_area(aice,"S"), col, ls="-." if ls=="-" else ls, alpha=0.6)
ax[0].set_ylabel("global mean SST [°C]"); ax[0].legend(fontsize=8); ax[0].set_title("C+KPP (dt1800) vs Fortran KPP — global integrals, monthly")
ax[1].set_ylabel("global mean SSS [PSU]")
ax[2].set_ylabel("sea-ice area [1e12 m²]\n(NH solid/SH faint)"); ax[2].set_xlabel("year")
for a in ax: a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT/"timeseries_5yr.png", dpi=130); plt.close(fig)
print(f"wrote {OUT}/timeseries_5yr.png")

# ---------- spatial maps ----------
def mean_field(path, suf, var, years):
    acc=[];
    for y in years:
        try: acc.append(np.nanmean(monthly(path,var,suf,y),0))
        except Exception: break
    return np.nanmean(np.array(acc),0) if acc else None

for var, unit, dlim in [("sst","°C",0.3), ("sss","PSU",0.15)]:
    cC  = mean_field(C,".monthly",var,YEARS)            # C 5yr mean
    cF  = mean_field(Fk5,"",var,YEARS)                  # Fortran KPP dt1800 5yr mean (exact)
    cF2 = mean_field(Fd18,"",var,range(1958,1960))      # Fortran KPP dt1800 2yr mean
    cC2 = mean_field(C,".monthly",var,range(1958,1960)) # C 2yr mean
    nan  = np.full_like(cC, np.nan)
    cFok = cF if cF is not None else nan
    lo = np.nanpercentile(np.concatenate([cC[np.isfinite(cC)], cFok[np.isfinite(cFok)]]), 1)
    hi = np.nanpercentile(np.concatenate([cC[np.isfinite(cC)], cFok[np.isfinite(cFok)]]), 99)
    dC  = np.where(np.isfinite(cC)&np.isfinite(cFok), cC-cFok, np.nan)
    dC2 = np.where(np.isfinite(cC2)&np.isfinite(cF2),  cC2-cF2, np.nan)
    panels = [(f"C+KPP (5yr mean)",               cC,  "viridis", lo, hi,   f"{var} [{unit}]"),
              (f"Fortran+KPP (5yr mean)",         cFok,"viridis", lo, hi,   f"{var} [{unit}]"),
              (f"Δ = C − Fortran (5yr mean)",     dC,  "RdBu_r", -dlim, dlim, f"Δ{var} [{unit}]"),
              (f"Δ = C − Fortran (2yr mean)",     dC2, "RdBu_r", -dlim, dlim, f"Δ{var} [{unit}]")]
    fig = plt.figure(figsize=(16,9))
    for i,(ttl,d,cmap,vmn,vmx,cl) in enumerate(panels):
        ax = fig.add_subplot(2,2,i+1, projection=ccrs.Robinson()); ax.set_global()
        nr.plot(d, lon, lat, projection="rob", cmap=cmap, vmin=vmn, vmax=vmx, ax=ax,
                colorbar=True, colorbar_label=cl, title=ttl)
    fig.suptitle(f"{var.upper()}: C+KPP vs Fortran+KPP, dt=1800, 5-yr (1958-1962)", fontsize=13)
    fig.tight_layout(); fig.savefig(OUT/f"maps_{var}_5yr.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {OUT}/maps_{var}_5yr.png")

# ---------- quantitative: per-year global RMS, C vs Fortran-KPP (dt1800) ----------
print("\n=== C+KPP vs Fortran+KPP (dt1800) global SST/SSS RMS+bias by year (area-weighted) ===")
for var, unit in [("sst","°C"), ("sss","PSU")]:
    for y in YEARS:
        try:
            cc = np.nanmean(monthly(C, var, ".monthly", y), 0)
            ff = np.nanmean(monthly(Fk5, var, "", y), 0)
        except Exception:
            continue
        m = np.isfinite(cc) & np.isfinite(ff)
        rms  = np.sqrt(np.sum((cc[m]-ff[m])**2 * na[m]) / np.sum(na[m]))
        bias = np.sum((cc[m]-ff[m]) * na[m]) / np.sum(na[m])
        print(f"  {var.upper():3s} {y}: RMS={rms:.4f} {unit}  bias={bias:+.4f}")
