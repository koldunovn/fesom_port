#!/usr/bin/env python3
"""Comparison plots: C-port vs FESOM-Fortran, one model year (monthly means).

Produces:
  ocean_seasonal.png   — temp/salt/ssh/sst/sss area-weighted seasonal cycle
  seaice_seasonal.png  — ice area / ice volume / snow volume / mean thickness
  scatter_agreement.png— C-vs-Fortran scatter for surface T,S at mid + end year
  maps_sst_diff.png    — spatial (C-Fortran) SST difference, Jan / Jul

Usage:
  python3 plot_c_vs_fortran.py <c_dir> <fortran_dir> <mesh_diag.nc> <out_dir> [--year 1958]
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

FILL = 1e30
MONTHS = np.arange(1, 13)
MNAMES = ["J","F","M","A","M","J","J","A","S","O","N","D"]


def valid(a):
    a = np.asarray(a)
    return (a != 0) & (np.abs(a) < FILL)


def masked_mean(a):
    a = np.asarray(a); m = valid(a)
    return a[m].mean() if m.sum() else np.nan


def area_wmean(x, area):
    """Area-weighted mean of a 1D node field over ocean nodes."""
    x = np.asarray(x); m = valid(x)
    return np.sum(x[m] * area[m]) / np.sum(area[m]) if m.sum() else np.nan


def load_year(path, var, n_months=12):
    ds = nc.Dataset(path)
    arr = np.asarray(ds[var][...])
    return arr[:n_months]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("c_dir", type=pathlib.Path)
    ap.add_argument("f_dir", type=pathlib.Path)
    ap.add_argument("mesh_diag", type=pathlib.Path)
    ap.add_argument("out_dir", type=pathlib.Path)
    ap.add_argument("--year", type=int, default=1958)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    yr = args.year

    # Surface node area + geographic coords from mesh.diag
    md = nc.Dataset(args.mesh_diag)
    area = np.asarray(md["nod_area"][0, :])           # [nod2]
    lon = np.degrees(np.asarray(md["lon"][:])) if "lon" in md.variables else None
    lat = np.degrees(np.asarray(md["lat"][:])) if "lat" in md.variables else None

    def cpath(v):  return args.c_dir / f"{v}.fesom.{yr}.monthly.nc"
    def fpath(v):  return args.f_dir / f"{v}.fesom.{yr}.nc"

    C = "#1f77b4"; F = "#d62728"   # blue C, red Fortran

    # ---- Figure 1: ocean seasonal cycle -----------------------------
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(f"Ocean seasonal cycle {yr}: C-port (dt=500) vs Fortran (dt=1800)", fontsize=13)

    def panel(a, var2d_or_3d, title, ylabel, weighted=True, is3d=False):
        c = load_year(cpath(var2d_or_3d[0]), var2d_or_3d[0])
        f = load_year(fpath(var2d_or_3d[1]), var2d_or_3d[1])
        cy, fy = [], []
        for m in range(12):
            if is3d:
                cy.append(masked_mean(c[m])); fy.append(masked_mean(f[m]))
            else:
                cy.append(area_wmean(c[m], area)); fy.append(area_wmean(f[m], area))
        a.plot(MONTHS, fy, "-o", color=F, label="Fortran", ms=4)
        a.plot(MONTHS, cy, "-s", color=C, label="C port", ms=4)
        a.set_title(title); a.set_ylabel(ylabel); a.set_xticks(MONTHS)
        a.set_xticklabels(MNAMES); a.grid(alpha=0.3); a.legend(fontsize=8)
        return np.array(cy), np.array(fy)

    panel(ax[0,0], ("temp","temp"),  "3D temperature (vol mean)", "°C", is3d=True)
    panel(ax[0,1], ("salt","salt"),  "3D salinity (vol mean)",    "PSU", is3d=True)
    panel(ax[0,2], ("ssh","ssh"),    "SSH (area mean)",           "m")
    panel(ax[1,0], ("sst","sst"),    "SST (area mean)",           "°C")
    panel(ax[1,1], ("sss","sss"),    "SSS (area mean)",           "PSU")
    # 6th: divergence of the temp mean over the year
    ct = np.array([masked_mean(load_year(cpath("temp"),"temp")[m]) for m in range(12)])
    ft = np.array([masked_mean(load_year(fpath("temp"),"temp")[m]) for m in range(12)])
    ax[1,2].plot(MONTHS, (ct-ft)*1000, "-^", color="k")
    ax[1,2].set_title("3D temp mean drift (C − Fortran)"); ax[1,2].set_ylabel("m°C")
    ax[1,2].set_xticks(MONTHS); ax[1,2].set_xticklabels(MNAMES); ax[1,2].grid(alpha=0.3)
    ax[1,2].axhline(0, color="gray", lw=0.8)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(args.out_dir / "ocean_seasonal.png", dpi=110); plt.close(fig)

    # ---- Figure 2: sea-ice seasonal cycle ---------------------------
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"Sea-ice seasonal cycle {yr}: C-port vs Fortran", fontsize=13)

    a_c = load_year(cpath("a_ice"), "a_ice"); a_f = load_year(fpath("a_ice"), "a_ice")
    mi_c = load_year(cpath("m_ice"), "m_ice"); mi_f = load_year(fpath("m_ice"), "m_ice")
    ms_c = load_year(cpath("m_snow"), "m_snow"); ms_f = load_year(fpath("m_snow"), "m_snow")

    def integ(field, m):  # area-integrated (sum field*area) over valid nodes
        x = np.asarray(field[m]); msk = valid(x)
        return np.sum(x[msk] * area[msk])

    ice_area_c = [integ(a_c, m)/1e12 for m in range(12)]     # million km²
    ice_area_f = [integ(a_f, m)/1e12 for m in range(12)]
    ice_vol_c  = [integ(mi_c, m)/1e12 for m in range(12)]    # thousand km³ (m³/1e12)
    ice_vol_f  = [integ(mi_f, m)/1e12 for m in range(12)]
    snow_vol_c = [integ(ms_c, m)/1e12 for m in range(12)]
    snow_vol_f = [integ(ms_f, m)/1e12 for m in range(12)]
    # mean thickness where a_ice>0.15
    def mthick(mi, ai, m):
        x = np.asarray(mi[m]); a = np.asarray(ai[m]); msk = valid(x) & (a > 0.15)
        return x[msk].mean() if msk.sum() else np.nan
    th_c = [mthick(mi_c, a_c, m) for m in range(12)]
    th_f = [mthick(mi_f, a_f, m) for m in range(12)]

    for a, (yc, yf, title, yl) in zip(ax.flat, [
        (ice_area_c, ice_area_f, "Ice area  Σ(a_ice·area)", "10⁶ km²"),
        (ice_vol_c,  ice_vol_f,  "Ice volume  Σ(m_ice·area)", "10³ km³"),
        (snow_vol_c, snow_vol_f, "Snow volume  Σ(m_snow·area)", "10³ km³"),
        (th_c, th_f, "Mean ice thickness (a_ice>0.15)", "m"),
    ]):
        a.plot(MONTHS, yf, "-o", color=F, label="Fortran", ms=4)
        a.plot(MONTHS, yc, "-s", color=C, label="C port", ms=4)
        a.set_title(title); a.set_ylabel(yl); a.set_xticks(MONTHS)
        a.set_xticklabels(MNAMES); a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(args.out_dir / "seaice_seasonal.png", dpi=110); plt.close(fig)

    # ---- Figure 3: scatter agreement (surface T,S, Jul & Dec) -------
    fig, ax = plt.subplots(2, 2, figsize=(11, 10))
    fig.suptitle(f"C vs Fortran spatial agreement {yr} (each dot = 1 node)", fontsize=13)
    for col, (cv, fv, name, unit) in enumerate([
        ("sst","sst","SST","°C"), ("sss","sss","SSS","PSU")]):
        c_all = load_year(cpath(cv), cv); f_all = load_year(fpath(fv), fv)
        for row, mon in enumerate([6, 11]):   # Jul (idx6), Dec (idx11)
            a = ax[row, col]
            x = np.asarray(f_all[mon]); y = np.asarray(c_all[mon])
            m = valid(x) & valid(y)
            a.scatter(x[m], y[m], s=1, alpha=0.15, color="#444")
            lo = min(x[m].min(), y[m].min()); hi = max(x[m].max(), y[m].max())
            a.plot([lo,hi],[lo,hi], "r-", lw=1)
            corr = np.corrcoef(x[m], y[m])[0,1]
            a.set_title(f"{name} month {mon+1}  (r={corr:.4f})")
            a.set_xlabel(f"Fortran [{unit}]"); a.set_ylabel(f"C port [{unit}]")
            a.grid(alpha=0.3)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(args.out_dir / "scatter_agreement.png", dpi=110); plt.close(fig)

    # ---- Figure 4: spatial SST difference maps (Jan, Jul) -----------
    if lon is not None and lat is not None:
        fig, ax = plt.subplots(1, 2, figsize=(16, 5.5))
        fig.suptitle(f"SST difference (C − Fortran) {yr}", fontsize=13)
        c_all = load_year(cpath("sst"), "sst"); f_all = load_year(fpath("sst"), "sst")
        for a, mon, mlabel in [(ax[0], 0, "January"), (ax[1], 6, "July")]:
            d = np.asarray(c_all[mon]) - np.asarray(f_all[mon])
            m = valid(c_all[mon]) & valid(f_all[mon])
            sc = a.scatter(lon[m], lat[m], c=d[m], s=2, cmap="RdBu_r", vmin=-2, vmax=2)
            a.set_title(f"{mlabel}  (mean Δ={d[m].mean():+.3f}, rms={np.sqrt((d[m]**2).mean()):.3f} °C)")
            a.set_xlabel("lon"); a.set_ylabel("lat"); a.set_xlim(-180,180); a.set_ylim(-90,90)
            plt.colorbar(sc, ax=a, label="ΔSST [°C]")
        fig.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(args.out_dir / "maps_sst_diff.png", dpi=110); plt.close(fig)

    print(f"Wrote plots to {args.out_dir}")
    for p in sorted(args.out_dir.glob("*.png")):
        print(f"  {p}")


if __name__ == "__main__":
    main()
