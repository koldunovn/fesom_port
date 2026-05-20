#!/usr/bin/env python3
"""Spatial comparison maps (C-port vs Fortran) using nereus for proper
unstructured-mesh regridding + projected maps.

Run with the nereus conda env python:
  /work/ab0995/a270088/mambaforge/envs/nereus/bin/python plot_maps_nereus.py
"""
import pathlib
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import nereus as nr

C_DIR = pathlib.Path("/work/ab0995/a270088/port/core2_864_2yr")
F_DIR = pathlib.Path("/scratch/a/a270088/fesom_compare_1y")
MESH_DIAG = F_DIR / "fesom.mesh.diag.nc"
OUT = pathlib.Path("/home/a/a270088/port2/fesom2_port/docs/compare_plots")
OUT.mkdir(parents=True, exist_ok=True)
YR = 1958
FILL = 1e30

mesh = nr.fesom.load_mesh(MESH_DIAG)
lon = mesh["lon"].values
lat = mesh["lat"].values
print(f"mesh: {len(lon)} nodes, lon[{lon.min():.1f},{lon.max():.1f}] lat[{lat.min():.1f},{lat.max():.1f}]")


def cvar(v, m):  # C-port monthly mean, surface (or 2D)
    a = np.asarray(nc.Dataset(C_DIR / f"{v}.fesom.{YR}.monthly.nc")[v][m, ...])
    return a[0] if a.ndim == 2 else a   # 3D -> surface row


def fvar(v, m):
    a = np.asarray(nc.Dataset(F_DIR / f"{v}.fesom.{YR}.nc")[v][m, ...])
    return a[0] if a.ndim == 2 else a


def clean(a, mask_zero=True):
    # ocean 3D fields: 0 = below bathymetry (mask). ice fields: 0 = open
    # water (keep). Always mask netCDF fill values.
    a = np.array(a, dtype=float)
    a[np.abs(a) > FILL] = np.nan
    if mask_zero:
        a[a == 0] = np.nan
    return a


def triptych(var, month, mlabel, cmap, vmin, vmax, dvlim, unit, proj, fname,
             mask_zero=True):
    """Three-panel: C / Fortran / (C-Fortran) for one field+month."""
    c = clean(cvar(var, month), mask_zero); f = clean(fvar(var, month), mask_zero)
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.5),
                           subplot_kw={} )
    plt.close(fig)  # we'll let nr.plot make its own axes via ax= reuse
    fig = plt.figure(figsize=(20, 6))
    import cartopy.crs as ccrs
    projd = {"rob": ccrs.Robinson(), "npstere": ccrs.NorthPolarStereo()}[proj]
    panels = [("C port", c, cmap, vmin, vmax),
              ("Fortran", f, cmap, vmin, vmax),
              ("C − Fortran", c - f, "RdBu_r", -dvlim, dvlim)]
    for i, (ttl, dat, cm, lo, hi) in enumerate(panels):
        a = fig.add_subplot(1, 3, i + 1, projection=projd)
        nr.plot(dat, lon, lat, projection=proj, cmap=cm, vmin=lo, vmax=hi,
                ax=a, colorbar=True, colorbar_label=unit,
                title=f"{var} {mlabel} — {ttl}")
    fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fname}")


# SST — July (largest ocean differences), global Robinson. Tight diff scale.
triptych("sst", 6, "Jul", "RdYlBu_r", -2, 30, 1.0, "°C", "rob", "map_sst_jul.png")
# SSS — July, global
triptych("sss", 6, "Jul", "viridis", 28, 38, 1.0, "PSU", "rob", "map_sss_jul.png")
# a_ice — March (NH winter max), Arctic stereographic. Keep zeros (open water).
triptych("a_ice", 2, "Mar", "Blues", 0, 1, 0.3, "frac", "npstere", "map_aice_mar.png",
         mask_zero=False)
# m_ice — March, Arctic stereographic. Keep zeros.
triptych("m_ice", 2, "Mar", "Blues", 0, 4, 1.0, "m", "npstere", "map_mice_mar.png",
         mask_zero=False)

print("done")
