#!/usr/bin/env python3
"""C-port (dt=1800, monthly) vs Fortran (dt=1800 2yr ref) comparison maps,
using nereus for unstructured-mesh projected maps.

The C dt=1800 run blew up at step 5386 (~day 112) so only Jan/Feb/Mar 1958
monthly means exist. We map March (month index 2) — NH ice maximum and the
last healthy C month before the central-Arctic blow-up.

Run:
  PYTHONPATH=/home/a/a270088/PYTHON \
  /work/ab0995/a270088/mambaforge/envs/nereus/bin/python plot_maps_dt1800.py
"""
import pathlib
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
import nereus as nr

C_DIR = pathlib.Path("/work/ab0995/a270088/port/core2_864_2yr_dt1800_monthly")
F_DIR = pathlib.Path("/scratch/a/a270088/fortran_2yr_dt1800")
MESH_DIAG = F_DIR / "fesom.mesh.diag.nc"
OUT = pathlib.Path("/home/a/a270088/port2/fesom2_port/docs/compare_plots_dt1800")
OUT.mkdir(parents=True, exist_ok=True)
YR = 1958
MON = 2          # March (0-based) — last healthy C month
MLABEL = "Mar"
FILL = 1e30

mesh = nr.fesom.load_mesh(MESH_DIAG)
lon = mesh["lon"].values
lat = mesh["lat"].values
print(f"mesh: {len(lon)} nodes, lon[{lon.min():.1f},{lon.max():.1f}] lat[{lat.min():.1f},{lat.max():.1f}]")


def cvar(v, m):
    a = np.asarray(nc.Dataset(C_DIR / f"{v}.fesom.{YR}.monthly.nc")[v][m, ...])
    return a[0] if a.ndim == 2 else a


def fvar(v, m):
    a = np.asarray(nc.Dataset(F_DIR / f"{v}.fesom.{YR}.nc")[v][m, ...])
    return a[0] if a.ndim == 2 else a


def clean(a, mask_zero=True):
    a = np.array(a, dtype=float)
    a[np.abs(a) > FILL] = np.nan
    if mask_zero:
        a[a == 0] = np.nan
    return a


def triptych(var, cmap, vmin, vmax, dvlim, unit, proj, fname, mask_zero=True):
    c = clean(cvar(var, MON), mask_zero)
    f = clean(fvar(var, MON), mask_zero)
    fig = plt.figure(figsize=(20, 6))
    projd = {"rob": ccrs.Robinson(), "npstere": ccrs.NorthPolarStereo()}[proj]
    panels = [("C port", c, cmap, vmin, vmax),
              ("Fortran", f, cmap, vmin, vmax),
              ("C - Fortran", c - f, "RdBu_r", -dvlim, dvlim)]
    for i, (ttl, dat, cm, lo, hi) in enumerate(panels):
        a = fig.add_subplot(1, 3, i + 1, projection=projd)
        nr.plot(dat, lon, lat, projection=proj, cmap=cm, vmin=lo, vmax=hi,
                ax=a, colorbar=True, colorbar_label=unit,
                title=f"{var} {MLABEL} {YR} - {ttl}")
    fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fname}")


# Ocean — global Robinson
triptych("sst", "RdYlBu_r", -2, 30, 1.0, "degC", "rob", "map_sst_mar.png")
triptych("sss", "viridis", 28, 38, 1.0, "PSU", "rob", "map_sss_mar.png")
triptych("ssh", "RdBu_r", -2, 1.5, 0.1, "m", "rob", "map_ssh_mar.png")
# Sea ice — Arctic stereographic (keep zeros = open water)
triptych("a_ice", "Blues", 0, 1, 0.3, "frac", "npstere", "map_aice_mar.png", mask_zero=False)
triptych("m_ice", "Blues", 0, 4, 1.0, "m", "npstere", "map_mice_mar.png", mask_zero=False)
# Vertical mixing — surface Kv. C caps at 0.1; Fortran has convective values.
triptych("Kv", "magma", 0, 0.5, 0.3, "m2/s", "npstere", "map_kv_mar.png")

print("done")
