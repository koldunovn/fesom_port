#!/usr/bin/env python3
"""Illustrate the central-Arctic cell-velocity 2Δx checkerboard in the C port.

Per-TRIANGLE (element) coloring of the surface velocity from C-port snapshots,
zoomed into the central Arctic. The checkerboard shows as adjacent triangles
flipping sign/magnitude. Three times show the growth: smooth -> developing ->
erupting.

Run:
  PYTHONPATH=/home/a/a270088/PYTHON \
  /work/ab0995/a270088/mambaforge/envs/nereus/bin/python plot_checkerboard.py
"""
import pathlib
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import netCDF4 as nc
import cartopy.crs as ccrs

SNAP_DIR = pathlib.Path("/work/ab0995/a270088/port/dt1800_snap")
OUT = pathlib.Path("/home/a/a270088/port2/fesom2_port/docs/compare_plots_dt1800")
OUT.mkdir(parents=True, exist_ok=True)
SNAPS = [(1340, "step 1340 (day 28) — smooth"),
         (4020, "step 4020 (day 84) — developing"),
         (5200, "step 5200 (day 108) — erupting")]
DT = 1800.0

proj = ccrs.NorthPolarStereo()
geo  = ccrs.PlateCarree()


def load(step):
    ds = nc.Dataset(SNAP_DIR / f"snap_{step:06d}.nc")
    lon = np.asarray(ds["lon"][:]); lat = np.asarray(ds["lat"][:])
    en  = np.asarray(ds["elem_nodes"][:])           # (elem, 3)
    if en.min() == 1:                                # 1-based -> 0-based
        en = en - 1
    ue = np.asarray(ds["u"][0, 0, :])                # surface element u
    ve = np.asarray(ds["v"][0, 0, :])
    ds.close()
    return lon, lat, en, ue, ve


# Use the last snapshot to locate the eruption (max surface speed) and frame.
lon, lat, en, ue, ve = load(SNAPS[-1][0])
xy = proj.transform_points(geo, lon, lat)
x, y = xy[:, 0], xy[:, 1]
# element centroids
cx = x[en].mean(1); cy = y[en].mean(1)
spd = np.hypot(ue, ve)
# focus element = max speed among elements north of 80N
clat = lat[en].mean(1)
cand = np.where(clat > 80.0)[0]
foc = cand[np.argmax(spd[cand])]
x0, y0 = cx[foc], cy[foc]
HALF = 110e3   # 110 km half-window — tight enough to resolve per-triangle 2dx
print(f"focus element {foc}: centroid lat={clat[foc]:.2f} speed={spd[foc]:.3f} m/s")

fig, axes = plt.subplots(1, 3, figsize=(19, 7))
for ax, (step, label) in zip(axes, SNAPS):
    lon, lat, en, ue, ve = load(step)
    xy = proj.transform_points(geo, lon, lat)
    x, y = xy[:, 0], xy[:, 1]
    triang = mtri.Triangulation(x, y, triangles=en)
    cx = x[en].mean(1); cy = y[en].mean(1)
    # keep only triangles inside the window (for speed + clean axes)
    inwin = (np.abs(cx - x0) < HALF) & (np.abs(cy - y0) < HALF)
    mask = ~inwin
    triang.set_mask(mask)
    uu = np.asarray(ue)
    tpc = ax.tripcolor(triang, facecolors=uu, cmap="RdBu_r",
                       vmin=-0.5, vmax=0.5, shading="flat",
                       edgecolors="0.25", linewidth=0.25)
    ax.set_xlim(x0 - HALF, x0 + HALF); ax.set_ylim(y0 - HALF, y0 + HALF)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{label}\nmax|u| in window = {np.abs(uu[inwin]).max():.2f} m/s",
                 fontsize=11)
    plt.colorbar(tpc, ax=ax, shrink=0.7, label="surface u  [m/s]  (per cell)")

fig.suptitle("C port dt=1800 — central-Arctic surface CELL velocity (per-triangle): "
             "2Δx checkerboard growth", fontsize=13)
fig.tight_layout()
fp = OUT / "checkerboard_cellvel.png"
fig.savefig(fp, dpi=130, bbox_inches="tight")
print("wrote", fp)
