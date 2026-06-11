#!/usr/bin/env python3
"""M5 end-to-end mEVP climate validation: does C+FESOM_WHICH_EVP=1 reproduce
the Fortran whichEVP=1 reference?

  Cm = c_mevp_2yr          (C, FESOM_WHICH_EVP=1, KPP+linfs, dt=1800, 1958-59)
  Fm = fortran_mevp_2yr    (Fortran whichEVP=1 — work_linfs_2yr_b clone,
                            frozen binary, the reference)
  Ce = c_evp_2yr           (C std-EVP default, SAME binary as Cm)
  Fe = fortran_linfs_2yr_b (Fortran whichEVP=0, SAME frozen binary as Fm)

Four gates (plan M5):
  1. Ocean: dSST/dSSS (Cm−Fm) area-weighted bias+RMS, global+regions, both
     years (drift check). Bar: TKE-class (~0.005/0.003, biases ~0).
  2. Ice: NH+SH extent (a_ice>0.15), volume, drift speed monthly cycles,
     C vs F (nan_to_num→0 BEFORE any time mean —
     feedback_nan_mask_temporal_average).
  3. VECTOR gate: per-node monthly-mean (uice,vice) magnitude ratio AND
     inter-vector angle Cm vs Fm (rotation-class bugs are invisible to
     speed/extent metrics — feedback_magnitude_invariant_masks_rotation).
  4. diff-of-diffs: (Cm−Ce) vs (Fm−Fe) pattern correlation for m_ice, drift
     speed, SST — must MATCH and be visibly NONZERO (option-live check).
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import nereus as nr

Cm = "/work/ab0995/a270088/port/mevp/c_mevp_2yr"
Fm = "/work/ab0995/a270088/port/mevp/fortran_mevp_2yr"
Ce = "/work/ab0995/a270088/port/mevp/c_evp_2yr"
Fe = "/work/ab0995/a270088/port/zstar/fortran_linfs_2yr_b"
DIAG = f"{Fm}/fesom.mesh.diag.nc"
FILL = 1e29

mesh = nr.fesom.load_mesh(DIAG)
lon = np.asarray(mesh["lon"].values); lat = np.asarray(mesh["lat"].values)
if lon.max() > 180: lon = np.where(lon > 180, lon - 360, lon)
ea = np.asarray(nc.Dataset(DIAG).variables["elem_area"][:])
fn = np.asarray(nc.Dataset(DIAG).variables["face_nodes"][:]); fn = fn.T if fn.shape[0] == 3 else fn
fn = fn - 1 if fn.min() == 1 else fn
na = np.zeros(len(lat))
for k in range(3): np.add.at(na, fn[:, k], ea / 3.0)

def monthly(path, var, suffix, year):
    """(12, nod) with NaN->0 BEFORE any averaging (the ice-edge open-water
    months trap). Ocean fields keep NaN (land mask) — caller chooses."""
    a = np.asarray(nc.Dataset(f"{path}/{var}.fesom.{year}{suffix}.nc").variables[var][:])
    return np.where(np.abs(a) < FILL, a, np.nan)

def ice_monthly(path, var, suffix, year):
    return np.nan_to_num(monthly(path, var, suffix, year), nan=0.0)

def wmean(d, m): return np.sum(d[m] * na[m]) / np.sum(na[m])
def wrms(d, m):  return np.sqrt(np.sum(d[m] ** 2 * na[m]) / np.sum(na[m]))

REGIONS = lambda wet: [
    ("GLOBAL",           wet),
    ("Tropics |lat|<23", wet & (np.abs(lat) < 23)),
    ("Midlat 40-60",     wet & (np.abs(lat) >= 40) & (np.abs(lat) < 60)),
    ("N high >60",       wet & (lat > 60)),
    ("S high <-60",      wet & (lat < -60)),
]

print("#" * 84)
print("# GATE 1 — ocean climate: C-mEVP vs Fortran-mEVP (and vs F-EVP contrast)")
print("#" * 84)
for year in (1958, 1959):
    print(f"\nYEAR {year}  dSST/dSSS area-weighted   mEVP = Cm−Fm (headline) | "
          f"xEVP = Cm−Fe (rheology contrast)")
    for var, unit in [("sst", "degC"), ("sss", "PSU")]:
        fF = np.nanmean(monthly(Fm, var, "", year), 0)
        fC = np.nanmean(monthly(Cm, var, ".monthly", year), 0)
        fE = np.nanmean(monthly(Fe, var, "", year), 0)
        wet = np.isfinite(fF) & np.isfinite(fC) & np.isfinite(fE)
        dM = np.where(wet, fC - fF, np.nan)
        dX = np.where(wet, fC - fE, np.nan)
        print(f"  d{var.upper()} [{unit}]  {'region':18s} {'mEVP bias':>10s} {'mEVP rms':>9s} | "
              f"{'xEVP bias':>10s} {'xEVP rms':>9s}")
        for name, m in REGIONS(wet):
            print(f"               {name:18s} {wmean(dM, m):+10.4f} {wrms(dM, m):9.4f} | "
                  f"{wmean(dX, m):+10.4f} {wrms(dX, m):9.4f}")

print()
print("#" * 84)
print("# GATE 2 — ice metrics: NH+SH extent / volume / drift-speed seasonal cycles")
print("#" * 84)
year = 1959
aC = ice_monthly(Cm, "a_ice", ".monthly", year); aF = ice_monthly(Fm, "a_ice", "", year)
mC = ice_monthly(Cm, "m_ice", ".monthly", year); mF = ice_monthly(Fm, "m_ice", "", year)
uC = ice_monthly(Cm, "uice",  ".monthly", year); uF = ice_monthly(Fm, "uice", "", year)
vC = ice_monthly(Cm, "vice",  ".monthly", year); vF = ice_monthly(Fm, "vice", "", year)
spC = np.sqrt(uC**2 + vC**2); spF = np.sqrt(uF**2 + vF**2)
for hemi, hm in [("NH", lat > 0), ("SH", lat < 0)]:
    print(f"\n  {hemi} {year}:  month   extent[1e6 km2]  C / F      volume[1e3 km3]  C / F"
          f"      drift[cm/s, a>0.15]  C / F")
    for mo in range(12):
        extC = np.sum(na[hm & (aC[mo] > 0.15)]) / 1e12
        extF = np.sum(na[hm & (aF[mo] > 0.15)]) / 1e12
        volC = np.sum(na[hm] * mC[mo][hm]) / 1e12
        volF = np.sum(na[hm] * mF[mo][hm]) / 1e12
        mskC = hm & (aC[mo] > 0.15); mskF = hm & (aF[mo] > 0.15)
        drC = 100 * np.mean(spC[mo][mskC]) if mskC.any() else 0.0
        drF = 100 * np.mean(spF[mo][mskF]) if mskF.any() else 0.0
        print(f"        {mo+1:5d}   {extC:7.2f} / {extF:7.2f}        "
              f"{volC:6.2f} / {volF:6.2f}          {drC:6.2f} / {drF:6.2f}")

print()
print("#" * 84)
print("# GATE 3 — ice-velocity VECTOR check Cm vs Fm (magnitude ratio + angle)")
print("#" * 84)
# FRAME ALIGNMENT (found 2026-06-11): the C IO writes uice/vice in the NATIVE
# ROTATED frame (fesom_io_stream.c resolve_uice accumulates ice->uice raw);
# the Fortran io_meandata sets do_rotation=.TRUE. for (uice,vice) and writes
# GEOGRAPHIC components. Comparing raw C vs F shows ratio==1 (rotation is an
# isometry) but position-dependent angles up to +-40deg — the exact
# magnitude-invariant signature feedback_magnitude_invariant_masks_rotation
# warns about. Rotate the C vectors r2g (gen_modules_rotate_grid.F90 port)
# before comparing.
al_e, be_e, ga_e = np.radians([50.0, 15.0, -90.0])   # alphaEuler/betaEuler/gammaEuler
RM = np.array([
 [np.cos(ga_e)*np.cos(al_e)-np.sin(ga_e)*np.cos(be_e)*np.sin(al_e),
  np.cos(ga_e)*np.sin(al_e)+np.sin(ga_e)*np.cos(be_e)*np.cos(al_e),
  np.sin(ga_e)*np.sin(be_e)],
 [-np.sin(ga_e)*np.cos(al_e)-np.cos(ga_e)*np.cos(be_e)*np.sin(al_e),
  -np.sin(ga_e)*np.sin(al_e)+np.cos(ga_e)*np.cos(be_e)*np.cos(al_e),
  np.cos(ga_e)*np.sin(be_e)],
 [np.sin(be_e)*np.sin(al_e), -np.sin(be_e)*np.cos(al_e), np.cos(be_e)]])
_glon = np.radians(lon); _glat = np.radians(lat)
_xg = np.cos(_glat)*np.cos(_glon); _yg = np.cos(_glat)*np.sin(_glon); _zg = np.sin(_glat)
_xr = RM[0,0]*_xg + RM[0,1]*_yg + RM[0,2]*_zg
_yr = RM[1,0]*_xg + RM[1,1]*_yg + RM[1,2]*_zg
_zr = RM[2,0]*_xg + RM[2,1]*_yg + RM[2,2]*_zg
_rlat = np.arcsin(np.clip(_zr, -1, 1)); _rlon = np.arctan2(_yr, _xr)

def vec_r2g(tlon, tlat):
    """Port of vector_r2g (gen_modules_rotate_grid.F90:164-202)."""
    txg = -tlat*np.sin(_rlat)*np.cos(_rlon) - tlon*np.sin(_rlon)
    tyg = -tlat*np.sin(_rlat)*np.sin(_rlon) + tlon*np.cos(_rlon)
    tzg =  tlat*np.cos(_rlat)
    txr = RM[0,0]*txg + RM[1,0]*tyg + RM[2,0]*tzg   # M^T (vector_r2g:195-197)
    tyr = RM[0,1]*txg + RM[1,1]*tyg + RM[2,1]*tzg
    tzr = RM[0,2]*txg + RM[1,2]*tyg + RM[2,2]*tzg
    vlat = -np.sin(_glat)*np.cos(_glon)*txr - np.sin(_glat)*np.sin(_glon)*tyr + np.cos(_glat)*tzr
    vlon = -np.sin(_glon)*txr + np.cos(_glon)*tyr
    return vlon, vlat

for mo, label in [(2, "Mar"), (8, "Sep")]:
    uCg, vCg = vec_r2g(uC[mo], vC[mo])
    spCg = np.sqrt(uCg**2 + vCg**2)
    for hemi, hm in [("NH", lat > 0), ("SH", lat < 0)]:
        m = hm & (aF[mo] > 0.15) & (spF[mo] > 0.005) & (spCg > 0.005)
        if m.sum() < 10:
            print(f"  {label} {hemi}: <10 valid nodes, skipped"); continue
        ratio = spCg[m] / spF[mo][m]
        dot   = uCg[m]*uF[mo][m] + vCg[m]*vF[mo][m]
        cross = uF[mo][m]*vCg[m] - vF[mo][m]*uCg[m]
        ang   = np.degrees(np.arctan2(cross, dot))
        wts   = na[m]
        print(f"  {label} 1959 {hemi} (n={m.sum():6d}): |u| ratio median={np.median(ratio):.4f} "
              f"[IQR {np.percentile(ratio,25):.4f}..{np.percentile(ratio,75):.4f}]  "
              f"angle median={np.median(ang):+.3f}deg "
              f"[IQR {np.percentile(ang,25):+.3f}..{np.percentile(ang,75):+.3f}]  "
              f"wmean={np.sum(ang*wts)/np.sum(wts):+.3f}deg")

print()
print("#" * 84)
print("# GATE 4 — diff-of-diffs: (Cm−Ce) vs (Fm−Fe), annual means 1959")
print("#" * 84)
def annual_ice(path, var, suffix):  return np.mean(ice_monthly(path, var, suffix, 1959), 0)
def annual_oce(path, var, suffix):  return np.nanmean(monthly(path, var, suffix, 1959), 0)
mCe = annual_ice(Ce, "m_ice", ".monthly"); mFe = annual_ice(Fe, "m_ice", "")
uCe = annual_ice(Ce, "uice", ".monthly"); vCe = annual_ice(Ce, "vice", ".monthly")
uFe = annual_ice(Fe, "uice", ""); vFe = annual_ice(Fe, "vice", "")
spCe = np.sqrt(uCe**2 + vCe**2); spFe = np.sqrt(uFe**2 + vFe**2)
mCm_a = np.mean(mC, 0); mFm_a = np.mean(mF, 0)
spCm_a = np.mean(spC, 0); spFm_a = np.mean(spF, 0)
cases = [
    ("m_ice [m]",       mCm_a - mCe,  mFm_a - mFe,  (np.mean(aC,0) > 0.15) | (np.mean(aF,0) > 0.15)),
    ("drift speed [m/s]", spCm_a - spCe, spFm_a - spFe, (np.mean(aC,0) > 0.15) | (np.mean(aF,0) > 0.15)),
]
sstCm = annual_oce(Cm, "sst", ".monthly"); sstCe = annual_oce(Ce, "sst", ".monthly")
sstFm = annual_oce(Fm, "sst", "");          sstFe = annual_oce(Fe, "sst", "")
wetsst = np.isfinite(sstCm) & np.isfinite(sstCe) & np.isfinite(sstFm) & np.isfinite(sstFe)
cases.append(("SST [degC]", np.where(wetsst, sstCm - sstCe, 0.0),
                            np.where(wetsst, sstFm - sstFe, 0.0), wetsst))
for name, dCfld, dFfld, m in cases:
    dC = np.nan_to_num(dCfld, nan=0.0)[m]; dF = np.nan_to_num(dFfld, nan=0.0)[m]; w = na[m]
    wm = lambda x: np.sum(x * w) / np.sum(w)
    cC = dC - wm(dC); cF = dF - wm(dF)
    corr = wm(cC * cF) / np.sqrt(wm(cC**2) * wm(cF**2))
    print(f"  {name:18s} pattern corr={corr:+.3f}  "
          f"wstd C={np.sqrt(wm(cC**2)):.4f} F={np.sqrt(wm(cF**2)):.4f}  "
          f"max|C|={np.abs(dC).max():.3f} max|F|={np.abs(dF).max():.3f}  "
          f"(match + visibly nonzero = option live & faithful)")
print("\ndone.")
