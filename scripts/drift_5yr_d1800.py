#!/usr/bin/env python3
"""5-year model-drift comparison: C port vs Fortran, IDENTICAL config (dt=1800,
PP, linfs, same mesh/forcing). Only the code differs → any growing divergence is
true model drift. nereus diagnostics.
 C = c_5yr_dt1200 (monthly), F = fortran_5yr_d1200.

Outputs to docs/drift_5yr/:
 - drift_timeseries.png : 60-month global SST/SSS (nr.surface_mean), NH/SH ice
                          area & volume (nr.ice_area/volume_nh/sh); C, F, and C−F.
 - drift_bands.png      : annual volume-mean T/S in depth bands vs year (drift of
                          the interior, nr.volume_mean).
 - drift_maps.png       : SST/SSS C−F at year 1 vs year 5 (spatial drift growth).

Run: PYTHONPATH=/home/a/a270088/PYTHON \
  /work/ab0995/a270088/mambaforge/envs/nereus/bin/python drift_5yr.py
"""
import pathlib, warnings; warnings.filterwarnings("ignore")
import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import nereus as nr

C=pathlib.Path("/work/ab0995/a270088/port/stab_5yr_dt1800")
F=pathlib.Path("/scratch/a/a270088/fortran_pp_5yr_d1800")
DIAG=pathlib.Path("/scratch/a/a270088/fortran_pp_2yr/fesom.mesh.diag.nc")  # same mesh
OUT=pathlib.Path("/home/a/a270088/port2/fesom2_port/docs/drift_5yr_d1800"); OUT.mkdir(parents=True,exist_ok=True)
YEARS=[1958,1959,1960,1961,1962]; FILL=1e20
mesh=nr.fesom.load_mesh(DIAG)
area=np.asarray(mesh["area"].values); lat=np.asarray(mesh["lat"].values); lon=np.asarray(mesh["lon"].values)
depth=np.asarray(mesh["depth"].values); thick=np.asarray(mesh["layer_thickness"].values); NL=len(depth)
def clean(a): return np.where(np.abs(a)>FILL,np.nan,a)
def path(d,var,yr): return d/(f"{var}.fesom.{yr}.monthly.nc" if d==C else f"{var}.fesom.{yr}.nc")
def have(d,var,yr): return path(d,var,yr).exists()
def m2d(d,var,yr,mo): return clean(np.asarray(nc.Dataset(path(d,var,yr))[var][mo]))

# ---------- monthly time series over available months ----------
def collect(d):
    sst=[]; sss=[]; aN=[]; aS=[]; vN=[]; vS=[]; tt=[]
    for yi,yr in enumerate(YEARS):
        if not have(d,"sst",yr): break
        ns=nc.Dataset(path(d,"sst",yr)).dimensions["time"].size
        for mo in range(ns):
            sst.append(nr.surface_mean(m2d(d,"sst",yr,mo),area,mask=np.isfinite(m2d(d,"sst",yr,mo))))
            sss.append(nr.surface_mean(m2d(d,"sss",yr,mo),area,mask=np.isfinite(m2d(d,"sss",yr,mo))))
            ai=m2d(d,"a_ice",yr,mo); mi=m2d(d,"m_ice",yr,mo)
            aN.append(nr.ice_area_nh(ai,area,lat)/1e12); aS.append(nr.ice_area_sh(ai,area,lat)/1e12)
            vN.append(nr.ice_volume_nh(mi,area,lat)/1e12); vS.append(nr.ice_volume_sh(mi,area,lat)/1e12)
            tt.append(yi+mo/12.0)
    return dict(t=np.array(tt),sst=np.array(sst),sss=np.array(sss),aN=np.array(aN),aS=np.array(aS),vN=np.array(vN),vS=np.array(vS))
cc=collect(C); ff=collect(F)
n=min(len(cc["t"]),len(ff["t"]))
print(f"months available: C={len(cc['t'])} F={len(ff['t'])} -> comparing {n}")
if n>0:
    fig,ax=plt.subplots(3,2,figsize=(15,12))
    def panel(a,key,ttl,ylab):
        a.plot(cc["t"][:n],cc[key][:n],'-',c='tab:red',label='C'); a.plot(ff["t"][:n],ff[key][:n],'-',c='k',label='Fortran')
        a2=a.twinx(); d=cc[key][:n]-ff[key][:n]; a2.plot(cc["t"][:n],d,':',c='tab:blue',lw=1)
        a2.set_ylabel('C−F',color='tab:blue'); a2.tick_params(axis='y',labelcolor='tab:blue'); a2.axhline(0,c='tab:blue',lw=.4,ls=':')
        a.set_title(f"{ttl}   (C−F end: {d[-1]:+.3g})"); a.set_ylabel(ylab); a.set_xlabel('year'); a.grid(alpha=.3); a.legend(fontsize=8,loc='upper left')
    panel(ax[0,0],"sst","global SST","°C"); panel(ax[0,1],"sss","global SSS","PSU")
    panel(ax[1,0],"aN","NH ice area","10⁶ km²"); panel(ax[1,1],"aS","SH ice area","10⁶ km²")
    panel(ax[2,0],"vN","NH ice volume","10³ km³"); panel(ax[2,1],"vS","SH ice volume","10³ km³")
    fig.suptitle("5-yr drift: C vs Fortran (identical dt=1800/PP) — solid=value, dotted=C−F",fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig(OUT/"drift_timeseries.png",dpi=120,bbox_inches="tight"); print("wrote drift_timeseries.png")

# ---------- interior drift: annual volume-mean T/S bands ----------
def ann3d(d,var,yr):
    ds=nc.Dataset(path(d,var,yr)); ns=ds.dimensions["time"].size
    acc=np.zeros((NL,len(lon))); cnt=np.zeros((NL,len(lon)))
    for m in range(ns):
        a=clean(np.asarray(ds[var][m])); a=np.where(a==0,np.nan,a); ok=np.isfinite(a); acc[ok]+=a[ok]; cnt[ok]+=1
    ds.close(); return np.where(cnt>0,acc/np.where(cnt>0,cnt,1),np.nan)
bands=[(0,100),(100,700),(700,2000),(2000,6500)]
print("\n=== annual volume-mean T (°C) per band, C vs Fortran ===")
print("  band        " + "  ".join(f"{y}" for y in YEARS))
yrs_ok=[y for y in YEARS if have(C,"temp",y) and have(F,"temp",y)]
driftT={b:([],[]) for b in bands}; driftS={b:([],[]) for b in bands}
for yr in yrs_ok:
    tC=ann3d(C,"temp",yr); tF=ann3d(F,"temp",yr); sC=ann3d(C,"salt",yr); sF=ann3d(F,"salt",yr)
    for b in bands:
        driftT[b][0].append(nr.volume_mean(tC,area,thick,depth,depth_min=b[0],depth_max=b[1]))
        driftT[b][1].append(nr.volume_mean(tF,area,thick,depth,depth_min=b[0],depth_max=b[1]))
        driftS[b][0].append(nr.volume_mean(sC,area,thick,depth,depth_min=b[0],depth_max=b[1]))
        driftS[b][1].append(nr.volume_mean(sF,area,thick,depth,depth_min=b[0],depth_max=b[1]))
for b in bands:
    dc=np.array(driftT[b][0])-np.array(driftT[b][1])
    print(f"  {b[0]:4d}-{b[1]:<5d} ΔT C−F: " + "  ".join(f"{x:+.3f}" for x in dc))
if yrs_ok:
    fig,ax=plt.subplots(1,2,figsize=(14,5))
    for b in bands:
        ax[0].plot(yrs_ok,np.array(driftT[b][0])-np.array(driftT[b][1]),'-o',label=f"{b[0]}-{b[1]}m")
        ax[1].plot(yrs_ok,np.array(driftS[b][0])-np.array(driftS[b][1]),'-o',label=f"{b[0]}-{b[1]}m")
    ax[0].set_title("ΔT (C−F) per depth band vs year"); ax[0].set_ylabel("ΔT [°C]")
    ax[1].set_title("ΔS (C−F) per depth band vs year"); ax[1].set_ylabel("ΔS [PSU]")
    for a in ax: a.set_xlabel("year"); a.axhline(0,c='gray',lw=.6); a.grid(alpha=.3); a.legend(fontsize=8)
    fig.suptitle("Interior drift: volume-mean T/S difference vs year (nr.volume_mean)")
    fig.tight_layout(); fig.savefig(OUT/"drift_bands.png",dpi=120,bbox_inches="tight"); print("wrote drift_bands.png")

# ---------- spatial drift: SST/SSS C−F at year1 vs last available year ----------
def annsurf(d,var,yr):
    a=clean(np.asarray(nc.Dataset(path(d,var,yr))[var][:])); return np.nanmean(a,axis=0)
if yrs_ok:
    y1=yrs_ok[0]; yN=yrs_ok[-1]
    fig=plt.figure(figsize=(15,9))
    for i,(var,u,dv) in enumerate([("sst","°C",2.0),("sss","PSU",0.5)]):
        for j,yr in enumerate([y1,yN]):
            d=annsurf(C,var,yr)-annsurf(F,var,yr)
            ax=fig.add_subplot(2,2,i*2+j+1,projection=ccrs.Robinson())
            nr.plot(d,lon,lat,projection="rob",cmap="RdBu_r",vmin=-dv,vmax=dv,ax=ax,colorbar=True,colorbar_label=f"Δ{var} [{u}]",title=f"Δ{var.upper()} (C−F) {yr}")
    fig.suptitle(f"Spatial drift: surface C−F, year {y1} vs {yN}",fontsize=13)
    fig.savefig(OUT/"drift_maps.png",dpi=110,bbox_inches="tight"); print("wrote drift_maps.png")
print("done →",OUT)
