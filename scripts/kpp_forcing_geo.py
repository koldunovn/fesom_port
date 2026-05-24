import glob, os, re
import numpy as np, netCDF4 as nc
FN = re.compile(r"kpp_dump_s(\d+)_(.+)_rank(\d+)\.txt$")
def load(d, tag):
    out={}
    for f in glob.glob(os.path.join(d,f"kpp_dump_s*_{tag}_rank*.txt")):
        if not FN.search(os.path.basename(f)): continue
        for line in open(f):
            if line.startswith("#") or not line.strip(): continue
            p=line.split(); out[int(p[0])]=np.array([float(x) for x in p[1:]])
    return out
R="/work/ab0995/a270088/port/kpp_fdump/dump"       # Fortran
T="/work/ab0995/a270088/port/kpp_k6_validate/kpp_dump"  # C (live, not replay)
prR,prT=load(R,"prestep"),load(T,"prestep")
g=sorted(set(prR)&set(prT))
ur=np.array([prR[x][0] for x in g]); ut=np.array([prT[x][0] for x in g])
br=np.array([prR[x][1] for x in g]); bt=np.array([prT[x][1] for x in g])
du=np.abs(ur-ut); db=np.abs(br-bt)
# mesh lat/lon (global, 0-based by gid-1)
DIAG="/scratch/a/a270088/fortran_2yr_dt1800/fesom.mesh.diag.nc"
import os.path
if not os.path.exists(DIAG): DIAG="/scratch/a/a270088/fortran_pp_2yr/fesom.mesh.diag.nc"
ds=nc.Dataset(DIAG)
latv=np.asarray(ds.variables["lat"][:]); lonv=np.asarray(ds.variables["lon"][:])
if latv.max()<3.2: latv=np.degrees(latv); lonv=np.degrees(lonv)   # radians?
gi=np.array(g)-1
lat=latv[gi]; lon=lonv[gi]
print(f"nodes={len(g)}  |du|>1e-9: {int((du>1e-9).sum())} ({100*(du>1e-9).sum()/len(g):.0f}%)   |db|>1e-9: {int((db>1e-9).sum())} ({100*(db>1e-9).sum()/len(g):.0f}%)")
print("\n=== |d ustar| by latitude band (n, frac with |du|>1e-9, max|du|) ===")
bands=[("Antarctic <-60",lat<-60),("S mid -60..-40",(lat>=-60)&(lat<-40)),("S sub -40..-23",(lat>=-40)&(lat<-23)),
       ("Tropics |lat|<23",np.abs(lat)<23),("N sub 23..40",(lat>=23)&(lat<40)),("N mid 40..60",(lat>=40)&(lat<60)),
       ("Arctic >60",lat>60)]
for nm,m in bands:
    if m.sum()>0:
        frac=100*(du[m]>1e-9).sum()/m.sum()
        print(f"  {nm:18s} n={int(m.sum()):6d}  frac(du>1e-9)={frac:5.1f}%  max|du|={du[m].max():.2e}  max|db|={db[m].max():.2e}")
# the big-ustar-diff nodes: where are they?
big=du>1e-4
print(f"\n=== nodes with |d ustar|>1e-4 (n={int(big.sum())}): latitude distribution ===")
if big.sum()>0:
    print(f"  lat range [{lat[big].min():.1f}, {lat[big].max():.1f}]  median |lat|={np.median(np.abs(lat[big])):.1f}")
    print(f"  frac |lat|>50: {100*(np.abs(lat[big])>50).sum()/big.sum():.0f}%   frac |lat|>60: {100*(np.abs(lat[big])>60).sum()/big.sum():.0f}%")
