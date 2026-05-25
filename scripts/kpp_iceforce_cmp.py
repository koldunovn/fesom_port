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
R="/work/ab0995/a270088/port/kpp_fdump/dump"
T="/work/ab0995/a270088/port/kpp_iceforce/dump"
fr,ft=load(R,"iceforce"),load(T,"iceforce")
g=sorted(set(fr)&set(ft)); n=len(g)
A=np.array([fr[x] for x in g]); B=np.array([ft[x] for x in g])   # (n,10)
DIAG="/scratch/a/a270088/fortran_2yr_dt1800/fesom.mesh.diag.nc"
ds=nc.Dataset(DIAG); latv=np.asarray(ds.variables["lat"][:])
if latv.max()<3.2: latv=np.degrees(latv)
lat=latv[np.array(g)-1]
polar = np.abs(lat)>50
names=["a_ice","m_ice","u_ice","v_ice","srfoce_u","srfoce_v","stress_ioc_x","stress_ioc_y","stress_ns_x","stress_ns_y"]
print(f"common nodes={n}  polar(|lat|>50)={int(polar.sum())}")
print(f"\n{'field':13s} {'glob max|d|':>12s} {'polar max|d|':>12s} {'polar max|val|':>13s} {'polar #(|d|>1e-12)':>18s}")
for c,nm in enumerate(names):
    d=np.abs(A[:,c]-B[:,c])
    pv=np.abs(B[polar,c])
    print(f"{nm:13s} {d.max():12.3e} {d[polar].max():12.3e} {pv.max():13.3e} {int((d[polar]>1e-12).sum()):18d}")
# decisive: where stress_node_surf differs, is it a_ice or u_ice?
dsns=np.maximum(np.abs(A[:,8]-B[:,8]),np.abs(A[:,9]-B[:,9]))
da=np.abs(A[:,0]-B[:,0]); duv=np.maximum(np.abs(A[:,2]-B[:,2]),np.abs(A[:,3]-B[:,3]))
big=dsns>1e-6
print(f"\nnodes with |d stress_node_surf|>1e-6: {int(big.sum())}")
if big.sum()>0:
    print(f"  among them: max|d a_ice|={da[big].max():.3e}   max|d (u_ice,v_ice)|={duv[big].max():.3e}")
    print(f"  a_ice bit-identical (|d|<1e-12) at {100*(da[big]<1e-12).sum()/big.sum():.0f}% of them;  u/v_ice differs (|d|>1e-9) at {100*(duv[big]>1e-9).sum()/big.sum():.0f}%")
print(f"\nsrfoce_u/v global max|val| (both codes): C={np.abs(B[:,4:6]).max():.3e}  F={np.abs(A[:,4:6]).max():.3e}  (≈0 confirms ocean-at-rest)")
