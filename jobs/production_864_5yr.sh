#!/bin/bash
#SBATCH --job-name=fesom_prod_kpp
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=03:00:00
#SBATCH -A ab0995
#SBATCH -o /work/ab0995/a270088/port/prod_kpp_5yr/slurm.%j.out
#SBATCH -e /work/ab0995/a270088/port/prod_kpp_5yr/slurm.%j.err
#
# Production CORE2 run — the validated reference configuration.
# KPP vertical mixing (the default), dt=1800, 864 ranks, 5 calendar years
# (1958-1962 JRA55-do) = 87600 steps. Monthly streams only, debug snapshots
# disabled (snap_every = -1). Needs the dist_864 partition under the mesh.

source /sw/etc/profile.levante
source /home/a/a270088/port2/fesom2_port/env.sh
ulimit -s 204800

REPO=/home/a/a270088/port2/fesom2_port
BIN=$REPO/build/fesom_port
MESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
PHC=/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc

# One output dir per job — never shared between concurrent jobs (jobs/README.md).
OUT=/work/ab0995/a270088/port/prod_kpp_5yr
mkdir -p "$OUT"
( cd "$REPO" && git rev-parse HEAD ) > "$OUT/SHA.txt" 2>/dev/null   # record the build

export FESOM_MIX_SCHEME=KPP        # default; set =PP to opt out (cf. production_pp.sh)
export FESOM_PRINT_EVERY=1000      # per-step stats every 1000 steps (~20.8 model-days)

# mesh  out   dt    nsteps snap_every  phc    jra_year
srun -l --ntasks=864 "$BIN" "$MESH" "$OUT" 1800 87600 -1 "$PHC" 1958 \
    > "$OUT/run.log" 2> "$OUT/run.err"
echo "exit: $?"
echo "--- tail run.log ---";   tail -3 "$OUT/run.log"
echo "--- blowup / uv>5 ? ---"; grep -iE "blowup|uv>5|nan" "$OUT/run.err" | tail -3
echo "--- monthly files ---";   ls "$OUT"/*.monthly.nc 2>/dev/null | wc -l
