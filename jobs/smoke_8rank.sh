#!/bin/bash
#SBATCH --job-name=fesom_smoke
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH -A ab0995
#SBATCH -o /work/ab0995/a270088/port/smoke_8rank/slurm.%j.out
#SBATCH -e /work/ab0995/a270088/port/smoke_8rank/slurm.%j.err
#
# Quick sanity run: does the binary build, start, integrate, and write output?
# 8 ranks on a single node, dt=1800, ~1 model month (1488 steps = 31 days from
# 1958-01-01), debug snapshots every 5 days. Finishes in a few minutes.
# Needs the dist_8 partition under the mesh. See jobs/README.md.

source /sw/etc/profile.levante
source /home/a/a270088/port2/fesom2_port/env.sh
ulimit -s 204800

REPO=/home/a/a270088/port2/fesom2_port
BIN=$REPO/build/fesom_port
MESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
PHC=/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc

# One output dir per job — never shared between concurrent jobs (jobs/README.md).
OUT=/work/ab0995/a270088/port/smoke_8rank
mkdir -p "$OUT"
( cd "$REPO" && git rev-parse HEAD ) > "$OUT/SHA.txt" 2>/dev/null   # record the build

# mesh  out   dt    nsteps snap_every  phc    jra_year
srun -l "$BIN" "$MESH" "$OUT" 1800 1488 240 "$PHC" 1958 \
    > "$OUT/run.log" 2> "$OUT/run.err"
echo "exit: $?"
echo "--- tail run.log ---"; tail -20 "$OUT/run.log"
echo "--- snapshots ---";    ls "$OUT"/snap_*.nc 2>/dev/null | wc -l
echo "--- monthly files ---"; ls "$OUT"/*.monthly.nc 2>/dev/null
