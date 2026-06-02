#!/bin/bash
#SBATCH --job-name=fesom_valid_2yr
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=01:30:00
#SBATCH -A ab0995
#SBATCH -o /work/ab0995/a270088/port/valid_2yr/slurm.%j.out
#SBATCH -e /work/ab0995/a270088/port/valid_2yr/slurm.%j.err
#
# Validation run for comparison against the Fortran reference: dt=1800, 864
# ranks, 2 calendar years (1958-1959) = 35040 steps, default KPP mixing.
# Writes BOTH daily and monthly streams on every default variable via the
# io.config override, so the output can be diffed field-by-field against a
# matched Fortran run. Records the build SHA + working-tree status for
# reproducibility. Needs the dist_864 partition.

source /sw/etc/profile.levante
source /home/a/a270088/port2/fesom2_port/env.sh
ulimit -s 204800

REPO=/home/a/a270088/port2/fesom2_port
BIN=$REPO/build/fesom_port
MESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
PHC=/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc

# One output dir per job — never shared between concurrent jobs (jobs/README.md).
OUT=/work/ab0995/a270088/port/valid_2yr
mkdir -p "$OUT"
( cd "$REPO" && git rev-parse HEAD ) > "$OUT/SHA.txt"
( cd "$REPO" && git status --short ) > "$OUT/status.txt"

# Daily + monthly cadence on every default variable.
export FESOM_IO_CONFIG="$REPO/io.config.daily_monthly"
export FESOM_PRINT_EVERY=1000

# mesh  out   dt    nsteps snap_every  phc    jra_year
srun -l --ntasks=864 "$BIN" "$MESH" "$OUT" 1800 35040 -1 "$PHC" 1958 \
    > "$OUT/run.log" 2> "$OUT/run.err"
echo "exit: $?"
echo "--- tail run.log ---";    tail -3 "$OUT/run.log"
echo "--- blowup / uv>5 ? ---";  grep -iE "blowup|uv>5|nan" "$OUT/run.err" | tail -3
echo "--- daily files ---";      ls "$OUT"/*.daily.nc   2>/dev/null | wc -l
echo "--- monthly files ---";    ls "$OUT"/*.monthly.nc 2>/dev/null | wc -l
