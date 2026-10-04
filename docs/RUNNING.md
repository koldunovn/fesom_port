# Running fesom_port outside Levante

This page shows how to build the C port on any Linux cluster and run the CORE2 configuration
used in the paper for one model year, using the input data archived on Zenodo.

## 1. Build

Requirements: a C99 compiler, MPI, netCDF-C (with `nc-config` on the `PATH`) and CMake 3.18 or
newer. `configure.sh` and `env.sh` load the Levante modules; on another machine, load your own
compiler, MPI and netCDF and call CMake directly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpicc
cmake --build build -j
# -> build/fesom_port
```

Use `Release` (`-O3`). Without a build type CMake builds `Debug` (`-O0`), which is correct but slow.
The paper's binaries were built with GCC 11.2, Open MPI 4.1.2 and netCDF-C 4.8.1.

## 2. Input data

Everything the CORE2 runs of the paper read for the year 1958 is archived as a single Zenodo record,
[doi:10.5281/zenodo.21324318](https://doi.org/10.5281/zenodo.21324318) (version 1.2,
[doi:10.5281/zenodo.21770779](https://doi.org/10.5281/zenodo.21770779)). Its files are
byte-identical to the ones the paper's runs read.

```bash
mkdir -p ~/fesom-data && cd ~/fesom-data
for f in core2_mesh_ic.zip core2_partitions.zip core2_forcing_1958.zip; do
    wget "https://zenodo.org/records/21770779/files/$f?download=1" -O "$f"
    unzip -q "$f"
done
```

`core2_forcing_1958.zip` is about 10 GB. The archives contain:

| Path | Content |
|---|---|
| `core2_mesh_ic/mesh_core2_raw/` | the CORE2 mesh in FESOM2 text format (`nod2d.out`, `elem2d.out`, `aux3d.out`, `nlvls.out`, `elvls.out`, `edges.out`, `edge_tri.out`, `edgenum.out`) |
| `core2_partitions/dist_N/` | domain decompositions for N = 2 to 864 ranks |
| `core2_mesh_ic/phc3.0_winter.nc` | PHC3.0 initial temperature and salinity |
| `core2_mesh_ic/PHC2_salx.nc` | sea-surface salinity for restoring |
| `core2_mesh_ic/CORE2_runoff.nc` | river runoff |
| `core2_mesh_ic/Sweeney_2005.nc` | chlorophyll climatology for shortwave penetration |
| `core2_forcing_1958/` | the eight JRA55-do v1.4.0 fields for 1958 (`uas vas tas huss rsds rlds prra prsn`), named `<var>.1958.nc` |

The model looks for its decomposition in `<mesh_dir>/dist_<number of MPI ranks>`, so link the
partitions into the mesh directory:

```bash
cd ~/fesom-data/core2_mesh_ic/mesh_core2_raw
for d in ~/fesom-data/core2_partitions/dist_*; do ln -s "$d" .; done
```

For years after 1958, add the JRA55-do files of those years to the forcing directory, named
`<var>.<year>.nc`. The paper used JRA55-do v1.4.0 for 1958 to 2018
([doi:10.22033/ESGF/input4MIPs.10842](https://doi.org/10.22033/ESGF/input4MIPs.10842)) and v1.5.0
for 2019 and 2020 ([doi:10.22033/ESGF/input4MIPs.15017](https://doi.org/10.22033/ESGF/input4MIPs.15017)),
both distributed through the Earth System Grid Federation.

## 3. Run one year of CORE2

The paper's configuration is z* vertical coordinate, cvmix_TKE vertical mixing, mEVP sea-ice
dynamics and the deterministic initial-condition fill, at a time step of 1800 s:

```bash
D=~/fesom-data
export FESOM_ALE=zstar FESOM_MIX_SCHEME=TKE FESOM_WHICH_EVP=1 FESOM_IC_EXTRAP=det
export FESOM_FORCING_DIR=$D/core2_forcing_1958
export FESOM_SSS_PATH=$D/core2_mesh_ic/PHC2_salx.nc
export FESOM_RUNOFF_PATH=$D/core2_mesh_ic/CORE2_runoff.nc
export FESOM_CHL_FILE=$D/core2_mesh_ic/Sweeney_2005.nc
mkdir -p out
mpirun -np 128 build/fesom_port $D/core2_mesh_ic/mesh_core2_raw out 1800 17520 -1 \
       $D/core2_mesh_ic/phc3.0_winter.nc 1958
```

The positional arguments are `<mesh_dir> <output_dir> <dt_s> <nsteps> <snapshot_every> <phc_file>
<first_forcing_year>`: 17520 steps of 1800 s are the 365 days of 1958, and `-1` switches the debug
snapshots off. The output is a set of monthly-mean files, `out/<variable>.fesom.1958.monthly.nc`.

The paper's C hindcast ran on 512 ranks (`dist_512`). Results depend on the rank count at round-off
level, because the order of global sums and halo exchanges follows the decomposition.

## 4. Restarts and longer runs

| Variable | Meaning |
|---|---|
| `FESOM_RESTART_OUT=<dir>` | write restart files into `<dir>` |
| `FESOM_RESTART_AT=<s1,s2,...>` | write a restart after these step numbers |
| `FESOM_RESTART_EVERY=<n>` | write a restart every `n` steps (0 = off) |
| `FESOM_RESTART_IN=<file>` | start from this restart file instead of the initial condition |

The paper's hindcast wrote a restart on every 1 January and ran in 21-year segments; the job scripts
are part of the paper's archive.

## 5. Other options

| Variable | Values |
|---|---|
| `FESOM_ALE` | `linfs` (default), `zstar` |
| `FESOM_MIX_SCHEME` | unset (KPP, the default), `PP`, `TKE` |
| `FESOM_WHICH_EVP` | `0` (EVP, default), `1` (mEVP) |
| `FESOM_IC_EXTRAP` | `legacy` (default), `det` |
| `FESOM_IO_CONFIG` | file that changes the output variables and their averaging period |

Section "Environment knobs" of the README lists the diagnostic switches.
