# SLURM job templates

Curated, ready-to-submit SLURM batch scripts for running the port on Levante.
Each is self-contained: it sources `env.sh`, points at the CORE2 mesh and PHC
initial condition, and calls `build/fesom_port` with the positional arguments
documented in the top-level [`README.md`](../README.md#running).

| Script | Ranks | dt | Length | Mixing | Output | Use |
|---|---|---|---|---|---|---|
| `smoke_8rank.sh` | 8 (1 node) | 1800 | ~1 month | KPP | snapshots + monthly | Quick "does it build, run, and write output?" check — a few minutes. |
| `production_864_5yr.sh` | 864 | 1800 | 5 yr (1958–1962) | KPP (default) | monthly streams | The validated reference production run. |
| `production_pp.sh` | 864 | 1800 | 5 yr | PP | monthly streams | Same as above with Pacanowski–Philander mixing (`FESOM_MIX_SCHEME=PP`). |
| `validation_2yr.sh` | 864 | 1800 | 2 yr (1958–1959) | KPP | daily + monthly | Field-by-field comparison against a matched Fortran run; records build SHA. |

Submit with `sbatch jobs/<script>` from the repo root.

## Before you submit

1. **Build first:** `bash -l configure.sh` (see the top-level README). The
   scripts run `build/fesom_port`; they do not build it.
2. **Partitions:** multi-rank runs need a pre-generated `dist_<NPES>` directory
   under the mesh for the exact `--ntasks`. The templates use 8 and 864; other
   counts seen in this tree: 1, 16, 32, 144, 256, 512.
3. **Create the output dir / log dir.** SLURM needs the `-o`/`-e` directory to
   exist at submit time, so create the `OUT` path before the first `sbatch`
   (the script also `mkdir -p`s it, which covers reruns):
   `mkdir -p /work/ab0995/a270088/port/<name>`.

## Adapting a template

Everything Levante-specific lives in a handful of variables near the top —
edit these and the `#SBATCH` lines to match your run:

| What | Where |
|---|---|
| Rank count | `#SBATCH --ntasks` / `--nodes` (+ `--ntasks=` on `srun`) and the `dist_<N>` partition |
| Output location | `OUT=...` **and** the `#SBATCH -o/-e` paths (keep them in sync) |
| Run length | the `nsteps` positional arg (1 yr ≈ 17520 steps at dt=1800) |
| Mixing scheme | `export FESOM_MIX_SCHEME=KPP` / `PP` |
| Forcing year | the trailing `jra_year` arg (model calendar starts here) |
| Paths (`MESH`, `PHC`, `BIN`) | the variable block — all hard-coded to Levante |

**Always use a unique `OUT` per job.** Concurrent jobs that share an output
directory clobber each other's logs and NetCDF and produce spurious "crashes"
that are really log corruption. Change `OUT`, the `#SBATCH -o/-e` paths, and
`--job-name` together for each new experiment.

For the full argument reference, environment knobs (diagnostics, physics-bisect
toggles), and end-of-run sanity checks, see the top-level
[`README.md`](../README.md).
