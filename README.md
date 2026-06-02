# FESOM2 C port

A literal C port of the FESOM2 (AWI) ocean+sea-ice model, built and
exercised on Levante. Around 20 kLoC across `src/*.{c,h}` (≈17 k in
`.c`), mirroring the Fortran modules with the same loop bounds and halo
discipline.

The port reproduces the FESOM2 CORE2 reference configuration and has
been validated against the Fortran model over multi-year runs
(1958–1962 JRA55-do, dt=1800 s, 864 ranks): a full clean 5-year
integration with global SST/SSS RMS differences vs the Fortran climate
at the O(10⁻³) °C / PSU level. See **Validation status** below.

Currently included physics (out of the FESOM2 superset):

- **Ocean**: linfs ALE, JM-EOS, hydrostatic PGF, AB2 Coriolis,
  momentum advection (`momadv_opt=2`, `momentum_adv_scalar`),
  biharmonic horizontal viscosity (`opt_visc=7`, `visc_filt_bidiff` —
  CORE2's real scheme; `opt_visc=5` backscatter is also ported but
  unused), MFCT 3rd-order horizontal + FCT vertical tracer advection,
  CG SSH solver, shortwave penetration (`use_sw_pene`, Sweeney
  chlorophyll), JRA55-do forcing with NCAR bulk formulae and geo→rotated
  wind rotation, SSS restoring (local reference salinity) + runoff.
- **Vertical mixing**: **KPP (K-Profile) is the default scheme** — the
  CORE2 production choice. Pacanowski–Philander (PP) is also fully
  ported and selectable via `FESOM_MIX_SCHEME=PP`. Both include
  convective adjustment and reproduce the Fortran climate.
- **Sea ice**: EVP dynamics (`whichEVP=0`) with atmosphere–ice wind
  stress, thermodynamics with T-clamping at freezing, ice–ocean
  coupling (ocean2ice / oce_fluxes / oce_fluxes_mom), FCT advection of
  `a_ice / m_ice / m_snow`.
- **GM/Redi**: Ferrari et al. (2010) bolus velocity + isoneutral Redi
  diffusion, including ODM95 slope tapering and `scaling_GMzexp`
  vertical scaling. Master off-switch `FESOM_NO_GMREDI`.
- **I/O**: FESOM-style time-mean output streams with a model calendar
  (the deliberate non-port redesign — see below). Per-variable cadence
  override via `FESOM_IO_CONFIG`. Debug snapshots still available.

Not yet included: restart I/O (cold start only), `zlevel` / `zstar` ALE
(linfs only), and the rest of the FESOM2 superset not exercised by
CORE2.

The model is **Levante-only**. The build script and `env.sh` source the
gcc-11 / openmpi-4.1.2 / netcdf-c modules and set the OpenMPI runtime
knobs that DKRZ's `fesom2/env/levante.dkrz.de/shell.gnu` uses. Other
clusters would need their own equivalents.

## Build

```
cd /home/a/a270088/port2/fesom2_port
bash -l configure.sh           # incremental
bash -l configure.sh --clean   # wipe build/ and start over
```

Output: `build/fesom_port`. The `bash -l` invocation is required so the
`module load` calls in `env.sh` see Levante's profile. Don't try to
build inside a mambaforge/conda environment — it will pick up a wrong
mpicc/netcdf and fail at link time. If you've activated one in your
shell, `module --force purge` in `env.sh` clears it.

## Data on Levante

Hard-coded in the job templates. To run elsewhere, change the paths in
the `.sh` files.

| Asset | Path |
|---|---|
| CORE2 mesh | `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2` |
| PHC3.0 winter IC | `/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc` |
| JRA55-do forcing | reader hard-codes the `/pool/data/.../JRA55-do/` layout for the requested year |
| Output | `/work/ab0995/a270088/port/...` (or `--out_dir`) |

Mesh layout the model expects under `<mesh_dir>`:

- `nod2d.out`, `elem2d.out`, `aux3d.out`, `nlvls.out`, `elvls.out`,
  `edges.out`, `edge_tri.out`, `depth.out` — global mesh files.
- `dist_<NPES>/my_list*.out` + `com_info*.out` — partitioning files
  per rank count. The model picks `dist_<npes>` based on the SLURM
  `--ntasks` value. **Multi-rank runs need pre-generated partitions
  for the exact rank count you intend to use** — rank counts seen so
  far in this tree: 1, 8, 16, 32, 144, 256, 512, 864.

## Running

The binary signature:

```
fesom_port <mesh_dir> [output_dir] [dt_seconds] [nsteps] [snap_every] [phc_nc_path] [jra55_year]
```

Positional arguments:

| arg | required | default | notes |
|---|---|---|---|
| `mesh_dir` | yes | — | path to the CORE2-style mesh directory |
| `output_dir` | no | `.` | where the monthly streams (and optional `snap_*.nc`) land |
| `dt_seconds` | no | compiled `FESOM_PHASE1_DT` (2400 s) | **production runs pass `1800`** — the validated CORE2 timestep |
| `nsteps` | no | 500 | total steps to run (1 yr ≈ 17 520 @ dt=1800) |
| `snap_every` | no | 25 | `> 0` snapshot every N steps; `0` → default 25; **`< 0` disables snapshots** (use for production — rely on monthly streams) |
| `phc_nc_path` | no | none → analytical forcing | PHC initial condition file |
| `jra55_year` | no | 0 → no JRA forcing | first year of JRA55-do forcing (1958 is the validation start) |

Empty-string args (`""`) act as "use default" — useful when you want
to override only later args. The compiled-in `dt` default (2400 s) is
only a fallback; **always pass `1800` explicitly** for a run you intend
to compare against the reference.

The model calendar starts at `jra55_year` (Gregorian) so the JRA reader
and the output stream filenames are dated correctly; with JRA disabled
it defaults to 1958 for CF time stamps.

## Job templates

Ready-to-submit SLURM scripts live in [`jobs/`](jobs/) — see
[`jobs/README.md`](jobs/README.md) for the full list and how to adapt them:

| Script | Ranks | Length | Notes |
|---|---|---|---|
| `smoke_8rank.sh` | 8 | ~1 month | quick "does it build, run, write output?" check |
| `production_864_5yr.sh` | 864 | 5 yr | validated reference run, KPP (default) |
| `production_pp.sh` | 864 | 5 yr | same, with `FESOM_MIX_SCHEME=PP` |
| `validation_2yr.sh` | 864 | 2 yr | daily+monthly output for Fortran comparison |

Submit from the repo root with `sbatch jobs/<script>`. The 864-rank,
5-year production template (`jobs/production_864_5yr.sh`) is shaped like
this:

```bash
#!/bin/bash
#SBATCH --job-name=fesom_port
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=03:00:00
#SBATCH -A ab0995
#SBATCH -o /work/ab0995/a270088/port/run_<stamp>/slurm.%j.out
#SBATCH -e /work/ab0995/a270088/port/run_<stamp>/slurm.%j.err

source /sw/etc/profile.levante
source /home/a/a270088/port2/fesom2_port/env.sh
ulimit -s 204800

REPO=/home/a/a270088/port2/fesom2_port
BIN=$REPO/build/fesom_port
MESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
PHC=/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc
OUT=/work/ab0995/a270088/port/run_<stamp>
mkdir -p "$OUT"
( cd "$REPO" && git rev-parse HEAD ) > "$OUT/SHA.txt"   # record the build

export FESOM_MIX_SCHEME=KPP        # default; =PP to opt out
export FESOM_PRINT_EVERY=1000

# dt=1800, 5 calendar yr (1958-1962) = 87600 steps, snapshots disabled (-1)
srun -l --ntasks=864 "$BIN" "$MESH" "$OUT" 1800 87600 -1 "$PHC" 1958 \
     > "$OUT/run.log" 2> "$OUT/run.err"
echo "exit: $?"
tail -3 "$OUT/run.log"
ls "$OUT"/*.monthly.nc | wc -l
```

**Always use a unique `OUT_DIR` per job.** Concurrent jobs that share an
`OUT_DIR` clobber each other's logs and output and can cause spurious
"crashes" that are actually log corruption. The `_<stamp>` suffix earns
its keep. Recording `git rev-parse HEAD` into the run dir makes a run
reproducible after the fact.

Shorter smoke / validation jobs (1–256 ranks, a month or a year) follow
the same shape — see `jobs/smoke_8rank.sh` and `jobs/validation_2yr.sh`.

## Output

Two output paths, both anchored to the model calendar:

### Time-mean streams (primary)

One stream per cadence owns a set of variables; on each calendar
boundary it gathers and appends one record per variable to a file named
`<var>.fesom.<date>.<cadence>.nc` (monthly/daily/yearly roll up per
year, e.g. `temp.fesom.1958.monthly.nc` holds 12 records). The
**compiled-in default is monthly-only, 17 variables** (12 ocean + 5
sea-ice; the ice ones are dropped when sea ice isn't initialised):

`temp, salt, ssh, sst, sss, u, v, w, Kv, Av, density, bvfreq`
(ocean) + `a_ice, m_ice, m_snow, uice, vice` (ice).

`Kv`, `Av`, `bvfreq` and `w` are **interface** quantities written on the
`nl` (=48) levels; `temp/salt/density/u/v` are mid-layer (`nl-1`=47).

Override the per-variable cadence with a config file:

```
export FESOM_IO_CONFIG=/path/to/io.config
```

Grammar (see `docs/plans/completed/20260425-io-system.md` and
`src/fesom_io_config.h`): one line per `varlist  WSEP  cadlist`, e.g.

```
temp,salt    monthly,daily      # both vars at monthly AND daily
ssh          step               # every step
```

`cadence := step | hourly | daily | monthly | yearly`. A variable named
in the config **replaces** its default cadence list. An example file
lives at `io.config.daily_monthly`.

This IO subsystem is a **deliberate, authorised non-port redesign** —
do not "fix" it into a literal port of FESOM2's `io_meandata.F90`. See
the banner in `src/fesom_io_stream.h`.

### Debug snapshots (optional)

When `output_dir` is set and `snap_every > 0`, also writes
`<output_dir>/snap_NNNNNN.nc` every `snap_every` steps (NNNNNN = step,
zero-padded). Each holds mesh, full ocean state (`T, S, eta_n, w, u, v,
density_m_rho0, bvfreq, pgf_x, pgf_y, Kv, Av`), and sea ice when
configured. Production runs disable these with `snap_every = -1`.

Per-step text stats land in `run.log` (one line per `FESOM_PRINT_EVERY`
steps). Per-rank stderr lands in `run.err`.

## Environment knobs

All read once and cached at first use unless noted otherwise.

### Scheme selection

| Knob | What it does |
|---|---|
| `FESOM_MIX_SCHEME=PP` | Use Pacanowski–Philander vertical mixing instead of the default KPP. Anything not starting with `P`/`p` keeps KPP. |
| `FESOM_IO_CONFIG=<path>` | Per-variable output cadence override file (see Output). Unset → compiled monthly-only default. |
| `FESOM_CHL_SOURCE=<name>` | Shortwave-penetration chlorophyll source (default `Sweeney`). |
| `FESOM_CHL_FILE=<path>` | Override the chlorophyll climatology file. |

### Master switches for physics phases

| Knob | What it does |
|---|---|
| `FESOM_NO_GMREDI=1` | Skip ALL GM/Redi physics. Treats `ctx->gm` as `NULL` → no sigma_xy, init_Redi_GM, fer_solve_Gamma, fer_w accumulator, bolus add/sub, Redi K33, G7a/G7b. With this set the binary was **byte-identical to the pre-G HEAD** — the gating correctness invariant. |
| `FESOM_NO_ICE_DYN=1` | Skip the EVP dynamics step inside `fesom_ice_step`. |
| `FESOM_NO_ICE_ADV=1` | Skip ice tracer advection (TG_rhs + fct_solve). cut_off still runs. |
| `FESOM_NO_ICE_THERMO=1` | Skip ice thermodynamics + ice→ocean flux update. |

### Diagnostics / instrumentation

| Knob | What it does |
|---|---|
| `FESOM_PRINT_EVERY=N` | Print per-step stats every N steps (defaults to `snap_every`, or 1000 if snapshots are off). |
| `FESOM_VERBOSE_CG=1` | Print SSH CG iteration details. |
| `FESOM_EVP_DUMP_DIR=/path` | Sea-ice EVP debug dumps (per-substep state, one dir per rank) — see `reference_evp_dump_diagnostic` in project memory. Always-on companion: an unconditional stderr trace whenever `ice_strength > 1e7`. |
| `FESOM_KPP_DUMP_DIR=/path`, `FESOM_KPP_DUMP_STEP=N` | KPP per-step dumps for validation against the Fortran K-routines. |
| `FESOM_KPP_REPLAY_DIR=/path` | Inject dumped KPP inputs to isolate the algebra from upstream input noise (the controlled-replay technique). |
| `FESOM_DIAG_MICE / _GID / _SPREAD / _SSHSLV / _VISCEDGE=1` | Targeted stderr probes kept from debug sessions (ice-mass monthly trace, per-GID, field spread, SSH-solver, viscosity-at-edge). |

### Physics-bisect toggles (kept in tree from MPI debug session)

These zero out forcing terms after they've been computed, useful when
hunting for which subsystem causes a multi-rank divergence.

| Knob | What it does |
|---|---|
| `FESOM_NO_WIND=1` | Zero `stress_surf` and `stress_node_surf` after bulk + runoff write them. Also skips `fesom_ice_oce_fluxes_mom` so ice-ocean drag isn't re-introduced. |
| `FESOM_NO_HFLUX=1` | Zero `heat_flux`, `water_flux`, `virtual_salt`, `relax_salt`. |
| `FESOM_FREEZE_TS=1` | Restore T,S to the IC values at the start of every step (so dynamics evolves but tracers don't drift). |
| `FESOM_NO_TRADV=1` | Skip ocean tracer advection (T then S). |
| `FESOM_NO_TRDIFF=1` | Skip ocean tracer implicit vertical diffusion. |
| `FESOM_NO_SFLOOR=1` | Skip the salinity-floor clamp (debug). |
| `FESOM_VISC_MULT=<x>` | Scale horizontal viscosity by `x` (timestep-stability sensitivity probe). |

These are independent — set as many as you like in one run.

## Validation status

The port is validated against the Fortran FESOM2 CORE2 reference at
`dt=1800` s, 864 ranks, JRA55-do 1958+:

- **dt=1800 stability**: a full clean 5-year run (1958–1962, 87 600
  steps) with `uv_max` peaking at ~4.65 m/s (the `uv>5` guard never
  fires). The earlier dt=1200 blow-up at ~year 4.9 is resolved.
- **Climate match (KPP, default)**: C+KPP reproduces the Fortran+KPP
  climate over the 5-year run; drift is tiny and non-growing.
- **Climate match (PP)**: after the bvfreq N² horizontal-smoothing fix,
  C+PP matches Fortran+PP to global SST RMS ≈ 0.004 °C, SSS ≈ 0.002 PSU,
  stable over 3 years.
- **GM/Redi gating**: `FESOM_NO_GMREDI=1` was byte-identical to the
  pre-GM baseline (commit `30b59b3`, 0.0e+00 max diff on every variable
  at 16 ranks).

Validation figures and scripts live under `docs/` (`drift_5yr_d1800/`,
`validation_2yr_d1800/`, `kpp_5yr_fix_figures/`, …) and `scripts/`.

### Cross-rank consistency

Sea-ice (Phase E) was bit-identical at printed precision across rank
counts. With GM/Redi active there is a small partition-dependent signal
at hot spots driven by PHC IC partition-dependence (the per-rank
bilin+extrap, mirroring the Fortran — see `feedback_phc_rank_dependent`
in project memory). Mean diffs stay O(1e-4) relative.

```bash
# Compare snapshots across rank counts
python3 -c "
import netCDF4 as nc, numpy as np
a = nc.Dataset('run_1r/snap_000200.nc')
b = nc.Dataset('run_8r/snap_000200.nc')
for v in ['T','S','u','v']:
    d = np.abs(a.variables[v][0,...] - b.variables[v][0,...])
    print(f'{v} max|d|={d.max():.2e}  mean|d|={d.mean():.2e}')"
```

### Things to check at the end of a run

```bash
# 1. SLURM thinks it succeeded
sacct -j <jobid> -o JobID,State,ExitCode,Elapsed --noheader -P | head

# 2. No HUGE-ice_strength fires (always-on sea-ice sanity trace)
grep -c "HUGE ice_strength" run.err   # expect 0

# 3. No blow-up / uv>5 guard trips
grep -iE "blowup|uv>5|nan" run.err     # expect empty

# 4. Step-N ocean stats look physical
grep -E "[0-9]+[ ]+it=" run.log | tail -5
#  -> uv, T-range, S-range, stress, CG iters; T should hold -2.06 °C
#     (ice freezing clamp) where there's sea ice

# 5. Monthly streams were written
ls *.monthly.nc | wc -l
```

## Source layout

```
src/
  fesom_main.c               CLI parsing, init, timestep loop, exit
  fesom_step.{c,h}           one ALE+adv+diff+ice timestep; mix-scheme dispatch
  fesom_constants.{c,h}      compiled physics params (dt, viscosity, switches)
  fesom_mesh.{c,h}           mesh I/O + derived geometry (areas, gradient_sca)
  fesom_partit.{c,h}         partition + halo communicator state
  fesom_halo.{c,h}           generic halo exchange primitives
  fesom_mpi.{c,h}            MPI bootstrap
  fesom_calendar.{c,h}       model calendar (Gregorian/noleap/360-day)

  fesom_aux.{c,h}            ocean aux: density, hpressure, bvfreq (N2 smoothing),
                             sw_alpha/sw_beta, Kv, Av, pgf, MLD1_ind
  fesom_dyn.{c,h}            dynamic state: uv, w, eta, ssh_rhs, fer_uv, fer_w
  fesom_tracers.{c,h}        T, S (and tracer ID infra)
  fesom_forcing.{c,h}        stress_node_surf, stress_surf, fluxes,
                             forcing fields shared with bulk + ice
  fesom_forcing_analytical.{c,h}   analytical forcing (zonal-wind toy)
  fesom_sss_runoff.{c,h}     SSS restoring (local ref salinity) + runoff

  fesom_eos.{c,h}            JM-EOS, pressure_bv, sw_alpha_beta
  fesom_ale.{c,h}            linfs ALE: thickness, vert_vel (with fer_w
                             accumulator), CFLz, wsplit, commit
  fesom_momentum.{c,h}       compute_vel_rhs, momentum_adv_scalar,
                             visc_filt_bidiff (opt_visc=7) / bcksct (5),
                             impl_vert_visc, update_vel
  fesom_pp.{c,h}             PP vertical mixing + convective adjustment
  fesom_kpp.{c,h}            KPP (K-Profile) vertical mixing — default scheme
  fesom_ssh.{c,h}            SSH stiffness build + parallel CG solver
  fesom_tracer_adv.{c,h}     MFCT 3rd-order horizontal + FCT vertical advection
  fesom_tracer_diff.{c,h}    implicit vert diffusion + Redi K33
  fesom_io.{c,h}             stream orchestrator, var registry, snap_*.nc writer
  fesom_io_stream.{c,h}      one stream = N vars at one cadence (accumulate + write)
  fesom_io_stream_dispatch.c pure sizing/routing helpers (unit-testable)
  fesom_io_config.{c,h}      FESOM_IO_CONFIG cadence-override file parser
  fesom_io_gather.h          gather-plan API used by the streams

  fesom_jra55.{c,h}          JRA55-do reader, bilin to mesh, g2r wind rotation
  fesom_bulk.{c,h}           NCAR bulk formulae (open-water fluxes)
  fesom_phc.{c,h}            PHC IC reader + cyclic-seam bilin + extrap_nod3D
  fesom_ic.{c,h}             hnode/helem from depth, T,S from PHC

  fesom_ice.{c,h}            sea-ice driver (ocean2ice + EVP + FCT +
                             cut_off + thermo + oce_fluxes)
  fesom_ice_types.h          T_ICE struct mirror
  fesom_ice_evp.{c,h}        EVP dynamics (whichEVP=0) + stress_atmice
  fesom_ice_thermo.{c,h}     thermodynamics, cut_off, runoff fold
  fesom_ice_coupling.{c,h}   ocean2ice, oce_fluxes (heat/water/virt/relax),
                             oce_fluxes_mom (ice-mediated stress)
  fesom_ice_fct.{c,h}        ice tracer FCT advection (mass_matrix_fill,
                             TG_rhs, fct_solve, solve_low/high, fem_fct)

  fesom_gm.{c,h}             GM/Redi: sigma_xy, neutral_slope, init_Redi_GM,
                             fer_solve_Gamma, fer_gamma2vel,
                             diff_ver_part_redi_expl, diff_part_hor_redi
```

## Documentation in tree

- `docs/plans/completed/` — finished phase plans (MPI port, GM/Redi,
  IO system, KPP vertical mixing).
- `docs/plans/<active>.md` — in-progress / recent plans (MFCT,
  shortwave penetration).
- `docs/PORTING_LESSONS.md`, `docs/PORT_EXPERIENCE_REPORT.md`,
  `docs/MPI_PORT_REPORT.md` — narrative writeups of the gotchas paid for.
- `docs/DT1800_HANDOFF.md`, `docs/FORCING_STEP1_DIFFERENCE.md` —
  deep-dive records on the dt=1800 stabilization and the (resolved)
  step-1 forcing diff.
- Project memory under `~/.claude/projects/-home-a-a270088-port2/memory/`
  has feedback files for the recurring bug patterns:
  - `feedback_write_loops_halo.md` — Fortran loops `myDim+eDim` on a
    write must port as `myDim+eDim` in C, not `myDim`.
  - `feedback_array_size_vs_reader_loop.md` — forcing arrays sized
    `myDim` only that get read at halo by another module are silent
    OOB reads (multi-rank-only bug).
  - `feedback_tracer_stride_nl.md` — tracer arrays are `[N * nl]`,
    not `[N * (nl-1)]`. Wrong stride inflates derived gradients ~1000×.
  - `feedback_bvfreq_smoothing_gap.md` — N² needs horizontal smoothing
    (`smooth_nod3D`); missing it dominated the old PP↔Fortran gap.
  - `feedback_mfct_gradient_from_values.md` — MFCT gradient is built
    from `values`, not `valuesAB`.
  - `feedback_magnitude_invariant_masks_rotation.md` — a wrong-angle
    wind-stress bug passes every scalar/|τ| check; compare the vector.
  - `feedback_cyclic_seam_interp.md` — the PHC longitude axis is
    periodic; the wrap halo must be padded correctly.
  - `feedback_bolus_divergence_balance.md` — never clamp `fer_uv` per
    cell without scaling `fer_gamma` to match.
  - `feedback_io_system_design.md` — the IO redesign is the authorised
    non-port exception; don't replace it with a literal port.
  - `feedback_unique_outdir_per_job.md` — concurrent SLURM jobs must
    not share `OUT_DIR`.
  - `feedback_phc_rank_dependent.md` — PHC IC differs across rank counts
    by design (extrap_nod3D Gauss-Seidel ordering); Fortran-faithful.
  - `feedback_levante_build.md` — never use mambaforge; build via
    `bash -l configure.sh`.

## Troubleshooting

**Build fails with linker errors against MPI / NetCDF.** You're not in
a clean shell. Quit the shell, log back in, and start again with
`bash -l configure.sh`. Check `which mpicc` returns
`/sw/spack-levante/openmpi-...`, not `~/mambaforge/bin/...`.

**Blow-up / `uv>5` guard trips, especially at dt=1800.** The dt=1800
stability rests on several fixes that must all be in place: the AB2
stabilization-offset epsilon = 0.1 (not 1e-9), biharmonic viscosity
(`opt_visc=7`), `use_wsplit=0`, and the bulk Ch/Ce halo coverage. If a
new change reintroduces a blow-up, bisect with `FESOM_VISC_MULT` and the
`FESOM_NO_*` toggles.

**Multi-rank crash early in step 2 with CG NaN.** Most likely a halo
coverage gap on a forcing or aux field that some kernel reads at halo.
Bisect with `FESOM_NO_TRADV=1`, then `_TRDIFF`, `_NO_WIND`, `_NO_HFLUX`
to localize the amplifier. See `feedback_array_size_vs_reader_loop` and
`feedback_write_loops_halo`.

**`run.log shows S.max > 41.12 PSU at step 1.`** Tracer advection
overshot — a divergence-balance issue in the bolus add
(`feedback_bolus_divergence_balance`) or a stride bug in a gradient
(`feedback_tracer_stride_nl`).

**`HUGE ice_strength` lines in stderr.** Sea-ice EVP detected
`ice_strength > 1e7` somewhere — almost always a halo node with bogus
ice mass from upstream OOB. Set
`FESOM_EVP_DUMP_DIR=/work/.../evp_$(date +%s)` to capture per-substep
dumps and use `scripts/evp_dump_diff.py` to compare against a 1-rank
reference.

**Job hangs at startup.** Wrong `dist_<npes>` partition for your
`--ntasks`. Match `--ntasks` to a partition directory that exists under
the mesh.

**`bash: module: command not found` in a job script.** Missing
`source /sw/etc/profile.levante` before `source env.sh`.

## Reference Fortran source

Used throughout for cross-checking: `/home/a/a270088/port2/fesom2/src/`.
Every function in this port has a `Fortran <file>:<line>` comment
pointing back at the producer routine. Search for `Mirror of` or
`Fortran` to navigate.
