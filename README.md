# FESOM2 C port

A literal C port of the FESOM2 (AWI) ocean+sea-ice model, built and
exercised on Levante. Around 15 kLoC across `src/*.{c,h}`, mirroring
the Fortran modules with the same loop bounds and halo discipline.

Currently included physics (out of the FESOM2 superset):

- Ocean: linfs ALE, JM-EOS, hydrostatic PGF, AB2 Coriolis, FCT tracer
  advection, opt_visc=5 backscatter, PP mixing + convective adjustment,
  CG SSH solver, JRA55 forcing with NCAR bulk formulae, SSS restoring
  + runoff.
- Sea ice (Phases A–E): EVP dynamics (whichEVP=0), thermodynamics with
  T-clamping at freezing, ice-ocean coupling (ocean2ice / oce_fluxes /
  oce_fluxes_mom), FCT advection of `a_ice / m_ice / m_snow`.
- GM/Redi (Phases G0–G9): Ferrari et al. (2010) bolus velocity +
  isoneutral Redi diffusion, including ODM95 slope tapering and
  scaling_GMzexp vertical scaling. Master off-switch `FESOM_NO_GMREDI`.

Not yet included: KPP mixing, GM/Redi at multi-decadal scales (zlevel /
zstar ALE), restart I/O, time-mean output streams.

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
  far in this tree: 1, 8, 16, 32, 144, 256.

## Running

The binary signature:

```
fesom_port <mesh_dir> [output_dir] [dt_seconds] [nsteps] [snap_every] [phc_nc_path] [jra55_year]
```

Positional arguments:

| arg | required | default | notes |
|---|---|---|---|
| `mesh_dir` | yes | — | path to the CORE2-style mesh directory |
| `output_dir` | no | none → no snapshot output | where `snap_*.nc` lands |
| `dt_seconds` | no | compiled FESOM_PHASE1_DT (500 s) | overrides timestep |
| `nsteps` | no | 500 | total steps to run |
| `snap_every` | no | 25 | step interval between NetCDF snapshots |
| `phc_nc_path` | no | none → analytical forcing | PHC initial condition file |
| `jra55_year` | no | 0 → no JRA forcing | year to use for JRA55-do (1958 was used for sea-ice/GM validation) |

Empty-string args (`""`) act as "use default" — useful when you want
to override only later args.

## Sample SLURM job (16 ranks, 200 steps)

Most copies in tree look like `job_core2_16_*`. Bare-bones template:

```bash
#!/bin/bash
#SBATCH --job-name=fesom_port_16
#SBATCH -p compute
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH -A ab0995
#SBATCH -o /work/ab0995/a270088/port/run_<stamp>/slurm.%j.out
#SBATCH -e /work/ab0995/a270088/port/run_<stamp>/slurm.%j.err

source /sw/etc/profile.levante
source /home/a/a270088/port2/fesom2_port/env.sh
ulimit -s 204800

OUT_DIR=/work/ab0995/a270088/port/run_<stamp>
mkdir -p "$OUT_DIR"

BIN=/home/a/a270088/port2/fesom2_port/build/fesom_port
MESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
PHC=/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc

DT=500
NSTEPS=200
SNAP_EVERY=50
export FESOM_PRINT_EVERY=50
JRA55_YEAR=1958

srun -l "$BIN" "$MESH" "$OUT_DIR" "$DT" "$NSTEPS" "$SNAP_EVERY" "$PHC" "$JRA55_YEAR" \
    > "$OUT_DIR/run.log" 2> "$OUT_DIR/run.err"
echo "Exit: $?"
tail -3 "$OUT_DIR/run.log"
```

**Always use a unique `OUT_DIR` per job.** Concurrent jobs that share an
`OUT_DIR` clobber each other's logs and can cause spurious "crashes"
that are actually log corruption. The `_<stamp>` suffix (e.g.
`20260425_204500`) earns its keep.

For 256-rank month/year runs see `job_core2_256_month_e` and
`job_core2_256_year_e` — same shape, longer wall, more snapshots.

## Environment knobs

All read once and cached at first use unless noted otherwise.

### Master switches for new physics phases

| Knob | What it does |
|---|---|
| `FESOM_NO_GMREDI=1` | Skip ALL GM/Redi physics. Treats `ctx->gm` as `NULL` for the rest of the step → no compute_sigma_xy, no init_Redi_GM, no fer_solve_Gamma, no fer_w accumulator, no bolus add/sub, no Redi K33, no G7a/G7b. With this set the binary is **byte-identical to the pre-G HEAD** — the gating correctness invariant. |
| `FESOM_NO_ICE_DYN=1` | Skip the EVP dynamics step inside `fesom_ice_step`. |
| `FESOM_NO_ICE_ADV=1` | Skip ice tracer advection (TG_rhs + fct_solve). cut_off still runs. |
| `FESOM_NO_ICE_THERMO=1` | Skip ice thermodynamics + ice→ocean flux update. |

### Diagnostics / instrumentation

| Knob | What it does |
|---|---|
| `FESOM_PRINT_EVERY=N` | Print per-step stats every N steps (defaults to `snap_every`). |
| `FESOM_VERBOSE_CG=1` | Print SSH CG iteration details. |
| `FESOM_EVP_DUMP_DIR=/path` | Sea-ice EVP debug dumps (per-substep state). One dir per rank — see `reference_evp_dump_diagnostic` in project memory. Always-on companion: an unconditional stderr trace whenever `ice_strength > 1e7`. |

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

These are independent — set as many as you like in one run.

## Output

`<output_dir>/snap_NNNNNN.nc` files, one per `snap_every` (NNNNNN is
the step number, zero-padded to 6 digits). Each contains:

- Mesh: `lon`, `lat`, `zbar`, `Z`, `elem_nodes`, `nlevels_nod2D`,
  `nlevels`.
- Ocean: `T(time, nz_1, nod2)`, `S(...)`, `eta_n(time, nod2)`,
  `w(time, nz, nod2)`, `u(time, nz_1, elem)`, `v(...)`,
  `density_m_rho0`, `bvfreq`, `pgf_x`, `pgf_y`, `Kv`, `Av`.
- Sea ice (when configured): `a_ice`, `m_ice`, `m_snow`, `uice`,
  `vice`, `h_ice`, `h_snow`.
- Globals: `:step`, `:dt`.

Per-step text stats land in `run.log` (one line per `FESOM_PRINT_EVERY`
steps). Per-rank stderr lands in `run.err`.

## Validating a run

Things to check at the end:

```bash
# 1. SLURM thinks it succeeded
sacct -j <jobid> -o JobID,State,ExitCode,Elapsed --noheader -P | head

# 2. No HUGE-ice_strength fires (always-on sea-ice sanity trace)
grep -c "HUGE ice_strength" run.err   # expect 0

# 3. Step-N ocean stats look physical
grep -E "[0-9]+[ ]+it=" run.log | tail -5
#  -> uv, T-range, S-range, stress, CG iters; T should hold -2.06 °C
#     (ice freezing clamp) once sea ice lands

# 4. Compare snapshots across rank counts (Python)
python3 -c "
import netCDF4 as nc, numpy as np
a = nc.Dataset('run_1r/snap_000200.nc')
b = nc.Dataset('run_8r/snap_000200.nc')
for v in ['T','S','u','v']:
    d = np.abs(a.variables[v][0,...] - b.variables[v][0,...])
    print(f'{v} max|d|={d.max():.2e}  mean|d|={d.mean():.2e}')"
```

For a fresh GM/Redi-on smoke at 16 ranks expect roughly:

- exit 0, ~3 min wall
- `it=` (CG iters) ~30–35 throughout
- T.min clamped at −2.06 °C, T.max ~30 °C, S in [0, 41.12]
- `stress` ~1–2 N/m², `pgf` ~1.4e-4 m/s²

Cross-rank variation: Phase E sea-ice was bit-identical at printed
precision; with GM/Redi active there is a small partition-dependent
signal at hot spots driven by PHC IC partition-dependence
(see `feedback_phc_rank_dependent` memory). Mean diffs stay O(1e-4)
relative.

## Bit-identity off-switch

The strongest correctness test for the GM/Redi gating:

```bash
export FESOM_NO_GMREDI=1
sbatch job_core2_16_g9_off
# After it runs, confirm byte-identical T/S/u/v/w vs the pre-GM HEAD
# baseline (G4 commit 704d352).
```

Tested at 16 ranks on commit `30b59b3` — 0.0e+00 max diff vs baseline
on every variable.

## Source layout

```
src/
  fesom_main.c               CLI parsing, init, timestep loop, exit
  fesom_step.{c,h}           one ALE+adv+diff+ice timestep
  fesom_mesh.{c,h}           mesh I/O + derived geometry (areas, gradient_sca)
  fesom_partit.{c,h}         partition + halo communicator state
  fesom_halo.{c,h}           generic halo exchange primitives
  fesom_mpi.{c,h}             MPI bootstrap

  fesom_aux.{c,h}            ocean aux: density, hpressure, bvfreq,
                              sw_alpha/sw_beta, Kv, Av, pgf, MLD1_ind
  fesom_dyn.{c,h}            dynamic state: uv, w, eta, ssh_rhs, fer_uv, fer_w
  fesom_tracers.{c,h}        T, S (and tracer ID infra)
  fesom_forcing.{c,h}        stress_node_surf, stress_surf, fluxes,
                              forcing fields shared with bulk + ice
  fesom_forcing_analytical.{c,h}   analytical forcing (zonal-wind toy)

  fesom_eos.{c,h}            JM-EOS, pressure_bv, sw_alpha_beta
  fesom_ale.{c,h}             linfs ALE: thickness, vert_vel (with fer_w
                              accumulator), CFLz, wsplit, commit
  fesom_momentum.{c,h}       compute_vel_rhs, visc_filt_bcksct,
                              impl_vert_visc, update_vel
  fesom_pp.{c,h}             PP vertical mixing + convective adjustment
  fesom_ssh.{c,h}             SSH stiffness build + parallel CG solver
  fesom_tracer_adv.{c,h}     FCT tracer advection (T then S)
  fesom_tracer_diff.{c,h}    implicit vert diffusion + Redi K33
  fesom_io.{c,h}              snap_*.nc gather + serial-NetCDF write

  fesom_jra55.{c,h}          JRA55-do reader, bilin to mesh, time interp
  fesom_bulk.{c,h}            NCAR bulk formulae (open-water fluxes)
  fesom_phc.{c,h}             PHC IC reader + bilin + extrap_nod3D
  fesom_ic.{c,h}              hnode/helem from depth, T,S from PHC

  fesom_ice.{c,h}             sea-ice driver (ocean2ice + EVP + FCT +
                              cut_off + thermo + oce_fluxes)
  fesom_ice_types.h           T_ICE struct mirror
  fesom_ice_evp.{c,h}        EVP dynamics (whichEVP=0)
  fesom_ice_thermo.{c,h}     thermodynamics, cut_off, runoff fold
  fesom_ice_coupling.{c,h}   ocean2ice, oce_fluxes (heat/water/virt/relax),
                              oce_fluxes_mom (Phase D5 ice-mediated stress)
  fesom_ice_fct.{c,h}        ice tracer FCT advection (mass_matrix_fill,
                              TG_rhs, fct_solve, solve_low/high, fem_fct)

  fesom_gm.{c,h}              GM/Redi: sigma_xy, neutral_slope, init_Redi_GM,
                              fer_solve_Gamma, fer_gamma2vel,
                              diff_ver_part_redi_expl, diff_part_hor_redi
```

## Documentation in tree

- `docs/plans/completed/` — finished phase plans
  (`20260425-sea-ice-port.md`, `20260425-gm-redi-port.md`, …).
- `docs/plans/<active>.md` — in-progress plan, if any.
- Project memory under `~/.claude/projects/-home-a-a270088-port2/memory/`
  has feedback files for the gotchas we've paid for:
  - `feedback_write_loops_halo.md` — Fortran loops `myDim+eDim` on a
    write must port as `myDim+eDim` in C, not `myDim`.
  - `feedback_array_size_vs_reader_loop.md` — forcing arrays sized
    `myDim` only that get read at halo by another module are silent
    OOB reads (multi-rank-only bug).
  - `feedback_tracer_stride_nl.md` — tracer arrays are `[N * nl]`,
    not `[N * (nl-1)]`. Wrong stride inflates derived gradients ~1000×.
  - `feedback_unported_consumer_gaps.md` — alloc + exchange ≠ field
    populated. New consumers expose unported producers.
  - `feedback_bolus_divergence_balance.md` — never clamp `fer_uv` per
    cell without scaling `fer_w` to match.
  - `feedback_unique_outdir_per_job.md` — concurrent SLURM jobs must
    not share `OUT_DIR`.
  - `feedback_phc_rank_dependent.md` — PHC IC differs O(1°C)/O(10 PSU)
    across rank counts due to extrap_nod3D's Gauss-Seidel ordering.
    Accepted as Fortran-faithful baseline.
  - `feedback_levante_build.md` — never use mambaforge; build via
    `bash -l configure.sh`.

## Troubleshooting

**Build fails with linker errors against MPI / NetCDF.** You're not in
a clean shell. Quit the shell, log back in, and start again with
`bash -l configure.sh`. Check `which mpicc` returns
`/sw/spack-levante/openmpi-...`, not `~/mambaforge/bin/...`.

**`Multi-rank crash early in step 2 with CG NaN.`** Most likely cause is
a halo coverage gap on a forcing or aux field that some kernel reads
at halo. Workflow: try the bisect knobs in order
(`FESOM_NO_TRADV=1`, then `_TRDIFF`, then `_NO_WIND`, then
`_NO_HFLUX`) to localize the amplifier. The
`feedback_array_size_vs_reader_loop` and `feedback_write_loops_halo`
memories capture the pattern.

**`run.log shows S.max > 41.12 PSU at step 1.`** Tracer FCT has
overshot. With GM/Redi this means either a divergence-balance issue
in the bolus add (`feedback_bolus_divergence_balance`) or a stride
bug in a gradient computation (`feedback_tracer_stride_nl`).

**`HUGE ice_strength` lines in stderr.** Sea-ice EVP detected
`ice_strength > 1e7` somewhere — almost always a halo node has bogus
ice mass from upstream OOB. Set `FESOM_EVP_DUMP_DIR=/work/.../evp_$(date +%s)`
to capture per-substep dumps and use `scripts/evp_dump_diff.py` to
compare against a 1-rank reference.

**Job hangs at startup.** Wrong `dist_<npes>` partition for your
`--ntasks`. Match `--ntasks` to one of the partition directories that
exists under the mesh.

**`bash: module: command not found` in a job script.** Missing
`source /sw/etc/profile.levante` before `source env.sh`.

## Reference Fortran source

Used throughout for cross-checking: `/home/a/a270088/port2/fesom2/src/`.
Every function in this port has a `Fortran <file>:<line>` comment
pointing back at the producer routine. Search for `Mirror of` or
`Fortran` to navigate.
