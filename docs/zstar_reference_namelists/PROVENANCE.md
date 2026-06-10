# zstar reference runs — provenance (Phase Z0, 2026-06-10)

Fortran tree: `/home/a/a270088/port2/fesom2` at HEAD `9271ae92`
("instrumentation: env-gated reference dumps for the FESOM2->C port comparison").
Intel oneapi build via `env/levante.dkrz.de/shell`.

## The three reference runs

| Job | Work dir | Output (purge-safe /work) | Config | Binary |
|---|---|---|---|---|
| 25492898 `ftn_zstar2y` | `fesom2/work_zstar_kpp` | `/work/ab0995/a270088/port/zstar/fortran_zstar_2yr/` | 2yr, dt=1800, 864r, KPP, **which_ALE='zstar'** | frozen copy (below) |
| 25492899 `ftn_linfs2yb` | `fesom2/work_linfs_2yr_b` | `/work/ab0995/a270088/port/zstar/fortran_linfs_2yr_b/` | 2yr, dt=1800, 864r, KPP, which_ALE='linfs' | frozen copy (below) |
| 25493084 `zstar_fdump` | `fesom2/work_zstar_dump` | `/work/ab0995/a270088/port/zstar/fdump/` (dumps in `fdump/dump/`) | 3 steps, dt=1800, 16r, KPP, zstar; `FESOM_ALE_DUMP_DIR` + `FESOM_ALE_DUMP_STEPS=3` | live `bin/fesom.x` with the (uncommitted) `ale_dump_mod` gate |

Namelists archived in this directory = `work_zstar_kpp`'s set (all 10 files).

## One-knob guarantee

`diff work_linfs_2yr_b/ work_zstar_kpp/` over all 10 namelists shows exactly two
lines: `ResultPath` and `which_ALE = 'linfs' → 'zstar'`. Everything else (KPP,
dt=1800, MFCT, GM/Redi, EVP, JRA55, PHC, opt_visc=7, use_wsplit=F,
use_ssh_se_subcycl=F) is identical between the pair — and identical to the
original `work_linfs_2yr` (v1.0-era reference) except ResultPath.

## Binary provenance

- **Both 2yr runs** execute a **frozen copy**
  (`fesom2/frozen_zstar_refs/bin/fesom.x` + `lib64/libfesom.so`,
  sha256 `700e2d19…629e`) taken from HEAD `9271ae92` BEFORE the ALE dump gate
  was added. Work-dir `fesom.x` symlinks point at the frozen copy, so later
  rebuilds of `bin/`/`lib64/` cannot perturb them (both jobs were observed
  RUNNING with the frozen lib before the dump-gate install at 18:51).
- **The dump run** uses the live `bin/fesom.x` (lib sha256 `3afe8268…1a09`,
  built 2026-06-10 18:51) = HEAD `9271ae92` + the **uncommitted** `ale_dump_mod`
  in `oce_ale.F90` (env-gated; no-op when `FESOM_ALE_DUMP_DIR` is unset),
  following the KPP/EVP uncommitted-dump-gate convention.

## Differences vs the ORIGINAL `work_linfs_2yr` run (May 21)

The original 2yr linfs output (`/scratch/a/a270088/fortran_2yr_dt1800`) was
**purged from scratch** — hence the rerun. The May 21 binary predates two
committed Fortran changes that are in the frozen pair binary:
1. `2682a9fb` — sbc: rotate initial wind/stress coefficients on rotated grids
   (the step-1 forcing cold-start fix, `project_forcing_step1_diff`; the C was
   always correct, so the new pair REMOVES a known C-vs-Fortran transient).
2. `49225df7`/`9271ae92` — env-gated JAX dump shims (no-op when env unset).

## Reference config values re-confirmed (plan table)

- `which_ALE='zstar'` (namelist.config:71, the one knob)
- `use_partial_cell=.false.` (namelist.config:75)
- `use_floatice=.false.` (namelist.config:107) → EVP/vel_rhs `use_pice=0`
- `which_pgf` NOT set in namelist.oce → module default `'shchepetkin'`
  (oce_modules.F90)
- `i_vert_diff=.true.` (namelist.tra:74) → implicit TDMA carries surface BCs
- `use_ssh_se_subcycl=.false.` (namelist.dyn:60) → CG path
- `use_wsplit=.false.`, `opt_visc=7` (namelist.dyn) — as in the linfs pair
- `alpha`/`theta` not in namelist.dyn → module defaults, same as v1.0 pair
- `use_virt_salt`/`is_nonlinfs` are DERIVED in oce_setup_step.F90 from
  which_ALE — zstar ⇒ use_virt_salt=.false., is_nonlinfs=1.0

## Fortran ALE dump gate (uncommitted, Intel tree)

`oce_ale.F90` top-of-file `module ale_dump_mod` + 6 one-line call sites in
`oce_timestep_ale`:

| Tag | After | Fields (ncomp) |
|---|---|---|
| `forcing` | step entry | water_flux, virtual_salt, relax_salt, real_salt_flux (4, nodes) |
| `pgf_x`/`pgf_y` | pressure_force_4_* | nl-1 levels (elems) |
| `sshsolve` | solve_ssh_ale | ssh_rhs, d_eta (2, nodes) |
| `hbar` + `dhe` | eta_n update | hbar, hbar_old, ssh_rhs_old, eta_n (4, nodes); dhe (1, elems) |
| `Wvel`, `hnode_new` | vert_vel_ale | nl / nl-1 levels (nodes) |
| `hnode`, `zbar_3d_n`, `Z_3d_n`, `helem` | update_thickness_ale | nl-1 / nl / nl-1 (nodes); nl-1 (elems) |

Format identical to the KPP dumps: `ale_dump_s<step>_<tag>_rank<R>.txt`,
header `# step= tag= rank= N= ncomp=`, gid-keyed rows, `es24.16`, OWNED
entities only. Reader: `scripts/ale_dump_diff.py`.
