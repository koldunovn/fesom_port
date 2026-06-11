# TKE reference runs — provenance (Phase T0, 2026-06-11)

Fortran tree: `/home/a/a270088/port2/fesom2` at HEAD `9271ae92`
("instrumentation: env-gated reference dumps for the FESOM2->C port comparison").
Intel oneapi build via `env/levante.dkrz.de/shell`. `CVMIX:BOOL=ON` in the build
cache; `nm` confirms `g_cvmix_tke_mp_{init,calc}_cvmix_tke_` in BOTH binaries
below.

## The two reference runs

| Job | Work dir | Output (purge-safe /work) | Config | Binary |
|---|---|---|---|---|
| 25500644 `ftn_linfstke` | `fesom2/work_linfs_tke` | `/work/ab0995/a270088/port/tke/fortran_linfs_tke/` | 2yr 1958-59, dt=1800, 864r, **mix_scheme='cvmix_TKE'**, which_ALE='linfs' | frozen copy (below) |
| 25500775 `tke_fdump` | `fesom2/work_tke_dump` | `/work/ab0995/a270088/port/tke/fdump/` (dumps in `fdump/dump/`) | 3 steps, dt=1800, 16r, cvmix_TKE, linfs; `FESOM_TKE_DUMP_DIR` + `FESOM_TKE_DUMP_STEPS=3` | live `bin/fesom.x` with the (uncommitted) `tke_dump_mod` gate |

Namelists archived in this directory = `work_linfs_tke`'s set (all 10 files).

## One-knob guarantee

`diff work_linfs_2yr_b/ work_linfs_tke/` over all 10 namelists shows exactly two
lines: `ResultPath` and `mix_scheme = 'KPP' → 'cvmix_TKE'`. Everything else
(linfs, dt=1800, MFCT, GM/Redi, EVP, JRA55, PHC, opt_visc=7, use_wsplit=F,
use_ssh_se_subcycl=F) is identical to the zstar-era paired baseline
`work_linfs_2yr_b` (job 25492899, the linfs+KPP 2yr).

`diff work_linfs_tke/ work_tke_dump/` shows exactly: `run_length 2y → 3s` and
`ResultPath`.

## Binary provenance

- **The 2yr run** executes the **frozen copy**
  (`fesom2/frozen_zstar_refs/bin/fesom.x` sha256 `69a35d75…`, +
  `lib64/libfesom.so` sha256 `700e2d19…`, RPATH `$ORIGIN/../lib64`) — the SAME
  binary as the zstar/linfs 2yr pair → C-TKE-vs-F-TKE and C-KPP-vs-F-KPP
  comparisons share one Fortran build. `work_linfs_tke/fesom.x` symlinks the
  frozen copy.
- **The dump run** uses live `bin/fesom.x` (lib sha256 `3b711e99…`, built
  2026-06-11) = HEAD `9271ae92` + the **uncommitted** `tke_dump_mod` in
  `gen_modules_cvmix_tke.F90` (env-gated; no-op when `FESOM_TKE_DUMP_DIR`
  unset) + the older uncommitted `ale_dump_mod` in `oce_ale.F90` (inactive,
  env unset).

## namelist.cvmix &param_tke values consumed (echo verified in BOTH run logs)

```
tke_only       = T          (default; nmb=5 keeps it true)
tke_c_k        = 0.1
tke_c_eps      = 0.7
tke_alpha      = 30.0
tke_cd         = 3.75       (NAMELIST value; module default is 1.0 — namelist wins)
tke_surf_min   = 1.0e-4
tke_min        = 1.0e-6
tke_kappaM_min = 0.0
tke_kappaM_max = 100.0
tke_mxl_choice = 2
tke_dolangmuir = F
use_ubound_dirichlet / use_lbound_dirichlet = F / F (defaults, not in namelist)
```

The `init_cvmix_tke` echo block (gen_modules_cvmix_tke.F90:228-241) appearing in
the run log proves (a) `namelist.cvmix` was actually read (silent cd=1.0
fallback otherwise) and (b) `__cvmix` is compiled in (otherwise nmb=5 silently
computes NO mixing — a plausible-looking broken reference). Verified for job
25500644 at launch (2026-06-11): `tke_cd = 3.75`, `tke_only = T`.

## TKE dump gate (uncommitted, Intel tree)

`gen_modules_cvmix_tke.F90` top-of-file `module tke_dump_mod` + capture
buffers in `calc_cvmix_tke`. Per TKE call (= ocean step) 1..3, OWNED entities:

| Tag | What | ncomp |
|---|---|---|
| `normstress` | √(stress_node_surf²)/ρ₀ per node | 1 |
| `vshear2` / `bvfreq2` / `dztrr` / `tkeold` | per-column INPUTS at the cvmix call | nl (nodes) |
| `tke` / `tkeav` / `tkekv` | post-loop outputs (post-zeroing, pre-exchange) | nl (nodes) |
| `tbpr tspr tdif tdis twin tiwf tbck ttot lmix pr` | the 10 budget diagnostics | nl (nodes) |
| `kv` | final Kv after `Kv = tke_Kv` | nl (nodes) |
| `av` | final element Av after the node→elem mean | nl (elems) |

Format identical to the KPP/ALE dumps: `tke_dump_s<step>_<tag>_rank<R>.txt`,
header `# step= tag= rank= N= ncomp=`, gid-keyed rows, `es24.16`. Reader:
`scripts/tke_dump_diff.py`.
