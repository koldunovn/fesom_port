# mEVP reference run provenance (M0, 2026-06-11)

These are the namelists of the Fortran 2yr mEVP reference run (`work_mevp`),
archived per plan `docs/plans/20260611-mevp-sea-ice-rheology.md` Phase M0.

## Run setup

- **Work dir:** `/home/a/a270088/port2/fesom2/work_mevp/`
- **SLURM job:** 25514027 (submitted 2026-06-11), 864 ranks (dist_864),
  `job_mevp_2yr_864`.
- **Gate:** run log echoes `EVP scheme option=           1`
  (fesom_module.F90:323) — confirmed at job start.
- **ResultPath:** `/work/ab0995/a270088/port/mevp/fortran_mevp_2yr/`
- **Config:** linfs + KPP + whichEVP=1, dt=1800 (step_per_day=48), 2yr
  (730 days = 35040 steps), cold start 1958, JRA55-do forcing.

## Base config lineage (the "single knob" claim)

`work_mevp` is a clone of `work_linfs_2yr_b` (the Fortran **standard-EVP**
linfs+KPP 2yr reference, run 2026-06-10 for the zstar/TKE matrix, output at
`/work/ab0995/a270088/port/zstar/fortran_linfs_2yr_b/`) with exactly TWO edits:

1. `namelist.ice`: `whichEVP = 0` → `1` — **the only physics change**;
2. `namelist.config`: `ResultPath` → unique output dir
   (`feedback_unique_outdir_per_job`).

`work_linfs_2yr_b`'s physics namelists (ice/forcing/tra/oce) were
diff-verified IDENTICAL to `fesom27/fesom2/work_core` (the plan's named
config source); its namelist.config pins linfs + run_length=2y. So this run
IS "work_core config with only whichEVP flipped" in the plan's sense, while
also being single-knob against the existing validated F-EVP leg — which is
the F-EVP side of the M5 diff-of-diffs.

mEVP parameters confirmed in `namelist.ice` (all equal to the MOD_ICE.F90
module defaults — no namelist-over-default trap): `alpha_evp=250`,
`beta_evp=250`, `evp_rheol_steps=120`, `Pstar=30000`, `ellipse=2.0`,
`c_pressure=20.0`, `delta_min=1.0e-11`, `Cd_oce_ice=0.0055`,
`ice_ave_steps=1` (→ ice_dt = dt = 1800 s; `rdt` = FULL step in mEVP).
`theta_io=0.0` is present but **not read** by the mEVP path.

## Binary provenance

`work_mevp/fesom.x` →
`/home/a/a270088/port2/fesom2/frozen_zstar_refs/bin/fesom.x`, whose RPATH
(`$ORIGIN/../lib64`) pins `frozen_zstar_refs/lib64/libfesom.so` (built
2026-05-25 16:46). This is **bit-identical to the binary+lib pair that ran
the F-EVP leg** (`work_linfs_2yr_b`) — zero binary skew in the M5
diff-of-diffs — and includes the May 25 sbc cold-start rotation fix
(2682a9fb) plus the committed KPP-era dump instrumentation (9271ae92,
env-gated, inactive here). The 2026-06-11 uncommitted `EVPdynamics_m` dump
patch (M0 instrumentation) is NOT in this frozen pair; it only affects
`build/lib/libfesom.so`, used by the separate 16r dump run
(`work_mevp_dump`, job 25514438).

NOTE: fesom.x is a thin launcher; the model code lives in `libfesom.so`
resolved via RPATH. Never point a reference run at `build/bin/fesom.x` —
it would silently track later rebuilds.

## Run completion (M0 close-out, 2026-06-11)

Job 25514027 **COMPLETED** clean: runtime 1573 s for 35040 steps, no error
strings, full 1958+1959 output set (62 files incl. uice/vice + a_ice/m_ice).
Option-liveness check vs the F-EVP leg (1959): m_ice max|diff| = 2.51 m,
mean|diff| = 0.006 m (visibly nonzero, near-identical global means 0.4274 vs
0.4280 m); Sept p99 |uice| 0.346 (mEVP) vs 0.161 m/s (EVP) — the rheology
flip is live and physically plausible. These are also the F-legs of the M5
diff-of-diffs.

The 16r dump run (job 25514438, instrumented build) COMPLETED with 448 dump
files (2 steps × 14 points × 16 ranks) at
`/work/ab0995/a270088/port/mevp/fdump_16r/dump/`; validated finite + plausible
(pressure_fac max 32.35 ≈ det2·pstar·msum·exp(·); u_aux 6 mm/s @substep 1 →
0.33 m/s @120; UF ≡ it120; initial-state NH/SH ice values exact).

## ⚠️ Stale `work_core_evp1` — do NOT use

`/home/a/a270088/fesom27/fesom2/work_core_evp1/` is a January experiment:
its log shows `EVP scheme option= 0` (the flip was never active there),
namelist.ice identical to work_core, and OLD forcing settings
(`nm_nc_iyear=1948`, `nm_nc_tmid=1`; current convention is 1900/0).
It was not reused. Also note `fesom27/fesom2/work_core/namelist.config`
itself is January-stale (which_ALE='zstar', run_length=20y) — only its
physics namelists are the live lineage; the run-control config comes from
`work_linfs_2yr_b` as described above.
