# Next-session prompt — start the zstar port (Phase Z0)

Copy the block below as the next session's opening prompt.
(The previous prompt's task — the step-1 forcing diff — was RESOLVED 2026-05-25: Fortran
cold-start bug in `nc_sbc_ini`, C correct; see memory `project_forcing_step1_diff`.)

---

> Continue the FESOM2→C port. Baseline: **tag `v1.0` on main@1a674fb** (linfs+KPP, complete,
> climate-validated, 5yr stable). **This session: start implementing the zstar vertical
> coordinate — Phase Z0 of `docs/plans/20260610-zstar-vertical-coordinate.md`** (plan-reviewed
> 2026-06-10, all findings folded in; companion TKE plan `20260610-tke-vertical-mixing.md`
> comes AFTER zstar per the approved order).
>
> 1. **Read the plan first**, then the auto-memory `project_zstar_tke_plans` (approved
>    decisions + the plan-review catches — esp. the `real_salt_flux` blocker: the C currently
>    DISCARDS `rsf` at `fesom_ice_thermo.c:616`; under zstar it replaces virtual_salt).
> 2. **Branch from main@v1.0** (the repo may still have `jax-mesh-export` checked out with
>    unrelated JAX-dump WIP — don't disturb it). Commit the two plan files (they're untracked).
> 3. **Z0 first moves:** add the `FESOM_ALE=linfs|zstar` switch + mesh arrays (Z_3d_n, dhe,
>    bottom_*_thickness, real_salt_flux) at linfs-identical values; **launch the two Fortran
>    reference runs immediately** (SLURM queue time): `work_zstar_dump` (16r, dt=1800, 3 steps,
>    ALE dump gate) and `work_zstar_kpp` (2yr, dt=1800; copy of `work_linfs_2yr` with ONLY
>    `which_ALE='zstar'` flipped). Unique OUT_DIRs. Archive namelists + PROVENANCE.
> 4. **Gate:** linfs (default) byte-identical to v1.0 before any zstar branch goes live.
>
> GUIDING PRINCIPLE: strict faithful port — namelist values over module defaults; port only
> what the zstar reference exercises (NO zlevel/local-zstar fallback); every phase ends with
> its Validate gate (dump-diff vs Fortran + linfs byte-gate).

---

## Quick reference (paths + commands)

- **Plans:** `docs/plans/20260610-zstar-vertical-coordinate.md` (Z0–Z10),
  `docs/plans/20260610-tke-vertical-mixing.md` (T0–T6, after zstar).
- **Approved run matrix:** zstar+KPP 2yr pair, linfs+TKE 2yr pair, zstar+TKE 2yr pair +
  5yr C stability (memory `project_zstar_tke_plans`).
- **Fortran tree:** `/home/a/a270088/port2/fesom2` (Intel, `CVMIX:BOOL=ON` already; carries
  uncommitted dump shims — line numbers drift, anchor by routine). Reference work dirs:
  `work_linfs_2yr` (the template), new `work_zstar_*` / `work_*_tke` to be created.
- Build C: `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`.
  Build Fortran: `source /sw/etc/profile.levante && source env/levante.dkrz.de/shell &&
  make -C build fesom.x -j8 && make -C build install` (install mandatory — `bin/fesom.x`
  loads `lib64/libfesom.so`). Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`
  (PYTHONPATH=/home/a/a270088/PYTHON).
- **Byte-gate tooling:** `scripts/kpp_byteident_check.py` (reuse); dump-diff scripts to be
  cloned from `scripts/kpp_dump_diff.py` → `scripts/ale_dump_diff.py`.
