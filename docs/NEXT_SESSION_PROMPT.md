# Next-session prompt — start the TKE port (Phase T0)

Copy the block below as the next session's opening prompt.
(The previous prompt's task — zstar Z7/Z8/Z9 verdicts → Z10 — was COMPLETED 2026-06-11:
**the zstar port is DONE, Z0–Z10**, tag `zstar-validated-2026-06-10`, plan moved to
`docs/plans/completed/`. User kept **linfs as the default**; `FESOM_ALE=zstar` opts in.)

---

> Continue the FESOM2→C port in the worktree **`/home/a/a270088/port2/fesom2_port_zstar`**
> (branch `zstar`; the main `fesom2_port/` checkout still holds jax-mesh-export WIP — don't
> disturb it). **This session: start the CVMix classical-TKE port — Phase T0 of
> `docs/plans/20260610-tke-vertical-mixing.md`** (plan-reviewed 2026-06-10; the second of
> the two approved plans, zstar-first order ✓ done).
>
> 1. **Read the plan first**, then auto-memories `project_zstar_tke_plans` (approved
>    decisions: two C files mirroring Fortran — `fesom_cvmix_tke.c` = cvmix_tke.F90 column
>    core, `fesom_tke.c` = gen_modules_cvmix_tke.F90 driver; `FESOM_MIX_SCHEME=TKE` third
>    option, KPP stays default; diagnostics computed per-column but stored only under
>    `FESOM_TKE_DIAG=1`; `tke_cd=3.75` from namelist.cvmix NOT the 1.0 module default;
>    `handle_old_vals`/old_KappaM/H are DEAD; driver loops OWNED-ONLY; `iw_diss` read
>    unconditionally → pass a real zero column) and `project_zstar_port_state`.
> 2. **T0 first moves:** `struct fesom_tke` + `FESOM_MIX_SCHEME=TKE` abort-stub dispatch +
>    `FESOM_TKE_DIAG` parse; **launch the Fortran reference runs immediately** (queue
>    time): `work_tke_dump` (16r, dt=1800, 3 steps, dump gate) and `work_linfs_tke`
>    (2yr; copy work_linfs_2yr, flip `mix_scheme='cvmix_TKE'`, copy namelist.cvmix in;
>    same Intel binary — CVMIX=ON already). Unique OUT_DIRs, purge-safe /work, archive
>    namelists + PROVENANCE (the zstar Z0 pattern; consider a frozen binary copy like
>    `frozen_zstar_refs/`).
> 3. **Provenance hardening (plan finding 6):** grep BOTH reference logs for the
>    `init_cvmix_tke` echo (`tke_cd = 3.75`, `tke_only = T`) — proves the namelist was
>    consumed AND `__cvmix` is compiled (else a plausible-looking broken reference).
> 4. **Gate:** TKE off → byte-identical to the current HEAD (16r/200/dt=500;
>    `jobs/job_zstar_z0_byteident` is the recipe — point run A at a HEAD-built binary).
>
> GUIDING PRINCIPLE unchanged: strict faithful port; namelist values over module defaults;
> port what the reference exercises (Langmuir/IDEMIX/Dirichlet gate-only); every phase ends
> with its Validate gate.

---

## Quick reference

- **zstar (DONE):** tag `zstar-validated-2026-06-10`; results
  `docs/validation_zstar_2yr/RESULTS.md`; references under
  `/work/ab0995/a270088/port/zstar/` (fortran_zstar_2yr, fortran_linfs_2yr_b, fdump/dump);
  switch `FESOM_ALE=linfs|zstar` (default linfs). One leftover formality: the
  linfs-baseline cross-rank job 25499700 (`z7_xrank_linfs/`) — paste its day-5 spread
  next to the zstar one in RESULTS.md if it looks as expected (same IC-noise scale).
- **Run matrix remainder (approved):** after TKE validates per-feature →
  Fortran `work_zstar_tke` 2yr (both knobs) + C zstar+TKE 2yr compare + 5yr C
  zstar+TKE stability (re-measure the autumn uv margin; linfs was 0.35).
- Build C: `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`
  (worktree). Fortran: `make -C build fesom.x -j8 && make -C build install` (install
  mandatory). Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`.
- Byte/dump tooling to clone: `jobs/job_zstar_z0_byteident`, `scripts/ale_dump_diff.py`
  → `scripts/tke_dump_diff.py`; Fortran dump-gate pattern = `ale_dump_mod` (uncommitted,
  top of `fesom2/src/oce_ale.F90`) / the KPP gates in `oce_ale_mixing_kpp.F90`.
