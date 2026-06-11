# Next-session prompt — TKE port DONE; run matrix CLOSED; open user decisions

(The previous prompt's task — the CVMix classical-TKE port, T0–T6 of
`docs/plans/completed/20260610-tke-vertical-mixing.md` — was COMPLETED
2026-06-11 in one session, tag **`tke-validated-2026-06-11`**; AND the
approved run matrix was fully closed the same session: zstar+KPP ✓ (Z9),
linfs+TKE ✓ (T5), zstar+TKE 2yr ✓, 5yr zstar+TKE stability ✓ (peak uv 2.93,
margin 2.07). See `docs/validation_tke_2yr/RESULTS.md` and auto-memory
`project_tke_port_state`.)

There is NO queued implementation work. What remains are USER DECISIONS:

1. **Merge** branch `zstar` (carrying both the zstar AND TKE ports, each
   opt-in and validated) into `main`? The main checkout `fesom2_port/` still
   holds jax-mesh-export WIP — merge order/timing is yours.
2. **Defaults:** KPP and linfs remain the defaults (`FESOM_MIX_SCHEME=TKE`,
   `FESOM_ALE=zstar` opt in). Flip either?
3. **IDEMIX** (`cvmix_TKE+cvmix_IDEMIX`): port later? The gate-only seams
   (iw-input plumbing) are staged in fesom_cvmix_tke.c/fesom_tke.c.
4. Manual science checks suggested by the plan (MLD maps vs obs in
   Labrador/Weddell/GIN with TKE vs KPP; TKE budget profiles vs ICON
   literature with FESOM_TKE_DIAG=1 + FESOM_IO_CONFIG).

---

## Quick reference

- **TKE (DONE 2026-06-11):** commits 45afc01/fcd152e/f5fbe58+, tag
  `tke-validated-2026-06-11`; results `docs/validation_tke_2yr/RESULTS.md`;
  references under `/work/ab0995/a270088/port/tke/` (fortran_linfs_tke,
  fortran_zstar_tke, fdump/dump, replay, c_tke_2yr, c_zstar_tke_2yr/5yr);
  switches: `FESOM_MIX_SCHEME=KPP(default)|PP|TKE`, `FESOM_TKE_DIAG=1`,
  `FESOM_TKE_DUMP_DIR`/`FESOM_TKE_REPLAY_DIR`/`FESOM_TKE_HALO_PROBE`.
- **zstar (DONE 2026-06-10):** tag `zstar-validated-2026-06-10`;
  `FESOM_ALE=linfs(default)|zstar`.
- Build C: `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`
  (worktree). Fortran: `make -C build fesom.x -j8 && make -C build install`.
  Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`.
- Port lessons added this session: `feedback_r8_literals_double` (-r8 makes
  ALL default-real literals double — replay caught a 6.6-as-single bug);
  job self-cleaning (never manual `rm` on submitted jobs' output dirs).
