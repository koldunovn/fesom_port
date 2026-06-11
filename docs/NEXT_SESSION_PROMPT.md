# Next-session prompt — run-matrix close-out (zstar+TKE combined legs)

Copy the block below as the next session's opening prompt.
(The previous prompt's task — the CVMix classical-TKE port, T0–T6 of
`docs/plans/20260610-tke-vertical-mixing.md` — was COMPLETED 2026-06-11:
**the TKE port is DONE**; see `docs/validation_tke_2yr/RESULTS.md` and
auto-memory `project_tke_port_state`. KPP stays the default;
`FESOM_MIX_SCHEME=TKE` opts in.)

---

> Continue the FESOM2→C port in the worktree **`/home/a/a270088/port2/fesom2_port_zstar`**
> (branch `zstar`). **This session: close the approved run matrix** (2026-06-10
> decisions, shared by the zstar and TKE plans — both features are now
> individually validated; the COMBINED zstar+TKE legs remain):
>
> 1. Check the Fortran combined reference `work_zstar_tke` (job 25503765,
>    two knobs: which_ALE='zstar' + mix_scheme='cvmix_TKE', frozen binary,
>    output `/work/ab0995/a270088/port/tke/fortran_zstar_tke/`) — if it was
>    still queued/running at session end, verify exit 0 + the init_cvmix_tke
>    echo (cd=3.75) + 730 days.
> 2. Submit `jobs/job_tke_m1_zstar_tke_2yr` (C zstar+TKE 2yr 864r) and
>    compare vs the Fortran combined reference (adapt
>    `scripts/tke_climate_compare.py` paths; bar: the same KPP-class RMS).
> 3. Submit `jobs/job_tke_m2_zstar_tke_5yr` (5yr stability; re-measure the
>    autumn uv-transient margin — linfs lore: peak 4.65, margin 0.35).
> 4. Record results in `docs/validation_tke_2yr/RESULTS.md` (matrix section)
>    + auto-memory; the matrix is then CLOSED:
>    zstar+KPP ✓ (Z9), linfs+TKE ✓ (T5), zstar+TKE ✓ (this session).
>
> OPEN DECISIONS for the user (don't decide unilaterally): merging branch
> `zstar` (zstar + TKE, both validated, both opt-in) into main; whether TKE
> or zstar ever become defaults; IDEMIX port (gate-only seams staged).

---

## Quick reference

- **TKE (DONE 2026-06-11):** commit `45afc01`+, tag `tke-validated-2026-06-11`;
  results `docs/validation_tke_2yr/RESULTS.md`; references under
  `/work/ab0995/a270088/port/tke/` (fortran_linfs_tke, fdump/dump, replay);
  switch `FESOM_MIX_SCHEME=KPP(default)|PP|TKE`; `FESOM_TKE_DIAG=1` stores
  the 13 budget slabs; `FESOM_TKE_DUMP_DIR`/`FESOM_TKE_REPLAY_DIR`/
  `FESOM_TKE_HALO_PROBE` for diagnostics.
- **zstar (DONE 2026-06-10):** tag `zstar-validated-2026-06-10`;
  `FESOM_ALE=linfs(default)|zstar`.
- Build C: `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`
  (worktree). Fortran: `make -C build fesom.x -j8 && make -C build install`.
  Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`.
- Lessons learned this port: `-r8` makes ALL default-real Fortran literals
  double (`feedback_r8_literals_double`); job outputs self-clean inside the
  job scripts — never manual `rm` against submitted jobs' dirs.
