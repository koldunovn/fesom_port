# Next-session prompt — zstar Z7/Z8/Z9 verdicts → Z10

Copy the block below as the next session's opening prompt.
(Previous prompt's task — "start the zstar port" — was executed 2026-06-10 FAR past Z0:
**Z0–Z7 are implemented, dump-validated, and committed** on branch `zstar`;
three validation runs were left in flight.)

---

> Continue the zstar port in the worktree **`/home/a/a270088/port2/fesom2_port_zstar`**
> (branch `zstar`, HEAD 89611d5; the main `fesom2_port/` checkout still holds the
> jax-mesh-export WIP — don't disturb it). Read auto-memory `project_zstar_port_state`
> and the plan `docs/plans/20260610-zstar-vertical-coordinate.md` (Z0–Z7 checked off with
> evidence; Z8/Z9 partially noted).
>
> **This session: evaluate the three in-flight validation jobs and close Z7→Z9:**
> 1. **25495447 `z7_xrank`** (`/work/ab0995/a270088/port/zstar/z7_xrank/slurm.*.out`):
>    zstar 5-day at 1/8/16 ranks + pairwise final-snapshot compare. Bar: agreement at the
>    PHC rank-noise scale (GM/sea-ice precedent) — NOT bit-identity. A large/growing
>    cross-rank gap = halo bug in the new zstar write loops (Wvel/hnode_new/helem/
>    ssh_rhs_old/hbar/dhe) → exchange-and-compare probe per feedback_write_loops_halo.
> 2. **25495448 `z8_month`** (`.../z8_month/`): 30-day dt=1800 16r zstar. Bar: clean exit,
>    bounded max|uv| (linfs margin lore: autumn peak 4.65), physical SSH/SST/SSS. Watch
>    Aleutian el 194724 + the 1-ring hbar-spread checkerboard.
> 3. **25495449 `z9_zstar2yr`** (`.../c_zstar_2yr/`): C zstar+KPP 2yr 864r. Run
>    `scripts/zstar_climate_compare.py` (paths pre-wired: F-zstar reference
>    `fortran_zstar_2yr`, contrast `fortran_linfs_2yr_b` — both COMPLETE, exit 0).
>    Bar: KPP-class (SST/SSS RMS ~0.005/0.002, biases ~0, non-drifting); the LIN contrast
>    column must be CLEARLY larger than ZS (the comparison resolves the coordinate).
> 4. Then Z9 checkbox + `docs/validation_zstar_2yr/` numbers, and **Z10**: handoff,
>    plan → completed/, memory update, tag suggestion `zstar-validated-<date>`,
>    user decision point on the default-coordinate flip. The TKE plan
>    (`docs/plans/20260610-tke-vertical-mixing.md`) comes next.
>
> Known/expected: steps 2-3 cross-code dump diffs are threshold-amplified input noise
> (kbl flips 7.9%, fully analyzed in the plan Z8 note + the Z2-Z7 commit message) — do
> NOT re-litigate; the climate run is the binding test. linfs byte-gates PASSED at
> Z0/Z1/Z2-Z6/Z2-Z7 (worst |Δ|=0 every time; run.log diffs are the 3 out_dir path
> strings only) — `jobs/job_zstar_z0_byteident` re-runs it in ~15 min if needed.

---

## Quick reference

- **Worktrees:** `fesom2_port_zstar` (branch zstar — THE working tree),
  `fesom2_port_v10` (detached v1.0 — byte-gate reference binary; keep).
- **Fortran refs (purge-safe):** `/work/ab0995/a270088/port/zstar/{fortran_zstar_2yr,
  fortran_linfs_2yr_b,fdump/dump,fdump_k2/kdump}`; frozen binary
  `fesom2/frozen_zstar_refs/` (sha 700e2d19…); namelists archived in
  `docs/zstar_reference_namelists/` + PROVENANCE.md. Fortran ALE dump gate =
  `ale_dump_mod` UNCOMMITTED in `fesom2/src/oce_ale.F90` (KPP-gate convention).
- **Dump tooling:** `FESOM_ALE_DUMP_DIR` + `FESOM_ALE_DUMP_STEPS` (both codes,
  identical format), `scripts/ale_dump_diff.py`; per-run C dumps under
  `/work/ab0995/a270088/port/zstar/{z1_smoke,z2_cdump,z7_kdump}`.
- Build C: `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`
  (from the worktree). Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`
  (PYTHONPATH=/home/a/a270088/PYTHON).
- Switch: `FESOM_ALE=linfs|zstar` (default linfs = byte-identical to v1.0).
