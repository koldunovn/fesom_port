# Next-session prompt — mEVP port: execute the plan, starting at M0

(The previous prompt's open items were resolved 2026-06-11: zstar+TKE were MERGED to
main — fast-forward 1a674fb→2714071, so now `main` == `zstar` == 2714071; defaults stay
KPP+linfs; IDEMIX stays unported. Same day the mEVP work was brainstormed and planned.)

**THE TASK:** implement the mEVP sea-ice rheology port (whichEVP=1) by executing
`docs/plans/20260611-mevp-sea-ice-rheology.md` — committed on branch **`mevp`** in THIS
worktree (`/home/a/a270088/port2/fesom2_port_zstar/`, checked out on `mevp`). The plan
was written from the 2026-06-11 brainstorm and revised after TWO plan-review agents —
**read it in full before coding**, especially:

- the **14-item fidelity-traps checklist** (headliners: `rdt = FULL ice_dt` AND drag
  carries rdt; the 0.5 lives in the sigma11/22 update NOT in pressure_fac; NO theta_io;
  non-ice nodes are SKIPPED not zeroed — do NOT copy the std-EVP `else u=0`;
  `uice_old` never touched);
- the **`bc_index_nod2D` rebuild** (plan-review catch): the C array EXISTS
  (fesom_mesh.h:141, built fesom_ice.c:199-209) but is edge_tri-built =
  multi-rank-unsafe, and mEVP is its FIRST reader → M1 rebuilds it with
  `myList_edge2D > edge2D_in`. If M3/M4 show seam-concentrated velocity diffs,
  suspect this first.

All design decisions are user-approved (2026-06-11) — do not re-ask: knob
`FESOM_WHICH_EVP=1` (numeric; unset/0 default, 2 aborts); run matrix = single-knob
Fortran reference + 5yr C stability (NO combined zstar/TKE/mEVP leg); new
`fesom_ice_maevp.{c,h}`; default byte-identity is a gate in every phase.

## M0 first (queue-time dominates — submit before touching C)

1. **Byte-gate baseline BEFORE any src change** (HEAD = 2714071 + docs-only commits →
   binary ≡ 2714071): build, run 16r dist_16 / 200 steps / dt=500 / knob unset, archive
   snapshots + run.log → `/work/ab0995/a270088/port/mevp/baseline_2714071/`.
2. **Fortran 2yr reference**: copy `/home/a/a270088/fesom27/fesom2/work_core` →
   `work_mevp` (⚠️ NOT the stale `work_core_evp1` — January experiment, whichEVP=0 in
   its log, old forcing); flip `whichEVP=1` (ONLY knob; alpha/beta=250 + steps=120
   already in namelist.ice); unique ResultPath; submit 2yr dt=1800.
   Gate: log echoes `EVP scheme option= 1` (fesom_module.F90:323).
3. **Fortran dump instrumentation**: env-gated calls in `EVPdynamics_m`
   (ice_maEVP.F90:429-882) reusing `module fesom_evp_dump` (ice_EVP.F90:20-84);
   uncommitted Intel-tree patch (KPP/TKE precedent); 16r 2-step dump run
   (`FESOM_EVP_DUMP_DIR`). Dump points per the plan's Validation Strategy.
4. Archive namelists → `docs/mevp_reference_namelists/` + PROVENANCE.md.

Then M1 (C skeleton) → M2 (the one-function port) → M3 (substep dump-diff) →
M4 (cross-rank) → M5 (2yr climate + ice metrics + VECTOR check + diff-of-diffs) →
M6 (5yr + tag + handoff). Phase gates and bars are in the plan.

## Quick reference

- **Plan:** `docs/plans/20260611-mevp-sea-ice-rheology.md` (this worktree, branch mevp).
- **Memory:** `project_mevp_port_state` (decisions, review catches, repo-state proofs).
- **Fortran source:** `/home/a/a270088/port2/fesom2/src/ice_maEVP.F90` (EVPdynamics_m
  429-882); MOD_ICE.F90:198/743-748/889-895; dispatch ice_setup_step.F90:209-226.
- **C side:** `src/fesom_ice.c` (init :68/:199-209, dispatch :327, atoi precedent
  :392-393), `src/fesom_ice_types.h`, `src/fesom_ice_evp.c` (std-EVP template + dump
  infra — beware the traps), `src/fesom_mesh.h:35-40,141`.
- **Build C:** `source /sw/etc/profile.levante && source env.sh && make -C build
  fesom_port` (fresh configure: `bash -l configure.sh`). Fortran: `make -C build
  fesom.x -j8 && make -C build install`.
  Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`.
- **Jobs:** clone `jobs/job_tke_t0_byteident` / `job_zstar_z7_crossrank` /
  `validation_2yr.sh` / `production_864_5yr.sh` → `job_mevp_*`; UNIQUE OUT_DIR per job.
- **Outputs root:** `/work/ab0995/a270088/port/mevp/`.
- **Switches today:** `FESOM_MIX_SCHEME=KPP(default)|PP|TKE`,
  `FESOM_ALE=linfs(default)|zstar`, `FESOM_EVP_DUMP_DIR`; this port adds
  `FESOM_WHICH_EVP=0(default)|1`.
