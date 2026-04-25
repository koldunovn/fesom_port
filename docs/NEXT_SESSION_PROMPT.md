# Next session: Phase D5 (oce_fluxes_mom) → Phase E (FCT advection)

## TL;DR

Phase D EVP is **done and validated** at 1, 8, 16 ranks (200 steps clean,
T.min=−2.06°C, bit-identical step-1 stats). The previous "16-rank crash at
step 2" was a JRA forcing OOB read fixed in commit `55f6b6c` (now on `main`
via `918b0d4`). Two phases remain to finish sea-ice port:

- **Phase D5**: port `oce_fluxes_mom` so ice mediates atmospheric stress on
  the ocean (today the ocean still sees raw bulk wind stress under ice cover).
- **Phase E**: port FCT advection of `a_ice`, `m_ice`, `m_snow`.

Both are scoped in `docs/plans/20260425-sea-ice-port.md` (Tasks D5 → D6 → E1…E8).

## Read these memory files first (in order)

```
~/.claude/projects/-home-a-a270088-port2/memory/MEMORY.md
~/.claude/projects/-home-a-a270088-port2/memory/project_sea_ice_port_state.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_array_size_vs_reader_loop.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_write_loops_halo.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_phc_rank_dependent.md
~/.claude/projects/-home-a-a270088-port2/memory/reference_evp_dump_diagnostic.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_unique_outdir_per_job.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_levante_build.md
```

`feedback_array_size_vs_reader_loop.md` and `feedback_write_loops_halo.md`
are critical for Phase E — FCT is exactly where this bug class hides.

## Repo state (start of next session)

```
main          918b0d4  Merge debug-evp: Phase D EVP + JRA halo fix (HEAD)
              55f6b6c  Fix multi-rank sea-ice EVP crash: JRA forcing halo coverage
              edc9e65  Phases A–C
              3da152b  Phase 4 wrap-up (ocean only, pre-sea-ice)
debug-evp     == main  (merged; safe to delete or reuse for next phase)
```

Always:
```
bash -l configure.sh --clean
```

Always use unique OUT_DIR per SLURM job (`feedback_unique_outdir_per_job.md`).

## What's already in tree

- Sea-ice Phases A–D fully wired in `fesom_ice.{c,h}`,
  `fesom_ice_thermo.{c,h}`, `fesom_ice_coupling.{c,h}`, `fesom_ice_evp.{c,h}`.
- EVP debug instrumentation env-gated by `FESOM_EVP_DUMP_DIR` plus an
  always-on `if (ice_strength[el] > 1e7) fprintf(stderr, ...)` in EVP Step 3.
  See `reference_evp_dump_diagnostic.md`. Use this for any multi-rank
  divergence hunt — it earned its keep finding the JRA bug.
- Diff scripts under `scripts/`: `evp_dump_diff.py`, `evp_dump_qcheck.py`,
  `phc_dump_diff.py`. Job templates `job_core2_{1,8,16}r_evp_dump`,
  `job_core2_8r_evp`.

## Phase D5 — port `oce_fluxes_mom`

**Source**: `~/port2/fesom2/src/ice_oce_coupling.F90:53` (subroutine `oce_fluxes_mom`).

**Goal**: blend bulk wind stress with ice-ocean drag from `cd_oce_ice` and
ice velocity, write into `dyn->stress_surf` (per-element). Currently
`dyn->stress_surf` is the unmodified bulk stress from
`fesom_bulk_compute` — under ice cover this is wrong.

**Design**:
- Add `fesom_ice_oce_fluxes_mom` to `fesom_ice_coupling.{c,h}`.
- Call from `fesom_step.c` BEFORE `compute_vel_rhs` (which consumes
  `stress_surf`). Place AFTER `fesom_ice_step` so `ice->uice/vice` are
  populated by the EVP just run.
- Per Fortran: element-based loop (`do el=1, myDim_elem2D`), reads
  `u_ice/v_ice` at the 3 vertices (halo OK from end-of-EVP exchange), reads
  `a_ice` per element (from vertex average), writes 2 components of
  `stress_surf[el]`.
- **Loop bound for stress_surf write**: check Fortran. If
  `do el = 1, myDim_elem2D + eDim_elem2D` then C must too. Apply the audit
  rule from `feedback_write_loops_halo.md` — this is the cheapest safety net.

**Validation** (Task D6 in plan):
1. Clean build.
2. Run `job_core2_16_fix` 200 steps; under ice cover the polar `|stress_surf|`
   should be visibly smaller than pre-D5 (because ice is now mediating).
3. Multi-rank cross-check: `|Δuice|_∞ < 1e-6` between 8/16/32 ranks at step 200.
4. Watch the always-on HUGE-ice_strength stderr — should never fire.

**Watchouts for D5**:
- The EVP halo exchange of u_ice/v_ice already happens at end of subcycle.
  After `fesom_ice_step`, halo of u_ice IS valid for stress_mom to read.
- Scattered prior attempt at `/tmp/sea_ice_modified.patch` may be stale —
  re-derive from Fortran.
- `oce_fluxes_mom` does NOT touch `stress_node_surf` (that's still bulk-direct
  for the ocean RHS at nodes). Read Fortran carefully.

## Phase E — FCT advection (highest-risk)

**Source**: `~/port2/fesom2/src/ice_fct.F90` (1624 LoC). Plan Tasks E1–E8.

**Why high-risk**: ocean FCT had two halo-write bugs (`init_tracers_AB` and
`fct_ttf_max/min`) that took days to find. The same bug class will recur.
**Apply the proactive audit**: grep every C `for.*<.*myDim_(nod|elem)2D;`
write-loop in `fesom_ice_fct.c` against the matching Fortran loop bound,
BEFORE running anything. The plan calls this Task E7 — do it as you write,
not after.

**Tasks (from plan)**:
- E1: `ice_mass_matrix_fill` (one-time setup)
- E2: `ice_TG_rhs` + `ice_TG_rhs_div`
- E3: `ice_solve_low_order` + `ice_solve_high_order`
- E4: `ice_fem_fct` (Zalesak limiter)
- E5: `ice_fct_solve` (driver)
- E6: wire into `fesom_ice_step`
- E7: proactive halo-write audit
- E8: validation

After E6: between EVP and thermo in `fesom_ice_step`, call
`fesom_ice_tg_rhs → fesom_ice_fct_solve → fesom_ice_cut_off` (mirroring
`ice_setup_step.F90:258-295`).

**Validation** (Task E8): full month at multi-rank; ice-cover area in
correct order of magnitude (~10% of ocean area); a_ice physically bounded
(0…1); compare T.min trajectory to Fortran reference if available.

## Workflow rules earned the hard way

1. **Always clean build** (`bash -l configure.sh --clean`).
2. **Unique OUT_DIR per SLURM job** (and unique `FESOM_EVP_DUMP_DIR` if
   using diagnostics) — `feedback_unique_outdir_per_job.md`.
3. **2× repeat tests** before declaring "validated". One PASS isn't enough.
4. **Commit at every working state** with a clear message.
5. **Follow Fortran literally** — no "simplifying" loop bounds, no
   "saving memory" by sizing arrays smaller. The forced halo coverage
   isn't accidental.
6. **PHC IC differs across rank counts** — accept it. Don't burn cycles
   diffing PHC at 1-rank vs N-rank. See `feedback_phc_rank_dependent.md`.
7. **Use the dump diagnostic** — set `FESOM_EVP_DUMP_DIR` and run
   `scripts/evp_dump_diff.py` whenever multi-rank diverges from 1-rank.
   The HUGE diagnostic is always on; check stderr first.

## Quick commands

```bash
# Dev cycle
cd ~/port2/fesom2_port
bash -l configure.sh --clean

# 16-rank smoke test (200 steps)
sbatch job_core2_16_fix

# 16-rank with full dumps (NSTEPS=2 for diagnostic)
sbatch job_core2_16r_evp_dump

# Diff against 1-rank reference
python3 scripts/evp_dump_diff.py <1r_dir> <16r_dir>
```

## Bug class to watch for in Phase E

The exact pattern that bit us in Phase D was: array sized myDim,
read at halo by a compute-over-halo kernel → garbage halo values →
downstream physics explodes. In Phase E, the analogues are arrays like
`fct_tmax`, `fct_tmin`, `fct_plus`, `fct_minus`, `fct_fluxes`,
`values_rhs`, `values_div_rhs`, `dvalues`, `valuesl`. Each:
1. allocated to size N (check what N is — `myDim` or `myDim+eDim`?)
2. written by which function with which loop bound?
3. read by which function with which loop bound?

If READER bound > WRITER bound (or > allocation size): bug. The Fortran
sizes them all `myDim+eDim` — match it in C.
