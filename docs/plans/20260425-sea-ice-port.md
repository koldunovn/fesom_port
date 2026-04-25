# Sea-Ice Port for FESOM2 C Port

## Overview

Bring the ocean-only C port (currently runs CORE2 stably 1–256 ranks) to feature
parity with the standard FESOM2 sea-ice path:

- `whichEVP = 0` (standard EVP only — `ice_EVP.F90`)
- `__icepack` not defined (non-coupled thermodynamics — `ice_thermo_oce.F90`)
- `use_meltponds = .false.`, no cavity, no icebergs, no `__oifs`/`__yac`
- `num_itracers = 3` (a_ice, m_ice, m_snow)
- Reference namelist: `/home/a/a270088/fesom27/fesom2/work_core/namelist.ice`

**Why now:** the ocean-only port has `T.min` drifting to ≈ −27 °C over a model
month because nothing insulates the ocean from cold air at high latitudes
(see `project_mpi_port_state.md`). Sea ice is the missing physics; resolving
that drift is the single best end-to-end validation that the port lands.

**Side-effect that must be fixed in lock-step:** `fesom_sss_runoff_step` currently
subtracts runoff from `water_flux` because we have no `therm_ice` to fold it
into `fresh_wa_flux`. The handover happens in **two paired tasks**: B4b folds
runoff into `flx_fw` inside `therm_ice`, and C3b removes the subtraction in
`fesom_sss_runoff` once C2 wires `flx_fw → water_flux` via `oce_fluxes`. Both
must land for freshwater to balance — neither alone is correct
(see `feedback_runoff_fold_when_ice.md`).

## Context (from discovery)

### Files involved

**Reference Fortran** (all under `/home/a/a270088/port2/fesom2/src/`):

| File | LoC | Role |
|---|---|---|
| `MOD_ICE.F90` | 950 | `T_ICE` type + `ice_init` (line 572) — state shape + namelist |
| `ice_setup_step.F90` | 563 | `ice_timestep` (96), `ice_initial_state` (358), `ice_setup` (51) |
| `ice_EVP.F90` | 833 | `stress_tensor` (46), `stress2rhs` (177), `EVPdynamics` (330) |
| `ice_fct.F90` | 1624 | `ice_TG_rhs` (91), `ice_fct_solve` (210), `ice_solve_low_order` (243), `ice_solve_high_order` (347), `ice_fem_fct` (498), `ice_mass_matrix_fill` (1145), `ice_TG_rhs_div` (1253), `ice_update_for_div` (1452) |
| `ice_thermo_oce.F90` | 1015 | `cut_off` (70), `thermodynamics` (148), `therm_ice` (367), `budget` (657), `obudget` (777), `flooding` (872), `TFrez` (899) |
| `ice_oce_coupling.F90` | 841 | `oce_fluxes_mom` (53), `ocean2ice` (154), `oce_fluxes` (249); non-ICEPACK branch at 393–398 |

**Out of scope (skipped, with loud aborts where needed):**
- `ice_maEVP.F90` (mEVP/aEVP — `whichEVP ∈ {1,2}`). Stub the EVP dispatcher to abort if `whichEVP != 0`.
- `ice_thermo_cpl.F90` (coupled-atmosphere thermodynamics).
- `ice_meltponds.F90` (`thermo%use_meltponds = .false.`).
- Cavity, icebergs, `__oifs`/`__yac`/`__oasis` paths.
- `ice_update_for_div` (`ice_fct.F90:1452`) — commented out at the call site (`ice_setup_step.F90:264`); not used in the standard path.
- ~~`ice_initial_state` (`ice_setup_step.F90:358`) — covered implicitly by calloc-zeroed fields + snapshot read in Task A5~~. **CORRECTED 2026-04-25**: actually needed. Cold-start (Fortran lines 500-521) sets m_ice=1.0 (NH) / 2.0 (SH), m_snow=0.1/0.5, a_ice=0.9 wherever surface T<0°C — calloc-zero is wrong. Ported as `fesom_ice_initial_state` in `fesom_ice.{c,h}`, called from `fesom_main.c` after PHC IC. The "from-file" branch (`ini_ice_from_file=true`) remains out of scope.
- CMIP6 dynamic-growth diagnostics `dyngr/dyngrsn/dyngra` (`MOD_ICE.F90:901-904`, `ice_setup_step.F90:203-205,307-309`) — output-only fields with no physics impact. Skip.

**C port target** (all under `/home/a/a270088/port2/fesom2_port/src/`):

| File | Mirror of | Notes |
|---|---|---|
| `fesom_ice_types.h` | parts of `MOD_ICE.F90` | `struct fesom_ice` mirrors `T_ICE`; nested `work` and `thermo` mirror `T_ICE_WORK`/`T_ICE_THERMO` |
| `fesom_ice.{c,h}` | `ice_setup_step.F90` + `MOD_ICE.F90:ice_init` | `fesom_ice_init`, `fesom_ice_initial_state`, `fesom_ice_setup`, `fesom_ice_step`; env knobs |
| `fesom_ice_evp.{c,h}` | `ice_EVP.F90` | `fesom_ice_stress_tensor`, `fesom_ice_stress2rhs`, `fesom_ice_evp_dynamics` |
| `fesom_ice_fct.{c,h}` | `ice_fct.F90` | `fesom_ice_tg_rhs`, `fesom_ice_fct_solve`, `fesom_ice_solve_low_order`, `fesom_ice_solve_high_order`, `fesom_ice_fem_fct`, `fesom_ice_mass_matrix_fill` |
| `fesom_ice_thermo.{c,h}` | `ice_thermo_oce.F90` | `fesom_ice_thermodynamics`, `fesom_therm_ice`, `fesom_ice_budget`, `fesom_ice_obudget`, `fesom_ice_flooding`, `fesom_ice_cut_off`, `fesom_ice_tfrez` |
| `fesom_ice_coupling.{c,h}` | `ice_oce_coupling.F90` | `fesom_ocean2ice`, `fesom_ice_oce_fluxes`, `fesom_ice_oce_fluxes_mom` |

### Patterns to follow (already established in the C port)

- **One C file per Fortran file.** Mirror subroutine names verbatim (with `fesom_` prefix). Per `feedback_port_fidelity.md`.
- **Halo-write audit on every C write-loop.** If Fortran loops `myDim+eDim` on a write, the C port must too. Per `feedback_write_loops_halo.md`.
- **Env-knob bisection.** Each new subsystem gets a `FESOM_NO_ICE_*=1` knob to disable it, mirroring `FESOM_NO_TRADV` etc.
- **Build via `bash -l configure.sh`.** Never mambaforge. Per `feedback_levante_build.md`.
- **Multi-rank validation on SLURM only.** Login node forbidden for CORE2.

### Key Fortran invariants that must be preserved

- `T_ICE.data(1:3)` = (a_ice, m_ice, m_snow). Allocated `myDim+eDim`.
- `whichEVP=0` path uses element-based stress (`sigma11/12/22`, `eps11/12/22` are `elem_size`) with node-based velocity rhs.
- `evp_rheol_steps = 120` subcycles per ice timestep.
- `ice_ave_steps = 1` → ice timestep = ocean timestep.
- `mesh%bc_index_nod2D` (allocated inside `ice_init`, line 889): set to 1 on interior nodes, 0 on open-boundary edge endpoints. Required even for `whichEVP=0`. Must be ported.
- `MOD_ICE.F90:743-754` allocates `uice_aux/vice_aux` (and `alpha_evp_array/beta_evp_array` for aEVP) ONLY when `whichEVP /= 0`. We do not allocate these.
- Non-ICEPACK `oce_fluxes` (line 393–398): `heat_flux = -net_heat_flux; water_flux = -fresh_wa_flux`. **No `-runoff` term** — runoff is already in `fresh_wa_flux` from `therm_ice`.

## Development Approach

### Testing approach: validation-driven, not TDD

This codebase does not have a unit-test framework (the existing 30+ `fesom_*.{c,h}`
files have no companion `_test.c` files). Validation is:

1. **Build cleanly** after each task (`bash -l configure.sh`).
2. **Multi-rank correctness** at end of each phase: run `job_core2_16_fix` (200 steps,
   16 ranks — fastest reproducer of partition-specific bugs) and confirm no blowup,
   no halo-bug signature, T/S within reasonable bounds.
3. **Halo-write audit** after every new `.c` file lands: grep
   `for.*<.*myDim_(nod|elem)2D;` against the matching Fortran loop bound.
4. **Env-knob bisection** is the debugging tool when validation fails. Add knobs
   per phase before they're needed.

The plan therefore has no `write tests for X` checkboxes — instead each phase
ends with explicit **build + halo audit + 16-rank validation run + check
acceptance criterion**. Treat that triplet with the same rigor as a passing
test suite: do not start the next phase until it passes.

### Other rules

- Complete each task fully before moving to the next.
- Make small, focused changes — one Fortran subroutine → one C function per task where reasonable.
- **CRITICAL: literal Fortran→C port** — no simplification. Mirror loop bounds, mirror per-element vs. per-node accounting, preserve area/areasvol ratios even when 1.0.
- **CRITICAL: halo-write audit before phase-end run** — this saves the days of debugging that the ocean port took to find three of these.
- **CRITICAL: update this plan file when scope changes** — add ➕ for newly discovered tasks, ⚠️ for blockers, mark `[x]` immediately when done.
- Run all CORE2 validation through SLURM. Login node only for build, snapshot reads, and pi-mesh smoke tests.

## Testing Strategy

### Per-task: build + audit
- `bash -l configure.sh` from `/home/a/a270088/port2/fesom2_port/` must succeed with no warnings introduced.
- For any task that adds a write-loop: grep the new loop bounds against the matching Fortran loop. Document the audit result in the task checkbox notes.

### Per-phase: 16-rank correctness
- `sbatch job_core2_16_fix` (200 steps, 16 ranks, ~5 min wallclock).
- Compare T/S min/max against the pre-phase baseline. Drift outside historical bounds = block on next phase until investigated.
- For phases B–E (where physics changes): the comparison is also against the Fortran reference run if available (`/home/a/a270088/fesom27/fesom2/work_core/`).

### Phase-E acceptance: full model month
- `sbatch job_core2_256_month` (5184 steps, 256 ranks, ~4 min compute).
- **Acceptance**: T.min stays above ≈ −2 °C (sea-ice insulation working), S.min and S.max within ±0.5 of pre-port baseline, no rank-dependent drift between 8/16/32/256.

## Progress Tracking
- mark completed items with `[x]` immediately when done
- add newly discovered tasks with ➕ prefix
- document issues/blockers with ⚠️ prefix
- update plan if implementation deviates from original scope
- keep plan in sync with actual work done

## What Goes Where
- **Implementation Steps** (`[ ]` checkboxes): code changes in `/home/a/a270088/port2/fesom2_port/`, halo audits, build runs, snapshot validation, multi-rank SLURM runs.
- **Post-Completion** (no checkboxes): physics validation against Fortran reference, memory updates, plan archival.

## Implementation Steps

### Phase A — State allocation, snapshot I/O, no-op driver

Goal: `fesom_ice` struct exists, allocates correctly at every rank count,
snapshot round-trips, `fesom_ice_step` is a no-op gated by both env knobs on.
No physics yet.

#### Task A1: Define `struct fesom_ice` and supporting types

**Files:**
- Create: `fesom2_port/src/fesom_ice_types.h`

- [x] declare `struct fesom_ice_data { double *values, *values_old, *values_rhs, *values_div_rhs, *dvalues, *valuesl; int id; }` mirroring `T_ICE_DATA`
- [x] declare `struct fesom_ice_work` with `fct_tmax`, `fct_tmin`, `fct_plus`, `fct_minus`, `fct_fluxes[3*elem_size]`, `fct_massmatrix`, `sigma11/12/22`, `eps11/12/22`, `ice_strength`, `inv_areamass`, `inv_mass`
- [x] declare `struct fesom_ice_thermo` with `t_skin`, `thdgr`, `thdgrsn`, `thdgra`, `thdgr_old`, `ustar` plus the scalar constants from `T_ICE_THERMO` lines 53–84 (rhoair, rhowat, rhoice, rhosno, cpair, cpice, cpsno, cc, cl, clhw, clhi, tmelt, boltzmann, etc.). **Skip** `dyngr/dyngrsn/dyngra` — CMIP6 output-only diagnostics, see Out-of-scope.
- [x] declare `struct fesom_ice` with `uice/vice/uice_rhs/vice_rhs/uice_old/vice_old`, stress fields, `srfoce_*` fields, `flx_fw`, `flx_h`, `h_ice`, `h_snow`, `data[3]`, embedded `work` and `thermo`, scalar params (pstar, ellipse, c_pressure, delta_min, evp_rheol_steps, ice_gamma_fct, ice_diff, theta_io, cd_oce_ice, ice_dt, Tevp_inv, whichEVP, ice_ave_steps)
- [x] hard-default `whichEVP = 0`, `evp_rheol_steps = 120`, `ice_ave_steps = 1`, `pstar = 30000.0`, `ellipse = 2.0`, `c_pressure = 20.0`, `delta_min = 1e-11`, `cd_oce_ice = 5.5e-3`, `ice_gamma_fct = 0.25`, `ice_diff = 10.0` (CORE2 namelist values)
- [x] **NO** `uice_aux/vice_aux/alpha_evp_array/beta_evp_array` fields (whichEVP=0 only)

#### Task A2: Implement `fesom_ice_init` (allocator + namelist defaults)

**Files:**
- Create: `fesom2_port/src/fesom_ice.c`
- Create: `fesom2_port/src/fesom_ice.h`

- [x] declare `fesom_ice_init(struct fesom_ice *ice, fesom_partit *p, struct fesom_mesh *m)` mirroring `MOD_ICE.F90:572` (skip namelist read — hardcode CORE2 defaults per FRESH_START.md §1)
- [x] allocate every field in `struct fesom_ice` with size `myDim+eDim` for nodes, `myDim+eDim+eXDim` for elements (matching `MOD_ICE.F90:712-754`); calloc-zero is fine since Fortran init sets everything to 0
- [x] set `ice->Tevp_inv = ice->evp_rheol_steps / ice->ice_dt` (computed in `ice_setup` line ~85 of ice_setup_step.F90 — verify exact source)
- [x] allocate `mesh->bc_index_nod2D[myDim+eDim]`, set to 1.0; loop `myDim_edge2D` and zero entries on open-boundary nodes per `MOD_ICE.F90:889-895`
- [x] declare `fesom_ice_free(struct fesom_ice *ice)` and free everything in reverse order
- [x] **halo audit**: bc_index_nod2D loop is per-edge in Fortran (`myDim_edge2D`); mirror exactly

#### Task A2b: Implement `fesom_ice_setup` (one-time post-init)

**Files:**
- Modify: `fesom2_port/src/fesom_ice.c`
- Modify: `fesom2_port/src/fesom_ice.h`

- [x] declare `fesom_ice_setup(struct fesom_ice *ice, fesom_partit *p, struct fesom_mesh *m)` mirroring `ice_setup_step.F90:51-91`
- [x] compute the three time constants in order (Fortran lines 74-78):
      `ice_dt = ice_ave_steps * fesom_phase1_dt;`
      `Tevp_inv = 3.0 / ice_dt;`
      `Clim_evp = Clim_evp * (evp_rheol_steps/ice_dt)^2 / Tevp_inv;`
- [x] **stub call** to `fesom_ice_mass_matrix_fill(ice, p, m)` — the actual implementation lands in Task E1; for Phase A leave as `/* TODO Phase E */` (no-op)
- [x] do **not** call `ice_initial_state` — covered implicitly by calloc-zero per Out-of-scope



**Files:**
- Modify: `fesom2_port/src/fesom_ice.c`
- Modify: `fesom2_port/src/fesom_ice.h`

- [x] declare `fesom_ice_step(int step, struct fesom_ice *ice, fesom_partit *p, struct fesom_mesh *m)` mirroring the structure of `ice_setup_step.F90:96` (ice_timestep)
- [x] read env vars `FESOM_NO_ICE_DYN`, `FESOM_NO_ICE_ADV`, `FESOM_NO_ICE_THERMO` once at startup (cache in static int)
- [x] driver body: `if (!FESOM_NO_ICE_DYN) /* TODO Phase D */ ; if (!FESOM_NO_ICE_ADV) /* TODO Phase E */ ; if (!FESOM_NO_ICE_THERMO) /* TODO Phase B */ ;` (all stubs for now)
- [x] post-step diagnostic: `h_ice[i] = data[1].values[i] / max(data[0].values[i], 1e-3)`, same for `h_snow` (matches `ice_setup_step.F90:319-322`)
- [x] **halo audit**: post-step diagnostic loop in Fortran is `1, myDim_nod2D+eDim_nod2D` — mirror

#### Task A4: Wire `fesom_ice` into main + step driver

**Files:**
- Modify: `fesom2_port/src/fesom_main.c` (allocate `struct fesom_ice` after mesh ready; free on shutdown)
- Modify: `fesom2_port/src/fesom_step.c` (call `fesom_ice_step` after the ocean step, mirroring placement in Fortran `oce_timestep_ale` which calls ice after the dynamics+tracer step)
- Modify: `fesom2_port/src/fesom_step.h` (extend `fesom_step_ctx` if needed to carry `struct fesom_ice *`)

- [x] add `struct fesom_ice *ice` to `fesom_step_ctx` (or pass as additional arg — match existing convention)
- [x] in `fesom_main`: after mesh + partit + aux + dyn allocated, call `fesom_ice_init` then `fesom_ice_setup`
- [x] in `fesom_main`: on shutdown, call `fesom_ice_free`
- [x] in `fesom_step.c`: call `fesom_ice_step(step_n, ctx->ice, ctx->partit, mesh)` at end of timestep (after final tracer commit)
- [x] confirm timestep order matches `oce_ale.F90` placement of ice driver

#### Task A5: Snapshot writer extended with ice fields

**Files:**
- Modify: `fesom2_port/src/fesom_io.c`
- Modify: `fesom2_port/src/fesom_io.h`

**Scope note**: the existing snapshot machinery is write-only — there is no
reader in the codebase yet, so a "round-trip" test isn't possible at Phase A.
Restart/read is its own future project. For Phase A we just need ice fields
visible in the output for diagnostic confirmation.

- [x] extend `fesom_io_write_snapshot` signature with `const struct fesom_ice *ice` (NULL = no ice fields)
- [x] gather + write `a_ice` (=`data[0].values`), `m_ice`, `m_snow`, `uice`, `vice`, `h_ice`, `h_snow` — all per-node, single-level
- [x] update the call site in `fesom_main.c` to pass `&ice`
- [x] verify with `ncdump -h` of a Phase-A snapshot: 7 new ice variables present, all zero (since no physics runs yet)

#### Task A6: Phase A validation

- [x] `bash -l configure.sh` clean build
- [x] **halo audit**: grep `for.*<.*myDim_(nod|elem)2D;` in `fesom_ice.c`; every write-loop matches the Fortran bound
- [x] login node smoke test: 1-rank pi mesh, 10 steps, with both `FESOM_NO_ICE_*` knobs UNSET (default off so the no-op stubs run); confirm no crash, ice fields all zero in snapshot
- [x] `sbatch job_core2_16_fix` — 200 steps at 16 ranks; T/S **bit-identical** to pre-phase baseline (SLURM job 24486031, 2:48 wallclock, exit 0; per-step rank-0 signature diff = 0 lines vs. `/tmp/baseline_pre_phaseA.log`)
- [x] mark Phase A complete in plan

---

### Phase B — Thermodynamics (`fesom_ice_thermo.c`)

Goal: port non-ICEPACK thermodynamics. Once `therm_ice` produces `fresh_wa_flux`,
remove the runoff subtraction in `fesom_sss_runoff_step`. Single-cell `therm_ice`
is the heart — `thermodynamics` is just the per-node loop over it.

#### Task B1: Port `cut_off`, `flooding`, `TFrez`

**Files:**
- Create: `fesom2_port/src/fesom_ice_thermo.c`
- Create: `fesom2_port/src/fesom_ice_thermo.h`

- [x] port `TFrez(S)` (line 899 of ice_thermo_oce.F90) — pure function, freezing-point polynomial
- [x] port `flooding(ithermp, h, hsn)` (line 872) — single-cell, snow→ice conversion when ice freeboard goes negative
- [x] port `cut_off(ice, partit, mesh)` (line 70) — per-node loop applying hmin/Armin cutoffs; confirmed near-zero-zero rescaling (Fortran zeros all three tracers when a_ice or m_ice drops below 1e-9)
- [x] **halo audit**: cut_off Fortran loop bound is `myDim+eDim` (line 102); C port uses `N = myDim+eDim` — matches

#### Task B2: Port `obudget` (open-water surface energy budget)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_thermo.c`
- Modify: `fesom2_port/src/fesom_ice_thermo.h`

- [x] port `obudget(ithermp, qa, fsh, flo, t, ug, ta, ch, ce, geolon, geolat, fh, evap, ...)` (line 777) — single-cell open-water heat budget; computes evap, sensible+latent fluxes, longwave. Asserts `open_water_albedo == 0` (CORE2 default; Taylor/Briegleb branches need clock + zenith routine that aren't ported yet).
- [x] this is a pure function on inputs — no halo concerns

#### Task B3: Port `budget` (ice surface energy budget)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_thermo.c`
- Modify: `fesom2_port/src/fesom_ice_thermo.h`

- [x] port `budget(ithermp, hice, hsn, t, ta, qa, fsh, flo, ug, S_oc, ch_i, ce_i, fh, subli)` (line 657) — surface temperature solver under ice. Albedo selection (freezing/melting × snow/no-snow) preserved; `q1=11637800.0`, `q2=-5897.8` are saturation-over-ice constants (different from obudget's saturation-over-water).
- [x] iterative skin-temperature solver — fixed 5-iter Newton-Raphson (Fortran does no convergence check; just 5 iterations). t clamped to ≤ 0 after iteration. Final fluxes use converged t + last-iter B.

#### Task B4a: Port `therm_ice` core (basic single-class growth/melt)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_thermo.c`
- Modify: `fesom2_port/src/fesom_ice_thermo.h`

- [x] port the core scaffold of `therm_ice(ithermp, h, hsn, A, fsh, flo, Ta, qa, rain, snow, runo, rsss, ...)` (`ice_thermo_oce.F90:367`). Implemented the FULL `iclasses` loop (default 7) since it's only 8 lines and required for CORE2 parity; only the `new_iclasses=true` (h_cutoff/hpdf) sub-branch is skipped (no hpdf in struct).
- [x] preserve `t_skin` solver implicitly via per-class `fesom_ice_budget` calls; ice salinity flux via `rsf = fwice * Sice`; mass exchange via dhgrowth/dhsngrowth in fwice/fwsnw
- [x] **runoff fold IS already in B4a** — `prec = rain + runo + snow*(1-A)` (line 526) → `fw = prec + evap + fwice + fwsnw` (line 616). No separate runoff-fold step needed; the contract change at C3b enables this to reach `water_flux`.
- [ ] hold off on flooding integration (lines 637-649) — Task B4b. Currently `iflice=0` and salt-conservation correction skipped.
- [ ] hold off on `use_virt_salt=true` branch (lines 619-622) — Task B4b. Asserts hardcoded false in B4a.

#### Task B4b: Add lateral melt, snow→ice conversion, runoff fold

**Files:**
- Modify: `fesom2_port/src/fesom_ice_thermo.c`
- Modify: `fesom2_port/src/fesom_ice_thermo.h`

- [x] (already in B4a) lateral-melt is the `c_melt * fmin(rh, 0)` term in compactness update; not a separate branch
- [x] snow→ice conversion via `fesom_ice_flooding` when freeboard goes negative (lines 638-641 in Fortran). Salt-conservation correction also added (lines 645-649).
- [x] (already in B4a) multi-class loop preserved; only `new_iclasses=true` (h_cutoff/hpdf path) is skipped — out-of-scope.
- [x] (already in B4a) runoff fold is structural: `prec=rain+runo+snow*(1-A)` → `fw=prec+evap+fwice+fwsnw`. C3b will remove the parallel subtraction in fesom_sss_runoff.
- [x] **use_virt_salt branch added** (lines 619-622) as a function parameter; the call site in `fesom_ice_thermodynamics` passes `sr->use_virt_salt`.

#### Task B5: Port `thermodynamics` (per-node driver)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_thermo.c`
- Modify: `fesom2_port/src/fesom_ice_thermo.h`

- [x] port `fesom_ice_thermodynamics(ice, partit, mesh, forcing, jra, sr)` (mirrors Fortran line 148). Two-pass: ustar friction-velocity loop (myDim only, then halo-exchange) + per-node therm_ice loop (myDim+eDim). Ice tracers + thdgr/thdgrsn/thdgra/t_skin + flx_fw/flx_h written.
- [x] write `flx_fw[n]` and `flx_h[n]` per node
- [x] update `data[0..2].values` + `values_old`, `thdgr*`, `t_skin`
- [x] **halo audit**: ustar loop is `myDim` only (matches Fortran line 230); followed by `fesom_exchange_nod2D(ustar)` (matches Fortran line 238). Therm_ice loop is `myDim+eDim` (matches Fortran line 244). Both correct.
- [x] `fesom_ice_cut_off` called from `fesom_ice_step` after thermodynamics (mirrors ice_setup_step.F90:295)
- [x] **dependency**: required `forcing->Ch_atm_oce` and `Ce_atm_oce` per-node arrays — added to fesom_forcing struct; populated in fesom_bulk.c by storing the per-node ch/ce already computed by ncar_ocean_fluxes_mode (was being discarded).

#### Task B6: Wire thermodynamics into driver

**Files:**
- Modify: `fesom2_port/src/fesom_ice.c`
- Modify: `fesom2_port/src/fesom_step.c` (if needed)

- [x] in `fesom_ice_step`: extended signature to take `dyn`, `tracers`, `forcing`, `jra`, `sr` (jra+sr added to `fesom_step_ctx`); thermo block now calls `fesom_ice_thermodynamics` + `fesom_ice_cut_off`. NULL forcing/jra/sr silently disables thermo (allows pi smoke test without JRA).
- [x] confirm forcing inputs reachable — `jra->{shortwave,longwave,Tair,shum,prec_rain,prec_snow,u_wind,v_wind}` and `forcing->{runoff,Ch_atm_oce,Ce_atm_oce}` all wired.
- [x] populate `srfoce_*` from ocean state at start of `fesom_ice_step` — surface T/S from `tracers->data[T/S].values[n*nl+0]`, surface u/v from `dyn->uvnode[n*nl*2+0..1]`, ssh from `dyn->eta_n[n]`. Loop covers myDim+eDim.

#### Task B7: Phase B validation

- [x] `bash -l configure.sh` clean build
- [x] **halo audit**: only myDim-only loop in fesom_ice_thermo.c is the ustar friction velocity (line 435), matches Fortran ice_thermo_oce.F90:230. Followed by `fesom_exchange_nod2D(ustar)` matching Fortran's `call exchange_nod`.
- [x] `sbatch job_core2_16_fix` (job 24486349, 2:48 wallclock, exit 0); 200 steps at 16 ranks, no NaN, no blowup
- [x] **Acceptance MET**: 23,310 / 126,858 nodes have nonzero ice at step 200 (Arctic + Antarctic polar caps); max a_ice=0.97, max m_ice=1.62 m (realistic mid-winter Arctic thickness)
- [x] **Note on freshwater confirmed**: Ocean T/S signature is bit-identical to the Phase A baseline at step 200, confirming `flx_fw`/`flx_h` are written into the ice struct but NOT yet plumbed into `water_flux`/`heat_flux` (that's Phase C2). The runoff subtraction in `fesom_sss_runoff` correctly remains in place through end of Phase B; C3b will remove it once C2 lands.
- [x] mark Phase B complete in plan

---

### Phase C — Ocean→ice + ice→ocean heat/water coupling (`fesom_ice_coupling.c`, no momentum yet)

Goal: get ocean state into the ice and get heat/water fluxes back into the
ocean. The momentum-side coupling (`oce_fluxes_mom`) is **deferred to Phase D**
because it reads `u_ice` to compute ice-ocean drag — meaningless until EVP runs.

#### Task C1: Port `ocean2ice`

**Files:**
- Create: `fesom2_port/src/fesom_ice_coupling.c`
- Create: `fesom2_port/src/fesom_ice_coupling.h`

- [x] port `fesom_ocean2ice(ice, dyn, tracers, partit, mesh)` mirroring `ice_oce_coupling.F90:154`. T_oc/S_oc/elevation from surface tracer + hbar (myDim+eDim loop). u_w/v_w via area-weighted mean over surrounding elements at level 0 (myDim only loop, then halo exchange). Cavity skip via ulevels_nod2D > 1. Replaces the simplified copy in B6.
- [x] **halo audit**: u_w/v_w computation is myDim only (matches Fortran line 219); followed by `fesom_exchange_nod2D(u_w/v_w)` (matches line 244)
- [x] halo exchanges of u_w, v_w added after the myDim-only computation

#### Task C2: Port `oce_fluxes` (non-ICEPACK branch — heat/water only)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_coupling.c`
- Modify: `fesom2_port/src/fesom_ice_coupling.h`

- [x] port `fesom_ice_oce_fluxes(ice, partit, mesh, tracers, forcing, sr)` non-ICEPACK branch: `heat_flux = -flx_h`, `water_flux = -flx_fw`, NO `-runoff` term
- [x] virtual_salt block (lines 436-457) ported with use_virt_salt branch + global mean balance
- [x] relax_salt block (lines 498-528) ported with global mean balance and SSS restoring; runs unconditionally now (mirrors Fortran)
- [x] cavity / wiso / iceberg branches skipped (out of scope)
- [x] **halo audit**: virtual_salt write loop is myDim+eDim, balance loop is myDim only (matches Fortran 452-456); same for relax_salt. Halo exchanges of virtual_salt, relax_salt, heat_flux, water_flux added.

#### Task C3: Wire ocean2ice + oce_fluxes into driver

**Files:**
- Modify: `fesom2_port/src/fesom_ice.c`
- Modify: `fesom2_port/src/fesom_step.c`

- [x] in `fesom_ice_step`: replaced simplified ocean2ice copy with `fesom_ocean2ice` call at start; `fesom_ice_oce_fluxes` call added at end (after thermo+cut_off)
- [x] **moved** `fesom_ice_step` call from inside `fesom_timestep` to `fesom_main.c` BEFORE `fesom_timestep`. This matches Fortran flow (ice runs before ocean step in oce_timestep_ale) and ensures the ice-aware heat_flux/water_flux are seen by impl_vert_diff_tracers.
- [x] do not wire `oce_fluxes_mom` yet (Phase D5 — needs uice from EVP)

#### Task C3b: Remove runoff subtraction from `fesom_sss_runoff` (the contract handover)

**Files:**
- Modify: `fesom2_port/src/fesom_sss_runoff.c`
- Update: `~/.claude/projects/-home-a-a270088-port2/memory/feedback_runoff_fold_when_ice.md`

- [x] contract verified: B4a confirmed runoff is structurally in `flx_fw` via `prec=rain+runo+snow*(1-A)`; C2 confirmed `water_flux = -flx_fw`
- [x] runoff subtraction block in `fesom_sss_runoff_step` removed (lines 360-368 → comment citing contract)
- [x] SSS-restoring + virtual_salt blocks left in place; oce_fluxes overwrites them anyway (harmless duplication; cleaner refactor deferred)
- [x] one-line citation comment added pointing to feedback_runoff_fold_when_ice.md
- [ ] memory file update — pending (will do at Phase G/end)

#### Task C4: Phase C validation

- [x] `bash -l configure.sh` clean build
- [x] **halo audit**: 4 myDim-only loops in fesom_ice_coupling.c (integrate, ocean2ice u_w/v_w, virtual_salt balance, relax_salt balance); all match Fortran loop bounds. Each followed by appropriate halo exchange.
- [x] `sbatch job_core2_16_fix` (job 24486565, 2:48 wallclock, exit 0); 200 steps at 16 ranks completed cleanly
- [x] **T.min jumped from −5.64°C to −2.06°C** — ocean now correctly clamped at the freezing point by sea-ice insulation (the headline expected result of ice→ocean coupling). T.max and S.max unchanged in tropics (no ice). S.min dropped 5.65→5.60 from ice melt.
- [x] **flux signature now ice-aware**: heat_flux max 1.44e+03 → 4.63e+03 W/m² (ice latent + radiation), water_flux max 2.47e−6 → 1.12e−5 m/s (runoff folded properly). No NaN.
- [x] ice growth bounded to physical levels (m_ice max 0.27 m vs. Phase B's unphysical 1.62 m — feedback works)
- [x] surface stress is unchanged from Phase B (still bulk-direct — Phase D5 fixes that)
- [x] mark Phase C complete

---

### Phase D — Standard EVP dynamics (`fesom_ice_evp.c`)

Goal: port the standard EVP solver. **`whichEVP=0` only** — stub the dispatcher
to abort if `whichEVP != 0`. mEVP/aEVP can be added later if needed.

#### Task D1: Port `stress_tensor`

**Files:**
- Create: `fesom2_port/src/fesom_ice_evp.c`
- Create: `fesom2_port/src/fesom_ice_evp.h`

- [ ] port `stress_tensor(ice, partit, mesh)` (`ice_EVP.F90:46`) — element loop computing `eps11/12/22` from velocity, then EVP stress update
- [ ] preserve `det1/det2/dte/zeta/delta/delta_inv/r1/r2/r3` exact algebra
- [ ] preserve `ulevels(el) > 1` cavity skip and `ice_strength(el) > 0` guard
- [ ] **halo audit**: Fortran loop is `1, myDim_elem2D` — verify and match

#### Task D2: Port `stress2rhs`

**Files:**
- Modify: `fesom2_port/src/fesom_ice_evp.c`
- Modify: `fesom2_port/src/fesom_ice_evp.h`

- [ ] port `stress2rhs(ice, partit, mesh)` (line 177) — divergence of stress tensor scattered to nodes; produces `uice_rhs`, `vice_rhs`
- [ ] **halo audit**: scatter loops — likely Fortran does element-loop with node-write; check whether final node loop is `myDim` or `myDim+eDim`

#### Task D3: Port `EVPdynamics` (subcycled solver)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_evp.c`
- Modify: `fesom2_port/src/fesom_ice_evp.h`

- [ ] port `EVPdynamics(ice, partit, mesh)` (line 330) — main subcycled loop, `evp_rheol_steps=120` iterations of `stress_tensor` + `stress2rhs` + velocity update
- [ ] preserve halo exchanges of `uice/vice` (or whatever Fortran exchanges) inside the subcycle; this is where MPI cost lives
- [ ] preserve ocean-ice drag formulation using `srfoce_u/v/ssh`
- [ ] **halo audit**: critical here — every write inside the subcycle must respect Fortran's halo bounds, since 120 iterations amplify any halo gap

#### Task D4: Stub the dispatcher with abort for whichEVP != 0

**Files:**
- Modify: `fesom2_port/src/fesom_ice.c`

- [ ] in `fesom_ice_step` dynamics block: `if (ice->whichEVP == 0) fesom_ice_evp_dynamics(ice, p, m); else { fprintf(stderr, "fesom_ice: whichEVP=%d not supported (only standard EVP)\n", ice->whichEVP); MPI_Abort(MPI_COMM_WORLD, 1); }`
- [ ] confirm `MOD_ICE.F90:743-754` aux allocations are NOT done in our `fesom_ice_init` (defensive — they'd leak otherwise)

#### Task D5: Port `oce_fluxes_mom` (ice-aware surface stress — was C2, deferred)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_coupling.c`
- Modify: `fesom2_port/src/fesom_ice_coupling.h`
- Modify: `fesom2_port/src/fesom_step.c`

- [ ] port `oce_fluxes_mom(ice, dynamics, partit, mesh)` (`ice_oce_coupling.F90:53`) — element-based wind stress now blended with ice-ocean drag from `cd_oce_ice` and ice velocity (`ice%uice/vice` populated by Task D3 EVP)
- [ ] this REPLACES the current direct write of bulk `stress_surf` into `dyn->stress_surf`
- [ ] **halo audit**: per-element loop bounds; Fortran reads `u_ice/v_ice` at element vertices, so it depends on those being halo-exchanged (which D3 should already do)
- [ ] in `fesom_step.c`: call `fesom_ice_oce_fluxes_mom` before the momentum RHS step (`compute_vel_rhs` consumes `stress_surf`)
- [ ] **placement rationale**: this lands AFTER D3 so EVP has run and `u_ice` is populated; landing it in Phase C against `u_ice ≡ 0` would silently produce wrong (open-water-only weighted) drag

#### Task D6: Phase D validation

- [ ] `bash -l configure.sh` clean build
- [ ] **halo audit**: extra-careful pass on `fesom_ice_evp.c` and on the `oce_fluxes_mom` write-loops; the 120-subcycle inner loop is the highest-risk place in the port for multi-rank divergence; apply the secondary `grep -v eDim` filter from the recipe in Technical Details
- [ ] `sbatch job_core2_16_fix` — confirm ice velocities are non-zero in polar ice cover (max < 0.5 m/s), no NaN, no blowup
- [ ] surface stress check: in Arctic ice cover, `|stress_surf|` should be visibly smaller than the Phase C baseline (because ice now mediates) — quantify with a Phase-C-vs-Phase-D printout at one polar element
- [ ] **multi-rank cross-check** (concrete tolerance): `job_core2_8_fix` vs `job_core2_32_fix` vs `job_core2_16_fix` at step 200 — `|Δuice|_∞ < 1e-8` (EVP is not bitwise reproducible across ranks because of MPI reductions in the subcycle, so a small floor is expected; > 1e-6 is a halo bug)
- [ ] add `FESOM_PRINT_EVERY=10` and watch ice velocity max — should be < 0.5 m/s in ice cover, exactly 0 in ice-free regions
- [ ] mark Phase D complete

---

### Phase E — Ice FCT advection (`fesom_ice_fct.c`)

Goal: port FCT advection of (a_ice, m_ice, m_snow). This is the highest-risk
phase for halo bugs — the ocean FCT had two of these (`init_tracers_AB` and
`fct_ttf_max/min`). Apply the audit *proactively* rather than reactively.

#### Task E1: Port `ice_mass_matrix_fill` (one-time setup)

**Files:**
- Create: `fesom2_port/src/fesom_ice_fct.c`
- Create: `fesom2_port/src/fesom_ice_fct.h`

- [ ] port `ice_mass_matrix_fill(ice, partit, mesh)` (`ice_fct.F90:1145`) — builds `ice%work%fct_massmatrix` once at startup
- [ ] call from `fesom_ice_setup` (one-time)
- [ ] **halo audit**: this is a setup routine — match Fortran exactly

#### Task E2: Port `ice_TG_rhs` and `ice_TG_rhs_div`

**Files:**
- Modify: `fesom2_port/src/fesom_ice_fct.c`
- Modify: `fesom2_port/src/fesom_ice_fct.h`

- [ ] port `ice_TG_rhs(ice, partit, mesh)` (line 91) — Taylor-Galerkin RHS for high-order solution
- [ ] port `ice_TG_rhs_div(ice, partit, mesh)` (line 1253) — divergence form RHS
- [ ] **halo audit — proactive**: grep Fortran for write-loops bounded `myDim+eDim` on `values_rhs`, `values_div_rhs`, `dvalues`, `valuesl` — any such loop in C must use `myDim+eDim`

#### Task E3: Port `ice_solve_low_order` and `ice_solve_high_order`

**Files:**
- Modify: `fesom2_port/src/fesom_ice_fct.c`
- Modify: `fesom2_port/src/fesom_ice_fct.h`

- [ ] port `ice_solve_low_order(ice, partit, mesh)` (line 243) — first-order upwind reference solution into `valuesl`
- [ ] port `ice_solve_high_order(ice, partit, mesh)` (line 347) — Taylor-Galerkin high-order solution
- [ ] **halo audit**: this is the equivalent of ocean's `oce_tra_adv_fct` step a1; the bug was that `fct_ttf_max/min` halo was not filled. Mirror Fortran loop bounds exactly. **Pre-emptively** check Fortran for `myDim+eDim` write-loops in this region

#### Task E4: Port `ice_fem_fct` (Zalesak limiter)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_fct.c`
- Modify: `fesom2_port/src/fesom_ice_fct.h`

- [ ] port `ice_fem_fct(tr_array_id, ice, partit, mesh)` (line 498) — per-tracer Zalesak limiter (called 3× per step for a_ice, m_ice, m_snow)
- [ ] **halo audit — proactive**: this is precisely where the ocean FCT bug 2 lived; pre-emptively check every loop bound in Fortran and mirror in C

#### Task E5: Port `ice_fct_solve` (driver)

**Files:**
- Modify: `fesom2_port/src/fesom_ice_fct.c`
- Modify: `fesom2_port/src/fesom_ice_fct.h`

- [ ] port `ice_fct_solve(ice, partit, mesh)` (line 210) — sequences `ice_solve_low_order` → `ice_solve_high_order` → 3× `ice_fem_fct` → write back to `data[1..3].values`
- [ ] confirm halo exchange placement matches Fortran

#### Task E6: Wire FCT into driver

**Files:**
- Modify: `fesom2_port/src/fesom_ice.c`

- [ ] in `fesom_ice_step` advection block (between dynamics and thermodynamics): call `fesom_ice_tg_rhs` then `fesom_ice_fct_solve` then `fesom_ice_cut_off`, mirroring `ice_setup_step.F90:258-295`

#### Task E7: Halo-write audit pass on the entire `fesom_ice_fct.c`

- [ ] grep `for.*<.*myDim_(nod|elem)2D;` in `fesom_ice_fct.c`; for each match, check matching Fortran loop bound
- [ ] document each write-loop in a comment: `/* matches ice_fct.F90:NNN — Fortran bound: myDim or myDim+eDim */`
- [ ] this is the most important task in Phase E — saves the days of debugging that the ocean FCT halo bugs took

#### Task E8: Phase E validation

- [ ] `bash -l configure.sh` clean build
- [ ] `sbatch job_core2_16_fix` — confirm ice tracers a_ice/m_ice/m_snow advect without blowup
- [ ] **multi-rank cross-check**: 8/16/32-rank runs at 200 steps must produce ice fields within numerical noise of each other; if 16-rank diverges from 8/32, that is the same partition-dependent halo signature we hit twice in ocean FCT — debug with `FESOM_NO_ICE_THERMO=1` then `FESOM_NO_ICE_DYN=1` to localize
- [ ] `sbatch job_core2_256_month` — full model month — **acceptance criterion below**
- [ ] mark Phase E complete

---

### Task F: Verify acceptance criteria

- [ ] T.min stays above ≈ −2 °C in 256-rank model month (was −27 °C without ice)
- [ ] S.min and S.max within ±0.5 of pre-port baseline at month-end
- [ ] No partition-dependent drift between 8/16/32/256 ranks at step 200
- [ ] Ice extent at end of model month qualitatively matches Fortran reference for January (Arctic ice cover, Antarctic Marginal ice zone)
- [ ] All `FESOM_NO_ICE_*` env knobs work (each can disable its subsystem cleanly)
- [ ] Snapshot round-trips ice state correctly
- [ ] Build still clean; no warnings introduced

### Task G: Update memory + archive plan

- [ ] update `feedback_runoff_fold_when_ice.md` to mark resolved (runoff now folded by `fesom_ice_thermodynamics`)
- [ ] add new `project_sea_ice_port_state.md` memory summarizing where the sea-ice port landed, residuals if any, and any new env knobs added
- [ ] update `MEMORY.md` index to point at the new memory
- [ ] move this plan to `fesom2_port/docs/plans/completed/`

## Technical Details

### `struct fesom_ice` shape (mirrors `T_ICE`)

```c
struct fesom_ice {
    /* tracers — a_ice, m_ice, m_snow */
    struct fesom_ice_data data[3];

    /* velocity (nodes) */
    double *uice, *uice_rhs, *uice_old;
    double *vice, *vice_rhs, *vice_old;

    /* surface stresses (nodes) */
    double *stress_atmice_x, *stress_iceoce_x;
    double *stress_atmice_y, *stress_iceoce_y;

    /* ocean surface state seen by ice (nodes) */
    double *srfoce_temp, *srfoce_salt, *srfoce_ssh;
    double *srfoce_u,    *srfoce_v;

    /* fluxes from thermodynamics (nodes) */
    double *flx_fw;  /* fresh_wa_flux — includes runoff */
    double *flx_h;   /* net_heat_flux */

    /* diagnostic ice/snow thicknesses (nodes) */
    double *h_ice, *h_snow;

    /* embedded work and thermo */
    struct fesom_ice_work   work;
    struct fesom_ice_thermo thermo;

    /* scalar parameters */
    double pstar, ellipse, c_pressure, delta_min, Clim_evp, zeta_min;
    int    evp_rheol_steps;
    double ice_gamma_fct, ice_diff, theta_io, cd_oce_ice;
    double ice_dt, Tevp_inv;
    int    whichEVP;       /* hard 0 — abort if changed */
    int    ice_ave_steps;  /* hard 1 */
};
```

### Driver order inside `fesom_ice_step` (mirrors `ice_setup_step.F90:96`)

```
if (!FESOM_NO_ICE_DYN)     fesom_ice_evp_dynamics(ice, p, m);
if (!FESOM_NO_ICE_ADV)   { fesom_ice_tg_rhs(ice, p, m); fesom_ice_fct_solve(ice, p, m); fesom_ice_cut_off(ice, p, m); }
if (!FESOM_NO_ICE_THERMO)  fesom_ice_thermodynamics(ice, p, m);
/* post-step diagnostic */
for (n = 0; n < myDim+eDim; n++) {
    h_ice [n] = data[1].values[n] / max(data[0].values[n], 1e-3);
    h_snow[n] = data[2].values[n] / max(data[0].values[n], 1e-3);
}
```

### Step-driver wiring (in `fesom_step.c`)

```
fesom_ice_oce_fluxes_mom(ice, ...);  /* Phase D5 — ice-aware stress, needs uice from prev EVP step */
... existing ocean step (compute_vel_rhs consumes stress_surf) ...
fesom_ice_step(step, ice, ...);      /* Phase A driver — runs ocean2ice, EVP, FCT, thermo */
fesom_ice_oce_fluxes(ice, ...);      /* Phase C — produces heat/water for next ocean step */
```

The `oce_fluxes_mom` call lives at the *top* of the next step (it consumes
`uice` from the previous `fesom_ice_step`), not in the middle. This matches
Fortran `oce_timestep_ale` order: ice runs at the end of the ocean step,
its outputs feed the next ocean step.

### Env knobs

| Knob | Effect |
|---|---|
| `FESOM_NO_ICE_DYN=1` | skip `fesom_ice_evp_dynamics` |
| `FESOM_NO_ICE_ADV=1` | skip `fesom_ice_tg_rhs` + `fesom_ice_fct_solve` + `fesom_ice_cut_off` |
| `FESOM_NO_ICE_THERMO=1` | skip `fesom_ice_thermodynamics` (and don't update `flx_fw/flx_h` — `oce_fluxes` will see zeros, ocean falls back to bulk-only) |

When all three are set, the C ice driver is a pure no-op; the ocean port should
produce bit-identical results to the pre-sea-ice baseline.

### whichEVP=0 enforcement

- Do not allocate `uice_aux/vice_aux/alpha_evp_array/beta_evp_array` (the Fortran
  guard at `MOD_ICE.F90:743-754` is mirrored by simple omission).
- Dispatcher in `fesom_ice_step` aborts loudly if `ice->whichEVP != 0`. The user
  asked specifically that whichEVP=0 land in the standard `EVPdynamics` path
  (file `fesom_ice_evp.c`), never in any mEVP/aEVP routine — there's no
  mEVP/aEVP code in the port at all, so this is structural rather than dispatch.

### Halo-write audit recipe (apply at end of every phase)

```bash
# step 1 — list all per-myDim write-loops in the new files
grep -nE 'for.*<.*myDim_(nod|elem)2D;' src/fesom_ice*.c

# step 2 — isolate the suspects (myDim-only, no +eDim extension)
grep -nE 'for.*<.*myDim_(nod|elem)2D;' src/fesom_ice*.c | grep -v eDim

# step 3 — for each suspect line, open the matching Fortran subroutine and
# find the same loop. If Fortran says `do n=1, myDim_*+eDim_*` (or
# `+eXDim_elem2D`), the C port MUST use the same extended bound.
```

Any divergence is a partition-dependent bug waiting to fire at some rank count.
Document each audited loop with an inline comment citing the Fortran source
line and the agreed bound, e.g.:

```c
/* matches ice_fct.F90:NNN — Fortran loop bound is myDim_nod2D, no halo write */
```

## Post-Completion

**Validation against Fortran reference** (manual, after Phase E):
- run Fortran reference for 1 model month at CORE2 (`/home/a/a270088/fesom27/fesom2/work_core/`)
- compare ice extent maps (a_ice ≥ 0.15 contour) — should match Fortran qualitatively
- compare T.min, S.min, S.max time series — within 1°C / 0.5 PSU
- compare integrated freshwater flux over a year — within 1%

**Memory updates** (Task G above): record port outcome and residuals.

**Possible follow-ups** (out of scope here, file as future plans):
- Port `ice_maEVP.F90` (mEVP, aEVP) — only when a run actually needs `whichEVP != 0`.
- Port `ice_meltponds.F90` — only when `use_meltponds = .true.` is needed.
- Port cavity ice clean-up routines (`cavity_ice_clean_vel`, `cavity_ice_clean_ma`) — only when cavity is enabled.
