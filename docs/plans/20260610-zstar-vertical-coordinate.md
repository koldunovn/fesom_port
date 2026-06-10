# Port the zstar (z*) vertical coordinate (which_ALE='zstar')

> Plan 1 of 2 from the 2026-06-10 brainstorm (companion: `20260610-tke-vertical-mixing.md`).
> **Revised 2026-06-10 after plan-review.** Corrections folded in: `real_salt_flux` is a LIVE
> producer in the standalone path under zstar (the C currently discards `rsf` — wiring it is a
> required fix, not optional plumbing); there is NO "non-linfs hpressure" term (under zstar
> Fortran computes no hpressure at all; the zxxxx shchepetkin PGF is self-contained); the S
> surface-BC is `+dt·(virtual_salt+relax_salt+real_salt_flux·is_nonlinfs)` (no sval·wf term);
> the freshwater/salt flip moved BEFORE the SSH/vert-vel phases (their dumps consume
> water_flux); the geometry audit greps the geometry arrays (not `linfs` comments) and adds
> mo_convect, vertical tracer advection (dual-geometry subtlety), and PP dz_inv.
> All design decisions are user-approved (2026-06-10): zstar-only scope behind a runtime
> switch, linfs byte-identical to tag `v1.0`, per-feature validation vs a NEW Fortran zstar
> reference, feeding the approved run matrix (per-feature + combined zstar+TKE 2yr + 5yr C).

## Overview

Port FESOM2's `which_ALE='zstar'` vertical coordinate from Fortran to C as a **strictly
faithful line-by-line port**, runtime-selected by `FESOM_ALE=linfs|zstar` (default `linfs`).
Under zstar the SSH change is distributed proportionally over the water column every step:
layer thicknesses `hnode`, element thicknesses `helem`, and the 3D level geometry
`zbar_3d_n`/`Z_3d_n` become time-varying; the SSH operator gains a per-step update; the
surface freshwater flux becomes a REAL volume flux (`use_virt_salt=.false.`) and the sea-ice
salt exchange becomes a REAL salt flux (`real_salt_flux`). This is the modern FESOM2
production coordinate and the prerequisite for the zstar+TKE combined target.

**GUIDING PRINCIPLE (non-negotiable):** THIS IS A PORT. Reproduce every formula, constant,
loop bound, branch, halo exchange, and ordering exactly. Where a value comes from a namelist,
port the **reference run's** value, not the module default (`feedback_namelist_over_codedefault`).
"Port what the reference uses" (`feedback_port_used_terms`): `zlevel` and its local-zstar
fallback (`lzstar_lev`/`min_hnode`/`CFL_z` guards, distribute/refill logic) are **NOT ported**
— no reference run exercises them. Cavity guards (`nzmin==1` / `ulevels>1 cycle`) are ported
as written even though C has no cavities.

**linfs stays byte-identical to `v1.0`** — verified by a byte-gate in EVERY phase that touches
shared code, not just at the end (the KPP K-phase pattern).

⚠️ **Line-number caveat:** `port2/fesom2` carries uncommitted dump-shim/debug lines (the
[FSSH] block, `dump_shim_record_*`, `FHBAR_DBG`), so oce_ale.F90 line numbers below can drift
by a few lines. Anchor by routine name + quoted code, not absolute line.

## Reference zstar configuration (VERIFIED 2026-06-10 — values to port)

The Fortran reference = `work_linfs_2yr` config with ONLY `which_ALE` flipped. Everything
else (KPP, dt=1800, MFCT, GM/Redi, EVP, JRA55, PHC) stays as in the v1.0-validated pair.

| Symbol | Value | Source | Consequence |
|---|---|---|---|
| `which_ALE` | `'zstar'` | namelist.config:71 (flipped from `'linfs'`) | the one knob |
| `use_partial_cell` | `.false.` | namelist.config:75 | full cells; `bottom_node/elem_thickness` = nominal `zbar` spacing |
| `which_pgf` | `'shchepetkin'` | **module default** `oce_modules.F90:190` — namelist.oce does NOT set it | zstar PGF = `pressure_force_4_zxxxx_shchepetkin`; (linfs ignores `which_pgf`: full-cell branch) |
| `use_virt_salt` | `.false.` (derived) | `gen_modules_forcing.F90:27`, set in `oce_setup_step.F90:115-125` | virtual_salt=0; real water-flux AND real salt-flux plumbing ON |
| `is_nonlinfs` | `1.0` (derived) | o_PARAM (set with use_virt_salt) | tracer surface-BC terms activate |
| `use_floatice` | `.false.` | namelist.ice (absent → default) | EVP `use_pice=0` → p_ice term ×0 (same in vel_rhs) |
| `alpha`, `theta` | as in namelist.oce (same as linfs runs) | namelist.oce | `factor=g·dt·alpha·theta` in stiff update |
| `i_vert_diff` | `.true.` | namelist.tra:74 | implicit TDMA path carries the surface BCs |
| `lzstar_lev`, `min_hnode` | irrelevant | zlevel-only machinery | NOT ported |
| `ldiag_DVD` | `.false.` | not enabled | `rescue_hnode_old`/hnode_new-restore in update_thickness NOT ported |
| `use_ssh_se_subcycl` | `.false.` | namelist.dyn | CG path only (as ported) |

➕ **Z0 sub-task:** re-confirm every value above against the namelists actually used for the
new `work_zstar_kpp` reference run; archive them in `docs/zstar_reference_namelists/` with a
PROVENANCE.md (KPP-plan pattern).

## Context (from discovery, 2026-06-10; spot-checked by plan-review)

### Fortran source (port reference) — all in `port2/fesom2/src/`
| File / lines (≈, see caveat) | Routine / branch | Role | Status |
|---|---|---|---|
| `oce_ale.F90:938-1035` | `init_thickness_ale`, zstar case | init `hnode=(zbar(nz)-zbar(nz+1))·(1+hbar/dd)` for nz=1..`nlevels_nod2D_min`-2, dd=zbar(1)-zbar(nlevels_nod2D_min-1); bottom-touching levels keep nominal; `hnode(nzmax)=bottom_node_thickness`; helem=mean(3 nodes), `dhe=mean(hbar)`, bottom helem=`bottom_elem_thickness`; then `hnode_new=hnode`, `exchange_elem(helem)` | port |
| `oce_ale.F90:1187-1247` | `update_thickness_ale`, zstar case | per-step COMMIT over **myDim+eDim**: for nz=`nlevels_nod2D_min`-2..1 (bottom→top): `hnode=hnode_new; zbar_3d_n(nz)=zbar_3d_n(nz+1)+hnode_new(nz); Z_3d_n(nz)=zbar_3d_n(nz+1)+hnode_new(nz)/2`; helem=mean over nz=1..nlevels(elem)-2; `exchange_elem(helem)` | port (skip ldiag_DVD rescue) |
| `oce_ale.F90:1701-1810` | `update_stiff_mat_ale` | **CUMULATIVE** edge-loop: `SSH_stiff%values(npos) += -dhe(elem)·(gradient_sca×edge_cross_dxdy)·g·dt·alpha·theta` with i/j sign flips; npos found by matching `elnodes(k)` in the row's column list (`nn_pos`) | port |
| `oce_ale.F90:~1931-1951` | `compute_ssh_rhs_ale`, non-linfs tail | `ssh_rhs(n) -= alpha·water_flux(n)·areasvol(nzmin,n)` then `+(1-alpha)·ssh_rhs_old(n)` (no-cavity branch) | port |
| `oce_ale.F90:~2071-2125` | `compute_hbar_ale`, non-linfs bits | `ssh_rhs_old(n) -= water_flux(n)·areasvol(1,n)` + `exchange_nod(ssh_rhs_old)` (exchange only in non-linfs!); then (shared) hbar_old=hbar (myDim+eDim), hbar update (myDim), exchange_nod(hbar); then `dhe(elem)=mean(hbar-hbar_old)` fill (UNCONDITIONAL in Fortran — runs under linfs too, simply unused) | port |
| `oce_ale.F90:2564-2640` | `vert_vel_ale`, zstar branch | over **myDim**: dd1=zbar_3d_n(nzmax,n), nzmax=`nlevels_nod2D_min`-1; dd=(hbar-hbar_old)/(zbar_3d_n(1)-dd1); dddt=dd/dt; for nz=1..nzmax-1: `Wvel(nz) -= (zbar_3d_n(nz)-dd1)·dddt; hnode_new(nz)=hnode(nz)+(zbar_3d_n(nz)-zbar_3d_n(nz+1))·dd`; then `Wvel(1) -= water_flux(n)`; NaN check | port |
| `oce_ale.F90:~2660-2690` | `vert_vel_ale` shared tail | negative-`hnode_new` fatal check (myDim+eDim) + **`exchange_nod(Wvel)` (:2679) + `exchange_nod(hnode_new)` (:2680)** | port |
| `oce_ale.F90:~3723` | timestep gate | `if (.not. linfs) call update_stiff_mat_ale` BEFORE compute_ssh_rhs_ale | port |
| `oce_ale_pressure_bv.F90:~369-401` | `hpressure` block | gated `if (which_ale=='linfs' .or. use_cavity)` — **under zstar (no cavity) Fortran computes NO hpressure at all**. The `-Z_3d_n(1)·ρ(1)·g` term is the LINFS no-cavity term the C already has (`fesom_eos.c:164`); the `0.5·(zbar_3d_n(1)-zbar_3d_n(2))·ρ·g` variant is the CAVITY arm (out of scope). `pressure_force_4_zxxxx_shchepetkin` has ZERO hpressure references — self-contained on `density_m_rho0` | gate C's hpressure to linfs (skip under zstar); nothing new to port here |
| `oce_ale_pressure_bv.F90:1873-1903` | `pressure_force_4_zxxxx` dispatcher | which_pgf='shchepetkin' → | port (gate-only others) |
| `oce_ale_pressure_bv.F90:2104-2339` | `pressure_force_4_zxxxx_shchepetkin` | density-Jacobian PGF on moving levels (~235 lines), reads `density_m_rho0` (computed unconditionally, C `fesom_eos.c:127`) + live `Z_3d_n`/`zbar_3d_n` | port |
| `ice_oce_coupling.F90:476-562` | virtual_salt / relax_salt block | `use_virt_salt=F` → virtual_salt stays 0 (computation + global balancing SKIPPED); `relax_salt=surf_relax_S·(Ssurf-salt(1,n))` + its `-net` balancing UNCHANGED | port the flip |
| `ice_oce_coupling.F90:~625-631` | water-flux assembly | non-virt-salt adds `flux -= thdgr·rhoice/rhowat + thdgrsn·rhosno/rhowat` (ice growth = real volume change) | port |
| `ice_thermo_oce.F90:~352, 617, 646` | **`rsf` → `real_salt_flux` producer (STANDALONE — plan-review blocker fix)** | `real_salt_flux(i)=rsf` is assigned **unconditionally** in the standalone `thermodynamics`; under `use_virt_salt=.false.`, `rsf = fwice·Sice` (:617) then `rsf -= iflice·rhoice·inv_rhowat·Sice` (:646) — **non-zero wherever ice grows/melts; it REPLACES virtual_salt**. C already computes `_rsf` (`fesom_ice_thermo.c:368,402`) but **discards it** (`:616 (void)rsf`), and `fesom_tracer_diff.c:66` hardcodes `real_salt_flux=0` | port (wire the producer + consumer) |
| `ice_thermo_oce.F90:613,645` | thermo fw-flux branches | the `(.not. use_virt_salt)` arms (C marks the linfs arm at `fesom_ice_thermo.c:370`) | port the flip |
| `oce_ale_tracer.F90:535` | explicit-diff surface flux | `flux = virtual_salt + relax_salt + real_salt_flux·is_nonlinfs` | port term |
| `oce_ale_tracer.F90:1517, 1523-1524` | `bc_surface` (implicit TDMA) | T (:1517): `bc_surface = -dt·(heat_flux/vcpw + sval·water_flux·is_nonlinfs)`; S (:1523-1524): `bc_surface = +dt·(virtual_salt + relax_salt + real_salt_flux·is_nonlinfs)` — **positive dt, NO sval·water_flux term in S** (under zstar: effectively `dt·(relax_salt + real_salt_flux)`) | port terms |
| `ice_EVP.F90:~576-727` | non-linfs rhs branch | `use_pice=0` (floatice off) → `p_ice·use_pice=0`; elevation-gradient rhs otherwise = linfs algebra; port the branch literally (cheap) | port |
| `oce_ale_vel_rhs.F90` | floatice pressure term | `use_virt_salt` is imported but **UNUSED** there; the zstar-relevant term is the floatice pressure `×use_pice` (=0 for our config) — byte-identical | audit-confirm only |

NOT writable per-step by anything: `areasvol`, `area`, `gradient_sca`, edge geometry — all
stay static under zstar (Fortran does not touch them; the vert_vel comment explains the
stretching is applied only where `area(nz,n)=area(1,n)`).

### C integration points (target = `fesom2_port/src/`, branch from main@v1.0)
| C file | Current state | Change |
|---|---|---|
| `fesom_ale.c/.h` | linfs-only (`fesom_ale_thickness_linfs`, `fesom_ale_vert_vel_linfs`) | + zstar init/update/vert-vel branches; dispatch on ale mode |
| `fesom_ssh.c` | stiffness built ONCE (`fesom_ssh_stiff_alloc_and_build`); header says "update_stiff_mat_ale NOT called (linfs gate)"; CG + preconditioner (`ssh_solve_preconditioner` port) | + `fesom_update_stiff_mat_ale` (cumulative += into CSR `values`, npos lookup); + non-linfs ssh_rhs term; **preconditioner/diag refresh per step** (check whether C builds `pr` per solve or caches — mirror Fortran `solver.F90` behavior) |
| `fesom_momentum.c:773-836` | hbar update, linfs skips water_flux (comments mark the sites) | + non-linfs `ssh_rhs_old -= wf·areasvol` + its exchange; + `dhe` fill after hbar update |
| `fesom_eos.c` | `hpressure` linfs branch (:160-164, already the `-Z_3d_n·ρ·g` form); `fesom_pressure_force_linfs_fullcell` (:279) | gate hpressure computation to linfs (Fortran skips it under zstar); + `fesom_pressure_force_zxxxx_shchepetkin`; dispatch on ale mode |
| `fesom_mesh.c/.h` | `hnode/hnode_new/helem/hbar/hbar_old` exist; `zbar_3d_n` static (:563-631, "constant in time"); **no per-node `Z_3d_n`, no `dhe`, no `bottom_node/elem_thickness`** | + `Z_3d_n[N·nl]` (init = mesh->Z), + `dhe[elem2D]`, + `bottom_node_thickness[N]`/`bottom_elem_thickness[E]` (nominal, full cells); zbar_3d_n becomes live under zstar |
| `fesom_step.c` | linfs call order | + `update_stiff_mat_ale` call before ssh_rhs when zstar; vert_vel/update_thickness dispatch |
| `fesom_ice_coupling.c` | linfs virtual-salt path (mirrors ice_oce_coupling.F90) | + use_virt_salt=F arms: virtual_salt=0, thdgr/thdgrsn water-flux terms, relax_salt unchanged |
| `fesom_ice_thermo.c` | computes `_rsf` (:368,402) then **discards** it (:616); virtual-salt branch marked at :370 | + the non-virt-salt arms; **wire `rsf` into the new `real_salt_flux[N]`** |
| `fesom_ice_evp.c:201-291` | linfs rhs branch (comments mark it) | + non-linfs branch with `use_pice=0` |
| `fesom_tracer_diff.c` | explicit flux + implicit `bc_surface`; `:66` hardcodes `real_salt_flux=0.0` | + `is_nonlinfs` terms with the CORRECTED S formula; consume the live `real_salt_flux` |
| `fesom_main.c` | env options block | + `FESOM_ALE` parsing → ale mode + use_virt_salt/is_nonlinfs |
| geometry consumers | see Phase Z7 (audit greps the geometry ARRAYS, not `linfs` comments) | re-point static-geometry assumptions to live arrays exactly where Fortran reads them |

### Out of scope (faithful = port only what the zstar reference exercises)
- **`zlevel` + local-zstar fallback** (`lzstar_lev`, `min_hnode`, `CFL_z` guards, distribute/
  refill) — zlevel-only machinery, no reference run.
- **Cavity bodies** — guards (`nzmin==1`, `ulevels>1 cycle`) ported; bodies unreachable.
- **`restart_thickness_ale`** — C has no restart system.
- **`ldiag_DVD` rescue of hnode_old** in update_thickness — diagnostic off.
- **`use_ssh_se_subcycl` branch** — CG path only.
- **floatice body** beyond the literal `p_ice·use_pice` (=0) factor (EVP and vel_rhs);
  icebergs/landice/age/CFC tracer arms; `use_cavity_fw2press`.
- **`pressure_force_4_zxxxx_cubicspline`/`easypgf`** — gate-only (abort if selected).
- **`gen_modules_diag` `is_nonlinfs` arms** (:2560,2562,3099,3101) + **io_meandata
  virtual-salt output gate** — diagnostics off in the reference; C stream IO unchanged.

### Key Fortran invariants (must hold in C)
1. **The stiffness update is CUMULATIVE**: `values += fy·factor` each step with
   `dhe = mean(hbar−hbar_old)` from the PREVIOUS step's `compute_hbar_ale`; the base matrix
   (init over unperturbed depth) is never rebuilt. Cold start: init_thickness sets
   `dhe=mean(hbar)=0` → step-1 update is a no-op. C must keep the same CSR object and
   increment it, replicating the i/j sign flips and the npos column lookup.
2. **Step order (non-subcycl)**: …momentum… → `update_stiff_mat_ale` → `compute_ssh_rhs_ale`
   (−α·wf·areasvol + (1−α)·ssh_rhs_old) → `solve_ssh_ale` → `update_vel` → `compute_hbar_ale`
   (ssh_rhs_old −wf·areasvol + exchange; hbar; dhe for the NEXT step) → `eta_n=α·hbar+(1−α)·hbar_old`
   (formula unchanged) → … → `vert_vel_ale` (zstar distribute + wf BC + exchanges) →
   tracers (generic ALE reconstruction — **zero `which_ale` branches in oce_ale_tracer.F90**,
   only the `is_nonlinfs` surface-BC terms) → `update_thickness_ale` (commit + geometry + helem).
3. **Bottom protection bounds**: thickness scaling/Wvel distribution run only to
   `nlevels_nod2D_min(n)−2` (init/commit geometry) / `nlevels_nod2D_min(n)−1` (vert_vel `nzmax`);
   levels intersecting bathymetry keep nominal spacing; bottom layer = `bottom_node_thickness`.
   Element helem updates stop at `nlevels(elem)−2`; element bottom keeps `bottom_elem_thickness`.
4. **Halo discipline** (`feedback_write_loops_halo`): `vert_vel_ale` writes `Wvel`/`hnode_new`
   over **myDim** then `exchange_nod(Wvel)` (:2679) + `exchange_nod(hnode_new)` (:2680);
   `update_thickness_ale` commits over **myDim+eDim** reading those halos. `compute_hbar_ale`
   non-linfs adds an `exchange_nod(ssh_rhs_old)` that linfs does NOT have. helem gets
   `exchange_elem` after both init and update. Replicate every bound and exchange exactly.
5. **Wvel correction is the vertically-integrated form**: `−(zbar_3d_n(nz)−dd1)·dddt`
   (NOT per-layer h·dd/dt) — it is the bottom-to-top sum of layer dh/dt; port the formula as
   written, including using the CURRENT (pre-update) `zbar_3d_n`.
6. **Freshwater/salt semantics flip**: `virtual_salt=0` (its computation + global balancing
   loop SKIPPED); `relax_salt` (SSS restoring, incl. its `−net` balancing) UNCHANGED;
   water_flux gains the ice growth terms (−thdgr·ρice/ρwat −thdgrsn·ρsno/ρwat); water_flux now
   enters `ssh_rhs`, `ssh_rhs_old`, and the surface `Wvel` BC. **`real_salt_flux` is a LIVE
   producer**: standalone `therm_ice` assigns it unconditionally (`ice_thermo_oce.F90:352`)
   and under zstar it carries the brine-rejection/melt salt flux that replaces virtual_salt —
   the C must stop discarding `rsf` and feed it to the S surface BC
   (`+dt·(virtual_salt + relax_salt + real_salt_flux·is_nonlinfs)`, no sval·wf term in S).
7. **PGF**: linfs ignores `which_pgf` (full-cell branch); zstar honors it → shchepetkin
   (self-contained on `density_m_rho0` + live geometry; NO hpressure involved — and under
   zstar the hpressure block itself is skipped, matching Fortran's `linfs .or. use_cavity` gate).
8. **EVP**: non-linfs branch = linfs algebra when `use_pice=0`; port the branch literally
   (p_ice computed, ×0) so the code matches the Fortran flow.
9. **No new AB terms** in the zstar path (α/θ are the θ-method weights, not AB) — but apply
   the AB-epsilon lesson (o_PARAM eps=0.1, `project_dt1800_state`) to anything AB-weighted
   encountered during the port.

## Development Approach
- Faithful port discipline: one Fortran routine per phase, dependency order, literal diff, no refactors.
- **Validation = dual-instrumentation vs Fortran** (dump-diff, NOT unit tests — the project's
  proven idiom; matches KPP/sea-ice/EVP/MFCT). Each phase ends with a Validate gate before the
  next phase starts.
- **linfs byte-identical to v1.0** is a discrete Validate item in EVERY phase that compiles
  changes into shared code.
- Update this plan (`[x]`, ➕, ⚠️) as you go. Unique OUT_DIR per SLURM job
  (`feedback_unique_outdir_per_job`). New node arrays sized `[N·nl]` stride nl
  (`feedback_tracer_stride_nl`), `myDim+eDim` (`feedback_array_size_vs_reader_loop`).

## Validation Strategy
- **Dual dump gate `FESOM_ALE_DUMP_DIR`** (both codes; Fortran gate uncommitted in the Intel
  tree like the KPP/EVP dumps): per-step dumps of `water_flux`, `virtual_salt`, `relax_salt`,
  `real_salt_flux`, `hbar`, `dhe`, `ssh_rhs`, `d_eta`, `Wvel`, `hnode_new`, `helem`,
  `zbar_3d_n`, `pgf_x/y` at steps 1–3, 16r, dt=1800, same partition → dump line i = local
  node i. Step-1 zstar from cold start: hbar=0 → thickness updates are near-no-ops; the chain
  becomes decisively non-trivial at steps 2–3 — dump those too.
- **Input-cross-check methodology** (KPP K3/K5): output diffs must be explained by input diffs
  (PHC rank-noise floor, `feedback_phc_rank_dependent`); 0 unexplained disagreements is the bar.
- **Byte-gates**: `FESOM_ALE=linfs` (default) runs bit-identical to v1.0 (16r, 200 steps,
  dt=500, snapshots + run.log; `scripts/kpp_byteident_check.py` reused).
- **Cross-rank**: 1/8/16r zstar short runs, exchange-and-compare probe on the new
  halo-written fields (Wvel, hnode_new, helem, ssh_rhs_old, hbar).
- **Climate**: 2yr C(zstar+KPP) 864r vs NEW Fortran `work_zstar_kpp` (clim-compare scripts
  reused). Acceptance bar: same order as KPP-vs-KPP (SST/SSS RMS ~0.005/0.002, biases ~0,
  non-drifting). The zstar-vs-linfs contrast run (C zstar vs Fortran linfs) should show a
  clearly LARGER coordinate signal — the comparison must resolve coordinates (KPP K9 pattern).
- **Stability**: 5-day dt=1800 smoke each phase-completion; full 5yr lands in the combined
  zstar+TKE matrix (see TKE plan / Post-Completion).

## Implementation Steps

### Phase Z0: Scaffolding — switch, mesh arrays, dump harness, reference runs (NO behavior change)
**Files:** Modify `src/fesom_main.c`, `src/fesom_mesh.{c,h}`, `src/fesom_ale.{c,h}`, `src/fesom_step.c`; Create `scripts/ale_dump_diff.py`, `jobs/job_zstar_*`; Fortran `oce_ale.F90` (dump gate, uncommitted)
- [x] Branch from main@v1.0. Add `FESOM_ALE=linfs|zstar` env parsing (default linfs) → step-ctx
      flags `ale_zstar`, `use_virt_salt=!ale_zstar`, `is_nonlinfs=ale_zstar?1.0:0.0`.
      ➕ 2026-06-10: branch `zstar` lives in git worktree `fesom2_port_zstar/` (the main
      checkout keeps jax-mesh-export WIP undisturbed); v1.0 reference binary built in
      detached worktree `fesom2_port_v10/`. Flags = globals in fesom_ale.{c,h} (mirror of
      the Fortran o_PARAM/g_config module globals) + copies in step-ctx; zstar-only
      startup print (linfs run.log stays byte-identical); unknown FESOM_ALE values abort.
- [x] Add mesh arrays, **initialized to linfs-identical values**: per-node `Z_3d_n[N·nl]`
      (=mesh->Z per level), `dhe[elem2D]`(=0), `bottom_node_thickness[N]`/`bottom_elem_thickness[E]`
      (nominal zbar spacing at the column bottom, full cells), forcing `real_salt_flux[N]`(=0).
      No consumer change yet. ➕ Z_3d_n/bottom_* built in fesom_mesh_compute_metrics (next
      to compute_zbar_3d_n); dhe in fesom_mesh_alloc_state; real_salt_flux in fesom_forcing.
- [x] Audit-confirm `oce_ale_vel_rhs.F90`: use_virt_salt imported-unused; floatice pressure
      term ×use_pice=0 → byte-identical for our config (record in plan).
      ➕ CONFIRMED: `use_virt_salt` only at :44 (import), zero body uses. `use_pice=0`
      (:93) flips to 1 only if `use_floatice .and. .not. linfs` (:94); use_floatice=.false.
      in the reference → the p_ice=0.0 arm (:177) runs under BOTH modes.
- [x] Launch the **Fortran reference runs NOW** (queue time): (a) `work_zstar_dump` 16r dt=1800
      3 steps with the ALE dump gate; (b) `work_zstar_kpp` 2yr dt=1800 (copy of work_linfs_2yr,
      flip `which_ALE='zstar'` only, same Intel binary). Unique OUT_DIRs. Archive namelists →
      `docs/zstar_reference_namelists/` + PROVENANCE.md; re-confirm the config table above.
      ➕ Launched 2026-06-10: jobs 25492898 (work_zstar_kpp) + 25492899 (work_linfs_2yr_b —
      RERUN of the linfs reference: the original /scratch/.../fortran_2yr_dt1800 was PURGED;
      rerun with the SAME frozen binary as the zstar job → perfectly-paired one-knob pair +
      restores the Z9 contrast baseline) + 25493084 (work_zstar_dump, COMPLETED exit 0:
      576 files = 12 tags × 3 steps × 16 ranks; step-1 virtual_salt=0 confirms the
      use_virt_salt=F derivation in the reference; hbar evolves s1→s3; no NaN). Outputs
      under /work/ab0995/a270088/port/zstar/ (purge-safe). Namelists archived +
      PROVENANCE.md (frozen-binary sha, one-knob diff proof, dump-gate site table).
      Config table re-confirmed against the actual run namelists (which_pgf absent ✓,
      i_vert_diff=T ✓, use_partial_cell=F ✓, use_floatice=F ✓, use_ssh_se_subcycl=F ✓).
- [x] Add `FESOM_ALE_DUMP_DIR` dump plumbing in C (writer skeleton; fields wired per phase)
      + the Fortran dump gate (uncommitted, Intel tree) + `scripts/ale_dump_diff.py`.
      ➕ C: src/fesom_ale_dump.{c,h} with ALL 12 tags wired at the 6 mirrored driver sites
      (not per-phase — every dumped field already exists in C); Fortran: `ale_dump_mod` at
      the top of oce_ale.F90 + 6 one-line call sites (uncommitted, KPP-gate convention);
      identical filename/format; FESOM_ALE_DUMP_STEPS default 3. ⚠️ The C vertvel dump
      sits BEFORE the GM bolus add (13a) which modifies dyn->w in place — matching the
      Fortran site (before solve_tracers_ale).
- [x] **Validate:** C with FESOM_ALE unset/linfs **byte-identical to v1.0** (16r/200/dt=500,
      snapshots + run.log). Fortran zstar dump run sane (hbar evolves, no NaN).
      ➕ Fortran dump sanity PASSED (above). Byte-gate job 25493375 PASSED
      (jobs/job_zstar_z0_byteident: v1.0 binary vs zstar-branch binary, both FESOM_ALE
      unset): all 9 snapshots BYTE-IDENTICAL (worst |Δ|=0); run.log content identical —
      the only diff lines are the out_dir PATH embedded in two log strings + srun rank
      interleaving order (sorted diff shows path-only lines). Z0 COMPLETE 2026-06-10.

### Phase Z1: Thickness machinery — init_thickness_ale + update_thickness_ale (zstar cases)
**Files:** Modify `src/fesom_ale.{c,h}`, `src/fesom_ic.c` (init call), `src/fesom_step.c` (commit call)
- [x] Port init_thickness_ale zstar: proportional hnode over nz=1..`nlevels_nod2D_min`-2
      (myDim+eDim loop), nominal bottom-touching levels, `hnode(nzmax)=bottom_node_thickness`,
      helem mean + `dhe=mean(hbar)` + `bottom_elem_thickness`, `hnode_new=hnode`,
      `exchange_elem(helem)`. (Cold start hbar=0 → values equal linfs init — by construction.)
      ➕ `fesom_ale_init_thickness_zstar` (fesom_ale.c) incl. the mode-shared pre-block
      (ssh_rhs_old/eta_n — REVERSED α weights, zeros at cold start); dispatched from
      fesom_main (zstar replaces fesom_ic_thickness). Cavity elem-arm = fatal guard
      (zbar_e_srf not ported, unreachable).
- [x] Port update_thickness_ale zstar: bottom→top commit loop (myDim+eDim) writing hnode,
      zbar_3d_n, Z_3d_n; helem mean to `nlevels(elem)-2`; `exchange_elem(helem)`. Wire as the
      step-end commit when zstar (linfs keeps `fesom_ale_thickness_linfs`).
      ⚠️ Note: until Z5 produces `hnode_new`, the commit is a NO-OP (hnode_new=hnode from
      init) — its bottom→top geometry writes are first meaningfully validated at Z5/Z8.
      ➕ `fesom_ale_update_thickness_zstar`; step-12 linfs hnode_new=hnode copy SKIPPED
      under zstar (Fortran has no per-step copy); step-14 dispatch (zstar path does NOT
      exchange hnode — Fortran doesn't; commit covers myDim+eDim from the exchanged
      hnode_new halo).
- [x] **Validate:** (a) zstar step-0 state (hnode/helem/zbar_3d_n/Z_3d_n) bit-equal to linfs at
      cold start (hbar=0 degenerate case); (b) linfs byte-gate vs v1.0; (c) dump-diff hnode/helem
      vs Fortran zstar init (16r) — bit-identical expected (pure geometry).
      ➕ Job 25493641 (jobs/job_zstar_z1_smoke, 16r dt=1800 3 steps zstar + linfs, dumps on):
      BOTH runs exit 0; ale_dump_diff zstar-vs-linfs over ALL 12 tags × 3 steps:
      **worst max|diff| = 0.000e+00 — the entire zstar run is BIT-IDENTICAL to linfs**
      (expected: no behavioral arm is live until Z2; even helem's mean-of-3 is exact
      because the CORE2 zbar values are round numbers). This subsumes (a); (b) covered by
      the Z0 byte-gate (25493375, snapshots worst |Δ|=0). (c) replaced by the stronger
      bit-identity: the Fortran fdump thickness tags reflect step-1..3 EVOLVED geometry
      (full chain live in Fortran from step 1), so the cross-code geometry compare lands
      at Z5/Z8 as the ⚠️ above says; the init-state equivalence is proven exactly here.
      Z1 COMPLETE 2026-06-10.

### Phase Z2: Freshwater/salt plumbing flip (use_virt_salt=.false.) + real_salt_flux producer + ice
**Files:** Modify `src/fesom_ice_coupling.c`, `src/fesom_ice_thermo.c`, `src/fesom_ice_evp.c`, `src/fesom_tracer_diff.c`, `src/fesom_sss_runoff.c` (if restoring touched)
*(Moved before the SSH/vert-vel phases — their dump validations consume `water_flux`, which
must already have the zstar composition; plan-review finding 5.)*
- [ ] Flip `fesom_ice_coupling.c` to the non-virt-salt arms: virtual_salt=0 (skip computation
      + balancing); add `flux -= thdgr·ρice/ρwat + thdgrsn·ρsno/ρwat`; relax_salt computation +
      `−net` balancing UNCHANGED. Map line-for-line against ice_oce_coupling.F90:476-660.
- [ ] **Wire the `real_salt_flux` producer (plan-review BLOCKER fix):** stop discarding `rsf`
      at `fesom_ice_thermo.c:616` — assign `real_salt_flux[n] = rsf` exactly as
      `ice_thermo_oce.F90:352` does (unconditional; the value itself is use_virt_salt-branched
      at :617/:646). Flip the `:370`-area fw-flux branches to the `(.not. use_virt_salt)` arms.
- [ ] Port the EVP non-linfs branch (`use_pice=0`, p_ice×0) in `fesom_ice_evp.c`.
- [ ] Tracer surface BCs with the CORRECTED formulas: explicit path (tracer:535
      `flux = virtual_salt + relax_salt + real_salt_flux·is_nonlinfs`); implicit `bc_surface`
      T (:1517) `−dt·(heat_flux/vcpw + sval·water_flux·is_nonlinfs)`; S (:1523-1524)
      `+dt·(virtual_salt + relax_salt + real_salt_flux·is_nonlinfs)` — replace the
      `fesom_tracer_diff.c:66` hardcoded zero with the live array.
- [ ] ⚠️ Audit what C's restoring/runoff currently folds into water_flux vs Fortran's
      relax_salt/water_flux split (`feedback_runoff_fold_when_ice`, ref_sss_local b8e1d36) —
      the zstar arms must match Fortran's composition exactly, not C's linfs shortcut if any.
- [ ] **Validate:** step-1 dump-diff vs Fortran zstar of `water_flux`, `virtual_salt`(=0),
      `relax_salt`, `real_salt_flux`, ice state (a_ice/m_ice/u_ice) — step 1 is clean
      (identical IC; ocean state hasn't diverged); later steps re-checked at Z8 once the full
      chain exists. Global salt-content accounting sanity (virtual_salt path off, real path on).
      linfs byte-gate (all flips gated on ale mode).

### Phase Z3: SSH plumbing — ssh_rhs + hbar water-flux terms + dhe
**Files:** Modify `src/fesom_ssh.c`, `src/fesom_momentum.c`
- [ ] Port the non-linfs `compute_ssh_rhs_ale` tail: `ssh_rhs(n) -= alpha·water_flux·areasvol(1,n)`
      (+ the existing `(1−alpha)·ssh_rhs_old` term in the same branch shape as Fortran).
- [ ] Port the non-linfs `compute_hbar_ale` bits: `ssh_rhs_old(n) -= water_flux·areasvol(1,n)`
      + **`exchange_nod(ssh_rhs_old)`** (non-linfs only!); `dhe(elem)=mean(hbar−hbar_old)` fill
      after the hbar exchange — the fill is UNCONDITIONAL in Fortran (runs under linfs too,
      unused there): mirror that, byte-gate stays green since dhe is a new array.
- [ ] **Validate:** dump-diff `ssh_rhs`, `hbar`, `dhe` steps 1–3 vs Fortran zstar (16r) —
      water_flux already zstar-composed (Z2), so diffs must be at the input-noise floor;
      explained-by-input bar; linfs byte-gate.

### Phase Z4: Per-step stiffness update + solver
**Files:** Modify `src/fesom_ssh.{c,h}`, `src/fesom_step.c`
- [ ] Port `update_stiff_mat_ale`: edge loop, `row>myDim cycle`, npos lookup of `elnodes(k)`
      in the CSR row, `fy=-dhe·(gradient_sca×edge_cross_dxdy)` with `if(i==2)fy=-fy; if(j==2)fy=-fy`,
      `values(npos) += fy·g·dt·alpha·theta`. Called before ssh_rhs when zstar (gate :3723).
- [ ] Mirror the solver's preconditioner behavior: if the C Jacobi `pr` is cached from init,
      refresh the diagonal after each update (Fortran `solver.F90` reads the live matrix).
      ⚠️ Check first what C does; mirror, don't invent.
- [ ] Add (debug-gated) the Fortran's commented row-sum sanity: `sum(row values)/area·dt ≈ 1`.
- [ ] **Validate:** dump-diff `d_eta` steps 1–3 vs Fortran zstar (step 1: dhe=0 → must equal
      the no-update solve bit-for-bit; steps 2–3 exercise the update); CG iteration counts
      comparable to Fortran log; linfs byte-gate (matrix path untouched when linfs).

### Phase Z5: Vertical velocity + hnode_new (vert_vel zstar branch)
**Files:** Modify `src/fesom_ale.c`, `src/fesom_step.c`
- [ ] Port the zstar branch: dd1/dd/dddt, Wvel integrated correction
      `−(zbar_3d_n(nz)−dd1)·dddt` + `hnode_new=hnode+(zbar_3d_n(nz)−zbar_3d_n(nz+1))·dd`
      over nz=1..`nlevels_nod2D_min`−2 (myDim), surface `Wvel(1) −= water_flux`, NaN check.
- [ ] Port the shared tail: negative-hnode_new fatal check (myDim+eDim) +
      **`exchange_nod(Wvel)` AND `exchange_nod(hnode_new)`** (:2679-2680; the latter feeds
      Z1's commit — write-loops-halo lesson).
- [ ] **Validate:** dump-diff `Wvel`, `hnode_new` steps 1–3 vs Fortran zstar; exchange-and-
      compare probe on Wvel/hnode_new; the Z1 commit now does real work — re-run the Z1
      geometry dump at step 2–3 (zbar_3d_n/Z_3d_n/helem evolve); linfs byte-gate.

### Phase Z6: Pressure — zxxxx shchepetkin PGF (+ hpressure gating)
**Files:** Modify `src/fesom_eos.{c,h}`, `src/fesom_step.c`
- [ ] Gate the C `hpressure` computation to linfs, mirroring Fortran's
      `if (which_ale=='linfs' .or. use_cavity)` — under zstar NO hpressure is computed
      (plan-review finding 2: there is no "non-linfs hpressure term"; nothing new to port).
- [ ] Port `pressure_force_4_zxxxx_shchepetkin` (:2104-2339) as
      `fesom_pressure_force_zxxxx_shchepetkin` (reads `density_m_rho0` + live
      `Z_3d_n`/`zbar_3d_n`); dispatch linfs→fullcell / zstar→shchepetkin; gate-only aborts
      for cubicspline/easypgf.
- [ ] **Validate:** dump-diff `pgf_x/pgf_y` steps 1–3 vs Fortran zstar (apples-to-apples:
      zstar-vs-zstar, NOT vs C linfs — different algorithm); linfs byte-gate (hpressure gating
      must be ale-conditional).

### Phase Z7: Geometry-consumer + halo audit
**Files:** Audit-driven small edits across `src/fesom_gm.c`, `src/fesom_kpp.c`, `src/fesom_pp.c`, `src/fesom_tracer_adv.c`, mo_convect site, sw-pene site, wsplit/CFL site
- [ ] **Audit method (plan-review finding 4): grep the geometry arrays themselves** —
      `mesh->Z\b`, `mesh->zbar\b`, `zbar_3d_n`, `Z_3d_n`, `hnode`, `helem` across `src/*.c` —
      and cross-reference each read against what the Fortran reads there. The `linfs`-comment
      grep is the secondary sweep (catches the annotated shortcuts: gm:486/515/812, eos:67/319,
      mesh:564, main:456/484/503, …).
- [ ] Explicitly covered consumers (from Fortran cross-check): **mo_convect**
      (`oce_mo_conv.F90:64` reads `zbar_3d_n` for the convection depth test); **vertical
      tracer advection** (`oce_adv_tra_ver.F90:417-422` reads `Z_3d_n`/`zbar_3d_n` directly,
      AND `:137-138` deliberately builds a LOCAL `zbar_n`/`Z_n` from **`hnode_new`**, not
      `zbar_3d_n` — a dual-geometry subtlety invisible under linfs; port each site's exact
      source); **PP `dz_inv`** (`oce_ale_mixing_pp.F90:49` uses `Z_3d_n`); KPP `zbar_3d_n`
      reads (bldepth); sw-pene decay depths; wsplit `CFL_z` from live hnode; GM slopes/bolus;
      MFCT/FCT `hnode_new` divisions.
- [ ] Halo audit of every NEW write loop (`feedback_write_loops_halo`): exchange-and-compare
      probe over Wvel/hnode_new/helem/ssh_rhs_old/hbar at 8r.
- [ ] **Validate:** 30-day zstar runs at 1/8/16r — cross-rank consistency at the PHC-noise
      scale (GM/sea-ice precedent); linfs byte-gate after the audit edits (they must be
      zstar-gated or value-identical under linfs).

### Phase Z8: Integration validation (short)
- [ ] Steps 1–3 FULL dump chain green (every Z1–Z6 field incl. the Z2 forcing fields at
      steps 2–3 now that the full chain exists; 0 unexplained); day-5 dt=1800 16r stable
      (max|uv|, SSH range, SST/SSS sane vs Fortran zstar day-5).
- [ ] ⚠️ Watch the known blowup locations (Aleutian trench el 194724) and the SSH-checkerboard
      indicator (1-ring hbar spread) — zstar changes the barotropic operator.
- [ ] **Validate:** dt=1800 30-day zstar run clean; linfs byte-gate final re-check.

### Phase Z9: 2yr climate validation vs work_zstar_kpp
**Files:** `scripts/zstar_climate_compare.py` (clone of kpp_climate_compare.py), `jobs/job_zstar_2yr`
- [ ] Run C zstar+KPP 2yr dt=1800 864r (unique OUT_DIR).
- [ ] Compare vs Fortran `work_zstar_kpp` output: global/regional SST/SSS RMS + bias maps +
      seasonal cycle. Bar: KPP-class agreement (RMS ~0.005/0.002, biases ~0, non-drifting).
- [ ] Contrast check: C(zstar) vs Fortran(linfs) must show a clearly larger coordinate signal
      than C(zstar) vs Fortran(zstar) — the comparison resolves the coordinate (K9 pattern).
- [ ] **Validate:** acceptance numbers recorded in `docs/validation_zstar_2yr/`.

### Phase Z10: Acceptance + handoff
- [ ] All Overview requirements met: FESOM_ALE switch works; linfs byte-identical to v1.0;
      every ported routine dump-validated; 2yr climate matches the zstar reference.
- [ ] Decision point (user): keep linfs default vs flip to zstar (KPP precedent: flip after
      validation). Tag suggestion: `zstar-validated-<date>`.
- [ ] Update docs/NEXT_SESSION_HANDOFF.md + auto-memory (`project_zstar_port_state`);
      commit per-phase; move this plan to `docs/plans/completed/`.

## Technical Details
- **New state**: mesh `Z_3d_n[N·nl]`, `dhe[E]`, `bottom_node_thickness[N]`,
  `bottom_elem_thickness[E]`; forcing `real_salt_flux[N]` (**a live producer under zstar** —
  wired from `rsf` in ice thermo, consumed by the S surface BC); step-ctx `ale_zstar`,
  `use_virt_salt`, `is_nonlinfs`. All node arrays `myDim+eDim`, stride nl.
- **zstar thickness math**: init `hnode=(zbar(nz)−zbar(nz+1))·(1+hbar/dd)`,
  `dd=zbar(1)−zbar(nlevels_nod2D_min−1)`; per-step `hnode_new=hnode+(zbar_3d_n(nz)−zbar_3d_n(nz+1))·dd`
  with `dd=(hbar−hbar_old)/(zbar_3d_n(1)−zbar_3d_n(nlevels_nod2D_min−1))`; Wvel correction
  `−(zbar_3d_n(nz)−dd1)·(dd/dt)` = vertically integrated dh/dt (bottom→top sum).
- **Stiffness update semantics**: cumulative `values(npos) += −dhe·(g·dt·α·θ)·(∇sca×edge_cross)`;
  the base matrix encodes unperturbed depth; Σdhe over steps ≡ hbar−hbar(init).
- **Dispatch**: single `ale_zstar` flag checked at each Fortran `which_ale` branch site —
  same shape as the Fortran `trim(which_ale)=='...'` chains (no function-pointer cleverness).
- **C step-order insertion points** mirror invariant 2 exactly; no reordering.

## Post-Completion
*Informational — external/manual:*

**Run matrix (approved 2026-06-10, shared with the TKE plan):**
- After BOTH features validate per-feature: Fortran `work_zstar_tke` 2yr (both knobs) + C
  zstar+TKE 2yr compare; then 5yr C zstar+TKE stability run (autumn uv-transient watch —
  linfs margin was 0.35 m/s; zstar changes the barotropic operator, re-measure the margin).

**Manual verification:**
- SSH/steric consistency: zstar global-mean SSH now reflects net freshwater volume input —
  compare its drift against Fortran zstar (linfs hides this).
- Conservation: global salt content drift under zstar (real-flux path: water_flux volume +
  real_salt_flux salt) ≈ Fortran's.

**Decisions deferred to the user:**
- Default coordinate flip (linfs → zstar) after validation.
- Whether to ever port zlevel + local-zstar fallback (only if a reference config needs it).
- Restart support for time-varying thickness state (no restart system in C today).
