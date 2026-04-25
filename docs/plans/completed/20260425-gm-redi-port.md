# GM/Redi port plan (FESOM2 → C)

Goal: port Gent–McWilliams (Ferrari et al. 2010 form) and Redi
isoneutral diffusion. Default-on, env-gated for easy disable.

> **v2 — incorporates plan-review feedback (2026-04-25)**:
>   - G2 split into G2a (`MLD1_ind` extension to existing pressure_bv,
>     unconditional) + G2b (sigma_xy + neutral_slope, gated).
>   - G6 split into G6a (`fer_w` inside the existing vert_vel_linfs
>     edge loop, validated against bit-identity when GM is off) + G6b
>     (bolus add/sub around tracer block).
>   - G7 was misfiled. Split into G7a (`diff_ver_part_redi_expl` —
>     vertical-explicit off-diagonal Redi) + G7b (`diff_part_hor_redi`
>     — horizontal Redi). Both live in tracer **diffusion**, not
>     advection.
>   - Mesh prerequisites added: `nlevels_nod2D_min`, `ulevels_nod2D_max`,
>     `mesh_resolution`, per-node `zbar_3d_n`, `tr_xy/tr_z/tr_xynodes`.
>   - LoC estimate revised to ~1400-1500 (was ~1050).
>   - Risk ranking corrected: G6a > G7b > G7a > G3 > G2b > rest.

## Scope

Two coupled physics modules that share the same prerequisites:

- **Redi**: along-isopycnal mixing — adds an off-diagonal tensor to
  the tracer diffusion. Acts in **three places** (this is the easy
  thing to miss; see G5/G7a/G7b):
    1. **G5**: vertical-implicit augmentation of K33 inside
       `diff_ver_part_impl_ale` (already plumbed in
       `fesom_tracer_diff.c:158` with hard zeros — swap them for the
       real expressions).
    2. **G7a**: vertical-explicit off-diagonal terms via
       `diff_ver_part_redi_expl` (writes into `del_ttf`).
    3. **G7b**: horizontal Redi flux via `diff_part_hor_redi`
       (writes into `del_ttf`).

- **GM** (Ferrari et al. 2010): bolus / eddy-induced transport
  velocity. Acts in tracer advection by adding `fer_uv`/`fer_w` to
  the velocity field, applied only during tracer transport then
  subtracted back.

Both use the same neutral-slope / tapering machinery and a shared
diffusivity scaling pipeline (`init_Redi_GM`).

### Source files (in Fortran)

| File | LoC | Subroutines we need |
|------|-----|---------------------|
| `oce_fer_gm.F90` | 552 | `init_Redi_GM` (215–510), `fer_solve_Gamma` (40–163), `fer_gamma2vel` (168–210); skip `init_RediGM_GINsea_mask` (default `scaling_GINsea=.false.`) |
| `oce_ale_pressure_bv.F90` | parts | `compute_sigma_xy` (2851), `compute_neutral_slope` (2949), `MLD1_ind` (~420-449 — inline inside `compute_pressure_bv`) |
| `oce_ale_tracer.F90` | parts | bolus add/sub around tracer block (~199–211, 284–292); Redi tensor in vertical-implicit TDMA (761–838); `diff_ver_part_redi_expl` (~1086–1169); `diff_part_hor_redi` (~1173+, 5 partial-cell regions A–E) |
| `oce_ale.F90` | parts | call `init_Redi_GM` + `fer_solve_Gamma` + `fer_gamma2vel` (3815–3827); **`fer_w` accumulation INSIDE the existing `vert_vel_ale` edge loop** (~2735–2860) — invasive modification, not a separate kernel |

### Skipped sub-features (default off in our config)

Per `namelist.oce` defaults:

- `scaling_Rossby = .false.` — no Rossby-radius cutoff
- `scaling_FESOM14 = .false.` — no NH boundary-layer mask
- `scaling_GINsea = .false.` — no GIN-sea upscale; skip `init_RediGM_GINsea_mask` and `point_in_polygon` / `smooth_nod` machinery
- `scaling_Ferreira = .false.` — no Ferreira vertical scaling
- `K_GM_Ktaper = .false.` — no GM-coefficient tapering (only Redi-tapering active)
- `K_GM_rampmax/rampmin = -1` — no resolution ramp
- `__oifs / __ifsinterface / SPP / __icepack` — not configured

Active sub-features (per default namelist):

- `Fer_GM = .true.`, `Redi = .true.`
- `K_GM_max = 1000`, `K_GM_min = 2`, `K_GM_cmin = 0.1`, `K_GM_cm = 3`
- `Redi_Kmax = 0` (auto-sync to `K_GM_max`), `Redi_Kmin = 100`
- `K_GM_resscalorder = 2` (area-based scaling), `scaling_resolution = .true.`
- `scaling_GMzexp = .true.`, `GMzexp_zref = 500`, `GMzexp_smin = 0.6`
- `Redi_Ktaper = .true.`, `K_GM_bvref = 1` (bottom-of-MLD reference — depends on `MLD1_ind`)

## Prerequisites in current C port

Plan-review revealed several gaps from v1. Updated table:

| Prerequisite | Status | Action |
|---|---|---|
| `mesh->bvfreq` | ✓ already in tree (`fesom_aux`) | none |
| `mesh->mesh_resolution` (Voronoi-radius per node) | ✗ missing | Phase G0 — port from `oce_mesh.F90`, compute once at setup |
| `mesh->nlevels_nod2D_min[n]` (min over surrounding elements) | ✗ missing | Phase G0 — derived from `nlevels` + `nod_in_elem2D`, compute once at setup |
| `mesh->ulevels_nod2D_max[n]` (max over surrounding elements) | ✗ missing | Phase G0 — analogous to above |
| `mesh->zbar_3d_n` (per-node 3D bottom depth, sized `nl × nod2D`) | ✗ missing | Phase G0 — derived per step from `hnode_new` (or one-time for linfs); used by `scaling_GMzexp` |
| `MLD1_ind[n]` (mixed-layer-depth index per node) | ✗ — explicitly deferred (`fesom_eos.h:35`) | Phase G2a — extend `fesom_pressure_bv` (~20 LoC, lives there in Fortran too); unconditional (free) |
| `tracers->work->tr_xy` / `tr_z` / `tr_xynodes` (per-tracer gradients on edges/nodes) | ✗ missing | Phase G7b — needed only for horizontal Redi; produce inside the per-tracer Redi pass |
| `compute_sigma_xy` → `aux->sigma_xy[2 * nl * nod2D]` | ✗ missing | Phase G2b |
| `compute_neutral_slope` → `aux->neutral_slope[3 * (nl-1) * nod2D]`, `aux->slope_tapered[3 * (nl-1) * nod2D]`, `aux->fer_tapfac[nl * nod2D]` | ✗ missing | Phase G2b |
| Tracer-diff `Ki` and `slope_tapered` plumbed | partial — currently hard-zeroed (`fesom_tracer_diff.c:158`) | Phase G5: swap `0` for actual array reads, gated on `redi_on` |
| Tracer-advection bolus velocity hook | ✗ missing | Phase G6b |

## Tasks

### Phase G0 — mesh prerequisites (one-time setup)

**Files:** `fesom_mesh.{c,h}`.

- [ ] Add `mesh->mesh_resolution[nod2D]` — Voronoi-cell-radius proxy per node, computed once at setup. Mirror Fortran `oce_mesh.F90` definition (sqrt(area / pi)? Check exact recipe).
- [ ] Add `mesh->nlevels_nod2D_min[n]` and `mesh->ulevels_nod2D_max[n]` — extrema across `nod_in_elem2D`. Compute once at setup.
- [ ] Add `mesh->zbar_3d_n[nod2D * nl]` — per-node 3D bottom-depth array. For linfs this is constant in time (built from `zbar` + per-node `nlevels`); when zlevel/zstar lands later it'll need per-step refresh.
- [ ] **Halo coverage**: all sized `myDim+eDim`; element extrema readers for `n>=myDim` rely on the `eDim` halo of `nod_in_elem2D` being valid (check).
- [ ] No exchange needed — derived from already-exchanged inputs.

### Phase G1 — state allocation

**Files:** `fesom_dyn.{c,h}`, `fesom_aux.{c,h}` (or new `fesom_gm.{c,h}` for GM-specific fields).

- [ ] In `fesom_dyn`: `fer_uv[(myDim+eDim+eXDim) * nl * 2]`, `fer_w[(myDim+eDim) * nl]`.
- [ ] In `fesom_aux` (or new `fesom_gm` struct): `sigma_xy[2 * nl * (myDim+eDim)]`, `neutral_slope[3 * (nl-1) * (myDim+eDim)]`, `slope_tapered[3 * (nl-1) * (myDim+eDim)]`, `fer_tapfac[nl * (myDim+eDim)]`, `fer_gamma[2 * nl * (myDim+eDim)]`, `fer_K[nl * (myDim+eDim)]`, `Ki[nl * (myDim+eDim)]`, `fer_C[(myDim+eDim)]`, `fer_scal[(myDim+eDim)]`.
- [ ] Initialise `fer_K = 1.0e3` (matches `oce_setup_step.F90:934`).
- [ ] Initialise `Ki = 0`, `fer_C = 0`, `slope_tapered = 0`, `fer_tapfac = 0`, etc.

### Phase G2a — MLD1_ind in pressure_bv

**Files:** `fesom_eos.{c,h}` (extend existing `fesom_pressure_bv`).

- [ ] Port `MLD1_ind` derivation from `oce_ale_pressure_bv.F90:420-449` directly into the existing `fesom_pressure_bv` C function. Lives where Fortran lives.
- [ ] Output: `aux->MLD1_ind[(myDim+eDim)]`, integer index per node.
- [ ] **Unconditional** — runs always (cheap, just an index, allows other diagnostics later).
- [ ] No new exchange needed; runs entirely from already-exchanged density / bvfreq.
- [ ] Update `fesom_eos.h` comment block to remove "MLD1_ind: deferred".

### Phase G2b — sigma_xy + neutral_slope + slope_tapered + fer_tapfac

**Files:** new `fesom_gm.{c,h}` (or extend `fesom_eos.c`); called from `fesom_step.c` after pressure_bv.

- [ ] Port `compute_sigma_xy` from `oce_ale_pressure_bv.F90:2851` — element-loop density gradient; also writes node-projected sigma.
- [ ] Port `compute_neutral_slope` from `oce_ale_pressure_bv.F90:2949` — fills `neutral_slope`, `slope_tapered`, `fer_tapfac` using ODM95 tapering (the active branch — confirm vs LDD97 namelist switch).
- [ ] **Halo audit AS WRITTEN**: every write-loop bound checked against Fortran. Per-node arrays cover `myDim+eDim` if Fortran does.
- [ ] Both gated together with the rest of GM/Redi under the `gmredi_on` master switch (G8). When off, the calls are skipped entirely → bit-identity is achievable.

### Phase G3 — init_Redi_GM (per-step coefficient builder)

**Files:** add to `fesom_gm.{c,h}`.

- [ ] Port `init_Redi_GM` (`oce_fer_gm.F90:215-506`).
  - Active sub-features only (see Scope).
  - F1(x,y): area-based resolution scaling for both K and Ki (ord=2).
  - F2(z): `exp(-z/z_ref)` downscaling with `smin=0.6` (`scaling_GMzexp`); reads `mesh->zbar_3d_n` (G0).
  - Redi tapering: yes (`Redi_Ktaper=.true.`), reads `fer_tapfac` (G2b).
  - K_GM_bvref=1 path reads `MLD1_ind` (G2a).
- [ ] Halo audit: lines 503-505 emit `exchange_nod` on `fer_c`, `fer_K`, `Ki`. Mirror in C.

### Phase G4 — fer_solve_Gamma + fer_gamma2vel

**Files:** add to `fesom_gm.{c,h}`.

- [ ] Port `fer_solve_Gamma` (`oce_fer_gm.F90:40-163`). Per-node 1D TDMA for streamfunction over `nzmin..nzmax` extrema (uses G0's `nlevels_nod2D_min` / `ulevels_nod2D_max`). Then `exchange_nod(fer_gamma)`.
- [ ] Port `fer_gamma2vel` (`oce_fer_gm.F90:168-210`). Element loop over `myDim_elem2D`, then `exchange_elem(fer_uv)`.

### Phase G5 — Redi in vertical-implicit tracer diffusion

**Files:** `fesom_tracer_diff.c`.

- [ ] Replace the hard-zeroed `Ki * slope_tapered**2` placeholder at lines 158-184 with the real expressions from `oce_ale_tracer.F90:761-838`. This is the **K33 augmentation** — the diagonal piece of the Redi tensor in z, projected from isopycnal.
- [ ] Gate by a single `redi_on` int (1 if Redi active, 0 otherwise).

### Phase G6a — fer_w inside vert_vel_linfs (HIGHEST RISK)

**Files:** `fesom_ale.c`.

- [ ] Port the GM block of Fortran `oce_ale.F90:2735-2860` — accumulates `fer_w` in the **same edge loop** that builds `dyn->w` and `dyn->w_e` from horizontal divergence. Two parallel accumulators per edge iteration.
- [ ] **Acceptance test for G6a**: with `gmredi_on=0`, the modified `fesom_ale_vert_vel_linfs` must produce `dyn->w` and `dyn->w_e` **bit-identical** to the pre-G6a kernel. Stand up a tiny diff harness on a 1-rank smoke and confirm before any wiring lands.
- [ ] Halo: `fer_w` must be exchange-coherent (Fortran does `exchange_nod(fer_Wvel, partit)` after vert_vel; mirror).

### Phase G6b — bolus velocity add/sub around tracer block

**Files:** `fesom_step.c`.

- [ ] Pre-tracer-advection: `dyn->uv += fer_uv` for elem in `myDim+eDim`; `dyn->w += fer_w` and `dyn->w_e += fer_w` for nod in `myDim+eDim`. Loop bounds match Fortran lines 199-211.
- [ ] Post-tracer-advection: subtract back (lines 284-292).
- [ ] **Important**: in current C, tracer T and S are advected separately with halo exchanges between. The bolus add must wrap **both** tracer iterations. Place the add/sub OUTSIDE the per-tracer loop, not inside it.
- [ ] No exchange needed for the add/sub itself — write happens to halo nodes/elements directly.

### Phase G7a — vertical-explicit Redi (`diff_ver_part_redi_expl`)

**Files:** `fesom_tracer_diff.c` (new function `fesom_diff_ver_part_redi_expl`).

- [ ] Port `oce_ale_tracer.F90:1086-1169`. Off-diagonal Redi tensor projected onto the vertical: `vd_flux(nz) ~ slope_tapered_x * tr_xynodes_x + slope_tapered_y * tr_xynodes_y`, then accumulated into `del_ttf`.
- [ ] Per-node loop (myDim_nod2D bound), reads `slope_tapered`, `Ki`, and the per-tracer horizontal gradients (G7b will produce these — order matters, but G7a doesn't need horizontal Redi flux yet, only the gradient inputs).
- [ ] Gate by `redi_on`. When off, skip the function call.

### Phase G7b — horizontal Redi (`diff_part_hor_redi`)

**Files:** `fesom_tracer_diff.c` (new function `fesom_diff_part_hor_redi`).

- [ ] Port `oce_ale_tracer.F90:~1173+`. Edge-based / element-based loop that adds the horizontal flux contribution to `del_ttf`.
- [ ] **Five partial-cell branches A/B/C/D/E** — handle level mismatch between two adjacent elements (top/bottom level offsets). Mirror each branch literally.
- [ ] Computes `tr_xy` (per-element horizontal tracer gradient) and `tr_xynodes` (per-node projection); these can be local scratch or stored on `tracers->work`.
- [ ] Gate by `redi_on`.
- [ ] **Halo audit**: this is the same code class as ocean FCT bug 2 (commit `20715c6`). Apply the audit AS the function is written, not after.

### Phase G8 — env knobs + driver wiring

**Files:** `fesom_main.c`, `fesom_step.c`.

- [ ] Knobs:
  - `FESOM_NO_GM=1` — set `gm_on = 0`. Skips `fer_solve_Gamma` + `fer_gamma2vel` and the bolus add/sub. `fer_uv`/`fer_w` stay zero so the add/sub is a no-op.
  - `FESOM_NO_REDI=1` — set `redi_on = 0`. Skips G5/G7a/G7b contributions.
  - `FESOM_NO_GMREDI=1` — set both AND set `gmredi_on = 0`. Master switch that ALSO skips `fesom_compute_sigma_xy`, `fesom_compute_neutral_slope`, and `init_Redi_GM`. Required for the bit-identity off-switch test.
- [ ] In `fesom_step.c`, wire the per-step order:
  1. EOS / pressure_bv (existing) — now also writes `aux->MLD1_ind` (G2a, unconditional).
  2. **NEW (gmredi_on)**: `fesom_compute_sigma_xy` → `fesom_compute_neutral_slope` (writes `slope_tapered`, `fer_tapfac`).
  3. **NEW (gmredi_on)**: `fesom_init_redi_gm`.
  4. **NEW (gm_on)**: `fer_solve_Gamma` + `fer_gamma2vel`.
  5. ALE thickness + vert_vel (existing); **G6a** modification builds `fer_w` inside this kernel when `gm_on`.
  6. **NEW (gm_on)**: pre-advection `uv += fer_uv`, `w += fer_w`, `w_e += fer_w`.
  7. Tracer advection (existing FCT) — sees the augmented velocity.
  8. **NEW (redi_on)** — invoke `fesom_diff_part_hor_redi` (G7b) and `fesom_diff_ver_part_redi_expl` (G7a) per tracer; both add to `del_ttf` before the implicit step.
  9. Tracer diffusion — vertical-implicit with G5 K33 augmentation (gated by `redi_on`).
  10. **NEW (gm_on)**: post-advection `uv -= fer_uv`, `w -= fer_w`, `w_e -= fer_w`.

### Phase G9 — Validation

**Files:** new `job_core2_*_gmredi` job templates.

- [ ] Clean build.
- [ ] **Off-switch bit-identity**: 200-step run with `FESOM_NO_GMREDI=1` must produce snapshots **bit-identical** to current pre-Phase-G HEAD (commit `6a13e9d` or its descendant). This is the strongest correctness test for the gating. *MLD1_ind always being computed is acceptable since it's an unread index when GM/Redi off.*
- [ ] **Smoke (200 steps, default-on)**: 1/8/16-rank smoke. Confirm clean exit, no NaNs/blowups, ice tracers still bounded.
- [ ] **Multi-rank consistency**: 8/16/32-rank step-200 ocean and ice statistics within numerical noise (matches Phase E standard).
- [ ] **Eddy-induced transport sanity**: `fer_uv` / `fer_w` magnitudes physical (~1e-3 m/s in mid-latitudes; ~0 in tropics with weak slopes).
- [ ] **Month run (256-rank)**: T/S range physical, no drift relative to no-GMRedi baseline beyond expected eddy mixing scale.
- [ ] **Year run with daily snapshots**: optional after month passes; serves as the second integration test.

### Phase G10 — Memory + plan archival

- [ ] Update or split `project_sea_ice_port_state.md` → `project_gmredi_port_state.md` with the GM/Redi commit hash, validation notes, env knobs, and any partition-dependence caveats.
- [ ] Update `MEMORY.md` index.
- [ ] Move this plan to `docs/plans/completed/`.

## Estimated effort (revised)

| Phase | New C LoC | Risk |
|---|---|---|
| G0 (mesh prereqs) | ~80 | low |
| G1 (state alloc) | ~60 | low |
| G2a (MLD1_ind in pressure_bv) | ~30 | low |
| G2b (sigma_xy + neutral_slope) | ~280 | medium — ODM95 tapering is dense |
| G3 (init_Redi_GM) | ~250 | medium — branch-heavy, easy to misread |
| G4 (fer_solve_Gamma + fer_gamma2vel) | ~150 | medium — TDMA-per-node, halo discipline |
| G5 (Redi K33 in tracer diff) | ~50 | low — slot already prepared |
| G6a (fer_w inside vert_vel_linfs) | ~80 | **HIGH** — invasive mod to validated kernel; bit-identity gate required |
| G6b (bolus add/sub) | ~50 | low |
| G7a (vertical-explicit Redi) | ~120 | medium |
| G7b (horizontal Redi + 5 partial-cell cases) | ~250 | medium-high — same class as ocean FCT halo bugs |
| G8 (knobs + wiring) | ~50 | low |
| G9 (validation jobs) | jobs | low |
| **Total** | ~1450 | medium overall |

## Risk ranking (revised)

1. **G6a** — `fer_w` accumulator inside the existing vert_vel_linfs edge loop. Modifies a validated kernel; must be byte-equivalent when GM is off. This is harder than any of the new function ports because we're changing a working function.
2. **G7b** — horizontal Redi flux on element edges with 5 partial-cell cases. Same code class as the two ocean FCT halo bugs we already paid for. Apply audit at write time.
3. **G7a** — vertical-explicit Redi. Less halo surface but unfamiliar projection of off-diagonal terms.
4. **G3** — many branches in `init_Redi_GM`. Easy to misread one of the active sub-feature gates.
5. **G2b** — ODM95 tapering recipe.
6. **All others** — standard porting risk.

## Halo-write gotchas to pre-empt

Pattern of bugs we've already paid for (sea-ice + ocean FCT):

- Per-node arrays sized `myDim` only that get read at halo by another module. Default to `myDim+eDim` for every fer_*, sigma_xy, etc.
- Write loops bounded `myDim_nod2D` whose downstream readers loop `myDim+eDim`. Apply the audit at write time.
- `exchange_nod` after every Fortran `exchange_nod` call. Lines to mirror:
  - `oce_fer_gm.F90:162` — `exchange_nod(fer_gamma)`
  - `oce_fer_gm.F90:208` — `exchange_elem(fer_uv)`
  - `oce_fer_gm.F90:503-505` — `exchange_nod(fer_c)`, `exchange_nod(fer_K)`, `exchange_nod(Ki)`
  - `oce_ale_pressure_bv.F90:~3070` — neutral_slope/slope_tapered exchanges
  - `oce_ale.F90:2681` — `exchange_nod(fer_Wvel)`

## Off-switch design — bit-identity criterion

`FESOM_NO_GMREDI=1` produces a binary that runs the path of pre-Phase-G
HEAD (commit `6a13e9d` or its descendant). Concretely:

- `gmredi_on=0` skips: `compute_sigma_xy`, `compute_neutral_slope`, `init_Redi_GM`, `fer_solve_Gamma`, `fer_gamma2vel`, the bolus add/sub block, all G7a/G7b/G5 contributions.
- `MLD1_ind` (G2a) runs unconditionally — it's just an integer index, no dependent caller when GM/Redi are off.
- The `fer_w` accumulator inside G6a is gated locally — when `gm_on=0`, the new code path becomes a no-op (and produces zero `fer_w`); the existing `dyn->w` / `dyn->w_e` math is unchanged.

The G9 bit-identity test on a 1-rank 200-step smoke is the strongest
verification that the gating is clean.

## What this gives us

- Mesoscale eddy parameterisation for climate-scale fidelity in T/S
  transport (mid-latitude western boundary currents, ACC).
- Better representation of isopycnal mixing — fixes spurious diapycnal
  flux from horizontal-grid-aligned diffusion.
- Closes one of the three remaining roadmap items (Phase 27 of
  `FRESH_START.md` §16).
