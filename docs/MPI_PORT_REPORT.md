# MPI Port — Open Issue Report

Date: 2026-04-24
Status: **partially working** — single rank stable, multi-rank fails progressively faster with more partitions.

## Problem

CORE2 multi-rank runs accumulate numerical drift at partition boundaries that
eventually produces NaN in the SSH conjugate-gradient solver. The earlier
fixes (slice 30d follow-up: PHC `extrap_nod3D` MPI port, preconditioner OOB
read, etc.) brought CORE2 16-rank to bit-identical step-1 match against the
1-rank reference, but longer runs reveal a slow accumulating drift specific
to multi-rank, with rate roughly proportional to partition count.

## Symptoms

### Failure trajectory by rank count (CORE2 mesh, dt=500s)

| ranks | step-1 `max(pgf)` | step-1 `max(rs)` | failure step | notes |
|-------|-------------------|-------------------|---------------|-------|
| 1     | 1.93e-04          | 3.49e-05          | — (100+ stable) | baseline |
| 8     | **1.93e-04**      | **3.49e-05**      | 147 → CG NaN   | bit-identical to 1-rank at step 1 |
| 16    | 2.56e-04          | 4.62e-05          | ~95 → CG NaN   | step-1 immediate divergence |
| 32    | (similar to 16)   | (similar to 16)   | 85 → CG NaN    | step-1 divergence |
| 144   | 2.20e-04          | 4.62e-05          | 58 → BLOWUP    | step-1 divergence |

### Visible artifacts in snapshots

- `eta` (sea surface height) shows clear discontinuities along partition boundaries
- `density` (= EOS(T,S)) shows the same boundaries less strongly
- `T` and `S` show faint partition seams visible only with adjusted colour bars

### Specific oddity: `rs` values cluster in two groups

- {1-rank, 8-rank}: max(relax_salt) = 3.49e-05
- {16-rank, 32-rank, 144-rank}: max(relax_salt) = 4.62e-05

Same value across vastly different partition counts (16/32/144) suggests an MPI
reduction whose order depends on rank count crossing some power-of-two
threshold, but the **magnitude** of the difference (~30%) is far larger than
expected from FP-rounding alone.

### `max(w)` divergence starts immediately

At step 18, with `max(uv)` bit-identical between 1-rank and 16-rank (0.855):
- 1-rank: max(w) = 5.54e-3
- 16-rank: max(w) = 8.17e-3 (47% bigger)

Same uv but different vertical-divergence integral. Inputs to vert_vel_ale (uv,
helem, edge_cross_dxdy) are halo-exchanged or computed once. Why w differs is
not yet explained.

## What we verified is working

### Mesh metadata bit-identical across ranks
Probe at 5 GIDs (78867, 9694, 23377, 12906, 20478) on 16-rank shows OWN/HALO
values for `area[0]`, `areasvol[0]`, surrounding-element count are
**bit-identical**. So mesh setup, scatter, and `compute_metrics` are correct.

### Single-rank baseline stable
1-rank with current code runs 100 steps with values matching the pre-MPI
1-rank reference: max(uv)=1.79, max(eta)=2.23, S=[5.65, 41.12].

### 8-rank step-1 bit-identical to 1-rank
All 16 instrumented fields (uv, eta, w, T-range, S-range, stress, hf, wf, vs,
rs, hp, pgf, rho, bv-range, Kv, Av) match 1-rank exactly through step 5.
Means dist_8 partition correctness is solid.

### Halo identity test passes
The halo exchange machinery itself is correct (round-trip GID test passes
positive + negative cases on 8 ranks).

### PHC IC `extrap_nod3D` MPI propagation
After porting Fortran's `exchange_nod` between sweeps + `MPI_Allreduce` on
glob_max, the IC produces step-1 values bit-identical to 1-rank.

## What we tried and ruled out

- **All known missing halo exchanges added**: see "Halo exchanges audit"
  section below — every per-step exchange Fortran does is in our code.
- **Preconditioner OOB read**: fixed (sized `diag_values` to N_alloc and
  halo-exchanged). Runtime values unchanged on 16-rank.
- **bulk_compute halo exchange**: added for `stress_node_surf`, `heat_flux`,
  `water_flux`. No observable effect *at the time* — but this fix was
  INCOMPLETE and the gap became the dt=1800 blocker. `bulk_compute` also
  writes `Ch_atm_oce`/`Ce_atm_oce` (NCAR exchange coefficients) and looped
  `myDim` only, while Fortran `ncar_ocean_fluxes_mode`
  (gen_bulk_formulae.F90:181) loops `myDim+eDim`. Exchanging the *fluxes* but
  not the *coefficients* left Ch/Ce stale at eDim; sea-ice `thermodynamics`
  reads them at eDim → divergent `a_ice` → divergent `stress_surf` on
  replicated boundary elements → non-conservative SSH RHS → dt-dependent
  blow-up. Real fix (2026-05-20): loop `bulk_compute` over `myDim+eDim` (one
  line). Lesson: exchanging *some* outputs of a routine is not the same as
  covering *every* field a downstream consumer reads at halo. See
  PORT_EXPERIENCE_REPORT.md §3.1 "most insidious instance".
- **`compute_hbar` / `compute_ssh_rhs_linfs` / `ale_vert_vel_linfs` edge
  loops**: changed `myDim+eDim_edge2D` → `myDim_edge2D` to match Fortran.
  No observable effect on 16-rank values.
- **`pp_mixing` / `mo_convect` node loops**: changed `myDim` → `myDim+eDim`
  (matches Fortran). No observable effect.
- **`fct_LO` halo exchange**: added before Zalesak limiter. No effect.
- **`visc_filt_bcksct` U_b/V_b and U_c/V_c exchanges**: added (Fortran has
  these). No effect.
- **`hbar_old` halo write**: extended to `myDim+eDim`. No effect.
- **Element-area exchange ordering**: refactored so `elem_area` halo-exchanged
  BEFORE per-node area accumulation. No effect.
- **`nod_in_elem2D` extended over halo elements**: now mirrors Fortran's
  full-extent adjacency. No effect.
- **`elem_nodes` for halo elements**: scatter populates myDim+eDim+eXDim with
  -1 sentinel for outer-halo unmappable vertices. No effect on dynamics.

## Halo exchanges currently in place (per timestep)

These are the fields we exchange and where, in execution order. **All are
present** — we audited against every `call exchange_*` in Fortran's
`oce_dyn.F90`, `oce_ale.F90`, `oce_ale_pressure_bv.F90`, `oce_ale_tracer.F90`,
`oce_adv_tra_*.F90`, `oce_ale_mixing_pp.F90`, `oce_mo_conv.F90`.

| stage | field(s) exchanged | Fortran ref |
|-------|--------------------|--------------|
| after pressure_bv | `density_m_rho0`, `hpressure`, `bvfreq`, `sw_alpha`, `sw_beta` | `oce_ale_pressure_bv.F90:2844-45`, `sw_alpha_beta` |
| after pressure_force | `pgf_x`, `pgf_y` | (we exchange; Fortran doesn't need to since downstream reads only myDim) |
| after compute_vel_nodes | `uvnode` | `oce_dyn.F90:225` |
| after pp_mixing | `Kv` (nod3D), `Av` (elem3D) | inside Fortran `oce_mixing_PP` |
| after mo_convect | `Kv` (nod3D), `Av` (elem3D) | inside Fortran `mo_convect` |
| after compute_vel_rhs | `uv_rhs` (elem3D vec), `uv_rhsAB` | (we exchange; Fortran doesn't need) |
| after visc_filt_bcksct (internal) | `u_b`, `v_b` (elem3D); then `u_c`, `v_c` (nod3D) | `oce_dyn.F90:367, 395` |
| after visc_filt_bcksct (external) | `uv_rhs` | (we exchange; Fortran doesn't need) |
| after impl_vert_visc | `uv_rhs` | (we exchange; Fortran doesn't need) |
| after compute_ssh_rhs_linfs | `ssh_rhs` | `oce_ale.F90:1954` |
| inside CG | `pp` (each iter, before SpMV); `rr` (after residual update); `X` (final) | `solver.F90` parallel CG |
| before CG (preconditioner) | `diag_values` (one-time at preconditioner build) | `solver.F90:73` |
| after CG | `d_eta` | `oce_ale.F90:3117` |
| after update_vel | `uv` | `oce_dyn.F90:172` |
| after compute_hbar | `ssh_rhs_old`, `hbar` | `oce_ale.F90:2078, 2102` |
| after ALE vert_vel | `w` | `oce_ale.F90:2679` |
| after ALE compute_cflz | `cfl_z` | (we exchange) |
| after ALE wvel_split | `w_e`, `w_i` | (we exchange) |
| after tracer FCT (each) | `T`, `S` (per tracer) | `oce_ale_tracer.F90:268` |
| inside tracer FCT | `fct_LO` (after compute), `fct_plus`/`fct_minus` (before limiter) | `oce_adv_tra_driver.F90:294`, `oce_adv_tra_fct.F90:401` |
| after impl_vert_diff | `T`, `S` | (we exchange) |
| after commit_thickness | `hnode`, `helem` | `oce_ale.F90:1027, 1249` |

### One-time halo exchanges at startup

- `mesh.elem_area` (ELEM2D and ELEM2D_FULL)
- `mesh.elem_cos`, `mesh.metric_factor`, `mesh.coriolis`
- `mesh.elem_center_x`, `mesh.elem_center_y`
- `MPI_Allreduce` on `mesh.ocean_area`

### MPI Allreduce sums (per timestep)

- CG: `pp·App`, `rr·zz`, `rr·rr` per CG iteration
- `integrate_nod_2D` for: virtual_salt mean, relax_salt mean, water_flux balance

## Hypotheses for the remaining drift

Ordered by what we still consider plausible:

1. **Allreduce summation order non-determinism amplified by chaos**
   - At 8 ranks, FP rounding from MPI tree reduction happens to land bit-equal
     to single-rank serial sum (lucky coincidence at small N).
   - At 16+ ranks, deeper tree → different rounding → 1e-15 perturbation per
     reduction.
   - Ocean dynamics is chaotic; small perturbations grow exponentially.
   - The fact that {1, 8} and {16, 32, 144} fall into TWO discrete groups
     for `max(rs)` is a smoking-gun signature of this hypothesis.
   - Counter-argument: the magnitude (4.62e-5 vs 3.49e-5) is too large for
     pure FP rounding (would expect ~1e-15 difference, not 1e-5).

2. **Replicated-element computation discrepancy**
   - For elements in BOTH ranks A's and B's `myDim_elem2D` (replication for
     boundary cells), gradient_sca is computed independently on each rank.
     If somehow the local computation differs, pgf at that element differs.
   - We verified mesh metadata at NODES is bit-identical across OWN/HALO,
     but did NOT directly probe `gradient_sca` or `metric_factor` at
     replicated ELEMENTS.

3. **Uninitialized array slot read in our diagnostic**
   - `pgf_x` allocated `myDim+eDim+eXDim`. We compute it for `myDim` only.
     If our diagnostic accidentally scans halo entries (we don't think we do,
     but it's worth re-verifying), we could pick up garbage.

4. **Unaudited Fortran exchange we missed**
   - Looking at `oce_ale_pressure_bv.F90`, the `sigma_xy` exchange (line 2937)
     is in `compute_sigma_xy` for Redi — we don't run Redi, so should be fine.
     `compute_neutral_slope` exchanges `neutral_slope`/`slope_tapered` — also
     Redi-only, skipped. We believe nothing per-step is missing.

## Suggested next debug steps

1. **Probe `gradient_sca[6*e + k]` at a replicated element** — if rank A and
   rank B compute different values for the same physical element, this
   pinpoints (2). Easy to add to the existing `mesh_metadata_consistency`
   probe in `fesom_main.c`.

2. **Probe `hpressure[n, 0]` at a boundary node** between 1-rank and
   16-rank — if hpressure differs at any node despite T/S being identical,
   the `pressure_bv` accumulation (sequential downward integral) is producing
   different rounding. Add a one-shot probe immediately after the
   `pressure_bv` call at step 1.

3. **Use `MPI_Allreduce(SUM, REPRODUCIBLE)`** — OpenMPI supports
   `MPI_REPRODUCIBLE` reductions that guarantee bit-identical results
   regardless of rank count. If switching to this collapses the 16-rank
   `max(rs)` to 3.49e-5 (matching 1-rank), the bug is purely Allreduce noise.

4. **Compare snapshot file contents** between 1-rank and 16-rank step 0 with
   `xarray.open_dataset(...).var.values` — should be bit-identical at step 0.
   If not, the gather IO writer has a bug.

5. **Bisect commit history** — git checkout each prior `Phase 4 step 30*`
   commit and run a 16-rank step-1 diagnostic. Find the commit where
   step-1 `max(pgf)` first diverged from 1-rank.

## Build & run quick reference

```bash
# Build (Levante, gcc-11.2.0 + openmpi/4.1.2 + serial netcdf-c)
cd /home/a/a270088/port2/fesom2_port && bash -l configure.sh

# Diagnostic: 16-rank, 100 steps, per-step prints, snapshots every 10
sbatch /home/a/a270088/port2/fesom2_port/job_core2_16_dbg

# Production attempt: 8-rank model month
sbatch /home/a/a270088/port2/fesom2_port/job_core2_8

# 1-rank baseline
sbatch /tmp/job_core2_1r   # see body of last 1-rank run script
```

## Currently usable

- 1-rank for full model month, bit-identical to pre-MPI port
- 8-rank for ~1 model day (147 steps × dt=500s ≈ 20 hours model time)
- 16+ ranks: short runs only (50–95 steps before NaN)
