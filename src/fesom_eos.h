#ifndef FESOM_EOS_H
#define FESOM_EOS_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_aux;
struct fesom_tracers;
struct fesom_partit;

/*
 * Jackett-McDougall (JM, EOS80) equation of state — split form. Returns the
 * three secant-bulk-modulus coefficients (bulk_0, bulk_pz, bulk_pz2) and the
 * potential density rhopot. In-situ density is then assembled by the caller as
 *
 *   bulk = bulk_0 + z*(bulk_pz + z*bulk_pz2)
 *   rho  = bulk * rhopot / (bulk + 0.1 * z * state_eq_int)
 *
 * Mirror of densityJM_components (oce_ale_pressure_bv.F90:2605-2669).
 */
void fesom_eos_jm_components(real_t t, real_t s,
                             real_t *bulk_0, real_t *bulk_pz, real_t *bulk_pz2,
                             real_t *rhopot);

/*
 * Phase 1 subset of pressure_bv (oce_ale_pressure_bv.F90:194-501):
 *   - JM-EOS at every wet level
 *   - density_m_rho0  = rho_in_situ - density_0  (use_density_ref=.false.)
 *   - hpressure       (linfs branch, full cells, no cavity)
 *   - bvfreq          (N², padded; horizontal smoothing done by fesom_smooth_nod3D,
 *                      called from fesom_step after the bvfreq halo exchange)
 *   - dbsfc           (buoyancy re surface for KPP bldepth; stored unconditionally,
 *                      PP never reads it — Fortran oce_ale_pressure_bv.F90:332,337,339)
 *   - MLD1_ind        (Large et al. 1997 / FESOM 1.4 — Phase G2a;
 *                      stored 0-based, used by GM init_Redi_GM)
 *
 * Deferred:
 *   - vertical N² smoothing (N2smth_v=.false. in CORE2)
 *   - MLD2/MLD3 (Levitus, Griffies)
 *   - density_dmoc (diagnostic)
 *   - cavity / partial-cell branches
 */
void fesom_pressure_bv(const struct fesom_tracers *tracers,
                       const struct fesom_mesh    *mesh,
                       struct fesom_aux           *aux);

/*
 * Area-weighted node-patch horizontal smoother — port of smooth_nod3D
 * (gen_support.F90:99-198). Smooths `arr` [(myDim+eDim)*nl] in place with
 * `n_smooth` sweeps, each an element-patch area-weighted average followed by a
 * halo exchange. The caller must supply a valid halo `arr` on entry. Used for N²
 * (N2smth_h=.true., N2smth_hidx=1).
 */
void fesom_smooth_nod3D(real_t *arr, int nl, int n_smooth,
                        const struct fesom_mesh *mesh, struct fesom_partit *p);

/*
 * Pressure gradient force at elements (Phase 1: linfs + full cells).
 *   pgf_x[e][nz] = Σ_i ∇N_i^x * hpressure[nz][V_i(e)] / density_0
 *   pgf_y[e][nz] = Σ_i ∇N_i^y * hpressure[nz][V_i(e)] / density_0
 *
 * Mirror of pressure_force_4_linfs_fullcell (oce_ale_pressure_bv.F90:575-614).
 * Other variants (Shchepetkin / cubic-spline / partial cells / cavity)
 * deferred to later phases — namelist default is 'shchepetkin' but the
 * dispatcher picks linfs_fullcell when use_partial_cell=.false. AND
 * use_cavity_partial_cell=.false. (both true in Phase 1).
 */
void fesom_pressure_force_linfs_fullcell(const struct fesom_mesh *mesh,
                                         struct fesom_aux        *aux);

/*
 * Compute thermal expansion (sw_alpha) and saline contraction (sw_beta)
 * coefficients per node per level, from McDougall (1987). Mirror of
 * oce_ale_pressure_bv.F90:2751-2846 sw_alpha_beta. Outputs into
 * aux->sw_alpha and aux->sw_beta.
 *
 * Halo-exchanged at end (the caller may add an explicit exchange too).
 */
void fesom_compute_sw_alpha_beta(const struct fesom_tracers *tracers,
                                 const struct fesom_mesh    *mesh,
                                 struct fesom_aux           *aux);

#endif /* FESOM_EOS_H */
