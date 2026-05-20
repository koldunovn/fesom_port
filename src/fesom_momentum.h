#ifndef FESOM_MOMENTUM_H
#define FESOM_MOMENTUM_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_aux;
struct fesom_dyn;
struct fesom_forcing;
struct fesom_partit;

/*
 * Phase 1 momentum substeps. Literal ports of:
 *   compute_vel_rhs           (oce_ale_vel_rhs.F90:35-328)
 *   impl_vert_visc_ale        (oce_ale.F90:3124-3335)
 *   update_vel                (oce_dyn.F90:73-173)
 *   compute_hbar_ale          (oce_ale.F90:1974-2116)
 *
 *   momentum_adv_scalar       (oce_ale_vel_rhs.F90:335-589, momadv_opt=2)
 *
 * Deferred from these subroutines (will be added in later steps/phases):
 *   - sea-ice (m_ice, m_snow): linfs has use_pice=0 anyway
 *   - sea-level pressure / atmospheric pressure / tides
 *   - kinetic-energy diagnostics (ldiag_ke=false)
 *   - split-explicit subcycling (use_ssh_se_subcycl=false)
 *   - cavity (use_cavity=false in pi/CORE2 namelists)
 *   - water_flux contribution to ssh_rhs_old (linfs branch in Fortran skips it
 *     when which_ALE='linfs' — see oce_ale.F90:2071-2080)
 *   - dhe (delta-hbar-on-element): only used by update_stiff_mat which we
 *     don't call for linfs
 */

/*
 * compute_vel_rhs — Coriolis (AB2) + SSH gradient + PGF.
 *   is_first_step: 1 on the very first call (no AB history yet, use ff=1.0
 *                  instead of ab2). Mirrors Fortran's `lfirst, save` (line 59).
 */
void fesom_compute_vel_rhs(const struct fesom_mesh *mesh,
                           const struct fesom_aux  *aux,
                           struct fesom_dyn        *dyn,
                           int                      is_first_step,
                           struct fesom_partit     *partit);

/*
 * momentum_adv_scalar (momadv_opt=2) — momentum advection on scalar control
 * volumes; adds ke_adv into uv_rhsAB. Called from compute_vel_rhs.
 */
void fesom_momentum_adv_scalar(const struct fesom_mesh *mesh,
                               struct fesom_dyn        *dyn,
                               struct fesom_partit     *partit);

/*
 * impl_vert_visc_ale — implicit vertical viscosity TDMA + wind stress + bottom
 * drag. Solves M * (u_new - u_old) = u_rhs for both u and v components on each
 * element column. Modifies dyn->uv_rhs in place to hold the new (u_new -
 * u_old) increment.
 */
void fesom_impl_vert_visc(const struct fesom_mesh    *mesh,
                          const struct fesom_aux     *aux,
                          const struct fesom_forcing *forcing,
                          struct fesom_dyn           *dyn);

/*
 * update_vel — apply SSH-gradient correction and add uv_rhs to uv:
 *   uv += uv_rhs - g*theta*dt * grad(d_eta)
 */
void fesom_update_vel(const struct fesom_mesh *mesh,
                      struct fesom_dyn        *dyn);

/*
 * compute_hbar_ale — write transport-divergence into ssh_rhs_old, save
 * hbar_old=hbar, then update hbar = hbar_old + ssh_rhs_old*dt/areasvol.
 *
 * Note for linfs: ssh_rhs_old enters the NEXT step's compute_ssh_rhs through
 * the (1-alpha)*ssh_rhs_old term — so even with alpha=1.0 it must be computed
 * correctly each step.
 */
void fesom_compute_hbar(const struct fesom_mesh *mesh,
                        struct fesom_dyn        *dyn);

/*
 * Horizontal harmonic viscosity with easy backscatter (opt_visc=5).
 * Mirror of visc_filt_bcksct (oce_dyn.F90:278-426).
 *   Per non-boundary edge:
 *     vi = dt · max(γ₀, max(γ₁·|du|, γ₂·|du|²)) · sqrt(Σ elem_area)
 *     u_b[el(1)] -= u1·vi/area(el1);  u_b[el(2)] += u1·vi/area(el2)
 *   Smooth u_b → u_c at nodes (area-weighted mean of surrounding cells).
 *   uv_rhs += u_b − easybsreturn · mean(u_c at 3 vertices).
 *
 * Inserted between compute_vel_rhs and impl_vert_visc (FRESH_START.md §5).
 * Skipped at non-boundary edges (Fortran tests myList_edge2D > edge2D_in;
 * we use edge_tri[2*e + 1] >= 0, equivalent for serial).
 */
struct fesom_partit;
void fesom_visc_filt_bcksct(const struct fesom_mesh *mesh,
                            struct fesom_dyn        *dyn,
                            struct fesom_partit     *partit);

/*
 * Biharmonic flow-aware viscosity (opt_visc=7) — visc_filt_bidiff
 * (oce_dyn.F90:592-744). This is what CORE2 runs (work_core opt_visc=7).
 * Two-stage ∇⁴ operator (k⁴ grid-scale damping); required for dt=1800
 * stability where the harmonic opt_visc=5 lets a 2Δx mode grow.
 */
void fesom_visc_filt_bidiff(const struct fesom_mesh *mesh,
                            struct fesom_dyn        *dyn,
                            struct fesom_partit     *partit);

#endif /* FESOM_MOMENTUM_H */
