#ifndef FESOM_MOMENTUM_H
#define FESOM_MOMENTUM_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_aux;
struct fesom_dyn;
struct fesom_forcing;

/*
 * Phase 1 momentum substeps. Literal ports of:
 *   compute_vel_rhs           (oce_ale_vel_rhs.F90:35-328)
 *   impl_vert_visc_ale        (oce_ale.F90:3124-3335)
 *   update_vel                (oce_dyn.F90:73-173)
 *   compute_hbar_ale          (oce_ale.F90:1974-2116)
 *
 * Deferred from these subroutines (will be added in later steps/phases):
 *   - momentum_adv_scalar (momadv_opt=2): no advection in Phase 1 first slice
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
                           int                      is_first_step);

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

#endif /* FESOM_MOMENTUM_H */
