#ifndef FESOM_DYN_H
#define FESOM_DYN_H

#include "fesom_types.h"

struct fesom_mesh;

/*
 * Phase 1 subset of Fortran t_dyn (MOD_DYN.F90:44+).
 *
 *   uv         [elem2D * nl * 2]   instant horizontal velocity (u, v)
 *   uv_rhs     [elem2D * nl * 2]   current momentum RHS
 *   uv_rhsAB   [elem2D * nl * 2]   AB2 history (single slot — see §14.4)
 *   w          [nod2D  * nl]       vertical velocity at scalar locations
 *   eta_n      [nod2D]             SSH at this step
 *   d_eta      [nod2D]             SSH increment
 *   ssh_rhs    [nod2D]             current SSH RHS
 *   ssh_rhs_old[nod2D]             previous SSH RHS (for AM2 blending)
 *
 * Deferred to later phases: fer_uv, fer_w, w_e, w_i, w_old, cfl_z, uvnode,
 * uvnode_rhs, all se_* split-explicit arrays, iceberg arrays.
 */
typedef struct fesom_dyn {
    real_t *uv;
    real_t *uv_rhs;
    real_t *uv_rhsAB;
    real_t *w;
    real_t *w_i;        /* implicit part of w (compute_Wvel_split). Zero in
                           Phase 1 because wsplit only kicks in when CFL_z >
                           wsplit_maxcfl; with no forcing that never triggers,
                           so w_i stays at the calloc'd zero — bit-identical to
                           Fortran behaviour for stable cases (compute_Wvel_split
                           sets w_i=0 in the else branch, oce_ale.F90:3065). */
    real_t *eta_n;
    real_t *d_eta;
    real_t *ssh_rhs;
    real_t *ssh_rhs_old;

    real_t *uvnode;     /* [nod2D * nl * 2] u,v at NODES — interpolated from
                           element-center UV via area-weighted mean. Required
                           by PP mixing (Fortran tracers%data%uvnode is held
                           in dynamics%uvnode). */

    int AB_order;   /* fixed at 2 in Phase 1 */
} fesom_dyn;

void fesom_dyn_alloc(fesom_dyn *d, const struct fesom_mesh *mesh);
void fesom_dyn_free (fesom_dyn *d);

#endif /* FESOM_DYN_H */
