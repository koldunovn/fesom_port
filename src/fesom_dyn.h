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
    real_t *w_e;        /* explicit part of w (for tracer advection)         */
    real_t *w_i;        /* implicit part of w (for impl_vert_visc)
                           When CFL_z[nz,n] ≤ wsplit_maxcfl: w_e = w, w_i = 0
                           (compute_Wvel_split, oce_ale.F90:3026-3074).       */
    real_t *cfl_z;      /* [nod2D * nl] vertical CFL at interfaces            */
    real_t *eta_n;
    real_t *d_eta;
    real_t *ssh_rhs;
    real_t *ssh_rhs_old;

    real_t *uvnode;     /* [nod2D * nl * 2] u,v at NODES — interpolated from
                           element-center UV via area-weighted mean. Required
                           by PP mixing (Fortran tracers%data%uvnode is held
                           in dynamics%uvnode). */

    /* Work arrays for visc_filt_bcksct (opt_visc=5). Mirror of
       dynamics%work%u_b/v_b/u_c/v_c (MOD_DYN.F90 t_dyn_work). */
    real_t *u_b;        /* [elem2D * nl] per-element raw harmonic-visc update */
    real_t *v_b;
    real_t *u_c;        /* [nod2D  * nl] node-smoothed visc-update            */
    real_t *v_c;

    int AB_order;   /* fixed at 2 in Phase 1 */
} fesom_dyn;

void fesom_dyn_alloc(fesom_dyn *d, const struct fesom_mesh *mesh);
void fesom_dyn_free (fesom_dyn *d);

#endif /* FESOM_DYN_H */
