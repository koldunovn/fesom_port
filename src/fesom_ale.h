#ifndef FESOM_ALE_H
#define FESOM_ALE_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;

/*
 * ALE (Arbitrary Lagrangian-Eulerian) step. Phase 1 implements only the
 * `which_ALE = 'linfs'` branch:
 *
 *   - hnode_new = hnode (no thickness motion)
 *   - w computed by accumulating horizontal edge fluxes over the two
 *     adjacent cells (Gauss-law div integration), then cumsum vertically
 *     from bottom up, then division by area(nz, n) to get m/s.
 *
 * Mirror of the linfs path of vert_vel_ale (oce_ale.F90:2132+) plus
 * update_thickness_ale (oce_ale.F90:1035+).
 *
 * Future: zlevel surface-flux correction, zstar layer scaling, Fer_GM
 * bolus contribution, partial cells.
 */
void fesom_ale_vert_vel_linfs(const struct fesom_mesh *mesh,
                              struct fesom_dyn        *dyn);

/* hnode_new := hnode (linfs invariant). Called at the start of the ALE step. */
void fesom_ale_thickness_linfs(struct fesom_mesh *mesh);

/* Commit at end of timestep: hnode := hnode_new and helem := mean of 3 vertex hnode. */
void fesom_ale_commit_thickness(struct fesom_mesh *mesh);

/*
 * Vertical CFL (compute_CFLz, oce_ale.F90:2935-3022):
 *   CFL_z[nz, n] gets contributions from W at the two interfaces above and
 *   below layer nz, scaled by the appropriate layer thickness.
 *
 * wsplit (compute_Wvel_split, oce_ale.F90:3026-3074):
 *   For each interface where CFL_z > wsplit_maxcfl, split W into explicit
 *   (w_e, used by tracer advection) and implicit (w_i, used by impl_vert_visc)
 *   parts so the tracer scheme stays CFL-bounded. Default behaviour
 *   (CFL ≤ maxcfl) is w_e = w, w_i = 0.
 */
void fesom_ale_compute_cflz(const struct fesom_mesh *mesh,
                            struct fesom_dyn        *dyn);

void fesom_ale_compute_wvel_split(const struct fesom_mesh *mesh,
                                  struct fesom_dyn        *dyn);

#endif /* FESOM_ALE_H */
