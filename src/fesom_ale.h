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

#endif /* FESOM_ALE_H */
