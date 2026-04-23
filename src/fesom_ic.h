#ifndef FESOM_IC_H
#define FESOM_IC_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;
struct fesom_tracers;

/*
 * Phase 1 initial conditions.
 *
 * fesom_ic_thickness:
 *   Initialise resting layer thicknesses from zbar in the linfs convention
 *   (no partial cells, no cavity yet — those branches in init_thickness_ale
 *   are deferred to Phase 3+):
 *     hnode[n][nz]    = zbar[nz] - zbar[nz+1]    for nz < nlevels_nod2D[n] - 1
 *     hnode_new[n][:] = hnode[n][:]              (linfs invariant per §6)
 *     helem[e][nz]    = zbar[nz] - zbar[nz+1]    for nz < nlevels[e]    - 1
 *     hbar = hbar_old = 0
 *     ssh_rhs_old = 0
 *     dyn->eta_n = 0
 *
 * fesom_ic_tracers_constant:
 *   Set T (FESOM_TRACER_T) and S (FESOM_TRACER_S) to constant values in
 *   every wet layer at every node. Below the column (nz >= nlevels_nod2D-1)
 *   the values stay zero.
 */
void fesom_ic_thickness(const struct fesom_mesh *mesh,
                        struct fesom_mesh       *mesh_state,
                        struct fesom_dyn        *dyn);

void fesom_ic_tracers_constant(const struct fesom_mesh *mesh,
                               struct fesom_tracers    *tracers,
                               real_t T_val, real_t S_val);

/*
 * fesom_ic_tracer_T_blob:
 *   Add a Gaussian temperature anomaly centred at geographic (lon0, lat0)
 *   with horizontal e-folding sigma_deg degrees and vertical e-folding
 *   sigma_z metres. Maximum amplitude amp_C added at the surface centre,
 *   tapering with horizontal distance and depth. Operates additively on
 *   every wet layer.
 */
void fesom_ic_tracer_T_blob(const struct fesom_mesh *mesh,
                            struct fesom_tracers    *tracers,
                            real_t lon0_deg, real_t lat0_deg,
                            real_t sigma_deg, real_t sigma_z,
                            real_t amp_C);

#endif /* FESOM_IC_H */
