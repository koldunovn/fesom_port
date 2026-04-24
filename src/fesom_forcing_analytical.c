#include "fesom_forcing_analytical.h"
#include "fesom_constants.h"
#include "fesom_forcing.h"
#include "fesom_mesh.h"

#include <math.h>
#include <string.h>

void fesom_forcing_set_analytical(const struct fesom_mesh *mesh,
                                  struct fesom_forcing    *forcing,
                                  real_t tau0,
                                  real_t Ly_factor)
{
    /* Zero all forcing arrays first (idempotent). */
    memset(forcing->heat_flux,        0, (size_t)(mesh->myDim_nod2D + mesh->eDim_nod2D) * sizeof(real_t));
    memset(forcing->water_flux,       0, (size_t)(mesh->myDim_nod2D + mesh->eDim_nod2D) * sizeof(real_t));
    memset(forcing->stress_node_surf, 0, (size_t)(mesh->myDim_nod2D + mesh->eDim_nod2D) * 2 * sizeof(real_t));
    memset(forcing->stress_surf,      0, (size_t)(mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D) * 2 * sizeof(real_t));

    const real_t inv_period = 2.0 / Ly_factor;   /* 2π / (Ly_factor * π) — see header */

    /* Wind stress at elements: zonal cosine pattern based on geographic
       centroid latitude. Use mesh->geo_coord_nod2D so the wind tracks the
       physical (un-rotated) latitude — otherwise a rotated mesh would
       produce a wind pattern that doesn't align with the equator. */
    for (int e = 0; e < (mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D); ++e) {
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];
        real_t lat = (mesh->geo_coord_nod2D[2*n0 + 1]
                    + mesh->geo_coord_nod2D[2*n1 + 1]
                    + mesh->geo_coord_nod2D[2*n2 + 1]) / 3.0;
        forcing->stress_surf[2*e + 0] = -tau0 * cos(inv_period * lat);
        forcing->stress_surf[2*e + 1] = 0.0;
    }

    /* stress_node_surf interpolated from element stress (area-weighted mean
       over surrounding elements). impl_vert_visc reads stress_surf at
       elements directly; stress_node_surf is for code paths we haven't
       wired (sea ice, etc.). For Phase 2 keep both consistent. */
    for (int n = 0; n < (mesh->myDim_nod2D + mesh->eDim_nod2D); ++n) {
        int o0 = mesh->nod_in_elem2D_offsets[n];
        int o1 = mesh->nod_in_elem2D_offsets[n + 1];
        real_t tot_a = 0.0, sx = 0.0, sy = 0.0;
        for (int k = o0; k < o1; ++k) {
            int e = mesh->nod_in_elem2D[k];
            real_t a = mesh->elem_area[e];
            tot_a += a;
            sx += a * forcing->stress_surf[2*e + 0];
            sy += a * forcing->stress_surf[2*e + 1];
        }
        if (tot_a > 0.0) {
            forcing->stress_node_surf[2*n + 0] = sx / tot_a;
            forcing->stress_node_surf[2*n + 1] = sy / tot_a;
        }
    }
}
