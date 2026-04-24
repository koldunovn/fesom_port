#include "fesom_ic.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_mesh.h"
#include "fesom_tracers.h"

#include <math.h>
#include <string.h>

void fesom_ic_thickness(const struct fesom_mesh *mesh_in,
                        struct fesom_mesh       *mesh,
                        struct fesom_dyn        *dyn)
{
    /* mesh_in vs mesh: same pointer, two roles. mesh_in carries read-only
       counts and zbar; mesh receives the writes into hnode/helem/hbar.
       Splitting them in the prototype reminds the call site that this
       routine *does* mutate the mesh struct. */
    (void)mesh_in;

    const int nl = mesh->nl;

    /* hnode at nodes — full-cell linfs. Loop over interior + halo so every
     * local node has thickness defined (each rank computes the same value
     * from its zbar / nlevels independently — no exchange needed). */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int E = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;   /* 1-based → 0-based */
        int nzmax = mesh->nlevels_nod2D[n] - 1;   /* exclusive bound  */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t h = mesh->zbar[nz] - mesh->zbar[nz + 1];
            mesh->hnode    [FESOM_NODE3D(n, nz, nl)] = h;
            mesh->hnode_new[FESOM_NODE3D(n, nz, nl)] = h;   /* linfs invariant */
        }
    }

    /* helem at cells — full-cell linfs.
       For Phase 1 (full cells, no partial), this is identical to the
       (1/3) Σ hnode formula since all 3 vertices have the same column thickness
       in the layers where the cell exists. We use the literal Fortran form
       (zbar diff) to mirror init_thickness_ale lines 849-852. */
    for (int e = 0; e < E; ++e) {
        int nzmin = mesh->ulevels[e] - 1;
        int nzmax = mesh->nlevels[e] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            mesh->helem[FESOM_ELEM3D(e, nz, nl)] =
                mesh->zbar[nz] - mesh->zbar[nz + 1];
        }
    }

    /* hbar / hbar_old / ssh_rhs_old / eta_n / d_eta — at rest. Sized
     * myDim+eDim per rank (allocations updated in fesom_dyn_alloc /
     * fesom_mesh_alloc_state). */
    memset(mesh->hbar,        0, (size_t)N * sizeof(real_t));
    memset(mesh->hbar_old,    0, (size_t)N * sizeof(real_t));
    memset(dyn->eta_n,        0, (size_t)N * sizeof(real_t));
    memset(dyn->d_eta,        0, (size_t)N * sizeof(real_t));
    memset(dyn->ssh_rhs,      0, (size_t)N * sizeof(real_t));
    memset(dyn->ssh_rhs_old,  0, (size_t)N * sizeof(real_t));
}

void fesom_ic_tracers_constant(const struct fesom_mesh *mesh,
                               struct fesom_tracers    *tracers,
                               real_t T_val, real_t S_val)
{
    const int nl = mesh->nl;
    real_t *T = tracers->data[FESOM_TRACER_T].values;
    real_t *S = tracers->data[FESOM_TRACER_S].values;

    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t idx = FESOM_NODE3D(n, nz, nl);
            T[idx] = T_val;
            S[idx] = S_val;
        }
    }
}

void fesom_ic_tracer_T_blob(const struct fesom_mesh *mesh,
                            struct fesom_tracers    *tracers,
                            real_t lon0_deg, real_t lat0_deg,
                            real_t sigma_deg, real_t sigma_z,
                            real_t amp_C)
{
    const int nl = mesh->nl;
    real_t *T = tracers->data[FESOM_TRACER_T].values;

    /* Wrap lon0 to [-180, 180) for distance compatibility. */
    while (lon0_deg >= 180.0)  lon0_deg -= 360.0;
    while (lon0_deg <  -180.0) lon0_deg += 360.0;

    const real_t inv_sig2_h = 1.0 / (sigma_deg * sigma_deg);
    const real_t inv_sig2_z = 1.0 / (sigma_z   * sigma_z);

    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    for (int n = 0; n < N; ++n) {
        real_t lon = (real_t)mesh->geo_coord_nod2D[2*n + 0] / FESOM_RAD;
        real_t lat = (real_t)mesh->geo_coord_nod2D[2*n + 1] / FESOM_RAD;
        /* Wrap node lon similarly */
        while (lon >= 180.0)  lon -= 360.0;
        while (lon <  -180.0) lon += 360.0;

        real_t dlon = lon - lon0_deg;
        if (dlon >  180.0) dlon -= 360.0;
        if (dlon < -180.0) dlon += 360.0;
        real_t dlat = lat - lat0_deg;

        /* Cheap small-circle correction so the patch isn't grossly
           stretched at high lats. */
        real_t cos_lat0 = cos(lat0_deg * FESOM_RAD);
        real_t r2_h = (dlon*cos_lat0)*(dlon*cos_lat0)*inv_sig2_h
                    + dlat*dlat*inv_sig2_h;

        if (r2_h > 16.0) continue;       /* skip beyond 4 sigma horizontally */

        real_t prof_h = exp(-0.5 * r2_h);

        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t z = (real_t)mesh->Z[nz];      /* Z negative downward */
            real_t prof_z = exp(-0.5 * z*z * inv_sig2_z);
            T[FESOM_NODE3D(n, nz, nl)] += amp_C * prof_h * prof_z;
        }
    }
}
