#include "fesom_ic.h"
#include "fesom_dyn.h"
#include "fesom_mesh.h"
#include "fesom_tracers.h"

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

    /* hnode at nodes — full-cell linfs */
    for (int n = 0; n < mesh->nod2D; ++n) {
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
    for (int e = 0; e < mesh->elem2D; ++e) {
        int nzmin = mesh->ulevels[e] - 1;
        int nzmax = mesh->nlevels[e] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            mesh->helem[FESOM_ELEM3D(e, nz, nl)] =
                mesh->zbar[nz] - mesh->zbar[nz + 1];
        }
    }

    /* hbar / hbar_old / ssh_rhs_old / eta_n / d_eta — at rest */
    memset(mesh->hbar,        0, (size_t)mesh->nod2D * sizeof(real_t));
    memset(mesh->hbar_old,    0, (size_t)mesh->nod2D * sizeof(real_t));
    memset(dyn->eta_n,        0, (size_t)mesh->nod2D * sizeof(real_t));
    memset(dyn->d_eta,        0, (size_t)mesh->nod2D * sizeof(real_t));
    memset(dyn->ssh_rhs,      0, (size_t)mesh->nod2D * sizeof(real_t));
    memset(dyn->ssh_rhs_old,  0, (size_t)mesh->nod2D * sizeof(real_t));
}

void fesom_ic_tracers_constant(const struct fesom_mesh *mesh,
                               struct fesom_tracers    *tracers,
                               real_t T_val, real_t S_val)
{
    const int nl = mesh->nl;
    real_t *T = tracers->data[FESOM_TRACER_T].values;
    real_t *S = tracers->data[FESOM_TRACER_S].values;

    for (int n = 0; n < mesh->nod2D; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t idx = FESOM_NODE3D(n, nz, nl);
            T[idx] = T_val;
            S[idx] = S_val;
        }
    }
}
