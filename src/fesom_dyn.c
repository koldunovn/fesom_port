#include "fesom_dyn.h"
#include "fesom_mesh.h"
#include "fesom_constants.h"

#include <stdlib.h>
#include <string.h>

void fesom_dyn_alloc(fesom_dyn *d, const struct fesom_mesh *mesh)
{
    memset(d, 0, sizeof(*d));
    d->AB_order = FESOM_PHASE1_AB_ORDER;

    /* MPI: size for full local extent (myDim+eDim for nodes,
     * myDim+eDim+eXDim for elements). */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int E = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    size_t e_nl_2 = (size_t)E * (size_t)mesh->nl * 2;
    size_t n_nl   = (size_t)N * (size_t)mesh->nl;
    size_t n      = (size_t)N;

    d->uv          = calloc(e_nl_2, sizeof(real_t));
    d->uv_rhs      = calloc(e_nl_2, sizeof(real_t));
    d->uv_rhsAB    = calloc(e_nl_2, sizeof(real_t));
    d->w           = calloc(n_nl,   sizeof(real_t));
    d->w_e         = calloc(n_nl,   sizeof(real_t));
    d->w_i         = calloc(n_nl,   sizeof(real_t));
    d->cfl_z       = calloc(n_nl,   sizeof(real_t));
    d->uvnode      = calloc(n_nl * 2, sizeof(real_t));
    d->uvnode_rhs  = calloc(n_nl * 2, sizeof(real_t));
    size_t e_nl = (size_t)E * (size_t)mesh->nl;
    d->u_b = calloc(e_nl, sizeof(real_t));
    d->v_b = calloc(e_nl, sizeof(real_t));
    d->u_c = calloc(n_nl, sizeof(real_t));
    d->v_c = calloc(n_nl, sizeof(real_t));
    d->eta_n       = calloc(n,      sizeof(real_t));
    d->d_eta       = calloc(n,      sizeof(real_t));
    d->ssh_rhs     = calloc(n,      sizeof(real_t));
    d->ssh_rhs_old = calloc(n,      sizeof(real_t));
    /* GM bolus velocity (Phase G1) — zero until G4/G6a write to them. */
    d->fer_uv      = calloc(e_nl_2, sizeof(real_t));
    d->fer_w       = calloc(n_nl,   sizeof(real_t));
    FESOM_CHECK(d->uv && d->uv_rhs && d->uv_rhsAB && d->w && d->w_e && d->w_i
             && d->cfl_z && d->uvnode && d->uvnode_rhs
             && d->u_b && d->v_b && d->u_c && d->v_c
             && d->eta_n && d->d_eta && d->ssh_rhs && d->ssh_rhs_old
             && d->fer_uv && d->fer_w,
             "fesom_dyn alloc: out of memory");
}

void fesom_dyn_free(fesom_dyn *d)
{
    free(d->uv);
    free(d->uv_rhs);
    free(d->uv_rhsAB);
    free(d->w);
    free(d->w_e);
    free(d->w_i);
    free(d->cfl_z);
    free(d->uvnode);
    free(d->uvnode_rhs);
    free(d->u_b);
    free(d->v_b);
    free(d->u_c);
    free(d->v_c);
    free(d->eta_n);
    free(d->d_eta);
    free(d->ssh_rhs);
    free(d->ssh_rhs_old);
    free(d->fer_uv);
    free(d->fer_w);
    memset(d, 0, sizeof(*d));
}
