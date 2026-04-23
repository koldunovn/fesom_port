#include "fesom_dyn.h"
#include "fesom_mesh.h"
#include "fesom_constants.h"

#include <stdlib.h>
#include <string.h>

void fesom_dyn_alloc(fesom_dyn *d, const struct fesom_mesh *mesh)
{
    memset(d, 0, sizeof(*d));
    d->AB_order = FESOM_PHASE1_AB_ORDER;

    size_t e_nl_2 = (size_t)mesh->elem2D * (size_t)mesh->nl * 2;
    size_t n_nl   = (size_t)mesh->nod2D  * (size_t)mesh->nl;
    size_t n      = (size_t)mesh->nod2D;

    d->uv          = calloc(e_nl_2, sizeof(real_t));
    d->uv_rhs      = calloc(e_nl_2, sizeof(real_t));
    d->uv_rhsAB    = calloc(e_nl_2, sizeof(real_t));
    d->w           = calloc(n_nl,   sizeof(real_t));
    d->w_i         = calloc(n_nl,   sizeof(real_t));
    d->uvnode      = calloc(n_nl * 2, sizeof(real_t));
    d->eta_n       = calloc(n,      sizeof(real_t));
    d->d_eta       = calloc(n,      sizeof(real_t));
    d->ssh_rhs     = calloc(n,      sizeof(real_t));
    d->ssh_rhs_old = calloc(n,      sizeof(real_t));
    FESOM_CHECK(d->uv && d->uv_rhs && d->uv_rhsAB && d->w && d->w_i && d->uvnode
             && d->eta_n && d->d_eta && d->ssh_rhs && d->ssh_rhs_old,
             "fesom_dyn alloc: out of memory");
}

void fesom_dyn_free(fesom_dyn *d)
{
    free(d->uv);
    free(d->uv_rhs);
    free(d->uv_rhsAB);
    free(d->w);
    free(d->w_i);
    free(d->uvnode);
    free(d->eta_n);
    free(d->d_eta);
    free(d->ssh_rhs);
    free(d->ssh_rhs_old);
    memset(d, 0, sizeof(*d));
}
