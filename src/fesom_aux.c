#include "fesom_aux.h"
#include "fesom_mesh.h"

#include <stdlib.h>
#include <string.h>

void fesom_aux_alloc(fesom_aux *a, const struct fesom_mesh *mesh)
{
    memset(a, 0, sizeof(*a));
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int E = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    size_t n_nl = (size_t)N * (size_t)mesh->nl;
    size_t e_nl = (size_t)E * (size_t)mesh->nl;

    a->density_m_rho0 = calloc(n_nl, sizeof(real_t));
    a->hpressure      = calloc(n_nl, sizeof(real_t));
    a->bvfreq         = calloc(n_nl, sizeof(real_t));
    a->sw_alpha       = calloc(n_nl, sizeof(real_t));
    a->sw_beta        = calloc(n_nl, sizeof(real_t));
    a->Kv             = calloc(n_nl, sizeof(real_t));
    a->dbsfc          = calloc(n_nl, sizeof(real_t));

    a->Av             = calloc(e_nl, sizeof(real_t));
    a->pgf_x          = calloc(e_nl, sizeof(real_t));
    a->pgf_y          = calloc(e_nl, sizeof(real_t));

    /* Phase G2a: MLD1_ind (filled by fesom_pressure_bv). */
    a->MLD1_ind       = calloc((size_t)N, sizeof(int));

    FESOM_CHECK(a->density_m_rho0 && a->hpressure && a->bvfreq
             && a->sw_alpha && a->sw_beta && a->Kv && a->dbsfc
             && a->Av && a->pgf_x && a->pgf_y
             && a->MLD1_ind,
             "fesom_aux alloc: out of memory");
}

void fesom_aux_free(fesom_aux *a)
{
    free(a->density_m_rho0);
    free(a->hpressure);
    free(a->bvfreq);
    free(a->sw_alpha);
    free(a->sw_beta);
    free(a->Kv);
    free(a->dbsfc);
    free(a->Av);
    free(a->pgf_x);
    free(a->pgf_y);
    free(a->MLD1_ind);
    memset(a, 0, sizeof(*a));
}
