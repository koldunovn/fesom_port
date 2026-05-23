#include "fesom_forcing.h"
#include "fesom_mesh.h"

#include <stdlib.h>
#include <string.h>

void fesom_forcing_alloc(fesom_forcing *f, const struct fesom_mesh *mesh)
{
    memset(f, 0, sizeof(*f));
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int E = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    size_t n  = (size_t)N;
    size_t e2 = (size_t)E * 2;

    f->heat_flux        = calloc(n,      sizeof(real_t));
    f->water_flux       = calloc(n,      sizeof(real_t));
    f->stress_node_surf = calloc(n * 2,  sizeof(real_t));
    f->stress_surf      = calloc(e2,     sizeof(real_t));
    f->runoff           = calloc(n,      sizeof(real_t));
    f->Ssurf            = calloc(n,      sizeof(real_t));
    f->virtual_salt     = calloc(n,      sizeof(real_t));
    f->relax_salt       = calloc(n,      sizeof(real_t));
    f->Ch_atm_oce       = calloc(n,      sizeof(real_t));
    f->Ce_atm_oce       = calloc(n,      sizeof(real_t));
    f->chl              = calloc(n,             sizeof(real_t));
    f->sw_3d            = calloc(n * (size_t)mesh->nl, sizeof(real_t));
    FESOM_CHECK(f->heat_flux && f->water_flux
             && f->stress_node_surf && f->stress_surf
             && f->runoff && f->Ssurf
             && f->virtual_salt && f->relax_salt
             && f->Ch_atm_oce && f->Ce_atm_oce
             && f->chl && f->sw_3d,
             "fesom_forcing alloc: out of memory");
}

void fesom_forcing_free(fesom_forcing *f)
{
    free(f->heat_flux);
    free(f->water_flux);
    free(f->stress_node_surf);
    free(f->stress_surf);
    free(f->runoff);
    free(f->Ssurf);
    free(f->virtual_salt);
    free(f->relax_salt);
    free(f->Ch_atm_oce);
    free(f->Ce_atm_oce);
    free(f->chl);
    free(f->sw_3d);
    memset(f, 0, sizeof(*f));
}
