#include "fesom_forcing.h"
#include "fesom_mesh.h"

#include <stdlib.h>
#include <string.h>

void fesom_forcing_alloc(fesom_forcing *f, const struct fesom_mesh *mesh)
{
    memset(f, 0, sizeof(*f));
    size_t n  = (size_t)mesh->nod2D;
    size_t e2 = (size_t)mesh->elem2D * 2;

    f->heat_flux        = calloc(n,      sizeof(real_t));
    f->water_flux       = calloc(n,      sizeof(real_t));
    f->stress_node_surf = calloc(n * 2,  sizeof(real_t));
    f->stress_surf      = calloc(e2,     sizeof(real_t));
    FESOM_CHECK(f->heat_flux && f->water_flux
             && f->stress_node_surf && f->stress_surf,
             "fesom_forcing alloc: out of memory");
}

void fesom_forcing_free(fesom_forcing *f)
{
    free(f->heat_flux);
    free(f->water_flux);
    free(f->stress_node_surf);
    free(f->stress_surf);
    memset(f, 0, sizeof(*f));
}
