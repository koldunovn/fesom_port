#include "fesom_tracers.h"
#include "fesom_mesh.h"

#include <stdlib.h>
#include <string.h>

void fesom_tracers_alloc(fesom_tracers *t, const struct fesom_mesh *mesh)
{
    memset(t, 0, sizeof(*t));
    t->num_tracers = FESOM_NUM_TRACERS;

    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    size_t n_nl = (size_t)N * (size_t)mesh->nl;
    for (int k = 0; k < FESOM_NUM_TRACERS; ++k) {
        t->data[k].values    = calloc(n_nl, sizeof(real_t));
        t->data[k].valuesAB  = calloc(n_nl, sizeof(real_t));
        t->data[k].valuesold = calloc(n_nl, sizeof(real_t));
        FESOM_CHECK(t->data[k].values && t->data[k].valuesAB && t->data[k].valuesold,
                    "tracer %d alloc: out of memory", k);
    }
    t->del_ttf = calloc(n_nl, sizeof(real_t));
    FESOM_CHECK(t->del_ttf, "del_ttf alloc: out of memory");
}

void fesom_tracers_free(fesom_tracers *t)
{
    for (int k = 0; k < FESOM_NUM_TRACERS; ++k) {
        free(t->data[k].values);
        free(t->data[k].valuesAB);
        free(t->data[k].valuesold);
    }
    free(t->del_ttf);
    memset(t, 0, sizeof(*t));
}
