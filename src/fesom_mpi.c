#include "fesom_mpi.h"

void fesom_mpi_init(fesom_mpi *p, int argc, char **argv)
{
    (void)argc; (void)argv;
    p->rank = 0;
    p->nranks = 1;
    p->dummy_comm = 0;
}

void fesom_mpi_finalize(fesom_mpi *p)
{
    (void)p;
}

void fesom_exchange_nod(real_t *field, int nlev, const struct fesom_mesh *mesh,
                        const fesom_mpi *p)
{
    (void)field; (void)nlev; (void)mesh; (void)p;
    /* Serial mode: no halo nodes, nothing to exchange. */
}
