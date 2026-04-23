#ifndef FESOM_MPI_H
#define FESOM_MPI_H

#include "fesom_types.h"

struct fesom_mesh;  /* forward decl; full type in fesom_mesh.h */

/*
 * Minimal MPI shim. Phase 1 operates as MPI-with-1-rank: rank=0,
 * nranks=1, fesom_exchange_nod() is a no-op. The struct shape mirrors
 * the Fortran t_partit so that switching to real MPI later changes
 * function bodies, not call sites.
 */
typedef struct fesom_mpi {
    int rank;          /* mype  */
    int nranks;        /* npes  */
    int dummy_comm;    /* placeholder for MPI_Comm */
} fesom_mpi;

void fesom_mpi_init(fesom_mpi *p, int argc, char **argv);
void fesom_mpi_finalize(fesom_mpi *p);

/* Halo exchange of a 3D node-column field. Serial: no-op. */
void fesom_exchange_nod(real_t *field, int nlev, const struct fesom_mesh *mesh,
                        const fesom_mpi *p);

#endif /* FESOM_MPI_H */
