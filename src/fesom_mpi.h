#ifndef FESOM_MPI_H
#define FESOM_MPI_H

/*
 * Phase 4: MPI is a hard dependency (matches Fortran). The real partition
 * state lives in fesom_partit.h; this file is kept as a thin compatibility
 * facade so existing call sites (`fesom_mpi`, `fesom_mpi_init`, etc.) still
 * resolve. New code should include and use fesom_partit.h directly.
 */

#include "fesom_partit.h"

typedef fesom_partit fesom_mpi;

/* Replaces the old serial `fesom_mpi_init(p, argc, argv)`. The mesh dir is
 * needed to find dist_<NPES>/. */
static inline void fesom_mpi_init(fesom_mpi *p, const char *mesh_dir,
                                  int argc, char **argv)
{
    fesom_partit_init(p, mesh_dir, argc, argv);
}

static inline void fesom_mpi_finalize(fesom_mpi *p)
{
    fesom_partit_finalize(p);
}

#endif /* FESOM_MPI_H */
