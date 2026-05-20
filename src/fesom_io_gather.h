/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 * See fesom_io.h banner for context.
 *
 * fesom_io_gather.h - the gather_plan abstraction shared between
 * fesom_io.c (snapshot writer) and fesom_io_stream.c (time-mean
 * streams). Lives in its own header so the two modules don't have
 * to round-trip through fesom_io.h.
 *
 * Internal struct fields are documented but not for direct manipulation
 * by callers; use the four public functions only.
 *
 * Important MPI semantics:
 *   - Nodes are partitioned exclusively across ranks → Σ myDim_n == nod2D.
 *   - Elements may be REPLICATED across rank boundaries (Σ myDim_e ≥ elem2D).
 *     Duplicate writes for the same global element id are harmless because
 *     each contributing rank produces the same value after halo exchange.
 */
#ifndef FESOM_IO_GATHER_H
#define FESOM_IO_GATHER_H

#include "fesom_types.h"

#include <mpi.h>

struct fesom_mesh;
struct fesom_partit;

typedef struct gather_plan {
    int  npes;
    int  mype;
    int  myDim_n;            /* this rank's interior node count        */
    int  myDim_e;            /* this rank's interior element count     */
    int *all_n;              /* [npes] per-rank myDim_nod2D (rank 0)   */
    int *all_e;              /* [npes] per-rank myDim_elem2D (rank 0)  */
    int *displ_n;            /* [npes] node-count displacements        */
    int *displ_e;            /* [npes] element-count displacements     */
    int  total_n;            /* sum of all_n  == mesh->nod2D           */
    int  total_e;            /* sum of all_e  >= mesh->elem2D          */
    int *gathered_myList_n;  /* [total_n] rank-0 view of myList_nod2D  */
    int *gathered_myList_e;  /* [total_e] rank-0 view of myList_elem2D */
} gather_plan;

void gather_plan_init(gather_plan *gp,
                      const struct fesom_mesh *mesh,
                      struct fesom_partit *partit);
void gather_plan_free(gather_plan *gp);

/* Gather a [myDim_n × stride] node-indexed buffer to rank 0's
 * [total_n × stride] global buffer indexed by global node id (0-based).
 * No-op on non-zero ranks (rank 0 is the only writer in this IO model). */
void gather_node(const real_t *local, int stride,
                 const gather_plan *gp,
                 real_t *global,
                 MPI_Comm comm);

/* Same for element-indexed [myDim_e × stride] → [total_e × stride]. */
void gather_elem(const real_t *local, int stride,
                 const gather_plan *gp,
                 real_t *global,
                 MPI_Comm comm);

#endif /* FESOM_IO_GATHER_H */
