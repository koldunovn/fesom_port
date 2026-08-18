/*
 * fesom_io_gather.c — the gather/scatter plan shared by the snapshot writer,
 * the time-mean streams and the restart layer.
 *
 * Split out of fesom_io.c so that code needing only the collective (the
 * restart reader and its unit test) does not have to link the whole output
 * subsystem.
 */
#include "fesom_io_gather.h"

#include "fesom_mesh.h"
#include "fesom_partit.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */

void gather_plan_init(gather_plan *gp,
                      const struct fesom_mesh *mesh,
                      struct fesom_partit *partit)
{
    memset(gp, 0, sizeof(*gp));
    gp->npes    = partit->npes;
    gp->mype    = partit->mype;
    gp->myDim_n = mesh->myDim_nod2D;
    gp->myDim_e = mesh->myDim_elem2D;
    if (gp->mype == 0) {
        gp->all_n   = malloc((size_t)gp->npes * sizeof(int));
        gp->all_e   = malloc((size_t)gp->npes * sizeof(int));
        gp->displ_n = malloc((size_t)gp->npes * sizeof(int));
        gp->displ_e = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(gp->all_n && gp->all_e && gp->displ_n && gp->displ_e,
                    "gather_plan: oom");
    }
    MPI_Gather(&gp->myDim_n, 1, MPI_INT, gp->all_n, 1, MPI_INT, 0, partit->MPI_COMM_FESOM);
    MPI_Gather(&gp->myDim_e, 1, MPI_INT, gp->all_e, 1, MPI_INT, 0, partit->MPI_COMM_FESOM);
    if (gp->mype == 0) {
        gp->displ_n[0] = 0;
        gp->displ_e[0] = 0;
        for (int r = 1; r < gp->npes; ++r) {
            gp->displ_n[r] = gp->displ_n[r-1] + gp->all_n[r-1];
            gp->displ_e[r] = gp->displ_e[r-1] + gp->all_e[r-1];
        }
        gp->total_n = gp->displ_n[gp->npes-1] + gp->all_n[gp->npes-1];
        gp->total_e = gp->displ_e[gp->npes-1] + gp->all_e[gp->npes-1];
        /* Nodes are partitioned exclusively (each global node has one owner)
         * so Σ myDim_n = nod2D. Elements may be REPLICATED across rank
         * boundaries (Fortran convention; sum can exceed elem2D). When
         * unpacking, duplicate writes for the same global element id are
         * harmless because each contributing rank produces the same value
         * after halo exchange. */
        FESOM_CHECK(gp->total_n == mesh->nod2D,
                    "gather_plan: Σmyn=%d != nod2D=%d", gp->total_n, mesh->nod2D);
        FESOM_CHECK(gp->total_e >= mesh->elem2D,
                    "gather_plan: Σmye=%d < elem2D=%d (cannot cover global)",
                    gp->total_e, mesh->elem2D);
        gp->gathered_myList_n = malloc((size_t)gp->total_n * sizeof(int));
        gp->gathered_myList_e = malloc((size_t)gp->total_e * sizeof(int));
        FESOM_CHECK(gp->gathered_myList_n && gp->gathered_myList_e, "gather_plan: oom");
    }
    MPI_Gatherv(partit->myList_nod2D, gp->myDim_n, MPI_INT,
                gp->gathered_myList_n, gp->all_n, gp->displ_n, MPI_INT,
                0, partit->MPI_COMM_FESOM);
    MPI_Gatherv(partit->myList_elem2D, gp->myDim_e, MPI_INT,
                gp->gathered_myList_e, gp->all_e, gp->displ_e, MPI_INT,
                0, partit->MPI_COMM_FESOM);
}

void gather_plan_free(gather_plan *gp)
{
    free(gp->all_n);   free(gp->all_e);
    free(gp->displ_n); free(gp->displ_e);
    free(gp->gathered_myList_n); free(gp->gathered_myList_e);
    memset(gp, 0, sizeof(*gp));
}

/* Gather an [my_n × stride] node-indexed buffer to rank 0's
 * [total_n × stride] global buffer indexed by global node id (0-based). */
void gather_node(const real_t *local, int stride,
                 const gather_plan *gp,
                 real_t *global,
                 MPI_Comm comm)
{
    real_t *recv = NULL;
    int *counts = NULL, *displs = NULL;
    if (gp->mype == 0) {
        recv   = malloc((size_t)gp->total_n * (size_t)stride * sizeof(real_t));
        counts = malloc((size_t)gp->npes * sizeof(int));
        displs = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(recv && counts && displs, "gather_node: oom");
        for (int r = 0; r < gp->npes; ++r) {
            counts[r] = gp->all_n[r] * stride;
            displs[r] = gp->displ_n[r] * stride;
        }
    }
    MPI_Gatherv((void *)local, gp->myDim_n * stride, MPI_DOUBLE,
                recv, counts, displs, MPI_DOUBLE, 0, comm);
    if (gp->mype == 0) {
        for (int r = 0; r < gp->npes; ++r) {
            for (int i = 0; i < gp->all_n[r]; ++i) {
                int gid = gp->gathered_myList_n[gp->displ_n[r] + i] - 1;
                memcpy(global + (size_t)gid * stride,
                       recv + ((size_t)gp->displ_n[r] + i) * stride,
                       (size_t)stride * sizeof(real_t));
            }
        }
        free(recv); free(counts); free(displs);
    }
}

/* Same for element-indexed [my_e × stride] → [total_e × stride]. */
void gather_elem(const real_t *local, int stride,
                 const gather_plan *gp,
                 real_t *global,
                 MPI_Comm comm)
{
    real_t *recv = NULL;
    int *counts = NULL, *displs = NULL;
    if (gp->mype == 0) {
        recv   = malloc((size_t)gp->total_e * (size_t)stride * sizeof(real_t));
        counts = malloc((size_t)gp->npes * sizeof(int));
        displs = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(recv && counts && displs, "gather_elem: oom");
        for (int r = 0; r < gp->npes; ++r) {
            counts[r] = gp->all_e[r] * stride;
            displs[r] = gp->displ_e[r] * stride;
        }
    }
    MPI_Gatherv((void *)local, gp->myDim_e * stride, MPI_DOUBLE,
                recv, counts, displs, MPI_DOUBLE, 0, comm);
    if (gp->mype == 0) {
        for (int r = 0; r < gp->npes; ++r) {
            for (int i = 0; i < gp->all_e[r]; ++i) {
                int gid = gp->gathered_myList_e[gp->displ_e[r] + i] - 1;
                memcpy(global + (size_t)gid * stride,
                       recv + ((size_t)gp->displ_e[r] + i) * stride,
                       (size_t)stride * sizeof(real_t));
            }
        }
        free(recv); free(counts); free(displs);
    }
}

/* ------------------------------------------------------------------ *
 * scatter_node / scatter_elem — the inverse of the two gathers above. *
 *                                                                     *
 * Rank 0 permutes its global buffer into rank order using the same    *
 * gathered_myList_* it uses on the way out, then MPI_Scatterv sends   *
 * each rank its interior block. Because the permutation is the        *
 * gather's read side run backwards, a field written on N ranks and    *
 * read on M ranks lands on the same global node ids in both cases.    *
 * ------------------------------------------------------------------ */
void scatter_node(const real_t *global, int stride,
                  const gather_plan *gp,
                  real_t *local,
                  MPI_Comm comm)
{
    real_t *send = NULL;
    int *counts = NULL, *displs = NULL;
    if (gp->mype == 0) {
        send   = malloc((size_t)gp->total_n * (size_t)stride * sizeof(real_t));
        counts = malloc((size_t)gp->npes * sizeof(int));
        displs = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(send && counts && displs, "scatter_node: oom");
        for (int r = 0; r < gp->npes; ++r) {
            counts[r] = gp->all_n[r] * stride;
            displs[r] = gp->displ_n[r] * stride;
            for (int i = 0; i < gp->all_n[r]; ++i) {
                int gid = gp->gathered_myList_n[gp->displ_n[r] + i] - 1;
                memcpy(send + ((size_t)gp->displ_n[r] + i) * stride,
                       global + (size_t)gid * stride,
                       (size_t)stride * sizeof(real_t));
            }
        }
    }
    MPI_Scatterv(send, counts, displs, MPI_DOUBLE,
                 local, gp->myDim_n * stride, MPI_DOUBLE, 0, comm);
    if (gp->mype == 0) { free(send); free(counts); free(displs); }
}

void scatter_elem(const real_t *global, int stride,
                  const gather_plan *gp,
                  real_t *local,
                  MPI_Comm comm)
{
    real_t *send = NULL;
    int *counts = NULL, *displs = NULL;
    if (gp->mype == 0) {
        send   = malloc((size_t)gp->total_e * (size_t)stride * sizeof(real_t));
        counts = malloc((size_t)gp->npes * sizeof(int));
        displs = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(send && counts && displs, "scatter_elem: oom");
        for (int r = 0; r < gp->npes; ++r) {
            counts[r] = gp->all_e[r] * stride;
            displs[r] = gp->displ_e[r] * stride;
            for (int i = 0; i < gp->all_e[r]; ++i) {
                int gid = gp->gathered_myList_e[gp->displ_e[r] + i] - 1;
                memcpy(send + ((size_t)gp->displ_e[r] + i) * stride,
                       global + (size_t)gid * stride,
                       (size_t)stride * sizeof(real_t));
            }
        }
    }
    MPI_Scatterv(send, counts, displs, MPI_DOUBLE,
                 local, gp->myDim_e * stride, MPI_DOUBLE, 0, comm);
    if (gp->mype == 0) { free(send); free(counts); free(displs); }
}

