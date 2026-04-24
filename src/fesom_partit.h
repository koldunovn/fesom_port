#ifndef FESOM_PARTIT_H
#define FESOM_PARTIT_H

#include <mpi.h>

#include "fesom_types.h"

/*
 * Partition / communicator state. Mirror of Fortran T_PARTIT (MOD_PARTIT.F90)
 * and com_struct, restricted to fields actually used in the model loop.
 *
 * Domain partitioning is precomputed by the Fortran mesh-partition tool and
 * lives on disk under <mesh_dir>/dist_<NPES>/. Files (ASCII):
 *
 *   rpart.out                  global partition vector (which rank owns each
 *                              global node, by node ID range)
 *   my_list<NNNNN>.out         per-rank: counts + global-id lists for nod2D,
 *                              elem2D (with extended halo eXDim), edge2D
 *   com_info<NNNNN>.out        per-rank: 3 com_struct blocks for nod2D,
 *                              elem2D (small halo), elem2D_full (extended halo)
 *
 * com_struct file layout (verified from pi/dist_2/com_info00000.out):
 *   rPEnum               - number of receiver-partner ranks
 *   rPE [1..rPEnum]      - their global rank ids
 *   rptr[1..rPEnum+1]    - cumulative offsets into rlist (1-based)
 *   rlist[1..eDim]       - LOCAL halo indices (1-based; range
 *                          (myDim..myDim+eDim])
 *   sPEnum               - sender-partner count
 *   sPE [1..sPEnum]      - their global rank ids
 *   sptr[1..sPEnum+1]    - cumulative offsets into slist (1-based)
 *   slist[1..sptr[end]-1]- LOCAL interior indices (1-based; range
 *                          [1..myDim]) — entries this rank sends to halos of
 *                          other ranks
 *
 * Index convention: file values are 1-based Fortran. We store them ALSO
 * 1-based in C (no offset shift) so the indexing arithmetic mirrors Fortran
 * line-for-line; halo exchange code subtracts 1 when reading/writing field[].
 */

#define FESOM_MAX_NEIGHBOR_PARTITIONS 32  /* Fortran MOD_PARTIT.F90:15 */

typedef struct fesom_com_struct {
    int  rPEnum;
    int  rPE [FESOM_MAX_NEIGHBOR_PARTITIONS];
    int  rptr[FESOM_MAX_NEIGHBOR_PARTITIONS + 1];
    int *rlist;                /* dynamic, size eDim */

    int  sPEnum;
    int  sPE [FESOM_MAX_NEIGHBOR_PARTITIONS];
    int  sptr[FESOM_MAX_NEIGHBOR_PARTITIONS + 1];
    int *slist;                /* dynamic, size sptr[sPEnum]-1 */
} fesom_com_struct;

typedef struct fesom_partit {
    int       npes;
    int       mype;
    MPI_Comm  MPI_COMM_FESOM;

    /* Three communication structures (Fortran T_PARTIT). */
    fesom_com_struct com_nod2D;
    fesom_com_struct com_elem2D;
    fesom_com_struct com_elem2D_full;

    /* Global partition vector (npes+1 entries; Fortran 1-based prefix-sum
     * over node IDs). part[k] for k in [1..npes+1]; rank k owns global node
     * ids in [part[k] .. part[k+1]-1]. */
    int *part;

    /* Per-rank counts and global-id lists (1-based global IDs).
     *   myDim_*  - this rank's interior count
     *   eDim_*   - small halo count (read from neighbours)
     *   eXDim_*  - extended halo count (additional, only for full-halo elem)
     *   myList_* - global IDs, [0..myDim+eDim[+eXDim])
     */
    int  myDim_nod2D, eDim_nod2D;
    int *myList_nod2D;

    int  myDim_elem2D, eDim_elem2D, eXDim_elem2D;
    int *myList_elem2D;

    int  myDim_edge2D, eDim_edge2D;
    int *myList_edge2D;
} fesom_partit;

/* Initialise MPI + read partition files for the supplied mesh directory.
 * If npes == 1, synthesise a trivial partit (myDim = global counts, eDim/eXDim
 * = 0, no neighbours) and skip file reads — this preserves bit-identical
 * behaviour with the existing serial run.
 * Aborts with a clear message if dist_<npes>/ is absent and npes > 1. */
void fesom_partit_init(fesom_partit *p,
                       const char   *mesh_dir,
                       int           argc,
                       char        **argv);

void fesom_partit_finalize(fesom_partit *p);

/* Set global counts so the synthesised npes==1 partit knows the mesh size.
 * Called once after fesom_mesh_read on rank 0 before any halo exchange.
 * For multi-rank, this is a no-op (partit already populated). */
void fesom_partit_set_global_counts_serial(fesom_partit *p,
                                           int nod2D,
                                           int elem2D,
                                           int edge2D);

#endif /* FESOM_PARTIT_H */
