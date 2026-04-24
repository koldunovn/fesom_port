/*
 * fesom_partit.c — partition + communicator state, ASCII partition file
 * readers. Mirrors gen_modules_partitioning.F90::par_init plus the partition
 * file readers in oce_mesh.F90:280-330 (my_list) and 800-879 (com_info).
 */
#include "fesom_partit.h"

#include "fesom_types.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MPI_CHECK(call) do {                                              \
    int _ec = (call);                                                     \
    if (_ec != MPI_SUCCESS) {                                             \
        char _s[MPI_MAX_ERROR_STRING]; int _n=0;                          \
        MPI_Error_string(_ec, _s, &_n);                                   \
        FESOM_DIE("MPI error: %s", _s);                                   \
    }                                                                     \
} while (0)

/* ------------------------------------------------------------------------ *
 * ASCII helpers — Fortran's free-format read(fileID,*) accepts whitespace- *
 * and newline-separated integers indistinguishably. fscanf("%d") matches.  *
 * ------------------------------------------------------------------------ */
static int read_int(FILE *f)
{
    int v;
    if (fscanf(f, " %d", &v) != 1) {
        FESOM_DIE("fesom_partit: expected integer at offset %ld", ftell(f));
    }
    return v;
}

static void read_int_array(FILE *f, int *dst, int n)
{
    for (int i = 0; i < n; ++i) dst[i] = read_int(f);
}

/* Build dist_<npes> directory path from mesh_dir. */
static void build_dist_dir(const char *mesh_dir, int npes, char *out, size_t out_sz)
{
    size_t mlen = strlen(mesh_dir);
    int has_slash = mlen > 0 && mesh_dir[mlen - 1] == '/';
    int n = snprintf(out, out_sz, "%s%sdist_%d", mesh_dir, has_slash ? "" : "/", npes);
    FESOM_CHECK(n > 0 && (size_t)n < out_sz, "fesom_partit: dist_dir buffer too small");
}

/* Path under dist_dir: e.g. dist_8/my_list00003.out */
static void build_rank_file(const char *dist_dir, const char *base, int rank,
                            char *out, size_t out_sz)
{
    int n = snprintf(out, out_sz, "%s/%s%05d.out", dist_dir, base, rank);
    FESOM_CHECK(n > 0 && (size_t)n < out_sz, "fesom_partit: rank-file buffer too small");
}

/* ------------------------------------------------------------------------ *
 * Read rpart.out — Fortran oce_mesh.F90:269-298                            *
 *                                                                          *
 * File layout (mype==0 reads, broadcasts):                                 *
 *   line 1: npes (must match runtime npes)                                 *
 *   line 2: per-rank counts (npes integers); converted to a prefix sum     *
 *           starting at part[0]=1, so part[k+1] - part[k] = rank-k count.  *
 * ------------------------------------------------------------------------ */
static void read_rpart(fesom_partit *p, const char *dist_dir)
{
    char path[2200];
    snprintf(path, sizeof(path), "%s/rpart.out", dist_dir);

    int file_npes = -1;
    p->part = malloc((size_t)(p->npes + 1) * sizeof(int));
    FESOM_CHECK(p->part, "fesom_partit: part alloc");

    if (p->mype == 0) {
        FILE *f = fopen(path, "r");
        if (!f) FESOM_DIE("fesom_partit: cannot open %s", path);
        file_npes = read_int(f);
        if (file_npes != p->npes) {
            FESOM_DIE("fesom_partit: rpart.out npes=%d does not match runtime npes=%d",
                      file_npes, p->npes);
        }
        p->part[0] = 1;          /* Fortran 1-based starting offset */
        for (int k = 1; k <= p->npes; ++k) {
            int cnt = read_int(f);
            p->part[k] = p->part[k - 1] + cnt;
        }
        fclose(f);
    }
    MPI_CHECK(MPI_Bcast(p->part, p->npes + 1, MPI_INT, 0, p->MPI_COMM_FESOM));
}

/* ------------------------------------------------------------------------ *
 * Read my_list<NNNNN>.out — Fortran oce_mesh.F90:305-327                   *
 *                                                                          *
 * File layout (each rank reads its own file):                              *
 *   line 1:  rank (sanity check, must equal mype)                          *
 *   line 2:  myDim_nod2D                                                   *
 *   line 3:  eDim_nod2D                                                    *
 *   ...:     myList_nod2D[1..myDim+eDim]   (whitespace-separated)          *
 *   line N:  myDim_elem2D                                                  *
 *   line N+1: eDim_elem2D                                                  *
 *   line N+2: eXDim_elem2D                                                 *
 *   ...:     myList_elem2D[1..myDim+eDim+eXDim]                            *
 *   line P:  myDim_edge2D                                                  *
 *   line P+1: eDim_edge2D                                                  *
 *   ...:     myList_edge2D[1..myDim+eDim]                                  *
 * ------------------------------------------------------------------------ */
static void read_my_list(fesom_partit *p, const char *dist_dir)
{
    char path[2200];
    build_rank_file(dist_dir, "my_list", p->mype, path, sizeof(path));

    FILE *f = fopen(path, "r");
    if (!f) FESOM_DIE("fesom_partit: cannot open %s", path);

    int rank = read_int(f);
    if (rank != p->mype) {
        FESOM_DIE("fesom_partit: %s says rank=%d, expected %d", path, rank, p->mype);
    }

    p->myDim_nod2D = read_int(f);
    p->eDim_nod2D  = read_int(f);
    p->myList_nod2D = malloc((size_t)(p->myDim_nod2D + p->eDim_nod2D) * sizeof(int));
    FESOM_CHECK(p->myList_nod2D, "fesom_partit: myList_nod2D alloc");
    read_int_array(f, p->myList_nod2D, p->myDim_nod2D + p->eDim_nod2D);

    p->myDim_elem2D  = read_int(f);
    p->eDim_elem2D   = read_int(f);
    p->eXDim_elem2D  = read_int(f);
    int n_elem = p->myDim_elem2D + p->eDim_elem2D + p->eXDim_elem2D;
    p->myList_elem2D = malloc((size_t)n_elem * sizeof(int));
    FESOM_CHECK(p->myList_elem2D, "fesom_partit: myList_elem2D alloc");
    read_int_array(f, p->myList_elem2D, n_elem);

    p->myDim_edge2D = read_int(f);
    p->eDim_edge2D  = read_int(f);
    p->myList_edge2D = malloc((size_t)(p->myDim_edge2D + p->eDim_edge2D) * sizeof(int));
    FESOM_CHECK(p->myList_edge2D, "fesom_partit: myList_edge2D alloc");
    read_int_array(f, p->myList_edge2D, p->myDim_edge2D + p->eDim_edge2D);

    fclose(f);
}

/* ------------------------------------------------------------------------ *
 * Read one com_struct block from com_info — Fortran oce_mesh.F90:802-879   *
 *                                                                          *
 * rlist_size is supplied by the caller because the file does not encode it *
 * directly (it follows from the partition convention: nod2D rlist=eDim,    *
 * elem2D rlist=eDim, elem2D_full rlist=eDim+eXDim).                        *
 * ------------------------------------------------------------------------ */
static void read_com_block(FILE *f, fesom_com_struct *cs, int rlist_size,
                           const char *kind)
{
    cs->rPEnum = read_int(f);
    if (cs->rPEnum > FESOM_MAX_NEIGHBOR_PARTITIONS) {
        FESOM_DIE("fesom_partit: %s rPEnum=%d exceeds FESOM_MAX_NEIGHBOR_PARTITIONS=%d",
                  kind, cs->rPEnum, FESOM_MAX_NEIGHBOR_PARTITIONS);
    }
    read_int_array(f, cs->rPE,  cs->rPEnum);
    read_int_array(f, cs->rptr, cs->rPEnum + 1);

    cs->rlist = malloc((size_t)rlist_size * sizeof(int));
    FESOM_CHECK(cs->rlist, "fesom_partit: %s rlist alloc", kind);
    read_int_array(f, cs->rlist, rlist_size);

    cs->sPEnum = read_int(f);
    if (cs->sPEnum > FESOM_MAX_NEIGHBOR_PARTITIONS) {
        FESOM_DIE("fesom_partit: %s sPEnum=%d exceeds FESOM_MAX_NEIGHBOR_PARTITIONS=%d",
                  kind, cs->sPEnum, FESOM_MAX_NEIGHBOR_PARTITIONS);
    }
    read_int_array(f, cs->sPE,  cs->sPEnum);
    read_int_array(f, cs->sptr, cs->sPEnum + 1);

    int slist_size = cs->sptr[cs->sPEnum] - cs->sptr[0];  /* sptr is 1-based */
    cs->slist = malloc((size_t)slist_size * sizeof(int));
    FESOM_CHECK(cs->slist, "fesom_partit: %s slist alloc", kind);
    read_int_array(f, cs->slist, slist_size);
}

/* ------------------------------------------------------------------------ *
 * Read com_info<NNNNN>.out — Fortran oce_mesh.F90:802-879                  *
 *                                                                          *
 * File layout: rank header, then 3 com_struct blocks (nod2D, elem2D,       *
 * elem2D_full). Verified format from pi/dist_2/com_info00000.out (47 lines,*
 * no trailing nreq — that field only appears in the binary restart format).*
 * ------------------------------------------------------------------------ */
static void read_com_info(fesom_partit *p, const char *dist_dir)
{
    char path[2200];
    build_rank_file(dist_dir, "com_info", p->mype, path, sizeof(path));

    FILE *f = fopen(path, "r");
    if (!f) FESOM_DIE("fesom_partit: cannot open %s", path);

    int rank = read_int(f);
    if (rank != p->mype) {
        FESOM_DIE("fesom_partit: %s says rank=%d, expected %d", path, rank, p->mype);
    }

    read_com_block(f, &p->com_nod2D,        p->eDim_nod2D,                       "nod2D");
    read_com_block(f, &p->com_elem2D,       p->eDim_elem2D,                      "elem2D");
    read_com_block(f, &p->com_elem2D_full,  p->eDim_elem2D + p->eXDim_elem2D,    "elem2D_full");

    fclose(f);
}

/* ------------------------------------------------------------------------ *
 * Public API                                                                *
 * ------------------------------------------------------------------------ */
void fesom_partit_init(fesom_partit *p,
                       const char   *mesh_dir,
                       int           argc,
                       char        **argv)
{
    memset(p, 0, sizeof(*p));

    int provided = 0;
    /* Fortran uses MPI_THREAD_SINGLE-equivalent (MPI_Init); we follow suit. */
    int already_initialised = 0;
    MPI_CHECK(MPI_Initialized(&already_initialised));
    if (!already_initialised) {
        MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided));
    }
    p->MPI_COMM_FESOM = MPI_COMM_WORLD;
    MPI_CHECK(MPI_Comm_size(p->MPI_COMM_FESOM, &p->npes));
    MPI_CHECK(MPI_Comm_rank(p->MPI_COMM_FESOM, &p->mype));

    if (p->npes == 1) {
        /* Synthesised partit — defer setting myDim/eDim until mesh dims known.
         * fesom_partit_set_global_counts_serial() fills them. */
        if (p->mype == 0) {
            printf("[fesom_partit] running on 1 rank (synthesised partit, no file reads)\n");
        }
        return;
    }

    char dist_dir[2048];
    build_dist_dir(mesh_dir, p->npes, dist_dir, sizeof(dist_dir));

    /* Quick existence check on rank 0. */
    if (p->mype == 0) {
        char probe[2200];
        snprintf(probe, sizeof(probe), "%s/rpart.out", dist_dir);
        FILE *f = fopen(probe, "r");
        if (!f) {
            FESOM_DIE("fesom_partit: %s missing — no partition for npes=%d under %s",
                      probe, p->npes, mesh_dir);
        }
        fclose(f);
        printf("[fesom_partit] running on %d ranks, partition: %s\n",
               p->npes, dist_dir);
    }

    read_rpart   (p, dist_dir);
    read_my_list (p, dist_dir);
    read_com_info(p, dist_dir);

    /* Sanity: per-rank Σ counts. */
    int total_my_n = 0, total_my_e = 0, total_my_ed = 0;
    MPI_CHECK(MPI_Reduce(&p->myDim_nod2D,  &total_my_n,  1, MPI_INT, MPI_SUM, 0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Reduce(&p->myDim_elem2D, &total_my_e,  1, MPI_INT, MPI_SUM, 0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Reduce(&p->myDim_edge2D, &total_my_ed, 1, MPI_INT, MPI_SUM, 0, p->MPI_COMM_FESOM));
    if (p->mype == 0) {
        printf("[fesom_partit] sum across ranks: myDim_nod2D=%d myDim_elem2D=%d myDim_edge2D=%d\n",
               total_my_n, total_my_e, total_my_ed);
    }
    /* Per-rank breakdown. */
    printf("[fesom_partit]   rank %3d: myDim_nod2D=%d eDim_nod2D=%d "
           "myDim_elem2D=%d eDim_elem2D=%d eXDim_elem2D=%d "
           "myDim_edge2D=%d eDim_edge2D=%d\n",
           p->mype,
           p->myDim_nod2D, p->eDim_nod2D,
           p->myDim_elem2D, p->eDim_elem2D, p->eXDim_elem2D,
           p->myDim_edge2D, p->eDim_edge2D);
}

void fesom_partit_set_global_counts_serial(fesom_partit *p,
                                           int nod2D, int elem2D, int edge2D)
{
    if (p->npes != 1) return;
    p->myDim_nod2D  = nod2D;   p->eDim_nod2D   = 0;
    p->myDim_elem2D = elem2D;  p->eDim_elem2D  = 0;  p->eXDim_elem2D = 0;
    p->myDim_edge2D = edge2D;  p->eDim_edge2D  = 0;

    /* Synthesise myList_* as identity (1-based global IDs == local indices+1). */
    p->myList_nod2D = malloc((size_t)nod2D * sizeof(int));
    FESOM_CHECK(p->myList_nod2D, "fesom_partit: serial myList_nod2D alloc");
    for (int i = 0; i < nod2D; ++i) p->myList_nod2D[i] = i + 1;

    p->myList_elem2D = malloc((size_t)elem2D * sizeof(int));
    FESOM_CHECK(p->myList_elem2D, "fesom_partit: serial myList_elem2D alloc");
    for (int i = 0; i < elem2D; ++i) p->myList_elem2D[i] = i + 1;

    p->myList_edge2D = malloc((size_t)edge2D * sizeof(int));
    FESOM_CHECK(p->myList_edge2D, "fesom_partit: serial myList_edge2D alloc");
    for (int i = 0; i < edge2D; ++i) p->myList_edge2D[i] = i + 1;

    p->part = malloc(2 * sizeof(int));
    FESOM_CHECK(p->part, "fesom_partit: serial part alloc");
    p->part[0] = 1;
    p->part[1] = nod2D + 1;

    /* Empty com_struct triples — no halo, no neighbours. */
    /* memset on init already zeroed rPEnum/sPEnum/rlist/slist. */
}

static void free_com_struct(fesom_com_struct *cs)
{
    free(cs->rlist); cs->rlist = NULL;
    free(cs->slist); cs->slist = NULL;
}

void fesom_partit_finalize(fesom_partit *p)
{
    free_com_struct(&p->com_nod2D);
    free_com_struct(&p->com_elem2D);
    free_com_struct(&p->com_elem2D_full);
    free(p->part);             p->part = NULL;
    free(p->myList_nod2D);     p->myList_nod2D = NULL;
    free(p->myList_elem2D);    p->myList_elem2D = NULL;
    free(p->myList_edge2D);    p->myList_edge2D = NULL;

    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Finalize();
}
