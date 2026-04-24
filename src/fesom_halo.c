/*
 * fesom_halo.c — generic halo exchange.
 *
 * Mirror of gen_comm.F90's exchange_nod / exchange_elem family, collapsed
 * into one parameterised path. Pack/unpack with explicit MPI buffers; the
 * Fortran builds MPI_TYPE_INDEXED over detected contiguous segments of
 * rlist/slist (gen_modules_partitioning.F90:237-281), but for C the simpler
 * pack-unpack loop is correct and equivalent.
 */
#include "fesom_halo.h"

#include "fesom_partit.h"
#include "fesom_types.h"

#include <math.h>
#include <stddef.h>
#include <mpi.h>
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
 * Per-com-struct scratch — one set per kind. Allocated lazily on first use,*
 * grown on demand. Indexed by enum value.                                  *
 * ------------------------------------------------------------------------ */
typedef struct halo_scratch {
    real_t      *send_buf;
    real_t      *recv_buf;
    size_t       send_cap;       /* in real_t entries */
    size_t       recv_cap;
    MPI_Request *reqs;
    int          reqs_cap;
} halo_scratch;

/* One scratch slot per fesom_halo_kind (5 entries). */
static halo_scratch g_scratch[5];

static void scratch_grow_buf(real_t **buf, size_t *cap, size_t need)
{
    if (need <= *cap) return;
    free(*buf);
    *buf = malloc(need * sizeof(real_t));
    FESOM_CHECK(*buf, "fesom_halo: scratch buffer alloc (need=%zu real_t)", need);
    *cap = need;
}

static void scratch_grow_reqs(halo_scratch *s, int need)
{
    if (need <= s->reqs_cap) return;
    free(s->reqs);
    s->reqs = malloc((size_t)need * sizeof(MPI_Request));
    FESOM_CHECK(s->reqs, "fesom_halo: reqs alloc");
    s->reqs_cap = need;
}

void fesom_halo_free_buffers(void)
{
    for (int i = 0; i < 5; ++i) {
        free(g_scratch[i].send_buf); g_scratch[i].send_buf = NULL;
        free(g_scratch[i].recv_buf); g_scratch[i].recv_buf = NULL;
        free(g_scratch[i].reqs);     g_scratch[i].reqs = NULL;
        g_scratch[i].send_cap = 0;
        g_scratch[i].recv_cap = 0;
        g_scratch[i].reqs_cap = 0;
    }
}

/* Resolve the right com_struct for a given kind. */
static fesom_com_struct *get_com(fesom_partit *p, fesom_halo_kind kind)
{
    switch (kind) {
    case FESOM_HALO_NOD2D:
    case FESOM_HALO_NOD3D:        return &p->com_nod2D;
    case FESOM_HALO_ELEM2D:
    case FESOM_HALO_ELEM3D:       return &p->com_elem2D;
    case FESOM_HALO_ELEM2D_FULL:  return &p->com_elem2D_full;
    default:
        FESOM_DIE("fesom_halo: unknown kind %d", (int)kind);
    }
}

/* MPI tag base — distinct per-kind so concurrent exchanges of different
 * kinds (we don't do this today, but cheap insurance) can't cross-talk. */
static int kind_tag(fesom_halo_kind kind) { return 2000 + (int)kind; }

/* ------------------------------------------------------------------------ *
 * Optional debug: nuke halo entries with NaN before exchange. Any kernel   *
 * that reads stale halo will produce NaN visible in the next snapshot.     *
 * Compile with -DFESOM_HALO_NUKE to enable.                                *
 * ------------------------------------------------------------------------ */
#ifdef FESOM_HALO_NUKE
static void halo_nuke(real_t *field, fesom_com_struct *cs,
                      int n_levels, int n_components)
{
    int stride = n_levels * n_components;
    for (int k = 0; k < cs->rPEnum; ++k) {
        int beg = cs->rptr[k] - 1;          /* 1-based → 0-based */
        int end = cs->rptr[k + 1] - 1;
        for (int j = beg; j < end; ++j) {
            int local_idx = cs->rlist[j] - 1;   /* 1-based → 0-based */
            real_t *p = field + (size_t)local_idx * stride;
            for (int s = 0; s < stride; ++s) p[s] = (real_t)NAN;
        }
    }
}
#else
static void halo_nuke(real_t *field, fesom_com_struct *cs,
                      int n_levels, int n_components)
{
    (void)field; (void)cs; (void)n_levels; (void)n_components;
}
#endif

/* ------------------------------------------------------------------------ *
 * The exchange itself.                                                      *
 *                                                                          *
 * Steps (mirror gen_comm.F90:exchange_nod3D pattern):                      *
 *   1. (debug) NaN-fill halo entries to catch missing exchanges.           *
 *   2. Allocate send + recv buffers sized to total payload bytes.          *
 *   3. Post Irecv for each receiver-PE into a contiguous slice of recv.    *
 *   4. Pack send buffer per sender-PE from slist[] indices.                *
 *   5. Post Isend for each sender-PE.                                      *
 *   6. Waitall on (rPEnum + sPEnum) requests.                              *
 *   7. Unpack recv buffer per receiver-PE into rlist[] slots in field[].   *
 * ------------------------------------------------------------------------ */
void fesom_halo_exchange(real_t          *field,
                         fesom_halo_kind  kind,
                         int              n_levels,
                         int              n_components,
                         fesom_partit    *p)
{
    if (p->npes == 1) return;          /* no halos */

    fesom_com_struct *cs = get_com(p, kind);
    if (cs->rPEnum == 0 && cs->sPEnum == 0) return;

    halo_scratch *sc = &g_scratch[kind];
    int    stride = n_levels * n_components;
    size_t recv_total = (size_t)(cs->rptr[cs->rPEnum] - cs->rptr[0]) * stride;
    size_t send_total = (size_t)(cs->sptr[cs->sPEnum] - cs->sptr[0]) * stride;

    halo_nuke(field, cs, n_levels, n_components);

    scratch_grow_buf(&sc->recv_buf, &sc->recv_cap, recv_total);
    scratch_grow_buf(&sc->send_buf, &sc->send_cap, send_total);
    scratch_grow_reqs(sc, cs->rPEnum + cs->sPEnum);

    int tag    = kind_tag(kind);
    int n_reqs = 0;

    /* 3. Post receives. Each receiver-PE k has rptr[k+1]-rptr[k] entries
     *    starting at recv_buf + (rptr[k]-rptr[0])*stride. */
    for (int k = 0; k < cs->rPEnum; ++k) {
        int    nseg  = cs->rptr[k + 1] - cs->rptr[k];
        size_t off   = (size_t)(cs->rptr[k] - cs->rptr[0]) * stride;
        int    nelem = nseg * stride;
        MPI_CHECK(MPI_Irecv(sc->recv_buf + off, nelem, MPI_DOUBLE,
                            cs->rPE[k], tag,
                            p->MPI_COMM_FESOM, &sc->reqs[n_reqs++]));
    }

    /* 4 + 5. Pack per sender-PE then post send. */
    for (int k = 0; k < cs->sPEnum; ++k) {
        int    nseg  = cs->sptr[k + 1] - cs->sptr[k];
        size_t off   = (size_t)(cs->sptr[k] - cs->sptr[0]) * stride;
        int    nelem = nseg * stride;
        real_t *dst  = sc->send_buf + off;
        for (int j = 0; j < nseg; ++j) {
            int local_idx = cs->slist[cs->sptr[k] - 1 + j] - 1;   /* 1-based → 0 */
            const real_t *src = field + (size_t)local_idx * stride;
            memcpy(dst + (size_t)j * stride, src, (size_t)stride * sizeof(real_t));
        }
        MPI_CHECK(MPI_Isend(sc->send_buf + off, nelem, MPI_DOUBLE,
                            cs->sPE[k], tag,
                            p->MPI_COMM_FESOM, &sc->reqs[n_reqs++]));
    }

    /* 6. Wait. */
    MPI_CHECK(MPI_Waitall(n_reqs, sc->reqs, MPI_STATUSES_IGNORE));

    /* 7. Unpack into halo region. */
    for (int k = 0; k < cs->rPEnum; ++k) {
        int    nseg = cs->rptr[k + 1] - cs->rptr[k];
        size_t off  = (size_t)(cs->rptr[k] - cs->rptr[0]) * stride;
        const real_t *src = sc->recv_buf + off;
        for (int j = 0; j < nseg; ++j) {
            int local_idx = cs->rlist[cs->rptr[k] - 1 + j] - 1;   /* 1-based → 0 */
            real_t *dst   = field + (size_t)local_idx * stride;
            memcpy(dst, src + (size_t)j * stride, (size_t)stride * sizeof(real_t));
        }
    }
}

/* ------------------------------------------------------------------------ *
 * Identity smoke test (slice 30b acceptance gate).                         *
 *                                                                          *
 * Strategy: pack each node's global ID into f[local_idx], halo-exchange,   *
 * then verify f[halo_idx] == myList_nod2D[halo_idx] (i.e., the halo entry  *
 * carries the global ID of the node it represents). This works regardless  *
 * of which rank owns the node — the contract is that exchange writes the   *
 * sender's value of f at the same global node into the receiver's halo.    *
 * ------------------------------------------------------------------------ */
void fesom_halo_identity_test(fesom_partit *p)
{
    if (p->npes == 1) {
        if (p->mype == 0) {
            printf("[fesom_halo] identity test: skipped (1 rank, no halo)\n");
        }
        return;
    }

    int N = p->myDim_nod2D + p->eDim_nod2D;
    real_t *f = malloc((size_t)N * sizeof(real_t));
    FESOM_CHECK(f, "fesom_halo: identity test field alloc");

    /* Round 1: positive test. */
    for (int i = 0; i < p->myDim_nod2D; ++i) {
        f[i] = (real_t)p->myList_nod2D[i];        /* gid of interior node */
    }
    for (int h = 0; h < p->eDim_nod2D; ++h) {
        f[p->myDim_nod2D + h] = (real_t)(-1);     /* sentinel */
    }

    fesom_exchange_nod2D(f, p);

    int local_failures = 0;
    for (int h = 0; h < p->eDim_nod2D; ++h) {
        real_t got      = f[p->myDim_nod2D + h];
        int    expected = p->myList_nod2D[p->myDim_nod2D + h];   /* expected gid */
        if ((int)got != expected) {
            if (local_failures < 5) {
                fprintf(stderr,
                        "[fesom_halo] rank %d FAIL: halo[%d] got=%g expected_gid=%d\n",
                        p->mype, h, (double)got, expected);
            }
            local_failures++;
        }
    }
    int global_failures = 0;
    MPI_CHECK(MPI_Allreduce(&local_failures, &global_failures, 1, MPI_INT, MPI_SUM,
                            p->MPI_COMM_FESOM));
    if (global_failures > 0) {
        FESOM_DIE("fesom_halo identity test: %d halo nodes mismatched (rank %d had %d locally)",
                  global_failures, p->mype, local_failures);
    }
    if (p->mype == 0) {
        printf("[fesom_halo] identity test (positive): all halo entries carry correct gid\n");
    }

    /* Round 2: negative test — corrupt halo[0] on rank 0, re-exchange, verify
     * the correct value is restored. */
    if (p->eDim_nod2D > 0) {
        if (p->mype == 0) {
            f[p->myDim_nod2D + 0] = (real_t)(-99);
        }
        fesom_exchange_nod2D(f, p);
        local_failures = 0;
        if (p->mype == 0) {
            real_t got      = f[p->myDim_nod2D + 0];
            int    expected = p->myList_nod2D[p->myDim_nod2D + 0];
            if ((int)got != expected) local_failures = 1;
        }
        MPI_CHECK(MPI_Allreduce(&local_failures, &global_failures, 1, MPI_INT, MPI_SUM,
                                p->MPI_COMM_FESOM));
        if (global_failures > 0) {
            FESOM_DIE("fesom_halo identity test: corruption-recovery FAILED — "
                      "exchange did not overwrite halo[0] on rank 0");
        }
        if (p->mype == 0) {
            printf("[fesom_halo] identity test (corruption recovery): exchange overwrote stale halo\n");
        }
    }

    free(f);
}
