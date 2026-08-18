/*
 * fesom_io_restart.c — see fesom_io_restart.h and
 * docs/plans/20260818-restart-io-port.md.
 *
 * The file carries every array whose value at the start of a step was not
 * produced by that step. Where a field is in principle derivable — helem,
 * zbar_3d_n, Z_3d_n from hnode — it is written anyway: the gate is a
 * bit-for-bit round trip, and storing the value is cheaper to be sure of than
 * an argument that a re-derivation rounds the same way. The one thing that is
 * re-derived on read is hnode_new, which the step overwrites before reading.
 */
#include "fesom_io_restart.h"

#include "fesom_dyn.h"
#include "fesom_halo.h"
#include "fesom_ice_types.h"
#include "fesom_io_gather.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_ssh.h"
#include "fesom_tke.h"
#include "fesom_tracers.h"

#include <netcdf.h>
#include <stdlib.h>
#include <string.h>

#define NC_CHECK(call) do {                                              \
    int _ec = (call);                                                    \
    if (_ec != NC_NOERR) {                                               \
        FESOM_DIE("restart: NetCDF error: %s", nc_strerror(_ec));        \
    }                                                                    \
} while (0)

/* ------------------------------------------------------------------ *
 * The field table                                                     *
 * ------------------------------------------------------------------ */

typedef struct {
    const char     *name[2];   /* one netCDF variable per component      */
    real_t         *base;      /* local array, NULL = field not present  */
    int             is_elem;   /* 0 = node-indexed, 1 = element-indexed  */
    int             levels;    /* 1 for a surface field, nl for a column */
    int             ncomp;     /* components interleaved innermost       */
    int             halo;      /* fesom_halo_kind, or -1 for none        */
} rst_var;

#define RST_MAX 40

/* halo == -1: the readers of this field loop owned entries only, so its halo
 * copies are never read and exchanging them would be a lie about what the
 * continuous run holds there. That is true of the EVP stresses (the mEVP
 * element loop is bounded by myDim_elem2D), of dhe (filled over owned
 * elements and consumed by update_stiff_mat_ale over the same range) and of
 * tke (see its entry below). */
static int rst_build(rst_var *v,
                     const struct fesom_mesh *mesh,
                     struct fesom_dyn        *dyn,
                     struct fesom_tracers    *tr,
                     struct fesom_ice        *ice,
                     struct fesom_tke        *tke)
{
    const int nl = mesh->nl;
    int n = 0;

    /* --- ocean, node, surface ---------------------------------------- */
    v[n++] = (rst_var){ {"ssh",         NULL}, dyn->eta_n,       0, 1, 1, FESOM_HALO_NOD2D };
    v[n++] = (rst_var){ {"ssh_rhs_old", NULL}, dyn->ssh_rhs_old, 0, 1, 1, FESOM_HALO_NOD2D };
    /* Not in the Fortran restart: the SSH solve warm-starts from the previous
     * step's d_eta (fesom_ssh_solve_cg takes dyn->d_eta as its initial guess
     * and never zeroes it). Starting from zero instead lands on a different
     * point inside the solver tolerance — a difference of order the tolerance,
     * every step, which is what the round-trip gate first measured. */
    v[n++] = (rst_var){ {"d_eta",       NULL}, dyn->d_eta,       0, 1, 1, FESOM_HALO_NOD2D };
    v[n++] = (rst_var){ {"hbar",        NULL}, ((struct fesom_mesh *)mesh)->hbar,
                                                                 0, 1, 1, FESOM_HALO_NOD2D };
    v[n++] = (rst_var){ {"hbar_old",    NULL}, ((struct fesom_mesh *)mesh)->hbar_old,
                                                                 0, 1, 1, FESOM_HALO_NOD2D };

    /* --- ocean, node, column ----------------------------------------- */
    v[n++] = (rst_var){ {"hnode",     NULL}, ((struct fesom_mesh *)mesh)->hnode,
                                                                 0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"zbar_3d_n", NULL}, ((struct fesom_mesh *)mesh)->zbar_3d_n,
                                                                 0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"Z_3d_n",    NULL}, ((struct fesom_mesh *)mesh)->Z_3d_n,
                                                                 0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"temp",    NULL}, tr->data[FESOM_TRACER_T].values,   0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"temp_AB", NULL}, tr->data[FESOM_TRACER_T].valuesAB, 0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"temp_M1", NULL}, tr->data[FESOM_TRACER_T].valuesold,0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"salt",    NULL}, tr->data[FESOM_TRACER_S].values,   0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"salt_AB", NULL}, tr->data[FESOM_TRACER_S].valuesAB, 0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"salt_M1", NULL}, tr->data[FESOM_TRACER_S].valuesold,0, nl, 1, FESOM_HALO_NOD3D };
    /* Not in the Fortran restart, and the subtlest of the additions: uvnode is
     * written in the middle of the step (fesom_step.c:136) but read at the TOP
     * of the next iteration by the bulk formula (fesom_bulk.c:261), so at the
     * moment a restart is written it holds the velocity from the step BEFORE
     * the one being saved. Recomputing it from the saved uv gives a value one
     * step too new, and the surface stress — hence TKE, hence everything — is
     * wrong from the first step of the resumed run. */
    v[n++] = (rst_var){ {"uvnode_u", "uvnode_v"}, dyn->uvnode, 0, nl, 2, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"w",       NULL}, dyn->w,   0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"w_expl",  NULL}, dyn->w_e, 0, nl, 1, FESOM_HALO_NOD3D };
    v[n++] = (rst_var){ {"w_impl",  NULL}, dyn->w_i, 0, nl, 1, FESOM_HALO_NOD3D };
    /* tke is present only under FESOM_MIX_SCHEME=TKE, exactly as the Fortran
     * restart registers it only for mix_scheme_nmb 5 or 56.
     * halo == -1: calc_cvmix_tke runs over myDim_nod2D and tke is never
     * exchanged (see the header of fesom_tke.c), so its halo is zero for the
     * life of an uninterrupted run. Exchanging it here would fill it with real
     * values and leave the resumed run holding something the continuous run
     * never held — harmless, since nothing reads it, but a difference for no
     * reason, and the full-extent probe reported it as one. */
    v[n++] = (rst_var){ {"tke", NULL}, tke ? tke->tke : NULL, 0, nl, 1, -1 };

    /* --- ocean, element ---------------------------------------------- */
    v[n++] = (rst_var){ {"helem", NULL}, ((struct fesom_mesh *)mesh)->helem,
                                                                 1, nl, 1, FESOM_HALO_ELEM3D };
    v[n++] = (rst_var){ {"dhe",   NULL}, ((struct fesom_mesh *)mesh)->dhe, 1, 1, 1, -1 };
    v[n++] = (rst_var){ {"u", "v"},              dyn->uv,       1, nl, 2, FESOM_HALO_ELEM3D };
    v[n++] = (rst_var){ {"urhs_AB", "vrhs_AB"},  dyn->uv_rhsAB, 1, nl, 2, FESOM_HALO_ELEM3D };

    /* --- sea ice ------------------------------------------------------ */
    if (ice) {
        v[n++] = (rst_var){ {"area",  NULL}, ice->data[FESOM_ICE_AICE].values,  0, 1, 1, FESOM_HALO_NOD2D };
        v[n++] = (rst_var){ {"hice",  NULL}, ice->data[FESOM_ICE_MICE].values,  0, 1, 1, FESOM_HALO_NOD2D };
        v[n++] = (rst_var){ {"hsnow", NULL}, ice->data[FESOM_ICE_MSNOW].values, 0, 1, 1, FESOM_HALO_NOD2D };
        v[n++] = (rst_var){ {"uice",  NULL}, ice->uice, 0, 1, 1, FESOM_HALO_NOD2D };
        v[n++] = (rst_var){ {"vice",  NULL}, ice->vice, 0, 1, 1, FESOM_HALO_NOD2D };
        /* Not in the Fortran restart: read before written in the
         * thermodynamics (fesom_ice_thermo.c:470 before :508). */
        v[n++] = (rst_var){ {"t_skin", NULL}, ice->thermo.t_skin, 0, 1, 1, FESOM_HALO_NOD2D };
        /* Not in the Fortran restart: the EVP stress tensor is updated
         * incrementally and never zeroed on entry. */
        v[n++] = (rst_var){ {"sigma11", NULL}, ice->work.sigma11, 1, 1, 1, -1 };
        v[n++] = (rst_var){ {"sigma12", NULL}, ice->work.sigma12, 1, 1, 1, -1 };
        v[n++] = (rst_var){ {"sigma22", NULL}, ice->work.sigma22, 1, 1, 1, -1 };
    }

    FESOM_CHECK(n <= RST_MAX, "restart: field table overflow (%d)", n);
    return n;
}

/* ------------------------------------------------------------------ *
 * Canonicalising the replicated elements (FESOM_RESTART_CANONICALIZE)  *
 *                                                                     *
 * Nodes are owned exclusively; elements are not. gen_comm.F90:264-270
 * puts an element in myList_elem2D of EVERY rank that owns one of its
 * nodes, and each of those ranks then computes it redundantly. For a
 * quantity assembled from its own three nodes that is harmless -- same
 * inputs, same result. It is not harmless for one assembled by an edge
 * loop: visc_filt_bcksct sums into U_b over myDim_edge2D+eDim_edge2D
 * (oce_dyn.F90:330-360, mirrored in fesom_momentum.c), and two ranks
 * reach the same element's three edges in a different order, so they get
 * the same sum with different rounding. exchange_elem does not repair
 * it: the element is owned on both ranks, so it is in neither one's halo.
 *
 * The disagreement is real and it grows. Measured on CORE2 at 128 ranks
 * (FESOM_IO_GATHER_DUPCHECK=1): uv agrees on all owners after one step,
 * differs on 167 of 15348 replicated slots after two (max 7e-18) and on
 * 3311 after twenty (max 4e-16); uv_rhsAB reaches 3e-11. It is inherited
 * from the Fortran, not introduced by the port.
 *
 * The consequence for a restart is unavoidable. The state of a running
 * model is not a function of the global fields alone -- it includes which
 * rank holds which rounding of those 15348 slots -- and a
 * partition-independent file cannot carry that. No file can reproduce an
 * uninterrupted run bit for bit, whatever is added to it.
 *
 * What CAN be exact is the case that matters operationally: a run that
 * checkpoints at step N and carries on, against a run resumed from that
 * checkpoint. Writing the restart pushes the gathered values back into
 * the live arrays, so the writer leaves the model in exactly the state
 * the reader will produce and every owner of a replicated element agrees
 * from then on. The cost is one extra scatter per element field per
 * checkpoint, and a perturbation of the trajectory the size of the
 * disagreement being removed.
 *
 * It does mean the trajectory depends on the checkpoint interval, at
 * roundoff. Set FESOM_RESTART_CANONICALIZE=0 to measure that, or to see
 * the raw disagreement; the round-trip gate then fails by design.
 * ------------------------------------------------------------------ */
static int rst_canonicalise(void)
{
    const char *e = getenv("FESOM_RESTART_CANONICALIZE");
    return e ? atoi(e) != 0 : 1;
}

/* ------------------------------------------------------------------ *
 * The ALE stiffness matrix                                            *
 *                                                                     *
 * fesom_update_stiff_mat_ale adds one dhe-weighted increment per step,
 * so the matrix is the integral of the run, not a function of its
 * current state. A resumed run that rebuilds the base matrix and then
 * applies only the last increment is missing every increment before
 * it — which is what made the first round-trip attempt disagree in the
 * third decimal place rather than in the last bit.
 *
 * The matrix is stored per owned row, columns sorted by GLOBAL node id,
 * so the file does not depend on how the mesh was partitioned. Column
 * ids travel with the values and are checked on read: a sparsity
 * pattern that does not match means the file belongs to another mesh.
 * ------------------------------------------------------------------ */

/* Positions of row `row`'s entries, ordered by global column id. Rows have a
 * node's neighbours in them — seven or so entries — so insertion sort is the
 * right tool. Returns the row length. */
static int stiff_row_order(const fesom_ssh_stiff *S, const fesom_partit *p,
                           int row, int *pos, int *gcol)
{
    const int a = S->rowptr[row], b = S->rowptr[row + 1];
    const int len = b - a;
    for (int i = 0; i < len; ++i) {
        pos[i]  = a + i;
        gcol[i] = p->myList_nod2D[S->colind[a + i]];
    }
    for (int i = 1; i < len; ++i) {
        int g = gcol[i], q = pos[i], j = i - 1;
        while (j >= 0 && gcol[j] > g) { gcol[j + 1] = gcol[j]; pos[j + 1] = pos[j]; --j; }
        gcol[j + 1] = g; pos[j + 1] = q;
    }
    return len;
}

static int stiff_max_row(const fesom_ssh_stiff *S)
{
    int mx = 0;
    for (int r = 0; r < S->dim; ++r) {
        int len = S->rowptr[r + 1] - S->rowptr[r];
        if (len > mx) mx = len;
    }
    return mx;
}

/* Flatten this rank's rows into (row length, global column, value) triples in
 * local row order. Caller frees. */
static void stiff_pack(const fesom_ssh_stiff *S, const fesom_partit *p,
                       int **len_out, int **gcol_out, real_t **val_out)
{
    const int rows = S->dim, nnz = S->rowptr[rows];
    int    *len  = malloc((size_t)(rows > 0 ? rows : 1) * sizeof(int));
    int    *gcol = malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(int));
    real_t *val  = malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(real_t));
    int    *pos  = malloc((size_t)(stiff_max_row(S) + 1) * sizeof(int));
    int    *gtmp = malloc((size_t)(stiff_max_row(S) + 1) * sizeof(int));
    FESOM_CHECK(len && gcol && val && pos && gtmp, "restart: stiff pack alloc");
    int k = 0;
    for (int r = 0; r < rows; ++r) {
        len[r] = stiff_row_order(S, p, r, pos, gtmp);
        for (int i = 0; i < len[r]; ++i) {
            gcol[k] = gtmp[i];
            val[k]  = S->values[pos[i]];
            ++k;
        }
    }
    FESOM_CHECK(k == nnz, "restart: stiff pack counted %d of %d entries", k, nnz);
    free(pos); free(gtmp);
    *len_out = len; *gcol_out = gcol; *val_out = val;
}

/* Rank 0 assembles the global CSR, ordered by global row id. Returns the
 * global nnz; the three out arrays are allocated on rank 0 only. */
static long long stiff_gather(const fesom_ssh_stiff *S, const gather_plan *gp,
                              const struct fesom_mesh *mesh, fesom_partit *p,
                              int **glen, int **gcol, real_t **gval)
{
    int *len = NULL, *col = NULL; real_t *val = NULL;
    stiff_pack(S, p, &len, &col, &val);
    const int mynnz = S->rowptr[S->dim];

    int *nnz_all = NULL, *nnz_displ = NULL, *len_ro = NULL;
    if (gp->mype == 0) {
        nnz_all   = malloc((size_t)gp->npes * sizeof(int));
        nnz_displ = malloc((size_t)gp->npes * sizeof(int));
        len_ro    = malloc((size_t)gp->total_n * sizeof(int));
        FESOM_CHECK(nnz_all && nnz_displ && len_ro, "restart: stiff gather alloc");
    }
    MPI_Gather(&mynnz, 1, MPI_INT, nnz_all, 1, MPI_INT, 0, p->MPI_COMM_FESOM);
    MPI_Gatherv(len, S->dim, MPI_INT, len_ro, gp->all_n, gp->displ_n, MPI_INT,
                0, p->MPI_COMM_FESOM);

    long long total_nnz = 0;
    int *col_ro = NULL; real_t *val_ro = NULL;
    if (gp->mype == 0) {
        nnz_displ[0] = 0;
        for (int r = 1; r < gp->npes; ++r) nnz_displ[r] = nnz_displ[r-1] + nnz_all[r-1];
        total_nnz = (long long)nnz_displ[gp->npes-1] + nnz_all[gp->npes-1];
        col_ro = malloc((size_t)total_nnz * sizeof(int));
        val_ro = malloc((size_t)total_nnz * sizeof(real_t));
        FESOM_CHECK(col_ro && val_ro, "restart: stiff gather alloc (%lld nnz)", total_nnz);
    }
    MPI_Gatherv(col, mynnz, MPI_INT,    col_ro, nnz_all, nnz_displ, MPI_INT,
                0, p->MPI_COMM_FESOM);
    MPI_Gatherv(val, mynnz, MPI_DOUBLE, val_ro, nnz_all, nnz_displ, MPI_DOUBLE,
                0, p->MPI_COMM_FESOM);
    free(len); free(col); free(val);

    if (gp->mype == 0) {
        /* rank order -> global row order */
        *glen = malloc((size_t)mesh->nod2D * sizeof(int));
        *gcol = malloc((size_t)total_nnz * sizeof(int));
        *gval = malloc((size_t)total_nnz * sizeof(real_t));
        FESOM_CHECK(*glen && *gcol && *gval, "restart: stiff assemble alloc");
        long long *start = malloc((size_t)mesh->nod2D * sizeof(long long));
        FESOM_CHECK(start, "restart: stiff assemble alloc");
        for (int i = 0; i < mesh->nod2D; ++i) (*glen)[i] = -1;

        long long here = 0;
        for (int r = 0; r < gp->npes; ++r) {
            for (int i = 0; i < gp->all_n[r]; ++i) {
                int gid = gp->gathered_myList_n[gp->displ_n[r] + i] - 1;
                (*glen)[gid] = len_ro[gp->displ_n[r] + i];
            }
        }
        for (int i = 0; i < mesh->nod2D; ++i) {
            FESOM_CHECK((*glen)[i] >= 0, "restart: global row %d has no owner", i);
            start[i] = here;
            here += (*glen)[i];
        }
        FESOM_CHECK(here == total_nnz, "restart: stiff row lengths sum to %lld, not %lld",
                    here, total_nnz);
        for (int r = 0; r < gp->npes; ++r) {
            long long off = nnz_displ[r];
            for (int i = 0; i < gp->all_n[r]; ++i) {
                int gid = gp->gathered_myList_n[gp->displ_n[r] + i] - 1;
                int len_i = len_ro[gp->displ_n[r] + i];
                memcpy(*gcol + start[gid], col_ro + off, (size_t)len_i * sizeof(int));
                memcpy(*gval + start[gid], val_ro + off, (size_t)len_i * sizeof(real_t));
                off += len_i;
            }
        }
        free(start); free(col_ro); free(val_ro);
        free(nnz_all); free(nnz_displ); free(len_ro);
    }
    return total_nnz;
}

/* The inverse: rank 0 holds the global CSR, every rank ends up with its own
 * rows' values in its own column order. */
static void stiff_scatter(fesom_ssh_stiff *S, const gather_plan *gp,
                          const struct fesom_mesh *mesh, fesom_partit *p,
                          const int *glen, const int *gcol, const real_t *gval)
{
    const int mynnz = S->rowptr[S->dim];
    int *nnz_all = NULL, *nnz_displ = NULL;
    int *col_ro = NULL; real_t *val_ro = NULL;
    if (gp->mype == 0) {
        nnz_all   = malloc((size_t)gp->npes * sizeof(int));
        nnz_displ = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(nnz_all && nnz_displ, "restart: stiff scatter alloc");
    }
    MPI_Gather(&mynnz, 1, MPI_INT, nnz_all, 1, MPI_INT, 0, p->MPI_COMM_FESOM);
    if (gp->mype == 0) {
        nnz_displ[0] = 0;
        for (int r = 1; r < gp->npes; ++r) nnz_displ[r] = nnz_displ[r-1] + nnz_all[r-1];
        long long total = (long long)nnz_displ[gp->npes-1] + nnz_all[gp->npes-1];
        col_ro = malloc((size_t)total * sizeof(int));
        val_ro = malloc((size_t)total * sizeof(real_t));
        FESOM_CHECK(col_ro && val_ro, "restart: stiff scatter alloc (%lld nnz)", total);

        long long *start = malloc((size_t)mesh->nod2D * sizeof(long long));
        FESOM_CHECK(start, "restart: stiff scatter alloc");
        long long here = 0;
        for (int i = 0; i < mesh->nod2D; ++i) { start[i] = here; here += glen[i]; }
        long long off = 0;
        for (int r = 0; r < gp->npes; ++r) {
            for (int i = 0; i < gp->all_n[r]; ++i) {
                int gid = gp->gathered_myList_n[gp->displ_n[r] + i] - 1;
                memcpy(col_ro + off, gcol + start[gid], (size_t)glen[gid] * sizeof(int));
                memcpy(val_ro + off, gval + start[gid], (size_t)glen[gid] * sizeof(real_t));
                off += glen[gid];
            }
        }
        free(start);
    }
    int *col = malloc((size_t)(mynnz > 0 ? mynnz : 1) * sizeof(int));
    real_t *val = malloc((size_t)(mynnz > 0 ? mynnz : 1) * sizeof(real_t));
    FESOM_CHECK(col && val, "restart: stiff scatter alloc");
    MPI_Scatterv(col_ro, nnz_all, nnz_displ, MPI_INT,    col, mynnz, MPI_INT,
                 0, p->MPI_COMM_FESOM);
    MPI_Scatterv(val_ro, nnz_all, nnz_displ, MPI_DOUBLE, val, mynnz, MPI_DOUBLE,
                 0, p->MPI_COMM_FESOM);
    if (gp->mype == 0) { free(col_ro); free(val_ro); free(nnz_all); free(nnz_displ); }

    int *pos  = malloc((size_t)(stiff_max_row(S) + 1) * sizeof(int));
    int *gtmp = malloc((size_t)(stiff_max_row(S) + 1) * sizeof(int));
    FESOM_CHECK(pos && gtmp, "restart: stiff unpack alloc");
    int k = 0;
    for (int r = 0; r < S->dim; ++r) {
        int len = stiff_row_order(S, p, r, pos, gtmp);
        for (int i = 0; i < len; ++i) {
            FESOM_CHECK(col[k] == gtmp[i],
                        "restart: stiffness sparsity differs at row %d (file column %d, "
                        "mesh column %d) — the file is for a different mesh",
                        r, col[k], gtmp[i]);
            S->values[pos[i]] = val[k];
            ++k;
        }
    }
    FESOM_CHECK(k == mynnz, "restart: stiff unpack counted %d of %d", k, mynnz);
    free(pos); free(gtmp); free(col); free(val);
}

/* ------------------------------------------------------------------ *
 * Configuration                                                       *
 * ------------------------------------------------------------------ */

void fesom_restart_config(fesom_restart_cfg *cfg)
{
    const char *out   = getenv("FESOM_RESTART_OUT");
    const char *in    = getenv("FESOM_RESTART_IN");
    const char *every = getenv("FESOM_RESTART_EVERY");
    const char *at    = getenv("FESOM_RESTART_AT");
    cfg->out_dir = (out && out[0]) ? out : NULL;
    cfg->in_path = (in  && in[0])  ? in  : NULL;
    cfg->every   = (every && every[0]) ? atoi(every) : 0;
    FESOM_CHECK(cfg->every >= 0, "FESOM_RESTART_EVERY must not be negative");
    cfg->n_at = 0;
    for (const char *p = at; p && *p; ) {
        char *end;
        long v = strtol(p, &end, 10);
        FESOM_CHECK(end != p, "FESOM_RESTART_AT: not a number at \"%s\"", p);
        FESOM_CHECK(v > 0, "FESOM_RESTART_AT: step %ld is not positive", v);
        FESOM_CHECK(cfg->n_at < FESOM_RESTART_MAX_AT,
                    "FESOM_RESTART_AT: more than %d steps", FESOM_RESTART_MAX_AT);
        cfg->at[cfg->n_at++] = (int)v;
        while (*end == ',' || *end == ' ') ++end;
        p = end;
    }
}

int fesom_restart_due(const fesom_restart_cfg *cfg, int step, int nsteps)
{
    if (!cfg->out_dir) return 0;
    if (step == nsteps) return 1;
    for (int i = 0; i < cfg->n_at; ++i)
        if (cfg->at[i] == step) return 1;
    return cfg->every > 0 && step % cfg->every == 0;
}

void fesom_restart_path(char *buf, size_t buflen,
                        const char *out_dir, const fesom_calendar_t *cal)
{
    int sod = cal->hour * 3600 + cal->minute * 60 + (int)cal->second;
    snprintf(buf, buflen, "%s/fesom.%04d%02d%02d_%05d.restart.nc",
             out_dir, cal->year, cal->month, cal->day, sod);
}

/* ------------------------------------------------------------------ *
 * Shared helpers                                                      *
 * ------------------------------------------------------------------ */

static size_t rst_count(const rst_var *v, const struct fesom_mesh *mesh)
{
    return (size_t)(v->is_elem ? mesh->elem2D : mesh->nod2D);
}

/* Global sum of nlevels_nod2D over owned nodes — the same number the campaign
 * uses to tell the pre-corruption CORE2 mesh (3 832 750) from the /pool copy
 * that was altered on 2026-07-03 (3 832 745). Stored in the file so a restart
 * cannot be replayed onto a different mesh without being caught. */
static long long rst_mesh_signature(const struct fesom_mesh *mesh, fesom_partit *p)
{
    long long local = 0;
    for (int n = 0; n < mesh->myDim_nod2D; ++n) local += mesh->nlevels_nod2D[n];
    long long global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_LONG_LONG, MPI_SUM, p->MPI_COMM_FESOM);
    return global;
}

/* ------------------------------------------------------------------ *
 * Write                                                               *
 * ------------------------------------------------------------------ */

void fesom_restart_write(const char                 *path,
                         int                         step,
                         double                      dt,
                         const fesom_calendar_t     *cal,
                         const struct fesom_mesh    *mesh,
                         const struct fesom_dyn     *dyn,
                         const struct fesom_tracers *tracers,
                         const struct fesom_ice     *ice,
                         const struct fesom_tke     *tke,
                         const struct fesom_ssh_stiff *stiff,
                         struct fesom_partit        *partit)
{
    rst_var vars[RST_MAX];
    /* The table holds writable pointers because the reader shares it; nothing
     * in this function writes through them. */
    int nvars = rst_build(vars, mesh,
                          (struct fesom_dyn *)dyn,
                          (struct fesom_tracers *)tracers,
                          (struct fesom_ice *)ice,
                          (struct fesom_tke *)tke);

    gather_plan gp;
    gather_plan_init(&gp, mesh, partit);
    const long long signature = rst_mesh_signature(mesh, partit);

    /* The stiffness matrix is gathered before the file is created because its
     * global nnz is one of the dimensions. */
    int *slen = NULL, *scol = NULL; real_t *sval = NULL;
    const long long stiff_nnz =
        stiff_gather((const fesom_ssh_stiff *)stiff, &gp, mesh, partit, &slen, &scol, &sval);

    int ncid = -1, dim_node = -1, dim_elem = -1, dim_nz = -1, dim_nnz = -1;
    int varid[RST_MAX][2];
    int vid_slen = -1, vid_scol = -1, vid_sval = -1;
    if (partit->mype == 0) {
        /* CDF-5, not netCDF-4/HDF5. The file is a flat set of large 1-D and
         * 2-D double arrays with no need for groups, chunking or compression,
         * and the HDF5 layer failed on the second 1.2 GB write within one run
         * ("NetCDF: HDF error" from nc_put_var_double). CDF-5 carries 64-bit
         * sizes and NC_INT64 attributes, which is all this format uses. */
        NC_CHECK(nc_create(path, NC_CLOBBER | NC_64BIT_DATA, &ncid));
        NC_CHECK(nc_def_dim(ncid, "node", (size_t)mesh->nod2D,  &dim_node));
        NC_CHECK(nc_def_dim(ncid, "elem", (size_t)mesh->elem2D, &dim_elem));
        NC_CHECK(nc_def_dim(ncid, "nz",   (size_t)mesh->nl,     &dim_nz));
        NC_CHECK(nc_def_dim(ncid, "stiff_nnz", (size_t)stiff_nnz, &dim_nnz));
        NC_CHECK(nc_def_var(ncid, "stiff_row_len", NC_INT,    1, &dim_node, &vid_slen));
        NC_CHECK(nc_def_var(ncid, "stiff_colind",  NC_INT,    1, &dim_nnz,  &vid_scol));
        NC_CHECK(nc_def_var(ncid, "stiff_values",  NC_DOUBLE, 1, &dim_nnz,  &vid_sval));
        for (int i = 0; i < nvars; ++i) {
            for (int c = 0; c < vars[i].ncomp; ++c) {
                varid[i][c] = -1;
                if (!vars[i].base) continue;
                int dims[2];
                int ndims = 0;
                if (vars[i].levels > 1) dims[ndims++] = dim_nz;
                dims[ndims++] = vars[i].is_elem ? dim_elem : dim_node;
                NC_CHECK(nc_def_var(ncid, vars[i].name[c], NC_DOUBLE, ndims, dims,
                                    &varid[i][c]));
            }
        }
        int fmt = FESOM_RESTART_FORMAT_VERSION;
        NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "restart_format_version", NC_INT, 1, &fmt));
        NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "step", NC_INT, 1, &step));
        NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "dt", NC_DOUBLE, 1, &dt));
        int cali[6] = { (int)cal->kind, cal->year, cal->month, cal->day,
                        cal->hour, cal->minute };
        NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "calendar", NC_INT, 6, cali));
        NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "second", NC_DOUBLE, 1, &cal->second));
        int origin[3] = { cal->origin_year, cal->origin_month, cal->origin_day };
        NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "calendar_origin", NC_INT, 3, origin));
        long long sig = signature;
        NC_CHECK(nc_put_att_longlong(ncid, NC_GLOBAL, "sum_nlevels_nod2D", NC_INT64, 1, &sig));
        if (ice) {
            int isu[2] = { ice->ice_steps_since_upd, ice->ice_update };
            NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "ice_step_state", NC_INT, 2, isu));
        }
        NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "source", 11, "fesom_port"));
        NC_CHECK(nc_enddef(ncid));
    }

    const int canon = rst_canonicalise();

    for (int i = 0; i < nvars; ++i) {
        if (!vars[i].base) continue;
        const size_t count  = rst_count(&vars[i], mesh);
        const int    stride = vars[i].levels * vars[i].ncomp;
        real_t *global = NULL, *plane = NULL;
        if (partit->mype == 0) {
            global = malloc(count * (size_t)stride * sizeof(real_t));
            plane  = malloc(count * (size_t)vars[i].levels * sizeof(real_t));
            FESOM_CHECK(global && plane, "restart: gather buffer alloc for %s",
                        vars[i].name[0]);
        }
        gather_set_label(vars[i].name[0]);
        if (vars[i].is_elem)
            gather_elem(vars[i].base, stride, &gp, global, partit->MPI_COMM_FESOM);
        else
            gather_node(vars[i].base, stride, &gp, global, partit->MPI_COMM_FESOM);

        /* Push the gathered element field straight back into the live arrays
         * (see rst_canonicalise above). Node fields are owned exclusively, so
         * for them this would be the identity; skipped. */
        if (vars[i].is_elem && canon) {
            scatter_elem(global, stride, &gp, vars[i].base, partit->MPI_COMM_FESOM);
            if (vars[i].halo >= 0)
                fesom_halo_exchange(vars[i].base, (fesom_halo_kind)vars[i].halo,
                                    vars[i].levels, vars[i].ncomp, partit);
        }

        if (partit->mype == 0) {
            for (int c = 0; c < vars[i].ncomp; ++c) {
                /* [entity][lev][comp] -> [lev][entity], the (nz, entity) order
                 * the Fortran restarts use. */
                for (size_t e = 0; e < count; ++e)
                    for (int lev = 0; lev < vars[i].levels; ++lev)
                        plane[(size_t)lev * count + e] =
                            global[e * (size_t)stride + (size_t)lev * vars[i].ncomp + c];
                NC_CHECK(nc_put_var_double(ncid, varid[i][c], plane));
            }
            free(global); free(plane);
        }
    }

    if (partit->mype == 0) {
        NC_CHECK(nc_put_var_int   (ncid, vid_slen, slen));
        NC_CHECK(nc_put_var_int   (ncid, vid_scol, scol));
        NC_CHECK(nc_put_var_double(ncid, vid_sval, sval));
        free(slen); free(scol); free(sval);
        NC_CHECK(nc_close(ncid));
    }
    gather_plan_free(&gp);
    MPI_Barrier(partit->MPI_COMM_FESOM);
}

/* ------------------------------------------------------------------ *
 * Read                                                                *
 * ------------------------------------------------------------------ */

void fesom_restart_read(const char           *path,
                        int                  *step,
                        double               *dt,
                        fesom_calendar_t     *cal,
                        struct fesom_mesh    *mesh,
                        struct fesom_dyn     *dyn,
                        struct fesom_tracers *tracers,
                        struct fesom_ice     *ice,
                        struct fesom_tke     *tke,
                        struct fesom_ssh_stiff *stiff,
                        struct fesom_partit  *partit)
{
    rst_var vars[RST_MAX];
    int nvars = rst_build(vars, mesh, dyn, tracers, ice, tke);

    gather_plan gp;
    gather_plan_init(&gp, mesh, partit);
    const long long signature = rst_mesh_signature(mesh, partit);

    int ncid = -1;
    int varid[RST_MAX][2];
    int header[9];              /* step, cal kind/y/m/d/h/min, origin y/m/d */
    double header_d[2];         /* dt, second */

    if (partit->mype == 0) {
        int ec = nc_open(path, NC_NOWRITE, &ncid);
        if (ec != NC_NOERR)
            FESOM_DIE("restart: cannot open %s: %s", path, nc_strerror(ec));

        int fmt = 0;
        NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "restart_format_version", &fmt));
        FESOM_CHECK(fmt == FESOM_RESTART_FORMAT_VERSION,
                    "restart: %s is format version %d, this build writes %d",
                    path, fmt, FESOM_RESTART_FORMAT_VERSION);

        /* Mesh identity. A restart replayed onto a different mesh produces a
         * plausible-looking run that is silently wrong, which is exactly the
         * failure mode the CORE2 bathymetry change created. */
        size_t len;
        int dimid;
        NC_CHECK(nc_inq_dimid(ncid, "node", &dimid));
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &len));
        FESOM_CHECK((int)len == mesh->nod2D,
                    "restart: %s has %zu nodes, this mesh has %d", path, len, mesh->nod2D);
        NC_CHECK(nc_inq_dimid(ncid, "elem", &dimid));
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &len));
        FESOM_CHECK((int)len == mesh->elem2D,
                    "restart: %s has %zu elements, this mesh has %d", path, len, mesh->elem2D);
        NC_CHECK(nc_inq_dimid(ncid, "nz", &dimid));
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &len));
        FESOM_CHECK((int)len == mesh->nl,
                    "restart: %s has %zu levels, this mesh has %d", path, len, mesh->nl);
        long long sig = 0;
        NC_CHECK(nc_get_att_longlong(ncid, NC_GLOBAL, "sum_nlevels_nod2D", &sig));
        FESOM_CHECK(sig == signature,
                    "restart: %s was written on a mesh with sum(nlevels_nod2D)=%lld, "
                    "this mesh sums to %lld — different bathymetry", path, sig, signature);

        NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "step", &header[0]));
        int cali[6];
        NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "calendar", cali));
        for (int i = 0; i < 6; ++i) header[1 + i] = cali[i];
        int origin[3];
        NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "calendar_origin", origin));
        header[7] = origin[0];
        NC_CHECK(nc_get_att_double(ncid, NC_GLOBAL, "dt", &header_d[0]));
        NC_CHECK(nc_get_att_double(ncid, NC_GLOBAL, "second", &header_d[1]));
        header[8] = origin[1] * 100 + origin[2];   /* month, day packed */

        for (int i = 0; i < nvars; ++i) {
            for (int c = 0; c < vars[i].ncomp; ++c) {
                varid[i][c] = -1;
                if (!vars[i].base) continue;
                int ec2 = nc_inq_varid(ncid, vars[i].name[c], &varid[i][c]);
                FESOM_CHECK(ec2 == NC_NOERR,
                            "restart: %s has no variable '%s' — this build needs it",
                            path, vars[i].name[c]);
            }
        }
        if (ice) {
            int isu[2];
            if (nc_get_att_int(ncid, NC_GLOBAL, "ice_step_state", isu) == NC_NOERR) {
                ice->ice_steps_since_upd = isu[0];
                ice->ice_update          = isu[1];
            }
        }
    }

    MPI_Bcast(header,   9, MPI_INT,    0, partit->MPI_COMM_FESOM);
    MPI_Bcast(header_d, 2, MPI_DOUBLE, 0, partit->MPI_COMM_FESOM);
    if (ice) {
        int isu[2] = { ice->ice_steps_since_upd, ice->ice_update };
        MPI_Bcast(isu, 2, MPI_INT, 0, partit->MPI_COMM_FESOM);
        ice->ice_steps_since_upd = isu[0];
        ice->ice_update          = isu[1];
    }
    *step = header[0];
    *dt   = header_d[0];
    cal->kind         = (fesom_calendar_kind_t)header[1];
    cal->year         = header[2];
    cal->month        = header[3];
    cal->day          = header[4];
    cal->hour         = header[5];
    cal->minute       = header[6];
    cal->second       = header_d[1];
    cal->origin_year  = header[7];
    cal->origin_month = header[8] / 100;
    cal->origin_day   = header[8] % 100;

    for (int i = 0; i < nvars; ++i) {
        if (!vars[i].base) continue;
        const size_t count  = rst_count(&vars[i], mesh);
        const int    stride = vars[i].levels * vars[i].ncomp;
        real_t *global = NULL, *plane = NULL;
        if (partit->mype == 0) {
            global = malloc(count * (size_t)stride * sizeof(real_t));
            plane  = malloc(count * (size_t)vars[i].levels * sizeof(real_t));
            FESOM_CHECK(global && plane, "restart: scatter buffer alloc for %s",
                        vars[i].name[0]);
            for (int c = 0; c < vars[i].ncomp; ++c) {
                NC_CHECK(nc_get_var_double(ncid, varid[i][c], plane));
                for (size_t e = 0; e < count; ++e)
                    for (int lev = 0; lev < vars[i].levels; ++lev)
                        global[e * (size_t)stride + (size_t)lev * vars[i].ncomp + c] =
                            plane[(size_t)lev * count + e];
            }
        }
        if (vars[i].is_elem)
            scatter_elem(global, stride, &gp, vars[i].base, partit->MPI_COMM_FESOM);
        else
            scatter_node(global, stride, &gp, vars[i].base, partit->MPI_COMM_FESOM);
        if (partit->mype == 0) { free(global); free(plane); }

        if (vars[i].halo >= 0)
            fesom_halo_exchange(vars[i].base, (fesom_halo_kind)vars[i].halo,
                                vars[i].levels, vars[i].ncomp, partit);
    }

    /* The accumulated free-surface stiffness matrix. */
    {
        int *slen = NULL, *scol = NULL; real_t *sval = NULL;
        if (partit->mype == 0) {
            size_t nnz;
            int dimid, vid;
            NC_CHECK(nc_inq_dimid(ncid, "stiff_nnz", &dimid));
            NC_CHECK(nc_inq_dimlen(ncid, dimid, &nnz));
            slen = malloc((size_t)mesh->nod2D * sizeof(int));
            scol = malloc(nnz * sizeof(int));
            sval = malloc(nnz * sizeof(real_t));
            FESOM_CHECK(slen && scol && sval, "restart: stiff read alloc");
            NC_CHECK(nc_inq_varid(ncid, "stiff_row_len", &vid));
            NC_CHECK(nc_get_var_int(ncid, vid, slen));
            NC_CHECK(nc_inq_varid(ncid, "stiff_colind", &vid));
            NC_CHECK(nc_get_var_int(ncid, vid, scol));
            NC_CHECK(nc_inq_varid(ncid, "stiff_values", &vid));
            NC_CHECK(nc_get_var_double(ncid, vid, sval));
        }
        stiff_scatter(stiff, &gp, mesh, partit, slen, scol, sval);
        if (partit->mype == 0) { free(slen); free(scol); free(sval); }
    }

    /* hnode_new is the one derived array not in the file: the ALE step opens by
     * setting it from hnode, so its incoming value is never read. Set it here
     * anyway so a mid-step inspection of a freshly restarted model is not
     * looking at the cold-start allocation. */
    memcpy(mesh->hnode_new, mesh->hnode,
           (size_t)(mesh->myDim_nod2D + mesh->eDim_nod2D) * (size_t)mesh->nl * sizeof(real_t));

    if (partit->mype == 0) NC_CHECK(nc_close(ncid));
    gather_plan_free(&gp);
    MPI_Barrier(partit->MPI_COMM_FESOM);
}
