/*
 * NetCDF snapshot writer with MPI gather to rank 0.
 *
 * Per snapshot: each rank packs its myDim slice of every output field; rank 0
 * gathers via MPI_Gatherv and unpacks into a global array indexed by global
 * node/element id (using each rank's myList_*). Rank 0 writes a single file
 * with `time` as an UNLIMITED dimension of length 1, so multiple snapshots
 * can be opened together via `xarray.open_mfdataset(...)`.
 */
#include "fesom_io.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_ice_types.h"
#include "fesom_io_config.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <mpi.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NC_CHECK(call) do {                                              \
    int _ec = (call);                                                    \
    if (_ec != NC_NOERR) {                                               \
        FESOM_DIE("NetCDF error: %s", nc_strerror(_ec));                 \
    }                                                                    \
} while (0)

/* ------------------------------------------------------------------ */
/* gather_plan implementation. Struct definition + public prototypes  */
/* live in fesom_io.h so fesom_io_stream.c can reuse the plumbing.    */
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

/* Integer variants for nlevels / elem_nodes. */
static void gather_node_int(const int *local, int stride,
                            const gather_plan *gp,
                            int *global,
                            MPI_Comm comm)
{
    int *recv = NULL, *counts = NULL, *displs = NULL;
    if (gp->mype == 0) {
        recv   = malloc((size_t)gp->total_n * (size_t)stride * sizeof(int));
        counts = malloc((size_t)gp->npes * sizeof(int));
        displs = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(recv && counts && displs, "gather_node_int: oom");
        for (int r = 0; r < gp->npes; ++r) {
            counts[r] = gp->all_n[r] * stride;
            displs[r] = gp->displ_n[r] * stride;
        }
    }
    MPI_Gatherv((void *)local, gp->myDim_n * stride, MPI_INT,
                recv, counts, displs, MPI_INT, 0, comm);
    if (gp->mype == 0) {
        for (int r = 0; r < gp->npes; ++r) {
            for (int i = 0; i < gp->all_n[r]; ++i) {
                int gid = gp->gathered_myList_n[gp->displ_n[r] + i] - 1;
                memcpy(global + (size_t)gid * stride,
                       recv + ((size_t)gp->displ_n[r] + i) * stride,
                       (size_t)stride * sizeof(int));
            }
        }
        free(recv); free(counts); free(displs);
    }
}

static void gather_elem_int(const int *local, int stride,
                            const gather_plan *gp,
                            int *global,
                            MPI_Comm comm)
{
    int *recv = NULL, *counts = NULL, *displs = NULL;
    if (gp->mype == 0) {
        recv   = malloc((size_t)gp->total_e * (size_t)stride * sizeof(int));
        counts = malloc((size_t)gp->npes * sizeof(int));
        displs = malloc((size_t)gp->npes * sizeof(int));
        FESOM_CHECK(recv && counts && displs, "gather_elem_int: oom");
        for (int r = 0; r < gp->npes; ++r) {
            counts[r] = gp->all_e[r] * stride;
            displs[r] = gp->displ_e[r] * stride;
        }
    }
    MPI_Gatherv((void *)local, gp->myDim_e * stride, MPI_INT,
                recv, counts, displs, MPI_INT, 0, comm);
    if (gp->mype == 0) {
        for (int r = 0; r < gp->npes; ++r) {
            for (int i = 0; i < gp->all_e[r]; ++i) {
                int gid = gp->gathered_myList_e[gp->displ_e[r] + i] - 1;
                memcpy(global + (size_t)gid * stride,
                       recv + ((size_t)gp->displ_e[r] + i) * stride,
                       (size_t)stride * sizeof(int));
            }
        }
        free(recv); free(counts); free(displs);
    }
}

/* ------------------------------------------------------------------ */
/* Transpose helpers (rank-0 only — operate on global arrays)         */
/* ------------------------------------------------------------------ */
static void transpose_node_columns(const real_t *src, real_t *dst,
                                   int n, int nl, int n_layers)
{
    for (int i = 0; i < n; ++i) {
        for (int nz = 0; nz < n_layers; ++nz) {
            dst[(size_t)nz * (size_t)n + (size_t)i] = src[(size_t)i * (size_t)nl + (size_t)nz];
        }
    }
}

static void extract_uv_component(const real_t *uv, real_t *dst,
                                 int n_elem, int nl, int n_layers, int comp)
{
    for (int e = 0; e < n_elem; ++e) {
        for (int nz = 0; nz < n_layers; ++nz) {
            size_t k = (size_t)e * (size_t)nl * 2 + (size_t)nz * 2 + (size_t)comp;
            dst[(size_t)nz * (size_t)n_elem + (size_t)e] = uv[k];
        }
    }
}

void fesom_io_write_snapshot(const char                  *path,
                             int                          step_n,
                             real_t                       dt,
                             const fesom_calendar_t      *cal,
                             const struct fesom_mesh     *mesh,
                             const struct fesom_dyn      *dyn,
                             const struct fesom_tracers  *tracers,
                             const struct fesom_aux      *aux,
                             const struct fesom_ice      *ice,
                             struct fesom_partit         *partit)
{
    const int nl     = mesh->nl;
    const int n_lay  = nl - 1;
    const int nod2D  = mesh->nod2D;
    const int elem2D = mesh->elem2D;

    /* Build gather plan once (collective). */
    gather_plan gp;
    gather_plan_init(&gp, mesh, partit);
    MPI_Comm comm = partit->MPI_COMM_FESOM;
    int mype = partit->mype;

    /* ---- Allocate global buffers on rank 0 -------------------------- */
    real_t *g_T = NULL, *g_S = NULL, *g_eta = NULL;
    real_t *g_w = NULL, *g_uv = NULL;
    real_t *g_dens = NULL, *g_bv = NULL, *g_pgfx = NULL, *g_pgfy = NULL;
    real_t *g_Kv = NULL, *g_Av = NULL;
    real_t *g_lon = NULL, *g_lat = NULL;
    int    *g_elem_nodes = NULL, *g_nlev_nod = NULL, *g_nlev_elem = NULL;
    /* sea-ice fields (NULL if ice == NULL) */
    real_t *g_aice = NULL, *g_mice = NULL, *g_msnow = NULL;
    real_t *g_uice = NULL, *g_vice = NULL;
    real_t *g_hice = NULL, *g_hsnow = NULL;
    if (mype == 0) {
        g_T    = malloc((size_t)nod2D  * (size_t)nl * sizeof(real_t));
        g_S    = malloc((size_t)nod2D  * (size_t)nl * sizeof(real_t));
        g_eta  = malloc((size_t)nod2D  * sizeof(real_t));
        g_w    = malloc((size_t)nod2D  * (size_t)nl * sizeof(real_t));
        g_uv   = malloc((size_t)elem2D * (size_t)nl * 2 * sizeof(real_t));
        g_dens = malloc((size_t)nod2D  * (size_t)nl * sizeof(real_t));
        g_bv   = malloc((size_t)nod2D  * (size_t)nl * sizeof(real_t));
        g_pgfx = malloc((size_t)elem2D * (size_t)nl * sizeof(real_t));
        g_pgfy = malloc((size_t)elem2D * (size_t)nl * sizeof(real_t));
        g_Kv   = malloc((size_t)nod2D  * (size_t)nl * sizeof(real_t));
        g_Av   = malloc((size_t)elem2D * (size_t)nl * sizeof(real_t));
        g_lon  = malloc((size_t)nod2D  * sizeof(real_t));
        g_lat  = malloc((size_t)nod2D  * sizeof(real_t));
        g_elem_nodes = malloc((size_t)elem2D * 3 * sizeof(int));
        g_nlev_nod   = malloc((size_t)nod2D  * sizeof(int));
        g_nlev_elem  = malloc((size_t)elem2D * sizeof(int));
        FESOM_CHECK(g_T && g_S && g_eta && g_w && g_uv && g_dens && g_bv &&
                    g_pgfx && g_pgfy && g_Kv && g_Av && g_lon && g_lat &&
                    g_elem_nodes && g_nlev_nod && g_nlev_elem,
                    "io: oom (global buffers)");
        if (ice) {
            g_aice  = malloc((size_t)nod2D * sizeof(real_t));
            g_mice  = malloc((size_t)nod2D * sizeof(real_t));
            g_msnow = malloc((size_t)nod2D * sizeof(real_t));
            g_uice  = malloc((size_t)nod2D * sizeof(real_t));
            g_vice  = malloc((size_t)nod2D * sizeof(real_t));
            g_hice  = malloc((size_t)nod2D * sizeof(real_t));
            g_hsnow = malloc((size_t)nod2D * sizeof(real_t));
            FESOM_CHECK(g_aice && g_mice && g_msnow && g_uice && g_vice
                     && g_hice && g_hsnow, "io: oom (ice global buffers)");
        }
    }

    /* ---- Gather per-vertex coordinates (geo radians → degrees on rank 0) -- */
    {
        real_t *lon_loc = malloc((size_t)gp.myDim_n * sizeof(real_t));
        real_t *lat_loc = malloc((size_t)gp.myDim_n * sizeof(real_t));
        FESOM_CHECK(lon_loc && lat_loc, "io: oom (lon/lat local)");
        for (int n = 0; n < gp.myDim_n; ++n) {
            lon_loc[n] = mesh->geo_coord_nod2D[2*n + 0];
            lat_loc[n] = mesh->geo_coord_nod2D[2*n + 1];
        }
        gather_node(lon_loc, 1, &gp, g_lon, comm);
        gather_node(lat_loc, 1, &gp, g_lat, comm);
        free(lon_loc); free(lat_loc);
        if (mype == 0) {
            for (int n = 0; n < nod2D; ++n) {
                g_lon[n] /= FESOM_RAD;
                g_lat[n] /= FESOM_RAD;
            }
        }
    }

    /* ---- Gather elem_nodes: translate local→global node IDs first --- */
    {
        int *eln_loc = malloc((size_t)gp.myDim_e * 3 * sizeof(int));
        FESOM_CHECK(eln_loc, "io: oom (elem_nodes local)");
        for (int e = 0; e < gp.myDim_e; ++e) {
            for (int k = 0; k < 3; ++k) {
                int loc = mesh->elem_nodes[3*e + k];   /* local node id */
                /* Translate to GLOBAL node id (0-based). myList_nod2D is 1-based. */
                eln_loc[3*e + k] = (loc >= 0)
                                   ? (partit->myList_nod2D[loc] - 1)
                                   : -1;
            }
        }
        gather_elem_int(eln_loc, 3, &gp, g_elem_nodes, comm);
        free(eln_loc);
    }

    /* ---- Gather nlevels (per-node and per-element) ------------------ */
    gather_node_int(mesh->nlevels_nod2D, 1, &gp, g_nlev_nod,  comm);
    gather_elem_int(mesh->nlevels,       1, &gp, g_nlev_elem, comm);

    /* ---- Gather state fields --------------------------------------- */
    gather_node(tracers->data[FESOM_TRACER_T].values, nl, &gp, g_T,  comm);
    gather_node(tracers->data[FESOM_TRACER_S].values, nl, &gp, g_S,  comm);
    gather_node(dyn->eta_n,                            1, &gp, g_eta, comm);
    gather_node(dyn->w,                               nl, &gp, g_w,  comm);
    gather_elem(dyn->uv,                            nl*2, &gp, g_uv, comm);
    gather_node(aux->density_m_rho0,                  nl, &gp, g_dens, comm);
    gather_node(aux->bvfreq,                          nl, &gp, g_bv,   comm);
    gather_elem(aux->pgf_x,                           nl, &gp, g_pgfx, comm);
    gather_elem(aux->pgf_y,                           nl, &gp, g_pgfy, comm);
    gather_node(aux->Kv,                              nl, &gp, g_Kv,   comm);
    gather_elem(aux->Av,                              nl, &gp, g_Av,   comm);

    if (ice) {
        gather_node(ice->data[FESOM_ICE_AICE].values,  1, &gp, g_aice,  comm);
        gather_node(ice->data[FESOM_ICE_MICE].values,  1, &gp, g_mice,  comm);
        gather_node(ice->data[FESOM_ICE_MSNOW].values, 1, &gp, g_msnow, comm);
        gather_node(ice->uice,   1, &gp, g_uice,  comm);
        gather_node(ice->vice,   1, &gp, g_vice,  comm);
        gather_node(ice->h_ice,  1, &gp, g_hice,  comm);
        gather_node(ice->h_snow, 1, &gp, g_hsnow, comm);
    }

    /* Non-rank-0 ranks are done. */
    if (mype != 0) {
        gather_plan_free(&gp);
        return;
    }

    /* ====================================================================
     * Rank 0: write the global file with time as UNLIMITED dim, length 1.
     * ==================================================================== */
    int ncid;
    NC_CHECK(nc_create(path, NC_CLOBBER | NC_NETCDF4, &ncid));

    int dim_time, dim_nod, dim_elem, dim_nz, dim_nz_1, dim_three;
    NC_CHECK(nc_def_dim(ncid, "time",  NC_UNLIMITED,    &dim_time));
    NC_CHECK(nc_def_dim(ncid, "nod2",  (size_t)nod2D,   &dim_nod));
    NC_CHECK(nc_def_dim(ncid, "elem",  (size_t)elem2D,  &dim_elem));
    NC_CHECK(nc_def_dim(ncid, "nz",    (size_t)nl,      &dim_nz));
    NC_CHECK(nc_def_dim(ncid, "nz_1",  (size_t)n_lay,   &dim_nz_1));
    NC_CHECK(nc_def_dim(ncid, "three", 3,               &dim_three));

    int var_time, var_lon, var_lat, var_zbar, var_Z;
    int var_elnodes, var_nlev_nod, var_nlev_elem;
    int var_T, var_S, var_eta, var_w, var_u, var_v;
    int var_dens, var_bv, var_pgfx, var_pgfy, var_Kv, var_Av;

    int dims_time[1]      = { dim_time };
    int dims_nod[1]       = { dim_nod };
    int dims_elem[1]      = { dim_elem };
    int dims_nz[1]        = { dim_nz };
    int dims_nz1[1]       = { dim_nz_1 };
    int dims_eln[2]       = { dim_elem, dim_three };
    /* time-varying layouts: outermost = time. */
    int dims_t_T[3]       = { dim_time, dim_nz_1, dim_nod };
    int dims_t_w[3]       = { dim_time, dim_nz,   dim_nod };
    int dims_t_uv[3]      = { dim_time, dim_nz_1, dim_elem };
    int dims_t_eta[2]     = { dim_time, dim_nod };

    NC_CHECK(nc_def_var(ncid, "time",           NC_DOUBLE, 1, dims_time, &var_time));
    NC_CHECK(nc_def_var(ncid, "lon",            NC_DOUBLE, 1, dims_nod,  &var_lon));
    NC_CHECK(nc_def_var(ncid, "lat",            NC_DOUBLE, 1, dims_nod,  &var_lat));
    NC_CHECK(nc_def_var(ncid, "zbar",           NC_DOUBLE, 1, dims_nz,   &var_zbar));
    NC_CHECK(nc_def_var(ncid, "Z",              NC_DOUBLE, 1, dims_nz1,  &var_Z));
    NC_CHECK(nc_def_var(ncid, "elem_nodes",     NC_INT,    2, dims_eln,  &var_elnodes));
    NC_CHECK(nc_def_var(ncid, "nlevels_nod2D",  NC_INT,    1, dims_nod,  &var_nlev_nod));
    NC_CHECK(nc_def_var(ncid, "nlevels",        NC_INT,    1, dims_elem, &var_nlev_elem));

    NC_CHECK(nc_def_var(ncid, "T",              NC_DOUBLE, 3, dims_t_T,   &var_T));
    NC_CHECK(nc_def_var(ncid, "S",              NC_DOUBLE, 3, dims_t_T,   &var_S));
    NC_CHECK(nc_def_var(ncid, "eta_n",          NC_DOUBLE, 2, dims_t_eta, &var_eta));
    NC_CHECK(nc_def_var(ncid, "w",              NC_DOUBLE, 3, dims_t_w,   &var_w));
    NC_CHECK(nc_def_var(ncid, "u",              NC_DOUBLE, 3, dims_t_uv,  &var_u));
    NC_CHECK(nc_def_var(ncid, "v",              NC_DOUBLE, 3, dims_t_uv,  &var_v));
    NC_CHECK(nc_def_var(ncid, "density_m_rho0", NC_DOUBLE, 3, dims_t_T,   &var_dens));
    NC_CHECK(nc_def_var(ncid, "bvfreq",         NC_DOUBLE, 3, dims_t_T,   &var_bv));
    NC_CHECK(nc_def_var(ncid, "pgf_x",          NC_DOUBLE, 3, dims_t_uv,  &var_pgfx));
    NC_CHECK(nc_def_var(ncid, "pgf_y",          NC_DOUBLE, 3, dims_t_uv,  &var_pgfy));
    NC_CHECK(nc_def_var(ncid, "Kv",             NC_DOUBLE, 3, dims_t_T,   &var_Kv));
    NC_CHECK(nc_def_var(ncid, "Av",             NC_DOUBLE, 3, dims_t_uv,  &var_Av));

    int var_aice = -1, var_mice = -1, var_msnow = -1;
    int var_uice = -1, var_vice = -1, var_hice = -1, var_hsnow = -1;
    if (ice) {
        NC_CHECK(nc_def_var(ncid, "a_ice",  NC_DOUBLE, 2, dims_t_eta, &var_aice));
        NC_CHECK(nc_def_var(ncid, "m_ice",  NC_DOUBLE, 2, dims_t_eta, &var_mice));
        NC_CHECK(nc_def_var(ncid, "m_snow", NC_DOUBLE, 2, dims_t_eta, &var_msnow));
        NC_CHECK(nc_def_var(ncid, "uice",   NC_DOUBLE, 2, dims_t_eta, &var_uice));
        NC_CHECK(nc_def_var(ncid, "vice",   NC_DOUBLE, 2, dims_t_eta, &var_vice));
        NC_CHECK(nc_def_var(ncid, "h_ice",  NC_DOUBLE, 2, dims_t_eta, &var_hice));
        NC_CHECK(nc_def_var(ncid, "h_snow", NC_DOUBLE, 2, dims_t_eta, &var_hsnow));
    }

    /* CF-compliant time stamp anchored to the model calendar. Replaces the
     * earlier "1900-01-01" placeholder so xarray.open_mfdataset on a stack
     * of snapshots produces dates in the run's actual year (e.g. 1958). */
    {
        char units_buf[64];
        fesom_calendar_units_string(cal, units_buf, sizeof units_buf);
        const char *cal_name = fesom_calendar_cf_name(cal->kind);
        NC_CHECK(nc_put_att_text(ncid, var_time, "units",     strlen(units_buf), units_buf));
        NC_CHECK(nc_put_att_text(ncid, var_time, "calendar",  strlen(cal_name),  cal_name));
        NC_CHECK(nc_put_att_text(ncid, var_time, "long_name", 10, "model time"));
    }
    NC_CHECK(nc_put_att_text(ncid, var_lon,  "units", 12, "degrees_east"));
    NC_CHECK(nc_put_att_text(ncid, var_lat,  "units", 13, "degrees_north"));
    NC_CHECK(nc_put_att_text(ncid, var_zbar, "units", 1, "m"));
    NC_CHECK(nc_put_att_text(ncid, var_Z,    "units", 1, "m"));
    NC_CHECK(nc_put_att_text(ncid, var_T,    "units", 9, "degrees_C"));
    NC_CHECK(nc_put_att_text(ncid, var_S,    "units", 3, "PSU"));
    NC_CHECK(nc_put_att_text(ncid, var_eta,  "units", 1, "m"));
    NC_CHECK(nc_put_att_text(ncid, var_w,    "units", 3, "m/s"));
    NC_CHECK(nc_put_att_text(ncid, var_u,    "units", 3, "m/s"));
    NC_CHECK(nc_put_att_text(ncid, var_v,    "units", 3, "m/s"));
    NC_CHECK(nc_put_att_text(ncid, var_dens, "units", 6, "kg/m^3"));
    NC_CHECK(nc_put_att_text(ncid, var_bv,   "units", 4, "s^-2"));
    NC_CHECK(nc_put_att_text(ncid, var_pgfx, "units", 5, "m/s^2"));
    NC_CHECK(nc_put_att_text(ncid, var_pgfy, "units", 5, "m/s^2"));
    NC_CHECK(nc_put_att_text(ncid, var_Kv,   "units", 5, "m^2/s"));
    NC_CHECK(nc_put_att_text(ncid, var_Av,   "units", 5, "m^2/s"));

    if (ice) {
        NC_CHECK(nc_put_att_text(ncid, var_aice,  "units", 1, "1"));
        NC_CHECK(nc_put_att_text(ncid, var_mice,  "units", 1, "m"));
        NC_CHECK(nc_put_att_text(ncid, var_msnow, "units", 1, "m"));
        NC_CHECK(nc_put_att_text(ncid, var_uice,  "units", 3, "m/s"));
        NC_CHECK(nc_put_att_text(ncid, var_vice,  "units", 3, "m/s"));
        NC_CHECK(nc_put_att_text(ncid, var_hice,  "units", 1, "m"));
        NC_CHECK(nc_put_att_text(ncid, var_hsnow, "units", 1, "m"));
    }

    int    step_attr = step_n;
    double dt_attr   = (double)dt;
    NC_CHECK(nc_put_att_int   (ncid, NC_GLOBAL, "step", NC_INT,    1, &step_attr));
    NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "dt",   NC_DOUBLE, 1, &dt_attr));

    NC_CHECK(nc_enddef(ncid));

    /* ---- time scalar (single record) -------------------------------
     * Pulled from the model calendar so the value lines up with the
     * "seconds since YYYY-MM-DD" units written above. Equals
     * step_n * dt for runs starting at YYYY-01-01 00:00:00, which is
     * what fesom_io_init enforces. step_n / dt remain available for the
     * "step" / "dt" global attributes set just above. */
    {
        size_t start[1] = {0}, count[1] = {1};
        double t_val = fesom_calendar_seconds_since_origin(cal);
        NC_CHECK(nc_put_vara_double(ncid, var_time, start, count, &t_val));
    }

    /* ---- static fields --------------------------------------------- */
    {
        double *lon_d = malloc((size_t)nod2D * sizeof(double));
        double *lat_d = malloc((size_t)nod2D * sizeof(double));
        FESOM_CHECK(lon_d && lat_d, "io: oom (lon/lat double cast)");
        for (int n = 0; n < nod2D; ++n) {
            lon_d[n] = (double)g_lon[n];
            lat_d[n] = (double)g_lat[n];
        }
        NC_CHECK(nc_put_var_double(ncid, var_lon, lon_d));
        NC_CHECK(nc_put_var_double(ncid, var_lat, lat_d));
        free(lon_d); free(lat_d);
    }
    NC_CHECK(nc_put_var_double(ncid, var_zbar, mesh->zbar));
    NC_CHECK(nc_put_var_double(ncid, var_Z,    mesh->Z));
    NC_CHECK(nc_put_var_int   (ncid, var_elnodes,   g_elem_nodes));
    NC_CHECK(nc_put_var_int   (ncid, var_nlev_nod,  g_nlev_nod));
    NC_CHECK(nc_put_var_int   (ncid, var_nlev_elem, g_nlev_elem));

    /* ---- transposed time-varying state at time index 0 ------------- */
    real_t *buf_nl   = malloc((size_t)n_lay  * (size_t)nod2D  * sizeof(real_t));
    real_t *buf_nl_e = malloc((size_t)n_lay  * (size_t)elem2D * sizeof(real_t));
    real_t *buf_w    = malloc((size_t)nl     * (size_t)nod2D  * sizeof(real_t));
    FESOM_CHECK(buf_nl && buf_nl_e && buf_w, "io: oom (transpose buffers)");

    /* Slab-write: start at time=0, count=1 along time, full extent for the rest. */
    size_t startT[3] = {0, 0, 0}, countT[3] = {1, (size_t)n_lay, (size_t)nod2D};
    size_t startW[3] = {0, 0, 0}, countW[3] = {1, (size_t)nl,    (size_t)nod2D};
    size_t startE[3] = {0, 0, 0}, countE[3] = {1, (size_t)n_lay, (size_t)elem2D};
    size_t start1[2] = {0, 0},    count1[2] = {1, (size_t)nod2D};

    transpose_node_columns(g_T,    buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_vara_double(ncid, var_T, startT, countT, buf_nl));
    transpose_node_columns(g_S,    buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_vara_double(ncid, var_S, startT, countT, buf_nl));

    NC_CHECK(nc_put_vara_double(ncid, var_eta, start1, count1, g_eta));

    transpose_node_columns(g_w,    buf_w, nod2D, nl, nl);
    NC_CHECK(nc_put_vara_double(ncid, var_w, startW, countW, buf_w));

    extract_uv_component(g_uv, buf_nl_e, elem2D, nl, n_lay, /*comp=*/0);
    NC_CHECK(nc_put_vara_double(ncid, var_u, startE, countE, buf_nl_e));
    extract_uv_component(g_uv, buf_nl_e, elem2D, nl, n_lay, /*comp=*/1);
    NC_CHECK(nc_put_vara_double(ncid, var_v, startE, countE, buf_nl_e));

    transpose_node_columns(g_dens, buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_vara_double(ncid, var_dens, startT, countT, buf_nl));
    transpose_node_columns(g_bv,   buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_vara_double(ncid, var_bv,   startT, countT, buf_nl));
    transpose_node_columns(g_Kv,   buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_vara_double(ncid, var_Kv,   startT, countT, buf_nl));

    /* Element-column transpose for pgf and Av: data already in [elem][nl] row-major. */
    for (int e = 0; e < elem2D; ++e) for (int nz = 0; nz < n_lay; ++nz)
        buf_nl_e[(size_t)nz * elem2D + e] = g_pgfx[(size_t)e * nl + nz];
    NC_CHECK(nc_put_vara_double(ncid, var_pgfx, startE, countE, buf_nl_e));
    for (int e = 0; e < elem2D; ++e) for (int nz = 0; nz < n_lay; ++nz)
        buf_nl_e[(size_t)nz * elem2D + e] = g_pgfy[(size_t)e * nl + nz];
    NC_CHECK(nc_put_vara_double(ncid, var_pgfy, startE, countE, buf_nl_e));
    for (int e = 0; e < elem2D; ++e) for (int nz = 0; nz < n_lay; ++nz)
        buf_nl_e[(size_t)nz * elem2D + e] = g_Av[(size_t)e * nl + nz];
    NC_CHECK(nc_put_vara_double(ncid, var_Av,   startE, countE, buf_nl_e));

    free(buf_nl); free(buf_nl_e); free(buf_w);

    if (ice) {
        NC_CHECK(nc_put_vara_double(ncid, var_aice,  start1, count1, g_aice));
        NC_CHECK(nc_put_vara_double(ncid, var_mice,  start1, count1, g_mice));
        NC_CHECK(nc_put_vara_double(ncid, var_msnow, start1, count1, g_msnow));
        NC_CHECK(nc_put_vara_double(ncid, var_uice,  start1, count1, g_uice));
        NC_CHECK(nc_put_vara_double(ncid, var_vice,  start1, count1, g_vice));
        NC_CHECK(nc_put_vara_double(ncid, var_hice,  start1, count1, g_hice));
        NC_CHECK(nc_put_vara_double(ncid, var_hsnow, start1, count1, g_hsnow));
    }

    NC_CHECK(nc_close(ncid));

    free(g_T); free(g_S); free(g_eta); free(g_w); free(g_uv);
    free(g_dens); free(g_bv); free(g_pgfx); free(g_pgfy);
    free(g_Kv); free(g_Av);
    free(g_lon); free(g_lat);
    free(g_elem_nodes); free(g_nlev_nod); free(g_nlev_elem);
    free(g_aice); free(g_mice); free(g_msnow);
    free(g_uice); free(g_vice); free(g_hice); free(g_hsnow);

    gather_plan_free(&gp);
}

/* ====================================================================== */
/* Task 6 — orchestrator + default monthly variable table.               */
/*                                                                        */
/* Each resolver ADDS its local-slice contribution to the accumulator    */
/* buffer. For 1:1 fields the loop is just `out[i] += src[i]`. Synthesised */
/* fields (sst from T-surface, u/v from uv components) do their own       */
/* indexing.                                                             */
/* ====================================================================== */

static void resolve_temp(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->tracers->data[FESOM_TRACER_T].values;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_salt(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->tracers->data[FESOM_TRACER_S].values;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_sst(const fesom_state *s, real_t *out, size_t n)
{
    /* Surface slice of T tracer; tracer storage stride is `nl`. */
    const real_t *T = s->tracers->data[FESOM_TRACER_T].values;
    const int nl = s->mesh->nl;
    for (size_t i = 0; i < n; ++i) out[i] += T[i * (size_t)nl + 0];
}
static void resolve_sss(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *S = s->tracers->data[FESOM_TRACER_S].values;
    const int nl = s->mesh->nl;
    for (size_t i = 0; i < n; ++i) out[i] += S[i * (size_t)nl + 0];
}
static void resolve_ssh(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->dyn->eta_n;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_u(const fesom_state *s, real_t *out, size_t n)
{
    /* uv layout: [elem][nl][2]. 3D_ELEM_MID accumulator: [elem][nl].
     * Halo-exchanged at fesom_step.c:127 — replicated-element values
     * agree across ranks at the moment of accumulation, so local-only
     * sum stays consistent under MPI. */
    const real_t *uv = s->dyn->uv;
    const int nl = s->mesh->nl;
    const size_t E = n / (size_t)nl;
    for (size_t e = 0; e < E; ++e) {
        for (int k = 0; k < nl; ++k) {
            out[e * (size_t)nl + (size_t)k] += uv[e * (size_t)nl * 2 + (size_t)k * 2 + 0];
        }
    }
}
static void resolve_v(const fesom_state *s, real_t *out, size_t n)
{
    /* Same exchange-status note as resolve_u; both components share uv. */
    const real_t *uv = s->dyn->uv;
    const int nl = s->mesh->nl;
    const size_t E = n / (size_t)nl;
    for (size_t e = 0; e < E; ++e) {
        for (int k = 0; k < nl; ++k) {
            out[e * (size_t)nl + (size_t)k] += uv[e * (size_t)nl * 2 + (size_t)k * 2 + 1];
        }
    }
}
static void resolve_w(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->dyn->w;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_Kv(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->aux->Kv;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_Av(const fesom_state *s, real_t *out, size_t n)
{
    /* aux->Av is element-replicated; halo-exchanged at fesom_step.c:92,96.
     * Local-only mean accumulation is correct (see resolve_u banner). */
    const real_t *src = s->aux->Av;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_density(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->aux->density_m_rho0;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_bvfreq(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->aux->bvfreq;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_a_ice(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->ice->data[FESOM_ICE_AICE].values;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_m_ice(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->ice->data[FESOM_ICE_MICE].values;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_m_snow(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->ice->data[FESOM_ICE_MSNOW].values;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_uice(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->ice->uice;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}
static void resolve_vice(const fesom_state *s, real_t *out, size_t n)
{
    const real_t *src = s->ice->vice;
    for (size_t i = 0; i < n; ++i) out[i] += src[i];
}

/* Default monthly variable table. Lifetime = run.
 * If sea ice isn't initialised, callers should pass nvars = 12 instead
 * of FESOM_DEFAULT_MONTHLY_NVARS to skip the trailing 5 ice entries. */
#define FESOM_DEFAULT_MONTHLY_NVARS 17

static const fesom_var_desc_t fesom_default_monthly_table[FESOM_DEFAULT_MONTHLY_NVARS] = {
    /* 3D ocean state */
    { "temp",    "potential temperature",          "degC",   FESOM_VAR_3D_NODE_MID,   resolve_temp    },
    { "salt",    "salinity",                       "PSU",    FESOM_VAR_3D_NODE_MID,   resolve_salt    },
    { "ssh",     "sea surface height",             "m",      FESOM_VAR_2D_NODE,       resolve_ssh     },
    { "sst",     "sea surface temperature",        "degC",   FESOM_VAR_2D_NODE,       resolve_sst     },
    { "sss",     "sea surface salinity",           "PSU",    FESOM_VAR_2D_NODE,       resolve_sss     },
    { "u",       "zonal velocity",                 "m/s",    FESOM_VAR_3D_ELEM_MID,   resolve_u       },
    { "v",       "meridional velocity",            "m/s",    FESOM_VAR_3D_ELEM_MID,   resolve_v       },
    { "w",       "vertical velocity at interfaces","m/s",    FESOM_VAR_3D_NODE_IFACE, resolve_w       },
    /* Kv, bvfreq (N^2) and Av are INTERFACE quantities (defined at the nl level
     * interfaces, like w), so they must be written on nz (nl=48) — matching the
     * Fortran. They were wrongly registered as MID (nz_1=47), which dropped the
     * bottom interface and mislabelled the depths as layer mid-points. Kv/bvfreq
     * are node-based (NODE_IFACE); Av is element-based (ELEM_IFACE). */
    { "Kv",      "vertical diffusivity",           "m^2/s",  FESOM_VAR_3D_NODE_IFACE, resolve_Kv      },
    { "Av",      "vertical viscosity",             "m^2/s",  FESOM_VAR_3D_ELEM_IFACE, resolve_Av      },
    { "density", "in-situ density minus rho0",     "kg/m^3", FESOM_VAR_3D_NODE_MID,   resolve_density },
    { "bvfreq",  "Brunt-Vaisala frequency squared","s^-2",   FESOM_VAR_3D_NODE_IFACE, resolve_bvfreq  },
    /* Ice fields (skipped when state->ice == NULL — see fesom_io_init) */
    { "a_ice",   "ice area fraction",              "1",      FESOM_VAR_2D_NODE,       resolve_a_ice   },
    { "m_ice",   "ice volume per area",            "m",      FESOM_VAR_2D_NODE,       resolve_m_ice   },
    { "m_snow",  "snow volume per area",           "m",      FESOM_VAR_2D_NODE,       resolve_m_snow  },
    { "uice",    "ice zonal velocity",             "m/s",    FESOM_VAR_2D_NODE,       resolve_uice    },
    { "vice",    "ice meridional velocity",        "m/s",    FESOM_VAR_2D_NODE,       resolve_vice    },
};

/* ====================================================================== */
/* Orchestrator                                                          */
/* ====================================================================== */

/* Default cadence-kind mapping: STEP -> INSTANT, others -> MEAN. */
static fesom_output_kind_t default_kind_for_period(fesom_period_kind_t p)
{
    return (p == FESOM_PERIOD_STEP) ? FESOM_OUT_INSTANT : FESOM_OUT_MEAN;
}

void fesom_io_init(fesom_io_t                  *io,
                   const char                  *out_dir,
                   fesom_calendar_kind_t        cal_kind,
                   int                          year,
                   int                          month,
                   int                          day,
                   const char                  *config_file_or_null,
                   const struct fesom_mesh     *mesh,
                   struct fesom_partit         *partit)
{
    memset(io, 0, sizeof(*io));
    fesom_calendar_init(&io->calendar, cal_kind, year, month, day);
    io->prev_calendar = io->calendar;

    /* Heap-owned out_dir copy. */
    {
        size_t n = strlen(out_dir) + 1;
        io->out_dir = malloc(n);
        FESOM_CHECK(io->out_dir, "fesom_io_init: oom (out_dir)");
        memcpy(io->out_dir, out_dir, n);
    }

    /* Optional override file. If parsing fails, abort — better to fail
     * loud than to silently produce unexpected output. */
    fesom_io_config_t cfg;
    int have_cfg = 0;
    if (config_file_or_null && config_file_or_null[0]) {
        if (fesom_io_config_parse(config_file_or_null, &cfg) != 0) {
            FESOM_DIE("fesom_io_init: failed to parse %s", config_file_or_null);
        }
        have_cfg = 1;
        if (partit->mype == 0) {
            printf("[fesom_io] config %s parsed: %d entries\n",
                   config_file_or_null, cfg.n_entries);
        }
    }

    /* For each variable in the default table, decide its cadence list:
     *   - if config has an entry: use its cadences (REPLACE the default)
     *   - else: keep the default [MONTHLY]
     * Then bin the variables per cadence so each stream gets its var list. */
    int per_cad_count[5] = {0};
    /* First pass: count vars per cadence to size buffers. */
    for (int v = 0; v < FESOM_DEFAULT_MONTHLY_NVARS; ++v) {
        const fesom_io_config_entry_t *e =
            have_cfg ? fesom_io_config_lookup(&cfg, fesom_default_monthly_table[v].name) : NULL;
        if (e) {
            for (int c = 0; c < e->n_cadences; ++c) per_cad_count[e->cadences[c]] += 1;
        } else {
            per_cad_count[FESOM_PERIOD_MONTHLY] += 1;
        }
    }
    /* Allocate the per-cadence var arrays. */
    int per_cad_idx[5] = {0};
    for (int p = 0; p < 5; ++p) {
        io->owned_nvars[p] = per_cad_count[p];
        if (per_cad_count[p] > 0) {
            io->owned_vars[p] = malloc((size_t)per_cad_count[p] * sizeof(fesom_var_desc_t));
            FESOM_CHECK(io->owned_vars[p], "fesom_io_init: oom (cadence %d)", p);
        }
    }
    /* Second pass: fill them. Var descriptors are copied from the static
     * table; resolver function pointers stay valid for the run. */
    for (int v = 0; v < FESOM_DEFAULT_MONTHLY_NVARS; ++v) {
        const fesom_io_config_entry_t *e =
            have_cfg ? fesom_io_config_lookup(&cfg, fesom_default_monthly_table[v].name) : NULL;
        if (e) {
            for (int c = 0; c < e->n_cadences; ++c) {
                fesom_period_kind_t p = e->cadences[c];
                io->owned_vars[p][per_cad_idx[p]++] = fesom_default_monthly_table[v];
            }
        } else {
            fesom_period_kind_t p = FESOM_PERIOD_MONTHLY;
            io->owned_vars[p][per_cad_idx[p]++] = fesom_default_monthly_table[v];
        }
    }
    if (have_cfg) fesom_io_config_free(&cfg);

    /* Init each cadence's stream. Empty cadences stay inactive. */
    for (int p = 0; p < 5; ++p) {
        if (io->owned_nvars[p] == 0) continue;
        fesom_io_stream_init(&io->streams[p],
                             (fesom_period_kind_t)p,
                             default_kind_for_period((fesom_period_kind_t)p),
                             io->owned_vars[p],
                             io->owned_nvars[p],
                             mesh, partit, io->out_dir);
        io->stream_active[p] = 1;
    }

    if (partit->mype == 0) {
        const char *names[5] = { "step", "hourly", "daily", "monthly", "yearly" };
        printf("[fesom_io] initialised: out_dir=%s, calendar=%s @ %04d-%02d-%02d\n",
               io->out_dir, fesom_calendar_cf_name(cal_kind), year, month, day);
        for (int p = 0; p < 5; ++p) {
            if (!io->stream_active[p]) continue;
            printf("[fesom_io]   %-7s : %d vars\n", names[p], io->owned_nvars[p]);
        }
    }
}

void fesom_io_step(fesom_io_t              *io,
                   double                   dt,
                   const fesom_state       *state,
                   struct fesom_partit     *partit)
{
    /* Save previous-step calendar BEFORE advancing, so streams can detect
     * period crossings via fesom_calendar_crossed(prev, curr, period). */
    io->prev_calendar = io->calendar;
    fesom_calendar_advance(&io->calendar, dt);

    /* Drive each active stream. (At Task 6, only MONTHLY is active; the
     * other slots are noops.) */
    for (int p = 0; p < 5; ++p) {
        if (!io->stream_active[p]) continue;
        fesom_io_stream_step(&io->streams[p],
                             state,
                             &io->prev_calendar,
                             &io->calendar,
                             partit);
    }
}

void fesom_io_finalize(fesom_io_t          *io,
                       struct fesom_partit *partit)
{
    for (int p = 0; p < 5; ++p) {
        if (io->stream_active[p]) {
            fesom_io_stream_finalize(&io->streams[p], &io->calendar, partit);
            fesom_io_stream_close(&io->streams[p]);
            io->stream_active[p] = 0;
        }
        free(io->owned_vars[p]);
        io->owned_vars[p]  = NULL;
        io->owned_nvars[p] = 0;
    }
    free(io->out_dir);
    io->out_dir = NULL;
}
