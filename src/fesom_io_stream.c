/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 * See fesom_io_stream.h banner for context.
 *
 * fesom_io_stream.c - per-cadence accumulator + writer. One stream owns N
 * variables sharing a calendar period (monthly, daily, ...) and writes one
 * netCDF file per (variable, year, cadence). The hot path (per-step) is a
 * local sum loop with no MPI; cold path (per-period flush) gathers via the
 * cached gather_plan and appends a record to the open file.
 *
 * Task 5 scope: MEAN path + MONTHLY (per-year rollover) only. INSTANT path
 * and other cadences/rollovers wired in Task 7.
 */
#include "fesom_io_stream.h"
#include "fesom_io.h"          /* full fesom_state typedef */
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_ice_types.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NC_CHECK(call) do {                                               \
    int _ec = (call);                                                     \
    if (_ec != NC_NOERR) {                                                \
        FESOM_DIE("io_stream: NetCDF error: %s", nc_strerror(_ec));       \
    }                                                                     \
} while (0)

/* Size dispatchers live in fesom_io_stream_dispatch.c so test_io_stream_unit
 * can link them without dragging netCDF/MPI in. */

static int kind_is_node(fesom_var_kind_t kind)
{
    return kind == FESOM_VAR_2D_NODE
        || kind == FESOM_VAR_3D_NODE_MID
        || kind == FESOM_VAR_3D_NODE_IFACE;
}

/* ------------------------------------------------------------------ */
/* Cache global node + element coords (geographic degrees) on rank 0  */
/* by gathering each rank's local view. Done once at stream init; the */
/* result is reused on every file open during the run. The local     */
/* element centroid is computed FROM each rank's own elem_nodes (with */
/* local node indices) BEFORE gather, so we never read off the end of */
/* a local-sized array.                                               */
/* ------------------------------------------------------------------ */
static void build_coord_cache(fesom_io_stream_t *s)
{
    const struct fesom_mesh *mesh = s->mesh;
    const double rad2deg = 180.0 / FESOM_PI;
    const int my_n = mesh->myDim_nod2D;
    const int my_e = mesh->myDim_elem2D;

    /* Local node coords in degrees (interior only). */
    real_t *loc_lon_n = malloc((size_t)my_n * sizeof(real_t));
    real_t *loc_lat_n = malloc((size_t)my_n * sizeof(real_t));
    FESOM_CHECK(loc_lon_n && loc_lat_n, "io_stream: oom (local node coord)");
    for (int n = 0; n < my_n; ++n) {
        loc_lon_n[n] = (real_t)((double)mesh->geo_coord_nod2D[2 * n + 0] * rad2deg);
        loc_lat_n[n] = (real_t)((double)mesh->geo_coord_nod2D[2 * n + 1] * rad2deg);
    }

    /* Local element centroids in degrees (interior only). Each elem's
     * 3 vertex IDs are LOCAL indices into geo_coord_nod2D, which is
     * sized myDim+eDim, so this access is always in-bounds. */
    real_t *loc_lon_e = malloc((size_t)my_e * sizeof(real_t));
    real_t *loc_lat_e = malloc((size_t)my_e * sizeof(real_t));
    FESOM_CHECK(loc_lon_e && loc_lat_e, "io_stream: oom (local elem centroid)");
    for (int e = 0; e < my_e; ++e) {
        double slon = 0.0, slat = 0.0;
        for (int k = 0; k < 3; ++k) {
            int n = mesh->elem_nodes[3 * e + k];
            slon += (double)mesh->geo_coord_nod2D[2 * n + 0];
            slat += (double)mesh->geo_coord_nod2D[2 * n + 1];
        }
        loc_lon_e[e] = (real_t)((slon / 3.0) * rad2deg);
        loc_lat_e[e] = (real_t)((slat / 3.0) * rad2deg);
    }

    /* Allocate global on rank 0 and gather. */
    if (s->mype == 0) {
        s->cached_lon_node = malloc((size_t)mesh->nod2D * sizeof(double));
        s->cached_lat_node = malloc((size_t)mesh->nod2D * sizeof(double));
        s->cached_lon_elem = malloc((size_t)mesh->elem2D * sizeof(double));
        s->cached_lat_elem = malloc((size_t)mesh->elem2D * sizeof(double));
        FESOM_CHECK(s->cached_lon_node && s->cached_lat_node
                    && s->cached_lon_elem && s->cached_lat_elem,
                    "io_stream: oom (cached coord)");
    }

    /* gather_node/gather_elem operate on real_t buffers; pack into
     * the global double caches via a temporary real_t global. */
    real_t *tmp_n = NULL, *tmp_e = NULL;
    if (s->mype == 0) {
        tmp_n = malloc((size_t)mesh->nod2D  * sizeof(real_t));
        tmp_e = malloc((size_t)mesh->elem2D * sizeof(real_t));
        FESOM_CHECK(tmp_n && tmp_e, "io_stream: oom (gather tmp)");
    }

    gather_node(loc_lon_n, 1, &s->gp, tmp_n, s->comm);
    if (s->mype == 0)
        for (int n = 0; n < mesh->nod2D; ++n) s->cached_lon_node[n] = (double)tmp_n[n];
    gather_node(loc_lat_n, 1, &s->gp, tmp_n, s->comm);
    if (s->mype == 0)
        for (int n = 0; n < mesh->nod2D; ++n) s->cached_lat_node[n] = (double)tmp_n[n];

    gather_elem(loc_lon_e, 1, &s->gp, tmp_e, s->comm);
    if (s->mype == 0)
        for (int e = 0; e < mesh->elem2D; ++e) s->cached_lon_elem[e] = (double)tmp_e[e];
    gather_elem(loc_lat_e, 1, &s->gp, tmp_e, s->comm);
    if (s->mype == 0)
        for (int e = 0; e < mesh->elem2D; ++e) s->cached_lat_elem[e] = (double)tmp_e[e];

    free(loc_lon_n); free(loc_lat_n);
    free(loc_lon_e); free(loc_lat_e);
    free(tmp_n); free(tmp_e);
}

/* ------------------------------------------------------------------ */
/* File naming + rollover-window keying                               */
/* ------------------------------------------------------------------ */

static void file_path_for(char *buf, size_t buflen,
                          const char *out_dir,
                          const char *var_name,
                          fesom_period_kind_t period,
                          int year, int month, int day)
{
    const char *cad = NULL;
    switch (period) {
        case FESOM_PERIOD_STEP:    cad = "step";    break;
        case FESOM_PERIOD_HOURLY:  cad = "hourly";  break;
        case FESOM_PERIOD_DAILY:   cad = "daily";   break;
        case FESOM_PERIOD_MONTHLY: cad = "monthly"; break;
        case FESOM_PERIOD_YEARLY:  cad = "yearly";  break;
    }
    /* Rollover defaults: STEP per-day, HOURLY per-month, others per-year. */
    if (period == FESOM_PERIOD_STEP) {
        snprintf(buf, buflen, "%s/%s.fesom.%04d-%02d-%02d.%s.nc",
                 out_dir, var_name, year, month, day, cad);
    } else if (period == FESOM_PERIOD_HOURLY) {
        snprintf(buf, buflen, "%s/%s.fesom.%04d-%02d.%s.nc",
                 out_dir, var_name, year, month, cad);
    } else {
        snprintf(buf, buflen, "%s/%s.fesom.%04d.%s.nc",
                 out_dir, var_name, year, cad);
    }
}

/* Window key: returns 1 if the file currently open for variable v
 * matches the period containing `cal`. */
static int window_matches(const fesom_io_stream_t *s, const fesom_calendar_t *cal)
{
    if (s->open_year < 0) return 0;
    if (s->period == FESOM_PERIOD_STEP) {
        return s->open_year  == cal->year
            && s->open_month == cal->month
            && s->open_day   == cal->day;
    }
    if (s->period == FESOM_PERIOD_HOURLY) {
        return s->open_year  == cal->year
            && s->open_month == cal->month;
    }
    return s->open_year == cal->year;
}

static void window_set(fesom_io_stream_t *s, const fesom_calendar_t *cal)
{
    s->open_year  = cal->year;
    s->open_month = (s->period == FESOM_PERIOD_STEP || s->period == FESOM_PERIOD_HOURLY)
                  ? cal->month : 0;
    s->open_day   = (s->period == FESOM_PERIOD_STEP) ? cal->day : 0;
}

/* ------------------------------------------------------------------ */
/* Open one variable's netCDF file for the current rollover window    */
/* ------------------------------------------------------------------ */

static void open_var_file(fesom_io_stream_t *s, int v, const fesom_calendar_t *cal)
{
    if (s->mype != 0) {
        s->ncid[v] = -1;       /* non-zero ranks never write */
        return;
    }
    char path[1024];
    file_path_for(path, sizeof path, s->out_dir, s->vars[v].name,
                  s->period, cal->year, cal->month, cal->day);

    int ncid;
    NC_CHECK(nc_create(path, NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL, &ncid));

    const struct fesom_mesh *mesh = s->mesh;
    const fesom_var_kind_t kind = s->vars[v].kind;
    const int nl = mesh->nl;

    /* Dimensions ---------------------------------------------------- */
    int dim_time, dim_nbnds;
    NC_CHECK(nc_def_dim(ncid, "time",  NC_UNLIMITED, &dim_time));
    NC_CHECK(nc_def_dim(ncid, "nbnds", 2,            &dim_nbnds));

    int dim_space = -1, dim_layers = -1, dim_two = -1;
    int n_layers = fesom_io_stream_n_layers_written(kind, nl);
    if (kind_is_node(kind)) {
        NC_CHECK(nc_def_dim(ncid, "nod2", mesh->nod2D, &dim_space));
    } else {
        NC_CHECK(nc_def_dim(ncid, "elem", mesh->elem2D, &dim_space));
    }
    if (n_layers > 0) {
        const char *lname = (kind == FESOM_VAR_3D_NODE_IFACE
                          || kind == FESOM_VAR_3D_ELEM_IFACE) ? "nz" : "nz_1";
        NC_CHECK(nc_def_dim(ncid, lname, n_layers, &dim_layers));
    }
    if (kind == FESOM_VAR_2D_ELEM_VEC) {
        NC_CHECK(nc_def_dim(ncid, "comp", 2, &dim_two));
    }

    /* time / time_bnds variables ------------------------------------ */
    int var_time, var_bnds;
    int time_dims[1] = {dim_time};
    int bnds_dims[2] = {dim_time, dim_nbnds};
    NC_CHECK(nc_def_var(ncid, "time",      NC_DOUBLE, 1, time_dims, &var_time));
    NC_CHECK(nc_def_var(ncid, "time_bnds", NC_DOUBLE, 2, bnds_dims, &var_bnds));

    char units_buf[64];
    fesom_calendar_units_string(cal, units_buf, sizeof units_buf);
    const char *cal_name = fesom_calendar_cf_name(cal->kind);
    NC_CHECK(nc_put_att_text(ncid, var_time, "units",     strlen(units_buf), units_buf));
    NC_CHECK(nc_put_att_text(ncid, var_time, "calendar",  strlen(cal_name),  cal_name));
    NC_CHECK(nc_put_att_text(ncid, var_time, "bounds",    9, "time_bnds"));
    NC_CHECK(nc_put_att_text(ncid, var_time, "long_name", 10, "model time"));

    /* Spatial coordinate variables ---------------------------------- */
    int var_lon, var_lat;
    int sp_dim[1] = {dim_space};
    if (kind_is_node(kind)) {
        NC_CHECK(nc_def_var(ncid, "lon", NC_DOUBLE, 1, sp_dim, &var_lon));
        NC_CHECK(nc_def_var(ncid, "lat", NC_DOUBLE, 1, sp_dim, &var_lat));
    } else {
        NC_CHECK(nc_def_var(ncid, "lon_elem", NC_DOUBLE, 1, sp_dim, &var_lon));
        NC_CHECK(nc_def_var(ncid, "lat_elem", NC_DOUBLE, 1, sp_dim, &var_lat));
    }
    NC_CHECK(nc_put_att_text(ncid, var_lon, "units", 12, "degrees_east"));
    NC_CHECK(nc_put_att_text(ncid, var_lat, "units", 13, "degrees_north"));

    /* Z-axis coordinate (only for 3D fields) ------------------------ */
    int var_z = -1;
    if (n_layers > 0) {
        int z_dim[1] = {dim_layers};
        const char *zname = (kind == FESOM_VAR_3D_NODE_IFACE
                          || kind == FESOM_VAR_3D_ELEM_IFACE) ? "zbar" : "Z";
        NC_CHECK(nc_def_var(ncid, zname, NC_DOUBLE, 1, z_dim, &var_z));
        NC_CHECK(nc_put_att_text(ncid, var_z, "units", 1, "m"));
        NC_CHECK(nc_put_att_text(ncid, var_z, "positive", 4, "down"));
    }

    /* Field variable ------------------------------------------------ */
    int var_data;
    int data_dims[4];
    int ndims = 0;
    data_dims[ndims++] = dim_time;
    if (n_layers > 0)            data_dims[ndims++] = dim_layers;
    if (kind == FESOM_VAR_2D_ELEM_VEC) data_dims[ndims++] = dim_two;
    data_dims[ndims++] = dim_space;
    NC_CHECK(nc_def_var(ncid, s->vars[v].name, NC_DOUBLE, ndims, data_dims, &var_data));
    if (s->vars[v].long_name)
        NC_CHECK(nc_put_att_text(ncid, var_data, "long_name",
                                 strlen(s->vars[v].long_name), s->vars[v].long_name));
    if (s->vars[v].units)
        NC_CHECK(nc_put_att_text(ncid, var_data, "units",
                                 strlen(s->vars[v].units), s->vars[v].units));
    {
        const char *cm = (s->output_kind == FESOM_OUT_MEAN) ? "time: mean" : "time: point";
        NC_CHECK(nc_put_att_text(ncid, var_data, "cell_methods", strlen(cm), cm));
    }

    /* Global attributes --------------------------------------------- */
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.8"));
    int origin[3] = {cal->origin_year, cal->origin_month, cal->origin_day};
    NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "time_origin", NC_INT, 3, origin));

    NC_CHECK(nc_enddef(ncid));

    /* Write static coords from the per-stream cache built at init. */
    if (kind_is_node(kind)) {
        NC_CHECK(nc_put_var_double(ncid, var_lon, s->cached_lon_node));
        NC_CHECK(nc_put_var_double(ncid, var_lat, s->cached_lat_node));
    } else {
        NC_CHECK(nc_put_var_double(ncid, var_lon, s->cached_lon_elem));
        NC_CHECK(nc_put_var_double(ncid, var_lat, s->cached_lat_elem));
    }
    if (n_layers > 0) {
        if (kind == FESOM_VAR_3D_NODE_IFACE || kind == FESOM_VAR_3D_ELEM_IFACE) {
            NC_CHECK(nc_put_var_double(ncid, var_z, mesh->zbar));
        } else {
            NC_CHECK(nc_put_var_double(ncid, var_z, mesh->Z));
        }
    }

    s->ncid[v]       = ncid;
    s->var_id[v]     = var_data;
    s->time_id[v]    = var_time;
    s->bnds_id[v]    = var_bnds;
    s->time_index[v] = 0;
}

static void close_open_files(fesom_io_stream_t *s)
{
    if (s->mype != 0) return;
    for (int v = 0; v < s->nvars; ++v) {
        if (s->ncid[v] >= 0) {
            NC_CHECK(nc_close(s->ncid[v]));
            s->ncid[v] = -1;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Flush: divide accum by n_accum, gather, write one record           */
/* ------------------------------------------------------------------ */

static void flush_one_var(fesom_io_stream_t *s, int v,
                          const fesom_calendar_t *cal,
                          const struct fesom_mesh *mesh)
{
    const fesom_var_kind_t kind = s->vars[v].kind;
    const int nl = mesh->nl;

    /* On non-zero ranks we still participate in the gather, but never
     * touch global buffers or the netCDF file. */
    int local_n = mesh->myDim_nod2D;
    int local_e = mesh->myDim_elem2D;
    size_t accum_n = s->accum_sz[v];
    real_t *accum = s->accum[v];

    /* Divide running sum by sample count to yield the mean. INSTANT
     * skips this branch (Task 7). */
    if (s->output_kind == FESOM_OUT_MEAN && s->n_accum > 0) {
        const real_t inv = (real_t)1.0 / (real_t)s->n_accum;
        for (size_t i = 0; i < accum_n; ++i) accum[i] *= inv;
    }

    /* Allocate global output (only meaningful on rank 0). */
    real_t *global = NULL;
    if (s->mype == 0) {
        size_t gsz = fesom_io_stream_global_size(kind, mesh->nod2D, mesh->elem2D, nl);
        global = malloc(gsz * sizeof(real_t));
        FESOM_CHECK(global, "io_stream: oom (global buf)");
    }

    /* Gather + write -------------------------------------------------*/
    int n_layers = fesom_io_stream_n_layers_written(kind, nl);
    switch (kind) {
        case FESOM_VAR_2D_NODE:
            gather_node(accum, 1, &s->gp, global, s->comm);
            break;
        case FESOM_VAR_3D_NODE_MID:
        case FESOM_VAR_3D_NODE_IFACE:
            gather_node(accum, nl, &s->gp, global, s->comm);
            break;
        case FESOM_VAR_2D_ELEM:
            gather_elem(accum, 1, &s->gp, global, s->comm);
            break;
        case FESOM_VAR_3D_ELEM_MID:
        case FESOM_VAR_3D_ELEM_IFACE:
            gather_elem(accum, nl, &s->gp, global, s->comm);
            break;
        case FESOM_VAR_2D_ELEM_VEC:
            gather_elem(accum, 2, &s->gp, global, s->comm);
            break;
    }
    (void)local_n; (void)local_e;

    if (s->mype == 0) {
        if (!window_matches(s, cal) || s->ncid[v] < 0) {
            /* Old window's file (if any for this var) is closed by
             * the caller (window-mismatch path). Open the new file. */
            open_var_file(s, v, cal);
        }
        size_t t_idx = (size_t)s->time_index[v];

        /* Write one record at time_index. Layout dispatched per kind. */
        if (kind == FESOM_VAR_2D_NODE || kind == FESOM_VAR_2D_ELEM) {
            size_t start[2] = {t_idx, 0};
            size_t count[2] = {1, (size_t)((kind == FESOM_VAR_2D_NODE) ? mesh->nod2D : mesh->elem2D)};
            NC_CHECK(nc_put_vara_double(s->ncid[v], s->var_id[v], start, count, global));
        } else if (kind == FESOM_VAR_2D_ELEM_VEC) {
            size_t start[3] = {t_idx, 0, 0};
            size_t count[3] = {1, 2, (size_t)mesh->elem2D};
            /* gather_elem laid out global as [elem][2]; transpose to [2][elem]. */
            real_t *xy = malloc((size_t)2 * (size_t)mesh->elem2D * sizeof(real_t));
            FESOM_CHECK(xy, "io_stream: oom (uv transpose)");
            for (int e = 0; e < mesh->elem2D; ++e) {
                xy[(size_t)0 * mesh->elem2D + e] = global[(size_t)e * 2 + 0];
                xy[(size_t)1 * mesh->elem2D + e] = global[(size_t)e * 2 + 1];
            }
            NC_CHECK(nc_put_vara_double(s->ncid[v], s->var_id[v], start, count, xy));
            free(xy);
        } else {
            /* 3D: layout in `global` is [N][nl] row-major.
             * netCDF dim order is (time, n_layers, N) — transpose to [n_layers][N]. */
            int N = kind_is_node(kind) ? mesh->nod2D : mesh->elem2D;
            size_t start[3] = {t_idx, 0, 0};
            size_t count[3] = {1, (size_t)n_layers, (size_t)N};
            real_t *t = malloc((size_t)n_layers * (size_t)N * sizeof(real_t));
            FESOM_CHECK(t, "io_stream: oom (3D transpose)");
            for (int n = 0; n < N; ++n) {
                for (int k = 0; k < n_layers; ++k) {
                    t[(size_t)k * N + n] = global[(size_t)n * nl + k];
                }
            }
            NC_CHECK(nc_put_vara_double(s->ncid[v], s->var_id[v], start, count, t));
            free(t);
        }

        /* Time scalar + bounds.
         *
         *   MEAN    → time = midpoint of the just-closed period,
         *             bnds = [period_start, period_end].
         *   INSTANT → time = exact instant being captured,
         *             bnds = [time, time]   (degenerate per CF practice).
         *
         * The INSTANT branch can't use period_window(curr, period) because
         * that would return the bounds of the NEXT period (curr lies at
         * its very start), which is misleading for an instantaneous
         * sample. */
        double t_start, t_mid, t_end;
        if (s->output_kind == FESOM_OUT_INSTANT) {
            t_mid   = fesom_calendar_seconds_since_origin(cal);
            t_start = t_end = t_mid;
        } else {
            fesom_calendar_period_window(cal, s->period, &t_start, &t_mid, &t_end);
        }
        size_t s1[1] = {t_idx}, c1[1] = {1};
        NC_CHECK(nc_put_vara_double(s->ncid[v], s->time_id[v], s1, c1, &t_mid));
        size_t s2[2] = {t_idx, 0}, c2[2] = {1, 2};
        double bnds[2] = {t_start, t_end};
        NC_CHECK(nc_put_vara_double(s->ncid[v], s->bnds_id[v], s2, c2, bnds));

        /* Force-flush so a crash mid-window doesn't lose written records. */
        NC_CHECK(nc_sync(s->ncid[v]));

        s->time_index[v] += 1;
        free(global);
    }
}

static void flush_all(fesom_io_stream_t *s,
                      const fesom_calendar_t *flush_anchor,
                      const struct fesom_mesh *mesh)
{
    if (s->n_accum == 0 && s->output_kind == FESOM_OUT_MEAN) return;

    /* Window rollover: if the file open for any variable is for a
     * different rollover window than `flush_anchor`, close all files
     * before flushing into the new window. */
    if (!window_matches(s, flush_anchor)) {
        close_open_files(s);
    }

    for (int v = 0; v < s->nvars; ++v) {
        flush_one_var(s, v, flush_anchor, mesh);
    }
    if (s->mype == 0) {
        window_set(s, flush_anchor);
    }
    /* Reset accumulators for the next window. */
    if (s->output_kind == FESOM_OUT_MEAN) {
        for (int v = 0; v < s->nvars; ++v) {
            memset(s->accum[v], 0, s->accum_sz[v] * sizeof(real_t));
        }
    }
    s->n_accum = 0;
}

/* ------------------------------------------------------------------ */
/* Public lifecycle                                                   */
/* ------------------------------------------------------------------ */

void fesom_io_stream_init(fesom_io_stream_t       *s,
                          fesom_period_kind_t      period,
                          fesom_output_kind_t      kind,
                          const fesom_var_desc_t  *vars,
                          int                      nvars,
                          const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          const char              *out_dir)
{
    memset(s, 0, sizeof(*s));
    s->period      = period;
    s->output_kind = kind;
    s->vars        = vars;
    s->nvars       = nvars;
    s->mesh        = mesh;
    s->comm        = partit->MPI_COMM_FESOM;
    s->mype        = partit->mype;
    s->open_year   = -1;
    s->open_month  = 0;
    s->open_day    = 0;
    s->n_accum     = 0;

    /* strdup is POSIX, not C99; do it by hand to keep -std=c99 clean. */
    {
        size_t n = strlen(out_dir) + 1;
        s->out_dir = malloc(n);
        FESOM_CHECK(s->out_dir, "io_stream: oom (out_dir)");
        memcpy(s->out_dir, out_dir, n);
    }

    /* Cached gather plan — built ONCE here; never per-flush. */
    gather_plan_init(&s->gp, mesh, partit);

    /* Cached coords on rank 0 — built ONCE here; reused on every file open. */
    build_coord_cache(s);

    s->ncid       = malloc((size_t)nvars * sizeof(int));
    s->var_id     = malloc((size_t)nvars * sizeof(int));
    s->time_id    = malloc((size_t)nvars * sizeof(int));
    s->bnds_id    = malloc((size_t)nvars * sizeof(int));
    s->time_index = malloc((size_t)nvars * sizeof(int));
    s->dim_time_id= malloc((size_t)nvars * sizeof(int));
    s->accum_sz   = malloc((size_t)nvars * sizeof(size_t));
    FESOM_CHECK(s->ncid && s->var_id && s->time_id && s->bnds_id
                && s->time_index && s->dim_time_id && s->accum_sz,
                "io_stream: oom (per-var arrays)");

    /* Always allocate per-variable accum buffers. MEAN streams accumulate
     * into them across the period; INSTANT streams reuse them as scratch
     * (zero-then-resolver-fill on each flush trigger). */
    s->accum = calloc((size_t)nvars, sizeof(real_t *));
    FESOM_CHECK(s->accum, "io_stream: oom (accum array)");

    for (int v = 0; v < nvars; ++v) {
        s->ncid[v]        = -1;
        s->var_id[v]      = -1;
        s->time_id[v]     = -1;
        s->bnds_id[v]     = -1;
        s->time_index[v]  = 0;
        s->dim_time_id[v] = -1;
        s->accum_sz[v]    = fesom_io_stream_accum_size(vars[v].kind,
                                                       mesh->myDim_nod2D,
                                                       mesh->myDim_elem2D,
                                                       mesh->nl);
        s->accum[v] = calloc(s->accum_sz[v], sizeof(real_t));
        FESOM_CHECK(s->accum[v], "io_stream: oom (accum[%d] %zu)",
                    v, s->accum_sz[v]);
    }
}

void fesom_io_stream_step(fesom_io_stream_t      *s,
                          const fesom_state      *state,
                          const fesom_calendar_t *prev,
                          const fesom_calendar_t *curr,
                          struct fesom_partit    *partit)
{
    (void)partit;  /* comm is cached in s->gp; partit unused on hot path */
    const struct fesom_mesh *mesh = s->mesh;

    /* MEAN hot path: per-variable resolver adds local slice to accumulator.
     * No MPI, no netCDF, no malloc. INSTANT streams skip per-step work
     * and only act on the boundary trigger. */
    if (s->output_kind == FESOM_OUT_MEAN) {
        for (int v = 0; v < s->nvars; ++v) {
            s->vars[v].resolve(state, s->accum[v], s->accum_sz[v]);
        }
        s->n_accum += 1;
    }

    if (fesom_calendar_crossed(prev, curr, s->period)) {
        if (s->output_kind == FESOM_OUT_INSTANT) {
            /* Capture current state into accum (zero, then resolver fills
             * via its `+=` contract) and flush. The flush anchor is `curr`
             * because INSTANT records carry the instant time, not a
             * period midpoint. */
            for (int v = 0; v < s->nvars; ++v) {
                memset(s->accum[v], 0, s->accum_sz[v] * sizeof(real_t));
                s->vars[v].resolve(state, s->accum[v], s->accum_sz[v]);
            }
            s->n_accum = 1;          /* division by 1 in flush is a no-op */
            flush_all(s, curr, mesh);
        } else {
            /* MEAN: `prev` lies inside the period that just closed, so
             * fesom_calendar_period_window(prev, period) yields the
             * correct (start, mid, end) for the just-closed window. */
            flush_all(s, prev, mesh);
        }
    }
}

void fesom_io_stream_finalize(fesom_io_stream_t   *s,
                              const fesom_calendar_t *cal,
                              struct fesom_partit *partit)
{
    (void)partit;
    /* Flush whatever is in the buffer against the current calendar
     * (mid-window data — better to write than to lose). */
    if (s->n_accum > 0) {
        flush_all(s, cal, s->mesh);
    }
    close_open_files(s);
}

void fesom_io_stream_close(fesom_io_stream_t *s)
{
    close_open_files(s);
    if (s->accum) {
        for (int v = 0; v < s->nvars; ++v) free(s->accum[v]);
        free(s->accum);
    }
    free(s->accum_sz);
    free(s->ncid); free(s->var_id);
    free(s->time_id); free(s->bnds_id);
    free(s->time_index); free(s->dim_time_id);
    free(s->out_dir);
    free(s->cached_lon_node); free(s->cached_lat_node);
    free(s->cached_lon_elem); free(s->cached_lat_elem);
    gather_plan_free(&s->gp);
    memset(s, 0, sizeof(*s));
}
