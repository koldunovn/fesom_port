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

/* ------------------------------------------------------------------ */
/* Geographic-frame vector output (FESOM_IO_VECTOR_FRAME)             */
/*                                                                    */
/* The model holds horizontal vectors in the ROTATED mesh frame.      */
/* Fortran io_meandata rotates recognised (x,y) pairs to geographic   */
/* components before writing (do_rotation, io_meandata.F90:2972-3004) */
/* — analyses expect geographic east/north. Mirror that here: at      */
/* flush time, rotate the accumulation buffers of paired vars with    */
/* vector_r2g. Linear op => rotating the sums (or instant values)     */
/* before the mean division equals accumulating rotated fields.       */
/*                                                                    */
/* Knob: FESOM_IO_VECTOR_FRAME = "geo" (default; rotate) | "rotated"  */
/* (write raw model-frame components — the pre-2026-06-11 behaviour,  */
/* kept for debugging and bitwise comparison against old runs).       */
/*                                                                    */
/* Scope: STREAM outputs only. The debug snapshot writer (fesom_io.c) */
/* and the env-gated reference dumps stay NATIVE-frame on purpose —   */
/* they are model-state debug media (Fortran-parallel: io_meandata    */
/* rotates, restarts stay native). This also keeps the pinned         */
/* byte-gate (snap_*.nc) lineage valid across this change.            */
/* ------------------------------------------------------------------ */
static const struct { const char *x, *y; } io_vec_pairs[] = {
    /* the Fortran do_rotation pair list, io_meandata.F90:2972-2979 */
    { "u",        "v"        },
    { "uice",     "vice"     },
    { "unod",     "vnod"     },
    { "tau_x",    "tau_y"    },
    { "atmice_x", "atmice_y" },
    { "atmoce_x", "atmoce_y" },
    { "iceoce_x", "iceoce_y" },
    { "utemp",    "vtemp"    },
};

static int io_vec_frame_geo(void)
{
    static int loaded = 0, geo = 1;
    if (!loaded) {
        const char *f = getenv("FESOM_IO_VECTOR_FRAME");
        if (!f || !*f || strcmp(f, "geo") == 0)      geo = 1;
        else if (strcmp(f, "rotated") == 0)          geo = 0;
        else
            FESOM_DIE("io_stream: FESOM_IO_VECTOR_FRAME='%s' "
                      "(expected 'geo' or 'rotated')", f);
        loaded = 1;
    }
    return geo;
}

/* Self-contained rotated->geographic vector transform — a deliberate
 * DUPLICATE of the static rotation helpers in fesom_mesh.c (mirrors Fortran
 * gen_modules_rotate_grid.F90: r2g_matrix / r2g / vector_r2g, flag_coord=0).
 *
 * WHY duplicated instead of exporting from fesom_mesh.c: adding ANY code to
 * fesom_mesh.c perturbs that TU's -O3/-ffp-contract code generation (inline
 * heuristics for its static helpers shift) -> last-ULP changes in the mesh
 * geometry (observed: pgf_y 3e-18 at step 0) -> bitwise byte-gate lineage vs
 * the pinned baseline breaks for ALL fields via chaos. The IO subsystem is
 * the documented redesign exception; it keeps a private copy so the mesh TU
 * stays bit-identical. Do NOT "clean this up" by moving it back. */
static void io_vector_r2g(real_t *u, real_t *v, real_t rlon, real_t rlat)
{
#if FESOM_FORCE_ROTATION
    static real_t M[9];
    static int built = 0;
    if (!built) {
        const real_t al = FESOM_ALPHA_EULER_DEG * FESOM_RAD;
        const real_t be = FESOM_BETA_EULER_DEG  * FESOM_RAD;
        const real_t ga = FESOM_GAMMA_EULER_DEG * FESOM_RAD;
        /* row-major flat 3x3 — same numeric layout as Fortran r2g_matrix(i,j) */
        M[0] = cos(ga)*cos(al) - sin(ga)*cos(be)*sin(al);
        M[1] = cos(ga)*sin(al) + sin(ga)*cos(be)*cos(al);
        M[2] = sin(ga)*sin(be);
        M[3] = -sin(ga)*cos(al) - cos(ga)*cos(be)*sin(al);
        M[4] = -sin(ga)*sin(al) + cos(ga)*cos(be)*cos(al);
        M[5] = cos(ga)*sin(be);
        M[6] = sin(be)*sin(al);
        M[7] = -sin(be)*cos(al);
        M[8] = cos(be);
        built = 1;
    }

    /* rotated coords -> geographic coords (Fortran r2g: M^T on Cartesian) */
    real_t xr = cos(rlat) * cos(rlon);
    real_t yr = cos(rlat) * sin(rlon);
    real_t zr = sin(rlat);
    real_t xg = M[0]*xr + M[3]*yr + M[6]*zr;
    real_t yg = M[1]*xr + M[4]*yr + M[7]*zr;
    real_t zg = M[2]*xr + M[5]*yr + M[8]*zr;
    real_t glat = asin(zg);
    real_t glon = (yg == 0.0 && xg == 0.0) ? 0.0 : atan2(yg, xg);

    real_t tlon = *u, tlat = *v;              /* rotated east/north comps */
    /* rotated vector -> Cartesian rotated */
    real_t txg = -tlat*sin(rlat)*cos(rlon) - tlon*sin(rlon);
    real_t tyg = -tlat*sin(rlat)*sin(rlon) + tlon*cos(rlon);
    real_t tzg =  tlat*cos(rlat);
    /* Cartesian rotated -> Cartesian geographic (M^T) */
    real_t txr = M[0]*txg + M[3]*tyg + M[6]*tzg;
    real_t tyr = M[1]*txg + M[4]*tyg + M[7]*tzg;
    real_t tzr = M[2]*txg + M[5]*tyg + M[8]*tzg;
    /* Cartesian geographic -> geographic east/north components */
    *v = -sin(glat)*cos(glon)*txr - sin(glat)*sin(glon)*tyr + cos(glat)*tzr;
    *u = -sin(glon)*txr + cos(glon)*tyr;
#else
    (void)u; (void)v; (void)rlon; (void)rlat;
#endif
}

/* Rotate one (x,y) buffer pair in place, all owned columns, all levels.
 * Per-column position: node vars use the node's rotated coords; element
 * vars use the PLAIN MEAN of the 3 vertex rotated coords — a literal
 * mirror of io_meandata.F90:2993-2994 (NOT cyclic-aware across the
 * rotated dateline; same convention as Fortran, kept deliberately). */
static void rotate_pair_buffers(real_t *bx, real_t *by,
                                fesom_var_kind_t kind,
                                const struct fesom_mesh *mesh)
{
    const int nl = mesh->nl;
    switch (kind) {
    case FESOM_VAR_2D_NODE:
    case FESOM_VAR_3D_NODE_MID:
    case FESOM_VAR_3D_NODE_IFACE: {
        const int lev = (kind == FESOM_VAR_2D_NODE) ? 1 : nl;
        for (int n = 0; n < mesh->myDim_nod2D; ++n) {
            const real_t rlon = mesh->coord_nod2D[2*n + 0];
            const real_t rlat = mesh->coord_nod2D[2*n + 1];
            for (int k = 0; k < lev; ++k)
                io_vector_r2g(&bx[(size_t)n*lev + k],
                              &by[(size_t)n*lev + k], rlon, rlat);
        }
        break;
    }
    case FESOM_VAR_2D_ELEM:
    case FESOM_VAR_3D_ELEM_MID:
    case FESOM_VAR_3D_ELEM_IFACE: {
        const int lev = (kind == FESOM_VAR_2D_ELEM) ? 1 : nl;
        for (int e = 0; e < mesh->myDim_elem2D; ++e) {
            const int n0 = mesh->elem_nodes[3*e + 0];
            const int n1 = mesh->elem_nodes[3*e + 1];
            const int n2 = mesh->elem_nodes[3*e + 2];
            const real_t rlon = (mesh->coord_nod2D[2*n0]   + mesh->coord_nod2D[2*n1]
                               + mesh->coord_nod2D[2*n2])   / 3.0;
            const real_t rlat = (mesh->coord_nod2D[2*n0+1] + mesh->coord_nod2D[2*n1+1]
                               + mesh->coord_nod2D[2*n2+1]) / 3.0;
            for (int k = 0; k < lev; ++k)
                io_vector_r2g(&bx[(size_t)e*lev + k],
                              &by[(size_t)e*lev + k], rlon, rlat);
        }
        break;
    }
    default:
        FESOM_DIE("io_stream: vector rotation unsupported for var kind %d",
                  (int)kind);
    }
}

static void rotate_vector_pairs(fesom_io_stream_t *s,
                                const struct fesom_mesh *mesh)
{
    if (!io_vec_frame_geo()) return;
    static int note_done = 0;
    for (size_t p = 0; p < sizeof(io_vec_pairs)/sizeof(io_vec_pairs[0]); ++p) {
        int ix = -1, iy = -1;
        for (int v = 0; v < s->nvars; ++v) {
            if (strcmp(s->vars[v].name, io_vec_pairs[p].x) == 0) ix = v;
            if (strcmp(s->vars[v].name, io_vec_pairs[p].y) == 0) iy = v;
        }
        if (ix < 0 && iy < 0) continue;
        FESOM_CHECK(ix >= 0 && iy >= 0,
            "io_stream: vector component '%s' in a stream without its partner "
            "'%s' — add the partner or set FESOM_IO_VECTOR_FRAME=rotated",
            ix >= 0 ? io_vec_pairs[p].x : io_vec_pairs[p].y,
            ix >= 0 ? io_vec_pairs[p].y : io_vec_pairs[p].x);
        FESOM_CHECK(s->vars[ix].kind == s->vars[iy].kind,
            "io_stream: vector pair (%s,%s) has mismatched var kinds",
            io_vec_pairs[p].x, io_vec_pairs[p].y);
        if (s->mype == 0 && !note_done)
            fprintf(stderr, "[fesom_io] %s and %s rotated to geographic "
                    "components at output (FESOM_IO_VECTOR_FRAME=geo)\n",
                    io_vec_pairs[p].x, io_vec_pairs[p].y);
        rotate_pair_buffers(s->accum[ix], s->accum[iy], s->vars[ix].kind, mesh);
    }
    if (s->mype == 0) note_done = 1;
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

    /* Geo-frame vector output (default): rotate paired vector buffers
     * BEFORE the per-var flush — both MEAN sums and INSTANT captures are
     * fully written by this point, and buffers reset after the flush, so
     * the rotation applies exactly once per window. */
    rotate_vector_pairs(s, mesh);

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
    io_vec_frame_geo();   /* validate FESOM_IO_VECTOR_FRAME at init (fail fast) */
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
