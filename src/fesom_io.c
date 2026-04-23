/*
 * NetCDF snapshot writer. One file per snapshot — time-outer-most layout
 * is deferred until we have a long run worth grouping.
 */
#include "fesom_io.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_mesh.h"
#include "fesom_tracers.h"

#include <math.h>
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

/* Transpose a node-column array from (nod2D, nl) C row-major (with the
 * trailing slot at level nl-1 unused) to (n_layers, nod2D) for NetCDF. */
static void transpose_node_columns(const real_t *src, real_t *dst,
                                   int nod2D, int nl, int n_layers)
{
    for (int n = 0; n < nod2D; ++n) {
        for (int nz = 0; nz < n_layers; ++nz) {
            dst[(size_t)nz * (size_t)nod2D + (size_t)n] = src[(size_t)n * (size_t)nl + (size_t)nz];
        }
    }
}

/* Same for element-column arrays. */
static void transpose_elem_columns(const real_t *src, real_t *dst,
                                   int elem2D, int nl, int n_layers)
{
    for (int e = 0; e < elem2D; ++e) {
        for (int nz = 0; nz < n_layers; ++nz) {
            dst[(size_t)nz * (size_t)elem2D + (size_t)e] = src[(size_t)e * (size_t)nl + (size_t)nz];
        }
    }
}

/* Extract the u (or v) component from element-vector storage (elem, nl, 2)
 * and transpose to (n_layers, elem). */
static void extract_uv_component(const real_t *uv, real_t *dst,
                                 int elem2D, int nl, int n_layers, int comp)
{
    for (int e = 0; e < elem2D; ++e) {
        for (int nz = 0; nz < n_layers; ++nz) {
            size_t k = (size_t)e * (size_t)nl * 2 + (size_t)nz * 2 + (size_t)comp;
            dst[(size_t)nz * (size_t)elem2D + (size_t)e] = uv[k];
        }
    }
}

void fesom_io_write_snapshot(const char                  *path,
                             int                          step_n,
                             real_t                       dt,
                             const struct fesom_mesh     *mesh,
                             const struct fesom_dyn      *dyn,
                             const struct fesom_tracers  *tracers,
                             const struct fesom_aux      *aux)
{
    const int nod2D  = mesh->nod2D;
    const int elem2D = mesh->elem2D;
    const int nl     = mesh->nl;
    const int n_lay  = nl - 1;

    int ncid;
    NC_CHECK(nc_create(path, NC_CLOBBER | NC_NETCDF4, &ncid));

    /* ---- Dimensions --------------------------------------------------- */
    int dim_nod, dim_elem, dim_nz, dim_nz_1, dim_three;
    NC_CHECK(nc_def_dim(ncid, "nod2",  (size_t)nod2D,  &dim_nod));
    NC_CHECK(nc_def_dim(ncid, "elem",  (size_t)elem2D, &dim_elem));
    NC_CHECK(nc_def_dim(ncid, "nz",    (size_t)nl,     &dim_nz));
    NC_CHECK(nc_def_dim(ncid, "nz_1",  (size_t)n_lay,  &dim_nz_1));
    NC_CHECK(nc_def_dim(ncid, "three", 3,              &dim_three));

    /* ---- Variable definitions --------------------------------------- */
    int var_lon, var_lat, var_zbar, var_Z;
    int var_elnodes, var_nlev_nod, var_nlev_elem;
    int var_T, var_S, var_eta, var_w, var_u, var_v;
    int var_dens, var_bv, var_pgfx, var_pgfy, var_Kv, var_Av;

    int dims_nod[1]    = { dim_nod };
    int dims_elem[1]   = { dim_elem };
    int dims_nz[1]     = { dim_nz };
    int dims_nz1[1]    = { dim_nz_1 };
    int dims_T[2]      = { dim_nz_1, dim_nod };
    int dims_w[2]      = { dim_nz,   dim_nod };
    int dims_uv[2]     = { dim_nz_1, dim_elem };
    int dims_eln[2]    = { dim_elem, dim_three };

    NC_CHECK(nc_def_var(ncid, "lon",            NC_DOUBLE, 1, dims_nod,  &var_lon));
    NC_CHECK(nc_def_var(ncid, "lat",            NC_DOUBLE, 1, dims_nod,  &var_lat));
    NC_CHECK(nc_def_var(ncid, "zbar",           NC_DOUBLE, 1, dims_nz,   &var_zbar));
    NC_CHECK(nc_def_var(ncid, "Z",              NC_DOUBLE, 1, dims_nz1,  &var_Z));
    NC_CHECK(nc_def_var(ncid, "elem_nodes",     NC_INT,    2, dims_eln,  &var_elnodes));
    NC_CHECK(nc_def_var(ncid, "nlevels_nod2D",  NC_INT,    1, dims_nod,  &var_nlev_nod));
    NC_CHECK(nc_def_var(ncid, "nlevels",        NC_INT,    1, dims_elem, &var_nlev_elem));

    NC_CHECK(nc_def_var(ncid, "T",              NC_DOUBLE, 2, dims_T,    &var_T));
    NC_CHECK(nc_def_var(ncid, "S",              NC_DOUBLE, 2, dims_T,    &var_S));
    NC_CHECK(nc_def_var(ncid, "eta_n",          NC_DOUBLE, 1, dims_nod,  &var_eta));
    NC_CHECK(nc_def_var(ncid, "w",              NC_DOUBLE, 2, dims_w,    &var_w));
    NC_CHECK(nc_def_var(ncid, "u",              NC_DOUBLE, 2, dims_uv,   &var_u));
    NC_CHECK(nc_def_var(ncid, "v",              NC_DOUBLE, 2, dims_uv,   &var_v));
    NC_CHECK(nc_def_var(ncid, "density_m_rho0", NC_DOUBLE, 2, dims_T,    &var_dens));
    NC_CHECK(nc_def_var(ncid, "bvfreq",         NC_DOUBLE, 2, dims_T,    &var_bv));
    NC_CHECK(nc_def_var(ncid, "pgf_x",          NC_DOUBLE, 2, dims_uv,   &var_pgfx));
    NC_CHECK(nc_def_var(ncid, "pgf_y",          NC_DOUBLE, 2, dims_uv,   &var_pgfy));
    NC_CHECK(nc_def_var(ncid, "Kv",             NC_DOUBLE, 2, dims_T,    &var_Kv));
    NC_CHECK(nc_def_var(ncid, "Av",             NC_DOUBLE, 2, dims_uv,   &var_Av));

    /* ---- Attributes -------------------------------------------------- */
    NC_CHECK(nc_put_att_text(ncid, var_lon, "units", strlen("degrees_east"),  "degrees_east"));
    NC_CHECK(nc_put_att_text(ncid, var_lat, "units", strlen("degrees_north"), "degrees_north"));
    NC_CHECK(nc_put_att_text(ncid, var_zbar, "units", 1, "m"));
    NC_CHECK(nc_put_att_text(ncid, var_Z,    "units", 1, "m"));
    NC_CHECK(nc_put_att_text(ncid, var_T, "units", strlen("degrees_C"), "degrees_C"));
    NC_CHECK(nc_put_att_text(ncid, var_S, "units", 3, "PSU"));
    NC_CHECK(nc_put_att_text(ncid, var_eta, "units", 1, "m"));
    NC_CHECK(nc_put_att_text(ncid, var_w, "units", 3, "m/s"));
    NC_CHECK(nc_put_att_text(ncid, var_u, "units", 3, "m/s"));
    NC_CHECK(nc_put_att_text(ncid, var_v, "units", 3, "m/s"));
    NC_CHECK(nc_put_att_text(ncid, var_dens, "units", strlen("kg/m^3"), "kg/m^3"));
    NC_CHECK(nc_put_att_text(ncid, var_bv,   "units", strlen("s^-2"),  "s^-2"));
    NC_CHECK(nc_put_att_text(ncid, var_pgfx, "units", strlen("m/s^2"), "m/s^2"));
    NC_CHECK(nc_put_att_text(ncid, var_pgfy, "units", strlen("m/s^2"), "m/s^2"));
    NC_CHECK(nc_put_att_text(ncid, var_Kv,   "units", strlen("m^2/s"), "m^2/s"));
    NC_CHECK(nc_put_att_text(ncid, var_Av,   "units", strlen("m^2/s"), "m^2/s"));

    int step_attr = step_n;
    double dt_attr = (double)dt;
    double model_time = step_n * (double)dt;
    NC_CHECK(nc_put_att_int   (ncid, NC_GLOBAL, "step",               NC_INT,    1, &step_attr));
    NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "dt",                 NC_DOUBLE, 1, &dt_attr));
    NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "model_time_seconds", NC_DOUBLE, 1, &model_time));

    NC_CHECK(nc_enddef(ncid));

    /* ---- Coordinates: convert geo radians → degrees ----------------- */
    double *lon_d = malloc((size_t)nod2D * sizeof(double));
    double *lat_d = malloc((size_t)nod2D * sizeof(double));
    FESOM_CHECK(lon_d && lat_d, "io: oom (lon/lat)");
    for (int n = 0; n < nod2D; ++n) {
        lon_d[n] = (double)mesh->geo_coord_nod2D[2*n + 0] / FESOM_RAD;
        lat_d[n] = (double)mesh->geo_coord_nod2D[2*n + 1] / FESOM_RAD;
    }
    NC_CHECK(nc_put_var_double(ncid, var_lon, lon_d));
    NC_CHECK(nc_put_var_double(ncid, var_lat, lat_d));
    free(lon_d); free(lat_d);

    NC_CHECK(nc_put_var_double(ncid, var_zbar, mesh->zbar));
    NC_CHECK(nc_put_var_double(ncid, var_Z,    mesh->Z));

    /* ---- Topology --------------------------------------------------- */
    NC_CHECK(nc_put_var_int(ncid, var_elnodes,    mesh->elem_nodes));
    NC_CHECK(nc_put_var_int(ncid, var_nlev_nod,   mesh->nlevels_nod2D));
    NC_CHECK(nc_put_var_int(ncid, var_nlev_elem,  mesh->nlevels));

    /* ---- State (transpose to NetCDF (nz, nod) / (nz, elem) layout) -- */
    real_t *buf_nl  = malloc((size_t)n_lay  * (size_t)nod2D  * sizeof(real_t));
    real_t *buf_nl_e= malloc((size_t)n_lay  * (size_t)elem2D * sizeof(real_t));
    real_t *buf_w   = malloc((size_t)nl     * (size_t)nod2D  * sizeof(real_t));
    FESOM_CHECK(buf_nl && buf_nl_e && buf_w, "io: oom (transpose buffers)");

    transpose_node_columns(tracers->data[FESOM_TRACER_T].values, buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_T, buf_nl));
    transpose_node_columns(tracers->data[FESOM_TRACER_S].values, buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_S, buf_nl));

    NC_CHECK(nc_put_var_double(ncid, var_eta, dyn->eta_n));

    transpose_node_columns(dyn->w, buf_w, nod2D, nl, nl);
    NC_CHECK(nc_put_var_double(ncid, var_w, buf_w));

    extract_uv_component(dyn->uv, buf_nl_e, elem2D, nl, n_lay, /*comp=*/0);
    NC_CHECK(nc_put_var_double(ncid, var_u, buf_nl_e));
    extract_uv_component(dyn->uv, buf_nl_e, elem2D, nl, n_lay, /*comp=*/1);
    NC_CHECK(nc_put_var_double(ncid, var_v, buf_nl_e));

    transpose_node_columns(aux->density_m_rho0, buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_dens, buf_nl));
    transpose_node_columns(aux->bvfreq, buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_bv, buf_nl));
    transpose_node_columns(aux->Kv, buf_nl, nod2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_Kv, buf_nl));

    transpose_elem_columns(aux->pgf_x, buf_nl_e, elem2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_pgfx, buf_nl_e));
    transpose_elem_columns(aux->pgf_y, buf_nl_e, elem2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_pgfy, buf_nl_e));
    transpose_elem_columns(aux->Av, buf_nl_e, elem2D, nl, n_lay);
    NC_CHECK(nc_put_var_double(ncid, var_Av, buf_nl_e));

    free(buf_nl); free(buf_nl_e); free(buf_w);

    NC_CHECK(nc_close(ncid));
}
