/*
 * SSS climatology + CORE2 runoff + per-timestep oce_fluxes salt+water balance.
 * Literal port of the Fortran routines listed in fesom_sss_runoff.h.
 */
#include "fesom_sss_runoff.h"

#include "fesom_constants.h"
#include "fesom_forcing.h"
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

/* ------------------------------------------------------------------------ *
 * interp_2d_field — gen_interpolation.F90:145-288                          *
 * (phase_flag=0 branch only — SSS / runoff are scalar fields).             *
 *                                                                          *
 * Bilinear interpolation from a regular (lon, lat) grid to a list of model *
 * nodes. Lon must be in [0, 360]; lat in [-90, 90] monotonically increasing.*
 * ------------------------------------------------------------------------ */
static void interp_2d_field(int            num_lon_reg,
                            int            num_lat_reg,
                            const real_t  *lon_reg,
                            const real_t  *lat_reg,
                            const real_t  *data_reg,    /* [num_lon × num_lat] row-major (lon-fastest) */
                            int            num_mod,
                            const real_t  *lon_mod,
                            const real_t  *lat_mod,
                            real_t        *data_mod)
{
    if (lon_reg[0] < 0.0 || lon_reg[num_lon_reg - 1] > 360.0) {
        FESOM_DIE("interp_2d_field: regular grid lon out of [0,360]");
    }

    for (int n = 0; n < num_mod; ++n) {
        real_t x = lon_mod[n];
        real_t y = lat_mod[n];
        int    ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l;
        real_t diff;
        real_t rt_lat1, rt_lat2, rt_lon1, rt_lon2;

        /* 1) north-south */
        if (y < lat_reg[0]) {
            ind_lat_h = 1;          /* Fortran 2 → C 1 */
            ind_lat_l = 0;          /* Fortran 1 → C 0 */
            y = lat_reg[0];
        } else if (y > lat_reg[num_lat_reg - 1]) {
            ind_lat_h = num_lat_reg - 1;
            ind_lat_l = num_lat_reg - 2;
            y = lat_reg[num_lat_reg - 1];
        } else {
            ind_lat_h = 1; ind_lat_l = 0;
            for (int i = 1; i < num_lat_reg; ++i) {
                if (lat_reg[i] >= y) {
                    ind_lat_h = i;
                    ind_lat_l = i - 1;
                    break;
                }
            }
        }
        diff = lat_reg[ind_lat_h] - lat_reg[ind_lat_l];
        rt_lat1 = (lat_reg[ind_lat_h] - y) / diff;
        rt_lat2 = 1.0 - rt_lat1;

        /* 2) east-west (with cyclic wrap when out of range) */
        if (x < lon_reg[0]) {
            ind_lon_h = 0;
            ind_lon_l = num_lon_reg - 1;
            diff      = lon_reg[ind_lon_h] + (360.0 - lon_reg[ind_lon_l]);
            rt_lon1   = (lon_reg[ind_lon_h] - x) / diff;
            rt_lon2   = 1.0 - rt_lon1;
        } else if (x > lon_reg[num_lon_reg - 1]) {
            ind_lon_h = 0;
            ind_lon_l = num_lon_reg - 1;
            diff      = lon_reg[ind_lon_h] + (360.0 - lon_reg[ind_lon_l]);
            rt_lon2   = (x - lon_reg[ind_lon_l]) / diff;
            rt_lon1   = 1.0 - rt_lon2;
        } else {
            ind_lon_h = 1; ind_lon_l = 0;
            for (int i = 1; i < num_lon_reg; ++i) {
                if (lon_reg[i] >= x) {
                    ind_lon_h = i;
                    ind_lon_l = i - 1;
                    break;
                }
            }
            diff = lon_reg[ind_lon_h] - lon_reg[ind_lon_l];
            rt_lon1 = (lon_reg[ind_lon_h] - x) / diff;
            rt_lon2 = 1.0 - rt_lon1;
        }

        real_t data_ll = data_reg[(size_t)ind_lat_l * num_lon_reg + ind_lon_l];
        real_t data_lh = data_reg[(size_t)ind_lat_h * num_lon_reg + ind_lon_l];
        real_t data_hl = data_reg[(size_t)ind_lat_l * num_lon_reg + ind_lon_h];
        real_t data_hh = data_reg[(size_t)ind_lat_h * num_lon_reg + ind_lon_h];

        data_mod[n] = (data_ll * rt_lon1 + data_hl * rt_lon2) * rt_lat1
                    + (data_lh * rt_lon1 + data_hh * rt_lon2) * rt_lat2;
    }
}

/* ------------------------------------------------------------------------ *
 * read_other_NetCDF — gen_modules_read_NetCDF.F90:6-184                    *
 *                                                                          *
 * Reads a 2D slice from a (lon, lat, time) NetCDF file, fills missing      *
 * values via 30-cell neighbourhood mean (when check_dummy), and bilinearly *
 * interpolates onto mesh nodes.                                            *
 * ------------------------------------------------------------------------ */
void fesom_read_other_NetCDF(const char *file,
                             const char *vari,
                             int         itime,
                             real_t     *model_2D,
                             int         check_dummy,
                             int         do_onvert,
                             const struct fesom_mesh *mesh)
{
    int ncid;
    if (nc_open(file, NC_NOWRITE, &ncid) != NC_NOERR) {
        FESOM_DIE("ERROR: CANNOT READ 2D netCDF FILE CORRECTLY: %s", file);
    }

    /* lat / lon dim+var ids */
    int latid, lonid, latvarid, lonvarid, varid;
    NC_CHECK(nc_inq_dimid(ncid, "lat", &latid));
    NC_CHECK(nc_inq_dimid(ncid, "lon", &lonid));
    size_t latlen, lonlen;
    NC_CHECK(nc_inq_dimlen(ncid, latid, &latlen));
    NC_CHECK(nc_inq_dimlen(ncid, lonid, &lonlen));

    NC_CHECK(nc_inq_varid(ncid, "lat", &latvarid));
    NC_CHECK(nc_inq_varid(ncid, "lon", &lonvarid));
    NC_CHECK(nc_inq_varid(ncid, vari,  &varid));

    real_t *lat = malloc(latlen * sizeof(real_t));
    real_t *lon = malloc(lonlen * sizeof(real_t));
    FESOM_CHECK(lat && lon, "read_other_NetCDF: lat/lon alloc");

    {
        size_t st = 0;
        double *tmp = malloc(latlen * sizeof(double));
        FESOM_CHECK(tmp, "read_other_NetCDF: tmp lat alloc");
        NC_CHECK(nc_get_vara_double(ncid, latvarid, &st, &latlen, tmp));
        for (size_t i = 0; i < latlen; ++i) lat[i] = (real_t)tmp[i];
        free(tmp);
        tmp = malloc(lonlen * sizeof(double));
        FESOM_CHECK(tmp, "read_other_NetCDF: tmp lon alloc");
        NC_CHECK(nc_get_vara_double(ncid, lonvarid, &st, &lonlen, tmp));
        for (size_t i = 0; i < lonlen; ++i) lon[i] = (real_t)tmp[i];
        free(tmp);
    }

    /* make sure range 0..360 (Fortran lines 88-92) */
    for (size_t i = 0; i < lonlen; ++i) {
        if (lon[i] < 0.0) lon[i] += 360.0;
    }

    real_t *ncdata      = malloc(lonlen * latlen * sizeof(real_t));
    real_t *ncdata_temp = malloc(lonlen * latlen * sizeof(real_t));
    FESOM_CHECK(ncdata && ncdata_temp, "read_other_NetCDF: ncdata alloc");
    memset(ncdata, 0, lonlen * latlen * sizeof(real_t));

    /* Read the time slice. NetCDF dim order is (time, lat, lon) for both files,
     * so the contiguous (lat, lon) block matches our ncdata[lat*lon] layout
     * directly. */
    {
        size_t start[3] = { (size_t)(itime - 1), 0, 0 };
        size_t count[3] = { 1, latlen, lonlen };
        float *tmp = malloc(latlen * lonlen * sizeof(float));
        FESOM_CHECK(tmp, "read_other_NetCDF: tmp data alloc");
        NC_CHECK(nc_get_vara_float(ncid, varid, start, count, tmp));
        /* Fortran stores ncdata(lon, lat) → row-major (lat, lon) for us. */
        for (size_t j = 0; j < latlen; ++j) {
            for (size_t i = 0; i < lonlen; ++i) {
                ncdata[j * lonlen + i] = (real_t)tmp[j * lonlen + i];
            }
        }
        free(tmp);
    }

    /* Fortran read_other_NetCDF (line 105) reads 'missing_value' but does NOT
       check the status, so an absent attribute is tolerated (the fill loop then
       matches only -99.0). The Sweeney chl field is pre-extrapolated and carries
       no missing_value attribute — so the read must be NON-FATAL. Sentinel
       -1e30 matches no real data when the attribute is absent; SSS/runoff files
       that DO carry missing_value overwrite it. */
    double miss_d = -1.0e30;
    nc_get_att_double(ncid, varid, "missing_value", &miss_d);
    real_t miss = (real_t)miss_d;

    NC_CHECK(nc_close(ncid));

    /* Fortran lines 113-139: missing-value fill.
     * Note Fortran also treats -99.0 as missing (line 116). */
    memcpy(ncdata_temp, ncdata, lonlen * latlen * sizeof(real_t));
    for (size_t i = 0; i < lonlen; ++i) {
        for (size_t j = 0; j < latlen; ++j) {
            real_t v = ncdata[j * lonlen + i];
            if (v == miss || v == (real_t)(-99.0)) {
                if (check_dummy) {
                    real_t aux = 0.0;
                    int    cnt = 0;
                    for (int k = 1; k <= 30; ++k) {
                        size_t i_lo = (i >= (size_t)k) ? i - k : 0;
                        size_t i_hi = (i + k < lonlen) ? i + k : lonlen - 1;
                        size_t j_lo = (j >= (size_t)k) ? j - k : 0;
                        size_t j_hi = (j + k < latlen) ? j + k : latlen - 1;
                        for (size_t ii = i_lo; ii <= i_hi; ++ii) {
                            for (size_t jj = j_lo; jj <= j_hi; ++jj) {
                                real_t w = ncdata_temp[jj * lonlen + ii];
                                if (w != miss && w != (real_t)(-99.0)) {
                                    aux += w;
                                    cnt += 1;
                                }
                            }
                        }
                        if (cnt > 0) {
                            ncdata[j * lonlen + i] = aux / (real_t)cnt;
                            break;
                        }
                    }
                } else {
                    ncdata[j * lonlen + i] = 0.0;
                }
            }
        }
    }

    /* Coordinate prep — Fortran 142-175 (do_onvert branch only). */
    int    num    = mesh->myDim_nod2D + mesh->eDim_nod2D;
    real_t *temp_x = malloc((size_t)num * sizeof(real_t));
    real_t *temp_y = malloc((size_t)num * sizeof(real_t));
    FESOM_CHECK(temp_x && temp_y, "read_other_NetCDF: temp_x/y alloc");
    if (do_onvert) {
        for (int n = 0; n < num; ++n) {
            real_t x = mesh->geo_coord_nod2D[2 * n + 0] / FESOM_RAD;
            real_t y = mesh->geo_coord_nod2D[2 * n + 1] / FESOM_RAD;
            if (x < 0.0) x += 360.0;
            temp_x[n] = x;
            temp_y[n] = y;
        }
    } else {
        FESOM_DIE("read_other_NetCDF: do_onvert=0 branch not ported");
    }

    interp_2d_field((int)lonlen, (int)latlen, lon, lat, ncdata,
                    num, temp_x, temp_y, model_2D);

    free(temp_x); free(temp_y);
    free(ncdata_temp); free(ncdata);
    free(lon); free(lat);
}

/* ------------------------------------------------------------------------ *
 * integrate_nod_2D — gen_support.F90:318-351                               *
 * Local sum over interior nodes only, then global Allreduce on multi-rank. *
 * ------------------------------------------------------------------------ */
static real_t integrate_nod_2D(const real_t *data,
                               const struct fesom_mesh *mesh,
                               struct fesom_partit     *partit)
{
    real_t lval = 0.0;
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        int ul = mesh->ulevels_nod2D[n] - 1;   /* 0-based */
        lval += data[n] * mesh->areasvol[FESOM_NODE3D(n, ul, mesh->nl)];
    }
    if (partit && partit->npes > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &lval, 1, MPI_DOUBLE, MPI_SUM,
                      partit->MPI_COMM_FESOM);
    }
    return lval;
}

/* ------------------------------------------------------------------------ *
 * Public                                                                    *
 * ------------------------------------------------------------------------ */
void fesom_sss_runoff_init(fesom_sss_runoff       *sr,
                           const struct fesom_mesh *mesh,
                           struct fesom_forcing    *forcing,
                           const char *sss_path,
                           const char *runoff_path)
{
    memset(sr, 0, sizeof(*sr));
    sr->sss_month_loaded = -1;

    /* Values from the CORE2 reference run's namelist.tra (NOT the oce_modules.F90
     * code defaults, which are ref_sss_local=.false./ref_sss=34.7). The reference
     * namelist sets ref_sss_local=.true. → the virtual-salt-flux reference salinity
     * rsss is the LOCAL surface salinity per node (ice_oce_coupling.F90:440), not a
     * global constant. Using the constant 34.7 over-strengthened the virtual salt
     * flux where local SSS is low (Arctic ~30) → systematic Arctic freshwater bias.
     * ref_sss (34.) is then unused, but set to the namelist value for fidelity. */
    sr->surf_relax_S  = 10.0 / (60.0 * 3600.0 * 24.0);   /* surf_relax_S = 1.929e-6 */
    sr->ref_sss       = 34.0;
    sr->ref_sss_local = 1;
    /* use_virt_salt set by which_ALE in Fortran oce_setup_step.F90:115-125;
     * for linfs (our case) this is .true. */
    sr->use_virt_salt = 1;

    if (sss_path && sss_path[0]) {
        strncpy(sr->sss_path, sss_path, sizeof(sr->sss_path) - 1);
    }
    if (runoff_path && runoff_path[0]) {
        strncpy(sr->runoff_path, runoff_path, sizeof(sr->runoff_path) - 1);

        /* Fortran sbc_ini lines 1267-1281 — CORE/CORE2 runoff is constant
         * (single time slice). Read once at init. */
        fesom_read_other_NetCDF(sr->runoff_path, "Foxx_o_roff", 1,
                                forcing->runoff,
                                /*check_dummy=*/0,
                                /*do_onvert=*/1,
                                mesh);
        /* Kg/s/m² → m/s (Fortran line 1281). read_other_NetCDF populated
         * runoff at myDim+eDim; apply the same unit conversion across halo so
         * any downstream read of runoff at halo sees m/s instead of the raw
         * kg/s/m² that's 1000× larger. */
        int N_full = mesh->myDim_nod2D + mesh->eDim_nod2D;
        for (int n = 0; n < N_full; ++n) {
            forcing->runoff[n] = forcing->runoff[n] / 1000.0;
        }
    }
}

void fesom_sss_runoff_free(fesom_sss_runoff *sr)
{
    memset(sr, 0, sizeof(*sr));
}

void fesom_sss_runoff_step(fesom_sss_runoff           *sr,
                           const struct fesom_mesh    *mesh,
                           const struct fesom_tracers *tracers,
                           struct fesom_forcing       *forcing,
                           struct fesom_partit        *partit,
                           int yearnew, int month_now,
                           int update_monthly_flag)
{
    (void)yearnew;

    /* sbc_do SSS section (Fortran 1571-1582) — read fresh slice each new month.
     * The Fortran fires update_monthly_flag on the LAST timestep of the month
     * (month still = ending month M) and then does `if (mstep>1) i=i+1` to load
     * the upcoming month M+1. OUR update_monthly_flag fires on the FIRST step of
     * the new month (fesom_calendar_crossed), where month_now is ALREADY M+1, so
     * we must NOT add the +1 — doing so double-advances and skips a month (the
     * old code loaded 1,3,4,..,12 and never read February → marginal seas absent
     * from PHC, e.g. Red Sea/Baltic, restored to the wrong month's extrapolated
     * target → systematic salty bias + missing seasonal cycle). Use month_now. */
    if (sr->sss_path[0] && update_monthly_flag) {
        int i = month_now;
        fesom_read_other_NetCDF(sr->sss_path, "SALT", i, forcing->Ssurf,
                                /*check_dummy=*/1, /*do_onvert=*/1, mesh);
        sr->sss_month_loaded = i;
    }

    /* ---------- oce_fluxes (ice_oce_coupling.F90:436-646) ---------- */

    /* All write-loops below run myDim+eDim to match Fortran ice_oce_coupling.F90
     * (lines 373-509) — halo entries of water_flux/virtual_salt/relax_salt must
     * stay in sync with the owner so any downstream reader at halo positions
     * sees correct values. Matches the port pattern documented in
     * `feedback_write_loops_halo.md`. */
    const int N_full = mesh->myDim_nod2D + mesh->eDim_nod2D;

    /* Phase C3b: removed runoff subtraction.
     * Runoff is now folded into ice->flx_fw inside fesom_therm_ice (via
     * `prec = rain + runo + snow*(1-A)`), then plumbed into water_flux by
     * fesom_ice_oce_fluxes. Subtracting again here would double-count.
     * See feedback_runoff_fold_when_ice.md (resolved 2026-04-25). */

    /* Fortran lines 436-457 — virtual_salt = rsss * water_flux, then
     * subtract global mean to enforce zero net salt flux. */
    if (sr->use_virt_salt) {
        real_t rsss = sr->ref_sss;
        for (int n = 0; n < N_full; ++n) {
            if (sr->ref_sss_local) {
                int ul = mesh->ulevels_nod2D[n] - 1;
                rsss = tracers->data[FESOM_TRACER_S].values[FESOM_NODE3D(n, ul, mesh->nl)];
            }
            forcing->virtual_salt[n] = rsss * forcing->water_flux[n];
        }
        real_t net = integrate_nod_2D(forcing->virtual_salt, mesh, partit) / mesh->ocean_area;
        for (int n = 0; n < N_full; ++n) {
            if (mesh->ulevels_nod2D[n] > 1) continue;     /* cavity */
            forcing->virtual_salt[n] -= net;
        }
    } else {
        for (int n = 0; n < N_full; ++n) forcing->virtual_salt[n] = 0.0;
    }

    /* Fortran lines 498-524 — relax_salt = surf_relax_S * (Ssurf - S_top),
     * then subtract global mean. Only meaningful where we actually loaded SSS. */
    if (sr->sss_path[0] && sr->sss_month_loaded > 0) {
        for (int n = 0; n < N_full; ++n) {
            int ul = mesh->ulevels_nod2D[n] - 1;
            real_t Stop = tracers->data[FESOM_TRACER_S].values[FESOM_NODE3D(n, ul, mesh->nl)];
            forcing->relax_salt[n] = sr->surf_relax_S * (forcing->Ssurf[n] - Stop);
        }
        real_t net = integrate_nod_2D(forcing->relax_salt, mesh, partit) / mesh->ocean_area;
        for (int n = 0; n < N_full; ++n) {
            if (mesh->ulevels_nod2D[n] > 1) continue;     /* cavity */
            forcing->relax_salt[n] -= net;
        }
    } else {
        for (int n = 0; n < N_full; ++n) forcing->relax_salt[n] = 0.0;
    }

    /* Fortran lines 544-563 — assemble local 'flux' = E + P + S + R, then
     * subtract its global mean from water_flux to enforce zero net mass.
     * (Without sea ice: a_ice_old = 0 → snow factor (1-a_ice_old)=1.) */
    {
        real_t *flux = malloc((size_t)N_full * sizeof(real_t));
        FESOM_CHECK(flux, "sss_runoff_step: flux alloc");
        for (int n = 0; n < N_full; ++n) {
            /* evaporation: in Fortran ice_thermo this is set to the bulk
             * evap (positive up). We store this inside water_flux via bulk
             * already, so the equivalent decomposition for our state is:
             *   flux ≡ -fresh_wa_flux + runoff = water_flux + runoff
             * (matches Fortran line 624 comment: "+ sign because we switched
             *  the sign of water_flux with water_flux=-fresh_wa_flux but evap,
             *  prec and runoff still have their original sign"). */
            flux[n] = forcing->water_flux[n] + forcing->runoff[n];
        }
        real_t net = integrate_nod_2D(flux, mesh, partit) / mesh->ocean_area;
        for (int n = 0; n < N_full; ++n) {
            forcing->water_flux[n] += net;
        }
        free(flux);
    }
}

/* ------------------------------------------------------------------ */
/* Task 4 of IO-system overhaul: calendar-driven adapter.
 *
 * Legacy caller built month_now from a static non-leap days_in_month[12]
 * table in fesom_main.c, and update_monthly_flag from a manual
 * month_prev != month_now comparison. Both come off the model calendar now:
 *   month_now           = cal->month
 *   update_monthly_flag = (n_step == 1) || crossed(prev, cal, MONTHLY)
 *
 * For runs starting on Jan 1 the n_step==1 trigger fires regardless of
 * `prev`, matching the legacy "fire on mstep==1" semantic exactly. */
void fesom_sss_runoff_step_cal(fesom_sss_runoff           *sr,
                               const struct fesom_mesh    *mesh,
                               const struct fesom_tracers *tracers,
                               struct fesom_forcing       *forcing,
                               struct fesom_partit        *partit,
                               int                         n_step,
                               const fesom_calendar_t     *prev,
                               const fesom_calendar_t     *cal)
{
    int yearnew  = cal->year;
    int month_now = cal->month;
    int update_monthly_flag = (n_step == 1)
                            || fesom_calendar_crossed(prev, cal, FESOM_PERIOD_MONTHLY);
    fesom_sss_runoff_step(sr, mesh, tracers, forcing, partit,
                          yearnew, month_now, update_monthly_flag);
}
