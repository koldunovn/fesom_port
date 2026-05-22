/*
 * JRA55-do daily forcing reader. Literal port of the routines listed in
 * fesom_jra55.h. No MPI yet — single rank reads all 8 files for the year.
 */
#include "fesom_jra55.h"
#include "fesom_constants.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"

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

/* ------------------------------------------------------------------------ *
 * julday — gen_surface_forcing.F90:1828                                    *
 * ------------------------------------------------------------------------ */
int fesom_jra_julday(int yyyy, int mm, int dd, const char *calendar)
{
    if (strcmp(calendar, "julian")              == 0 ||
        strcmp(calendar, "gregorian")           == 0 ||
        strcmp(calendar, "proleptic_gregorian") == 0 ||
        strcmp(calendar, "standard")            == 0) {

        const int IGREG = 15 + 31 * (10 + 12 * 1582);  /* Oct 15, 1582 */
        int jy = yyyy;
        if (jy == 0) {
            FESOM_DIE("julday: there is no year zero");
        }
        if (jy < 0) jy += 1;

        int jm;
        if (mm > 2) {
            jm = mm + 1;
        } else {
            jy -= 1;
            jm = mm + 13;
        }

        int julday = (int)(365.25 * jy) + (int)(30.6001 * jm) + dd + 1720995;
        if (dd + 31 * (mm + 12 * yyyy) >= IGREG) {
            int ja = (int)(0.01 * jy);
            julday = julday + 2 - ja + (int)(0.25 * ja);
        }
        return julday;
    }

    /* Fortran else branch: 365*yyyy */
    return 365 * yyyy;
}

/* ------------------------------------------------------------------------ *
 * binarysearch — gen_surface_forcing.F90:1940                              *
 *                                                                          *
 * Returns 1-based index of the largest array element ≤ value (within d=1e-9). *
 * Returns 0  if value < array[0]                                           *
 * Returns length if value > array[length-1]                                *
 * (Mirrors Fortran semantics: ind starts as `right` after loop exhausts;   *
 * when value > all, `right == length` so ind == length).                   *
 * ------------------------------------------------------------------------ */
int fesom_jra_binarysearch(int length, const real_t *array, real_t value)
{
    const real_t d = 1e-9;
    int left  = 1;        /* 1-based */
    int right = length;
    int ind   = right;
    while (1) {
        if (left > right) break;
        int middle = (int)floor((left + right) / 2.0 + 0.5);  /* Fortran nint() */
        real_t am = array[middle - 1];
        if (fabs(am - value) <= d) {
            return middle;
        }
        if (am > value) {
            right = middle - 1;
        } else {
            left = middle + 1;
        }
    }
    ind = right;
    return ind;  /* 0 if value < array[0]; length if value > array[length-1] */
}

/* ------------------------------------------------------------------------ *
 * Static helpers for the file/var lookups (Fortran tries a few names)      *
 * ------------------------------------------------------------------------ */
static int try_nf90_inq_dimid(int ncid, const char *const *names, int n)
{
    int id = -1;
    for (int i = 0; i < n; ++i) {
        if (nc_inq_dimid(ncid, names[i], &id) == NC_NOERR) return id;
    }
    FESOM_DIE("nc_inq_dimid: none of the candidate names matched (first='%s')",
              n > 0 ? names[0] : "(none)");
}

static int try_nf90_inq_varid(int ncid, const char *const *names, int n)
{
    int id = -1;
    for (int i = 0; i < n; ++i) {
        if (nc_inq_varid(ncid, names[i], &id) == NC_NOERR) return id;
    }
    FESOM_DIE("nc_inq_varid: none of the candidate names matched (first='%s')",
              n > 0 ? names[0] : "(none)");
}

/* ------------------------------------------------------------------------ *
 * nc_readTimeGrid — gen_surface_forcing.F90:223-538                        *
 *                                                                          *
 * Opens flf->path_prefix + yearnew + ".nc", reads lat/lon/time + calendar  *
 * attribute, allocates and fills nc_lon/nc_lat/nc_time. Closes the file    *
 * (getcoeffld will reopen with nc_get_vara on each refresh).               *
 * ------------------------------------------------------------------------ */
static void nc_read_time_grid(fesom_jra55 *jra,
                              fesom_jra55_field *flf,
                              int yearnew)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s%04d.nc", flf->path_prefix, yearnew);

    int ncid;
    if (nc_open(path, NC_NOWRITE, &ncid) != NC_NOERR) {
        FESOM_DIE("nc_open failed: %s", path);
    }

    /* Dimensions: try LAT/lat/latitude/LAT1, etc. */
    static const char *lat_names [] = {"LAT",  "lat",  "latitude",  "LAT1"};
    static const char *lon_names [] = {"LON",  "lon",  "longitude", "LON1"};
    static const char *time_names[] = {"TIME", "time",              "TIME1"};

    int id_latd = try_nf90_inq_dimid(ncid, lat_names,  4);
    int id_lond = try_nf90_inq_dimid(ncid, lon_names,  4);
    int id_timed= try_nf90_inq_dimid(ncid, time_names, 3);

    int id_lat  = try_nf90_inq_varid(ncid, lat_names,  4);
    int id_lon  = try_nf90_inq_varid(ncid, lon_names,  4);
    int id_time = try_nf90_inq_varid(ncid, time_names, 3);

    size_t Nlat, Nlon, Ntime;
    NC_CHECK(nc_inq_dimlen(ncid, id_latd,  &Nlat));
    NC_CHECK(nc_inq_dimlen(ncid, id_lond,  &Nlon));
    NC_CHECK(nc_inq_dimlen(ncid, id_timed, &Ntime));

    /* +2 for cyclic halo at lon, like Fortran flf%nc_Nlon=flf%nc_Nlon+2 (line 353) */
    flf->Nlon  = (int)Nlon + 2;
    flf->Nlat  = (int)Nlat;
    flf->Ntime = (int)Ntime;

    /* Allocate (or reallocate just the time axis on year change). */
    if (!flf->nc_lat) {
        flf->nc_lon = malloc((size_t)flf->Nlon * sizeof(real_t));
        flf->nc_lat = malloc((size_t)flf->Nlat * sizeof(real_t));
    }
    if (flf->nc_time) free(flf->nc_time);
    flf->nc_time = malloc((size_t)flf->Ntime * sizeof(real_t));
    FESOM_CHECK(flf->nc_lon && flf->nc_lat && flf->nc_time,
                "fesom_jra55: nc_lon/lat/time alloc failed");

    /* Read lat (full extent) */
    {
        size_t st = 0, ct = (size_t)flf->Nlat;
        double *tmp = malloc(ct * sizeof(double));
        FESOM_CHECK(tmp, "fesom_jra55: tmp lat alloc");
        NC_CHECK(nc_get_vara_double(ncid, id_lat, &st, &ct, tmp));
        for (int i = 0; i < flf->Nlat; ++i) flf->nc_lat[i] = (real_t)tmp[i];
        free(tmp);
    }

    /* Read lon into [1..Nlon-2] (Fortran indices 2..Nlon-1) — leave halo for later */
    {
        size_t st = 0, ct = (size_t)(flf->Nlon - 2);
        double *tmp = malloc(ct * sizeof(double));
        FESOM_CHECK(tmp, "fesom_jra55: tmp lon alloc");
        NC_CHECK(nc_get_vara_double(ncid, id_lon, &st, &ct, tmp));
        for (size_t i = 0; i < ct; ++i) flf->nc_lon[i + 1] = (real_t)tmp[i];
        free(tmp);
        /* Fortran lines 378-379: nc_lon[1] := nc_lon[Nlon-1]; nc_lon[Nlon] := nc_lon[2] */
        flf->nc_lon[0]            = flf->nc_lon[flf->Nlon - 2];
        flf->nc_lon[flf->Nlon - 1] = flf->nc_lon[1];
    }

    /* Read time */
    {
        size_t st = 0, ct = (size_t)flf->Ntime;
        double *tmp = malloc(ct * sizeof(double));
        FESOM_CHECK(tmp, "fesom_jra55: tmp time alloc");
        NC_CHECK(nc_get_vara_double(ncid, id_time, &st, &ct, tmp));
        for (int i = 0; i < flf->Ntime; ++i) flf->nc_time[i] = (real_t)tmp[i];
        free(tmp);
    }

    /* calendar attribute */
    size_t att_len = 0;
    if (nc_inq_attlen(ncid, id_time, "calendar", &att_len) == NC_NOERR
        && att_len < sizeof(flf->calendar)) {
        char buf[64] = {0};
        NC_CHECK(nc_get_att_text(ncid, id_time, "calendar", buf));
        buf[att_len] = '\0';
        /* lowercase */
        for (size_t i = 0; i < att_len; ++i) {
            if (buf[i] >= 'A' && buf[i] <= 'Z') buf[i] += 32;
        }
        strncpy(flf->calendar, buf, sizeof(flf->calendar) - 1);
        flf->calendar[sizeof(flf->calendar) - 1] = '\0';
    } else {
        strcpy(flf->calendar, "none");
        fprintf(stderr, "[jra55] could not read calendar attribute, defaulting to 'none'\n");
    }

    /* Time-axis transform — Fortran lines 473-503. */
    /* nc_time = nc_time / nm_nc_freq + julday(nm_nc_iyear, nm_nc_imm, nm_nc_idd) */
    int jd_origin = fesom_jra_julday(jra->nm_nc_iyear, jra->nm_nc_imm,
                                     jra->nm_nc_idd, flf->calendar);
    for (int i = 0; i < flf->Ntime; ++i) {
        flf->nc_time[i] = flf->nc_time[i] / (real_t)jra->nm_nc_freq + (real_t)jd_origin;
    }
    /* Determine year_orig from nc_time[0]. */
    int yyyy = -1;
    {
        int julian = (int)flf->nc_time[0];
        const char *cal = flf->calendar;
        if (strcmp(cal, "julian")              == 0 ||
            strcmp(cal, "gregorian")           == 0 ||
            strcmp(cal, "proleptic_gregorian") == 0 ||
            strcmp(cal, "standard")            == 0) {
            const int IGREG = 2299161;
            int ja;
            if (julian >= IGREG) {
                double x = ((julian - 1867216) - 0.25) / 36524.25;
                ja = julian + 1 + (int)x - (int)(0.25 * x);
            } else {
                ja = julian;
            }
            int jb = ja + 1524;
            int jc = (int)(6680 + ((jb - 2439870) - 122.1) / 365.25);
            int jd = (int)(365 * jc + (0.25 * jc));
            int je = (int)((jb - jd) / 30.6001);
            int mm_o = je - 1;
            if (mm_o > 12) mm_o -= 12;
            yyyy = jc - 4715;
            if (mm_o > 2) yyyy -= 1;
            if (yyyy <= 0) yyyy -= 1;
        } else {
            yyyy = (int)((flf->nc_time[0] + 1.0e-12) / 365.0);
        }
    }
    int jd_year_orig = fesom_jra_julday(yyyy, 1, 1, flf->calendar);
    int jd_yearnew   = fesom_jra_julday(yearnew, 1, 1, flf->calendar);
    for (int i = 0; i < flf->Ntime; ++i) {
        flf->nc_time[i] = flf->nc_time[i] - (real_t)jd_year_orig + (real_t)jd_yearnew;
    }
    /* nm_nc_tmid != 1 branch — Fortran 505-512: shift to mid-interval. */
    if (jra->nm_nc_tmid != 1 && flf->Ntime > 1) {
        for (int i = 0; i < flf->Ntime - 1; ++i) {
            flf->nc_time[i] = (flf->nc_time[i+1] + flf->nc_time[i]) * 0.5;
        }
        flf->nc_time[flf->Ntime - 1] += (flf->nc_time[flf->Ntime - 1]
                                         - flf->nc_time[flf->Ntime - 2]) * 0.5;
    }

    /* flip_lat — Fortran 519-526: if file stores -90→90, flip lat axis. */
    flf->flip_lat = 0;
    if (flf->Nlat > 1 && flf->nc_lat[0] > flf->nc_lat[flf->Nlat - 1]) {
        flf->flip_lat = 1;
        for (int i = 0; i < flf->Nlat / 2; ++i) {
            real_t t = flf->nc_lat[i];
            flf->nc_lat[i] = flf->nc_lat[flf->Nlat - 1 - i];
            flf->nc_lat[flf->Nlat - 1 - i] = t;
        }
    }

    /* ic_cyclic — Fortran 534-537: extend cyclic halo. */
    flf->nc_lon[0]            = flf->nc_lon[0]            - 360.0;
    flf->nc_lon[flf->Nlon - 1] = flf->nc_lon[flf->Nlon - 1] + 360.0;

    /* Keep the file open: getcoeffld reads time slices repeatedly. */
    flf->ncid        = ncid;
    flf->year_loaded = yearnew;
}

/* ------------------------------------------------------------------------ *
 * Build bilinear interpolation indices — first half of nc_sbc_ini          *
 * (gen_surface_forcing.F90:646-687). 1-based, with -1/0 sentinels for      *
 * out-of-domain (north/south of grid).                                     *
 * ------------------------------------------------------------------------ */
static void build_bilin_indices(fesom_jra55_field *flf,
                                const struct fesom_mesh *mesh)
{
    int N = mesh->myDim_nod2D;
    if (!flf->bilin_i) {
        flf->bilin_i = malloc((size_t)N * sizeof(int));
        flf->bilin_j = malloc((size_t)N * sizeof(int));
        FESOM_CHECK(flf->bilin_i && flf->bilin_j,
                    "fesom_jra55: bilin_i/j alloc");
    }
    for (int i = 0; i < N; ++i) {
        real_t x = mesh->geo_coord_nod2D[2*i + 0] / FESOM_RAD;
        if (x < 0.0) x += 360.0;
        real_t y = mesh->geo_coord_nod2D[2*i + 1] / FESOM_RAD;

        /* longitude */
        if (x < flf->nc_lon[flf->Nlon - 1] && x >= flf->nc_lon[0]) {
            flf->bilin_i[i] = fesom_jra_binarysearch(flf->Nlon, flf->nc_lon, x);
        } else if (x < flf->nc_lon[0]) {
            flf->bilin_i[i] = -1;
        } else {
            flf->bilin_i[i] = 0;
        }

        /* latitude */
        if (y < flf->nc_lat[flf->Nlat - 1] && y >= flf->nc_lat[0]) {
            flf->bilin_j[i] = fesom_jra_binarysearch(flf->Nlat, flf->nc_lat, y);
        } else if (y < flf->nc_lat[0]) {
            flf->bilin_j[i] = -1;
        } else {
            flf->bilin_j[i] = 0;
        }
    }
}

/* ------------------------------------------------------------------------ *
 * Read one (Nlon-2)×Nlat time slice into sbcdata (with cyclic halo fill).  *
 * Mirrors Fortran lines 884-924 — except the variable's lat axis may be    *
 * stored -90→90 and we have to flip rows on read.                          *
 * sbcdata layout: row-major [Nlat][Nlon]                                   *
 * ------------------------------------------------------------------------ */
static void read_one_time_slice(const fesom_jra55_field *flf,
                                int t_indx_1based,
                                real_t *sbcdata)
{
    int varid;
    NC_CHECK(nc_inq_varid(flf->ncid, flf->var_name, &varid));

    /* Read interior columns [1 .. Nlon-2] (Fortran 2 .. Nlon-1). */
    size_t start[3] = { (size_t)(t_indx_1based - 1), 0, 0 };
    size_t count[3] = { 1, (size_t)flf->Nlat, (size_t)(flf->Nlon - 2) };
    /* JRA55 stores float; allocate scratch and convert. */
    float *tmp = malloc(count[1] * count[2] * sizeof(float));
    FESOM_CHECK(tmp, "fesom_jra55: time-slice scratch alloc");
    NC_CHECK(nc_get_vara_float(flf->ncid, varid, start, count, tmp));

    /* Place into sbcdata with optional lat-flip. */
    int Nlon = flf->Nlon, Nlat = flf->Nlat;
    if (flf->flip_lat) {
        for (int j = 0; j < Nlat; ++j) {
            int j_src = Nlat - 1 - j;
            for (int i = 0; i < Nlon - 2; ++i) {
                sbcdata[(size_t)j * Nlon + (i + 1)] =
                    (real_t)tmp[(size_t)j_src * (Nlon - 2) + i];
            }
        }
    } else {
        for (int j = 0; j < Nlat; ++j) {
            for (int i = 0; i < Nlon - 2; ++i) {
                sbcdata[(size_t)j * Nlon + (i + 1)] =
                    (real_t)tmp[(size_t)j * (Nlon - 2) + i];
            }
        }
    }
    free(tmp);

    /* Cyclic halo — Fortran 900-902 / 921-923:
     *   sbcdata[*][1] := sbcdata[*][Nlon-1]
     *   sbcdata[*][Nlon] := sbcdata[*][2]   (1-based; 0/Nlon-1 zero-based) */
    for (int j = 0; j < Nlat; ++j) {
        sbcdata[(size_t)j * Nlon + 0]        = sbcdata[(size_t)j * Nlon + (Nlon - 2)];
        sbcdata[(size_t)j * Nlon + (Nlon-1)] = sbcdata[(size_t)j * Nlon + 1];
    }
}

/* ------------------------------------------------------------------------ *
 * getcoeffld — gen_surface_forcing.F90:702-998                             *
 *                                                                          *
 * Locate t_indx, refresh sbcdata1/sbcdata2 (with cache), bilinear interp   *
 * to mesh nodes, set coef_a/coef_b for time interpolation.                 *
 * ------------------------------------------------------------------------ */
static void getcoeffld(fesom_jra55_field *flf,
                       const struct fesom_mesh *mesh,
                       real_t rdate)
{
    int Nlon = flf->Nlon;
    int Nlat = flf->Nlat;
    int Ntime = flf->Ntime;

    /* Allocate sbcdata1/2 lazily. */
    if (!flf->sbcdata1) {
        flf->sbcdata1 = malloc((size_t)Nlon * Nlat * sizeof(real_t));
        flf->sbcdata2 = malloc((size_t)Nlon * Nlat * sizeof(real_t));
        FESOM_CHECK(flf->sbcdata1 && flf->sbcdata2,
                    "fesom_jra55: sbcdata alloc");
        flf->sbcdata1_t_index = -1;
        flf->sbcdata2_t_index = -1;
    }

    /* Find t_indx via binary search. */
    int t_indx    = fesom_jra_binarysearch(Ntime, flf->nc_time, rdate);
    int t_indx_p1;
    real_t delta_t;
    if (t_indx < Ntime && t_indx > 0) {
        t_indx_p1 = t_indx + 1;
        delta_t   = flf->nc_time[t_indx_p1 - 1] - flf->nc_time[t_indx - 1];
        /* (We skip the leap-day jump-over branch — include_fleapyear=true.) */
    } else if (t_indx > 0) {
        /* No extrapolation into future. */
        t_indx    = Ntime;
        t_indx_p1 = t_indx;
        delta_t   = 1.0;
    } else {
        /* No extrapolation back in time. */
        t_indx    = 1;
        t_indx_p1 = 1;
        delta_t   = 1.0;
    }
    flf->t_indx    = t_indx;
    flf->t_indx_p1 = t_indx_p1;

    /* Determine cache reuse — Fortran lines 845-876.
     * The simple/correct strategy: if sbcdata2 currently holds t_indx, swap
     * it into sbcdata1 (no read), then refresh sbcdata2. */
    real_t *data1 = flf->sbcdata1;
    real_t *data2 = flf->sbcdata2;
    int     read_data1 = 1;
    int     read_data2 = 1;

    if (flf->sbcdata1_t_index == t_indx) {
        read_data1 = 0;
    } else if (flf->sbcdata2_t_index == t_indx) {
        /* Swap pointers: sbcdata2 ↔ sbcdata1 */
        real_t *tmp = flf->sbcdata1; flf->sbcdata1 = flf->sbcdata2; flf->sbcdata2 = tmp;
        int     ti  = flf->sbcdata1_t_index;
        flf->sbcdata1_t_index = flf->sbcdata2_t_index;
        flf->sbcdata2_t_index = ti;
        data1 = flf->sbcdata1; data2 = flf->sbcdata2;
        read_data1 = 0;
    }
    if (read_data1) {
        read_one_time_slice(flf, t_indx, data1);
        flf->sbcdata1_t_index = t_indx;
    }
    if (flf->sbcdata2_t_index == t_indx_p1 && t_indx_p1 != t_indx) {
        read_data2 = 0;
    }
    if (read_data2) {
        read_one_time_slice(flf, t_indx_p1, data2);
        flf->sbcdata2_t_index = t_indx_p1;
    }

    /* Bilinear horiz interp + build coef_a/coef_b — Fortran 939-997. */
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        int i = flf->bilin_i[n];
        int j = flf->bilin_j[n];
        int ip1 = i + 1;
        int jp1 = j + 1;
        real_t x = mesh->geo_coord_nod2D[2*n + 0] / FESOM_RAD;
        if (x < 0.0) x += 360.0;
        real_t y = mesh->geo_coord_nod2D[2*n + 1] / FESOM_RAD;

        int extrp = 0;
        if (i == 0)  { i = Nlon; ip1 = i; extrp += 1; }
        if (i == -1) { i = 1;    ip1 = i; extrp += 1; }
        if (j == 0)  { j = Nlat; jp1 = j; extrp += 2; }
        if (j == -1) { j = 1;    jp1 = j; extrp += 2; }
        /* 1-based → 0-based */
        int i0 = i - 1, ip10 = ip1 - 1;
        int j0 = j - 1, jp10 = jp1 - 1;

        real_t x1 = flf->nc_lon[i0],   x2 = flf->nc_lon[ip10];
        real_t y1 = flf->nc_lat[j0],   y2 = flf->nc_lat[jp10];

        real_t d1, d2;
        if (extrp == 0) {
            real_t denom = (x2 - x1) * (y2 - y1);
            real_t s_ij   = data1[(size_t)j0   * Nlon + i0];
            real_t s_p1j  = data1[(size_t)j0   * Nlon + ip10];
            real_t s_ijp1 = data1[(size_t)jp10 * Nlon + i0];
            real_t s_p1p1 = data1[(size_t)jp10 * Nlon + ip10];
            d1 = ( s_ij   * (x2-x)*(y2-y)
                 + s_p1j  * (x-x1)*(y2-y)
                 + s_ijp1 * (x2-x)*(y-y1)
                 + s_p1p1 * (x-x1)*(y-y1) ) / denom;
            real_t t_ij   = data2[(size_t)j0   * Nlon + i0];
            real_t t_p1j  = data2[(size_t)j0   * Nlon + ip10];
            real_t t_ijp1 = data2[(size_t)jp10 * Nlon + i0];
            real_t t_p1p1 = data2[(size_t)jp10 * Nlon + ip10];
            d2 = ( t_ij   * (x2-x)*(y2-y)
                 + t_p1j  * (x-x1)*(y2-y)
                 + t_ijp1 * (x2-x)*(y-y1)
                 + t_p1p1 * (x-x1)*(y-y1) ) / denom;
        } else if (extrp == 1) {
            real_t denom = (y2 - y1);
            d1 = ( data1[(size_t)j0   * Nlon + i0] * (y2-y)
                 + data1[(size_t)jp10 * Nlon + ip10] * (y-y1) ) / denom;
            d2 = ( data2[(size_t)j0   * Nlon + i0] * (y2-y)
                 + data2[(size_t)jp10 * Nlon + ip10] * (y-y1) ) / denom;
        } else if (extrp == 2) {
            real_t denom = (x2 - x1);
            d1 = ( data1[(size_t)j0   * Nlon + i0]   * (x2-x)
                 + data1[(size_t)jp10 * Nlon + ip10] * (x-x1) ) / denom;
            d2 = ( data2[(size_t)j0   * Nlon + i0]   * (x2-x)
                 + data2[(size_t)jp10 * Nlon + ip10] * (x-x1) ) / denom;
        } else {
            d1 = data1[(size_t)j0 * Nlon + i0];
            d2 = data2[(size_t)j0 * Nlon + i0];
        }

        flf->coef_a[n] = (d2 - d1) / delta_t;
        flf->coef_b[n] = d1 - flf->coef_a[n] * flf->nc_time[t_indx - 1];
    }
}

/* ------------------------------------------------------------------------ *
 * Public API                                                                *
 * ------------------------------------------------------------------------ */

void fesom_jra55_init(fesom_jra55 *jra, const struct fesom_mesh *mesh)
{
    memset(jra, 0, sizeof(*jra));

    /* Namelist defaults from work_core/namelist.forcing */
    static const char *prefixes[FESOM_JRA_NFLD] = {
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/uas.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/vas.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/huss.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/rsds.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/rlds.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/tas.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/prra.",
        "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/prsn.",
    };
    static const char *vars[FESOM_JRA_NFLD] = {
        "uas", "vas", "huss", "rsds", "rlds", "tas", "prra", "prsn"
    };

    int N_my = mesh->myDim_nod2D;
    int N    = mesh->myDim_nod2D + mesh->eDim_nod2D;   /* physics arrays cover halo */
    for (int f = 0; f < FESOM_JRA_NFLD; ++f) {
        fesom_jra55_field *flf = &jra->fld[f];
        strncpy(flf->path_prefix, prefixes[f], sizeof(flf->path_prefix) - 1);
        flf->path_prefix[sizeof(flf->path_prefix) - 1] = '\0';
        strncpy(flf->var_name, vars[f], sizeof(flf->var_name) - 1);
        flf->var_name[sizeof(flf->var_name) - 1] = '\0';
        flf->year_loaded = -1;
        flf->ncid        = -1;
        /* coef_a/b are computed only at myDim — the halo of jra->* is filled
         * by exchange after the step compute, so coefs need not extend. */
        flf->coef_a = calloc((size_t)N_my, sizeof(real_t));
        flf->coef_b = calloc((size_t)N_my, sizeof(real_t));
        flf->sbcdata1_t_index = -1;
        flf->sbcdata2_t_index = -1;
        FESOM_CHECK(flf->coef_a && flf->coef_b,
                    "fesom_jra55_init: coef_a/b alloc");
    }

    /* Physics arrays sized myDim+eDim: thermo's compute-over-halo loop in
     * fesom_ice_thermo.c reads jra->shortwave[halo_n] etc.; this used to
     * be an out-of-bounds read when arrays were sized myDim only. */
    jra->u_wind    = calloc((size_t)N, sizeof(real_t));
    jra->v_wind    = calloc((size_t)N, sizeof(real_t));
    jra->shum      = calloc((size_t)N, sizeof(real_t));
    jra->shortwave = calloc((size_t)N, sizeof(real_t));
    jra->longwave  = calloc((size_t)N, sizeof(real_t));
    jra->Tair      = calloc((size_t)N, sizeof(real_t));
    jra->prec_rain = calloc((size_t)N, sizeof(real_t));
    jra->prec_snow = calloc((size_t)N, sizeof(real_t));
    FESOM_CHECK(jra->u_wind && jra->v_wind && jra->shum && jra->shortwave
                && jra->longwave && jra->Tair && jra->prec_rain && jra->prec_snow,
                "fesom_jra55_init: physics arrays alloc");

    jra->nm_nc_iyear   = 1900;
    jra->nm_nc_imm     = 1;
    jra->nm_nc_idd     = 1;
    jra->nm_nc_freq    = 1;
    jra->nm_nc_tmid    = 0;
    jra->include_fleapyear = 1;
    jra->z_wind = 10.0;
    jra->z_tair = 10.0;
    jra->z_shum = 10.0;
}

void fesom_jra55_free(fesom_jra55 *jra)
{
    for (int f = 0; f < FESOM_JRA_NFLD; ++f) {
        fesom_jra55_field *flf = &jra->fld[f];
        if (flf->ncid >= 0) nc_close(flf->ncid);
        free(flf->nc_lon);
        free(flf->nc_lat);
        free(flf->nc_time);
        free(flf->bilin_i);
        free(flf->bilin_j);
        free(flf->coef_a);
        free(flf->coef_b);
        free(flf->sbcdata1);
        free(flf->sbcdata2);
    }
    free(jra->u_wind);    free(jra->v_wind);
    free(jra->shum);
    free(jra->shortwave); free(jra->longwave);
    free(jra->Tair);
    free(jra->prec_rain); free(jra->prec_snow);
    memset(jra, 0, sizeof(*jra));
}

void fesom_jra55_open_year(fesom_jra55 *jra,
                           const struct fesom_mesh *mesh,
                           int yearnew)
{
    for (int f = 0; f < FESOM_JRA_NFLD; ++f) {
        fesom_jra55_field *flf = &jra->fld[f];
        if (flf->ncid >= 0 && flf->year_loaded != yearnew) {
            nc_close(flf->ncid);
            flf->ncid = -1;
        }
        if (flf->ncid < 0) {
            nc_read_time_grid(jra, flf, yearnew);
        }
        /* Bilinear indices need rebuilding only on first init (mesh doesn't change),
         * but cheap to redo per year. */
        if (!flf->bilin_i) {
            build_bilin_indices(flf, mesh);
        }
    }
    /* After year change, prefetch first slice + coef_a/coef_b at start-of-year. */
    real_t rdate = (real_t)fesom_jra_julday(yearnew, 1, 1, jra->fld[0].calendar);
    for (int f = 0; f < FESOM_JRA_NFLD; ++f) {
        getcoeffld(&jra->fld[f], mesh, rdate);
    }
}

void fesom_jra55_step(fesom_jra55 *jra,
                      const struct fesom_mesh *mesh,
                      struct fesom_partit     *partit,
                      int yearnew, int daynew, real_t timenew)
{
    /* Fortran sbc_do (line 1500-ish): rdate = julday(yearnew,1,1) + (daynew-1) + timenew/86400 */
    real_t rdate = (real_t)fesom_jra_julday(yearnew, 1, 1, jra->fld[0].calendar)
                 + (real_t)(daynew - 1)
                 + timenew / 86400.0;

    /* Refresh data + coefs whenever rdate has crossed t_indx_p1 in any field. */
    for (int f = 0; f < FESOM_JRA_NFLD; ++f) {
        fesom_jra55_field *flf = &jra->fld[f];
        int need_refresh = 0;
        if (flf->t_indx <= 0) need_refresh = 1;
        else {
            real_t lo = flf->nc_time[flf->t_indx    - 1];
            real_t hi = flf->nc_time[flf->t_indx_p1 - 1];
            if (rdate < lo || rdate > hi) need_refresh = 1;
        }
        if (need_refresh) {
            getcoeffld(flf, mesh, rdate);
        }
    }

    /* data_timeinterp + distribution to physics arrays. */
    int N = mesh->myDim_nod2D;
    for (int n = 0; n < N; ++n) {
        real_t a_xw  = jra->fld[FESOM_JRA_XWIND].coef_a[n];
        real_t b_xw  = jra->fld[FESOM_JRA_XWIND].coef_b[n];
        real_t a_yw  = jra->fld[FESOM_JRA_YWIND].coef_a[n];
        real_t b_yw  = jra->fld[FESOM_JRA_YWIND].coef_b[n];
        real_t a_sh  = jra->fld[FESOM_JRA_HUMI ].coef_a[n];
        real_t b_sh  = jra->fld[FESOM_JRA_HUMI ].coef_b[n];
        real_t a_sw  = jra->fld[FESOM_JRA_QSR  ].coef_a[n];
        real_t b_sw  = jra->fld[FESOM_JRA_QSR  ].coef_b[n];
        real_t a_lw  = jra->fld[FESOM_JRA_QLW  ].coef_a[n];
        real_t b_lw  = jra->fld[FESOM_JRA_QLW  ].coef_b[n];
        real_t a_ta  = jra->fld[FESOM_JRA_TAIR ].coef_a[n];
        real_t b_ta  = jra->fld[FESOM_JRA_TAIR ].coef_b[n];
        real_t a_pr  = jra->fld[FESOM_JRA_PREC ].coef_a[n];
        real_t b_pr  = jra->fld[FESOM_JRA_PREC ].coef_b[n];
        real_t a_sn  = jra->fld[FESOM_JRA_SNOW ].coef_a[n];
        real_t b_sn  = jra->fld[FESOM_JRA_SNOW ].coef_b[n];

        jra->u_wind   [n] = rdate * a_xw + b_xw;
        jra->v_wind   [n] = rdate * a_yw + b_yw;
        /* Rotate JRA wind (geographic east/north) into the model's rotated frame
         * — mirror of Fortran gen_surface_forcing.F90:1550-1551 (vector_g2r on the
         * wind coefficients). Rotation is linear & time-independent, so rotating
         * the time-interpolated wind == rotating coef_a/coef_b. Without it the
         * Arctic wind stress (on BOTH ocean and ice) is mis-directed by the large
         * g2r angle near the geographic pole -> too-slow / misdirected drift.
         * Magnitude-preserving, so |wind|-based thermodynamics is unaffected. */
        fesom_vector_g2r(&jra->u_wind[n], &jra->v_wind[n],
                         mesh->geo_coord_nod2D[2*n + 0], mesh->geo_coord_nod2D[2*n + 1],
                         mesh->coord_nod2D[2*n + 0],     mesh->coord_nod2D[2*n + 1]);
        jra->shum     [n] = rdate * a_sh + b_sh;
        jra->shortwave[n] = rdate * a_sw + b_sw;
        jra->longwave [n] = rdate * a_lw + b_lw;
        jra->Tair     [n] = (rdate * a_ta + b_ta) - 273.15;
        jra->prec_rain[n] = (rdate * a_pr + b_pr) / 1000.0;
        jra->prec_snow[n] = (rdate * a_sn + b_sn) / 1000.0;
    }

    /* Halo exchange so thermo's compute-over-halo loop sees consistent
     * forcing at halo nodes. Without this, halo of jra->* is zero / stale,
     * which surfaces as huge-m_ice-at-halo bug in fesom_ice_thermo. */
    if (partit) {
        fesom_exchange_nod2D(jra->u_wind,    partit);
        fesom_exchange_nod2D(jra->v_wind,    partit);
        fesom_exchange_nod2D(jra->shum,      partit);
        fesom_exchange_nod2D(jra->shortwave, partit);
        fesom_exchange_nod2D(jra->longwave,  partit);
        fesom_exchange_nod2D(jra->Tair,      partit);
        fesom_exchange_nod2D(jra->prec_rain, partit);
        fesom_exchange_nod2D(jra->prec_snow, partit);
    }
}

/* ------------------------------------------------------------------ */
/* Task 3 of IO-system overhaul: calendar-driven adapter.
 *
 * Maps cal -> (yearnew, daynew, timenew) for the existing fesom_jra55_step.
 * The legacy caller in fesom_main.c computed daynew = floor(t_sec/86400)+1
 * from a fixed run-start year; the adapter just reads year/day-of-year/
 * seconds-in-day off the model calendar. The fesom_jra55_step internal
 * "rdate = julday(year,1,1) + (daynew-1) + timenew/86400" formula is
 * monotonic, so for a date crossing a year boundary the legacy
 * (yearnew=Y, daynew=366) and the calendar-form (yearnew=Y+1, daynew=1)
 * produce identical rdate values, hence identical interpolation. */
void fesom_jra55_step_cal(fesom_jra55 *jra,
                          const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          const fesom_calendar_t  *cal)
{
    int yearnew = cal->year;

    /* Multi-year rollover: when the calendar crosses Jan 1, swap in the
     * next year's NetCDF files. Invalidate the per-year sbcdata cache
     * indices first so the new year's getcoeffld can't false-match a
     * t_indx left over from the previous year (the indices live on
     * fesom_jra55_field across open_year calls). */
    if (jra->fld[0].year_loaded != yearnew) {
        for (int f = 0; f < FESOM_JRA_NFLD; ++f) {
            jra->fld[f].sbcdata1_t_index = -1;
            jra->fld[f].sbcdata2_t_index = -1;
        }
        fesom_jra55_open_year(jra, mesh, yearnew);
    }

    int daynew  = fesom_calendar_day_of_year(cal);
    real_t timenew = (real_t)((double)cal->hour * 3600.0
                            + (double)cal->minute * 60.0
                            + cal->second);
    fesom_jra55_step(jra, mesh, partit, yearnew, daynew, timenew);
}
