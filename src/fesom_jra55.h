#ifndef FESOM_JRA55_H
#define FESOM_JRA55_H

#include "fesom_types.h"

#include "fesom_calendar.h"

struct fesom_partit;
struct fesom_mesh;

/*
 * JRA55-do daily forcing reader — literal port of:
 *   nc_readTimeGrid     (gen_surface_forcing.F90:223-538)
 *   nc_sbc_ini_fillnames(gen_surface_forcing.F90:540-588)
 *   nc_sbc_ini          (gen_surface_forcing.F90:603-700)
 *   getcoeffld          (gen_surface_forcing.F90:702-998)
 *   data_timeinterp     (gen_surface_forcing.F90:1000-1023)
 *   sbc_ini             (gen_surface_forcing.F90:1025-1300)
 *   sbc_do              (gen_surface_forcing.F90:1448-1825)
 *   julday              (gen_surface_forcing.F90:1828-1867)
 *   calendar_date       (gen_surface_forcing.F90:1870-1913)
 *   binarysearch        (gen_surface_forcing.F90:1940-1984)
 *
 * Plus the atmdata→physics-arrays distribution from gen_forcing_couple.F90:681-696.
 *
 * Field indices follow Fortran (i_xwind=0, etc.) and the namelist enables all of:
 *     uas, vas, huss, rsds, rlds, tas, prra, prsn  (8 fields, l_mslp=false)
 *
 * Time axis: NetCDF time is "days since nm_nc_iyear-nm_nc_imm-nm_nc_idd" (1900-01-01
 * for JRA55). Internal storage shifts to "julian days since year 0001":
 *     flf->nc_time = nc_time/nm_nc_freq + julday(nm_nc_iyear, ...) [Fortran lines 473, 477]
 *     flf->nc_time -= julday(flf->year_orig, 1, 1) [line 499]
 *     flf->nc_time += julday(yearnew, 1, 1)        [line 503]
 *
 * Spatial: bilinear interp from regular (lon, lat) grid to mesh nodes.
 *          Cyclic boundary handled by extending nc_lon by 2 halo cells:
 *          nc_lon[0]      := nc_lon[Nlon-1] - 360   (Fortran ic_cyclic block)
 *          nc_lon[Nlon-1] := nc_lon[Nlon-2] + 360
 *
 * Lat is flipped if file stores -90→90 (Fortran flip_lat==1, lines 519-526).
 *
 * Lifecycle:
 *   fesom_jra55_init()       — allocates state, parses namelist defaults
 *   fesom_jra55_open_year()  — opens NetCDF for given year (called at sbc_ini and on year change)
 *   fesom_jra55_step()       — reads/refreshes if t_indx changed, time-interpolates atmdata,
 *                              distributes to physics arrays (u_wind, shum, ...)
 *   fesom_jra55_free()       — closes files, frees memory
 */

enum {
    FESOM_JRA_XWIND = 0,
    FESOM_JRA_YWIND,
    FESOM_JRA_HUMI,
    FESOM_JRA_QSR,
    FESOM_JRA_QLW,
    FESOM_JRA_TAIR,
    FESOM_JRA_PREC,
    FESOM_JRA_SNOW,
    FESOM_JRA_NFLD
};

/*
 * Per-field NetCDF state — one for each of the 8 fields. Mirrors the Fortran
 * type flfi_type (gen_surface_forcing.F90:189-220), but only the fields we use.
 *
 * sbcdata1[Nlon*Nlat] and sbcdata2[Nlon*Nlat] are the two time-bracketing snapshots.
 * The cache logic (Fortran lines 845-876) avoids re-reading data1 on the next
 * step when it equals the previous data2; we replicate by tracking which t_indx
 * each cache slot holds.
 */
typedef struct fesom_jra55_field {
    char       path_prefix[512]; /* e.g. ".../uas." → opens .../uas.YYYY.nc */
    char       var_name   [40];  /* e.g. "uas" */
    int        year_loaded;      /* year of currently open NetCDF (-1 = none) */
    int        ncid;             /* NetCDF file handle */
    int        Nlon;             /* nc_Nlon (with +2 halo for cyclic wrap) */
    int        Nlat;             /* nc_Nlat */
    int        Ntime;            /* nc_Ntime */
    real_t    *nc_lon;           /* [Nlon] (degrees east, halo at [0] and [Nlon-1]) */
    real_t    *nc_lat;           /* [Nlat] (degrees north, possibly flipped) */
    real_t    *nc_time;          /* [Ntime] (julian days since year 0001 of yearnew) */
    int        flip_lat;         /* 1 if file stores -90→90 (need to flip data rows) */
    char       calendar[24];     /* e.g. "gregorian" */
    int       *bilin_i;          /* [nod2D] longitude bracket index (1-based, like Fortran) */
    int       *bilin_j;          /* [nod2D] latitude bracket index */
    real_t    *coef_a;           /* [nod2D] linear-interp slope */
    real_t    *coef_b;           /* [nod2D] linear-interp intercept */
    real_t    *sbcdata1;         /* [Nlon*Nlat] snapshot at t_indx (with cyclic halo) */
    real_t    *sbcdata2;         /* [Nlon*Nlat] snapshot at t_indx_p1 */
    int        sbcdata1_t_index; /* -1 if not yet filled */
    int        sbcdata2_t_index;
    int        t_indx;           /* 1-based, like Fortran */
    int        t_indx_p1;
} fesom_jra55_field;

/*
 * Top-level JRA55 forcing state — owns the 8 fields plus the
 * time-interpolated arrays the physics layer reads.
 *
 * atmdata[FESOM_JRA_NFLD][nod2D] is the raw output of data_timeinterp.
 * Physics-layer aliases follow gen_forcing_couple.F90:681-696:
 *   u_wind     = atmdata[XWIND]
 *   v_wind     = atmdata[YWIND]
 *   shum       = atmdata[HUMI]
 *   shortwave  = atmdata[QSR]
 *   longwave   = atmdata[QLW]
 *   Tair (°C)  = atmdata[TAIR] - 273.15
 *   prec_rain  = atmdata[PREC] / 1000   (kg/m²/s → m/s)
 *   prec_snow  = atmdata[SNOW] / 1000
 */
typedef struct fesom_jra55 {
    fesom_jra55_field fld[FESOM_JRA_NFLD];

    /* Per-node time-interpolated quantities (in physics units). */
    real_t *u_wind;     /* [nod2D]  m/s  */
    real_t *v_wind;     /* [nod2D]  m/s  */
    real_t *shum;       /* [nod2D]  kg/kg */
    real_t *shortwave;  /* [nod2D]  W/m² (downward) */
    real_t *longwave;   /* [nod2D]  W/m² (downward) */
    real_t *Tair;       /* [nod2D]  °C  */
    real_t *prec_rain;  /* [nod2D]  m/s */
    real_t *prec_snow;  /* [nod2D]  m/s */

    /* Namelist scalars (defaults from work_core/namelist.forcing). */
    int     nm_nc_iyear;   /* 1900 */
    int     nm_nc_imm;     /* 1 */
    int     nm_nc_idd;     /* 1 */
    int     nm_nc_freq;    /* 1 (data points per day in the file's time units) */
    int     nm_nc_tmid;    /* 0 */
    int     include_fleapyear; /* 1 */

    /* Reference height for bulk formulae (default 10 m). */
    real_t  z_wind;
    real_t  z_tair;
    real_t  z_shum;
} fesom_jra55;

/* julday — port of gen_surface_forcing.F90:1828.
 * Returns julian day at noon for (yyyy,mm,dd) under the named calendar.
 * Currently supports {julian, gregorian, proleptic_gregorian, standard}; any
 * other string falls back to Fortran's else branch (365*yyyy). */
int fesom_jra_julday(int yyyy, int mm, int dd, const char *calendar);

/* binarysearch — port of gen_surface_forcing.F90:1940. Returns 1-based index
 * (matching Fortran). Caller treats `ind == 0` as "before array start" and
 * `ind == length` as "exact-match-at-end" or "past-end" (Fortran getcoeffld
 * branches on this distinction). */
int fesom_jra_binarysearch(int length, const real_t *array, real_t value);

void fesom_jra55_init(fesom_jra55 *jra, const struct fesom_mesh *mesh);
void fesom_jra55_free(fesom_jra55 *jra);

/* Open all 8 NetCDF files for the given year, read coordinate arrays + time
 * axis, build bilinear indices, prefetch first time slice & coef_a/coef_b.
 * Mirrors nc_sbc_ini (Fortran 603-700): once per year change. */
void fesom_jra55_open_year(fesom_jra55 *jra,
                           const struct fesom_mesh *mesh,
                           int yearnew);

/* Per-timestep update — mirrors sbc_do (Fortran 1448-1825):
 *   1. Compute rdate = julday(yearnew,1,1) + (daynew-1) + timenew/86400
 *   2. For each field: if rdate moved past current t_indx_p1, refresh
 *      sbcdata + coef_a/coef_b (cached read avoids re-fetch when possible).
 *   3. data_timeinterp: atmdata = rdate*coef_a + coef_b for each node.
 *   4. Distribute atmdata to physics-unit arrays (Tair-273.15, prec/1000, ...).
 *
 * yearnew/daynew/timenew follow g_clock semantics:
 *   daynew  = day-of-year, 1-based
 *   timenew = seconds since start of day */
void fesom_jra55_step(fesom_jra55 *jra,
                      const struct fesom_mesh *mesh,
                      struct fesom_partit     *partit,
                      int yearnew, int daynew, real_t timenew);

/* Calendar-driven adapter for fesom_jra55_step. Extracts year/day-of-year/
 * seconds-in-day from the model calendar and delegates to fesom_jra55_step.
 * Bit-identical to the legacy caller for runs <= 1 calendar year — see
 * Task 3 of the IO-system plan for the rdate-equivalence argument that
 * extends this to multi-year runs once the JRA55 reader gains year-rollover
 * support.
 *
 * NOTE: callers must still arrange for fesom_jra55_open_year() to be called
 * when the calendar crosses a year boundary; this adapter does not do it.
 * For single-year runs (Task 3 gate), no rollover happens. */
void fesom_jra55_step_cal(fesom_jra55 *jra,
                          const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          const fesom_calendar_t  *cal);

#endif /* FESOM_JRA55_H */
