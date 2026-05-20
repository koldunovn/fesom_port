#ifndef FESOM_SSS_RUNOFF_H
#define FESOM_SSS_RUNOFF_H

#include "fesom_calendar.h"
#include "fesom_types.h"

struct fesom_mesh;
struct fesom_tracers;
struct fesom_forcing;
struct fesom_partit;

/*
 * Phase 3 step 25 — SSS restoring (PHC2_salx) + CORE2 runoff.
 *
 * Literal port of:
 *   read_other_NetCDF      (gen_modules_read_NetCDF.F90:6-184)
 *   interp_2d_field        (gen_interpolation.F90:145-288)
 *   sbc_ini runoff section (gen_surface_forcing.F90:1262-1331)
 *   sbc_do  SSS+runoff     (gen_surface_forcing.F90:1565-1652)
 *   oce_fluxes salt+water  (ice_oce_coupling.F90:436-646)
 *   integrate_nod_2D       (gen_support.F90:318-351)
 *
 * Structure is: a self-contained module that owns the SSS climatology raw
 * file handle (kept open so we can re-read each month) and the runoff time-1
 * snapshot (read once at init for CORE2). Per timestep it:
 *   1. If month boundary or first step → read SSS slice for current month
 *   2. Subtracts runoff from water_flux
 *   3. Computes virtual_salt = rsss * water_flux, balances net to zero
 *   4. Computes relax_salt = surf_relax_S * (Ssurf - S_top), balances net to zero
 *   5. Globally balances water_flux mean
 *
 * Defaults (Fortran oce_modules.F90:30-37):
 *   surf_relax_S = 10/(60·3600·24) = 1.929e-6 s⁻¹
 *   ref_sss      = 34.7 PSU
 *   use_virt_salt = true (linfs)
 */
typedef struct fesom_sss_runoff {
    char    sss_path  [512];      /* PHC2_salx.nc absolute path */
    char    runoff_path[512];     /* CORE2_runoff.nc absolute path */
    int     sss_month_loaded;     /* 1..12, -1 if not yet loaded */

    /* Tunables — Fortran namelist /tracer_phys/. */
    real_t  surf_relax_S;         /* 1/s, default 10/(60·3600·24) */
    real_t  ref_sss;              /* PSU, default 34.7 */
    int     use_virt_salt;        /* 1 for linfs (we are linfs), 0 otherwise */
    int     ref_sss_local;        /* 0 (use ref_sss); 1 → rsss = S_top per node */
} fesom_sss_runoff;

/* Initialise. Reads runoff once (CORE2 source). Does NOT read SSS yet — that
 * happens in fesom_sss_runoff_step on the first call.
 *   sss_path    : path to PHC2_salx.nc (or NULL to disable SSS restoring)
 *   runoff_path : path to CORE2_runoff.nc (or NULL to disable runoff) */
void fesom_sss_runoff_init(fesom_sss_runoff       *sr,
                           const struct fesom_mesh *mesh,
                           struct fesom_forcing    *forcing,
                           const char *sss_path,
                           const char *runoff_path);

void fesom_sss_runoff_free(fesom_sss_runoff *sr);

/* Per-timestep entry — to be called after fesom_bulk_compute, before
 * fesom_timestep. Replicates the relevant sections of sbc_do (monthly SSS
 * refresh) and oce_fluxes (virtual_salt + relax_salt + water_flux balance). */
void fesom_sss_runoff_step(fesom_sss_runoff           *sr,
                           const struct fesom_mesh    *mesh,
                           const struct fesom_tracers *tracers,
                           struct fesom_forcing       *forcing,
                           struct fesom_partit        *partit,
                           int yearnew, int month_now,
                           int update_monthly_flag);

/* Calendar-driven adapter for fesom_sss_runoff_step. Computes month_now
 * from `cal->month` and update_monthly_flag from `n_step == 1` OR a
 * MONTHLY-period crossing between `prev` and `cal`. Caller maintains
 * `prev` as the calendar state from the previous timestep — i.e. saves
 * `cal` before calling `fesom_calendar_advance`. */
void fesom_sss_runoff_step_cal(fesom_sss_runoff           *sr,
                               const struct fesom_mesh    *mesh,
                               const struct fesom_tracers *tracers,
                               struct fesom_forcing       *forcing,
                               struct fesom_partit        *partit,
                               int                         n_step,
                               const fesom_calendar_t     *prev,
                               const fesom_calendar_t     *cal);

/* Helper used by both SSS reader and the runoff reader.
 * Literal port of read_other_NetCDF (gen_modules_read_NetCDF.F90:6-184).
 *   file        : NetCDF path
 *   vari        : variable name (e.g. "SALT", "Foxx_o_roff")
 *   itime       : time slice (1-based, like Fortran)
 *   model_2D    : output [nod2D] interpolated to mesh
 *   check_dummy : 1 = fill missing values from neighbours; 0 = replace with 0
 *   do_onvert   : 1 = interpolate to mesh vertices; 0 = element centroids (we
 *                 only use 1 here) */
void fesom_read_other_NetCDF(const char *file,
                             const char *vari,
                             int         itime,
                             real_t     *model_2D,
                             int         check_dummy,
                             int         do_onvert,
                             const struct fesom_mesh *mesh);

#endif /* FESOM_SSS_RUNOFF_H */
