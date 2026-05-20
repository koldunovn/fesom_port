/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 *
 * The port-fidelity rule (CLAUDE.md / feedback_port_fidelity.md) requires
 * line-by-line porting from Fortran. The IO subsystem is the deliberate
 * exception, with explicit user authorisation. See:
 *   docs/plans/20260425-io-system.md
 *   feedback_io_system_design.md (once written)
 *
 * Do NOT "fix" this file by replacing it with a literal port of FESOM2's
 * io_meandata.F90 / io_netcdf_module.F90 / io_restart.F90 etc. - that is
 * specifically what we chose not to do.
 *
 * fesom_calendar.h - the model's single source of truth for "what date is
 * it now". Replaces the three scattered date implementations: JRA55's
 * julday, the static days_in_month[12] non-leap array in fesom_main.c,
 * and the month_now/update_monthly_flag detection.
 */
#ifndef FESOM_CALENDAR_H
#define FESOM_CALENDAR_H

#include <stddef.h>

typedef enum {
    FESOM_CAL_GREGORIAN  = 0,   /* proleptic gregorian, leap years */
    FESOM_CAL_NOLEAP_365 = 1,   /* every Feb is 28 days            */
    FESOM_CAL_360_DAY    = 2,   /* every month is 30 days          */
} fesom_calendar_kind_t;

typedef enum {
    FESOM_PERIOD_STEP    = 0,
    FESOM_PERIOD_HOURLY  = 1,
    FESOM_PERIOD_DAILY   = 2,
    FESOM_PERIOD_MONTHLY = 3,
    FESOM_PERIOD_YEARLY  = 4,
} fesom_period_kind_t;

typedef struct {
    fesom_calendar_kind_t kind;
    int    year;        /* full year, e.g. 1948                          */
    int    month;       /* 1..12                                          */
    int    day;         /* 1..days_in_month(kind, year, month)            */
    int    hour;        /* 0..23                                          */
    int    minute;      /* 0..59                                          */
    double second;      /* 0.0 .. <60.0; fractional for sub-second dt     */

    /* CF time-axis anchor. "seconds since YYYY-MM-DD 00:00:00".          */
    int    origin_year;
    int    origin_month;
    int    origin_day;
} fesom_calendar_t;

/* Initialise to (start_year, start_month, start_day) at 00:00:00.0.
 * The CF origin is set to the same date so seconds_since_origin starts
 * at 0 and time stamps grow monotonically from run start. */
void fesom_calendar_init(fesom_calendar_t *cal,
                         fesom_calendar_kind_t kind,
                         int year, int month, int day);

/* Advance the wall clock by dt_seconds. Robust to non-integer dt
 * (e.g. dt=2700, 500, fractional). Internally accumulates a 64-bit
 * stable seconds count then re-derives y/m/d/h/m/s, so drift over
 * long runs is bounded only by double-precision rounding on dt. */
void fesom_calendar_advance(fesom_calendar_t *cal, double dt_seconds);

/* Time since CF origin in seconds. */
double fesom_calendar_seconds_since_origin(const fesom_calendar_t *cal);

/* Days in month under this calendar kind. (year,month) 1-based;
 * year ignored for 365-day and 360-day kinds. */
int fesom_calendar_days_in_month(fesom_calendar_kind_t kind, int year, int month);

/* Days in year under this calendar kind. */
int fesom_calendar_days_in_year(fesom_calendar_kind_t kind, int year);

/* Day-of-year (1..365|366) for the current calendar moment. */
int fesom_calendar_day_of_year(const fesom_calendar_t *cal);

/* Return non-zero iff the last advance() crossed at least one boundary
 * of the given period (i.e., prev and curr lie in different windows).
 *
 * Comparing prev vs curr instances rather than modular arithmetic is
 * what makes variable-length months work cleanly; we just ask "did the
 * (year,month) tuple change?" instead of computing 30-day windows. */
int fesom_calendar_crossed(const fesom_calendar_t *prev,
                           const fesom_calendar_t *curr,
                           fesom_period_kind_t period);

/* Period bounds and midpoint (seconds since CF origin) for the period
 * window that `anchor` lies in. STEP returns a degenerate
 * (t_anchor, t_anchor, t_anchor). HOURLY/DAILY/MONTHLY/YEARLY return
 * the [start, end) of the calendar period containing `anchor`, with
 * t_mid = 0.5*(start+end) for the CF time variable. */
void fesom_calendar_period_window(const fesom_calendar_t *anchor,
                                  fesom_period_kind_t period,
                                  double *t_start, double *t_mid, double *t_end);

/* CF "units" attribute string, e.g. "seconds since 1948-01-01 00:00:00".
 * Writes up to buflen-1 chars + NUL; returns chars written excluding NUL. */
size_t fesom_calendar_units_string(const fesom_calendar_t *cal,
                                   char *buf, size_t buflen);

/* CF "calendar" attribute string. Pointer to static storage; do not free. */
const char *fesom_calendar_cf_name(fesom_calendar_kind_t kind);

#endif /* FESOM_CALENDAR_H */
