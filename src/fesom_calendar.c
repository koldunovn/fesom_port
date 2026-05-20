/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 * See fesom_calendar.h banner for context.
 *
 * fesom_calendar.c - calendar arithmetic for proleptic Gregorian, 365-day
 * no-leap, and 360-day calendars. All time math goes through a serial-day
 * conversion (days since YYYY-01-01 of year 1 in the chosen calendar) so
 * advance() over arbitrary dt is stable to fractional-second precision.
 */
#include "fesom_calendar.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* days_in_month / days_in_year                                       */
/* ------------------------------------------------------------------ */

static int is_gregorian_leap(int year)
{
    if ((year % 4) != 0)   return 0;
    if ((year % 100) != 0) return 1;
    if ((year % 400) != 0) return 0;
    return 1;
}

int fesom_calendar_days_in_month(fesom_calendar_kind_t kind, int year, int month)
{
    static const int dbm_31_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month < 1 || month > 12) return 0;
    switch (kind) {
        case FESOM_CAL_GREGORIAN:
            if (month == 2 && is_gregorian_leap(year)) return 29;
            return dbm_31_days[month - 1];
        case FESOM_CAL_NOLEAP_365:
            return dbm_31_days[month - 1];
        case FESOM_CAL_360_DAY:
            return 30;
    }
    return 0;
}

int fesom_calendar_days_in_year(fesom_calendar_kind_t kind, int year)
{
    switch (kind) {
        case FESOM_CAL_GREGORIAN:    return is_gregorian_leap(year) ? 366 : 365;
        case FESOM_CAL_NOLEAP_365:   return 365;
        case FESOM_CAL_360_DAY:      return 360;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* serial-day conversions                                             */
/*                                                                    */
/* For each kind, ymd <-> serial_day where serial_day = 0 means       */
/* (year=1, month=1, day=1) of that calendar. Strictly monotonic.     */
/* ------------------------------------------------------------------ */

/* Howard-Hinnant proleptic-Gregorian formula. Returns days from
 * 1970-01-01; pre-1970 years yield negative values. The public API only
 * ever takes differences between two same-kind serials, so the absolute
 * value of the epoch doesn't matter. */
static int64_t gregorian_serial_day(int y, int m, int d)
{
    y -= (m <= 2);
    int64_t era = (y >= 0 ? y : y - 399) / 400;
    int64_t yoe = y - era * 400;
    int64_t doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
    int64_t doe = yoe * 365 + yoe/4 - yoe/100 + doy;
    return era * 146097 + doe - 719468;
}

static void gregorian_civil_from_serial(int64_t s, int *y, int *m, int *d)
{
    int64_t z = s + 719468;
    int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    int64_t doe = z - era * 146097;
    int64_t yoe = (doe - doe/1460 + doe/36524 - doe/146096) / 365;
    int yy  = (int)(yoe + era * 400);
    int doy = (int)(doe - (365*yoe + yoe/4 - yoe/100));
    int mp  = (5*doy + 2) / 153;
    *d = doy - (153*mp + 2)/5 + 1;
    *m = mp < 10 ? mp + 3 : mp - 9;
    *y = yy + (*m <= 2 ? 1 : 0);
}

static const int dbm_31_days_offset[13] = {
    0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
};

static int64_t noleap_serial_day(int y, int m, int d)
{
    return (int64_t)(y - 1) * 365 + dbm_31_days_offset[m - 1] + (d - 1);
}

static void noleap_civil_from_serial(int64_t s, int *y, int *m, int *d)
{
    /* Floor-div to handle negatives correctly even though we only
     * expect positive years in practice. */
    int64_t yy = s / 365;
    int64_t rem = s - yy * 365;
    if (rem < 0) { rem += 365; yy -= 1; }
    *y = (int)(yy + 1);
    int month = 1;
    while (month < 12 && rem >= dbm_31_days_offset[month]) month++;
    *m = month;
    *d = (int)(rem - dbm_31_days_offset[month - 1] + 1);
}

static int64_t day360_serial_day(int y, int m, int d)
{
    return (int64_t)(y - 1) * 360 + (int64_t)(m - 1) * 30 + (d - 1);
}

static void day360_civil_from_serial(int64_t s, int *y, int *m, int *d)
{
    int64_t yy = s / 360;
    int64_t rem = s - yy * 360;
    if (rem < 0) { rem += 360; yy -= 1; }
    *y = (int)(yy + 1);
    *m = (int)(rem / 30) + 1;
    *d = (int)(rem % 30) + 1;
}

static int64_t serial_day(fesom_calendar_kind_t k, int y, int m, int d)
{
    switch (k) {
        case FESOM_CAL_GREGORIAN:    return gregorian_serial_day(y, m, d);
        case FESOM_CAL_NOLEAP_365:   return noleap_serial_day(y, m, d);
        case FESOM_CAL_360_DAY:      return day360_serial_day(y, m, d);
    }
    return 0;
}

static void civil_from_serial(fesom_calendar_kind_t k, int64_t s,
                              int *y, int *m, int *d)
{
    switch (k) {
        case FESOM_CAL_GREGORIAN:    gregorian_civil_from_serial(s, y, m, d); return;
        case FESOM_CAL_NOLEAP_365:   noleap_civil_from_serial(s, y, m, d);    return;
        case FESOM_CAL_360_DAY:      day360_civil_from_serial(s, y, m, d);    return;
    }
}

/* ------------------------------------------------------------------ */
/* Public API                                                         */
/* ------------------------------------------------------------------ */

void fesom_calendar_init(fesom_calendar_t *cal,
                         fesom_calendar_kind_t kind,
                         int year, int month, int day)
{
    cal->kind   = kind;
    cal->year   = year;
    cal->month  = month;
    cal->day    = day;
    cal->hour   = 0;
    cal->minute = 0;
    cal->second = 0.0;
    cal->origin_year  = year;
    cal->origin_month = month;
    cal->origin_day   = day;
}

void fesom_calendar_advance(fesom_calendar_t *cal, double dt_seconds)
{
    int64_t sd = serial_day(cal->kind, cal->year, cal->month, cal->day);
    double total = (double)sd * 86400.0
                 + (double)cal->hour * 3600.0
                 + (double)cal->minute * 60.0
                 + cal->second
                 + dt_seconds;

    /* Split into integer-seconds + sub-second fraction.
     * floor() ensures negative-dt cases would still bucket correctly
     * even though we only call with dt >= 0 in practice. */
    double tot_floor = floor(total);
    double frac = total - tot_floor;
    int64_t total_int = (int64_t)tot_floor;

    int64_t new_sd = total_int / 86400;
    int64_t rem    = total_int - new_sd * 86400;
    if (rem < 0) { rem += 86400; new_sd -= 1; }

    cal->hour   = (int)(rem / 3600);
    cal->minute = (int)((rem % 3600) / 60);
    cal->second = (double)(rem % 60) + frac;
    civil_from_serial(cal->kind, new_sd, &cal->year, &cal->month, &cal->day);
}

double fesom_calendar_seconds_since_origin(const fesom_calendar_t *cal)
{
    int64_t sd_now = serial_day(cal->kind, cal->year, cal->month, cal->day);
    int64_t sd_org = serial_day(cal->kind, cal->origin_year,
                                cal->origin_month, cal->origin_day);
    int64_t day_diff = sd_now - sd_org;
    return (double)day_diff * 86400.0
         + (double)cal->hour * 3600.0
         + (double)cal->minute * 60.0
         + cal->second;
}

int fesom_calendar_day_of_year(const fesom_calendar_t *cal)
{
    int64_t sd_jan1 = serial_day(cal->kind, cal->year, 1, 1);
    int64_t sd_now  = serial_day(cal->kind, cal->year, cal->month, cal->day);
    return (int)(sd_now - sd_jan1) + 1;
}

int fesom_calendar_crossed(const fesom_calendar_t *prev,
                           const fesom_calendar_t *curr,
                           fesom_period_kind_t period)
{
    switch (period) {
        case FESOM_PERIOD_STEP:
            return 1;
        case FESOM_PERIOD_HOURLY:
            return prev->year   != curr->year   ||
                   prev->month  != curr->month  ||
                   prev->day    != curr->day    ||
                   prev->hour   != curr->hour;
        case FESOM_PERIOD_DAILY:
            return prev->year   != curr->year   ||
                   prev->month  != curr->month  ||
                   prev->day    != curr->day;
        case FESOM_PERIOD_MONTHLY:
            return prev->year   != curr->year   ||
                   prev->month  != curr->month;
        case FESOM_PERIOD_YEARLY:
            return prev->year   != curr->year;
    }
    return 0;
}

void fesom_calendar_period_window(const fesom_calendar_t *anchor,
                                  fesom_period_kind_t period,
                                  double *t_start, double *t_mid, double *t_end)
{
    int64_t sd_org = serial_day(anchor->kind, anchor->origin_year,
                                anchor->origin_month, anchor->origin_day);

    /* t_start = seconds-since-origin at the start of the period containing
     * `anchor`; t_end = start of next period; t_mid = arithmetic midpoint.
     * STEP collapses to a degenerate point at the anchor itself. */
    int64_t sd_anchor = serial_day(anchor->kind, anchor->year, anchor->month, anchor->day);

    double seconds_in_day_anchor = (double)anchor->hour * 3600.0
                                 + (double)anchor->minute * 60.0
                                 + anchor->second;
    double tA = (double)(sd_anchor - sd_org) * 86400.0 + seconds_in_day_anchor;

    double tS = 0.0, tE = 0.0;

    switch (period) {
        case FESOM_PERIOD_STEP:
            tS = tA;
            tE = tA;
            break;
        case FESOM_PERIOD_HOURLY: {
            double seconds_into_hour = (double)anchor->minute * 60.0 + anchor->second;
            tS = tA - seconds_into_hour;
            tE = tS + 3600.0;
            break;
        }
        case FESOM_PERIOD_DAILY:
            tS = tA - seconds_in_day_anchor;
            tE = tS + 86400.0;
            break;
        case FESOM_PERIOD_MONTHLY: {
            int dim = fesom_calendar_days_in_month(anchor->kind,
                                                   anchor->year, anchor->month);
            int64_t sd_first_of_month = serial_day(anchor->kind,
                                                   anchor->year, anchor->month, 1);
            tS = (double)(sd_first_of_month - sd_org) * 86400.0;
            tE = tS + (double)dim * 86400.0;
            break;
        }
        case FESOM_PERIOD_YEARLY: {
            int diy = fesom_calendar_days_in_year(anchor->kind, anchor->year);
            int64_t sd_first_of_year = serial_day(anchor->kind, anchor->year, 1, 1);
            tS = (double)(sd_first_of_year - sd_org) * 86400.0;
            tE = tS + (double)diy * 86400.0;
            break;
        }
    }

    if (t_start) *t_start = tS;
    if (t_end)   *t_end   = tE;
    if (t_mid)   *t_mid   = 0.5 * (tS + tE);
}

size_t fesom_calendar_units_string(const fesom_calendar_t *cal,
                                   char *buf, size_t buflen)
{
    int n = snprintf(buf, buflen, "seconds since %04d-%02d-%02d 00:00:00",
                     cal->origin_year, cal->origin_month, cal->origin_day);
    if (n < 0) {
        if (buflen > 0) buf[0] = '\0';
        return 0;
    }
    return (size_t)n < buflen ? (size_t)n : buflen - 1;
}

const char *fesom_calendar_cf_name(fesom_calendar_kind_t kind)
{
    switch (kind) {
        case FESOM_CAL_GREGORIAN:    return "proleptic_gregorian";
        case FESOM_CAL_NOLEAP_365:   return "365_day";
        case FESOM_CAL_360_DAY:      return "360_day";
    }
    return "unknown";
}
