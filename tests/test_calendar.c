/*
 * Standalone unit tests for fesom_calendar. Build via the CMake target
 * `test_calendar`; run via `ctest -V`.
 *
 * Failure mode: assert() abort with file:line — easier to bisect than
 * a framework-formatted dump in an environment without a unit-test
 * framework set up.
 */
#include "fesom_calendar.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define APPROX_EQ(a, b, tol) (fabs((a) - (b)) <= (tol))

static int n_passed = 0;
static int n_failed = 0;

#define CHECK(expr)                                                          \
    do {                                                                     \
        if (expr) { ++n_passed; }                                            \
        else      {                                                          \
            ++n_failed;                                                      \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr);  \
        }                                                                    \
    } while (0)

/* ------------------------------------------------------------------ */
static void test_days_in_month_gregorian(void)
{
    /* Leap year rule: divisible by 4, except centuries not divisible by 400. */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 2000, 2) == 29); /* 400 → leap */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 2004, 2) == 29); /* /4 → leap */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 2100, 2) == 28); /* century, not /400 */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 1900, 2) == 28); /* century, not /400 */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 1948, 2) == 29); /* /4 leap */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 1948, 1) == 31);
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 1948, 4) == 30);
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_GREGORIAN, 1948, 12) == 31);
    CHECK(fesom_calendar_days_in_year(FESOM_CAL_GREGORIAN, 2000) == 366);
    CHECK(fesom_calendar_days_in_year(FESOM_CAL_GREGORIAN, 1900) == 365);
    CHECK(fesom_calendar_days_in_year(FESOM_CAL_GREGORIAN, 1948) == 366);
}

static void test_days_in_month_365_360(void)
{
    /* 365-day: Feb is always 28 regardless of year. */
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_NOLEAP_365, 2000, 2) == 28);
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_NOLEAP_365, 2004, 2) == 28);
    CHECK(fesom_calendar_days_in_month(FESOM_CAL_NOLEAP_365, 1948, 2) == 28);
    CHECK(fesom_calendar_days_in_year(FESOM_CAL_NOLEAP_365, 2000) == 365);
    CHECK(fesom_calendar_days_in_year(FESOM_CAL_NOLEAP_365, 1948) == 365);
    /* 360-day: every month is 30 regardless. */
    for (int m = 1; m <= 12; ++m) {
        CHECK(fesom_calendar_days_in_month(FESOM_CAL_360_DAY, 2000, m) == 30);
    }
    CHECK(fesom_calendar_days_in_year(FESOM_CAL_360_DAY, 2000) == 360);
}

/* ------------------------------------------------------------------ */
static void test_advance_day_boundary(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    /* Advance to 86399.0 = end of day. */
    fesom_calendar_advance(&cal, 86399.0);
    CHECK(cal.year == 1948 && cal.month == 1 && cal.day == 1);
    CHECK(cal.hour == 23 && cal.minute == 59 && APPROX_EQ(cal.second, 59.0, 1e-9));
    /* Cross midnight. */
    fesom_calendar_advance(&cal, 1.0);
    CHECK(cal.year == 1948 && cal.month == 1 && cal.day == 2);
    CHECK(cal.hour == 0 && cal.minute == 0 && APPROX_EQ(cal.second, 0.0, 1e-9));
}

static void test_advance_month_boundary_gregorian(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 31);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1948 && cal.month == 2 && cal.day == 1);
    /* Feb 1948 has 29 days (leap). */
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 2, 28);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1948 && cal.month == 2 && cal.day == 29);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1948 && cal.month == 3 && cal.day == 1);
    /* Feb 1949 has 28 days (not leap). */
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1949, 2, 28);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1949 && cal.month == 3 && cal.day == 1);
}

static void test_advance_year_boundary(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 12, 31);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1949 && cal.month == 1 && cal.day == 1);
}

static void test_advance_noleap_skips_feb29(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_NOLEAP_365, 1948, 2, 28);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1948 && cal.month == 3 && cal.day == 1);  /* No Feb 29 */
}

static void test_advance_360_all_30(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_360_DAY, 1948, 2, 30);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1948 && cal.month == 3 && cal.day == 1);
    fesom_calendar_init(&cal, FESOM_CAL_360_DAY, 1948, 12, 30);
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(cal.year == 1949 && cal.month == 1 && cal.day == 1);
}

static void test_advance_non_integer_dt(void)
{
    /* dt = 500s: typical CORE2 ocean step. */
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    for (int i = 0; i < 173; ++i) fesom_calendar_advance(&cal, 500.0);
    /* 173 * 500 = 86500s = 24h 1m 40s past 1948-01-01 00:00:00. */
    CHECK(cal.year == 1948 && cal.month == 1 && cal.day == 2);
    CHECK(cal.hour == 0 && cal.minute == 1 && APPROX_EQ(cal.second, 40.0, 1e-6));

    /* dt = 2700s: Fortran CORE2 with sea ice. */
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    for (int i = 0; i < 32; ++i) fesom_calendar_advance(&cal, 2700.0);
    /* 32 * 2700 = 86400s = exactly one day. */
    CHECK(cal.year == 1948 && cal.month == 1 && cal.day == 2);
    CHECK(cal.hour == 0 && cal.minute == 0 && APPROX_EQ(cal.second, 0.0, 1e-6));

    /* Fractional dt. */
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    fesom_calendar_advance(&cal, 0.5);
    CHECK(APPROX_EQ(cal.second, 0.5, 1e-9));
    fesom_calendar_advance(&cal, 0.25);
    CHECK(APPROX_EQ(cal.second, 0.75, 1e-9));
}

/* ------------------------------------------------------------------ */
static void test_crossed_per_period(void)
{
    fesom_calendar_t prev, curr;
    /* Same instant: nothing crossed except STEP. */
    fesom_calendar_init(&prev, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    curr = prev;
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_STEP)    == 1);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_HOURLY)  == 0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_DAILY)   == 0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_MONTHLY) == 0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_YEARLY)  == 0);

    /* Hour boundary. */
    prev = curr;
    fesom_calendar_advance(&curr, 3600.0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_HOURLY)  == 1);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_DAILY)   == 0);

    /* Day boundary (start of 1948-01-02). */
    fesom_calendar_init(&prev, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    fesom_calendar_advance(&prev, 86399.0);
    curr = prev;
    fesom_calendar_advance(&curr, 1.0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_DAILY)   == 1);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_MONTHLY) == 0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_YEARLY)  == 0);

    /* Month boundary. */
    fesom_calendar_init(&prev, FESOM_CAL_GREGORIAN, 1948, 1, 31);
    fesom_calendar_advance(&prev, 86399.0);
    curr = prev;
    fesom_calendar_advance(&curr, 1.0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_MONTHLY) == 1);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_YEARLY)  == 0);

    /* Year boundary. */
    fesom_calendar_init(&prev, FESOM_CAL_GREGORIAN, 1948, 12, 31);
    fesom_calendar_advance(&prev, 86399.0);
    curr = prev;
    fesom_calendar_advance(&curr, 1.0);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_YEARLY)  == 1);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_MONTHLY) == 1);
    CHECK(fesom_calendar_crossed(&prev, &curr, FESOM_PERIOD_DAILY)   == 1);
}

/* ------------------------------------------------------------------ */
static void test_period_window_monthly(void)
{
    /* Anchor mid-month: window covers full calendar month, not partial. */
    fesom_calendar_t anchor;
    fesom_calendar_init(&anchor, FESOM_CAL_GREGORIAN, 1948, 3, 15);
    fesom_calendar_advance(&anchor, 12.0 * 3600.0);  /* 12:00 noon */
    double tS, tM, tE;
    fesom_calendar_period_window(&anchor, FESOM_PERIOD_MONTHLY, &tS, &tM, &tE);
    /* origin = 1948-03-15; tS = start of March 1948 = -14 days. */
    CHECK(APPROX_EQ(tS, -14.0 * 86400.0, 1e-6));
    /* March has 31 days → tE = -14 + 31 = +17 days. */
    CHECK(APPROX_EQ(tE, 17.0 * 86400.0, 1e-6));
    /* Midpoint = 0.5*(tS+tE). */
    CHECK(APPROX_EQ(tM, 0.5 * (tS + tE), 1e-6));
}

static void test_period_window_anchor_on_boundary(void)
{
    /* Anchor exactly on period boundary (1948-04-01 00:00:00).
     * The window for MONTHLY containing this anchor is April 1948. */
    fesom_calendar_t anchor;
    fesom_calendar_init(&anchor, FESOM_CAL_GREGORIAN, 1948, 4, 1);  /* origin */
    double tS, tM, tE;
    fesom_calendar_period_window(&anchor, FESOM_PERIOD_MONTHLY, &tS, &tM, &tE);
    CHECK(APPROX_EQ(tS, 0.0, 1e-6));
    CHECK(APPROX_EQ(tE, 30.0 * 86400.0, 1e-6));  /* April has 30 days */
}

static void test_period_window_360_feb_crossing(void)
{
    /* 360-day calendar: February still has 30 days, no leap-day issues. */
    fesom_calendar_t anchor;
    fesom_calendar_init(&anchor, FESOM_CAL_360_DAY, 1948, 2, 15);
    double tS, tM, tE;
    fesom_calendar_period_window(&anchor, FESOM_PERIOD_MONTHLY, &tS, &tM, &tE);
    CHECK(APPROX_EQ(tS, -14.0 * 86400.0, 1e-6));
    CHECK(APPROX_EQ(tE, 16.0 * 86400.0, 1e-6));  /* +30 days */
}

static void test_period_window_daily_hourly(void)
{
    fesom_calendar_t anchor;
    fesom_calendar_init(&anchor, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    fesom_calendar_advance(&anchor, 13.0 * 3600.0 + 27.0 * 60.0 + 4.5);
    double tS, tM, tE;

    fesom_calendar_period_window(&anchor, FESOM_PERIOD_DAILY, &tS, &tM, &tE);
    CHECK(APPROX_EQ(tS, 0.0, 1e-6));
    CHECK(APPROX_EQ(tE, 86400.0, 1e-6));

    fesom_calendar_period_window(&anchor, FESOM_PERIOD_HOURLY, &tS, &tM, &tE);
    /* 13:27:04.5 → start of 13:00, end of 14:00. */
    CHECK(APPROX_EQ(tS, 13.0 * 3600.0, 1e-6));
    CHECK(APPROX_EQ(tE, 14.0 * 3600.0, 1e-6));
    CHECK(APPROX_EQ(tM, 13.5 * 3600.0, 1e-6));
}

static void test_period_window_yearly(void)
{
    fesom_calendar_t anchor;
    fesom_calendar_init(&anchor, FESOM_CAL_GREGORIAN, 1948, 7, 4);
    double tS, tM, tE;
    fesom_calendar_period_window(&anchor, FESOM_PERIOD_YEARLY, &tS, &tM, &tE);
    /* origin = 1948-07-04; start of 1948 = -185 days (Jan 1 + Jan(31) + Feb(29) + Mar(31) + Apr(30) + May(31) + Jun(30) + 4 = 186, so 1948-07-04 is day 186, start of year is 185 days back). */
    CHECK(APPROX_EQ(tS, -185.0 * 86400.0, 1e-6));
    CHECK(APPROX_EQ(tE, 181.0 * 86400.0, 1e-6));  /* leap year: 366 days total */
}

static void test_period_window_step_degenerate(void)
{
    fesom_calendar_t anchor;
    fesom_calendar_init(&anchor, FESOM_CAL_GREGORIAN, 1948, 6, 15);
    fesom_calendar_advance(&anchor, 17.5);
    double tS, tM, tE;
    fesom_calendar_period_window(&anchor, FESOM_PERIOD_STEP, &tS, &tM, &tE);
    CHECK(APPROX_EQ(tS, tM, 1e-9));
    CHECK(APPROX_EQ(tM, tE, 1e-9));
    /* origin = 1948-06-15; tA = days(0) + 17.5 seconds = 17.5. */
    CHECK(APPROX_EQ(tS, 17.5, 1e-9));
}

/* ------------------------------------------------------------------ */
static void test_seconds_since_origin(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    CHECK(APPROX_EQ(fesom_calendar_seconds_since_origin(&cal), 0.0, 1e-9));
    fesom_calendar_advance(&cal, 86400.0);
    CHECK(APPROX_EQ(fesom_calendar_seconds_since_origin(&cal), 86400.0, 1e-9));
    fesom_calendar_advance(&cal, 86400.0 * 365.0);  /* a leap year is 366 → cross 1949 */
    /* After 366 days (the Jan 2nd 1948 advance + 365 more): we should be at
     * 1949-01-01 because 1948 is leap (366 days from Jan 1 → Dec 31 + 1 → next-year Jan 1). */
    CHECK(cal.year == 1949 && cal.month == 1 && cal.day == 1);
    CHECK(APPROX_EQ(fesom_calendar_seconds_since_origin(&cal), 366.0 * 86400.0, 1e-9));
}

static void test_day_of_year(void)
{
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    CHECK(fesom_calendar_day_of_year(&cal) == 1);
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 12, 31);
    CHECK(fesom_calendar_day_of_year(&cal) == 366);  /* leap */
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1949, 12, 31);
    CHECK(fesom_calendar_day_of_year(&cal) == 365);
    fesom_calendar_init(&cal, FESOM_CAL_360_DAY, 1948, 12, 30);
    CHECK(fesom_calendar_day_of_year(&cal) == 360);
}

/* ------------------------------------------------------------------ */
static void test_cf_strings(void)
{
    char buf[64];
    fesom_calendar_t cal;
    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 1948, 1, 1);
    size_t n = fesom_calendar_units_string(&cal, buf, sizeof buf);
    CHECK(n > 0);
    CHECK(strcmp(buf, "seconds since 1948-01-01 00:00:00") == 0);

    fesom_calendar_init(&cal, FESOM_CAL_GREGORIAN, 2026, 4, 25);
    fesom_calendar_units_string(&cal, buf, sizeof buf);
    CHECK(strcmp(buf, "seconds since 2026-04-25 00:00:00") == 0);

    CHECK(strcmp(fesom_calendar_cf_name(FESOM_CAL_GREGORIAN),  "proleptic_gregorian") == 0);
    CHECK(strcmp(fesom_calendar_cf_name(FESOM_CAL_NOLEAP_365), "365_day")             == 0);
    CHECK(strcmp(fesom_calendar_cf_name(FESOM_CAL_360_DAY),    "360_day")             == 0);
}

/* ------------------------------------------------------------------ */
int main(void)
{
    test_days_in_month_gregorian();
    test_days_in_month_365_360();
    test_advance_day_boundary();
    test_advance_month_boundary_gregorian();
    test_advance_year_boundary();
    test_advance_noleap_skips_feb29();
    test_advance_360_all_30();
    test_advance_non_integer_dt();
    test_crossed_per_period();
    test_period_window_monthly();
    test_period_window_anchor_on_boundary();
    test_period_window_360_feb_crossing();
    test_period_window_daily_hourly();
    test_period_window_yearly();
    test_period_window_step_degenerate();
    test_seconds_since_origin();
    test_day_of_year();
    test_cf_strings();

    fprintf(stderr, "test_calendar: %d passed, %d failed\n", n_passed, n_failed);
    return n_failed == 0 ? 0 : 1;
}
