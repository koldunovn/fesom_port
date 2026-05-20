/*
 * Unit tests for fesom_io_config parser. Builds against
 * fesom_io_config.c + fesom_calendar.c (no MPI, no netCDF).
 */
#include "fesom_io_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int n_passed = 0, n_failed = 0;

#define CHECK(expr)                                                          \
    do {                                                                     \
        if (expr) ++n_passed;                                                \
        else      {                                                          \
            ++n_failed;                                                      \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr);  \
        }                                                                    \
    } while (0)

/* Write `content` to a fresh temp file and return its path (heap, caller frees). */
static char *write_temp_config(const char *content)
{
    char path[64];
    snprintf(path, sizeof path, "/tmp/test_io_config_%d.txt", getpid());
    FILE *f = fopen(path, "w");
    if (!f) return NULL;
    fwrite(content, 1, strlen(content), f);
    fclose(f);
    return strdup(path);
}

/* ------------------------------------------------------------------ */
static void test_valid_simple(void)
{
    const char *src =
        "# This is a comment\n"
        "\n"
        "sst    step,daily,monthly\n"
        "ssh    monthly,yearly\n"
        "temp,salt   daily,monthly,yearly\n";
    char *path = write_temp_config(src);
    fesom_io_config_t cfg;
    CHECK(fesom_io_config_parse(path, &cfg) == 0);
    CHECK(cfg.n_entries == 4);   /* sst, ssh, temp, salt */

    const fesom_io_config_entry_t *e = fesom_io_config_lookup(&cfg, "sst");
    CHECK(e != NULL);
    if (e) {
        CHECK(e->n_cadences == 3);
        CHECK(e->cadences[0] == FESOM_PERIOD_STEP);
        CHECK(e->cadences[1] == FESOM_PERIOD_DAILY);
        CHECK(e->cadences[2] == FESOM_PERIOD_MONTHLY);
    }
    e = fesom_io_config_lookup(&cfg, "temp");
    CHECK(e != NULL && e->n_cadences == 3);
    e = fesom_io_config_lookup(&cfg, "salt");
    CHECK(e != NULL && e->n_cadences == 3);
    CHECK(fesom_io_config_lookup(&cfg, "u") == NULL);  /* not in config */

    fesom_io_config_free(&cfg);
    unlink(path); free(path);
}

static void test_kind_suffix_parsed_then_ignored(void)
{
    /* :instant and :mean are accepted by the grammar; the orchestrator
     * (separately) ignores the suffix and uses cadence-default kinds. */
    const char *src = "ssh    monthly:instant,yearly:mean\n";
    char *path = write_temp_config(src);
    fesom_io_config_t cfg;
    CHECK(fesom_io_config_parse(path, &cfg) == 0);
    CHECK(cfg.n_entries == 1);
    const fesom_io_config_entry_t *e = fesom_io_config_lookup(&cfg, "ssh");
    CHECK(e != NULL);
    if (e) {
        CHECK(e->n_cadences == 2);
        CHECK(e->cadences[0] == FESOM_PERIOD_MONTHLY);
        CHECK(e->cadences[1] == FESOM_PERIOD_YEARLY);
    }
    fesom_io_config_free(&cfg);
    unlink(path); free(path);
}

static void test_blanks_and_comments_only(void)
{
    const char *src = "# all comments\n\n   \n# another\n";
    char *path = write_temp_config(src);
    fesom_io_config_t cfg;
    CHECK(fesom_io_config_parse(path, &cfg) == 0);
    CHECK(cfg.n_entries == 0);
    fesom_io_config_free(&cfg);
    unlink(path); free(path);
}

static void test_inline_comment(void)
{
    const char *src = "sst   monthly  # daily comes later\n";
    char *path = write_temp_config(src);
    fesom_io_config_t cfg;
    CHECK(fesom_io_config_parse(path, &cfg) == 0);
    const fesom_io_config_entry_t *e = fesom_io_config_lookup(&cfg, "sst");
    CHECK(e != NULL && e->n_cadences == 1);
    if (e) CHECK(e->cadences[0] == FESOM_PERIOD_MONTHLY);
    fesom_io_config_free(&cfg);
    unlink(path); free(path);
}

/* ------------------------------------------------------------------ */
/* Negative cases — must reject and the line number must appear in    */
/* stderr. We only check the return code here; visual stderr is OK    */
/* for now. (Test runs under ctest --output-on-failure.)              */
/* ------------------------------------------------------------------ */
static void expect_reject(const char *src, const char *label)
{
    char *path = write_temp_config(src);
    fesom_io_config_t cfg;
    fprintf(stderr, "  [expected error from '%s']:\n  ", label);
    int rc = fesom_io_config_parse(path, &cfg);
    if (rc == 0) {
        ++n_failed;
        fprintf(stderr, "FAIL: input '%s' was accepted, expected reject\n", label);
        fesom_io_config_free(&cfg);
    } else {
        ++n_passed;
    }
    unlink(path); free(path);
}

static void test_negative_cases(void)
{
    expect_reject("temp, salt  daily\n",        "whitespace inside LHS list");
    expect_reject(",foo  daily\n",              "leading comma in LHS");
    expect_reject("foo,  daily\n",              "trailing comma in LHS");
    expect_reject("foo  daily,\n",              "trailing comma in RHS");
    expect_reject("foo  daily, monthly\n",      "whitespace inside RHS list");
    expect_reject("foo\n",                      "missing RHS");
    expect_reject("foo  hourly:bogus\n",        "invalid kind");
    expect_reject("foo  bogus_period\n",        "invalid period");
    expect_reject("123foo  daily\n",            "varname starts with digit");
    expect_reject("foo,,bar  daily\n",          "empty field in LHS");
}

/* ------------------------------------------------------------------ */
int main(void)
{
    test_valid_simple();
    test_kind_suffix_parsed_then_ignored();
    test_blanks_and_comments_only();
    test_inline_comment();
    test_negative_cases();

    fprintf(stderr, "test_io_config: %d passed, %d failed\n", n_passed, n_failed);
    return n_failed == 0 ? 0 : 1;
}
