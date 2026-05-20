/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 * See fesom_io_config.h banner for context.
 *
 * fesom_io_config.c - line-by-line parser for the override config file.
 * No external deps beyond libc; intended to be linked into both
 * fesom_port and the standalone test_io_config executable.
 */
#include "fesom_io_config.h"

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Local helpers                                                      */
/* ------------------------------------------------------------------ */

static int is_varname_start(int c) { return (isalpha(c) || c == '_'); }
static int is_varname_cont (int c) { return (isalnum(c) || c == '_'); }

static int parse_period(const char *s, fesom_period_kind_t *out)
{
    if (strcmp(s, "step")    == 0) { *out = FESOM_PERIOD_STEP;    return 0; }
    if (strcmp(s, "hourly")  == 0) { *out = FESOM_PERIOD_HOURLY;  return 0; }
    if (strcmp(s, "daily")   == 0) { *out = FESOM_PERIOD_DAILY;   return 0; }
    if (strcmp(s, "monthly") == 0) { *out = FESOM_PERIOD_MONTHLY; return 0; }
    if (strcmp(s, "yearly")  == 0) { *out = FESOM_PERIOD_YEARLY;  return 0; }
    return -1;
}

static int parse_kind(const char *s)
{
    if (strcmp(s, "instant") == 0) return 0;
    if (strcmp(s, "mean")    == 0) return 0;
    return -1;
}

/* Validate that `s` (length len) matches [A-Za-z_][A-Za-z0-9_]*.
 * Length 0 is rejected. */
static int valid_varname(const char *s, size_t len)
{
    if (len == 0) return 0;
    if (!is_varname_start((unsigned char)s[0])) return 0;
    for (size_t i = 1; i < len; ++i)
        if (!is_varname_cont((unsigned char)s[i])) return 0;
    return 1;
}

static char *strndup_local(const char *s, size_t n)
{
    char *out = malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

/* Free an array of heap-allocated strings + the array itself. */
static void free_str_array(char **arr, int n)
{
    if (!arr) return;
    for (int k = 0; k < n; ++k) free(arr[k]);
    free(arr);
}

/* Print a parse error and return -1. */
static int parse_error(const char *path, int lineno, const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "%s:%d: ", path, lineno);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    return -1;
}

/* Comma-split `s` (no nulls inside; len bytes long), reject any
 * whitespace inside. On success, returns the number of fields and writes
 * heap-allocated strings to `out_fields[*]`. Caller frees on success and
 * also on partial failure (function frees its own partial allocations on
 * failure). Returns -1 on syntax error.
 *
 * Grammar enforced:
 *   - No leading or trailing comma  (',foo' / 'foo,')
 *   - No empty field between commas  ('foo,,bar')
 *   - No whitespace anywhere inside  ('foo, bar')
 */
static int comma_split_strict(const char *s, size_t len,
                              char ***out_fields,
                              int    *out_count,
                              const char *path, int lineno,
                              const char *what)
{
    /* Empty input is an error — the LHS or RHS must be non-empty. */
    if (len == 0) {
        return parse_error(path, lineno, "empty %s", what);
    }
    /* Reject comma at the very ends. */
    if (s[0] == ',') {
        return parse_error(path, lineno, "leading ',' in %s", what);
    }
    if (s[len-1] == ',') {
        return parse_error(path, lineno, "trailing ',' in %s", what);
    }

    /* Count commas to size the fields array. */
    int n = 1;
    for (size_t i = 0; i < len; ++i) if (s[i] == ',') ++n;
    char **fields = calloc((size_t)n, sizeof(char *));
    if (!fields) {
        return parse_error(path, lineno, "out of memory in %s", what);
    }

    int idx = 0;
    size_t start = 0;
    for (size_t i = 0; i <= len; ++i) {
        if (i == len || s[i] == ',') {
            size_t flen = i - start;
            if (flen == 0) {
                free_str_array(fields, idx);
                return parse_error(path, lineno, "empty field in %s", what);
            }
            for (size_t j = start; j < i; ++j) {
                if (isspace((unsigned char)s[j])) {
                    for (int k = 0; k < idx; ++k) free(fields[k]);
                    free(fields);
                    return parse_error(path, lineno,
                        "whitespace inside %s ('%.*s'); commas only, no spaces",
                        what, (int)flen, s + start);
                }
            }
            fields[idx] = strndup_local(s + start, flen);
            if (!fields[idx]) {
                free_str_array(fields, idx);
                return parse_error(path, lineno, "out of memory in %s", what);
            }
            ++idx;
            start = i + 1;
        }
    }
    *out_fields = fields;
    *out_count  = idx;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public API                                                         */
/* ------------------------------------------------------------------ */

int fesom_io_config_parse(const char *path, fesom_io_config_t *out)
{
    memset(out, 0, sizeof(*out));

    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "%s: cannot open: %s\n", path, strerror(errno));
        return -1;
    }

    int cap = 16;
    out->entries = calloc((size_t)cap, sizeof(*out->entries));
    if (!out->entries) {
        fclose(f);
        fprintf(stderr, "%s: out of memory\n", path);
        return -1;
    }
    out->n_entries = 0;

    char line[1024];
    int lineno = 0;
    while (fgets(line, sizeof line, f)) {
        ++lineno;

        /* Strip trailing newline + comments. Line is mutated in place. */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = '\0';
        for (size_t i = 0; i < len; ++i) {
            if (line[i] == '#') { line[i] = '\0'; len = i; break; }
        }

        /* Find LHS start (after leading whitespace). */
        size_t i = 0;
        while (i < len && isspace((unsigned char)line[i])) ++i;
        if (i == len) continue;  /* blank or comment-only */

        /* LHS runs to the first whitespace. */
        size_t lhs_start = i;
        while (i < len && !isspace((unsigned char)line[i])) ++i;
        size_t lhs_end = i;
        size_t lhs_len = lhs_end - lhs_start;

        /* Skip the whitespace separator. */
        while (i < len && isspace((unsigned char)line[i])) ++i;
        if (i == len) {
            fclose(f); fesom_io_config_free(out);
            parse_error(path, lineno, "missing cadence list");
            return -1;
        }

        /* RHS runs from `i` to end-of-line; whitespace inside it would
         * have been "trailing whitespace" the previous skip already
         * passed, so any internal whitespace here is a separator that
         * shouldn't exist. comma_split_strict catches it. */
        size_t rhs_start = i;
        size_t rhs_end   = len;
        /* Trim trailing whitespace from RHS. */
        while (rhs_end > rhs_start && isspace((unsigned char)line[rhs_end-1])) --rhs_end;
        size_t rhs_len = rhs_end - rhs_start;

        /* Parse LHS into varnames. */
        char **vars = NULL; int n_vars = 0;
        if (comma_split_strict(line + lhs_start, lhs_len,
                               &vars, &n_vars, path, lineno, "varlist") != 0) {
            fclose(f); fesom_io_config_free(out);
            return -1;
        }
        for (int v = 0; v < n_vars; ++v) {
            if (!valid_varname(vars[v], strlen(vars[v]))) {
                parse_error(path, lineno, "invalid varname '%s'", vars[v]);
                free_str_array(vars, n_vars);
                fclose(f); fesom_io_config_free(out);
                return -1;
            }
        }

        /* Parse RHS into cadence tokens. */
        char **cads = NULL; int n_cads = 0;
        if (comma_split_strict(line + rhs_start, rhs_len,
                               &cads, &n_cads, path, lineno, "cadence list") != 0) {
            free_str_array(vars, n_vars);
            fclose(f); fesom_io_config_free(out);
            return -1;
        }

        /* Convert each cadence token "period[:kind]" into an enum. */
        fesom_period_kind_t *cad_enums = malloc((size_t)n_cads * sizeof(*cad_enums));
        if (!cad_enums) {
            parse_error(path, lineno, "out of memory");
            free_str_array(vars, n_vars);
            free_str_array(cads, n_cads);
            fclose(f); fesom_io_config_free(out);
            return -1;
        }
        for (int c = 0; c < n_cads; ++c) {
            char *colon = strchr(cads[c], ':');
            if (colon) {
                *colon = '\0';
                if (parse_kind(colon + 1) != 0) {
                    parse_error(path, lineno,
                        "unknown kind '%s' (expected 'instant' or 'mean')",
                        colon + 1);
                    free(cad_enums);
                    free_str_array(vars, n_vars);
                    free_str_array(cads, n_cads);
                    fclose(f); fesom_io_config_free(out);
                    return -1;
                }
                /* :kind parsed but ignored — cadence default kinds apply.
                 * See fesom_io_config.h banner. */
            }
            if (parse_period(cads[c], &cad_enums[c]) != 0) {
                parse_error(path, lineno,
                    "unknown period '%s' (expected one of step/hourly/daily/monthly/yearly)",
                    cads[c]);
                free(cad_enums);
                free_str_array(vars, n_vars);
                free_str_array(cads, n_cads);
                fclose(f); fesom_io_config_free(out);
                return -1;
            }
        }
        free_str_array(cads, n_cads);

        /* Emit one entry per (varname, cadence-list). All vars on this
         * line share the same cadence list, but each gets its own copy. */
        for (int v = 0; v < n_vars; ++v) {
            if (out->n_entries == cap) {
                cap *= 2;
                fesom_io_config_entry_t *resized =
                    realloc(out->entries, (size_t)cap * sizeof(*out->entries));
                if (!resized) {
                    parse_error(path, lineno, "out of memory");
                    free(cad_enums);
                    /* vars[0..v-1] already moved into entries; release the
                     * still-owned tail and the array itself. */
                    for (int k = v; k < n_vars; ++k) free(vars[k]);
                    free(vars);
                    fclose(f); fesom_io_config_free(out);
                    return -1;
                }
                out->entries = resized;
            }
            fesom_io_config_entry_t *e = &out->entries[out->n_entries++];
            e->varname    = vars[v];          /* take ownership */
            e->n_cadences = n_cads;
            e->cadences   = malloc((size_t)n_cads * sizeof(fesom_period_kind_t));
            if (!e->cadences) {
                parse_error(path, lineno, "out of memory");
                free(e->varname); e->varname = NULL;
                --out->n_entries;
                free(cad_enums);
                for (int k = v + 1; k < n_vars; ++k) free(vars[k]);
                free(vars);
                fclose(f); fesom_io_config_free(out);
                return -1;
            }
            memcpy(e->cadences, cad_enums, (size_t)n_cads * sizeof(*cad_enums));
        }
        free(cad_enums);
        free(vars);   /* names transferred to entries; the array itself is owned here */
    }
    fclose(f);
    return 0;
}

void fesom_io_config_free(fesom_io_config_t *cfg)
{
    if (!cfg) return;
    for (int i = 0; i < cfg->n_entries; ++i) {
        free(cfg->entries[i].varname);
        free(cfg->entries[i].cadences);
    }
    free(cfg->entries);
    cfg->entries = NULL;
    cfg->n_entries = 0;
}

const fesom_io_config_entry_t *
fesom_io_config_lookup(const fesom_io_config_t *cfg, const char *varname)
{
    if (!cfg) return NULL;
    for (int i = 0; i < cfg->n_entries; ++i) {
        if (strcmp(cfg->entries[i].varname, varname) == 0) {
            return &cfg->entries[i];
        }
    }
    return NULL;
}
