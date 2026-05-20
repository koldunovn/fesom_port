/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 * See fesom_io_stream.h banner for context.
 *
 * fesom_io_config.h - tiny parser for the override config file. Loaded
 * once at fesom_io_init when FESOM_IO_CONFIG points to a path. Each
 * non-comment line names one or more variables and the cadences they
 * should be written at; the orchestrator uses this to override the
 * compiled-in monthly-only default per variable.
 *
 * Grammar (also reproduced in docs/plans/20260425-io-system.md):
 *
 *     line     := comment | blank | entry
 *     comment  := '#' [^\n]*
 *     blank    := [ \t]*
 *     entry    := varlist  WSEP  cadlist
 *     varlist  := varname (',' varname)*       ; NO whitespace inside
 *     cadlist  := cadence (',' cadence)*       ; NO whitespace inside
 *     varname  := [A-Za-z_][A-Za-z0-9_]*
 *     cadence  := period (':' kind)?
 *     period   := step | hourly | daily | monthly | yearly
 *     kind     := instant | mean        (parsed but currently ignored —
 *                                        cadence-default kinds apply)
 *     WSEP     := [ \t]+                ; one or more spaces/tabs
 */
#ifndef FESOM_IO_CONFIG_H
#define FESOM_IO_CONFIG_H

#include "fesom_calendar.h"

#include <stddef.h>

typedef struct fesom_io_config_entry {
    char                *varname;          /* heap, one per LHS varname */
    fesom_period_kind_t *cadences;         /* heap, n_cadences entries  */
    int                  n_cadences;
} fesom_io_config_entry_t;

typedef struct fesom_io_config {
    fesom_io_config_entry_t *entries;
    int                      n_entries;
} fesom_io_config_t;

/* Parse a config file. Returns 0 on success, non-zero on syntax error.
 * On error, prints `<path>:<lineno>: <reason>` to stderr and the partial
 * `*out` is freed for the caller. On success, fills `*out` (caller must
 * call fesom_io_config_free). */
int fesom_io_config_parse(const char *path, fesom_io_config_t *out);

void fesom_io_config_free(fesom_io_config_t *cfg);

/* Return the entry for the named variable, or NULL if not present.
 * Linear scan — config files are tiny (≤ a few dozen entries). */
const fesom_io_config_entry_t *
fesom_io_config_lookup(const fesom_io_config_t *cfg, const char *varname);

#endif /* FESOM_IO_CONFIG_H */
