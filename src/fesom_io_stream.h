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
 * fesom_io_stream.h - one stream owns N variables at the same cadence
 * (period). On each fesom_io_step it accumulates each variable's local
 * (myDim) slice; on calendar boundary it gathers via the public
 * gather_plan API and writes one record to a per-(variable, year, cadence)
 * netCDF file. FESOM-style file layout: temp.fesom.1948.monthly.nc, etc.
 */
#ifndef FESOM_IO_STREAM_H
#define FESOM_IO_STREAM_H

#include "fesom_calendar.h"
#include "fesom_io_gather.h"
#include "fesom_types.h"

#include <mpi.h>

struct fesom_mesh;
struct fesom_partit;
struct fesom_state;          /* full def in fesom_io.h */

/* ------------------------------------------------------------------ */
/* Variable shape — drives buffer sizing, gather routing, and netCDF  */
/* dimension binding. The stride invariant for FESOM_VAR_3D_NODE_MID  */
/* is `nl`, NOT `nl-1`, because the underlying tracer storage uses    */
/* stride nl (see feedback_tracer_stride_nl.md — cost a session of    */
/* debug time previously). The FILE only writes nl-1 mid-layer values */
/* per column; the unused last layer is simply not read at write.    */
/* ------------------------------------------------------------------ */
typedef enum {
    FESOM_VAR_2D_NODE       = 0,   /* [local_n]                              */
    FESOM_VAR_3D_NODE_MID   = 1,   /* [local_n × nl]; nl-1 layers written    */
    FESOM_VAR_3D_NODE_IFACE = 2,   /* [local_n × nl]; nl interfaces written  */
    FESOM_VAR_2D_ELEM       = 3,   /* [local_e]                              */
    FESOM_VAR_3D_ELEM_MID   = 4,   /* [local_e × nl]; nl-1 layers written    */
    FESOM_VAR_2D_ELEM_VEC   = 5,   /* [local_e × 2]                          */
    FESOM_VAR_3D_ELEM_IFACE = 6,   /* [local_e × nl]; nl interfaces written  */
} fesom_var_kind_t;

typedef enum {
    FESOM_OUT_INSTANT = 0,         /* sample at flush time                   */
    FESOM_OUT_MEAN    = 1,         /* sum / n_accum at flush                 */
} fesom_output_kind_t;

/* Variable resolver: fills `n` real_t elements at `out` from the model
 * state, ADDING to whatever's already there. For MEAN this means
 * `accum[i] += value[i]` per call; for INSTANT (Task 7) the caller
 * zeros the buffer first, so the same function works as a pure write.
 *
 * Why "add" rather than "return pointer": some logical variables don't
 * exist as contiguous buffers in the model state — e.g. SST is the
 * surface slice of the 3D T tracer (strided at nl), and `u` / `v` are
 * components of the [E × nl × 2] uv array. The resolver does the
 * extraction loop. For 1:1 fields it's just `for i: out[i] += src[i]`. */
typedef void (*fesom_resolve_fn)(const struct fesom_state *state, real_t *out, size_t n);

/* Variable descriptor. Lifetime = run; held by reference, never copied. */
typedef struct fesom_var_desc {
    const char         *name;       /* file basename: "temp", "ssh", "u", ... */
    const char         *long_name;  /* CF "long_name" attribute               */
    const char         *units;      /* CF "units" attribute                   */
    fesom_var_kind_t    kind;
    fesom_resolve_fn    resolve;
} fesom_var_desc_t;

/* Compute the size (in real_t elements) of the local accumulator buffer
 * for a given variable kind. Exposed so test_io_stream_unit can probe it
 * without instantiating a stream. */
size_t fesom_io_stream_accum_size(fesom_var_kind_t kind, int local_n, int local_e, int nl);

/* Compute the global-side write count for one record (in real_t elements)
 * and the netCDF file's per-variable spatial dimension extent. */
size_t fesom_io_stream_global_size(fesom_var_kind_t kind, int nod2D, int elem2D, int nl);
int    fesom_io_stream_n_layers_written(fesom_var_kind_t kind, int nl);

/* ------------------------------------------------------------------ */
/* Stream state. One per cadence. All fields private to the           */
/* implementation; declared here only so the orchestrator can store   */
/* streams by value inside its own struct.                            */
/* ------------------------------------------------------------------ */
typedef struct fesom_io_stream {
    fesom_period_kind_t      period;
    fesom_output_kind_t      output_kind;
    const fesom_var_desc_t  *vars;
    int                      nvars;

    /* Accumulator buffers, one per variable.
     * For MEAN: holds running sum; flush divides by n_accum.
     * For INSTANT: NULL (per-step path samples + writes immediately). */
    real_t                 **accum;
    size_t                  *accum_sz;          /* [nvars] elements per accum */
    int                      n_accum;           /* # of step() calls in window */

    /* netCDF handles per variable. ncid == -1 means file not open yet
     * (lazy open on first flush of the rollover window). */
    int                     *ncid;
    int                     *var_id;
    int                     *time_id;
    int                     *bnds_id;
    int                     *time_index;        /* next record index */
    int                     *dim_time_id;       /* per-file time dim ID, cached */

    /* Rollover-window key. -1 means no file open for any window.
     * Active window depends on period:
     *   STEP    -> per-day:    open_year, open_month, open_day
     *   HOURLY  -> per-month:  open_year, open_month (open_day = 0)
     *   DAILY/MONTHLY/YEARLY -> per-year: open_year (open_month/day = 0) */
    int                      open_year;
    int                      open_month;
    int                      open_day;

    /* Cached MPI plumbing — built once at init, freed at close.
     * Per-step gather happens against this plan; never rebuilt. */
    gather_plan              gp;
    MPI_Comm                 comm;
    int                      mype;

    /* Static metadata kept for write-time. */
    const struct fesom_mesh *mesh;
    char                    *out_dir;            /* owned copy */

    /* Cached global node + element coords (geographic degrees), built
     * once at init via gather. NULL on non-zero ranks (rank 0 is the
     * only writer). Element coords are the centroid (mean of 3 vertex
     * coords) — computed locally on each rank then gathered. */
    double                  *cached_lon_node;    /* [nod2D]  on rank 0 */
    double                  *cached_lat_node;
    double                  *cached_lon_elem;    /* [elem2D] on rank 0 */
    double                  *cached_lat_elem;
} fesom_io_stream_t;

/* ------------------------------------------------------------------ */
/* Lifecycle.                                                         */
/* ------------------------------------------------------------------ */

/* Initialise a stream. `vars` array lifetime = run (typically the
 * compiled-in default table). MEAN streams allocate per-variable accum
 * buffers; INSTANT streams skip that. The cached gather_plan is built
 * here (one MPI_Gather + two MPI_Gatherv). */
void fesom_io_stream_init(fesom_io_stream_t       *s,
                          fesom_period_kind_t      period,
                          fesom_output_kind_t      kind,
                          const fesom_var_desc_t  *vars,
                          int                      nvars,
                          const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          const char              *out_dir);

/* Per-step entry. Cheap: per variable, one local sum loop; no MPI
 * (for MEAN). Calls fesom_io_stream_flush internally when
 * fesom_calendar_crossed(prev, curr, period) returns non-zero. */
void fesom_io_stream_step(fesom_io_stream_t        *s,
                          const struct fesom_state *state,
                          const fesom_calendar_t   *prev,
                          const fesom_calendar_t   *curr,
                          struct fesom_partit      *partit);

/* Force a flush + close any open files. Called from fesom_io_finalize
 * to make sure mid-window data isn't lost when a run ends inside a
 * period. */
void fesom_io_stream_finalize(fesom_io_stream_t   *s,
                              const fesom_calendar_t *cal,
                              struct fesom_partit *partit);

/* Free all owned memory. Closes any still-open netCDF files. */
void fesom_io_stream_close(fesom_io_stream_t *s);

#endif /* FESOM_IO_STREAM_H */
