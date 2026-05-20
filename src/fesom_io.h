#ifndef FESOM_IO_H
#define FESOM_IO_H

#include "fesom_calendar.h"
#include "fesom_io_gather.h"
#include "fesom_io_stream.h"
#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;
struct fesom_tracers;
struct fesom_aux;
struct fesom_partit;
struct fesom_ice;
struct fesom_forcing;

/* ------------------------------------------------------------------ */
/* fesom_state — POD bundle of pointers into the model's per-step    */
/* state. Built once per timestep by the orchestrator and passed by  */
/* reference to each stream's resolver. No data is copied; lifetimes */
/* are tied to the live model state for the duration of one          */
/* fesom_io_step call.                                                */
/* ------------------------------------------------------------------ */
typedef struct fesom_state {
    const struct fesom_mesh    *mesh;
    const struct fesom_dyn     *dyn;
    const struct fesom_tracers *tracers;
    const struct fesom_aux     *aux;
    const struct fesom_ice     *ice;       /* may be NULL (no sea ice) */
    const struct fesom_forcing *forcing;
} fesom_state;

/* gather_plan API lives in fesom_io_gather.h (included above).      */

/* ------------------------------------------------------------------ */
/* Orchestrator. Owns the model-internal calendar and one stream per */
/* cadence. Per-step, the main loop calls only fesom_io_step().      */
/* ------------------------------------------------------------------ */
typedef struct fesom_io {
    fesom_calendar_t       calendar;
    fesom_calendar_t       prev_calendar;
    fesom_io_stream_t      streams[5];        /* indexed by fesom_period_kind_t */
    int                    stream_active[5];  /* 1 = stream init'd, 0 = inactive */
    /* Per-cadence var arrays, heap-owned and freed at finalize. Streams
     * hold a const pointer into one of these. NULL for inactive cadences. */
    struct fesom_var_desc *owned_vars[5];
    int                    owned_nvars[5];
    char                  *out_dir;           /* heap-owned copy */
} fesom_io_t;

/* Initialise the orchestrator. Sets up the calendar at (year, month, day),
 * and configures the monthly stream with the compiled-in default variable
 * table. Other four cadences (step/hourly/daily/yearly) are inactive at
 * Task 6; Task 7 wires them on. `config_file` (override-file path) is
 * accepted but parsed only in Task 8 — pass NULL for now. */
void fesom_io_init(fesom_io_t                  *io,
                   const char                  *out_dir,
                   fesom_calendar_kind_t        cal_kind,
                   int                          year,
                   int                          month,
                   int                          day,
                   const char                  *config_file_or_null,
                   const struct fesom_mesh     *mesh,
                   struct fesom_partit         *partit);

/* Per-step entry. Saves prev_calendar = calendar, advances calendar by dt,
 * dispatches to each active stream's step. The only thing the main loop
 * has to call to drive output. */
void fesom_io_step(fesom_io_t              *io,
                   double                   dt,
                   const fesom_state       *state,
                   struct fesom_partit     *partit);

/* Mid-window flush + cleanup. Call at the end of the run. */
void fesom_io_finalize(fesom_io_t          *io,
                       struct fesom_partit *partit);

/*
 * Write one snapshot of the model state to a NetCDF file.
 *
 * Layout (chosen for direct readability with xarray / ncview / Paraview):
 *
 *   dimensions:
 *     nod2 = mesh->nod2D
 *     elem = mesh->elem2D
 *     nz   = mesh->nl              (interfaces — for w, zbar)
 *     nz_1 = mesh->nl - 1          (mid-layers — for T, S, u, v, Z)
 *     three = 3                    (vertices per element)
 *
 *   coordinates (1D, time-independent):
 *     lon (nod2)        geographic longitude, degrees
 *     lat (nod2)        geographic latitude,  degrees
 *     zbar(nz)          interface depths,     m  (negative downward)
 *     Z   (nz_1)        mid-layer depths,     m
 *
 *   topology:
 *     elem_nodes (elem, three)     node IDs (0-based)
 *     nlevels_nod2D (nod2)         per-node level count
 *     nlevels       (elem)         per-element level count
 *
 *   state:
 *     T (nz_1, nod2)               potential temperature, °C
 *     S (nz_1, nod2)               salinity, PSU
 *     eta_n (nod2)                 SSH, m
 *     w (nz, nod2)                 vertical velocity at interfaces, m/s
 *     u (nz_1, elem)               zonal velocity at element mid-layers, m/s
 *     v (nz_1, elem)               meridional velocity, m/s
 *     density_m_rho0 (nz_1, nod2)  in-situ density minus 1030, kg/m³
 *     bvfreq (nz_1, nod2)          N², s⁻²
 *     pgf_x, pgf_y (nz_1, elem)    pressure-gradient force, m/s²
 *     Kv (nz_1, nod2), Av (nz_1, elem)   vertical diffusivity / viscosity
 *
 *   sea-ice state (only when ice != NULL; all per-node, single level):
 *     a_ice (nod2)                 ice area fraction
 *     m_ice (nod2)                 ice volume per area, m
 *     m_snow(nod2)                 snow volume per area, m
 *     uice (nod2), vice (nod2)     ice velocity, m/s
 *     h_ice (nod2), h_snow (nod2)  diagnostic ice/snow thickness, m
 *
 *   time dimension:
 *     time = 1 (UNLIMITED, model_time_seconds for this snapshot).
 *     All time-varying state variables get `time` as the outermost axis so
 *     `xarray.open_mfdataset("snap_*.nc", concat_dim="time")` produces a
 *     single dataset spanning all snapshot files.
 *
 *   global attributes:
 *     step (int)                   timestep counter
 *     dt (double)                  timestep size, s
 *
 * Multi-rank: when `partit->npes > 1`, every rank participates in MPI gathers
 * that ship local-myDim slices to rank 0; only rank 0 opens the file. Single-
 * rank or NULL partit: rank 0 writes its own data directly.
 */
void fesom_io_write_snapshot(const char                  *path,
                             int                          step_n,
                             real_t                       dt,
                             const fesom_calendar_t      *cal,
                             const struct fesom_mesh     *mesh,
                             const struct fesom_dyn      *dyn,
                             const struct fesom_tracers  *tracers,
                             const struct fesom_aux      *aux,
                             const struct fesom_ice      *ice,
                             struct fesom_partit         *partit);

#endif /* FESOM_IO_H */
