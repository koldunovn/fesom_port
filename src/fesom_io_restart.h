/*
 * fesom_io_restart.h — restart write and read.
 *
 * io_restart.F90 was a deliberate non-port. This is the replacement, and it is
 * NOT a literal port: see docs/plans/20260818-restart-io-port.md. Two ways it
 * departs from the Fortran, both on purpose:
 *
 *   1. It saves a STRICT SUPERSET of the Fortran field list. The Fortran
 *      restart omits the EVP stress tensor and the ice skin temperature, both
 *      of which survive a timestep in memory, so a Fortran restart is not the
 *      bitwise continuation of an uninterrupted run. Ours is.
 *   2. It carries the ALE free-surface stiffness matrix. That matrix is
 *      accumulated one increment per step (fesom_update_stiff_mat_ale), so it
 *      depends on the whole history, not on the current state. The Fortran
 *      rebuilds it after a restart by applying a single increment equal to the
 *      absolute elemental sea surface height, which telescopes to the right
 *      answer in exact arithmetic but not in floating point.
 *   3. One netCDF file per restart instead of one directory of one file per
 *      field, and no in-file time dimension.
 *
 * The file is partition-independent — written through the gather path, read
 * through the scatter path — so the rank count may change across a restart and
 * a file written by one backend can be read by another.
 *
 * Off unless asked for:
 *   FESOM_RESTART_OUT=<dir>    write restarts here
 *   FESOM_RESTART_EVERY=<n>    ... every n steps (0 = only at the last step)
 *   FESOM_RESTART_IN=<file>    start from this file instead of the PHC IC
 */
#ifndef FESOM_IO_RESTART_H
#define FESOM_IO_RESTART_H

#include "fesom_calendar.h"
#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;
struct fesom_tracers;
struct fesom_ice;
struct fesom_tke;
struct fesom_partit;
struct fesom_ssh_stiff;

#define FESOM_RESTART_FORMAT_VERSION 1

typedef struct fesom_restart_cfg {
    const char *out_dir;    /* NULL when writing is off  */
    const char *in_path;    /* NULL for a cold start     */
    int         every;      /* steps between restarts; 0 = last step only */
} fesom_restart_cfg;

/* Read the three environment variables. Never fails; an unset knob is off. */
void fesom_restart_config(fesom_restart_cfg *cfg);

/* True when step `step` of `nsteps` should produce a restart. */
int fesom_restart_due(const fesom_restart_cfg *cfg, int step, int nsteps);

/* The path a restart for this calendar date goes to:
 * <out_dir>/fesom.<year><mm><dd>_<sssss>.restart.nc */
void fesom_restart_path(char *buf, size_t buflen,
                        const char *out_dir, const fesom_calendar_t *cal);

/* Write every piece of persistent model state. Collective; rank 0 writes.
 * `dt` and `step` ride along as attributes for provenance and for the
 * continuity checks in the reader. */
void fesom_restart_write(const char                 *path,
                         int                         step,
                         double                      dt,
                         const fesom_calendar_t     *cal,
                         const struct fesom_mesh    *mesh,
                         const struct fesom_dyn     *dyn,
                         const struct fesom_tracers *tracers,
                         const struct fesom_ice     *ice,
                         const struct fesom_tke     *tke,
                         const struct fesom_ssh_stiff *stiff,
                         struct fesom_partit        *partit);

/* Read a restart into live model state, then re-derive what the file does not
 * carry (zbar_3d_n, Z_3d_n, helem, hnode_new, dhe), mirroring
 * restart_thickness_ale (oce_ale.F90:1449). Halo entries are exchanged here,
 * so the caller gets a state indistinguishable from the end of the step that
 * wrote it. `step` and `cal` are filled from the file.
 *
 * Aborts on a file that does not match this mesh (node/element/level counts),
 * on a missing field, and on a format version it does not know. A restart that
 * silently half-loads is worse than no restart. */
void fesom_restart_read(const char           *path,
                        int                  *step,
                        double               *dt,
                        fesom_calendar_t     *cal,
                        struct fesom_mesh    *mesh,
                        struct fesom_dyn     *dyn,
                        struct fesom_tracers *tracers,
                        struct fesom_ice     *ice,
                        struct fesom_tke     *tke,
                        struct fesom_ssh_stiff *stiff,
                        struct fesom_partit  *partit);

#endif /* FESOM_IO_RESTART_H */
