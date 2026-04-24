#ifndef FESOM_PHC_H
#define FESOM_PHC_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_tracers;
struct fesom_partit;

/*
 * Read PHC3.0 (or any FESOM-compatible 3D-IC NetCDF) and apply to the
 * tracer fields. Literal port of:
 *   nc_readGrid    (gen_ic3d.F90:68-217)
 *   nc_ic3d_ini    (gen_ic3d.F90:220-301)
 *   getcoeffld     (gen_ic3d.F90:303-491)
 *   do_ic3d        (gen_ic3d.F90:493-...) — top-level driver
 *   binarysearch   (gen_ic3d.F90:666-703)
 *   extrap_nod3D   (gen_support.F90:400-507) — horizontal then vertical fill
 *   insitu2pot     (oce_ale_pressure_bv.F90:3074-3118) — calls ptheta (Bryden 1973)
 *
 * Behaviour matches Fortran with `ic_cyclic = .true.` (default), `t_insitu`
 * controlled via the function arg. Loads T (var "temp") and S (var "salt")
 * from the same file. Land marker = `dummy = 1e10` (matches Fortran g_config).
 *
 * Note: Fortran nc_readGrid (lines 213-216) applies a cyclic wrap to nc_lon
 * (lon[0] -= 360, lon[Nlon-1] += 360). Together with the ncdata wrap
 * (getcoeffld lines 363-371: read [1..Nlon-2], copy [Nlon-2]→[0] and
 * [1]→[Nlon-1]), this implements periodic boundary at the dateline.
 * For PHC (lon = 0.5..359.5, no ghost cells in the file), the file's
 * data at lon[0]=0.5° and lon[Nlon-1]=359.5° is silently REPLACED by the
 * wrap-ghost values (358.5° and 1.5° data respectively). This is a
 * Fortran-side quirk we mirror exactly — produces a 2°-wide band of
 * smeared data near 0°/360° longitude that is invisible at ocean scale.
 */
void fesom_phc_load_ic(const char                  *path,
                       const struct fesom_mesh     *mesh,
                       struct fesom_tracers        *tracers,
                       int                          t_insitu,
                       struct fesom_partit         *partit);

#endif /* FESOM_PHC_H */
