#ifndef FESOM_IO_H
#define FESOM_IO_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;
struct fesom_tracers;
struct fesom_aux;
struct fesom_partit;

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
                             const struct fesom_mesh     *mesh,
                             const struct fesom_dyn      *dyn,
                             const struct fesom_tracers  *tracers,
                             const struct fesom_aux      *aux,
                             struct fesom_partit         *partit);

#endif /* FESOM_IO_H */
