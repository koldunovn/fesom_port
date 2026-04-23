#ifndef FESOM_FORCING_ANALYTICAL_H
#define FESOM_FORCING_ANALYTICAL_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_forcing;

/*
 * Phase 2 step 17 — analytical wind + (optional) SST restoring.
 *
 * The Fortran reads JRA55 NetCDF for real forcing (gen_surface_forcing.F90).
 * For Phase 2 we don't need NetCDF yet — we synthesise a known-good wind
 * pattern and (optionally) a latitudinally varying climatology to relax T to.
 *
 * fesom_forcing_set_analytical(mesh, forcing, tau0, Ly_factor):
 *   stress_surf[e][0] = -tau0 * cos(2π * lat_centroid / (Ly_factor * π/2))
 *   stress_surf[e][1] = 0
 *   heat_flux  [n]    = 0   (Slice A: no thermal forcing)
 *   water_flux [n]    = 0
 *
 * For the simplest first test, pass tau0 = 0.05 N/m² and Ly_factor = 2.0
 * so the wind switches sign every 90° latitude — a single-cell zonal
 * pattern, easterly trades + mid-latitude westerlies. This drives an
 * Ekman drift in the surface and gyre-scale circulation underneath.
 *
 * Centroid latitude is taken from mesh->geo_coord_nod2D averaged over
 * the 3 vertices (geographic radians, NOT rotated).
 */
void fesom_forcing_set_analytical(const struct fesom_mesh *mesh,
                                  struct fesom_forcing    *forcing,
                                  real_t tau0,
                                  real_t Ly_factor);

#endif /* FESOM_FORCING_ANALYTICAL_H */
