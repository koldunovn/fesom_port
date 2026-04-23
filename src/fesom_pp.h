#ifndef FESOM_PP_H
#define FESOM_PP_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_aux;
struct fesom_dyn;

/*
 * Pacanowski-Philander vertical mixing + convective adjustment.
 * Literal ports of:
 *   compute_vel_nodes  (oce_dyn.F90:177-226)
 *   oce_mixing_pp      (oce_ale_mixing_pp.F90:2-93)
 *   mo_convect         (oce_mo_conv.F90:78-120, use_instabmix branch only)
 *
 * Phase 1 deviations from Fortran defaults:
 *   - use_momix = false   (Fortran default true; momix needs ice + forcing)
 *   - use_windmix = false (matches default)
 *   - Kv0_const = true    (matches default; we don't port the lat/depth Kv0)
 *
 * compute_vel_nodes interpolates the cell-centred UV onto nodes via an
 * area-weighted mean over surrounding elements. PP uses the resulting
 * UVnode to compute vertical shear at scalar locations.
 */
void fesom_compute_vel_nodes(const struct fesom_mesh *mesh,
                             struct fesom_dyn        *dyn);

void fesom_pp_mixing(const struct fesom_mesh *mesh,
                     const struct fesom_dyn  *dyn,
                     struct fesom_aux        *aux);

/* Convective adjustment: where N² < 0, raise Kv (and Av at adjacent cells)
   to instabmix_kv. Mirror of the use_instabmix branch in mo_convect. */
void fesom_mo_convect(const struct fesom_mesh *mesh,
                      struct fesom_aux        *aux);

#endif /* FESOM_PP_H */
