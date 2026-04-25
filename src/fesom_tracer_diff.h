#ifndef FESOM_TRACER_DIFF_H
#define FESOM_TRACER_DIFF_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_aux;
struct fesom_tracers;
struct fesom_dyn;
struct fesom_forcing;
struct fesom_gm;

/*
 * Implicit vertical diffusion of tracers + surface heat/water flux BC.
 * Literal port of diff_ver_part_impl_ale (oce_ale_tracer.F90:562-1082)
 * for the FCT path (do_wimpl=false), no KPP nonlocal.
 *
 * Phase G5 plumbed the Redi K33 augmentation: when `gm != NULL`,
 * Ki and slope_tapered are read to extend the diagonal vertical
 * diffusivity (Fortran `Ty / Ty1` projection of horizontal-Redi onto
 * the vertical axis). With gm=NULL the routine reduces to its pre-G5
 * behaviour exactly (Ki and st² terms drop to zero).
 *
 * For each tracer (T, S):
 *   1. Per-node TDMA on the column with K_v from PP mixing (+ Redi K33
 *      if gm is non-NULL).
 *   2. Surface BC via bc_surface (oce_ale_tracer.F90:1517-1524):
 *        T:   bc_surface = -dt*(heat_flux/vcpw)        (linfs → is_nonlinfs=0)
 *        S:   bc_surface =  dt*(virtual_salt + relax_salt)   (zero in step 24)
 *   3. tracer += dT (the TDMA solution).
 *
 * Without this stage, surface heat_flux/water_flux are never applied — T/S
 * only respond to advection. Phase 3 step 24 needs this to make JRA55 visible.
 */
void fesom_impl_vert_diff_tracers(const struct fesom_mesh    *mesh,
                                  const struct fesom_aux     *aux,
                                  const struct fesom_forcing *forcing,
                                  struct fesom_tracers       *tracers,
                                  const struct fesom_gm      *gm);

#endif /* FESOM_TRACER_DIFF_H */
