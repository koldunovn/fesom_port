#ifndef FESOM_ICE_COUPLING_H
#define FESOM_ICE_COUPLING_H

#include "fesom_ice_types.h"

struct fesom_mesh;
struct fesom_partit;
struct fesom_dyn;
struct fesom_tracers;
struct fesom_forcing;
struct fesom_sss_runoff;

/*
 * Ocean → ice direction. Mirror of Fortran ocean2ice (ice_oce_coupling.F90:154).
 *
 * Populates ice->srfoce_temp/salt/u/v/ssh from the surface-layer ocean state.
 * Cavity skip via ulevels_nod2D > 1. The u_w/v_w computation is an area-weighted
 * mean of element-center UV at the surface, then halo-exchanged.
 *
 * `ice_update`-branch path (when ice steps less often than ocean) is mirrored
 * but inactive in our case (ice_ave_steps = 1).
 */
void fesom_ocean2ice(fesom_ice                 *ice,
                     const struct fesom_dyn    *dyn,
                     const struct fesom_tracers *tracers,
                     struct fesom_partit        *partit,
                     struct fesom_mesh          *mesh);

/*
 * Ice → ocean heat / freshwater fluxes. Mirror of Fortran oce_fluxes
 * (ice_oce_coupling.F90:249) — non-ICEPACK branch only (lines 393-398):
 *
 *   heat_flux  = -flx_h        (ice%flx_h is +down; ocean uses +up)
 *   water_flux = -flx_fw       (no -runoff term; runoff is in flx_fw via therm_ice)
 *
 * Then virtual_salt = rsss · water_flux (use_virt_salt branch) and
 * relax_salt = surf_relax_S · (Ssurf - S_top), both balanced to net zero.
 *
 * sr (sss_runoff state) provides ref_sss / ref_sss_local / use_virt_salt /
 * surf_relax_S. After this lands, the runoff subtraction inside
 * fesom_sss_runoff_step must be removed (Task C3b) to avoid double-counting.
 */
void fesom_ice_oce_fluxes(fesom_ice                     *ice,
                          struct fesom_partit           *partit,
                          struct fesom_mesh             *mesh,
                          const struct fesom_tracers    *tracers,
                          struct fesom_forcing          *forcing,
                          const struct fesom_sss_runoff *sr);

#endif /* FESOM_ICE_COUPLING_H */
