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

/*
 * Ice-mediated surface momentum flux. Mirror of Fortran oce_fluxes_mom
 * (ice_oce_coupling.F90:53). Must be called AFTER fesom_ice_step (so
 * ice->uice/vice halo is populated by the closing EVP exchange) and BEFORE
 * the ocean step's compute_vel_rhs / impl_vert_visc (which consume
 * forcing->stress_surf).
 *
 * The atmosphere-ocean stress is taken in place from forcing->stress_node_surf
 * as written by fesom_bulk_compute. This is faithful in our (non-OASIS,
 * non-transit) config: no other code reads the atmoce field separately.
 *
 * Writes:
 *   ice->stress_iceoce_x/y[n]               for n in [0, myDim+eDim)
 *   forcing->stress_node_surf[2n+0,1]       for n in [0, myDim+eDim)  (in place)
 *   forcing->stress_surf[2*elem+0,1]        for elem in [0, myDim_elem2D)
 *
 * Halo coverage matches Fortran: per-node loop covers myDim+eDim so no
 * exchange is needed before the element loop reads stress_node_surf at
 * elnodes. The element loop matches Fortran's myDim-only bound (the only
 * downstream reader, impl_vert_visc, also iterates myDim).
 */
void fesom_ice_oce_fluxes_mom(fesom_ice              *ice,
                              struct fesom_partit    *partit,
                              struct fesom_mesh      *mesh,
                              struct fesom_forcing   *forcing);

#endif /* FESOM_ICE_COUPLING_H */
