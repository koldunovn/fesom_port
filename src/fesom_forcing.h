#ifndef FESOM_FORCING_H
#define FESOM_FORCING_H

#include "fesom_types.h"

struct fesom_mesh;

/*
 * Surface forcing (FRESH_START.md §9). Phase 1 starts with zero forcing;
 * later steps wire JRA55-do bulk formulae and SSS restoring.
 *
 *   heat_flux         [nod2D]      net surface heat flux (W/m², +ve = ocean loses)
 *   water_flux        [nod2D]      E - P virtual freshwater (m/s, +ve = ocean loses)
 *   stress_node_surf  [nod2D * 2]  wind stress at nodes (N/m²)
 *   stress_surf       [elem2D * 2] wind stress at elements (N/m²)
 *
 * stress_node_surf and stress_surf are kept in sync by interpolation in the
 * Fortran code. We mirror both because impl_vert_visc reads stress at nodes
 * while compute_vel_rhs uses stress at elements.
 */
typedef struct fesom_forcing {
    real_t *heat_flux;
    real_t *water_flux;
    real_t *stress_node_surf;
    real_t *stress_surf;

    /* Phase 3 step 25 — SSS restoring + CORE2 runoff (o_ARRAYS in Fortran). */
    real_t *runoff;          /* [nod2D]  freshwater input from rivers, m/s   */
    real_t *Ssurf;           /* [nod2D]  monthly SSS climatology, PSU       */
    real_t *virtual_salt;    /* [nod2D]  rsss·water_flux (linfs only)        */
    real_t *relax_salt;      /* [nod2D]  surf_relax_S·(Ssurf - S_top)        */

    /* Per-node bulk transfer coefficients for the open-water (atm-ocean) path.
       Computed by fesom_bulk_compute via L&Y09; consumed by sea-ice
       thermodynamics (fesom_ice_thermodynamics). Mirrors Fortran
       Ch_atm_oce_arr / Ce_atm_oce_arr in g_forcing_arrays. */
    real_t *Ch_atm_oce;      /* [nod2D]  sensible-heat transfer coeff (open water) */
    real_t *Ce_atm_oce;      /* [nod2D]  evaporation transfer coeff   (open water) */

    /* Shortwave penetration (use_sw_pene). Mirrors Fortran g_forcing_arrays
       chl / sw_3d. chl: chlorophyll [mg/m³] (const or Sweeney monthly climatology).
       sw_3d: visible-band shortwave TEMPERATURE flux [K m/s] at level interfaces,
       built by fesom_cal_shortwave_rad (oce_shortwave_pene.F90), consumed by the
       temperature tracer equation (oce_ale_tracer.F90:990). */
    real_t *chl;             /* [nod2D]        chlorophyll concentration            */
    real_t *sw_3d;           /* [nod2D * nl]   shortwave T-flux at interfaces        */
} fesom_forcing;

void fesom_forcing_alloc(fesom_forcing *f, const struct fesom_mesh *mesh);
void fesom_forcing_free (fesom_forcing *f);

#endif /* FESOM_FORCING_H */
