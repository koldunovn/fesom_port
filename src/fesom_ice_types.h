#ifndef FESOM_ICE_TYPES_H
#define FESOM_ICE_TYPES_H

#include "fesom_types.h"

/*
 * Mirror of Fortran T_ICE / T_ICE_DATA / T_ICE_WORK / T_ICE_THERMO from
 * MOD_ICE.F90. Scope cut for the C port (per docs/plans/20260425-sea-ice-port.md):
 *
 *   - whichEVP = 0 only (standard EVP). uice_aux/vice_aux/alpha_evp_array/
 *     beta_evp_array NOT mirrored — Fortran allocates them only when whichEVP != 0.
 *   - non-ICEPACK build. No __oasis/__oifs/__yac/__ifsinterface fields.
 *   - num_itracers = 3 (a_ice, m_ice, m_snow). No ice_temp tracer (oifs only).
 *   - No icebergs (uice_ib/vice_ib).
 *   - No melt-pond fields beyond apnd/hpnd/ipnd defaults (use_meltponds = false).
 *   - No CMIP6 dyngr/dyngrsn/dyngra diagnostics (output-only, no physics impact).
 *   - No coupled-atmosphere atmcoupl substruct.
 */

/* Tracer indices into ice->data[]. Mirrors Fortran ice%data(1..3). */
#define FESOM_NUM_ICE_TRACERS 3
#define FESOM_ICE_AICE  0     /* ice area fraction      a_ice  */
#define FESOM_ICE_MICE  1     /* ice volume per area    m_ice  */
#define FESOM_ICE_MSNOW 2     /* snow volume per area   m_snow */

/* Mirror of T_ICE_DATA (MOD_ICE.F90:14-23). All sized [myDim+eDim] (nodes). */
typedef struct fesom_ice_data {
    real_t *values;          /* current value          */
    real_t *values_old;      /* previous timestep      */
    real_t *values_rhs;      /* high-order FCT RHS     */
    real_t *values_div_rhs;  /* divergence-form RHS    */
    real_t *dvalues;         /* increment / scratch    */
    real_t *valuesl;         /* low-order FCT solution */
    int     id;              /* 0..2, set at init      */
} fesom_ice_data;

/*
 * Mirror of T_ICE_WORK (MOD_ICE.F90:28-41). Mix of node and element arrays.
 *
 *   fct_tmax/tmin/plus/minus  [myDim+eDim]                  (nodes)
 *   fct_fluxes                [elem_size * 3]               (3 edge fluxes per element)
 *   fct_massmatrix            [sum(nn_num(1:myDim_nod2D))]  (CSR-style; sized at alloc)
 *   sigma11/12/22, eps11/12/22, ice_strength
 *                             [elem_size]                    (elements)
 *   inv_areamass, inv_mass    [myDim+eDim]                  (nodes)
 */
typedef struct fesom_ice_work {
    real_t *fct_tmax;
    real_t *fct_tmin;
    real_t *fct_plus;
    real_t *fct_minus;
    real_t *fct_fluxes;
    real_t *fct_massmatrix;
    real_t *sigma11, *sigma12, *sigma22;
    real_t *eps11,   *eps12,   *eps22;
    real_t *ice_strength;
    real_t *inv_areamass;
    real_t *inv_mass;
} fesom_ice_work;

/*
 * Mirror of T_ICE_THERMO (MOD_ICE.F90:46-101). Per-node arrays + scalar
 * physical constants and namelist parameters. Defaults match the Fortran
 * type defaults (T_ICE_THERMO lines 53-97), which match work_core/namelist.ice.
 *
 * Skipped:
 *   - dyngr/dyngrsn/dyngra (CMIP6 output-only, no physics impact — see Out-of-scope)
 *   - hpdf[15], h_cutoff, new_iclasses (active only when use_meltponds or new_iclasses=true)
 */
typedef struct fesom_ice_thermo {
    /* per-node arrays */
    real_t *t_skin;
    real_t *thdgr, *thdgrsn, *thdgra;
    real_t *thdgr_old;
    real_t *ustar;
    real_t *apnd, *hpnd, *ipnd;          /* meltpond fields, zeroed; use_meltponds=false */

    /* density and inverse (Fortran T_ICE_THERMO:53-57) */
    real_t rhoair,  inv_rhoair;
    real_t rhowat,  inv_rhowat;
    real_t rhofwt,  inv_rhofwt;
    real_t rhoice,  inv_rhoice;
    real_t rhosno,  inv_rhosno;

    /* specific heats J/(kg*K) */
    real_t cpair, cpice, cpsno;

    /* derived (computed in init): cc = rhowat*4190, cl = rhoice*3.34e5 */
    real_t cc, cl;

    /* latent heat (J/kg) */
    real_t clhw;          /* water -> water vapor */
    real_t clhi;          /* sea ice -> water vapor */

    real_t tmelt;         /* 273.15 K */
    real_t boltzmann;     /* sigma * emissivity */

    int    iclasses;      /* number of ice thickness gradations (default 7) */
    real_t hmin;          /* cutoff ice thickness, m */
    real_t armin;         /* minimum ice concentration */

    /* namelist /ice_therm/ */
    real_t con, consn;    /* thermal conductivities ice, snow (W/m/K) */
    real_t Sice;          /* ice salinity (ppt) */
    real_t h0, h0_s;      /* lead closing N hemi / S hemi (m) */
    real_t emiss_ice, emiss_wat;
    real_t albsn, albsnm; /* snow albedo: frozen, melting */
    real_t albi,  albim;  /* ice  albedo: frozen, melting */
    real_t albw;          /* open water albedo (LY2004) */
    real_t h_ml;          /* upper-layer thickness for heat available */

    /* additional namelist params */
    int    snowdist;
    int    open_water_albedo;   /* 0=standard 1=taylor 2=briegleb */
    real_t c_melt;              /* concentration eq. coefficient on melt */

    int    use_meltponds;       /* 0 in CORE2 */
} fesom_ice_thermo;

/*
 * Mirror of T_ICE (MOD_ICE.F90:128-222) — the top-level sea-ice state struct.
 *
 * Sizing convention (per fesom_dyn.c):
 *   node arrays:    myDim_nod2D  + eDim_nod2D
 *   element arrays: myDim_elem2D + eDim_elem2D + eXDim_elem2D
 */
typedef struct fesom_ice {
    /* tracers — a_ice, m_ice, m_snow (indices FESOM_ICE_AICE/MICE/MSNOW) */
    fesom_ice_data data[FESOM_NUM_ICE_TRACERS];

    /* zonal & meridional ice velocity (nodes) */
    real_t *uice, *uice_rhs, *uice_old;
    real_t *vice, *vice_rhs, *vice_old;

    /* surface stresses atm<->ice and oce<->ice (nodes) */
    real_t *stress_atmice_x, *stress_iceoce_x;
    real_t *stress_atmice_y, *stress_iceoce_y;

    /* surface ocean state seen by ice (populated by ocean2ice; nodes) */
    real_t *srfoce_temp, *srfoce_salt, *srfoce_ssh;
    real_t *srfoce_u,    *srfoce_v;

    /* fluxes from thermodynamics (nodes)
       flx_fw is fresh_wa_flux (m/s, includes runoff per non-ICEPACK contract)
       flx_h  is net_heat_flux (W/m^2)                                          */
    real_t *flx_fw;
    real_t *flx_h;

    /* diagnostic ice/snow thicknesses h = m / max(a, 1e-3) (nodes) */
    real_t *h_ice, *h_snow;

    /* embedded substructs */
    fesom_ice_work   work;
    fesom_ice_thermo thermo;

    /* --- RHEOLOGY scalar parameters (defaults from T_ICE:187-202 = namelist.ice) --- */
    real_t pstar;            /* 30000.0 N/m^2 */
    real_t ellipse;          /* 2.0           */
    real_t c_pressure;       /* 20.0          */
    real_t delta_min;        /* 1.0e-11 s^-1  */
    real_t Clim_evp;         /* 615 kg/m^2    */
    real_t zeta_min;         /* 4.0e+8 kg/s   */

    int    evp_rheol_steps;  /* 120  EVP subcycles per ice timestep */
    real_t ice_gamma_fct;    /* 0.5 smoothing in ice FCT advection (CORE2 namelist) */
    real_t ice_diff;         /* 10.0 stabilising diffusion          */
    real_t theta_io;         /* 0.0  ice-ocean rotation angle       */
    real_t cd_oce_ice;       /* 5.5e-3 ocean-ice drag coefficient   */
    int    ice_free_slip;    /* 0     */

    /* --- whichEVP locked at 0; aux/alpha_evp/beta_evp arrays not allocated --- */
    int    whichEVP;         /* MUST be 0; dispatcher aborts otherwise */

    /* --- timestep --- */
    int    ice_ave_steps;    /* 1 — ice timestep = ice_ave_steps * ocean step */
    real_t ice_dt;           /* set in fesom_ice_setup from ocean dt          */
    real_t Tevp_inv;         /* set in fesom_ice_setup = evp_rheol_steps/ice_dt */

    /* misc state */
    int    ice_steps_since_upd;  /* 0 */
    int    ice_update;           /* 1 */
} fesom_ice;

#endif /* FESOM_ICE_TYPES_H */
