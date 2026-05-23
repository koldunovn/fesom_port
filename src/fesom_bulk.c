/*
 * NCAR L&Y09 bulk formulae + open-water flux assembly + wind-stress assembly.
 * No sea ice — Phase 3 ocean-only.
 */
#include "fesom_bulk.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_halo.h"
#include "fesom_ice_types.h"
#include "fesom_jra55.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <string.h>

/* MOD_ICE.F90 ice%thermo defaults (lines 53-79). */
#define BULK_RHOAIR      1.3
#define BULK_INV_RHOAIR  (1.0 / 1.3)
#define BULK_CPAIR       1005.0
#define BULK_CLHW        2.501e6     /* J/kg, water → vapor */
#define BULK_TMELT       273.15
#define BULK_BOLTZMANN   5.67e-8
#define BULK_EMISS_WAT   0.97
#define BULK_ALBW        0.1         /* open-water albedo — CORE2 namelist.ice albw=0.1
                                        (open_water_albedo=0, constant). NOT the LY2004
                                        0.066: the reference run overrides it to 0.1, and
                                        the C used 0.066 -> ~3.8% excess SW absorption ->
                                        tropics-max SST warm bias. Keep == fesom_ice.c th->albw. */
#define BULK_INV_RHOWAT  (1.0 / 1000.0)
#define BULK_GRAV        9.80
#define BULK_VONKARM     0.40
#define BULK_Q1          640380.0
#define BULK_Q2          (-5107.4)
#define BULK_U10MIN      0.3         /* m/s wind floor */

#define BULK_N_ITTS      5
#define BULK_INC_RATIO   1.0e-4

/* core_coeff_2z helper not needed — we use the L&Y09 path only. */

/* ------------------------------------------------------------------------ *
 * NCAR L&Y04/09 ocean exchange coefficients (cd, ce, ch).                  *
 * Literal port of ncar_ocean_fluxes_mode (gen_bulk_formulae.F90:126-341).  *
 * Returns by reference; only called once per node per step.                *
 * ------------------------------------------------------------------------ */
static void ncar_ocean_fluxes_mode(real_t tair_C,    /* °C */
                                   real_t shum,      /* kg/kg */
                                   real_t u_wind,    /* m/s */
                                   real_t v_wind,    /* m/s */
                                   real_t T_oc_C,    /* °C */
                                   real_t u_w,       /* m/s */
                                   real_t v_w,       /* m/s */
                                   real_t z_wind,
                                   real_t z_tair,
                                   real_t z_shum,
                                   real_t *cd_out,
                                   real_t *ce_out,
                                   real_t *ch_out)
{
    /* Convert to Kelvin (Fortran lines 184-185). */
    real_t t  = tair_C  + BULK_TMELT;
    real_t ts = T_oc_C  + BULK_TMELT;
    real_t q  = shum;
    real_t qs = 0.98 * BULK_Q1 * BULK_INV_RHOAIR * exp(BULK_Q2 / ts); /* L-Y eqn 5 */
    real_t tv = t * (1.0 + 0.608 * q);

    real_t dux = u_wind - u_w;
    real_t dvy = v_wind - v_w;
    real_t u   = sqrt(dux*dux + dvy*dvy);
    if (u < BULK_U10MIN) u = BULK_U10MIN;

    real_t u10 = u, t10 = t, q10 = q;

    /* Initial cd_n10 — LY2009 eqn 11a/b. */
    real_t hl1 = (2.7/u10 + 0.142 + 0.0764*u10 - 3.14807e-10 * pow(u10, 6)) / 1.0e3;
    real_t cd_n10 = (0.5 - copysign(0.5, u10 - 33.0)) * hl1
                  + (0.5 + copysign(0.5, u10 - 33.0)) * 2.34e-3;
    real_t cd_n10_rt = sqrt(cd_n10);
    real_t ce_n10    = 34.6 * cd_n10_rt * 1.0e-3;
    real_t stab      = 0.5 + copysign(0.5, t - ts);
    real_t ch_n10    = (18.0*stab + 32.7*(1.0 - stab)) * cd_n10_rt * 1.0e-3;

    real_t cd = cd_n10, ch = ch_n10, ce = ce_n10;
    real_t cd_prev = cd;

    for (int it = 0; it < BULK_N_ITTS; ++it) {
        real_t cd_rt = sqrt(cd);
        real_t ustar = cd_rt * u;                     /* L-Y eqn 7a */
        real_t tstar = (ch / cd_rt) * (t10 - ts);     /* L-Y eqn 7b */
        real_t qstar = (ce / cd_rt) * (q10 - qs);     /* L-Y eqn 7c */
        real_t bstar = BULK_GRAV * (tstar/tv + qstar / (q10 + 1.0/0.608));

        /* zeta_u */
        real_t zeta_u = BULK_VONKARM * bstar * z_wind / (ustar * ustar);
        if (fabs(zeta_u) > 10.0) zeta_u = copysign(10.0, zeta_u);
        real_t x2 = sqrt(fabs(1.0 - 16.0 * zeta_u));
        if (x2 < 1.0) x2 = 1.0;
        real_t x  = sqrt(x2);
        real_t psi_m_u, psi_h_u;
        if (zeta_u > 0.0) {
            psi_m_u = -5.0 * zeta_u;
            psi_h_u = -5.0 * zeta_u;
        } else {
            psi_m_u = log((1.0 + 2.0*x + x2) * (1.0 + x2) / 8.0)
                    - 2.0 * (atan(x) - atan(1.0));
            psi_h_u = 2.0 * log((1.0 + x2) / 2.0);
        }

        /* zeta_t */
        real_t zeta_t = BULK_VONKARM * bstar * z_tair / (ustar * ustar);
        if (fabs(zeta_t) > 10.0) zeta_t = copysign(10.0, zeta_t);
        x2 = sqrt(fabs(1.0 - 16.0 * zeta_t));
        if (x2 < 1.0) x2 = 1.0;
        x = sqrt(x2);
        real_t psi_m_t, psi_h_t;
        (void)psi_m_t;
        if (zeta_t > 0.0) {
            psi_m_t = -5.0 * zeta_t;
            psi_h_t = -5.0 * zeta_t;
        } else {
            psi_m_t = log((1.0 + 2.0*x + x2) * (1.0 + x2) / 8.0)
                    - 2.0 * (atan(x) - atan(1.0));
            psi_h_t = 2.0 * log((1.0 + x2) / 2.0);
        }

        /* zeta_q */
        real_t zeta_q = BULK_VONKARM * bstar * z_shum / (ustar * ustar);
        if (fabs(zeta_q) > 10.0) zeta_q = copysign(10.0, zeta_q);
        x2 = sqrt(fabs(1.0 - 16.0 * zeta_q));
        if (x2 < 1.0) x2 = 1.0;
        x = sqrt(x2);
        real_t psi_m_q, psi_h_q;
        (void)psi_m_q;
        if (zeta_q > 0.0) {
            psi_m_q = -5.0 * zeta_q;
            psi_h_q = -5.0 * zeta_q;
        } else {
            psi_m_q = log((1.0 + 2.0*x + x2) * (1.0 + x2) / 8.0)
                    - 2.0 * (atan(x) - atan(1.0));
            psi_h_q = 2.0 * log((1.0 + x2) / 2.0);
        }

        /* Shift wind/temperature/humidity to 10 m and reference levels. */
        u10 = u / (1.0 + cd_n10_rt * (log(z_wind / 10.0) - psi_m_u) / BULK_VONKARM);
        if (u10 < BULK_U10MIN) u10 = BULK_U10MIN;
        t10 = t - tstar / BULK_VONKARM * (log(z_tair / z_wind) + psi_h_u - psi_h_t);
        q10 = q - qstar / BULK_VONKARM * (log(z_shum / z_wind) + psi_h_u - psi_h_q);
        tv  = t10 * (1.0 + 0.608 * q10);

        /* Update neutral 10 m coefs (LY2009 eqn 11a/b again). */
        hl1 = (2.7/u10 + 0.142 + 0.0764*u10 - 3.14807e-10 * pow(u10, 6)) / 1.0e3;
        cd_n10 = (0.5 - copysign(0.5, u10 - 33.0)) * hl1
               + (0.5 + copysign(0.5, u10 - 33.0)) * 2.34e-3;
        cd_n10_rt = sqrt(cd_n10);
        ce_n10    = 34.6 * cd_n10_rt * 1.0e-3;
        stab      = 0.5 + copysign(0.5, zeta_u);
        ch_n10    = (18.0*stab + 32.7*(1.0 - stab)) * cd_n10_rt * 1.0e-3;

        /* Shift cd/ch/ce back to measurement height. */
        real_t xx = (log(z_wind / 10.0) - psi_m_u) / BULK_VONKARM;
        cd = cd_n10 / pow(1.0 + cd_n10_rt * xx, 2);
        xx = (log(z_wind / 10.0) - psi_h_u) / BULK_VONKARM;
        ch = ch_n10 / (1.0 + ch_n10 * xx / cd_n10_rt) * sqrt(cd / cd_n10);
        ce = ce_n10 / (1.0 + ce_n10 * xx / cd_n10_rt) * sqrt(cd / cd_n10);

        real_t test = fabs(cd - cd_prev) / (cd + 1.0e-8);
        cd_prev = cd;
        if (test < BULK_INC_RATIO) break;
    }

    *cd_out = cd;
    *ce_out = ce;
    *ch_out = ch;
}

/* ------------------------------------------------------------------------ *
 * Open-water heat & freshwater fluxes.                                     *
 * Literal port of obudget (ice_thermo_oce.F90:777-868).                    *
 * Standard saturation_shum_formula = .true. branch.                        *
 * Returns: qsr (downward shortwave to ocean, +ve down)                     *
 *          qns (non-solar surface heat, +ve up = ocean loses)              *
 *          evap (kg/m²/s, +ve = ocean loses)                               *
 * ------------------------------------------------------------------------ */
static void obudget(real_t qa,        /* shum [kg/kg] */
                    real_t fsh,       /* shortwave down [W/m²] */
                    real_t flo,       /* longwave  down [W/m²] */
                    real_t t,         /* surface ocean T [°C] */
                    real_t ug,        /* wind speed scalar [m/s] */
                    real_t ta,        /* air temp [°C] */
                    real_t ch,
                    real_t ce,
                    real_t *qsr_out,  /* downward shortwave (W/m², +ve down) */
                    real_t *qns_out,  /* non-solar surface heat (W/m², +ve up) */
                    real_t *evap_out) /* m/s, +ve up */
{
    /* Standard saturation specific humidity formula (Fortran 840-844). */
    const real_t c1 = 3.8e-3;
    const real_t c4 = 17.27;
    const real_t c5 = 237.3;
    real_t b = c1 * exp(c4 * t / (t + c5));

    real_t hfswrow  = (1.0 - BULK_ALBW) * fsh;
    real_t hflwrow  = flo;
    real_t hflwrdout= -BULK_EMISS_WAT * BULK_BOLTZMANN * pow(t + BULK_TMELT, 4);
    real_t hfsenow  = BULK_RHOAIR * BULK_CPAIR * ch * ug * (ta - t);
    real_t evap     = BULK_RHOAIR * ce * ug * (qa - b);   /* kg/m²/s */
    real_t hflatow  = BULK_CLHW * evap;

    /* qsr: downward shortwave to ocean. qns: non-solar surface heat (positive
     * up = ocean loses). The Fortran ice path collects all heat into ehf
     * (positive down) then negates in oce_fluxes(). For our forcing struct
     * convention (positive up = ocean loses), we put qsr separately and
     * make qns = -(LW_in + LW_out + sensible + latent). */
    real_t qns = -(hflwrow + hflwrdout + hfsenow + hflatow);
    *qsr_out  = hfswrow;
    *qns_out  = qns;
    *evap_out = evap * BULK_INV_RHOWAT;   /* m/s, +ve up */
}

/* ------------------------------------------------------------------------ *
 * Drive bulk formulae over all nodes; assemble forcing arrays.             *
 * ------------------------------------------------------------------------ */
void fesom_bulk_compute(const struct fesom_jra55  *jra,
                        const struct fesom_mesh   *mesh,
                        const struct fesom_dyn    *dyn,
                        const struct fesom_tracers *tracers,
                        struct fesom_forcing       *forcing,
                        struct fesom_ice           *ice,
                        struct fesom_partit        *partit)
{
    /* Loop bound is myDim+eDim — matches Fortran ncar_ocean_fluxes_mode
     * (gen_bulk_formulae.F90:181) and the wind-stress loop
     * (gen_forcing_couple.F90:738), both `do i=1,myDim_nod2D+eDim_nod2D`.
     * Computing the exchange coefficients (Cd/Ch/Ce) and fluxes at eDim from
     * the (already halo-consistent) JRA + ocean-surface inputs keeps
     * Ch_atm_oce/Ce_atm_oce consistent at eDim. The sea-ice thermodynamics
     * reads Ch/Ce directly over myDim+eDim; computing them on myDim only left
     * eDim stale -> divergent a_ice across ranks -> divergent stress_surf on
     * replicated boundary elements -> non-conservative SSH RHS -> dt blow-up. */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int E = mesh->myDim_elem2D;

    /* For each node compute cd/ce/ch then qsr/qns/emp. */
    for (int n = 0; n < N; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) {
            /* Cavity node — fluxes are 0. */
            forcing->stress_node_surf[2*n + 0] = 0.0;
            forcing->stress_node_surf[2*n + 1] = 0.0;
            forcing->heat_flux [n] = 0.0;
            forcing->water_flux[n] = 0.0;
            ice->stress_atmice_x[n] = 0.0;
            ice->stress_atmice_y[n] = 0.0;
            continue;
        }

        real_t T_oc = tracers->data[FESOM_TRACER_T].values[FESOM_NODE3D(n, 0, mesh->nl)];
        /* Surface ocean current at NODE — use uvnode's surface layer. */
        real_t u_w = dyn->uvnode[(size_t)(2 * (FESOM_NODE3D(n, 0, mesh->nl))) + 0];
        real_t v_w = dyn->uvnode[(size_t)(2 * (FESOM_NODE3D(n, 0, mesh->nl))) + 1];

        real_t ua = jra->u_wind[n];
        real_t va = jra->v_wind[n];
        real_t ta = jra->Tair[n];        /* °C */
        real_t qa = jra->shum[n];        /* kg/kg */
        real_t fsh= jra->shortwave[n];
        real_t flo= jra->longwave[n];
        real_t pr = jra->prec_rain[n];   /* m/s */
        real_t ps = jra->prec_snow[n];   /* m/s */

        real_t cd, ce, ch;
        ncar_ocean_fluxes_mode(ta, qa, ua, va, T_oc, u_w, v_w,
                               jra->z_wind, jra->z_tair, jra->z_shum,
                               &cd, &ce, &ch);
        /* save per-node ch/ce for the sea-ice thermodynamics path
         * (Fortran: Ch_atm_oce_arr / Ce_atm_oce_arr in g_forcing_arrays) */
        forcing->Ch_atm_oce[n] = ch;
        forcing->Ce_atm_oce[n] = ce;

        /* Wind speed for obudget — Fortran ice_thermo_oce uses sqrt(u_wind² + v_wind²). */
        real_t ug = sqrt(ua*ua + va*va);
        real_t qsr, qns, evap;
        obudget(qa, fsh, flo, T_oc, ug, ta, ch, ce, &qsr, &qns, &evap);

        /* Heat flux (positive up = ocean loses). qsr is downward → enters as
         * -qsr in ocean's surface heat budget. We sum into a single
         * heat_flux per analytical-forcing convention; if/when shortwave
         * penetration is added we'll separate. */
        forcing->heat_flux[n]  = qns - qsr;
        /* Freshwater (E - P). evap is m/s positive up (ocean loses);
         * prec is m/s positive down (ocean gains). */
        forcing->water_flux[n] = evap - pr - ps;

        /* Wind stress at node — Fortran gen_forcing_couple.F90 lines 749-757
         * (Swind=0 default in our namelist). */
        real_t dux = ua - u_w;
        real_t dvy = va - v_w;
        real_t mag = sqrt(dux*dux + dvy*dvy) * BULK_RHOAIR;
        forcing->stress_node_surf[2*n + 0] = cd * mag * dux;
        forcing->stress_node_surf[2*n + 1] = cd * mag * dvy;

        /* Atmosphere-ICE momentum stress — Fortran gen_forcing_couple.F90:
         * 759-763. Same bulk form but wind RELATIVE TO ICE and the constant
         * atm-ice drag Cd_atm_ice (AOMIP_drag_coeff=.false. on CORE2). This
         * is the EVP's direct wind forcing (inv_mass·stress_atmice in
         * fesom_ice_evp.c); without it the sea ice never feels the wind. */
        real_t dux_i = ua - ice->uice[n];
        real_t dvy_i = va - ice->vice[n];
        real_t mag_i = sqrt(dux_i*dux_i + dvy_i*dvy_i) * BULK_RHOAIR;
        ice->stress_atmice_x[n] = FESOM_CD_ATM_ICE * mag_i * dux_i;
        ice->stress_atmice_y[n] = FESOM_CD_ATM_ICE * mag_i * dvy_i;
    }

    /* Halo-exchange node-side fields written above so the element
     * interpolation below sees correct values at halo vertices. Without this,
     * partition-boundary element stress is interpolated against stale or
     * uninitialised halo entries → spurious wind torque on boundary triangles
     * → blow-up during momentum step. heat_flux/water_flux are read by
     * tracer_diff for myDim only, but exchange them too for any future kernel
     * that touches them at halo. */
    if (partit && partit->npes > 1) {
        fesom_halo_exchange(forcing->stress_node_surf, FESOM_HALO_NOD2D, 1, 2, partit);
        fesom_halo_exchange(forcing->heat_flux,        FESOM_HALO_NOD2D, 1, 1, partit);
        fesom_halo_exchange(forcing->water_flux,       FESOM_HALO_NOD2D, 1, 1, partit);
    }

    /* Interpolate node stress → element stress (mean of 3 vertices). The
     * existing analytical-forcing path writes stress_surf directly at elements;
     * the bulk formulae assemble at nodes (matches Fortran), then average. */
    memset(forcing->stress_surf, 0, (size_t)E * 2 * sizeof(real_t));
    for (int e = 0; e < E; ++e) {
        real_t sx = 0.0, sy = 0.0;
        for (int k = 0; k < 3; ++k) {
            int v = mesh->elem_nodes[3*e + k];
            sx += forcing->stress_node_surf[2*v + 0];
            sy += forcing->stress_node_surf[2*v + 1];
        }
        forcing->stress_surf[2*e + 0] = sx / 3.0;
        forcing->stress_surf[2*e + 1] = sy / 3.0;
    }
}

/*===========================================================================
 * cal_shortwave_rad — literal port of oce_shortwave_pene.F90.
 *
 * Shortwave penetration into the ocean (Morel & Antoine 1994, Sweeney 2005).
 * Builds the visible-band shortwave TEMPERATURE flux sw_3d through the column
 * (two-band exponential with chl-dependent coefficients) and removes the
 * visible fraction (0.54) from heat_flux. No penetration under cavity or sea
 * ice. MUST be called AFTER ice→ocean coupling and BEFORE the temperature
 * tracer equation consumes heat_flux + sw_3d (oce_ale_tracer.F90:990).
 *===========================================================================*/
void fesom_cal_shortwave_rad(const struct fesom_mesh  *mesh,
                             const struct fesom_jra55 *jra,
                             const struct fesom_ice   *ice,
                             struct fesom_forcing     *forcing)
{
    if (!FESOM_PHASE1_USE_SW_PENE) return;

    const int    N    = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int    nl   = mesh->nl;
    const real_t albw = (real_t)BULK_ALBW;     /* open-water albedo (CORE2 namelist.ice albw=0.1) */
    const real_t vcpw = (real_t)FESOM_VCPW;
    const real_t *a_ice = ice->data[FESOM_ICE_AICE].values;
    real_t *sw_3d = forcing->sw_3d;
    real_t *chl   = forcing->chl;

    /* zero sw_3d over all local nodes/levels (Fortran 39-43) */
    memset(sw_3d, 0, (size_t)N * (size_t)nl * sizeof(real_t));

    for (int n2 = 0; n2 < N; ++n2) {
        if (mesh->ulevels_nod2D[n2] > 1) continue;   /* cavity: no penetration (F:51) */
        if (a_ice[n2] > 0.0)             continue;   /* under ice: none      (F:52) */

        /* visible shortwave into ocean [W/m²]; '+'-up: add back to heat_flux (F:56-60) */
        real_t swsurf = (1.0 - albw) * jra->shortwave[n2];   /* = qsr */
        swsurf *= 0.54;                                      /* visible part (300-750nm) */
        forcing->heat_flux[n2] += swsurf;

        /* Sweeney 2005 (Appendix A) two-band coefficients from chl (F:66-79) */
        real_t cc = chl[n2];
        if (cc < 0.02) cc = 0.02;                            /* limit from below */
        real_t c  = log10(cc);
        real_t c2 = c*c, c3 = c2*c, c4 = c3*c, c5 = c4*c;
        real_t v1  = 0.008*c + 0.132*c2 + 0.038*c3 - 0.017*c4 - 0.007*c5;
        real_t v2  = 0.679 - v1;
        v1         = 0.321 + v1;
        real_t sc1 = 1.54  - 0.197*c + 0.166*c2 - 0.252*c3 - 0.055*c4 + 0.042*c5;
        real_t sc2 = 7.925 - 6.644*c + 3.662*c2 - 1.815*c3 - 0.218*c4 + 0.502*c5;

        swsurf /= vcpw;                                      /* W/m² → K m/s (F:81) */

        const int nzmin = mesh->ulevels_nod2D[n2] - 1;       /* 0-based top interface */
        const int nzmax = mesh->nlevels_nod2D[n2] - 1;       /* 0-based bottom interface */
        sw_3d[FESOM_NODE3D(n2, nzmin, nl)] = swsurf;         /* Fortran nzmin (F:85) */
        for (int k = nzmin + 1; k <= nzmax; ++k) {           /* F:86-93 */
            real_t z   = mesh->zbar_3d_n[FESOM_NODE3D(n2, k, nl)];
            real_t aux = v1 * exp(z / sc1) + v2 * exp(z / sc2);
            sw_3d[FESOM_NODE3D(n2, k, nl)] = swsurf * aux;
            if (aux < 1.0e-5 || k == nzmax) {
                sw_3d[FESOM_NODE3D(n2, k, nl)] = 0.0;
                break;
            }
        }
    }
}
