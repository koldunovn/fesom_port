#include "fesom_ice_thermo.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_halo.h"
#include "fesom_jra55.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_sss_runoff.h"
#include "fesom_tracers.h"
#include "fesom_types.h"

#include <math.h>

/*
 * Mirror of Fortran TFrez(S) at ice_thermo_oce.F90:899.
 *   TFrez = -0.0575*S + 1.7105e-3 * sqrt(S**3) - 2.155e-4 * S*S
 * Millero (1978) / UNESCO; see Gill (1982).
 */
real_t fesom_ice_tfrez(real_t salinity)
{
    real_t s = salinity;
    return -0.0575 * s + 1.7105e-3 * sqrt(s * s * s) - 2.155e-4 * s * s;
}

/*
 * Mirror of Fortran flooding at ice_thermo_oce.F90:872.
 * Archimedes freeboard test; converts snow to ice when freeboard < 0.
 */
void fesom_ice_flooding(const fesom_ice_thermo *th,
                        real_t                 *h,
                        real_t                 *hsn)
{
    /* hdraft: displaced water (m) = (rhosno*hsn + h*rhoice) / rhowat */
    real_t hdraft = (th->rhosno * (*hsn) + (*h) * th->rhoice) * th->inv_rhowat;
    /* hflood: increase in ice thickness due to flooding */
    real_t hmin = (hdraft < *h) ? hdraft : *h;
    real_t hflood = hdraft - hmin;
    /* convert: add flooded volume to ice, remove mass-equivalent from snow */
    *h   += hflood;
    *hsn -= hflood * th->rhoice * th->inv_rhosno;
}

/*
 * Mirror of Fortran cut_off at ice_thermo_oce.F90:70.
 * Fortran loop bound: 1..myDim_nod2D+eDim_nod2D — halo included.
 */
void fesom_ice_cut_off(fesom_ice            *ice,
                       struct fesom_partit  *partit,
                       struct fesom_mesh    *mesh)
{
    (void)partit;
    real_t *a_ice  = ice->data[FESOM_ICE_AICE].values;
    real_t *m_ice  = ice->data[FESOM_ICE_MICE].values;
    real_t *m_snow = ice->data[FESOM_ICE_MSNOW].values;

    /* matches ice_thermo_oce.F90:102 — loop bound is myDim+eDim, halo included */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    for (int n = 0; n < N; ++n) {
        /* upper cutoff on a_ice */
        if (a_ice[n] > 1.0) a_ice[n] = 1.0;
        /* lower cutoff on a_ice: zero all three if near-zero */
        if (a_ice[n] < 1.0e-9) {
            a_ice [n] = 0.0;
            m_ice [n] = 0.0;
            m_snow[n] = 0.0;
        }
        /* lower cutoff on m_ice: zero all three if near-zero */
        if (m_ice[n] < 1.0e-9) {
            m_ice [n] = 0.0;
            m_snow[n] = 0.0;
            a_ice [n] = 0.0;
        }
    }
}

/*
 * Mirror of Fortran obudget at ice_thermo_oce.F90:777.
 *
 * Flow:
 *   hfswrow  = (1 - albw) * fsh                          net shortwave
 *   hflwrow  = flo                                       incoming longwave
 *   hflwrdout= -emiss_wat * sigma * (t + tmelt)^4        outgoing longwave
 *   hfradow  = sum of the above                          total radiation
 *   hfsenow  = rhoair * cpair * ch * ug * (ta - t)       sensible heat
 *   evap_kgm2s = rhoair * ce * ug * (qa - b)             evaporation mass flux
 *   hflatow  = clhw * evap_kgm2s                         latent heat
 *   fh       = -(hfradow + hfsenow + hflatow) / cl       growth rate in m ice/s
 *   evap     = evap_kgm2s * inv_rhowat                   convert to m water/s
 *
 * where b is the saturated surface specific humidity
 *   b = c1 * exp(c4 * t / (t + c5))       (standard formula; Fortran default)
 * with c1=3.8e-3, c4=17.27, c5=237.3.
 *
 * The `open_water_albedo > 0` branches (time-dependent Taylor / Briegleb
 * albedo) are out of scope for this port — CORE2 default is 0.
 */
void fesom_ice_obudget(const fesom_ice_thermo *th,
                       real_t qa, real_t fsh, real_t flo,
                       real_t t,  real_t ug,  real_t ta,
                       real_t ch, real_t ce,
                       real_t *fh_out,      real_t *evap_out,
                       real_t *hflatow,     real_t *hfsenow,
                       real_t *hflwrdout,   real_t *hfswrow,
                       real_t *hflwrow,     real_t *hfradow)
{
    FESOM_CHECK(th->open_water_albedo == 0,
                "fesom_ice_obudget: open_water_albedo=%d not supported "
                "(only 0 = fixed albw is ported)", th->open_water_albedo);

    /* saturated surface specific humidity — standard formula branch */
    const real_t c1 = 3.8e-3;
    const real_t c4 = 17.27;
    const real_t c5 = 237.3;
    real_t b = c1 * exp(c4 * t / (t + c5));

    /* radiation heat fluxes, W/m² */
    real_t _hfswr   = (1.0 - th->albw) * fsh;
    real_t _hflwr   = flo;
    real_t _hflwrd  = -th->emiss_wat * th->boltzmann * pow(t + th->tmelt, 4.0);
    real_t _hfrad   = _hfswr + _hflwr + _hflwrd;

    /* sensible heat, W/m² */
    real_t _hfsen = th->rhoair * th->cpair * ch * ug * (ta - t);

    /* evaporation in kg/(m²·s), then latent heat */
    real_t _evap_kgm2s = th->rhoair * ce * ug * (qa - b);
    real_t _hflat = th->clhw * _evap_kgm2s;

    real_t _hftot = _hfrad + _hfsen + _hflat;

    /* outputs */
    *fh_out   = -_hftot / th->cl;          /* m ice / s */
    *evap_out = _evap_kgm2s * th->inv_rhowat;  /* m water / s, negative up */

    if (hfswrow)   *hfswrow   = _hfswr;
    if (hflwrow)   *hflwrow   = _hflwr;
    if (hflwrdout) *hflwrdout = _hflwrd;
    if (hfradow)   *hfradow   = _hfrad;
    if (hfsenow)   *hfsenow   = _hfsen;
    if (hflatow)   *hflatow   = _hflat;
}

/*
 * Mirror of Fortran budget at ice_thermo_oce.F90:657.
 *
 * Albedo selection (line 718-730):
 *   t < 0:  freezing → albsn (snow present) | albi  (no snow)
 *   t ≥ 0:  melting  → albsnm                | albim
 *
 * 5-iteration Newton-Raphson on t (line 741-752):
 *   Residual A1+A2+C, derivative A3, update t += residual/derivative.
 *   B = q1/rhoair * exp(q2/(t+tmelt)) is the saturated specific humidity over ice.
 *   q1 = 11637800.0; q2 = -5897.8 (these are different from obudget's q1,q2
 *   because they parameterise saturation over ice, not water).
 *
 * Final fluxes (line 757-770) match obudget structure but use ice albedo and
 * ice latent heat (clhi instead of clhw) and write `subli` instead of `evap`.
 */
void fesom_ice_budget(const fesom_ice_thermo *th,
                      real_t hice, real_t hsn,
                      real_t *t,
                      real_t ta, real_t qa,
                      real_t fsh, real_t flo,
                      real_t ug, real_t S_oc,
                      real_t ch_i, real_t ce_i,
                      real_t *fh, real_t *subli)
{
    /* albedo selection */
    real_t alb;
    if (*t < 0.0) {
        alb = (hsn > 0.0) ? th->albsn  : th->albi;
    } else {
        alb = (hsn > 0.0) ? th->albsnm : th->albim;
    }

    const real_t q1 = 11637800.0;
    const real_t q2 = -5897.8;
    const int    imax = 5;

    real_t d1 = th->rhoair * th->cpair * ch_i;
    real_t d2 = th->rhoair * ce_i;
    real_t d3 = d2 * th->clhi;

    /* total incoming atmospheric heat flux (constant across iterations) */
    real_t A1 = (1.0 - alb) * fsh + flo + d1 * ug * ta + d3 * ug * qa;

    /* Newton-Raphson on t */
    real_t B = 0.0;  /* outlives loop for the final flux pass */
    for (int iter = 0; iter < imax; ++iter) {
        real_t tk = *t + th->tmelt;
        B = q1 * th->inv_rhoair * exp(q2 / tk);
        real_t A2 = -d1 * ug * (*t)
                    -d3 * ug * B
                    -th->emiss_ice * th->boltzmann * (tk * tk * tk * tk);
        real_t A3 = -d3 * ug * B * q2 / (tk * tk);
        real_t C  = th->con / hice;
        A3 += C + d1 * ug
              + 4.0 * th->emiss_ice * th->boltzmann * (tk * tk * tk);
        C *= (fesom_ice_tfrez(S_oc) - *t);
        *t += (A1 + A2 + C) / A3;
    }
    if (*t > 0.0) *t = 0.0;

    /* final fluxes (W/m²) using the converged t and last-iteration B */
    real_t tk = *t + th->tmelt;
    real_t hfrad = (1.0 - alb) * fsh
                   + flo
                   - th->emiss_ice * th->boltzmann * (tk * tk * tk * tk);
    real_t hfsen = d1 * ug * (ta - *t);
    real_t _subli = d2 * ug * (qa - B);
    real_t hflat  = th->clhi * _subli;
    real_t hftot  = hfrad + hfsen + hflat;

    *fh    = -hftot / th->cl;
    *subli = _subli * th->inv_rhowat;   /* negative upward */
}

/*
 * Mirror of Fortran therm_ice at ice_thermo_oce.F90:367.
 *
 * Phase B4a scope (per .h):
 *   - Skip flooding integration (lines 637-649)         → iflice = 0
 *   - Only use_virt_salt = false branch (lines 614-617)
 *   - Skip new_iclasses=true branch (no hpdf in struct)
 *   - Always assume snowdist = true (Fortran default)
 *
 * The body follows the Fortran line-by-line. Variable names match Fortran
 * to keep diff visible; obvious in/out semantics use C pointer style.
 */
static inline real_t fmin_w(real_t a, real_t b) { return a < b ? a : b; }
static inline real_t fmax_w(real_t a, real_t b) { return a > b ? a : b; }

void fesom_therm_ice(const fesom_ice_thermo *th,
                     int use_virt_salt,
                     real_t *h, real_t *hsn, real_t *A,
                     real_t fsh, real_t flo, real_t Ta, real_t qa,
                     real_t rain, real_t snow, real_t runo, real_t rsss,
                     real_t ug, real_t ustar, real_t T_oc, real_t S_oc,
                     real_t H_ML, real_t *t, real_t ice_dt,
                     real_t ch, real_t ce, real_t ch_i, real_t ce_i,
                     real_t evap_in,
                     real_t *fw, real_t *fwice, real_t *fwsnw,
                     real_t *ehf, real_t *evap, real_t *rsf,
                     real_t *dhgrowth, real_t *dhsngrowth, real_t *dAgrowth,
                     real_t *iflice,
                     real_t *hflatow, real_t *hfsenow, real_t *hflwrdout,
                     real_t *hfswrow, real_t *hflwrow, real_t *hfradow,
                     real_t lid_clo, real_t geolon, real_t geolat,
                     real_t *subli)
{
    (void)evap_in;   /* Fortran ignores the value too — see line 472 */

    /* B4a: assert snowdist=true (Fortran default for CORE2) */
    FESOM_CHECK(th->snowdist == 1, "fesom_therm_ice: snowdist=false not ported");

    /* dhgrowth holds the pre-step ice mass for end-of-step diff */
    real_t _dhgrowth = *h;

    /* Effective snow thickness (Semtner 1976 0-layer); always added under snowdist=true */
    real_t snthick = (*hsn) * (th->con / th->consn) / fmax_w(*A, th->armin);
    real_t thick   = (*h) / fmax_w(*A, th->armin);
    thick = snthick + thick;

    /* --- Open-ocean growth rate (single obudget call) --- */
    real_t rhow = 0.0, _evap = 0.0;
    real_t _hflatow = 0.0, _hfsenow = 0.0, _hflwrdout = 0.0;
    real_t _hfswrow = 0.0, _hflwrow = 0.0, _hfradow   = 0.0;
    fesom_ice_obudget(th, qa, fsh, flo, T_oc, ug, Ta, ch, ce,
                      &rhow, &_evap,
                      &_hflatow, &_hfsenow, &_hflwrdout,
                      &_hfswrow, &_hflwrow, &_hfradow);
    /* Scale radiative components by open-water fraction */
    real_t one_minus_A = 1.0 - (*A);
    _hflatow   *= one_minus_A;
    _hfsenow   *= one_minus_A;
    _hfradow   *= one_minus_A;
    _hfswrow   *= one_minus_A;
    _hflwrow   *= one_minus_A;
    _hflwrdout *= one_minus_A;
    if (hflatow)   *hflatow   = _hflatow;
    if (hfsenow)   *hfsenow   = _hfsenow;
    if (hflwrdout) *hflwrdout = _hflwrdout;
    if (hfswrow)   *hfswrow   = _hfswrow;
    if (hflwrow)   *hflwrow   = _hflwrow;
    if (hfradow)   *hfradow   = _hfradow;

    /* --- Ice-covered growth rate (iclasses loop, Hibler 1984) --- */
    real_t rhice = 0.0;
    real_t _subli = 0.0;
    if (thick > th->hmin) {
        for (int k = 1; k <= th->iclasses; ++k) {
            real_t thact = (real_t)(2*k - 1) * thick / (real_t)th->iclasses;
            /* snowdist=true → no per-class snthick adjustment */
            real_t shice = 0.0, subli_i = 0.0;
            fesom_ice_budget(th, thact, *hsn, t, Ta, qa, fsh, flo,
                             ug, S_oc, ch_i, ce_i, &shice, &subli_i);
            rhice  += shice;
            _subli += subli_i;
        }
        rhice  /= (real_t)th->iclasses;
        _subli /= (real_t)th->iclasses;
    }

    /* Convert growth rates [m/s] -> per-step [m/DT] */
    rhow  *= ice_dt;
    rhice *= ice_dt;

    /* Areal weighting */
    real_t show  = rhow  * one_minus_A;
    real_t shice = rhice * (*A);
    real_t sh    = show + shice;

    /* Atmospheric heat flux averaged over grid cell */
    real_t ahf = -th->cl * sh / ice_dt;

    /* Precipitation into ocean (m water/s); includes runoff */
    real_t prec = rain + runo + snow * one_minus_A;

    /* Snow fall above ice — accumulate */
    *hsn += snow * ice_dt * (*A) * 1000.0 * th->inv_rhosno;
    real_t _dhsngrowth = *hsn;

    _evap   *= one_minus_A;          /* open-water evap weighted */
    _subli  *= (*A);                  /* ice sublimation weighted */

    /* Snow melt by negative atmospheric heat flux (sh<0 → melting) */
    real_t hsntmp = -fmin_w(sh, 0.0) * th->rhoice * th->inv_rhosno;
    hsntmp = fmin_w(hsntmp, *hsn);
    *hsn -= hsntmp;

    /* Negative atmospheric heat flux remaining after snow melt */
    real_t rh = sh + hsntmp * th->rhosno / th->rhoice;
    *h = fmax_w(*h, 0.0);

    /* Ocean-to-ice heat flux (ustar parameterisation) */
    real_t Tfrez_S = fesom_ice_tfrez(S_oc);
    real_t o2ihf = (T_oc - Tfrez_S) * 0.006 * ustar * th->cc * (*A)
                 + (T_oc - Tfrez_S) * H_ML / ice_dt * th->cc * one_minus_A;
    rh -= o2ihf * ice_dt / th->cl;
    real_t qhst = *h + rh;

    /* Melt remaining snow if ML heat content left over */
    real_t sn = *hsn + fmin_w(qhst, 0.0) * th->rhoice * th->inv_rhosno;
    sn = fmax_w(sn, 0.0);
    *hsn = sn;
    *h   = fmax_w(qhst, 0.0);
    if (*h < 1.0e-6) *h = 0.0;

    /* Compute changes per second */
    _dhgrowth   = (*h   - _dhgrowth)   / ice_dt;
    _dhsngrowth = (*hsn - _dhsngrowth) / ice_dt;
    if (dhgrowth)   *dhgrowth   = _dhgrowth;
    if (dhsngrowth) *dhsngrowth = _dhsngrowth;

    /* Net surface heat flux to ocean */
    *ehf = ahf + th->cl * (_dhgrowth + (th->rhosno / th->rhoice) * _dhsngrowth);

    /* Freshwater fluxes — Fortran lines 612-623 */
    real_t _fwice, _fwsnw, _fw, _rsf;
    if (!use_virt_salt) {
        /* real-salt path (use_virt_salt = false): salt leaves with melted ice */
        _fwice = -_dhgrowth   * th->rhoice * th->inv_rhowat;
        _fwsnw = -_dhsngrowth * th->rhosno * th->inv_rhowat;
        _fw    = prec + _evap + _fwice + _fwsnw;
        _rsf   = _fwice * th->Sice;
    } else {
        /* virtual-salt path (linfs default): scale fwice by (rsss-Sice)/rsss */
        _fwice = -_dhgrowth   * th->rhoice * th->inv_rhowat * (rsss - th->Sice) / rsss;
        _fwsnw = -_dhsngrowth * th->rhosno * th->inv_rhowat;
        _fw    = prec + _evap + _fwice + _fwsnw;
        _rsf   = 0.0;
    }
    if (fwice) *fwice = _fwice;
    if (fwsnw) *fwsnw = _fwsnw;

    /* Compactness change (Hibler 1979 eq. 16) */
    rh = -fmin_w(*h, -rh);
    real_t rA  = rhow - o2ihf * ice_dt / th->cl;
    real_t Aold = *A;
    *A = *A + th->c_melt * fmin_w(rh, 0.0) * (*A) / fmax_w(*h, th->hmin)
            + fmax_w(rA, 0.0) * (1.0 - *A) / lid_clo;
    *A = fmin_w(*A, *h * 1.0e6);
    *A = fmin_w(fmax_w(*A, 0.0), 1.0);
    if (dAgrowth) *dAgrowth = (*A - Aold) / ice_dt;

    /* Flooding (snow-to-ice conversion) — Fortran lines 638-649 */
    real_t _iflice = *h;
    fesom_ice_flooding(th, h, hsn);
    _iflice = (*h - _iflice) / ice_dt;
    if (iflice) *iflice = _iflice;

    /*
     * Salt conservation correction for the flooding step.
     * Snow-derived ice carries no salt; without correction the system would
     * "produce" salt as snow becomes ice. The correction adjusts rsf (real
     * salt path) or fw (virtual salt path) to compensate.
     */
    if (!use_virt_salt) {
        _rsf -= _iflice * th->rhoice * th->inv_rhowat * th->Sice;
    } else {
        _fw  += _iflice * th->rhoice * th->inv_rhowat * th->Sice / rsss;
    }
    if (rsf) *rsf = _rsf;
    if (fw)  *fw  = _fw;

    if (evap)  *evap  = _evap + _subli;
    if (subli) *subli = _subli;

    (void)geolon; (void)geolat;  /* reserved for Taylor/Briegleb branch (out of scope) */
}

/*
 * Per-node sea-ice thermodynamics driver.
 * Mirror of Fortran thermodynamics at ice_thermo_oce.F90:148.
 *
 * Cavity skip uses ulevels_nod2D > 1.
 * Phase B: l_snow=true assumed (CORE2 has prec_snow). The l_snow=false branch
 * (synth snow from Tair) lands when JRA55 doesn't provide snow.
 */
void fesom_ice_thermodynamics(fesom_ice                     *ice,
                              struct fesom_partit           *partit,
                              struct fesom_mesh             *mesh,
                              const struct fesom_forcing    *forcing,
                              const struct fesom_jra55      *jra,
                              const struct fesom_sss_runoff *sr)
{
    int nl = mesh->nl;
    real_t  rsss_default = sr->ref_sss;
    int     ref_sss_local = sr->ref_sss_local;
    real_t  h_ml = ice->thermo.h_ml;

    /* ---- Loop 1 (myDim only): friction velocity ---- */
    /* Fortran ice_thermo_oce.F90:230-235 — myDim only, then exchange */
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) {
            ice->thermo.ustar[n] = 0.0;
            continue;
        }
        real_t du = ice->uice[n] - ice->srfoce_u[n];
        real_t dv = ice->vice[n] - ice->srfoce_v[n];
        ice->thermo.ustar[n] = sqrt((du*du + dv*dv) * ice->cd_oce_ice);
    }
    /* Halo exchange — Fortran line 238 (call exchange_nod) */
    fesom_exchange_nod2D(ice->thermo.ustar, partit);

    /* ---- Loop 2 (myDim+eDim, halo included): per-node therm_ice ---- */
    /* Fortran line 244 — loop bound includes halo, Fortran writes halo entries */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    for (int n = 0; n < N; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) continue;   /* cavity skip */

        /* Inputs */
        real_t h    = ice->data[FESOM_ICE_MICE].values[n];
        real_t hsn  = ice->data[FESOM_ICE_MSNOW].values[n];
        real_t A    = ice->data[FESOM_ICE_AICE].values[n];
        real_t fsh  = jra->shortwave[n];
        real_t flo  = jra->longwave[n];
        real_t Ta   = jra->Tair[n];
        real_t qa   = jra->shum[n];
        real_t rain = jra->prec_rain[n];
        real_t snow = jra->prec_snow[n];
        real_t runo = forcing->runoff[n];
        real_t ug   = sqrt(jra->u_wind[n]*jra->u_wind[n]
                         + jra->v_wind[n]*jra->v_wind[n]);
        real_t ustar = ice->thermo.ustar[n];
        real_t T_oc  = ice->srfoce_temp[n];
        real_t S_oc  = ice->srfoce_salt[n];
        real_t rsss  = ref_sss_local ? S_oc : rsss_default;
        real_t t     = ice->thermo.t_skin[n];
        real_t ch    = forcing->Ch_atm_oce[n];
        real_t ce    = forcing->Ce_atm_oce[n];
        real_t lid_clo = (mesh->geo_coord_nod2D[2*n + 1] > 0.0)
                            ? ice->thermo.h0
                            : ice->thermo.h0_s;
        real_t geolon = mesh->geo_coord_nod2D[2*n + 0];
        real_t geolat = mesh->geo_coord_nod2D[2*n + 1];

        real_t fw = 0, fwice = 0, fwsnw = 0, ehf = 0, evap = 0, rsf = 0;
        real_t ithdgr = 0, ithdgrsn = 0, ithdgra = 0, iflice = 0;
        real_t hflatow = 0, hfsenow = 0, hflwrdout = 0;
        real_t hfswrow = 0, hflwrow = 0, hfradow = 0, subli = 0;

        fesom_therm_ice(&ice->thermo,
                        sr->use_virt_salt,
                        &h, &hsn, &A,
                        fsh, flo, Ta, qa, rain, snow, runo, rsss,
                        ug, ustar, T_oc, S_oc, h_ml, &t, ice->ice_dt,
                        ch, ce, FESOM_CH_ATM_ICE, FESOM_CE_ATM_ICE,
                        /*evap_in=*/0.0,
                        &fw, &fwice, &fwsnw, &ehf, &evap, &rsf,
                        &ithdgr, &ithdgrsn, &ithdgra, &iflice,
                        &hflatow, &hfsenow, &hflwrdout,
                        &hfswrow, &hflwrow, &hfradow,
                        lid_clo, geolon, geolat, &subli);

        /* Backup of old values (Fortran lines 311-314) */
        ice->data[FESOM_ICE_MICE].values_old [n] = ice->data[FESOM_ICE_MICE].values [n];
        ice->data[FESOM_ICE_MSNOW].values_old[n] = ice->data[FESOM_ICE_MSNOW].values[n];
        ice->data[FESOM_ICE_AICE].values_old [n] = ice->data[FESOM_ICE_AICE].values [n];
        ice->thermo.thdgr_old[n]                 = ice->thermo.thdgr[n];

        /* New values */
        ice->data[FESOM_ICE_MICE].values [n] = h;
        ice->data[FESOM_ICE_MSNOW].values[n] = hsn;
        ice->data[FESOM_ICE_AICE].values [n] = A;

        ice->thermo.t_skin[n]   = t;
        ice->flx_fw[n]          = fw;        /* + down (will be negated in oce_fluxes) */
        ice->flx_h [n]          = ehf;       /* + down */
        ice->thermo.thdgr  [n]  = ithdgr;
        ice->thermo.thdgrsn[n]  = ithdgrsn;
        ice->thermo.thdgra [n]  = ithdgra;

        /* Diagnostics intentionally not stored in this Phase — we have no
         * arrays for evap/ice_sublimation/real_salt_flux/fw_ice/fw_snw/hf_Q*
         * yet. They land when oce_fluxes (Phase C) wires them up. */
        (void)evap; (void)rsf; (void)iflice;
        (void)fwice; (void)fwsnw;
        (void)hflatow; (void)hfsenow; (void)hflwrdout;
        (void)hfswrow; (void)hflwrow; (void)hfradow;
        (void)subli;
    }

    (void)nl;  /* unused for now; reserved if surface-level extraction is moved here */
}
