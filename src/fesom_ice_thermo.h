#ifndef FESOM_ICE_THERMO_H
#define FESOM_ICE_THERMO_H

#include "fesom_ice_types.h"

struct fesom_mesh;
struct fesom_partit;
struct fesom_dyn;
struct fesom_tracers;
struct fesom_forcing;
struct fesom_jra55;
struct fesom_sss_runoff;

/*
 * Sea-ice thermodynamics — non-ICEPACK (standard) path.
 * Mirrors Fortran ice_thermo_oce.F90. Scope per docs/plans/20260425-sea-ice-port.md:
 *   - no __icepack, no __oifs/__ifsinterface (so no ice_temp as 4th tracer)
 *   - no cavity (use_cavity=false)
 *   - no meltponds (use_meltponds=false)
 *   - no wiso tracers
 *
 * Phase B1: these three small building blocks. The per-node thermodynamics
 * driver and therm_ice land in later tasks.
 */

/*
 * Freezing-point temperature of seawater, Millero (1978) / UNESCO.
 * Mirror of Fortran TFrez(S) at ice_thermo_oce.F90:899.
 * Returns °C (or K offset) given salinity in PSU.
 */
real_t fesom_ice_tfrez(real_t salinity);

/*
 * Snow→ice conversion via flooding. Single-cell, in-place update of h and hsn.
 * Mirror of Fortran flooding at ice_thermo_oce.F90:872.
 *
 * Applies Archimedes: if the ice freeboard (h - hdraft where hdraft is the
 * displaced water from ice+snow mass) goes negative, convert enough snow to
 * ice that the freeboard is zero. Mass is conserved (converted snow volume
 * is scaled by rhoice/rhosno to preserve mass).
 */
void fesom_ice_flooding(const fesom_ice_thermo *th,
                        real_t                 *h,
                        real_t                 *hsn);

/*
 * Per-node sanity cutoffs for the ice tracers.
 * Mirror of Fortran cut_off at ice_thermo_oce.F90:70.
 *   - clamp a_ice to [0, 1]
 *   - if a_ice or m_ice falls below 1e-9, zero all three tracers at that node
 * Fortran loop bound: myDim_nod2D + eDim_nod2D (halo included).
 */
void fesom_ice_cut_off(fesom_ice            *ice,
                       struct fesom_partit  *partit,
                       struct fesom_mesh    *mesh);

/*
 * Open-water surface energy budget (single cell).
 * Mirror of Fortran obudget at ice_thermo_oce.F90:777.
 *
 * Inputs:
 *   qa  — specific humidity (kg/kg)
 *   fsh — incoming shortwave radiation (W/m²)
 *   flo — incoming longwave radiation  (W/m²)
 *   t   — sea surface temperature (°C)
 *   ug  — wind speed (m/s)
 *   ta  — air temperature (°C)
 *   ch  — sensible-heat transfer coefficient
 *   ce  — evaporation transfer coefficient
 *
 * Outputs (all may be NULL except fh and evap):
 *   fh        — ice growth rate (m ice / s; + = melt, − = grow)
 *   evap      — evaporation rate (m water / s, negative up)
 *   hflatow   — latent heat flux   (W/m²)
 *   hfsenow   — sensible heat flux (W/m²)
 *   hflwrdout — outgoing longwave  (W/m², negative)
 *   hfswrow   — net shortwave      (W/m²)
 *   hflwrow   — incoming longwave  (W/m²)
 *   hfradow   — total radiative    (W/m²)
 *
 * Scope cut: the `open_water_albedo > 0` branches (Taylor 1996 / Briegleb et al.)
 * need the model clock + solar-zenith routine, neither of which is ported yet.
 * CORE2 default is 0 (fixed albw from namelist), so we assert that.
 */
void fesom_ice_obudget(const fesom_ice_thermo *th,
                       real_t qa, real_t fsh, real_t flo,
                       real_t t,  real_t ug,  real_t ta,
                       real_t ch, real_t ce,
                       real_t *fh,        real_t *evap,
                       real_t *hflatow,   real_t *hfsenow,
                       real_t *hflwrdout, real_t *hfswrow,
                       real_t *hflwrow,   real_t *hfradow);

/*
 * Ice-surface skin-temperature solver + thick-ice growth/melt budget.
 * Mirror of Fortran budget at ice_thermo_oce.F90:657.
 *
 * Iterative Newton–Raphson (5 iterations) on the surface temperature t
 * to balance: net atmospheric heat in − sensible − latent − outgoing LW
 * + downward conduction (con/hice * (TFrez(S_oc) − t)) = 0.
 *
 * In/out:
 *   t — surface temperature in °C, updated by the iteration; clamped to ≤ 0
 *
 * Outputs:
 *   fh    — ice growth rate (m ice / s; + = melt, − = grow)
 *   subli — sublimation rate (m water / s, negative up)
 *
 * Inputs:
 *   hice — ice thickness, m
 *   hsn  — snow thickness, m (used only to pick albedo)
 *   ta   — air temperature, °C
 *   qa   — specific humidity, kg/kg
 *   fsh  — incoming shortwave, W/m²
 *   flo  — incoming longwave, W/m²
 *   ug   — wind speed, m/s
 *   S_oc — ocean salinity, ppt (used in TFrez for the ice-base temperature)
 *   ch_i, ce_i — sensible / evaporation transfer coefficients for ice
 */
void fesom_ice_budget(const fesom_ice_thermo *th,
                      real_t hice, real_t hsn,
                      real_t *t,
                      real_t ta, real_t qa,
                      real_t fsh, real_t flo,
                      real_t ug, real_t S_oc,
                      real_t ch_i, real_t ce_i,
                      real_t *fh, real_t *subli);

/*
 * Single-cell sea-ice thermodynamic growth/melt step.
 * Mirror of Fortran therm_ice at ice_thermo_oce.F90:367.
 *
 * Phase B4a scope: full body EXCEPT the flooding integration (lines 637-649)
 * and the use_virt_salt=true branch (lines 619-622). Both land in B4b.
 *
 * Note on runoff: included automatically via `prec = rain + runo + snow*(1-A)`
 * (line 526) and `fw = prec + evap + fwice + fwsnw` (line 616). No separate
 * "runoff fold" step — it's structural to therm_ice.
 *
 * Note on iclasses: the multi-class loop is preserved (CORE2 default = 7).
 * The `new_iclasses=true` branch (h_cutoff scaling, hpdf weighting) is NOT
 * ported — that field is not in struct fesom_ice_thermo.
 *
 * In/out:
 *   h    — ice mass per area, m
 *   hsn  — snow mass per area, m
 *   A    — ice compactness (area fraction), 0..1
 *   t    — temperature of snow/ice top surface, °C
 *
 * Outputs (must all be non-NULL except optional radiative flux outputs):
 *   fw      — total freshwater flux to ocean (m water / ice_dt; +up)
 *   fwice   — freshwater flux from ice melting (m water / ice_dt)
 *   fwsnw   — freshwater flux from snow melting (m water / ice_dt)
 *   ehf     — net surface heat flux to ocean (W/m², +down)
 *   evap    — evaporation + sublimation, m water / s
 *   rsf     — real salt flux (use_virt_salt=false branch only)
 *   dhgrowth, dhsngrowth, dAgrowth — per-second changes
 *   iflice  — flooding contribution to ice growth (0 in B4a; B4b populates)
 *   hflatow, hfsenow, hflwrdout, hfswrow, hflwrow, hfradow — open-water radiative components
 *   subli   — sublimation rate (m water / s)
 *
 * Inputs:
 *   fsh, flo — incoming SW, LW (W/m²)
 *   Ta, qa   — air temperature (°C), specific humidity (kg/kg)
 *   rain, snow, runo — precipitation rain/snow + river runoff (m water / s)
 *   rsss     — reference sea surface salinity (used in use_virt_salt branch — B4b)
 *   ug, ustar — wind speed and friction velocity (m/s)
 *   T_oc, S_oc, H_ML — mixed-layer T (°C), S (PSU), depth (m)
 *   ice_dt   — ice timestep (s)
 *   ch, ce, ch_i, ce_i — transfer coefficients (open ocean / ice)
 *   evap_in — open-water evap from prior call (currently unused — Fortran also uses 0 default)
 *   lid_clo — lead-closing parameter h0 (or h0_s in southern hemisphere)
 *   geolon, geolat — geographic position (only used by Taylor/Briegleb albedo, not in scope)
 */
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
                     real_t *subli);

/*
 * Per-node sea-ice thermodynamics driver.
 * Mirror of Fortran thermodynamics at ice_thermo_oce.F90:148.
 *
 * Two-pass:
 *   1. Compute friction velocity ustar from (u_ice − u_w) and (v_ice − v_w)
 *      over owned nodes only (myDim), then halo-exchange ustar to fill eDim.
 *   2. For every node (myDim+eDim), call fesom_therm_ice. Cavity nodes are
 *      skipped via ulevels_nod2D > 1.
 *
 * Inputs read:
 *   ice->uice/vice/srfoce_*       — ice and ocean-surface state
 *   ice->thermo (constants)
 *   forcing->Ch_atm_oce/Ce_atm_oce/runoff
 *   jra->shortwave/longwave/Tair/shum/u_wind/v_wind/prec_rain/prec_snow
 *   sr->ref_sss + sr->ref_sss_local
 *
 * Outputs written:
 *   ice->data[0..2].values + values_old (a_ice, m_ice, m_snow)
 *   ice->thermo->t_skin/thdgr/thdgrsn/thdgra/thdgr_old/ustar
 *   ice->flx_fw, ice->flx_h
 */
void fesom_ice_thermodynamics(fesom_ice                     *ice,
                              struct fesom_partit           *partit,
                              struct fesom_mesh             *mesh,
                              const struct fesom_forcing    *forcing,
                              const struct fesom_jra55      *jra,
                              const struct fesom_sss_runoff *sr);

#endif /* FESOM_ICE_THERMO_H */
