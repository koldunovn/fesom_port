#ifndef FESOM_ICE_H
#define FESOM_ICE_H

#include "fesom_ice_types.h"

struct fesom_mesh;
struct fesom_partit;

/*
 * Top-level lifecycle for the sea-ice subsystem.
 *
 * Mirrors Fortran ice_init (MOD_ICE.F90:572) and ice_setup
 * (ice_setup_step.F90:51). Init does the allocation and namelist-style
 * defaults; setup wires in the time-step constants (Tevp_inv) and the
 * one-time mass matrix.
 *
 * Call order from fesom_main:
 *     fesom_ice_init (ice, partit, mesh);
 *     fesom_ice_setup(ice, partit, mesh);   // Phase A2b — Tevp_inv, mass matrix stub
 *     ... model loop calls fesom_ice_step every step ...
 *     fesom_ice_free (ice);
 */

void fesom_ice_init(fesom_ice           *ice,
                    struct fesom_partit *partit,
                    struct fesom_mesh   *mesh);

/* One-time post-init: derives ice_dt/Tevp_inv/Clim_evp from the ocean dt
 * (fesom_phase1_dt) and, in Phase E, calls fesom_ice_mass_matrix_fill.
 * Mirrors ice_setup_step.F90:51-91. */
void fesom_ice_setup(fesom_ice           *ice,
                     struct fesom_partit *partit,
                     struct fesom_mesh   *mesh);

/* Cold-start ice IC: where surface T < 0°C, set ice/snow thicknesses and
 * concentration. Mirror of ice_setup_step.F90:500-521 (cold-start branch;
 * "initialize from file" branch out of scope).
 *
 * Per Fortran:
 *   if SST < 0:
 *     NH (geo lat > 0):  m_ice=1.0  m_snow=0.1
 *     SH:                m_ice=2.0  m_snow=0.5
 *     a_ice=0.9, u_ice=v_ice=0
 *   else: zero (open water)
 *
 * Must run AFTER tracer IC (PHC or constant) so surface T is realistic. */
struct fesom_tracers;
void fesom_ice_initial_state(fesom_ice                  *ice,
                             const struct fesom_tracers *tracers,
                             struct fesom_partit        *partit,
                             struct fesom_mesh          *mesh);

/* Per-step driver. Mirrors ice_setup_step.F90:96 (ice_timestep).
 *   FESOM_NO_ICE_DYN    skips EVP dynamics    (Phase D — currently stubbed)
 *   FESOM_NO_ICE_ADV    skips FCT advection   (Phase E — currently stubbed)
 *   FESOM_NO_ICE_THERMO skips thermodynamics  (Phase B — wired)
 *
 * Forcing/jra/sr may be NULL ONLY if FESOM_NO_ICE_THERMO is set (the
 * thermodynamics path needs all three). dyn/tracers are required to populate
 * srfoce_* (the simplified ocean2ice copy; full ocean2ice lands in Phase C). */
struct fesom_dyn;
struct fesom_tracers;
struct fesom_forcing;
struct fesom_jra55;
struct fesom_sss_runoff;
struct fesom_ssh_stiff;

void fesom_ice_step(int                            step,
                    fesom_ice                     *ice,
                    struct fesom_partit           *partit,
                    struct fesom_mesh             *mesh,
                    const struct fesom_dyn        *dyn,
                    const struct fesom_tracers    *tracers,
                    struct fesom_forcing          *forcing,
                    const struct fesom_jra55      *jra,
                    const struct fesom_sss_runoff *sr,
                    const struct fesom_ssh_stiff  *stiff);

void fesom_ice_free(fesom_ice *ice);

#endif /* FESOM_ICE_H */
