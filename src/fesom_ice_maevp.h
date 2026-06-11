#ifndef FESOM_ICE_MAEVP_H
#define FESOM_ICE_MAEVP_H

#include "fesom_ice_types.h"

struct fesom_mesh;
struct fesom_partit;

/*
 * Modified EVP (mEVP) rheology — Bouillon et al., Ocean Modelling 2013:
 * pseudo-time iteration with implicit alpha/beta stabilization replacing
 * the standard-EVP subcycling. Mirror of Fortran EVPdynamics_m
 * (ice_maEVP.F90:429-882) — the NR-optimized version with ssh2rhs,
 * stress_tensor_m and stress2rhs_m INLINED; those standalones are dead or
 * aEVP-only and are not ported.
 *
 * Scope per docs/plans/20260611-mevp-sea-ice-rheology.md:
 *   - selected by whichEVP == 1 (FESOM_WHICH_EVP=1); standard EVP
 *     (fesom_ice_evp.c) stays the default and is untouched.
 *   - levitating-ice ssh2rhs branch only (use_floatice=false — covers
 *     linfs AND zstar); cavity blocks compiled as written but no-op
 *     (use_cavity=false); aEVP (whichEVP=2) not ported.
 *   - NOTE the mEVP-vs-stdEVP asymmetries (plan fidelity-traps checklist):
 *     rdt = FULL ice_dt and the drag term carries rdt; no theta_io
 *     rotation; non-ice nodes are SKIPPED (velocity retained), not zeroed;
 *     uice_old/vice_old never touched; sigma persists across calls.
 *
 * Reads ice->uice_aux/vice_aux + ice->work.mevp_* scratch (allocated by
 * fesom_ice_init only when whichEVP != 0) and mesh->bc_index_nod2D.
 */
void fesom_ice_evp_dynamics_m(fesom_ice            *ice,
                              struct fesom_partit  *partit,
                              struct fesom_mesh    *mesh);

#endif /* FESOM_ICE_MAEVP_H */
