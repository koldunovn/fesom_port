#ifndef FESOM_ICE_EVP_H
#define FESOM_ICE_EVP_H

#include "fesom_ice_types.h"

struct fesom_mesh;
struct fesom_partit;

/*
 * Standard EVP rheology (Hibler/Hunke). Mirror of Fortran ice_EVP.F90.
 * Scope per docs/plans/20260425-sea-ice-port.md:
 *   - whichEVP = 0 only (this file). mEVP/aEVP (ice_maEVP.F90) are
 *     out of scope; the dispatcher in fesom_ice_step aborts otherwise.
 *
 * Three subroutines:
 *   stress_tensor — element loop: compute strain-rate tensor + update sigma
 *   stress2rhs    — divergence of sigma scattered to nodes (uice_rhs/vice_rhs)
 *   EVPdynamics   — 120-subcycle solver wrapping stress_tensor + stress2rhs
 *                   plus velocity update + halo exchanges
 */

/* Compute strain rates eps11/12/22 from u_ice/v_ice over each element,
 * then update sigma11/12/22 via the EVP rheology (Hibler 1979 with the
 * Hunke EVP modification). Element loop bound is myDim_elem2D (matches
 * Fortran ice_EVP.F90:94). Ice-free elements (ice_strength≤0) and cavity
 * elements (ulevels>1) are skipped. */
void fesom_ice_stress_tensor(fesom_ice                *ice,
                             struct fesom_partit      *partit,
                             struct fesom_mesh        *mesh);

/* Divergence of sigma scattered to nodes — produces uice_rhs/vice_rhs.
 * Mirror of stress2rhs at ice_EVP.F90:177. */
void fesom_ice_stress2rhs(fesom_ice            *ice,
                          struct fesom_partit  *partit,
                          struct fesom_mesh    *mesh);

/* Subcycled EVP solver (120 iterations of stress_tensor + stress2rhs
 * + velocity update + halo exchange of uice/vice).
 * Mirror of EVPdynamics at ice_EVP.F90:330. */
void fesom_ice_evp_dynamics(fesom_ice            *ice,
                            struct fesom_partit  *partit,
                            struct fesom_mesh    *mesh);

#endif /* FESOM_ICE_EVP_H */
