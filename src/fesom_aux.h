#ifndef FESOM_AUX_H
#define FESOM_AUX_H

#include "fesom_types.h"

struct fesom_mesh;

/*
 * Phase 1 subset of the Fortran o_ARRAYS module globals (oce_modules.F90).
 *
 * Per-node 3D fields:
 *   density_m_rho0  [nod2D * nl]   in-situ density minus density_ref (kg/m³)
 *   hpressure       [nod2D * nl]   hydrostatic pressure (Pa) — written by
 *                                  pressure_bv (linfs) or pressure_force (zstar)
 *   bvfreq          [nod2D * nl]   N² (s⁻²)
 *   sw_alpha        [nod2D * nl]   thermal expansion coeff (1/K)
 *   sw_beta         [nod2D * nl]   haline contraction coeff (1/(g/kg))
 *   Kv              [nod2D * nl]   vertical tracer diffusivity (m²/s)
 *
 * Per-element 3D fields:
 *   Av              [elem2D * nl]  vertical momentum viscosity (m²/s)
 *   pgf_x           [elem2D * nl]  pressure-gradient force, x-component (m/s²)
 *   pgf_y           [elem2D * nl]  pressure-gradient force, y-component (m/s²)
 *
 * KPP additions (Kv per tracer, blmc, ghats, OBL depth) are deferred to Phase 4.
 */
typedef struct fesom_aux {
    real_t *density_m_rho0;
    real_t *hpressure;
    real_t *bvfreq;
    real_t *sw_alpha;
    real_t *sw_beta;
    real_t *Kv;

    real_t *Av;
    real_t *pgf_x;
    real_t *pgf_y;
} fesom_aux;

void fesom_aux_alloc(fesom_aux *a, const struct fesom_mesh *mesh);
void fesom_aux_free (fesom_aux *a);

#endif /* FESOM_AUX_H */
