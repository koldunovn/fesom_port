#ifndef FESOM_TRACER_ADV_H
#define FESOM_TRACER_ADV_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;
struct fesom_tracers;

/*
 * Phase 1 simple-upwind tracer advection. Literal port of:
 *   init_tracers_AB                (oce_tracer_mod.F90:13-145, AB2 path)
 *   adv_tra_hor_upw1               (oce_adv_tra_hor.F90:64-257)
 *   adv_tra_ver_upw1               (oce_adv_tra_ver.F90:244-328)
 *   oce_tra_adv_flux2dtracer       (oce_adv_tra_driver.F90:423-575, non-FCT path)
 *   ALE reconstruction             (diff_tracers_ale lines 481-484)
 *   valuesold update               (oce_tracer_mod.F90:100-102, AB2)
 *
 * Deferred to Phase 2: FCT (high-order MUSCL/MFCT + Zalesak limiter),
 * tracer_gradient_elements, fill_up_dn_grad.
 *
 * Per-tracer scratch:
 *   adv_flux_hor   [edge2D * (nl-1)]
 *   adv_flux_ver   [nod2D  * nl    ]
 *   del_ttf_advhoriz, del_ttf_advvert  [nod2D * nl]
 *
 * fesom_tracer_adv_init must be called once after the mesh is read; it
 * allocates the scratch arrays. fesom_tracer_adv_free frees them.
 */
typedef struct fesom_tracer_adv_scratch {
    real_t *adv_flux_hor;        /* [edge2D * (nl-1)] */
    real_t *adv_flux_ver;        /* [nod2D  * nl]     */
    real_t *del_ttf_advhoriz;    /* [nod2D  * (nl-1)] (we allocate nl for symmetry) */
    real_t *del_ttf_advvert;     /* [nod2D  * (nl-1)] */
} fesom_tracer_adv_scratch;

void fesom_tracer_adv_init(fesom_tracer_adv_scratch *sc,
                           const struct fesom_mesh  *mesh);
void fesom_tracer_adv_free(fesom_tracer_adv_scratch *sc);

/*
 * Run one full upwind-advection step for one tracer:
 *   1. init_tracers_AB (zero del_ttf*, fill valuesAB, save valuesold)
 *   2. compute horizontal upwind flux into adv_flux_hor
 *   3. compute vertical upwind flux into adv_flux_ver
 *   4. accumulate fluxes into del_ttf_advhoriz, del_ttf_advvert
 *   5. del_ttf := del_ttf_advhoriz + del_ttf_advvert
 *   6. ALE reconstruction:  T_new = (T*hnode + del_ttf) / hnode_new
 *   7. update valuesold = (the values BEFORE ALE; matches Fortran ordering at
 *      line 100-102 inside init_tracers_AB which runs at start of next step
 *      but uses the value snapshot from before adv on the current step)
 *
 * Caller chooses tr_idx ∈ {FESOM_TRACER_T, FESOM_TRACER_S}.
 */
void fesom_tracer_advect_one(fesom_tracer_adv_scratch *sc,
                             int                       tr_idx,
                             const struct fesom_mesh  *mesh,
                             const struct fesom_dyn   *dyn,
                             struct fesom_tracers     *tracers);

#endif /* FESOM_TRACER_ADV_H */
