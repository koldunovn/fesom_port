#ifndef FESOM_GM_H
#define FESOM_GM_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_partit;
struct fesom_aux;
struct fesom_tracers;
struct fesom_dyn;

/*
 * GM / Redi state. Holds the auxiliary fields the parameterisation needs.
 * Mirror of the relevant globals in Fortran's o_ARRAYS module
 * (oce_modules.F90:296+) plus per-step scratch.
 *
 * Memory layout (flat, row-major; Fortran is column-major (comp, nz, n)):
 *   sigma_xy        [N * nl * 2]      (node, level, comp)
 *   neutral_slope   [N * (nl-1) * 3]  (node, level, {x, y, |s|})
 *   slope_tapered   [N * (nl-1) * 3]  (same as above)
 *   fer_tapfac      [N * nl]          (node, level)
 *   fer_gamma       [N * nl * 2]      (node, level, comp)
 *   fer_K           [N * nl]          (node, level) — GM thickness diffusivity
 *   Ki              [N * nl]          (node, level) — Redi diffusivity
 *   fer_C           [N]
 *   fer_scal        [N]
 *
 * All per-node arrays are sized myDim+eDim. Halo coverage is required
 * because downstream readers (init_Redi_GM, fer_solve_Gamma, fer_gamma2vel,
 * tracer-side Redi terms) all loop into the halo.
 *
 * fer_K is initialised to K_GM_max=1000 m²/s at allocation, mirroring
 * Fortran oce_setup_step.F90:934.
 */
typedef struct fesom_gm {
    real_t *sigma_xy;
    real_t *neutral_slope;
    real_t *slope_tapered;
    real_t *fer_tapfac;
    real_t *fer_gamma;
    real_t *fer_K;
    real_t *Ki;
    real_t *fer_C;
    real_t *fer_scal;
    /* G7 scratch — per-element horizontal gradient of the active tracer
     * (Fortran tr_xy from oce_tracer_mod.F90:178-181). One tracer at a
     * time; fesom_diff_ver_part_redi_expl rebuilds and halo-exchanges.
     * Layout: [E_full * (nl-1) * 2], comp innermost. */
    real_t *tr_xy;
    /* G7b scratch — per-node 3D vertical gradient of the active tracer
     * (Fortran tr_z from oce_tracer_mod.F90:222-227). Sized [N * nl]; nz=0
     * (top) and nz=nlevels-1 (bottom) are 0 boundaries. */
    real_t *tr_z;
} fesom_gm;

void fesom_gm_alloc(fesom_gm *g, const struct fesom_mesh *mesh);
void fesom_gm_free (fesom_gm *g);

/* --- Phase G2b: density gradients + neutral slope (per step) ----------- */
/* Forward declarations — implementations land in Phase G2b. */
void fesom_compute_sigma_xy      (struct fesom_aux         *aux,
                                  const struct fesom_tracers *tracers,
                                  const struct fesom_mesh  *mesh,
                                  fesom_gm                 *gm,
                                  struct fesom_partit      *partit);

void fesom_compute_neutral_slope (struct fesom_aux         *aux,
                                  const struct fesom_mesh  *mesh,
                                  fesom_gm                 *gm,
                                  struct fesom_partit      *partit);

/* --- Phase G3: per-step GM/Redi coefficient builder -------------------- */
void fesom_init_redi_gm          (struct fesom_aux         *aux,
                                  const struct fesom_mesh  *mesh,
                                  fesom_gm                 *gm,
                                  struct fesom_partit      *partit);

/* --- Phase G4: streamfunction solve and bolus velocity reconstruction -- */
void fesom_fer_solve_gamma       (const struct fesom_aux   *aux,
                                  const struct fesom_mesh  *mesh,
                                  fesom_gm                 *gm,
                                  struct fesom_partit      *partit);

void fesom_fer_gamma2vel         (struct fesom_dyn         *dyn,
                                  const struct fesom_mesh  *mesh,
                                  const fesom_gm           *gm,
                                  struct fesom_partit      *partit);

/* --- Phase G7a: vertical-explicit Redi (oce_ale_tracer.F90:1086-1169)
 *
 * Adds the off-diagonal Redi tensor's vertical projection to T values
 * for the given tracer. Reads valuesold (= pre-step T, matches Fortran
 * `values` at tr_xy build time). Builds gm->tr_xy as scratch and
 * halo-exchanges it. No-op if gm == NULL. */
void fesom_diff_ver_part_redi_expl(int                          tr_idx,
                                   fesom_gm                    *gm,
                                   const struct fesom_aux      *aux,
                                   const struct fesom_mesh     *mesh,
                                   struct fesom_tracers        *tracers,
                                   struct fesom_partit         *partit);

/* --- Phase G7b: horizontal Redi (oce_ale_tracer.F90:1173-1336)
 *
 * Edge-based horizontal Redi/eddy diffusivity flux with 5 partial-cell
 * branches A/B/C/D/E for level mismatches between the two adjacent
 * elements of an edge. Builds gm->tr_z (per-node vertical gradient)
 * and reuses gm->tr_xy from a preceding fesom_diff_ver_part_redi_expl
 * call (same tracer). No-op if gm == NULL.
 *
 * Must be called AFTER fesom_diff_ver_part_redi_expl in the same
 * per-tracer iteration so gm->tr_xy is populated. */
void fesom_diff_part_hor_redi    (int                          tr_idx,
                                  fesom_gm                    *gm,
                                  const struct fesom_aux      *aux,
                                  const struct fesom_mesh     *mesh,
                                  struct fesom_tracers        *tracers,
                                  struct fesom_partit         *partit);

#endif /* FESOM_GM_H */
