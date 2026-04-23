#ifndef FESOM_TRACERS_H
#define FESOM_TRACERS_H

#include "fesom_types.h"

struct fesom_mesh;

/*
 * Mirror of Fortran t_tracer_data + t_tracer (MOD_TRACER.F90).
 *
 * Phase 1 carries 2 tracers: T (index 0), S (index 1). Adams-Bashforth
 * order 2 → values + valuesAB only (no valuesold).
 *
 *   data[i].values   [nod2D * nl]   instant T or S
 *   data[i].valuesAB [nod2D * nl]   AB2 history for tracer i
 *   del_ttf          [nod2D * nl]   diffusion + ALE-reconstruction accumulator,
 *                                   zeroed at the start of each timestep
 *                                   (FRESH_START.md §5)
 *
 * Note: per FRESH_START.md §14.1, horizontal/explicit-vertical diffusion
 * accumulate into del_ttf; the implicit TDMA later operates on the
 * ALE-reconstructed tracer. Do NOT apply diffusion directly to values.
 */
typedef struct fesom_tracer_data {
    real_t *values;
    real_t *valuesAB;
    real_t *valuesold;   /* previous-timestep values; AB2 single slot
                            (Fortran: tracers%data(k)%valuesold(1,:,:) only) */
} fesom_tracer_data;

#define FESOM_NUM_TRACERS 2
#define FESOM_TRACER_T 0
#define FESOM_TRACER_S 1

typedef struct fesom_tracers {
    int                num_tracers;
    fesom_tracer_data  data[FESOM_NUM_TRACERS];
    real_t            *del_ttf;
} fesom_tracers;

void fesom_tracers_alloc(fesom_tracers *t, const struct fesom_mesh *mesh);
void fesom_tracers_free (fesom_tracers *t);

#endif /* FESOM_TRACERS_H */
