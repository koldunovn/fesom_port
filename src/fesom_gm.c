/*
 * GM / Redi state lifecycle.
 *
 * The physics functions (compute_sigma_xy, compute_neutral_slope,
 * init_Redi_GM, fer_solve_Gamma, fer_gamma2vel) are declared in
 * fesom_gm.h but implemented in later phases (G2b, G3, G4). Phase G1
 * lands the storage only.
 */
#include "fesom_gm.h"
#include "fesom_constants.h"
#include "fesom_mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GM_K_GM_MAX 1000.0   /* matches namelist.oce K_GM_max default
                              * and oce_setup_step.F90:934 fer_K init   */

#define GM_CHECK(p, what) do {                                       \
    if (!(p)) {                                                      \
        fprintf(stderr, "fesom_gm_alloc: out of memory (%s)\n",      \
                (what));                                             \
        abort();                                                     \
    }                                                                \
} while (0)

void fesom_gm_alloc(fesom_gm *g, const struct fesom_mesh *mesh)
{
    memset(g, 0, sizeof(*g));

    int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int nl = mesh->nl;

    size_t n_nl   = (size_t)N * (size_t)nl;
    size_t n_nlx2 = n_nl * 2;
    size_t n_nlm1 = (size_t)N * (size_t)(nl - 1);
    size_t n_nlm1x3 = n_nlm1 * 3;

    g->sigma_xy      = calloc(n_nlx2,    sizeof(real_t));
    g->neutral_slope = calloc(n_nlm1x3,  sizeof(real_t));
    g->slope_tapered = calloc(n_nlm1x3,  sizeof(real_t));
    g->fer_tapfac    = calloc(n_nl,      sizeof(real_t));
    g->fer_gamma     = calloc(n_nlx2,    sizeof(real_t));
    g->fer_K         = calloc(n_nl,      sizeof(real_t));
    g->Ki            = calloc(n_nl,      sizeof(real_t));
    g->fer_C         = calloc((size_t)N, sizeof(real_t));
    g->fer_scal      = calloc((size_t)N, sizeof(real_t));

    GM_CHECK(g->sigma_xy,      "sigma_xy");
    GM_CHECK(g->neutral_slope, "neutral_slope");
    GM_CHECK(g->slope_tapered, "slope_tapered");
    GM_CHECK(g->fer_tapfac,    "fer_tapfac");
    GM_CHECK(g->fer_gamma,     "fer_gamma");
    GM_CHECK(g->fer_K,         "fer_K");
    GM_CHECK(g->Ki,            "Ki");
    GM_CHECK(g->fer_C,         "fer_C");
    GM_CHECK(g->fer_scal,      "fer_scal");

    /* Initial value matches Fortran oce_setup_step.F90:934.
     * init_Redi_GM (G3) overwrites this every step, but pre-G3 readers
     * (none yet) should see the namelist default rather than 0. */
    for (size_t i = 0; i < n_nl; ++i) {
        g->fer_K[i] = GM_K_GM_MAX;
    }
}

void fesom_gm_free(fesom_gm *g)
{
    free(g->sigma_xy);
    free(g->neutral_slope);
    free(g->slope_tapered);
    free(g->fer_tapfac);
    free(g->fer_gamma);
    free(g->fer_K);
    free(g->Ki);
    free(g->fer_C);
    free(g->fer_scal);
    memset(g, 0, sizeof(*g));
}

/* Physics-function stubs land in later phases. Each stub aborts so that
 * a premature wire-up is caught loudly rather than silently producing
 * wrong physics. */

void fesom_compute_sigma_xy(struct fesom_aux *aux,
                            const struct fesom_tracers *tracers,
                            const struct fesom_mesh *mesh,
                            fesom_gm *gm,
                            struct fesom_partit *partit)
{
    (void)aux; (void)tracers; (void)mesh; (void)gm; (void)partit;
    fprintf(stderr, "fesom_compute_sigma_xy: not yet ported (Phase G2b)\n");
    abort();
}

void fesom_compute_neutral_slope(struct fesom_aux *aux,
                                 const struct fesom_mesh *mesh,
                                 fesom_gm *gm,
                                 struct fesom_partit *partit)
{
    (void)aux; (void)mesh; (void)gm; (void)partit;
    fprintf(stderr, "fesom_compute_neutral_slope: not yet ported (Phase G2b)\n");
    abort();
}

void fesom_init_redi_gm(struct fesom_aux *aux,
                        const struct fesom_mesh *mesh,
                        fesom_gm *gm,
                        struct fesom_partit *partit)
{
    (void)aux; (void)mesh; (void)gm; (void)partit;
    fprintf(stderr, "fesom_init_redi_gm: not yet ported (Phase G3)\n");
    abort();
}

void fesom_fer_solve_gamma(const struct fesom_aux *aux,
                           const struct fesom_mesh *mesh,
                           fesom_gm *gm,
                           struct fesom_partit *partit)
{
    (void)aux; (void)mesh; (void)gm; (void)partit;
    fprintf(stderr, "fesom_fer_solve_gamma: not yet ported (Phase G4)\n");
    abort();
}

void fesom_fer_gamma2vel(struct fesom_dyn *dyn,
                         const struct fesom_mesh *mesh,
                         const fesom_gm *gm,
                         struct fesom_partit *partit)
{
    (void)dyn; (void)mesh; (void)gm; (void)partit;
    fprintf(stderr, "fesom_fer_gamma2vel: not yet ported (Phase G4)\n");
    abort();
}
