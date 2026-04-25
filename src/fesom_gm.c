/*
 * GM / Redi state lifecycle and per-step physics.
 *
 * Phase G1 — state lifecycle (fesom_gm_alloc / fesom_gm_free).
 * Phase G2b — fesom_compute_sigma_xy / fesom_compute_neutral_slope.
 * Phase G3 — fesom_init_redi_gm  (deferred; abort stub).
 * Phase G4 — fesom_fer_solve_gamma / fesom_fer_gamma2vel (deferred).
 */
#include "fesom_gm.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { NL_MAX = 64 };          /* matches the cap used in fesom_eos / momentum */

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

/* ============================================================ */
/* G2b — compute_sigma_xy (oce_ale_pressure_bv.F90:2851)        */
/* ============================================================ */
/*
 * Density gradient on neutral surfaces — per node, per level, both
 * horizontal components (x, y).
 *
 * Recipe (per node n, per level nz):
 *   sigma_xy(comp, nz, n) = (-sw_alpha(nz,n) * dT/d{x,y}
 *                            +  sw_beta(nz,n) * dS/d{x,y}) * density_0
 *
 * Horizontal gradients of T and S are area-weighted means of the per-
 * element ∇T and ∇S over surrounding elements (gradient_sca is the
 * per-element ∂N_i/∂{x,y}).
 *
 * Halo audit: Fortran writes only for `n=1, myDim_nod2D` then halo-
 * exchanges sigma_xy (one component at a time, lines 2937-2942). We
 * mirror.
 *
 * Layout: sigma_xy[n * nl * 2 + nz * 2 + c]  (c = 0 → x, 1 → y).
 */
void fesom_compute_sigma_xy(struct fesom_aux *aux,
                            const struct fesom_tracers *tracers,
                            const struct fesom_mesh *mesh,
                            fesom_gm *gm,
                            struct fesom_partit *partit)
{
    const int  nl       = mesh->nl;
    const int  nl1      = nl - 1;             /* mid-layer count */
    const int  myDim    = mesh->myDim_nod2D;
    const real_t rho_ref = (real_t)FESOM_DENSITY_0;

    const real_t *T = tracers->data[FESOM_TRACER_T].values;
    const real_t *S = tracers->data[FESOM_TRACER_S].values;
    const real_t *alpha = aux->sw_alpha;     /* [N * nl] */
    const real_t *beta  = aux->sw_beta;      /* [N * nl] */
    real_t       *sxy   = gm->sigma_xy;      /* [N * nl * 2] */

    /* Stack scratch — Fortran allocates per-thread `tx[nl-1]` etc. */
    FESOM_CHECK(nl <= NL_MAX, "sigma_xy: nl %d > NL_MAX", nl);

    for (int n = 0; n < myDim; ++n) {
        int nle_node = mesh->nlevels_nod2D[n] - 1;     /* exclusive */
        int ule_node = mesh->ulevels_nod2D[n] - 1;     /* 0-based   */
        if (nle_node <= ule_node) continue;            /* dry */

        real_t tx[NL_MAX], ty[NL_MAX], sx[NL_MAX], sy[NL_MAX], vol[NL_MAX];
        for (int nz = ule_node; nz < nle_node; ++nz) {
            tx[nz] = ty[nz] = sx[nz] = sy[nz] = vol[nz] = 0.0;
        }

        int o0 = mesh->nod_in_elem2D_offsets[n];
        int o1 = mesh->nod_in_elem2D_offsets[n + 1];
        for (int k = o0; k < o1; ++k) {
            int el = mesh->nod_in_elem2D[k];
            int nle = mesh->nlevels[el] - 1;            /* exclusive  */
            int ule = mesh->ulevels[el] - 1;            /* 0-based    */
            const real_t *g = &mesh->gradient_sca[6 * el];   /* dN/dx, dN/dy */
            int v0 = mesh->elem_nodes[3*el + 0];
            int v1 = mesh->elem_nodes[3*el + 1];
            int v2 = mesh->elem_nodes[3*el + 2];
            real_t a = mesh->elem_area[el];

            for (int nz = ule; nz < nle; ++nz) {
                size_t i0 = (size_t)v0 * nl1 + nz;
                size_t i1 = (size_t)v1 * nl1 + nz;
                size_t i2 = (size_t)v2 * nl1 + nz;
                /* Note: Fortran loop bound is `ule_el..nle_el-1` per
                 * element; per-node accumulator only valid in
                 * [ule_node, nle_node-1] but that's a superset — for
                 * partial cells the per-element bound is correct. */
                if (nz < ule_node || nz >= nle_node) continue;

                vol[nz] += a;
                tx[nz] += (g[0]*T[i0] + g[1]*T[i1] + g[2]*T[i2]) * a;
                ty[nz] += (g[3]*T[i0] + g[4]*T[i1] + g[5]*T[i2]) * a;
                sx[nz] += (g[0]*S[i0] + g[1]*S[i1] + g[2]*S[i2]) * a;
                sy[nz] += (g[3]*S[i0] + g[4]*S[i1] + g[5]*S[i2]) * a;
            }
        }

        for (int nz = ule_node; nz < nle_node; ++nz) {
            real_t inv_vol = (vol[nz] > 0.0) ? 1.0 / vol[nz] : 0.0;
            real_t a = alpha[(size_t)n * nl + nz];
            real_t b = beta [(size_t)n * nl + nz];
            sxy[(size_t)n * nl * 2 + nz * 2 + 0] =
                (-a * tx[nz] + b * sx[nz]) * inv_vol * rho_ref;
            sxy[(size_t)n * nl * 2 + nz * 2 + 1] =
                (-a * ty[nz] + b * sy[nz]) * inv_vol * rho_ref;
        }
    }

    /* Halo exchange — Fortran lines 2937-2942 do it component-by-
     * component through an aux array. With our [n, nz, c] layout the
     * two components are interleaved, so a single 2-stride exchange
     * covers both. */
    if (partit && partit->npes > 1) {
        fesom_halo_exchange(sxy, FESOM_HALO_NOD2D, nl, 2, partit);
    }
}

/* ============================================================ */
/* G2b — compute_neutral_slope (oce_ale_pressure_bv.F90:2949)   */
/* ============================================================ */
/*
 * Active config:
 *   scaling_ODM95 = .true.  (c1 hyperbolic tanh tapering)
 *   scaling_LDD97 = .false. (c2 ≡ 1)
 *   Fer_GM = .true., Redi = .true., Redi_Ktaper = .true.
 * → fer_tapfac = c1 * c2 = c1;
 *   slope_tapered = neutral_slope * sqrt(c1*c2) = neutral_slope * sqrt(c1).
 *
 * Layout reminders:
 *   neutral_slope[n * (nl-1) * 3 + nz * 3 + comp]   (comp 0=x, 1=y, 2=|s|)
 *   slope_tapered[same]
 *   fer_tapfac[n * nl + nz]
 *
 * Halo: Fortran exchange_nod on neutral_slope and slope_tapered (lines
 * 3065-3066). We mirror.
 */
void fesom_compute_neutral_slope(struct fesom_aux *aux,
                                 const struct fesom_mesh *mesh,
                                 fesom_gm *gm,
                                 struct fesom_partit *partit)
{
    const int    nl       = mesh->nl;
    const int    nl1      = nl - 1;
    const int    myDim    = mesh->myDim_nod2D;
    const real_t g        = (real_t)FESOM_G;
    const real_t rho_ref  = (real_t)FESOM_DENSITY_0;
    const real_t eps      = 5.0e-6;
    const real_t eps_sq   = eps * eps;

    /* Active namelist params (active branches only). */
    const real_t ODM95_Scr = 0.2e-2;
    const real_t ODM95_Sd  = 1.0e-3;

    real_t *ns  = gm->neutral_slope;     /* [N * nl1 * 3] */
    real_t *st  = gm->slope_tapered;     /* [N * nl1 * 3] */
    real_t *tf  = gm->fer_tapfac;        /* [N * nl     ] */

    for (int n = 0; n < myDim; ++n) {
        int nle = mesh->nlevels_nod2D[n] - 1;       /* exclusive index for slopes */
        int ule = mesh->ulevels_nod2D[n] - 1;       /* 0-based */

        /* Fortran zeroes slope_tapered(:, :, n) = 0 for all nz at top of
         * the loop (line 2979). Mirror: zero the full nz range so dry
         * cells / gaps stay clean. */
        for (int nz = 0; nz < nl1; ++nz) {
            st[(size_t)n * nl1 * 3 + nz * 3 + 0] = 0.0;
            st[(size_t)n * nl1 * 3 + nz * 3 + 1] = 0.0;
            st[(size_t)n * nl1 * 3 + nz * 3 + 2] = 0.0;
        }

        if (nle <= ule) continue;

        /* Pass 1: neutral_slope from sigma_xy / N². Fortran lines 2997-3005. */
        for (int nz = ule; nz < nle; ++nz) {
            real_t bv_sum = aux->bvfreq[(size_t)n * nl + nz]
                          + aux->bvfreq[(size_t)n * nl + (nz + 1)];
            real_t denom  = (bv_sum > eps_sq) ? bv_sum : eps_sq;
            real_t ro_z_inv = 2.0 * g / rho_ref / denom;

            real_t sx = gm->sigma_xy[(size_t)n * nl * 2 + nz * 2 + 0] * ro_z_inv;
            real_t sy = gm->sigma_xy[(size_t)n * nl * 2 + nz * 2 + 1] * ro_z_inv;
            real_t sm = sqrt(sx*sx + sy*sy);
            ns[(size_t)n * nl1 * 3 + nz * 3 + 0] = sx;
            ns[(size_t)n * nl1 * 3 + nz * 3 + 1] = sy;
            ns[(size_t)n * nl1 * 3 + nz * 3 + 2] = sm;
        }

        /* Pass 2: ODM95 tapering c1. Fortran lines 3010-3016. */
        real_t c1[NL_MAX];
        for (int nz = 0; nz < nl1; ++nz) c1[nz] = 1.0;
        for (int nz = ule; nz < nle; ++nz) {
            real_t sm = ns[(size_t)n * nl1 * 3 + nz * 3 + 2];
            real_t v = 0.5 * (1.0 + tanh((ODM95_Scr - sm) / ODM95_Sd));
            real_t bv0 = aux->bvfreq[(size_t)n * nl + nz];
            real_t bv1 = aux->bvfreq[(size_t)n * nl + (nz + 1)];
            if (bv0 <= 0.0 || bv1 <= 0.0) v = 0.0;
            c1[nz] = v;
        }

        /* Pass 3: combine c1*c2 with c2≡1 in our config. Active branch
         * is `Fer_GM .and. Redi .and. Redi_Ktaper` (all true) → use
         * sqrt(c1*c2) for slope_tapered and store fer_tapfac=c1*c2.
         * Fortran lines 3049-3054. */
        for (int nz = ule; nz < nle; ++nz) {
            real_t cc = c1[nz];                          /* c1 * c2 = c1 */
            tf[(size_t)n * nl + nz] = cc;
            real_t s = sqrt(cc);
            real_t sx = ns[(size_t)n * nl1 * 3 + nz * 3 + 0];
            real_t sy = ns[(size_t)n * nl1 * 3 + nz * 3 + 1];
            real_t sm = ns[(size_t)n * nl1 * 3 + nz * 3 + 2];
            st[(size_t)n * nl1 * 3 + nz * 3 + 0] = sx * s;
            st[(size_t)n * nl1 * 3 + nz * 3 + 1] = sy * s;
            st[(size_t)n * nl1 * 3 + nz * 3 + 2] = sm * s;
        }
    }

    /* Halo exchange — Fortran lines 3065-3066. 3 components per (n,nz). */
    if (partit && partit->npes > 1) {
        fesom_halo_exchange(ns, FESOM_HALO_NOD2D, nl1, 3, partit);
        fesom_halo_exchange(st, FESOM_HALO_NOD2D, nl1, 3, partit);
        /* fer_tapfac is consumed in init_Redi_GM at myDim only — no
         * exchange needed (matches Fortran which doesn't exchange it). */
    }
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
