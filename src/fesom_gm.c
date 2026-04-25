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
#include "fesom_dyn.h"
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
    const int  myDim    = mesh->myDim_nod2D;
    const real_t rho_ref = (real_t)FESOM_DENSITY_0;

    const real_t *T = tracers->data[FESOM_TRACER_T].values;   /* [N * nl] stride */
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
                /* Tracer T/S are stored with stride nl (per fesom_tracers_alloc) — NOT nl-1. */
                size_t i0 = (size_t)v0 * nl + nz;
                size_t i1 = (size_t)v1 * nl + nz;
                size_t i2 = (size_t)v2 * nl + nz;
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

/* ============================================================ */
/* G3 — init_Redi_GM (oce_fer_gm.F90:215-506)                   */
/* ============================================================ */
/*
 * Per-step coefficient builder for GM (fer_K, fer_C, fer_scal) and
 * Redi (Ki). Three logical phases:
 *
 *   Pass 1 (F1 = horizontal scaling) — per-node loop using the
 *     conservative level extrema nlevels_nod2D_min / ulevels_nod2D_max
 *     (which are exactly what Fortran derives inline at lines 252-257).
 *     Computes fer_scal, fer_K[nzmin,n], fer_C[n] (all GM); Ki[nzmin,n]
 *     is filled here too but immediately overwritten in the next pass
 *     since Fer_GM ∧ Redi is true in our config.
 *
 *   Pass 2 (Redi=GM sync) — Ki[nzmin,n] = max(fer_scal·Redi_Kmax, K_GM_min)
 *     (Fortran lines 351-359; the FESOM 1.4 carry-over).
 *
 *   Pass 3 (F2 = vertical scaling) — per-node loop using regular
 *     per-node nlevels_nod2D / ulevels_nod2D. Applies the
 *     scaling_GMzexp depth-exp scaling to both fer_K and Ki. Redi
 *     also gets the Redi_Ktaper sqrt(fer_tapfac) tapering.
 *
 * Active sub-features only (default namelist.oce):
 *   Fer_GM=T, Redi=T, K_GM_resscalorder=2 (real, not int → 1/2=0.5),
 *   scaling_resolution=T, scaling_GMzexp=T, Redi_Ktaper=T,
 *   K_GM_max=1000, K_GM_min=2, Redi_Kmin=100, K_GM_cmin=0.1, K_GM_cm=3,
 *   GMzexp_zref=500, GMzexp_smin=0.6, refscalresol=100km.
 *   Skipped: scaling_Rossby, FESOM14, GINsea, Ferreira (and so the
 *   K_GM_bvref=1 MLD-reference branch never runs); also K_GM_Ktaper,
 *   K_GM_rampmax/min, scaling_LDD97. Redi_Kmax=0 auto-syncs to K_GM_max.
 *
 * Halo: Fortran exchanges fer_c, fer_K, Ki (lines 503-505). Mirror.
 */
void fesom_init_redi_gm(struct fesom_aux *aux,
                        const struct fesom_mesh *mesh,
                        fesom_gm *gm,
                        struct fesom_partit *partit)
{
    (void)aux;   /* unused in active branches (would be needed for Ferreira/MLD ref) */

    const int    nl       = mesh->nl;
    const int    myDim    = mesh->myDim_nod2D;
    const real_t pi       = 3.14159265358979323846;

    /* Active-config namelist constants. */
    const real_t K_GM_max     = 1000.0;
    const real_t K_GM_min     = 2.0;
    const real_t K_GM_cmin    = 0.1;
    const real_t K_GM_cm      = 3.0;
    const real_t Redi_Kmin    = 100.0;
    const real_t Redi_Kmax    = K_GM_max;     /* auto-sync (Fortran 240-242) */
    const real_t GMzexp_zref  = 500.0;
    const real_t GMzexp_smin  = 0.6;
    const real_t refscalresol = 100000.0;     /* 100 km */
    const real_t inv_refscalresol_sq = 1.0 / (refscalresol * refscalresol);

    /* --- Pass 1: F1(x,y) per-node scalar ------------------------- */
    /* Bounds nzmin/nzmax mirror Fortran's inline min(nlevels)/max(ulevels)
     * over surrounding elements (= our G0 mesh prereqs). 0-based. */
    for (int n = 0; n < myDim; ++n) {
        int nzmax = mesh->nlevels_nod2D_min[n] - 1;   /* exclusive */
        int nzmin = mesh->ulevels_nod2D_max[n] - 1;   /* 0-based */

        /* baroclinic wave speed cm — depth integral of |N|. */
        real_t cm_sum = 0.0;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t bv0 = aux->bvfreq[(size_t)n * nl + nz];
            real_t bv1 = aux->bvfreq[(size_t)n * nl + (nz + 1)];
            real_t v0 = (bv0 > 0.0) ? sqrt(bv0) : 0.0;
            real_t v1 = (bv1 > 0.0) ? sqrt(bv1) : 0.0;
            cm_sum += mesh->hnode_new[(size_t)n * nl + nz] * 0.5 * (v0 + v1);
        }
        /* Limited from below by K_GM_cmin, and divided by π·K_GM_cm
         * (Fortran 268-269). c1 itself is only used by scaling_Rossby
         * which is off — we still compute cm. */
        real_t cm = cm_sum / pi / K_GM_cm;
        if (cm < K_GM_cmin) cm = K_GM_cmin;

        /* scaling — start at 1, no Rossby-radius cutoff, area-based
         * resolution scaling for K_GM_resscalorder=2. */
        real_t area_n = mesh->area[(size_t)n * nl + 0];   /* surface CV area */
        real_t scaling = sqrt(area_n * inv_refscalresol_sq * 2.0);
        if (scaling > 1.0) scaling = 1.0;          /* Fortran: fer_scal=min(scaling,1) */

        gm->fer_scal[n] = scaling;

        /* fer_K[nzmin,n] = max(fer_scal*K_GM_max, K_GM_min). Fortran 328-330.
         * Stored 0-based — `nzmin` is the 0-based level index. */
        real_t k_top = scaling * K_GM_max;
        if (k_top < K_GM_min) k_top = K_GM_min;
        gm->fer_K[(size_t)n * nl + nzmin] = k_top;
        gm->fer_C[n] = cm * cm;

        /* Redi placeholder (always overwritten in pass 2 below) — the
         * Fortran K_hor**x branch from lines 340-344 is dead in our
         * config because Fer_GM ∧ Redi triggers the override below. */

        /* --- Pass 2 (inline): Redi = GM sync (Fortran 351-359) ---
         * Ki[nzmin,n] = max(fer_scal*Redi_Kmax, K_GM_min). */
        real_t ki_top = scaling * Redi_Kmax;
        if (ki_top < K_GM_min) ki_top = K_GM_min;
        gm->Ki[(size_t)n * nl + nzmin] = ki_top;
    }

    /* --- Pass 3: F2(z) per-node, per-level vertical scaling ------ */
    real_t zscaling[NL_MAX];
    for (int n = 0; n < myDim; ++n) {
        int nzmax = mesh->nlevels_nod2D[n] - 1;       /* exclusive */
        int nzmin = mesh->ulevels_nod2D[n] - 1;       /* 0-based */

        /* zscaling: GMzexp_smin + (1-GMzexp_smin)*exp(-|z|/zref);
         * clamp to [GMzexp_smin, 1]. Fortran 437-441. */
        for (int nz = 0; nz < nl; ++nz) zscaling[nz] = 1.0;
        for (int nz = nzmin; nz <= nzmax; ++nz) {
            real_t z = mesh->zbar_3d_n[(size_t)n * nl + nz];
            real_t v = GMzexp_smin
                     + (1.0 - GMzexp_smin) * exp(-fabs(z) / GMzexp_zref);
            if (v > 1.0)          v = 1.0;
            if (v < GMzexp_smin)  v = GMzexp_smin;
            zscaling[nz] = v;
        }

        /* Apply F2 to fer_K (Fer_GM). Fortran 459-464. */
        real_t k_top = gm->fer_K[(size_t)n * nl + nzmin];
        for (int nz = nzmin + 1; nz <= nzmax; ++nz) {
            gm->fer_K[(size_t)n * nl + nz] = k_top * zscaling[nz];
        }
        gm->fer_K[(size_t)n * nl + nzmin] = k_top * zscaling[nzmin];
        /* K_GM_Ktaper = false → skip the GM tapering branch. */

        /* Apply F2 to Ki (Redi). Fortran 482-487. */
        real_t ki_top = gm->Ki[(size_t)n * nl + nzmin];
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            gm->Ki[(size_t)n * nl + nz] =
                ki_top * 0.5 * (zscaling[nz] + zscaling[nz + 1]);
        }
        gm->Ki[(size_t)n * nl + nzmin] =
            ki_top * 0.5 * (zscaling[nzmin] + zscaling[nzmin + 1]);

        /* Redi_Ktaper = true (Fortran 490-497):
         *   Ki[nz,n] = Ki[nz,n] * sqrt(tapfac) + Redi_Kmin * |sqrt(tapfac)-1|
         */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t tf = gm->fer_tapfac[(size_t)n * nl + nz];
            real_t s  = sqrt(tf);
            real_t k  = gm->Ki[(size_t)n * nl + nz];
            gm->Ki[(size_t)n * nl + nz] = k * s + Redi_Kmin * fabs(s - 1.0);
        }
    }

    /* --- Halo exchange (Fortran lines 503-505) -------------------- */
    if (partit && partit->npes > 1) {
        fesom_exchange_nod2D(gm->fer_C, partit);
        fesom_halo_exchange(gm->fer_K, FESOM_HALO_NOD2D, nl, 1, partit);
        fesom_halo_exchange(gm->Ki,    FESOM_HALO_NOD2D, nl, 1, partit);
    }
}

/* ============================================================ */
/* G4 — fer_solve_Gamma (oce_fer_gm.F90:40-163)                 */
/* ============================================================ */
/*
 * Per-node 1D TDMA solving the GM streamfunction equation.
 *
 *   ∂z(C·∂z Γ) − N²·Γ = (g/ρ₀)·∇σ·K_GM(z)
 *
 * The two horizontal components (x, y) share the same matrix and are
 * solved together (Fortran does it with a (2, nl) RHS; we use a
 * 2-component sweep per node).
 *
 * Note the two distinct level bounds:
 *   - outer (zbar_n / Z_n setup):  nlevels_nod2D / ulevels_nod2D
 *   - inner (TDMA solve):          nlevels_nod2D_min / ulevels_nod2D_max
 *                                  (= conservative bounds)
 * In our linfs / no-cavity / no-partial-cells config
 * nzmin_outer == nzmin_inner = 0; nzmax_outer ≥ nzmax_inner.
 *
 * Halo: exchange_nod(fer_gamma) — Fortran line 162. We use stride 2
 * to cover both components in a single exchange.
 */
void fesom_fer_solve_gamma(const struct fesom_aux *aux,
                           const struct fesom_mesh *mesh,
                           fesom_gm *gm,
                           struct fesom_partit *partit)
{
    const int    nl       = mesh->nl;
    const int    myDim    = mesh->myDim_nod2D;
    const real_t g        = (real_t)FESOM_G;
    const real_t rho_ref  = (real_t)FESOM_DENSITY_0;

    real_t zbar_n[NL_MAX], Z_n[NL_MAX];
    real_t a[NL_MAX], b[NL_MAX], c[NL_MAX];
    real_t cp[NL_MAX], tp_x[NL_MAX], tp_y[NL_MAX];
    real_t tr_x[NL_MAX], tr_y[NL_MAX];
    FESOM_CHECK(nl <= NL_MAX, "fer_solve_gamma: nl %d > NL_MAX", nl);

    for (int n = 0; n < myDim; ++n) {
        /* Outer bounds (Fortran lines 71-72). 0-based throughout. */
        int nzmax_o = mesh->nlevels_nod2D[n] - 1;     /* index of bottom interface */
        int nzmin_o = mesh->ulevels_nod2D[n] - 1;     /* = 0 in our config       */
        if (nzmax_o <= nzmin_o + 1) continue;          /* degenerate column */

        /* Build zbar_n, Z_n on [nzmin_o, nzmax_o] (Fortran 76-86).
         * For linfs/no-cavity: zbar_n[nz] == mesh->zbar[nz]. */
        for (int nz = 0; nz < nl; ++nz) { zbar_n[nz] = 0.0; Z_n[nz] = 0.0; }
        zbar_n[nzmax_o] = mesh->zbar[nzmax_o];
        Z_n[nzmax_o - 1] = zbar_n[nzmax_o]
                         + mesh->hnode_new[(size_t)n * nl + (nzmax_o - 1)] * 0.5;
        for (int nz = nzmax_o - 1; nz >= nzmin_o + 1; --nz) {
            zbar_n[nz]    = zbar_n[nz + 1]
                          + mesh->hnode_new[(size_t)n * nl + nz];
            Z_n[nz - 1]   = zbar_n[nz]
                          + mesh->hnode_new[(size_t)n * nl + (nz - 1)] * 0.5;
        }
        zbar_n[nzmin_o] = zbar_n[nzmin_o + 1]
                        + mesh->hnode_new[(size_t)n * nl + nzmin_o];

        /* Inner bounds for TDMA (Fortran 92-93). */
        int nzmax = mesh->nlevels_nod2D_min[n] - 1;
        int nzmin = mesh->ulevels_nod2D_max[n] - 1;
        if (nzmax <= nzmin + 1) {
            /* No interior layers — Γ stays zero everywhere. */
            for (int nz = 0; nz < nl; ++nz) {
                gm->fer_gamma[(size_t)n * nl * 2 + nz * 2 + 0] = 0.0;
                gm->fer_gamma[(size_t)n * nl * 2 + nz * 2 + 1] = 0.0;
            }
            continue;
        }

        /* Top boundary (Dirichlet): a=c=0, b=1. */
        a[nzmin] = 0.0;
        c[nzmin] = 0.0;
        b[nzmin] = 1.0;

        /* zinv2 init for body (Fortran 104). */
        real_t zinv2 = 1.0 / (zbar_n[nzmin] - zbar_n[nzmin + 1]);

        /* Body: nz = nzmin+1 .. nzmax-1 (Fortran 107-114). */
        const real_t fc = gm->fer_C[n];
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t zinv1 = zinv2;
            zinv2 = 1.0 / (zbar_n[nz] - zbar_n[nz + 1]);
            real_t zinv  = 1.0 / (Z_n[nz - 1] - Z_n[nz]);
            a[nz] = fc * zinv1 * zinv;
            c[nz] = fc * zinv2 * zinv;
            real_t bv = aux->bvfreq[(size_t)n * nl + nz];
            if (bv < 1.0e-8) bv = 1.0e-8;
            b[nz] = -a[nz] - c[nz] - bv;
        }

        /* Bottom boundary (Dirichlet). */
        a[nzmax] = 0.0;
        c[nzmax] = 0.0;
        b[nzmax] = 1.0;

        /* RHS (Fortran 124-132). r = g/ρ₀; sigma_xy avg over two
         * adjacent mid-layers (nz-1, nz) at the level interface. */
        const real_t r = g / rho_ref;
        for (int nz = 0; nz < nl; ++nz) { tr_x[nz] = 0.0; tr_y[nz] = 0.0; }
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t sx_up = gm->sigma_xy[(size_t)n * nl * 2 + (nz - 1) * 2 + 0];
            real_t sx_dn = gm->sigma_xy[(size_t)n * nl * 2 + nz       * 2 + 0];
            real_t sy_up = gm->sigma_xy[(size_t)n * nl * 2 + (nz - 1) * 2 + 1];
            real_t sy_dn = gm->sigma_xy[(size_t)n * nl * 2 + nz       * 2 + 1];
            real_t k = gm->fer_K[(size_t)n * nl + nz];
            tr_x[nz] = r * 0.5 * (sx_up + sx_dn) * k;
            tr_y[nz] = r * 0.5 * (sy_up + sy_dn) * k;
        }

        /* Thomas sweep (Fortran 139-148). */
        cp[nzmin]   = c[nzmin]   / b[nzmin];
        tp_x[nzmin] = tr_x[nzmin] / b[nzmin];
        tp_y[nzmin] = tr_y[nzmin] / b[nzmin];
        for (int nz = nzmin + 1; nz <= nzmax; ++nz) {
            real_t m = b[nz] - cp[nz - 1] * a[nz];
            cp  [nz] = c[nz] / m;
            tp_x[nz] = (tr_x[nz] - tp_x[nz - 1] * a[nz]) / m;
            tp_y[nz] = (tr_y[nz] - tp_y[nz - 1] * a[nz]) / m;
        }

        /* Back-substitution into fer_gamma (Fortran 151-157). Zero
         * outside [nzmin, nzmax]. */
        for (int nz = 0; nz < nl; ++nz) {
            gm->fer_gamma[(size_t)n * nl * 2 + nz * 2 + 0] = 0.0;
            gm->fer_gamma[(size_t)n * nl * 2 + nz * 2 + 1] = 0.0;
        }
        gm->fer_gamma[(size_t)n * nl * 2 + nzmax * 2 + 0] = tp_x[nzmax];
        gm->fer_gamma[(size_t)n * nl * 2 + nzmax * 2 + 1] = tp_y[nzmax];
        for (int nz = nzmax - 1; nz >= nzmin; --nz) {
            real_t gx_below = gm->fer_gamma[(size_t)n * nl * 2 + (nz + 1) * 2 + 0];
            real_t gy_below = gm->fer_gamma[(size_t)n * nl * 2 + (nz + 1) * 2 + 1];
            gm->fer_gamma[(size_t)n * nl * 2 + nz * 2 + 0] = tp_x[nz] - cp[nz] * gx_below;
            gm->fer_gamma[(size_t)n * nl * 2 + nz * 2 + 1] = tp_y[nz] - cp[nz] * gy_below;
        }
    }

    /* Halo exchange fer_gamma (Fortran 162). 2 components per (n, nz). */
    if (partit && partit->npes > 1) {
        fesom_halo_exchange(gm->fer_gamma, FESOM_HALO_NOD2D, nl, 2, partit);
    }
}

/* ============================================================ */
/* G4 — fer_gamma2vel (oce_fer_gm.F90:168-210)                  */
/* ============================================================ */
/*
 * Element-loop reconstruction of bolus velocity from the streamfunction
 * differences across vertical interfaces:
 *   fer_uv(comp, nz, el) = (1/3) · Σ_v (Γ(comp, nz, v) − Γ(comp, nz+1, v))
 *                          / helem(nz, el)
 *
 * Halo: exchange_elem(fer_uv) — Fortran 208. Stride 2.
 */
void fesom_fer_gamma2vel(struct fesom_dyn *dyn,
                         const struct fesom_mesh *mesh,
                         const fesom_gm *gm,
                         struct fesom_partit *partit)
{
    const int  nl      = mesh->nl;
    const int  myDim_e = mesh->myDim_elem2D;
    const real_t onethird = 1.0 / 3.0;

    for (int el = 0; el < myDim_e; ++el) {
        int v0 = mesh->elem_nodes[3*el + 0];
        int v1 = mesh->elem_nodes[3*el + 1];
        int v2 = mesh->elem_nodes[3*el + 2];
        int nzmax = mesh->nlevels[el] - 1;       /* nl1 = layer count */
        int nzmin = mesh->ulevels[el] - 1;       /* 0-based */

        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t h = mesh->helem[(size_t)el * nl + nz];
            if (!(h > 0.0)) continue;            /* dry / sentinel layer */
            real_t zinv = onethird / h;

            real_t gx_top = gm->fer_gamma[(size_t)v0 * nl * 2 + nz       * 2 + 0]
                          + gm->fer_gamma[(size_t)v1 * nl * 2 + nz       * 2 + 0]
                          + gm->fer_gamma[(size_t)v2 * nl * 2 + nz       * 2 + 0];
            real_t gx_bot = gm->fer_gamma[(size_t)v0 * nl * 2 + (nz + 1) * 2 + 0]
                          + gm->fer_gamma[(size_t)v1 * nl * 2 + (nz + 1) * 2 + 0]
                          + gm->fer_gamma[(size_t)v2 * nl * 2 + (nz + 1) * 2 + 0];
            real_t gy_top = gm->fer_gamma[(size_t)v0 * nl * 2 + nz       * 2 + 1]
                          + gm->fer_gamma[(size_t)v1 * nl * 2 + nz       * 2 + 1]
                          + gm->fer_gamma[(size_t)v2 * nl * 2 + nz       * 2 + 1];
            real_t gy_bot = gm->fer_gamma[(size_t)v0 * nl * 2 + (nz + 1) * 2 + 1]
                          + gm->fer_gamma[(size_t)v1 * nl * 2 + (nz + 1) * 2 + 1]
                          + gm->fer_gamma[(size_t)v2 * nl * 2 + (nz + 1) * 2 + 1];

            dyn->fer_uv[(size_t)el * nl * 2 + nz * 2 + 0] = (gx_top - gx_bot) * zinv;
            dyn->fer_uv[(size_t)el * nl * 2 + nz * 2 + 1] = (gy_top - gy_bot) * zinv;
        }
    }

    if (partit && partit->npes > 1) {
        fesom_halo_exchange(dyn->fer_uv, FESOM_HALO_ELEM2D_FULL, nl, 2, partit);
    }
}
