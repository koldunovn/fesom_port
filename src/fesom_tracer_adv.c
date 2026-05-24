/*
 * Phase 1 simple-upwind tracer advection. Literal port — see header for
 * Fortran source line references.
 *
 * Fortran's upwind flux convention (oce_adv_tra_hor.F90:174-177):
 *
 *     vflux  = (-v*dx + u*dy) * helem               (dx,dy = el → edge_mid)
 *     flux  := -0.5 * [ T(n1) * (vflux + |vflux|)
 *                     + T(n2) * (vflux - |vflux|) ] - flux
 *
 * "- flux" in the writeback means the routine subtracts the previous content
 * (so calling with l_init_zero=true gives a fresh write; l_init_zero=false
 * gives an antidiffusive HO-LO difference). Here we are always called with
 * effective l_init_zero=true (UPW1 path) — we explicitly zero the flux array
 * up-front and the writeback is then `flux := -0.5*(...) - 0`.
 *
 * Sign convention sanity (FRESH_START.md §14.5 — "edge_vflux sign was reversed
 * in a previous port"). The vertical flux formula is:
 *
 *     flux(nz, n) := -0.5 * [ T(nz)   * (W(nz) + |W(nz)|)
 *                           + T(nz-1) * (W(nz) - |W(nz)|) ] * area(nz, n) - flux
 *
 * Note: vertical W is positive UPWARD. So at an interface where W>0 (upward
 * flow), the upwind-source is the cell BELOW the interface (T(nz)) — the
 * (W+|W|) factor doubles to 2W in that case, picking T(nz). Get this wrong
 * and tracer advection runs in reverse. Verified against the Fortran code.
 *
 * Surface vertical flux uses T(nz=top) directly (line 301):
 *     flux(top, n) := -W(top, n) * T(top, n) * area(top, n) - flux
 * Bottom vertical flux is set to 0 (line 306).
 *
 * STORAGE: all node-column arrays use stride `nl` and FESOM_NODE3D indexing —
 * matches mesh->hnode and tracers->del_ttf (allocated as N*nl). The unused
 * slot at level (nl-1) stays 0; this trades a few bytes for consistency.
 */
#include "fesom_tracer_adv.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void fesom_tracer_adv_init(fesom_tracer_adv_scratch *sc,
                           const struct fesom_mesh  *mesh)
{
    memset(sc, 0, sizeof(*sc));
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int E = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    int EG = mesh->myDim_edge2D + mesh->eDim_edge2D;
    size_t e_full = (size_t)EG * (size_t)mesh->nl;
    size_t n_full = (size_t)N  * (size_t)mesh->nl;

    sc->adv_flux_hor     = calloc(e_full,         sizeof(real_t));
    sc->adv_flux_ver     = calloc(n_full,         sizeof(real_t));
    sc->del_ttf_advhoriz = calloc(n_full,         sizeof(real_t));
    sc->del_ttf_advvert  = calloc(n_full,         sizeof(real_t));
    sc->fct_LO           = calloc(n_full,         sizeof(real_t));
    sc->fct_ttf_min      = calloc(n_full,         sizeof(real_t));
    sc->fct_ttf_max      = calloc(n_full,         sizeof(real_t));
    sc->fct_plus         = calloc(n_full,         sizeof(real_t));
    sc->fct_minus        = calloc(n_full,         sizeof(real_t));
    sc->fct_aux          = calloc((size_t)E * (size_t)mesh->nl * 2,
                                  sizeof(real_t));
    sc->tr_xy            = calloc((size_t)E * (size_t)mesh->nl * 2, sizeof(real_t));
    sc->edge_up_dn_grad  = calloc((size_t)mesh->myDim_edge2D * (size_t)mesh->nl * 4,
                                  sizeof(real_t));
    FESOM_CHECK(sc->adv_flux_hor && sc->adv_flux_ver
             && sc->del_ttf_advhoriz && sc->del_ttf_advvert
             && sc->fct_LO && sc->fct_ttf_min && sc->fct_ttf_max
             && sc->fct_plus && sc->fct_minus && sc->fct_aux
             && sc->tr_xy && sc->edge_up_dn_grad,
             "tracer adv scratch: out of memory");
}

void fesom_tracer_adv_free(fesom_tracer_adv_scratch *sc)
{
    free(sc->adv_flux_hor);
    free(sc->adv_flux_ver);
    free(sc->del_ttf_advhoriz);
    free(sc->del_ttf_advvert);
    free(sc->fct_LO);
    free(sc->fct_ttf_min);
    free(sc->fct_ttf_max);
    free(sc->fct_plus);
    free(sc->fct_minus);
    free(sc->fct_aux);
    free(sc->tr_xy);
    free(sc->edge_up_dn_grad);
    memset(sc, 0, sizeof(*sc));
}

/*--- compute_fct_LO --------------------------------------------------------
 * Mirror of oce_adv_tra_driver.F90:118-236 (FCT branch).
 * Pre: adv_flux_hor / adv_flux_ver hold UPWIND fluxes for the chosen tracer
 *      (i.e., adv_tra_hor_upw1 + adv_tra_ver_upw1 just ran with init_zero=true).
 *
 * Step 1: zero fct_LO, then accumulate horizontal flux divergence per node:
 *           fct_LO[nz, n1] += adv_flux_hor[nz, e]
 *           fct_LO[nz, n2] -= adv_flux_hor[nz, e]
 * Step 2: per-node finalisation:
 *           fct_LO[nz, n] = (T*hnode + (fct_LO + (vflux[nz] - vflux[nz+1]))
 *                                    * dt/areasvol) / hnode_new
 */
void fesom_tracer_compute_fct_LO(fesom_tracer_adv_scratch *sc,
                                 int                       tr_idx,
                                 const struct fesom_mesh  *mesh,
                                 const struct fesom_tracers *tracers)
{
    const int N    = mesh->myDim_nod2D;
    const int E    = mesh->myDim_edge2D;   /* Fortran adv loops are myDim only */
    const int nl   = mesh->nl;
    const real_t dt = (real_t)FESOM_PHASE1_DT;
    const real_t *T = tracers->data[tr_idx].values;

    /* Step 1a: zero fct_LO */
    memset(sc->fct_LO, 0, (size_t)N * (size_t)nl * sizeof(real_t));

    /* Step 1b: accumulate horizontal flux divergence per node.
       Range matches oce_adv_tra_driver.F90:170 — nu12 to nl12. */
    for (int e = 0; e < E; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;

        int nl12 = (nl1 > nl2) ? nl1 : nl2;
        int nu12 = nu1;
        if (nu2 > nu12) nu12 = nu2;

        for (int nz = nu12; nz < nl12; ++nz) {
            real_t f = sc->adv_flux_hor[FESOM_NODE3D(e, nz, nl)];
            sc->fct_LO[FESOM_NODE3D(n1, nz, nl)] += f;
            sc->fct_LO[FESOM_NODE3D(n2, nz, nl)] -= f;
        }
    }

    /* Step 2: per-node finalisation */
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t a = mesh->areasvol[k];
            real_t hnod_old = mesh->hnode    [k];
            real_t hnod_new = mesh->hnode_new[k];
            if (hnod_new <= 0.0 || a <= 0.0) {
                sc->fct_LO[k] = 0.0;
                continue;
            }
            real_t f_top = sc->adv_flux_ver[FESOM_NODE3D(n, nz,     nl)];
            real_t f_bot = sc->adv_flux_ver[FESOM_NODE3D(n, nz + 1, nl)];
            real_t numer = T[k] * hnod_old
                         + (sc->fct_LO[k] + (f_top - f_bot)) * dt / a;
            sc->fct_LO[k] = numer / hnod_new;
        }
    }
}

/*--- init_tracers_AB (oce_tracer_mod.F90:13-145, AB2 path) ----------------
 * All loops cover myDim_nod2D + eDim_nod2D (literal Fortran port).
 * Using myDim alone leaves halo entries of valuesAB / valuesold stale, which
 * makes the high-order horizontal advection read wrong tracer values on
 * halo endpoints of boundary edges → partition-dependent tracer drift. */
static void init_tracers_AB_one(int tr_idx,
                                const struct fesom_mesh *mesh,
                                struct fesom_tracers    *tracers,
                                fesom_tracer_adv_scratch *sc)
{
    const int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl = mesh->nl;
    /* AB2 stabilization offset (same as momentum, fesom_momentum.c): Fortran
       o_PARAM epsilon=0.1 (oce_modules.F90:92), used in the tracer AB2 at
       oce_tracer_mod.F90:53 as -(0.5+ε)*old+(1.5+ε)*new. The prior 1e-9 was the
       same mis-port corrected in the momentum RHS. */
    const real_t eps = 0.1;
    const real_t c_old = -(0.5 + eps);
    const real_t c_new =  (1.5 + eps);

    /* Zero del_ttf and the partial accumulators (lines 32-39). */
    size_t n_full = (size_t)N * (size_t)nl * sizeof(real_t);
    memset(tracers->del_ttf,        0, n_full);
    memset(sc->del_ttf_advhoriz,    0, n_full);
    memset(sc->del_ttf_advvert,     0, n_full);

    /* AB2 valuesAB = -(0.5+ε)*valuesold + (1.5+ε)*values  (lines 47-53). */
    real_t *vals    = tracers->data[tr_idx].values;
    real_t *valsAB  = tracers->data[tr_idx].valuesAB;
    real_t *valsold = tracers->data[tr_idx].valuesold;
    for (int n = 0; n < N; ++n) {
        for (int nz = 0; nz < nl; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            valsAB[k] = c_old * valsold[k] + c_new * vals[k];
        }
    }

    /* Save valuesold = values (lines 100-102, AB2). Runs after the valuesAB
       computation (Fortran sequence at lines 94-107). */
    memcpy(valsold, vals, n_full);
}

/*--- adv_tra_hor_upw1 (oce_adv_tra_hor.F90:64-257) ------------------------*/
static void adv_tra_hor_upw1(const struct fesom_mesh *mesh,
                             const struct fesom_dyn  *dyn,
                             const real_t            *ttf,
                             real_t                  *flux)
{
    const int E    = mesh->myDim_edge2D;   /* Fortran adv loops are myDim only */
    const int nl   = mesh->nl;

    /* Zero the flux array (l_init_zero=.true.) — lines 91-107 */
    memset(flux, 0, (size_t)E * (size_t)nl * sizeof(real_t));

    for (int e = 0; e < E; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;

        real_t dx1 = (el1 >= 0) ? mesh->edge_cross_dxdy[4*e + 0] : 0.0;
        real_t dy1 = (el1 >= 0) ? mesh->edge_cross_dxdy[4*e + 1] : 0.0;
        real_t dx2 = (el2 >= 0) ? mesh->edge_cross_dxdy[4*e + 2] : 0.0;
        real_t dy2 = (el2 >= 0) ? mesh->edge_cross_dxdy[4*e + 3] : 0.0;

        int nl12 = (el2 >= 0) ? ((nl1 < nl2) ? nl1 : nl2) : 0;
        int nu12 = nu1;
        if (nu2 > nu12) nu12 = nu2;

        /* (A) el1-only zone above shared layers (Fortran lines 167-178) */
        for (int nz = nu1; nz < nu12; ++nz) {
            real_t u  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
            real_t v  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
            real_t h  = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
            real_t vflux = (-v*dx1 + u*dy1) * h;
            real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            flux[FESOM_NODE3D(e, nz, nl)] = -0.5 * (T1*(vflux + fabs(vflux))
                                                  + T2*(vflux - fabs(vflux)));
        }

        /* (B) el2-only zone above shared layers (lines 185-199) */
        if (el2 >= 0) {
            for (int nz = nu2; nz < nu12; ++nz) {
                real_t u  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h  = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                real_t vflux = (v*dx2 - u*dy2) * h;
                real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
                real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
                flux[FESOM_NODE3D(e, nz, nl)] = -0.5 * (T1*(vflux + fabs(vflux))
                                                      + T2*(vflux - fabs(vflux)));
            }
        }

        /* (C) Both segments (lines 207-217) */
        for (int nz = nu12; nz < nl12; ++nz) {
            real_t u1 = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
            real_t v1 = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
            real_t h1 = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
            real_t u2 = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
            real_t v2 = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
            real_t h2 = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
            real_t vflux = (-v1*dx1 + u1*dy1) * h1
                         + ( v2*dx2 - u2*dy2) * h2;
            real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            flux[FESOM_NODE3D(e, nz, nl)] = -0.5 * (T1*(vflux + fabs(vflux))
                                                  + T2*(vflux - fabs(vflux)));
        }

        /* (D) el1-only remaining (lines 223-233) */
        for (int nz = nl12; nz < nl1; ++nz) {
            real_t u  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
            real_t v  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
            real_t h  = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
            real_t vflux = (-v*dx1 + u*dy1) * h;
            real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            flux[FESOM_NODE3D(e, nz, nl)] = -0.5 * (T1*(vflux + fabs(vflux))
                                                  + T2*(vflux - fabs(vflux)));
        }

        /* (E) el2-only remaining (lines 239-248) */
        if (el2 >= 0) {
            for (int nz = nl12; nz < nl2; ++nz) {
                real_t u  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h  = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                real_t vflux = (v*dx2 - u*dy2) * h;
                real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
                real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
                flux[FESOM_NODE3D(e, nz, nl)] = -0.5 * (T1*(vflux + fabs(vflux))
                                                      + T2*(vflux - fabs(vflux)));
            }
        }
    }
}

/*--- adv_tra_hor_central (interim Phase 2 high-order) --------------------
 * Central-difference horizontal flux. Two-sided face value Tmean = 0.5(T1+T2).
 * Equivalent to FESOM's MFCT (oce_adv_tra_hor.F90:546-834) with the
 * edge_up_dn_grad correction zeroed out — i.e., zero-gradient interpretation.
 *
 * Phase 2 interim choice: full MFCT pulls in tracer_gradient_elements +
 * fill_up_dn_grad + a per-edge upwind/downwind triangle index that we
 * haven't ported yet. Central is 2nd-order vs MFCT's 3rd-order with
 * upwind blend, but it is a valid HO scheme for FCT — the antidiffusive
 * flux (HO − LO) it produces is correctly limited by Zalesak.
 *
 * When init_zero is FALSE, this writes (HO − flux_existing); for FCT
 * use, callers run upwind with init_zero=true first, then central with
 * init_zero=false to obtain the antidiffusive horizontal flux.
 */
/* Retained on purpose (currently unused): the 2nd-order central HO flux is kept
 * available as an alternate FCT high-order scheme now that MFCT is the default.
 * Do NOT remove. __attribute__((unused)) silences the no-caller warning. */
__attribute__((unused))
static void adv_tra_hor_central(const struct fesom_mesh *mesh,
                                const struct fesom_dyn  *dyn,
                                const real_t            *ttf,
                                int                      init_zero,
                                real_t                  *flux)
{
    const int E    = mesh->myDim_edge2D;   /* Fortran adv loops are myDim only */
    const int nl   = mesh->nl;

    if (init_zero) {
        memset(flux, 0, (size_t)E * (size_t)nl * sizeof(real_t));
    }

    for (int e = 0; e < E; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;

        real_t dx1 = (el1 >= 0) ? mesh->edge_cross_dxdy[4*e + 0] : 0.0;
        real_t dy1 = (el1 >= 0) ? mesh->edge_cross_dxdy[4*e + 1] : 0.0;
        real_t dx2 = (el2 >= 0) ? mesh->edge_cross_dxdy[4*e + 2] : 0.0;
        real_t dy2 = (el2 >= 0) ? mesh->edge_cross_dxdy[4*e + 3] : 0.0;

        int nl12 = (el2 >= 0) ? ((nl1 < nl2) ? nl1 : nl2) : 0;
        int nu12 = (nu1 > nu2) ? nu1 : nu2;

        /* Same 5-zone structure as upwind (oce_adv_tra_hor.F90:167-249).
         * Only the per-layer flux formula differs:
         *   upwind   : flux = -0.5 [T1·(v+|v|) + T2·(v-|v|)]
         *   central  : flux = -v · 0.5(T1+T2)
         * The "− flux" writeback for init_zero=false yields the antidiffusive
         * flux when called after the upwind pass. */

        /* (A) el1-only zone above shared layers */
        for (int nz = nu1; nz < nu12; ++nz) {
            real_t u  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
            real_t v  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
            real_t h  = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
            real_t vflux = (-v*dx1 + u*dy1) * h;
            real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            real_t hi = -vflux * 0.5 * (T1 + T2);
            flux[FESOM_NODE3D(e, nz, nl)] = hi - flux[FESOM_NODE3D(e, nz, nl)];
        }
        /* (B) el2-only zone above shared */
        if (el2 >= 0) {
            for (int nz = nu2; nz < nu12; ++nz) {
                real_t u  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h  = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                real_t vflux = (v*dx2 - u*dy2) * h;
                real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
                real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
                real_t hi = -vflux * 0.5 * (T1 + T2);
                flux[FESOM_NODE3D(e, nz, nl)] = hi - flux[FESOM_NODE3D(e, nz, nl)];
            }
        }
        /* (C) Both segments — interior */
        for (int nz = nu12; nz < nl12; ++nz) {
            real_t u1 = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
            real_t v1 = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
            real_t h1 = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
            real_t u2 = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
            real_t v2 = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
            real_t h2 = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
            real_t vflux = (-v1*dx1 + u1*dy1) * h1
                         + ( v2*dx2 - u2*dy2) * h2;
            real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            real_t hi = -vflux * 0.5 * (T1 + T2);
            flux[FESOM_NODE3D(e, nz, nl)] = hi - flux[FESOM_NODE3D(e, nz, nl)];
        }
        /* (D) el1-only remaining */
        for (int nz = nl12; nz < nl1; ++nz) {
            real_t u  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
            real_t v  = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
            real_t h  = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
            real_t vflux = (-v*dx1 + u*dy1) * h;
            real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            real_t hi = -vflux * 0.5 * (T1 + T2);
            flux[FESOM_NODE3D(e, nz, nl)] = hi - flux[FESOM_NODE3D(e, nz, nl)];
        }
        /* (E) el2-only remaining */
        if (el2 >= 0) {
            for (int nz = nl12; nz < nl2; ++nz) {
                real_t u  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v  = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h  = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                real_t vflux = (v*dx2 - u*dy2) * h;
                real_t T1 = ttf[FESOM_NODE3D(n1, nz, nl)];
                real_t T2 = ttf[FESOM_NODE3D(n2, nz, nl)];
                real_t hi = -vflux * 0.5 * (T1 + T2);
                flux[FESOM_NODE3D(e, nz, nl)] = hi - flux[FESOM_NODE3D(e, nz, nl)];
            }
        }
    }
}

/*===========================================================================
 * MFCT 3rd-order horizontal tracer advection (replaces adv_tra_hor_central).
 *   tracer_gradient_elements (oce_tracer_mod.F90:175-185)
 *   fill_up_dn_grad          (oce_muscl_adv.F90:356-525)
 *   adv_tra_hor_mfct         (oce_adv_tra_hor.F90:546-834)
 *===========================================================================*/

/* tr_xy = gradient_sca·ttf per element per layer; on myDim then wide-exchanged. */
static void tracer_gradient_elements(const struct fesom_mesh *mesh,
                                     const real_t *ttf, real_t *tr_xy,
                                     struct fesom_partit *partit)
{
    const int nl = mesh->nl;
    const int E  = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    memset(tr_xy, 0, (size_t)E * (size_t)nl * 2 * sizeof(real_t));
    for (int el = 0; el < mesh->myDim_elem2D; ++el) {
        int n0 = mesh->elem_nodes[3*el + 0];
        int n1 = mesh->elem_nodes[3*el + 1];
        int n2 = mesh->elem_nodes[3*el + 2];
        const real_t *gs = &mesh->gradient_sca[6*el];
        int nu = mesh->ulevels[el] - 1, nle = mesh->nlevels[el] - 1;
        for (int nz = nu; nz < nle; ++nz) {
            real_t t0 = ttf[FESOM_NODE3D(n0, nz, nl)];
            real_t t1 = ttf[FESOM_NODE3D(n1, nz, nl)];
            real_t t2 = ttf[FESOM_NODE3D(n2, nz, nl)];
            tr_xy[FESOM_ELEMVEC(el, nz, nl) + 0] = gs[0]*t0 + gs[1]*t1 + gs[2]*t2;
            tr_xy[FESOM_ELEMVEC(el, nz, nl) + 1] = gs[3]*t0 + gs[4]*t1 + gs[5]*t2;
        }
    }
    if (partit && partit->npes > 1)
        fesom_halo_exchange(tr_xy, FESOM_HALO_ELEM2D_FULL, nl, 2, partit);
}

/* area-weighted mean of tr_xy over node's surrounding elements that own layer nz,
 * written into g[nz*4+cx], g[nz*4+cy], for nz in [nz0, nz1). (F:448-457 etc.) */
static void node_avg_grad(const struct fesom_mesh *mesh, const real_t *tr_xy,
                          int node, int nz0, int nz1, real_t *g, int cx, int cy)
{
    const int nl = mesh->nl;
    int b0 = mesh->nod_in_elem2D_offsets[node];
    int b1 = mesh->nod_in_elem2D_offsets[node + 1];
    for (int nz = nz0; nz < nz1; ++nz) {
        real_t tvol = 0.0, tx = 0.0, ty = 0.0;
        for (int k = b0; k < b1; ++k) {
            int el = mesh->nod_in_elem2D[k];
            if (nz >= mesh->nlevels[el] - 1 || nz < mesh->ulevels[el] - 1) continue;
            real_t a = mesh->elem_area[el];
            tvol += a;
            tx   += tr_xy[FESOM_ELEMVEC(el, nz, nl) + 0] * a;
            ty   += tr_xy[FESOM_ELEMVEC(el, nz, nl) + 1] * a;
        }
        if (tvol > 0.0) { g[nz*4 + cx] = tx / tvol; g[nz*4 + cy] = ty / tvol; }
    }
}

/* edge_up_dn_grad[4 per (edge,level)]: shared levels use the up/down-triangle
 * gradient; non-shared (cavity/bathy) + boundary edges use node-averaged. */
static void fill_up_dn_grad(const struct fesom_mesh *mesh, const real_t *tr_xy,
                            real_t *eud)
{
    const int nl  = mesh->nl;
    const int Eme = mesh->myDim_edge2D;
    memset(eud, 0, (size_t)Eme * (size_t)nl * 4 * sizeof(real_t));
    for (int ed = 0; ed < Eme; ++ed) {
        int n1 = mesh->edges[2*ed + 0];
        int n2 = mesh->edges[2*ed + 1];
        int up = mesh->edge_up_dn_tri[2*ed + 0];
        int dn = mesh->edge_up_dn_tri[2*ed + 1];
        real_t *g = &eud[(size_t)ed * nl * 4];
        if (up >= 0 && dn >= 0) {
            int a1 = mesh->ulevels_nod2D_max[n1], a2 = mesh->ulevels_nod2D_max[n2];
            int nzmin = ((a1 > a2) ? a1 : a2) - 1;                   /* 0-based */
            int b1 = mesh->nlevels_nod2D_min[n1], b2 = mesh->nlevels_nod2D_min[n2];
            int nzmax = ((b1 < b2) ? b1 : b2) - 1;
            node_avg_grad(mesh, tr_xy, n1, mesh->ulevels_nod2D[n1]-1, nzmin, g, 0, 2);
            node_avg_grad(mesh, tr_xy, n2, mesh->ulevels_nod2D[n2]-1, nzmin, g, 1, 3);
            for (int nz = nzmin; nz < nzmax; ++nz) {
                g[nz*4 + 0] = tr_xy[FESOM_ELEMVEC(up, nz, nl) + 0];
                g[nz*4 + 1] = tr_xy[FESOM_ELEMVEC(dn, nz, nl) + 0];
                g[nz*4 + 2] = tr_xy[FESOM_ELEMVEC(up, nz, nl) + 1];
                g[nz*4 + 3] = tr_xy[FESOM_ELEMVEC(dn, nz, nl) + 1];
            }
            node_avg_grad(mesh, tr_xy, n1, nzmax, mesh->nlevels_nod2D[n1]-1, g, 0, 2);
            node_avg_grad(mesh, tr_xy, n2, nzmax, mesh->nlevels_nod2D[n2]-1, g, 1, 3);
        } else {
            node_avg_grad(mesh, tr_xy, n1, mesh->ulevels_nod2D[n1]-1,
                          mesh->nlevels_nod2D[n1]-1, g, 0, 2);
            node_avg_grad(mesh, tr_xy, n2, mesh->ulevels_nod2D[n2]-1,
                          mesh->nlevels_nod2D[n2]-1, g, 1, 3);
        }
    }
}

/* 3rd-order reconstructed antidiffusive flux at one (edge, level). */
static inline void mfct_edge_flux(const real_t *ttf, real_t *flux, const real_t *eud,
                                  int n1, int n2, int e, int nz, int nl, real_t vflux,
                                  real_t a, real_t edx, real_t edy, real_t RE, real_t num_ord)
{
    real_t T1   = ttf[FESOM_NODE3D(n1, nz, nl)];
    real_t T2   = ttf[FESOM_NODE3D(n2, nz, nl)];
    real_t diff = T2 - T1;
    const real_t *g = &eud[((size_t)e * nl + nz) * 4];
    real_t Tmean1 = T1 + (2.0*diff + edx*a*g[0] + edy*RE*g[2]) / 6.0;
    real_t Tmean2 = T2 - (2.0*diff + edx*a*g[1] + edy*RE*g[3]) / 6.0;
    real_t av  = fabs(vflux);
    real_t cHO = (vflux + av)*Tmean1 + (vflux - av)*Tmean2;
    real_t hi  = -0.5*(1.0 - num_ord)*cHO - vflux*num_ord*0.5*(Tmean1 + Tmean2);
    flux[FESOM_NODE3D(e, nz, nl)] = hi - flux[FESOM_NODE3D(e, nz, nl)];
}

static void adv_tra_hor_mfct(const struct fesom_mesh *mesh, const struct fesom_dyn *dyn,
                             const real_t *ttf, const real_t *eud, real_t num_ord,
                             int init_zero, real_t *flux)
{
    const int E  = mesh->myDim_edge2D;
    const int nl = mesh->nl;
    const real_t RE = (real_t)FESOM_R_EARTH;
    if (init_zero) memset(flux, 0, (size_t)E * (size_t)nl * sizeof(real_t));

    for (int e = 0; e < E; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];
        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;
        real_t dx1 = (el1 >= 0) ? mesh->edge_cross_dxdy[4*e + 0] : 0.0;
        real_t dy1 = (el1 >= 0) ? mesh->edge_cross_dxdy[4*e + 1] : 0.0;
        real_t dx2 = (el2 >= 0) ? mesh->edge_cross_dxdy[4*e + 2] : 0.0;
        real_t dy2 = (el2 >= 0) ? mesh->edge_cross_dxdy[4*e + 3] : 0.0;
        int nl12 = (el2 >= 0) ? ((nl1 < nl2) ? nl1 : nl2) : 0;
        int nu12 = (nu1 > nu2) ? nu1 : nu2;
        real_t edx = mesh->edge_dxdy[2*e + 0];
        real_t edy = mesh->edge_dxdy[2*e + 1];
        real_t a = (el1 >= 0) ? RE * mesh->elem_cos[el1] : 0.0;
        if (el2 >= 0) a = 0.5 * (a + RE * mesh->elem_cos[el2]);

        for (int nz = nu1; nz < nu12; ++nz) {                /* (A) el1-only above */
            real_t u = dyn->uv[FESOM_ELEMVEC(el1,nz,nl)+0], v = dyn->uv[FESOM_ELEMVEC(el1,nz,nl)+1];
            real_t vflux = (-v*dx1 + u*dy1) * mesh->helem[FESOM_ELEM3D(el1,nz,nl)];
            mfct_edge_flux(ttf, flux, eud, n1, n2, e, nz, nl, vflux, a, edx, edy, RE, num_ord);
        }
        if (el2 >= 0) for (int nz = nu2; nz < nu12; ++nz) {  /* (B) el2-only above */
            real_t u = dyn->uv[FESOM_ELEMVEC(el2,nz,nl)+0], v = dyn->uv[FESOM_ELEMVEC(el2,nz,nl)+1];
            real_t vflux = (v*dx2 - u*dy2) * mesh->helem[FESOM_ELEM3D(el2,nz,nl)];
            mfct_edge_flux(ttf, flux, eud, n1, n2, e, nz, nl, vflux, a, edx, edy, RE, num_ord);
        }
        for (int nz = nu12; nz < nl12; ++nz) {               /* (C) both */
            real_t u1 = dyn->uv[FESOM_ELEMVEC(el1,nz,nl)+0], v1 = dyn->uv[FESOM_ELEMVEC(el1,nz,nl)+1];
            real_t u2 = dyn->uv[FESOM_ELEMVEC(el2,nz,nl)+0], v2 = dyn->uv[FESOM_ELEMVEC(el2,nz,nl)+1];
            real_t vflux = (-v1*dx1 + u1*dy1) * mesh->helem[FESOM_ELEM3D(el1,nz,nl)]
                         + ( v2*dx2 - u2*dy2) * mesh->helem[FESOM_ELEM3D(el2,nz,nl)];
            mfct_edge_flux(ttf, flux, eud, n1, n2, e, nz, nl, vflux, a, edx, edy, RE, num_ord);
        }
        for (int nz = nl12; nz < nl1; ++nz) {                /* (D) el1-only below */
            real_t u = dyn->uv[FESOM_ELEMVEC(el1,nz,nl)+0], v = dyn->uv[FESOM_ELEMVEC(el1,nz,nl)+1];
            real_t vflux = (-v*dx1 + u*dy1) * mesh->helem[FESOM_ELEM3D(el1,nz,nl)];
            mfct_edge_flux(ttf, flux, eud, n1, n2, e, nz, nl, vflux, a, edx, edy, RE, num_ord);
        }
        if (el2 >= 0) for (int nz = nl12; nz < nl2; ++nz) {  /* (E) el2-only below */
            real_t u = dyn->uv[FESOM_ELEMVEC(el2,nz,nl)+0], v = dyn->uv[FESOM_ELEMVEC(el2,nz,nl)+1];
            real_t vflux = (v*dx2 - u*dy2) * mesh->helem[FESOM_ELEM3D(el2,nz,nl)];
            mfct_edge_flux(ttf, flux, eud, n1, n2, e, nz, nl, vflux, a, edx, edy, RE, num_ord);
        }
    }
}

/*--- adv_tra_ver_qr4c (oce_adv_tra_ver.F90:332-434) -----------------------
 * Vertical 4th-order quadratic reconstruction with upwind blend.
 *   num_ord ∈ [0, 1]: fraction of pure 4th-order (1.0 = full 4th-order,
 *                     0.0 = full upwind-blend). Pi config sets it to 1.
 *
 * For Phase 1 (no partial cells, no cavity), Z_3d_n[nz, n] = mesh->Z[nz] and
 * zbar_3d_n[nz, n] = mesh->zbar[nz] for all n.
 *
 * Sign convention (cf. FRESH_START.md §14.5): vertical W is positive UPWARD.
 * Surface flux = -ttf·W·area  (W>0 advects ttf out of the surface layer).
 *
 * When init_zero is FALSE, this writes (HO_flux - flux_existing); for FCT
 * use, callers run upwind first with init_zero=true, then QR4C with
 * init_zero=false to obtain the antidiffusive vertical flux.
 */
static void adv_tra_ver_qr4c(const struct fesom_mesh *mesh,
                             const struct fesom_dyn  *dyn,
                             const real_t            *ttf,
                             real_t                   num_ord,
                             int                      init_zero,
                             real_t                  *flux)
{
    const int N  = mesh->myDim_nod2D;
    const int nl = mesh->nl;
    /* Tracer advection uses the EXPLICIT part of w (w_e). When wsplit is
       not triggered (CFL ≤ maxcfl) compute_Wvel_split sets w_e = w. */
    const real_t *W = dyn->w_e;

    if (init_zero) {
        memset(flux, 0, (size_t)N * (size_t)nl * sizeof(real_t));
    }

    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        if (nzmax <= nzmin) continue;

        /* Surface (Fortran nz=nzmin): -ttf·W·area - flux_existing */
        {
            int nz = nzmin;
            size_t k = FESOM_NODE3D(n, nz, nl);
            flux[k] = -ttf[k] * W[k] * mesh->area[k] - flux[k];
        }
        /* 2nd layer (nzmin+1): centred differences. Guard for shallow cols. */
        if (nzmax - nzmin >= 2) {
            int nz = nzmin + 1;
            real_t Tup = ttf[FESOM_NODE3D(n, nz - 1, nl)];
            real_t Tdn = ttf[FESOM_NODE3D(n, nz,     nl)];
            size_t k = FESOM_NODE3D(n, nz, nl);
            flux[k] = -0.5 * (Tup + Tdn) * W[k] * mesh->area[k] - flux[k];
        }
        /* Bottom - 1 (nzmax-1): centred differences. Same as second-layer. */
        if (nzmax - nzmin >= 2) {
            int nz = nzmax - 1;
            real_t Tup = ttf[FESOM_NODE3D(n, nz - 1, nl)];
            real_t Tdn = ttf[FESOM_NODE3D(n, nz,     nl)];
            size_t k = FESOM_NODE3D(n, nz, nl);
            flux[k] = -0.5 * (Tup + Tdn) * W[k] * mesh->area[k] - flux[k];
        }
        /* Bottom (nzmax): zero flux, less existing. */
        {
            int nz = nzmax;
            size_t k = FESOM_NODE3D(n, nz, nl);
            flux[k] = 0.0 - flux[k];
        }
        /* Interior (nzmin+2 .. nzmax-2 inclusive): 4th-order quadratic. */
        for (int nz = nzmin + 2; nz <= nzmax - 2; ++nz) {
            real_t Z_um1 = mesh->Z[nz - 1];
            real_t Z_u   = mesh->Z[nz];
            real_t Z_dn  = mesh->Z[nz + 1];
            real_t Z_um2 = mesh->Z[nz - 2];
            real_t Tum1 = ttf[FESOM_NODE3D(n, nz - 1, nl)];
            real_t Tu   = ttf[FESOM_NODE3D(n, nz,     nl)];
            real_t Tdn  = ttf[FESOM_NODE3D(n, nz + 1, nl)];
            real_t Tum2 = ttf[FESOM_NODE3D(n, nz - 2, nl)];

            real_t qc = (Tum1 - Tu)  / (Z_um1 - Z_u);
            real_t qu = (Tu   - Tdn) / (Z_u   - Z_dn);
            real_t qd = (Tum2 - Tum1)/ (Z_um2 - Z_um1);

            real_t zb = mesh->zbar[nz];
            real_t Tmean1 = Tu  + (2.0*qc + qu) * (zb - Z_u  ) / 3.0;
            real_t Tmean2 = Tum1 + (2.0*qc + qd) * (zb - Z_um1) / 3.0;
            real_t w_iface = W[FESOM_NODE3D(n, nz, nl)];
            real_t Tmean = (w_iface + fabs(w_iface)) * Tmean1
                         + (w_iface - fabs(w_iface)) * Tmean2;
            real_t a = mesh->area[FESOM_NODE3D(n, nz, nl)];
            real_t hi = (-0.5 * (1.0 - num_ord) * Tmean
                       - num_ord * 0.5 * (Tmean1 + Tmean2) * w_iface) * a;
            flux[FESOM_NODE3D(n, nz, nl)] = hi - flux[FESOM_NODE3D(n, nz, nl)];
        }
    }
}

/*--- adv_tra_ver_upw1 (oce_adv_tra_ver.F90:244-328) -----------------------*/
static void adv_tra_ver_upw1(const struct fesom_mesh *mesh,
                             const struct fesom_dyn  *dyn,
                             const real_t            *ttf,
                             real_t                  *flux)
{
    const int N  = mesh->myDim_nod2D;
    const int nl = mesh->nl;
    /* Use explicit-W (w_e). compute_Wvel_split sets w_e = w when CFL ≤ maxcfl. */
    const real_t *W = dyn->w_e;

    memset(flux, 0, (size_t)N * (size_t)nl * sizeof(real_t));

    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        if (nzmax <= nzmin) continue;

        /* Surface flux (line 301) */
        flux[FESOM_NODE3D(n, nzmin, nl)] =
            -W[FESOM_NODE3D(n, nzmin, nl)]
            * ttf[FESOM_NODE3D(n, nzmin, nl)]
            * mesh->area[FESOM_NODE3D(n, nzmin, nl)];

        /* Bottom flux = 0 (line 306) — already zeroed by memset. */

        /* Interior interfaces (lines 315-319) */
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t w_iface = W[FESOM_NODE3D(n, nz, nl)];
            real_t T_below = ttf[FESOM_NODE3D(n, nz,     nl)];   /* layer index nz */
            real_t T_above = ttf[FESOM_NODE3D(n, nz - 1, nl)];   /* layer index nz-1 */
            real_t a       = mesh->area[FESOM_NODE3D(n, nz, nl)];
            flux[FESOM_NODE3D(n, nz, nl)] =
                -0.5 * (T_below * (w_iface + fabs(w_iface))
                      + T_above * (w_iface - fabs(w_iface))) * a;
        }
    }
}

/*--- oce_tra_adv_flux2dtracer (oce_adv_tra_driver.F90:423-575, non-FCT) ---*/
static void flux2dtracer_upwind(const struct fesom_mesh *mesh,
                                const real_t            *flux_h,
                                const real_t            *flux_v,
                                real_t                  *dttf_h,
                                real_t                  *dttf_v)
{
    const int N  = mesh->myDim_nod2D;
    const int E  = mesh->myDim_edge2D;   /* Fortran adv loops are myDim only */
    const int nl = mesh->nl;
    const real_t dt = (real_t)FESOM_PHASE1_DT;

    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            real_t a = mesh->areasvol[FESOM_NODE3D(n, nz, nl)];
            if (a <= 0.0) continue;
            real_t f_top = flux_v[FESOM_NODE3D(n, nz,     nl)];
            real_t f_bot = flux_v[FESOM_NODE3D(n, nz + 1, nl)];
            dttf_v[FESOM_NODE3D(n, nz, nl)] += (f_top - f_bot) * dt / a;
        }
    }

    for (int e = 0; e < E; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;

        int nl12 = (nl1 > nl2) ? nl1 : nl2;
        int nu12 = nu1;
        if (el2 >= 0 && nu2 > 0 && nu2 < nu12) nu12 = nu2;

        for (int nz = nu12; nz < nl12; ++nz) {
            real_t f  = flux_h[FESOM_NODE3D(e, nz, nl)];
            real_t a1 = mesh->areasvol[FESOM_NODE3D(n1, nz, nl)];
            real_t a2 = mesh->areasvol[FESOM_NODE3D(n2, nz, nl)];
            if (a1 > 0.0) dttf_h[FESOM_NODE3D(n1, nz, nl)] += f * dt / a1;
            if (a2 > 0.0) dttf_h[FESOM_NODE3D(n2, nz, nl)] -= f * dt / a2;
        }
    }
}

/*--- ALE reconstruction (diff_tracers_ale lines 481-484) ------------------
 *   del_ttf += T * (hnode - hnode_new)
 *   T       += del_ttf / hnode_new
 */
static void ale_reconstruct(const struct fesom_mesh *mesh,
                            real_t                  *T,
                            real_t                  *del_ttf)
{
    const int N  = mesh->myDim_nod2D;
    const int nl = mesh->nl;

    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t hnode_old = mesh->hnode    [k];
            real_t hnode_new = mesh->hnode_new[k];
            del_ttf[k] += T[k] * (hnode_old - hnode_new);
            if (hnode_new > 0.0) {
                T[k] += del_ttf[k] / hnode_new;
            }
        }
    }
}

/*--- Zalesak limiter (oce_adv_tra_fct.F90:72-500) -------------------------
 * Inputs:
 *   ttf:    tracer at step n (=values BEFORE advection)
 *   lo:     fct_LO (low-order upwind solution, advanced one step)
 *   adf_h:  HO−LO antidiffusive horizontal flux (in/out — in: anti-diff;
 *                                                       out: limited)
 *   adf_v:  HO−LO antidiffusive vertical flux   (in/out)
 *
 * Workspace (in scratch):
 *   fct_ttf_min/max  (writes per-node admissible-increment bounds)
 *   fct_plus/minus   (limiter factors)
 *   fct_aux          (per-element max/min, reused as Fortran AUX(:,:,:))
 *
 * Algorithm:
 *   a1. fct_ttf_max = max(LO, ttf), fct_ttf_min = min(LO, ttf)   per node
 *   a2. AUX[1] = max over 3 vertices of fct_ttf_max,
 *       AUX[2] = min                                              per element
 *   a3. tvert_max = max over surrounding cells of AUX[1],
 *       tvert_min = min                                            per node
 *   a4. fct_ttf_max[surface] = tvert_max[surface] − LO,
 *       fct_ttf_max[middle ] = max over (nz-1, nz, nz+1) − LO,
 *       fct_ttf_max[bottom ] = tvert_max[bottom] − LO              (and min)
 *   b1. fct_plus  = Σ (max(0, +adf_v[nz]) + max(0, −adf_v[nz+1]))  per node
 *                 + Σ over edges of max(0, ±adf_h)
 *       fct_minus = same with min instead of max
 *   b2. flux  = fct_plus  · dt / areasvol / hnode_new + flux_eps
 *       fct_plus  = min(1, fct_ttf_max / flux)
 *       (similar for fct_minus, with −flux_eps)
 *   b3. Limiting (vertical):
 *         surface: ae = (adf_v ≥ 0 ? fct_plus[nz] : fct_minus[nz])
 *         interior: ae = adf_v ≥ 0 ? min(fct_minus[nz-1], fct_plus[nz])
 *                                   : min(fct_plus[nz-1],  fct_minus[nz])
 *         adf_v *= ae
 *       Limiting (horizontal): adf_h ≥ 0 ? min(fct_plus[n1],fct_minus[n2])
 *                                         : min(fct_minus[n1],fct_plus[n2])
 *         adf_h *= ae
 */
static void oce_tra_adv_fct(fesom_tracer_adv_scratch *sc,
                            const struct fesom_mesh  *mesh,
                            const real_t             *ttf,
                            real_t                   *adf_h,
                            real_t                   *adf_v)
{
    const int N    = mesh->myDim_nod2D;
    const int E    = mesh->myDim_elem2D;
    const int Eedg = mesh->myDim_edge2D;   /* Fortran FCT loops are myDim only */
    const int nl   = mesh->nl;
    const real_t dt        = (real_t)FESOM_PHASE1_DT;
    const real_t flux_eps  = 1e-16;
    const real_t bignumber = 1e3;

    real_t *fct_ttf_min = sc->fct_ttf_min;
    real_t *fct_ttf_max = sc->fct_ttf_max;
    real_t *fct_plus    = sc->fct_plus;
    real_t *fct_minus   = sc->fct_minus;
    real_t *AUX         = sc->fct_aux;     /* [E * nl * 2] : (max=0, min=1) */
    const real_t *LO    = sc->fct_LO;

    /* a1. Per-node max/min between LO and ttf.
     * Fortran oce_adv_tra_fct.F90:124 loops myDim+eDim so the halo entries
     * of fct_ttf_{min,max} match the owner — step a2 reads these values at
     * halo vertices of owned elements. Using myDim alone leaves halo entries
     * stale (partition-dependent tracer drift; this is bug 2 behind the
     * init_tracers_AB fix). */
    const int N_full_a1 = mesh->myDim_nod2D + mesh->eDim_nod2D;
    for (int n = 0; n < N_full_a1; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t a = LO[k], b = ttf[k];
            fct_ttf_max[k] = (a > b) ? a : b;
            fct_ttf_min[k] = (a < b) ? a : b;
        }
    }

    /* a2. Per-element max/min over 3 vertices. Pad below cell's nlevels with
       ±bignumber so the per-node tvert step never picks a stale layer. */
    for (int e = 0; e < E; ++e) {
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];
        int nu1 = mesh->ulevels[e] - 1;
        int nle = mesh->nlevels[e] - 1;
        for (int nz = nu1; nz < nle; ++nz) {
            real_t a0 = fct_ttf_max[FESOM_NODE3D(n0, nz, nl)];
            real_t a1 = fct_ttf_max[FESOM_NODE3D(n1, nz, nl)];
            real_t a2 = fct_ttf_max[FESOM_NODE3D(n2, nz, nl)];
            real_t mx = a0; if (a1 > mx) mx = a1; if (a2 > mx) mx = a2;
            real_t b0 = fct_ttf_min[FESOM_NODE3D(n0, nz, nl)];
            real_t b1 = fct_ttf_min[FESOM_NODE3D(n1, nz, nl)];
            real_t b2 = fct_ttf_min[FESOM_NODE3D(n2, nz, nl)];
            real_t mn = b0; if (b1 < mn) mn = b1; if (b2 < mn) mn = b2;
            AUX[FESOM_ELEM3D(e, nz, nl)*2 + 0] = mx;
            AUX[FESOM_ELEM3D(e, nz, nl)*2 + 1] = mn;
        }
        for (int nz = nle; nz < nl - 1; ++nz) {
            AUX[FESOM_ELEM3D(e, nz, nl)*2 + 0] = -bignumber;
            AUX[FESOM_ELEM3D(e, nz, nl)*2 + 1] = +bignumber;
        }
    }

    /* a3. Per-node max/min over surrounding cells (cluster). Use a per-node
       small buffer for tvert_max/min — allocated once, reused per node. */
    real_t *tvert_max = malloc((size_t)N * (size_t)nl * sizeof(real_t));
    real_t *tvert_min = malloc((size_t)N * (size_t)nl * sizeof(real_t));
    FESOM_CHECK(tvert_max && tvert_min, "FCT: oom (tvert)");
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        int o0 = mesh->nod_in_elem2D_offsets[n];
        int o1 = mesh->nod_in_elem2D_offsets[n + 1];
        for (int nz = nu1; nz < nl1; ++nz) {
            int e0 = mesh->nod_in_elem2D[o0];
            real_t mx = AUX[FESOM_ELEM3D(e0, nz, nl)*2 + 0];
            real_t mn = AUX[FESOM_ELEM3D(e0, nz, nl)*2 + 1];
            for (int k = o0 + 1; k < o1; ++k) {
                int e = mesh->nod_in_elem2D[k];
                real_t a = AUX[FESOM_ELEM3D(e, nz, nl)*2 + 0];
                real_t b = AUX[FESOM_ELEM3D(e, nz, nl)*2 + 1];
                if (a > mx) mx = a;
                if (b < mn) mn = b;
            }
            tvert_max[FESOM_NODE3D(n, nz, nl)] = mx;
            tvert_min[FESOM_NODE3D(n, nz, nl)] = mn;
        }
    }

    /* a4. Per-node admissible increment relative to LO.
          Includes vertical 3-layer cluster max/min (nz-1, nz, nz+1). */
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        if (nu1 >= nl1) continue;

        /* Surface */
        {
            int nz = nu1;
            size_t k = FESOM_NODE3D(n, nz, nl);
            fct_ttf_max[k] = tvert_max[k] - LO[k];
            fct_ttf_min[k] = tvert_min[k] - LO[k];
        }
        /* Interior */
        for (int nz = nu1 + 1; nz < nl1 - 1; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t mxA = tvert_max[FESOM_NODE3D(n, nz - 1, nl)];
            real_t mxB = tvert_max[FESOM_NODE3D(n, nz,     nl)];
            real_t mxC = tvert_max[FESOM_NODE3D(n, nz + 1, nl)];
            real_t mnA = tvert_min[FESOM_NODE3D(n, nz - 1, nl)];
            real_t mnB = tvert_min[FESOM_NODE3D(n, nz,     nl)];
            real_t mnC = tvert_min[FESOM_NODE3D(n, nz + 1, nl)];
            real_t mx = mxA; if (mxB > mx) mx = mxB; if (mxC > mx) mx = mxC;
            real_t mn = mnA; if (mnB < mn) mn = mnB; if (mnC < mn) mn = mnC;
            fct_ttf_max[k] = mx - LO[k];
            fct_ttf_min[k] = mn - LO[k];
        }
        /* Bottom layer (nl1 - 1) */
        if (nl1 - 1 > nu1) {
            int nz = nl1 - 1;
            size_t k = FESOM_NODE3D(n, nz, nl);
            fct_ttf_max[k] = tvert_max[k] - LO[k];
            fct_ttf_min[k] = tvert_min[k] - LO[k];
        }
    }
    free(tvert_max);
    free(tvert_min);

    /* b1. Sum positive/negative antidiffusive contributions per node.
       Vertical: from adf_v[nz] (in) and adf_v[nz+1] (out — flipped sign). */
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            fct_plus [k] = 0.0;
            fct_minus[k] = 0.0;
        }
    }
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            real_t fv_top  = adf_v[FESOM_NODE3D(n, nz,     nl)];
            real_t fv_bot  = adf_v[FESOM_NODE3D(n, nz + 1, nl)];
            real_t pos = (fv_top > 0.0 ? fv_top : 0.0)
                       + (-fv_bot > 0.0 ? -fv_bot : 0.0);
            real_t neg = (fv_top < 0.0 ? fv_top : 0.0)
                       + (-fv_bot < 0.0 ? -fv_bot : 0.0);
            size_t k = FESOM_NODE3D(n, nz, nl);
            fct_plus [k] += pos;
            fct_minus[k] += neg;
        }
    }
    /* Horizontal — per edge */
    for (int e = 0; e < Eedg; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;

        int nl12 = (nl1 > nl2) ? nl1 : nl2;
        int nu12 = nu1;
        if (el2 >= 0 && nu2 < nu12) nu12 = nu2;

        for (int nz = nu12; nz < nl12; ++nz) {
            real_t f = adf_h[FESOM_NODE3D(e, nz, nl)];
            size_t k1 = FESOM_NODE3D(n1, nz, nl);
            size_t k2 = FESOM_NODE3D(n2, nz, nl);
            fct_plus [k1] += (f > 0.0 ? f : 0.0);
            fct_minus[k1] += (f < 0.0 ? f : 0.0);
            fct_plus [k2] += (-f > 0.0 ? -f : 0.0);
            fct_minus[k2] += (-f < 0.0 ? -f : 0.0);
        }
    }

    /* b2. Per-node limiter factors. */
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t a = mesh->areasvol[k];
            real_t hnod_new = mesh->hnode_new[k];
            if (a <= 0.0 || hnod_new <= 0.0) {
                fct_plus[k]  = 1.0;
                fct_minus[k] = 1.0;
                continue;
            }
            real_t flux_pos = fct_plus[k]  * dt / a / hnod_new + flux_eps;
            real_t ratio_pos = fct_ttf_max[k] / flux_pos;
            fct_plus[k] = (ratio_pos < 1.0 ? ratio_pos : 1.0);

            real_t flux_neg = fct_minus[k] * dt / a / hnod_new - flux_eps;
            real_t ratio_neg = fct_ttf_min[k] / flux_neg;
            fct_minus[k] = (ratio_neg < 1.0 ? ratio_neg : 1.0);
        }
    }
    /* MPI: halo exchange of fct_plus / fct_minus so the limiter ratios are
     * consistent across rank boundaries. Mirrors Fortran
     *   call exchange_nod(fct_plus, fct_minus, partit, luse_g2g=.true.)
     * (oce_adv_tra_fct.F90:401). Without this, the per-edge horizontal
     * limiter at b3 picks an inconsistent (fct_plus, fct_minus) pair on
     * each side of a partition boundary, breaking flux limiting and causing
     * tracer overshoot → density gradient → momentum blowup within a few
     * steps. Two single-field calls (no 2-field variant in our halo). */
    if (sc->partit && sc->partit->npes > 1) {
        fesom_exchange_nod3D(fct_plus,  nl, sc->partit);
        fesom_exchange_nod3D(fct_minus, nl, sc->partit);
    }

    /* b3. Apply limits — vertical first, then horizontal. */
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        if (nu1 >= nl1) continue;
        /* surface interface (nz = nu1) */
        {
            int nz = nu1;
            real_t f = adf_v[FESOM_NODE3D(n, nz, nl)];
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t ae = (f >= 0.0) ? fct_plus[k] : fct_minus[k];
            adf_v[FESOM_NODE3D(n, nz, nl)] = ae * f;
        }
        /* interior interfaces */
        for (int nz = nu1 + 1; nz < nl1; ++nz) {
            real_t f = adf_v[FESOM_NODE3D(n, nz, nl)];
            real_t ae;
            if (f >= 0.0) {
                real_t a = fct_minus[FESOM_NODE3D(n, nz - 1, nl)];
                real_t b = fct_plus [FESOM_NODE3D(n, nz,     nl)];
                ae = (a < b) ? a : b;
            } else {
                real_t a = fct_plus [FESOM_NODE3D(n, nz - 1, nl)];
                real_t b = fct_minus[FESOM_NODE3D(n, nz,     nl)];
                ae = (a < b) ? a : b;
            }
            adf_v[FESOM_NODE3D(n, nz, nl)] = ae * f;
        }
        /* bottom flux is 0 by construction */
    }
    /* Horizontal */
    for (int e = 0; e < Eedg; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;

        int nl12 = (nl1 > nl2) ? nl1 : nl2;
        int nu12 = nu1;
        if (el2 >= 0 && nu2 < nu12) nu12 = nu2;

        for (int nz = nu12; nz < nl12; ++nz) {
            real_t f = adf_h[FESOM_NODE3D(e, nz, nl)];
            real_t ae;
            if (f >= 0.0) {
                real_t a = fct_plus [FESOM_NODE3D(n1, nz, nl)];
                real_t b = fct_minus[FESOM_NODE3D(n2, nz, nl)];
                ae = (a < b) ? a : b;
            } else {
                real_t a = fct_minus[FESOM_NODE3D(n1, nz, nl)];
                real_t b = fct_plus [FESOM_NODE3D(n2, nz, nl)];
                ae = (a < b) ? a : b;
            }
            adf_h[FESOM_NODE3D(e, nz, nl)] = ae * f;
        }
    }
}

/*--- flux2dtracer with use_lo=true (oce_adv_tra_driver.F90:451-572) -------
 * Adds the LO contribution to dttf_v then accumulates the (now-limited)
 * antidiffusive fluxes into dttf_h and dttf_v.
 *
 *   dttf_v[nz, n] += −ttf·hnode + LO·hnode_new          (LO transition)
 *   dttf_v[nz, n] += (flux_v[nz] − flux_v[nz+1]) · dt / areasvol
 *   dttf_h[nz, n1] += flux_h · dt / areasvol[n1]
 *   dttf_h[nz, n2] −= flux_h · dt / areasvol[n2]
 */
static void flux2dtracer_fct(const struct fesom_mesh *mesh,
                             const real_t            *ttf,
                             const real_t            *lo,
                             const real_t            *flux_h,
                             const real_t            *flux_v,
                             real_t                  *dttf_h,
                             real_t                  *dttf_v)
{
    const int N  = mesh->myDim_nod2D;
    const int E  = mesh->myDim_edge2D;   /* Fortran adv loops are myDim only */
    const int nl = mesh->nl;
    const real_t dt = (real_t)FESOM_PHASE1_DT;

    /* Vertical: LO transition + antidiffusive accumulation. */
    for (int n = 0; n < N; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nu1; nz < nl1; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            dttf_v[k] += -ttf[k] * mesh->hnode    [k]
                       +  lo[k] * mesh->hnode_new[k];
        }
        for (int nz = nu1; nz < nl1; ++nz) {
            real_t a = mesh->areasvol[FESOM_NODE3D(n, nz, nl)];
            if (a <= 0.0) continue;
            real_t f_top = flux_v[FESOM_NODE3D(n, nz,     nl)];
            real_t f_bot = flux_v[FESOM_NODE3D(n, nz + 1, nl)];
            dttf_v[FESOM_NODE3D(n, nz, nl)] += (f_top - f_bot) * dt / a;
        }
    }

    /* Horizontal */
    for (int e = 0; e < E; ++e) {
        int n1  = mesh->edges[2*e + 0];
        int n2  = mesh->edges[2*e + 1];
        int el1 = mesh->edge_tri[2*e + 0];
        int el2 = mesh->edge_tri[2*e + 1];

        int nl1 = (el1 >= 0) ? mesh->nlevels[el1] - 1 : 0;
        int nu1 = (el1 >= 0) ? mesh->ulevels[el1] - 1 : 0;
        int nl2 = (el2 >= 0) ? mesh->nlevels[el2] - 1 : 0;
        int nu2 = (el2 >= 0) ? mesh->ulevels[el2] - 1 : 0;

        int nl12 = (nl1 > nl2) ? nl1 : nl2;
        int nu12 = nu1;
        if (el2 >= 0 && nu2 > 0 && nu2 < nu12) nu12 = nu2;

        for (int nz = nu12; nz < nl12; ++nz) {
            real_t f  = flux_h[FESOM_NODE3D(e, nz, nl)];
            real_t a1 = mesh->areasvol[FESOM_NODE3D(n1, nz, nl)];
            real_t a2 = mesh->areasvol[FESOM_NODE3D(n2, nz, nl)];
            if (a1 > 0.0) dttf_h[FESOM_NODE3D(n1, nz, nl)] += f * dt / a1;
            if (a2 > 0.0) dttf_h[FESOM_NODE3D(n2, nz, nl)] -= f * dt / a2;
        }
    }
}

/*--- Public FCT entry: full FCT step for one tracer ------------------------*/
void fesom_tracer_advect_one_fct(fesom_tracer_adv_scratch *sc,
                                 int                       tr_idx,
                                 const struct fesom_mesh  *mesh,
                                 const struct fesom_dyn   *dyn,
                                 struct fesom_tracers     *tracers)
{
    const int N    = mesh->myDim_nod2D;
    const int nl   = mesh->nl;

    /* 1. init_tracers_AB (zero del_ttf*, fill valuesAB, save valuesold) */
    init_tracers_AB_one(tr_idx, mesh, tracers, sc);

    real_t *vals  = tracers->data[tr_idx].values;
    real_t *ttfAB = tracers->data[tr_idx].valuesAB;

    /* 2. LO upwind fluxes from current `values` (NOT valuesAB).
       Match Fortran do_oce_adv_tra:115 which calls hor_upw1 with `ttf`. */
    adv_tra_hor_upw1(mesh, dyn, vals, sc->adv_flux_hor);
    adv_tra_ver_upw1(mesh, dyn, vals, sc->adv_flux_ver);

    /* 3. fct_LO from upwind (advance one step of upwind solution). */
    fesom_tracer_compute_fct_LO(sc, tr_idx, mesh, tracers);

    /* Halo-exchange fct_LO — Fortran oce_adv_tra_driver.F90:294
     *   call exchange_nod(fct_LO, partit, luse_g2g=.true.)
     * Step 5 (oce_tra_adv_fct) reads LO at halo nodes when computing
     * T_max/T_min for the FCT bounds; without this exchange, halo LO is
     * stale → wrong fct_plus/fct_minus → boundary tracer drift. */
    if (sc->partit && sc->partit->npes > 1) {
        fesom_exchange_nod3D(sc->fct_LO, mesh->nl, sc->partit);
    }

    /* 4. High-order flux on top of LO (init_zero=false → adv_flux := HO − LO).
       Horizontal: MFCT 3rd-order (num_ord=0, CORE2 namelist.tra).
       FIDELITY: the MFCT flux mixes two fields exactly as the Fortran does —
         * node values + central diff  → valuesAB (ttfAB), passed to adv_tra_hor_mfct
           (oce_adv_tra_driver.F90:307 calls it with ttfAB).
         * up/down gradient correction → edge_up_dn_grad, built ONCE in
           init_tracers_AB from `values` (current step), NOT valuesAB — the
           valuesAB variant at oce_tracer_mod.F90:126 is commented out, :127 uses
           values. So tr_xy/edge_up_dn_grad must be built from vals, not ttfAB. */
    tracer_gradient_elements(mesh, vals, sc->tr_xy, sc->partit);
    fill_up_dn_grad         (mesh, sc->tr_xy, sc->edge_up_dn_grad);
    adv_tra_hor_mfct(mesh, dyn, ttfAB, sc->edge_up_dn_grad, /*num_ord=*/0.0,
                     /*init_zero=*/0, sc->adv_flux_hor);
    adv_tra_ver_qr4c   (mesh, dyn, ttfAB, /*num_ord=*/1.0, /*init_zero=*/0, sc->adv_flux_ver);

    /* 5. Zalesak limit the antidiffusive fluxes. */
    oce_tra_adv_fct(sc, mesh, vals, sc->adv_flux_hor, sc->adv_flux_ver);

    /* 6. flux2dtracer with LO (include the LO transition + limited HO). */
    flux2dtracer_fct(mesh, vals, sc->fct_LO,
                     sc->adv_flux_hor, sc->adv_flux_ver,
                     sc->del_ttf_advhoriz, sc->del_ttf_advvert);

    /* 7. del_ttf := del_ttf_advhoriz + del_ttf_advvert */
    for (int n = 0; n < N; ++n) {
        for (int nz = 0; nz < nl; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            tracers->del_ttf[k] = sc->del_ttf_advhoriz[k]
                                + sc->del_ttf_advvert [k];
        }
    }

    /* 8. ALE reconstruction. With FCT use_lo=true the math reduces to
       T_new = LO + limited_antidiff/areasvol/hnode_new — see fesom_tracer_adv.h */
    ale_reconstruct(mesh, vals, tracers->del_ttf);
}

/*--- Public entry: one full upwind step for one tracer --------------------*/
void fesom_tracer_advect_one(fesom_tracer_adv_scratch *sc,
                             int                       tr_idx,
                             const struct fesom_mesh  *mesh,
                             const struct fesom_dyn   *dyn,
                             struct fesom_tracers     *tracers)
{
    const int N    = mesh->myDim_nod2D;
    const int nl   = mesh->nl;

    /* 1. init_tracers_AB */
    init_tracers_AB_one(tr_idx, mesh, tracers, sc);

    /* 2-3. compute upwind fluxes from valuesAB (Fortran calls UPW1 with ttfAB) */
    real_t *ttfAB = tracers->data[tr_idx].valuesAB;
    adv_tra_hor_upw1(mesh, dyn, ttfAB, sc->adv_flux_hor);
    adv_tra_ver_upw1(mesh, dyn, ttfAB, sc->adv_flux_ver);

    /* 4. accumulate fluxes into del_ttf_advhoriz, del_ttf_advvert */
    flux2dtracer_upwind(mesh, sc->adv_flux_hor, sc->adv_flux_ver,
                        sc->del_ttf_advhoriz, sc->del_ttf_advvert);

    /* 5. del_ttf := del_ttf_advhoriz + del_ttf_advvert (lines 239-242) */
    for (int n = 0; n < N; ++n) {
        for (int nz = 0; nz < nl; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            tracers->del_ttf[k] = sc->del_ttf_advhoriz[k]
                                + sc->del_ttf_advvert [k];
        }
    }

    /* 6. ALE reconstruction */
    ale_reconstruct(mesh, tracers->data[tr_idx].values, tracers->del_ttf);
}
