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
#include "fesom_mesh.h"
#include "fesom_tracers.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void fesom_tracer_adv_init(fesom_tracer_adv_scratch *sc,
                           const struct fesom_mesh  *mesh)
{
    memset(sc, 0, sizeof(*sc));
    size_t e_full = (size_t)mesh->edge2D * (size_t)mesh->nl;   /* edge × nl (stride nl, last slot unused) */
    size_t n_full = (size_t)mesh->nod2D  * (size_t)mesh->nl;

    sc->adv_flux_hor     = calloc(e_full,         sizeof(real_t));
    sc->adv_flux_ver     = calloc(n_full,         sizeof(real_t));
    sc->del_ttf_advhoriz = calloc(n_full,         sizeof(real_t));
    sc->del_ttf_advvert  = calloc(n_full,         sizeof(real_t));
    sc->fct_LO           = calloc(n_full,         sizeof(real_t));
    sc->fct_ttf_min      = calloc(n_full,         sizeof(real_t));
    sc->fct_ttf_max      = calloc(n_full,         sizeof(real_t));
    sc->fct_plus         = calloc(n_full,         sizeof(real_t));
    sc->fct_minus        = calloc(n_full,         sizeof(real_t));
    sc->fct_aux          = calloc((size_t)mesh->elem2D * (size_t)mesh->nl * 2,
                                  sizeof(real_t));
    FESOM_CHECK(sc->adv_flux_hor && sc->adv_flux_ver
             && sc->del_ttf_advhoriz && sc->del_ttf_advvert
             && sc->fct_LO && sc->fct_ttf_min && sc->fct_ttf_max
             && sc->fct_plus && sc->fct_minus && sc->fct_aux,
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
    const int N    = mesh->nod2D;
    const int E    = mesh->edge2D;
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

/*--- init_tracers_AB (oce_tracer_mod.F90:13-145, AB2 path) ----------------*/
static void init_tracers_AB_one(int tr_idx,
                                const struct fesom_mesh *mesh,
                                struct fesom_tracers    *tracers,
                                fesom_tracer_adv_scratch *sc)
{
    const int N  = mesh->nod2D;
    const int nl = mesh->nl;
    const real_t eps = 1.0e-9;
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
    const int E    = mesh->edge2D;
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

/*--- adv_tra_ver_upw1 (oce_adv_tra_ver.F90:244-328) -----------------------*/
static void adv_tra_ver_upw1(const struct fesom_mesh *mesh,
                             const struct fesom_dyn  *dyn,
                             const real_t            *ttf,
                             real_t                  *flux)
{
    const int N  = mesh->nod2D;
    const int nl = mesh->nl;
    /* w_e = w in Phase 1 (no wsplit firing — see fesom_dyn.h). */
    const real_t *W = dyn->w;

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
    const int N  = mesh->nod2D;
    const int E  = mesh->edge2D;
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
    const int N  = mesh->nod2D;
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

/*--- Public entry: one full upwind step for one tracer --------------------*/
void fesom_tracer_advect_one(fesom_tracer_adv_scratch *sc,
                             int                       tr_idx,
                             const struct fesom_mesh  *mesh,
                             const struct fesom_dyn   *dyn,
                             struct fesom_tracers     *tracers)
{
    const int N    = mesh->nod2D;
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
