/*
 * Phase 1 momentum: literal port of compute_vel_rhs, impl_vert_visc_ale,
 * update_vel, compute_hbar_ale. Each routine follows the Fortran source
 * line-by-line; deferred branches are explicitly listed in fesom_momentum.h.
 */
#include "fesom_momentum.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_mesh.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================
 * compute_vel_rhs (oce_ale_vel_rhs.F90:35-328)
 *
 * Three sub-loops (we collapse the second AB-shift loop into the first since
 * Phase 1 uses AB_order=2, single-slot history):
 *
 * Per-element loop:
 *   ab1 = -(0.5 + epsilon),  ab2 = (1.5 + epsilon)        (lines 97-100)
 *   For each layer:
 *       UV_rhs[u] = ab1 * UV_rhsAB[u]      ← previous step's Coriolis
 *       UV_rhs[v] = ab1 * UV_rhsAB[v]
 *
 *   ff = coriolis(elem) * elem_area(elem)                  (line 157)
 *   p_eta = g * eta_n[3 vertices]                          (line 155)
 *   pre   = -p_eta                                         (Phase 1: no ice/air/tides)
 *   Fx, Fy = sum( gradient_sca[1..3 / 4..6] * pre )        (lines 188-189)
 *
 *   For each layer:
 *       UV_rhs[u] += (Fx - pgf_x) * elem_area              (lines 197-198)
 *       UV_rhs[v] += (Fy - pgf_y) * elem_area
 *       UV_rhsAB[u] =  v * ff      ← new Coriolis stored for AB next step
 *       UV_rhsAB[v] = -u * ff                              (lines 203-204)
 *
 * Final assembly loop (after the per-element loop completes):
 *   ff_step = ab2     (1.0 on the very first step, to match Fortran's lfirst)
 *   For each layer:
 *       UV_rhs[u] = dt * (UV_rhs[u] + UV_rhsAB[u] * ff_step) / elem_area
 *       UV_rhs[v] = dt * (UV_rhs[v] + UV_rhsAB[v] * ff_step) / elem_area
 *===========================================================================*/

void fesom_compute_vel_rhs(const struct fesom_mesh *mesh,
                           const struct fesom_aux  *aux,
                           struct fesom_dyn        *dyn,
                           int                      is_first_step)
{
    const int E       = mesh->myDim_elem2D;     /* interior only — halo gets exchange */
    const int nl      = mesh->nl;
    const real_t g    = (real_t)FESOM_G;
    const real_t dt   = (real_t)FESOM_PHASE1_DT;

    /* Adams-Bashforth coefficients (oce_ale_vel_rhs.F90:97-100). The
       small `epsilon` matches Fortran's `epsilon` in MOD_PARAM (=1e-9 by
       convention; FESOM uses it to keep the scheme stable at the first
       step). For literal byte-parity we use the same numeric value. */
    const real_t eps = 1.0e-9;
    const real_t ab1 = -(0.5 + eps);
    const real_t ab2 =  (1.5 + eps);

    /* Per-element computation: (i) shift AB history into UV_rhs,
       (ii) add SSH-gradient + PGF, (iii) overwrite UV_rhsAB with this step's
       Coriolis term. */
    for (int e = 0; e < E; ++e) {
        int nzmin = mesh->ulevels[e]   - 1;       /* 1-based → 0-based */
        int nzmax = mesh->nlevels[e]   - 1;
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];

        /* (i) AB history shift */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t k = FESOM_ELEMVEC(e, nz, nl);
            dyn->uv_rhs[k + 0] = ab1 * dyn->uv_rhsAB[k + 0];
            dyn->uv_rhs[k + 1] = ab1 * dyn->uv_rhsAB[k + 1];
        }

        /* SSH gradient (Phase 1: no ice / air / tides — pre = -p_eta) */
        real_t pre0 = -g * dyn->eta_n[n0];
        real_t pre1 = -g * dyn->eta_n[n1];
        real_t pre2 = -g * dyn->eta_n[n2];
        const real_t *gs = &mesh->gradient_sca[6 * e];
        real_t Fx = gs[0]*pre0 + gs[1]*pre1 + gs[2]*pre2;
        real_t Fy = gs[3]*pre0 + gs[4]*pre1 + gs[5]*pre2;

        real_t area = mesh->elem_area[e];
        real_t ff   = mesh->coriolis[e] * area;       /* Coriolis * area */

        /* (ii) and (iii) */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t k    = FESOM_ELEMVEC(e, nz, nl);
            size_t kel  = FESOM_ELEM3D(e, nz, nl);
            real_t pgfx = aux->pgf_x[kel];
            real_t pgfy = aux->pgf_y[kel];
            dyn->uv_rhs[k + 0] += (Fx - pgfx) * area;
            dyn->uv_rhs[k + 1] += (Fy - pgfy) * area;
            real_t u = dyn->uv[k + 0];
            real_t v = dyn->uv[k + 1];
            dyn->uv_rhsAB[k + 0] =  v * ff;
            dyn->uv_rhsAB[k + 1] = -u * ff;
        }
    }

    /* Final assembly with ff_step. Fortran tracks lfirst as a `save` flag —
       caller passes is_first_step here. */
    real_t ff_step = is_first_step ? 1.0 : ab2;
    for (int e = 0; e < E; ++e) {
        int nzmin = mesh->ulevels[e]   - 1;
        int nzmax = mesh->nlevels[e]   - 1;
        real_t inv_area = 1.0 / mesh->elem_area[e];
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t k = FESOM_ELEMVEC(e, nz, nl);
            dyn->uv_rhs[k + 0] = dt * (dyn->uv_rhs[k + 0]
                                     + dyn->uv_rhsAB[k + 0] * ff_step) * inv_area;
            dyn->uv_rhs[k + 1] = dt * (dyn->uv_rhs[k + 1]
                                     + dyn->uv_rhsAB[k + 1] * ff_step) * inv_area;
        }
    }
}

/*===========================================================================
 * impl_vert_visc_ale (oce_ale.F90:3124-3335)
 *
 * Per-element TDMA on (u, v) columns:
 *   - Build tridiagonal coefficients a, b, c using Av and helem-derived dz
 *   - Add advective contribution from Wvel_i (zero in Phase 1; see fesom_dyn.h)
 *   - Surface row: add wind stress / density_0 * zinv to RHS
 *   - Bottom row: add bottom drag friction * UV * zinv to RHS
 *   - Convert RHS from "u_new" formulation to "du = u_new - u_old"
 *   - Thomas (sweep) algorithm
 *   - Write result back into UV_rhs
 *
 * For Phase 1 (no partial cells), zbar_n[nz] = zbar[nz] for all elements,
 * Z_n[nz] = Z[nz] = 0.5*(zbar[nz]+zbar[nz+1]). But Fortran rebuilds zbar_n /
 * Z_n per element from helem (so partial-cell-aware). We do the same to keep
 * the code structurally identical.
 *===========================================================================*/

void fesom_impl_vert_visc(const struct fesom_mesh    *mesh,
                          const struct fesom_aux     *aux,
                          const struct fesom_forcing *forcing,
                          struct fesom_dyn           *dyn)
{
    const int E      = mesh->myDim_elem2D;
    const int nl     = mesh->nl;
    const real_t dt  = (real_t)FESOM_PHASE1_DT;
    const real_t Cd  = (real_t)FESOM_PHASE1_C_D;
    const real_t inv_density_0 = 1.0 / (real_t)FESOM_DENSITY_0;

    /* Per-column scratch — nl ≤ 48 so static cap of 64 is safe. */
    enum { NL_MAX = 64 };

    for (int e = 0; e < E; ++e) {
        int nzmin = mesh->ulevels[e]   - 1;
        int nzmax = mesh->nlevels[e]   - 1;
        if (nzmax - nzmin < 1) continue;
        FESOM_CHECK(nzmax + 1 < NL_MAX, "impl_vert_visc: nzmax %d > NL_MAX", nzmax);

        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];

        real_t zbar_n[NL_MAX], Z_n[NL_MAX];
        real_t a[NL_MAX], b[NL_MAX], c[NL_MAX];
        real_t ur[NL_MAX], vr[NL_MAX];
        real_t cp[NL_MAX], up[NL_MAX], vp[NL_MAX];

        /* Build zbar_n and Z_n from helem upward — mirror of lines 3167-3179.
           In 0-based C: nzmax is the bottom-interface index (Fortran nlevels-1).
           helem is layer thickness; zbar_n[nz] is the interface above layer nz. */
        for (int k = 0; k <= nzmax; ++k) zbar_n[k] = 0.0;
        for (int k = 0; k < nzmax; ++k)  Z_n   [k] = 0.0;
        zbar_n[nzmax]   = mesh->zbar[nzmax];   /* in Phase 1 = zbar_e_bot[e] */
        Z_n   [nzmax-1] = zbar_n[nzmax]
                        + mesh->helem[FESOM_ELEM3D(e, nzmax - 1, nl)] / 2.0;
        for (int nz = nzmax - 1; nz >= nzmin + 1; --nz) {
            zbar_n[nz]   = zbar_n[nz + 1] + mesh->helem[FESOM_ELEM3D(e, nz,     nl)];
            Z_n   [nz-1] = zbar_n[nz]     + mesh->helem[FESOM_ELEM3D(e, nz - 1, nl)] / 2.0;
        }
        zbar_n[nzmin] = zbar_n[nzmin + 1] + mesh->helem[FESOM_ELEM3D(e, nzmin, nl)];

        /* Tridiagonal coefficients: regular interior (lines 3185-3198) */
        for (int nz = nzmin + 1; nz <= nzmax - 2; ++nz) {
            real_t zinv = dt / (zbar_n[nz] - zbar_n[nz + 1]);
            a[nz] = -aux->Av[FESOM_ELEM3D(e, nz,     nl)] / (Z_n[nz - 1] - Z_n[nz]) * zinv;
            c[nz] = -aux->Av[FESOM_ELEM3D(e, nz + 1, nl)] / (Z_n[nz] - Z_n[nz + 1]) * zinv;
            b[nz] = -a[nz] - c[nz] + 1.0;
            /* Vertical advection update from Wvel_i (zero in Phase 1) */
            real_t wu = (dyn->w_i[FESOM_NODE3D(n0, nz,     nl)]
                       + dyn->w_i[FESOM_NODE3D(n1, nz,     nl)]
                       + dyn->w_i[FESOM_NODE3D(n2, nz,     nl)]) / 3.0;
            real_t wd = (dyn->w_i[FESOM_NODE3D(n0, nz + 1, nl)]
                       + dyn->w_i[FESOM_NODE3D(n1, nz + 1, nl)]
                       + dyn->w_i[FESOM_NODE3D(n2, nz + 1, nl)]) / 3.0;
            a[nz] = a[nz] + (wu < 0.0 ? wu : 0.0) * zinv;
            b[nz] = b[nz] + (wu > 0.0 ? wu : 0.0) * zinv;
            b[nz] = b[nz] - (wd < 0.0 ? wd : 0.0) * zinv;
            c[nz] = c[nz] - (wd > 0.0 ? wd : 0.0) * zinv;
        }

        /* Bottom row (lines 3200-3208) */
        {
            int nz = nzmax - 1;
            real_t zinv = dt / (zbar_n[nz] - zbar_n[nz + 1]);
            a[nz] = -aux->Av[FESOM_ELEM3D(e, nz, nl)] / (Z_n[nz - 1] - Z_n[nz]) * zinv;
            b[nz] = -a[nz] + 1.0;
            c[nz] = 0.0;
            real_t wu = (dyn->w_i[FESOM_NODE3D(n0, nz, nl)]
                       + dyn->w_i[FESOM_NODE3D(n1, nz, nl)]
                       + dyn->w_i[FESOM_NODE3D(n2, nz, nl)]) / 3.0;
            a[nz] = a[nz] + (wu < 0.0 ? wu : 0.0) * zinv;
            b[nz] = b[nz] + (wu > 0.0 ? wu : 0.0) * zinv;
        }

        /* Surface row (lines 3215-3231); zinv kept for the wind-stress RHS use */
        real_t zinv_top = dt / (zbar_n[nzmin] - zbar_n[nzmin + 1]);
        {
            int nz = nzmin;
            c[nz] = -aux->Av[FESOM_ELEM3D(e, nz + 1, nl)] / (Z_n[nz] - Z_n[nz + 1]) * zinv_top;
            a[nz] = 0.0;
            b[nz] = -c[nz] + 1.0;
            real_t wu = (dyn->w_i[FESOM_NODE3D(n0, nz,     nl)]
                       + dyn->w_i[FESOM_NODE3D(n1, nz,     nl)]
                       + dyn->w_i[FESOM_NODE3D(n2, nz,     nl)]) / 3.0;
            real_t wd = (dyn->w_i[FESOM_NODE3D(n0, nz + 1, nl)]
                       + dyn->w_i[FESOM_NODE3D(n1, nz + 1, nl)]
                       + dyn->w_i[FESOM_NODE3D(n2, nz + 1, nl)]) / 3.0;
            b[nz] = b[nz] + wu * zinv_top;
            b[nz] = b[nz] - (wd < 0.0 ? wd : 0.0) * zinv_top;
            c[nz] = c[nz] - (wd > 0.0 ? wd : 0.0) * zinv_top;
        }

        /* Build RHS from current uv_rhs (lines 3238-3245) */
        for (int nz = nzmin; nz <= nzmax - 1; ++nz) {
            ur[nz] = dyn->uv_rhs[FESOM_ELEMVEC(e, nz, nl) + 0];
            vr[nz] = dyn->uv_rhs[FESOM_ELEMVEC(e, nz, nl) + 1];
        }

        /* Surface forcing (wind stress); element wind stress at index 2*e */
        ur[nzmin] += zinv_top * forcing->stress_surf[2*e + 0] * inv_density_0;
        vr[nzmin] += zinv_top * forcing->stress_surf[2*e + 1] * inv_density_0;

        /* Bottom drag (lines 3267-3272). Phase 1: not toy_ocean, so quadratic
           drag from the bottom-layer velocity magnitude. */
        real_t zinv_bot = dt / (zbar_n[nzmax - 1] - zbar_n[nzmax]);
        {
            int nz = nzmax - 1;
            real_t u_bot = dyn->uv[FESOM_ELEMVEC(e, nz, nl) + 0];
            real_t v_bot = dyn->uv[FESOM_ELEMVEC(e, nz, nl) + 1];
            real_t spd   = sqrt(u_bot * u_bot + v_bot * v_bot);
            real_t friction = -Cd * spd;
            ur[nz] += zinv_bot * friction * u_bot;
            vr[nz] += zinv_bot * friction * v_bot;
        }

        /* Convert "solve for u_new" to "solve for du = u_new - u_old"
           (lines 3282-3292). */
        for (int nz = nzmin + 1; nz <= nzmax - 2; ++nz) {
            real_t u_up = dyn->uv[FESOM_ELEMVEC(e, nz - 1, nl) + 0];
            real_t v_up = dyn->uv[FESOM_ELEMVEC(e, nz - 1, nl) + 1];
            real_t u_th = dyn->uv[FESOM_ELEMVEC(e, nz,     nl) + 0];
            real_t v_th = dyn->uv[FESOM_ELEMVEC(e, nz,     nl) + 1];
            real_t u_dn = dyn->uv[FESOM_ELEMVEC(e, nz + 1, nl) + 0];
            real_t v_dn = dyn->uv[FESOM_ELEMVEC(e, nz + 1, nl) + 1];
            ur[nz] -= a[nz]*u_up + (b[nz] - 1.0)*u_th + c[nz]*u_dn;
            vr[nz] -= a[nz]*v_up + (b[nz] - 1.0)*v_th + c[nz]*v_dn;
        }
        {
            int nz = nzmin;
            real_t u_th = dyn->uv[FESOM_ELEMVEC(e, nz,     nl) + 0];
            real_t v_th = dyn->uv[FESOM_ELEMVEC(e, nz,     nl) + 1];
            real_t u_dn = dyn->uv[FESOM_ELEMVEC(e, nz + 1, nl) + 0];
            real_t v_dn = dyn->uv[FESOM_ELEMVEC(e, nz + 1, nl) + 1];
            ur[nz] -= (b[nz] - 1.0)*u_th + c[nz]*u_dn;
            vr[nz] -= (b[nz] - 1.0)*v_th + c[nz]*v_dn;
        }
        {
            int nz = nzmax - 1;
            real_t u_up = dyn->uv[FESOM_ELEMVEC(e, nz - 1, nl) + 0];
            real_t v_up = dyn->uv[FESOM_ELEMVEC(e, nz - 1, nl) + 1];
            real_t u_th = dyn->uv[FESOM_ELEMVEC(e, nz,     nl) + 0];
            real_t v_th = dyn->uv[FESOM_ELEMVEC(e, nz,     nl) + 1];
            ur[nz] -= a[nz]*u_up + (b[nz] - 1.0)*u_th;
            vr[nz] -= a[nz]*v_up + (b[nz] - 1.0)*v_th;
        }

        /* Thomas algorithm forward sweep (lines 3301-3312) */
        cp[nzmin] = c[nzmin] / b[nzmin];
        up[nzmin] = ur[nzmin] / b[nzmin];
        vp[nzmin] = vr[nzmin] / b[nzmin];
        for (int nz = nzmin + 1; nz <= nzmax - 1; ++nz) {
            real_t m = b[nz] - cp[nz - 1] * a[nz];
            cp[nz] = c[nz] / m;
            up[nz] = (ur[nz] - up[nz - 1] * a[nz]) / m;
            vp[nz] = (vr[nz] - vp[nz - 1] * a[nz]) / m;
        }
        /* Back substitution (lines 3314-3322) */
        ur[nzmax - 1] = up[nzmax - 1];
        vr[nzmax - 1] = vp[nzmax - 1];
        for (int nz = nzmax - 2; nz >= nzmin; --nz) {
            ur[nz] = up[nz] - cp[nz] * ur[nz + 1];
            vr[nz] = vp[nz] - cp[nz] * vr[nz + 1];
        }

        /* Write back into uv_rhs (lines 3328-3331) */
        for (int nz = nzmin; nz <= nzmax - 1; ++nz) {
            dyn->uv_rhs[FESOM_ELEMVEC(e, nz, nl) + 0] = ur[nz];
            dyn->uv_rhs[FESOM_ELEMVEC(e, nz, nl) + 1] = vr[nz];
        }
    }
}

/*===========================================================================
 * update_vel (oce_dyn.F90:73-173)
 *
 *   eta = -g * theta * dt * d_eta[3 vertices]
 *   Fx = sum( gradient_sca[1..3] * eta )
 *   Fy = sum( gradient_sca[4..6] * eta )
 *   uv += uv_rhs + (Fx, Fy)
 *===========================================================================*/

void fesom_update_vel(const struct fesom_mesh *mesh,
                      struct fesom_dyn        *dyn)
{
    const int E = mesh->myDim_elem2D;       /* interior — exchange UV after */
    const int nl = mesh->nl;
    const real_t coef = -(real_t)FESOM_G * (real_t)FESOM_PHASE1_THETA
                       * (real_t)FESOM_PHASE1_DT;

    for (int e = 0; e < E; ++e) {
        int nzmin = mesh->ulevels[e]   - 1;
        int nzmax = mesh->nlevels[e]   - 1;
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];

        real_t e0 = coef * dyn->d_eta[n0];
        real_t e1 = coef * dyn->d_eta[n1];
        real_t e2 = coef * dyn->d_eta[n2];
        const real_t *g = &mesh->gradient_sca[6 * e];
        real_t Fx = g[0]*e0 + g[1]*e1 + g[2]*e2;
        real_t Fy = g[3]*e0 + g[4]*e1 + g[5]*e2;

        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t k = FESOM_ELEMVEC(e, nz, nl);
            dyn->uv[k + 0] += dyn->uv_rhs[k + 0] + Fx;
            dyn->uv[k + 1] += dyn->uv_rhs[k + 1] + Fy;
        }
    }
}

/*===========================================================================
 * visc_filt_bcksct (oce_dyn.F90:278-426) — opt_visc=5 with easy backscatter
 *===========================================================================*/

void fesom_visc_filt_bcksct(const struct fesom_mesh *mesh,
                            struct fesom_dyn        *dyn)
{
    const int E    = mesh->myDim_elem2D;
    const int Eedg = mesh->myDim_edge2D + mesh->eDim_edge2D;
    const int nl   = mesh->nl;
    /* Allocations sized for full local extent — touch only the elements we
     * compute on (exchange is done in fesom_step after this returns). */
    int E_alloc = mesh->myDim_elem2D + mesh->eDim_elem2D + mesh->eXDim_elem2D;
    int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;
    (void)E_alloc; (void)N_alloc;
    const real_t dt = (real_t)FESOM_PHASE1_DT;
    const real_t g0 = (real_t)FESOM_PHASE1_VISC_GAMMA0;
    const real_t g1 = (real_t)FESOM_PHASE1_VISC_GAMMA1;
    const real_t g2 = (real_t)FESOM_PHASE1_VISC_GAMMA2;
    const real_t bsret = (real_t)FESOM_PHASE1_VISC_EASYBSRETURN;

    /* Zero the four work arrays (lines 316-322) over full local extent. */
    memset(dyn->u_b, 0, (size_t)E_alloc * (size_t)nl * sizeof(real_t));
    memset(dyn->v_b, 0, (size_t)E_alloc * (size_t)nl * sizeof(real_t));
    memset(dyn->u_c, 0, (size_t)N_alloc * (size_t)nl * sizeof(real_t));
    memset(dyn->v_c, 0, (size_t)N_alloc * (size_t)nl * sizeof(real_t));

    /* Edge loop — accumulate per-element viscous updates (lines 327-361).
       Skip boundary edges (Fortran: myList_edge2D > edge2D_in; serial: el2 < 0). */
    for (int ed = 0; ed < Eedg; ++ed) {
        int el1 = mesh->edge_tri[2*ed + 0];
        int el2 = mesh->edge_tri[2*ed + 1];
        if (el1 < 0 || el2 < 0) continue;       /* boundary edge */

        real_t a1 = mesh->elem_area[el1];
        real_t a2 = mesh->elem_area[el2];
        real_t len = sqrt(a1 + a2);

        int nu1 = mesh->ulevels[el1] - 1;
        int nu2 = mesh->ulevels[el2] - 1;
        int nzmin = (nu1 > nu2) ? nu1 : nu2;       /* maxval(ulevels) */
        int nl1   = mesh->nlevels[el1] - 1;
        int nl2   = mesh->nlevels[el2] - 1;
        int nzmax = (nl1 < nl2) ? nl1 : nl2;       /* minval(nlevels) - 1 */

        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t u1 = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0]
                       - dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
            real_t v1 = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1]
                       - dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
            real_t mag2 = u1*u1 + v1*v1;
            real_t mag  = sqrt(mag2);
            real_t vi_inner = (g1 * mag > g2 * mag2) ? g1 * mag : g2 * mag2;
            real_t vi = dt * (g0 > vi_inner ? g0 : vi_inner) * len;
            real_t du = u1 * vi;
            real_t dv = v1 * vi;
            dyn->u_b[FESOM_ELEM3D(el1, nz, nl)] -= du / a1;
            dyn->v_b[FESOM_ELEM3D(el1, nz, nl)] -= dv / a1;
            dyn->u_b[FESOM_ELEM3D(el2, nz, nl)] += du / a2;
            dyn->v_b[FESOM_ELEM3D(el2, nz, nl)] += dv / a2;
        }
    }

    /* Smooth U_b → U_c at nodes (area-weighted mean of surrounding cells).
       Interior nodes only; halo gets values via exchange afterward.
       Same neighbour traversal as compute_vel_nodes. */
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        int nu1 = mesh->ulevels_nod2D[n] - 1;
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        int o0 = mesh->nod_in_elem2D_offsets[n];
        int o1 = mesh->nod_in_elem2D_offsets[n + 1];
        for (int nz = nu1; nz < nl1; ++nz) {
            real_t tot = 0.0, tu = 0.0, tv = 0.0;
            for (int k = o0; k < o1; ++k) {
                int e = mesh->nod_in_elem2D[k];
                real_t a = mesh->elem_area[e];
                tot += a;
                tu  += dyn->u_b[FESOM_ELEM3D(e, nz, nl)] * a;
                tv  += dyn->v_b[FESOM_ELEM3D(e, nz, nl)] * a;
            }
            if (tot > 0.0) {
                dyn->u_c[FESOM_NODE3D(n, nz, nl)] = tu / tot;
                dyn->v_c[FESOM_NODE3D(n, nz, nl)] = tv / tot;
            }
        }
    }

    /* Apply: uv_rhs += u_b − easybsreturn · mean(u_c at 3 vertices)
       (lines 418-422 — non-split-explicit branch, what we use in Phase 2). */
    for (int e = 0; e < E; ++e) {
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];
        int nu1 = mesh->ulevels[e] - 1;
        int nle = mesh->nlevels[e] - 1;
        for (int nz = nu1; nz < nle; ++nz) {
            real_t uc_mean = (dyn->u_c[FESOM_NODE3D(n0, nz, nl)]
                            + dyn->u_c[FESOM_NODE3D(n1, nz, nl)]
                            + dyn->u_c[FESOM_NODE3D(n2, nz, nl)]) / 3.0;
            real_t vc_mean = (dyn->v_c[FESOM_NODE3D(n0, nz, nl)]
                            + dyn->v_c[FESOM_NODE3D(n1, nz, nl)]
                            + dyn->v_c[FESOM_NODE3D(n2, nz, nl)]) / 3.0;
            dyn->uv_rhs[FESOM_ELEMVEC(e, nz, nl) + 0]
                += dyn->u_b[FESOM_ELEM3D(e, nz, nl)] - bsret * uc_mean;
            dyn->uv_rhs[FESOM_ELEMVEC(e, nz, nl) + 1]
                += dyn->v_b[FESOM_ELEM3D(e, nz, nl)] - bsret * vc_mean;
        }
    }
}

/*===========================================================================
 * compute_hbar_ale (oce_ale.F90:1974-2116)
 *
 *   ssh_rhs_old = 0
 *   for each edge:
 *       c1 = Σ_nz (v[el1]*dx1 - u[el1]*dy1) * helem[el1, nz]
 *       c2 = -Σ_nz (v[el2]*dx2 - u[el2]*dy2) * helem[el2, nz]   (if el2 ≥ 0)
 *       ssh_rhs_old[n1] += c1 + c2
 *       ssh_rhs_old[n2] -= c1 + c2
 *   (linfs branch skips the water_flux subtraction at lines 2071-2080.)
 *
 *   hbar_old = hbar
 *   hbar = hbar_old + ssh_rhs_old * dt / areasvol[ulevels_nod2D[n], n]
 *===========================================================================*/

void fesom_compute_hbar(const struct fesom_mesh *mesh,
                        struct fesom_dyn        *dyn)
{
    const int N  = mesh->myDim_nod2D;
    const int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;
    /* Fortran oce_ale.F90:2014 uses `do ed=1, myDim_edge2D` only.
     * Iterating eDim_edge2D too would double-count contributions to interior
     * boundary nodes (the same physical edge is processed by the owner's
     * myDim loop too) → small per-step drift → CG NaN after ~85-95 steps. */
    const int E  = mesh->myDim_edge2D;
    const int nl = mesh->nl;
    const real_t dt = (real_t)FESOM_PHASE1_DT;

    /* Reset ssh_rhs_old (lines 2005-2009) over full local extent */
    memset(dyn->ssh_rhs_old, 0, (size_t)N_alloc * sizeof(real_t));

    /* Edge loop — accumulate transport divergence (lines 2014-2065) */
    for (int ed = 0; ed < E; ++ed) {
        int n1 = mesh->edges[2*ed + 0];
        int n2 = mesh->edges[2*ed + 1];
        int el1 = mesh->edge_tri[2*ed + 0];
        int el2 = mesh->edge_tri[2*ed + 1];

        real_t c1 = 0.0;
        if (el1 >= 0) {
            real_t dx1 = mesh->edge_cross_dxdy[4*ed + 0];
            real_t dy1 = mesh->edge_cross_dxdy[4*ed + 1];
            int nzmin = mesh->ulevels[el1] - 1;
            int nzmax = mesh->nlevels[el1] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                real_t u = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
                real_t v = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
                real_t h = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
                c1 += (v * dx1 - u * dy1) * h;
            }
        }
        real_t c2 = 0.0;
        if (el2 >= 0) {
            real_t dx2 = mesh->edge_cross_dxdy[4*ed + 2];
            real_t dy2 = mesh->edge_cross_dxdy[4*ed + 3];
            int nzmin = mesh->ulevels[el2] - 1;
            int nzmax = mesh->nlevels[el2] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                real_t u = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                c2 -= (v * dx2 - u * dy2) * h;
            }
        }
        dyn->ssh_rhs_old[n1] += (c1 + c2);
        dyn->ssh_rhs_old[n2] -= (c1 + c2);
    }

    /* hbar update (lines 2090-2102): linfs skips the water_flux term. */
    for (int n = 0; n < N; ++n) {
        mesh->hbar_old[n] = mesh->hbar[n];
    }
    for (int n = 0; n < N; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) continue;        /* cavity guard */
        int top_layer = mesh->ulevels_nod2D[n] - 1;       /* = 0 in Phase 1 */
        real_t area = mesh->areasvol[FESOM_NODE3D(n, top_layer, nl)];
        mesh->hbar[n] = mesh->hbar_old[n]
                       + dyn->ssh_rhs_old[n] * dt / area;
    }
    /* (Serial: no exchange_nod for ssh_rhs_old or hbar.) */
}
