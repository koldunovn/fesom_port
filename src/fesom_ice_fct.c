/*
 * Sea-ice FCT advection. Literal port of ice_fct.F90.
 *
 * The Fortran nn_pos(k,row) / nn_num(row) lookup is reused via the SSH
 * stiffness CSR (stiff->rowptr / stiff->colind):
 *
 *   nn_pos(k, row)   ↔  stiff->colind[stiff->rowptr[row] + (k-1)]   (Fortran 1-based)
 *   nn_num(row)      ↔  stiff->rowptr[row+1] - stiff->rowptr[row]
 *   ssh_stiff%rowptr(row) - ssh_stiff%rowptr(1) + 1   ↔  stiff->rowptr[row]   (0-based start)
 *
 * fct_massmatrix shares the ssh_stiff sparsity pattern; allocated to
 * stiff->nnz here.
 *
 * Halo audit (per feedback_write_loops_halo.md and
 * feedback_array_size_vs_reader_loop.md). Each write-loop annotates its
 * Fortran source line and bound. The bound is preserved literally: do NOT
 * "fix" a Fortran myDim_nod2D loop into myDim+eDim — Fortran's downstream
 * readers stay at myDim too, and changing the bound would diverge.
 */

#include "fesom_ice_fct.h"
#include "fesom_constants.h"
#include "fesom_halo.h"
#include "fesom_ice.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_ssh.h"
#include "fesom_types.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FESOM_ICE_SCALE_AREA
#define FESOM_ICE_SCALE_AREA 2.0e8   /* oce_modules.F90:29 — global, used by ice_diff scaling */
#endif

/* Logical tracer index used by ice_fem_fct — mirrors Fortran tr_array_id
 * but 0-based here. The Fortran/C mapping:
 *   logical 0 (Fortran 1) → m_ice  (data[FESOM_ICE_MICE])
 *   logical 1 (Fortran 2) → a_ice  (data[FESOM_ICE_AICE])
 *   logical 2 (Fortran 3) → m_snow (data[FESOM_ICE_MSNOW])
 * Caller (fct_solve) calls in this exact order — same as Fortran. */
enum { FCT_LOG_M = 0, FCT_LOG_A = 1, FCT_LOG_MS = 2 };

static int data_idx_for_logical(int log_id)
{
    switch (log_id) {
        case FCT_LOG_M:  return FESOM_ICE_MICE;
        case FCT_LOG_A:  return FESOM_ICE_AICE;
        case FCT_LOG_MS: return FESOM_ICE_MSNOW;
        default:
            fprintf(stderr, "fesom_ice_fct: bad logical tracer id %d\n", log_id);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
    }
}

/* ============================================================ */
/* E1 — ice_mass_matrix_fill (ice_fct.F90:1145)                 */
/* ============================================================ */

void fesom_ice_mass_matrix_fill(fesom_ice                    *ice,
                                const struct fesom_ssh_stiff *stiff,
                                struct fesom_partit          *partit,
                                const struct fesom_mesh      *mesh)
{
    (void)partit;

    /* (Re-)allocate to ssh_stiff->nnz, zero-initialised. Fortran allocates
     * with the same shape via mass_matrix=0 then accumulates. */
    if (ice->work.fct_massmatrix == NULL) {
        ice->work.fct_massmatrix = calloc((size_t)stiff->nnz, sizeof(real_t));
        if (!ice->work.fct_massmatrix) {
            fprintf(stderr, "fesom_ice_mass_matrix_fill: out of memory (nnz=%d)\n",
                    stiff->nnz);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        memset(ice->work.fct_massmatrix, 0, (size_t)stiff->nnz * sizeof(real_t));
    }

    real_t *mm = ice->work.fct_massmatrix;
    const int N_own = mesh->myDim_nod2D;

    /* --- a) accumulate (ice_fct.F90:1170-1210) ------------------------------
     * Fortran outer DO elem=1, myDim_elem2D; inner DO n=1,3 (vertices);
     * skip if row > myDim_nod2D (halo); inner DO q=1,3 with cavity skip.
     * Bound: myDim_elem2D — matches Fortran exactly.                       */
    for (int elem = 0; elem < mesh->myDim_elem2D; ++elem) {
        const int en[3] = { mesh->elem_nodes[3*elem + 0],
                            mesh->elem_nodes[3*elem + 1],
                            mesh->elem_nodes[3*elem + 2] };
        for (int nl = 0; nl < 3; ++nl) {
            int row = en[nl];
            if (row >= N_own) continue;     /* skip halo row (Fortran line 1177) */
            int rstart = stiff->rowptr[row];
            int rend   = stiff->rowptr[row + 1];

            for (int q = 0; q < 3; ++q) {
                if (mesh->ulevels[elem] > 1) continue;   /* cavity (line 1186) */

                /* find ipos: position of en[q] in row's sparsity pattern */
                int ipos = -1;
                for (int k = rstart; k < rend; ++k) {
                    if (stiff->colind[k] == en[q]) { ipos = k; break; }
                }
                if (ipos < 0) {
                    fprintf(stderr,
                            "fesom_ice_mass_matrix_fill: FATAL — row=%d en[q=%d]=%d not in CSR\n",
                            row, q, en[q]);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                mm[ipos] += mesh->elem_area[elem] / 12.0;
                if (q == nl) {
                    mm[ipos] += mesh->elem_area[elem] / 12.0;
                }
            }
        }
    }

    /* --- b) consistency check (ice_fct.F90:1213-1248) -----------------------
     * sum of row entries = area(ulevels_nod2D(q), q). Print warning if any
     * row mismatches by > 0.1. Fortran prints rank/index then continues. */
    for (int q = 0; q < N_own; ++q) {
        int ul = mesh->ulevels_nod2D[q];   /* 1-based level */
        if (ul > 1) continue;              /* cavity */
        int rstart = stiff->rowptr[q];
        int rend   = stiff->rowptr[q + 1];
        real_t aa = 0.0;
        for (int k = rstart; k < rend; ++k) aa += mm[k];
        real_t a_node = mesh->area[q * mesh->nl + (ul - 1)];   /* area(ul, q) */
        if (fabs(a_node - aa) > 0.1) {
            fprintf(stderr,
                    "#### MASS MATRIX PROBLEM rank=%d q=%d aa=%g area=%g ul=%d\n",
                    partit ? partit->mype : 0, q, (double)aa, (double)a_node, ul);
        }
    }
}

/* ============================================================ */
/* E2 — ice_TG_rhs (ice_fct.F90:91)                             */
/* ============================================================ */

void fesom_ice_tg_rhs(fesom_ice                *ice,
                      struct fesom_partit      *partit,
                      const struct fesom_mesh  *mesh)
{
    (void)partit;

    real_t *u_ice  = ice->uice;
    real_t *v_ice  = ice->vice;
    real_t *a_ice  = ice->data[FESOM_ICE_AICE].values;
    real_t *m_ice  = ice->data[FESOM_ICE_MICE].values;
    real_t *m_snow = ice->data[FESOM_ICE_MSNOW].values;
    real_t *rhs_a  = ice->data[FESOM_ICE_AICE].values_rhs;
    real_t *rhs_m  = ice->data[FESOM_ICE_MICE].values_rhs;
    real_t *rhs_ms = ice->data[FESOM_ICE_MSNOW].values_rhs;

    const real_t dt   = ice->ice_dt;
    const real_t idiff = ice->ice_diff;
    const real_t scale = (real_t)FESOM_ICE_SCALE_AREA;

    /* zero rhs (ice_fct.F90:137-144) — bound: myDim_nod2D only.
     * Halo rhs entries are intentionally NOT zeroed — Fortran does the same;
     * downstream consumers (solve_low/high_order) read rhs only at myDim. */
    for (int row = 0; row < mesh->myDim_nod2D; ++row) {
        rhs_m [row] = 0.0;
        rhs_a [row] = 0.0;
        rhs_ms[row] = 0.0;
    }

    /* assemble (ice_fct.F90:159-198) — bound: myDim_elem2D; cavity skip.
     * Writes to rhs_*[row] where row is an elnode (can be halo); halo
     * accumulations are stale-summed but harmless (consumers read myDim). */
    for (int elem = 0; elem < mesh->myDim_elem2D; ++elem) {
        if (mesh->ulevels[elem] > 1) continue;          /* cavity (line 163) */

        const int en[3] = { mesh->elem_nodes[3*elem + 0],
                            mesh->elem_nodes[3*elem + 1],
                            mesh->elem_nodes[3*elem + 2] };
        const real_t *g  = &mesh->gradient_sca[6 * elem];
        const real_t dx[3] = { g[0], g[1], g[2] };
        const real_t dy[3] = { g[3], g[4], g[5] };
        const real_t vol  = mesh->elem_area[elem];
        /* Fortran: um = sum(U_ice(elnodes)) — NOT divided by 3 (line 171) */
        const real_t um = u_ice[en[0]] + u_ice[en[1]] + u_ice[en[2]];
        const real_t vm = v_ice[en[0]] + v_ice[en[1]] + v_ice[en[2]];
        const real_t diff = idiff * sqrt(vol / scale);

        for (int n = 0; n < 3; ++n) {
            int row = en[n];
            real_t entries[3];
            for (int q = 0; q < 3; ++q) {
                /* Fortran line 184-187 */
                real_t a = (dx[n] * (um + u_ice[en[q]])
                          + dy[n] * (vm + v_ice[en[q]])) / 12.0;
                real_t b = diff * (dx[n] * dx[q] + dy[n] * dy[q]);
                real_t c = 0.5 * dt
                         * (um * dx[n] + vm * dy[n])
                         * (um * dx[q] + vm * dy[q]) / 9.0;
                entries[q] = vol * dt * (a - b - c);
            }
            real_t sm  = entries[0] * m_ice [en[0]]
                       + entries[1] * m_ice [en[1]]
                       + entries[2] * m_ice [en[2]];
            real_t sa  = entries[0] * a_ice [en[0]]
                       + entries[1] * a_ice [en[1]]
                       + entries[2] * a_ice [en[2]];
            real_t sms = entries[0] * m_snow[en[0]]
                       + entries[1] * m_snow[en[1]]
                       + entries[2] * m_snow[en[2]];
            rhs_m [row] += sm;
            rhs_a [row] += sa;
            rhs_ms[row] += sms;
        }
    }
}

/* ============================================================ */
/* E3 — ice_solve_low_order (ice_fct.F90:243)                   */
/* ============================================================ */

static void ice_solve_low_order(fesom_ice                    *ice,
                                const struct fesom_ssh_stiff *stiff,
                                struct fesom_partit          *partit,
                                const struct fesom_mesh      *mesh)
{
    real_t *a_ice  = ice->data[FESOM_ICE_AICE].values;
    real_t *m_ice  = ice->data[FESOM_ICE_MICE].values;
    real_t *m_snow = ice->data[FESOM_ICE_MSNOW].values;
    real_t *rhs_a  = ice->data[FESOM_ICE_AICE].values_rhs;
    real_t *rhs_m  = ice->data[FESOM_ICE_MICE].values_rhs;
    real_t *rhs_ms = ice->data[FESOM_ICE_MSNOW].values_rhs;
    real_t *a_l    = ice->data[FESOM_ICE_AICE].valuesl;
    real_t *m_l    = ice->data[FESOM_ICE_MICE].valuesl;
    real_t *ms_l   = ice->data[FESOM_ICE_MSNOW].valuesl;
    const real_t  *mm = ice->work.fct_massmatrix;

    const real_t gamma = ice->ice_gamma_fct;
    const int    nl    = mesh->nl;

    /* Fortran line 304: do row=1, myDim_nod2D — myDim only. */
    for (int row = 0; row < mesh->myDim_nod2D; ++row) {
        if (mesh->ulevels_nod2D[row] > 1) continue;        /* cavity (line 307) */

        int rs = stiff->rowptr[row];
        int re = stiff->rowptr[row + 1];

        real_t sm = 0.0, sa = 0.0, sms = 0.0;
        for (int k = rs; k < re; ++k) {
            int j = stiff->colind[k];
            sm  += mm[k] * m_ice [j];
            sa  += mm[k] * a_ice [j];
            sms += mm[k] * m_snow[j];
        }
        /* area(1, row) — Fortran is (1, row) i.e. surface level (1-based);
         * our area[] is flat [n*nl + level], 0-based. Surface = level 0. */
        real_t inv_area = 1.0 / mesh->area[row * nl + 0];
        m_l [row] = (rhs_m [row] + gamma * sm ) * inv_area + (1.0 - gamma) * m_ice [row];
        a_l [row] = (rhs_a [row] + gamma * sa ) * inv_area + (1.0 - gamma) * a_ice [row];
        ms_l[row] = (rhs_ms[row] + gamma * sms) * inv_area + (1.0 - gamma) * m_snow[row];
    }

    /* Fortran line 335: exchange_nod(m_icel, a_icel, m_snowl) */
    fesom_exchange_nod2D(m_l,  partit);
    fesom_exchange_nod2D(a_l,  partit);
    fesom_exchange_nod2D(ms_l, partit);
}

/* ============================================================ */
/* E3 — ice_solve_high_order (ice_fct.F90:347)                  */
/* ============================================================ */

static void ice_solve_high_order(fesom_ice                    *ice,
                                 const struct fesom_ssh_stiff *stiff,
                                 struct fesom_partit          *partit,
                                 const struct fesom_mesh      *mesh)
{
    real_t *rhs_a  = ice->data[FESOM_ICE_AICE].values_rhs;
    real_t *rhs_m  = ice->data[FESOM_ICE_MICE].values_rhs;
    real_t *rhs_ms = ice->data[FESOM_ICE_MSNOW].values_rhs;
    real_t *a_l    = ice->data[FESOM_ICE_AICE].valuesl;
    real_t *m_l    = ice->data[FESOM_ICE_MICE].valuesl;
    real_t *ms_l   = ice->data[FESOM_ICE_MSNOW].valuesl;
    real_t *da     = ice->data[FESOM_ICE_AICE].dvalues;
    real_t *dm     = ice->data[FESOM_ICE_MICE].dvalues;
    real_t *dms    = ice->data[FESOM_ICE_MSNOW].dvalues;
    const real_t  *mm = ice->work.fct_massmatrix;

    const int nl = mesh->nl;
    const int num_iter_solve = 3;   /* Fortran line 362 */

    /* First approximation (lines 400-410): myDim_nod2D bound */
    for (int row = 0; row < mesh->myDim_nod2D; ++row) {
        if (mesh->ulevels_nod2D[row] > 1) continue;
        real_t inv_area = 1.0 / mesh->area[row * nl + 0];
        dm [row] = rhs_m [row] * inv_area;
        da [row] = rhs_a [row] * inv_area;
        dms[row] = rhs_ms[row] * inv_area;
    }
    fesom_exchange_nod2D(dm,  partit);
    fesom_exchange_nod2D(da,  partit);
    fesom_exchange_nod2D(dms, partit);

    /* Iterate (lines 426-488): num_iter_solve - 1 = 2 passes */
    for (int it = 0; it < num_iter_solve - 1; ++it) {
        /* (a) update valuesl  (lines 433-452) — myDim_nod2D bound */
        for (int row = 0; row < mesh->myDim_nod2D; ++row) {
            if (mesh->ulevels_nod2D[row] > 1) continue;
            int rs = stiff->rowptr[row];
            int re = stiff->rowptr[row + 1];

            real_t sm = 0.0, sa = 0.0, sms = 0.0;
            for (int k = rs; k < re; ++k) {
                int j = stiff->colind[k];
                sm  += mm[k] * dm [j];
                sa  += mm[k] * da [j];
                sms += mm[k] * dms[j];
            }
            real_t inv_area = 1.0 / mesh->area[row * nl + 0];
            real_t rhs_new = rhs_m[row] - sm;
            m_l [row] = dm[row] + rhs_new * inv_area;
            rhs_new = rhs_a[row] - sa;
            a_l [row] = da[row] + rhs_new * inv_area;
            rhs_new = rhs_ms[row] - sms;
            ms_l[row] = dms[row] + rhs_new * inv_area;
        }

        /* (b) copy valuesl back to dvalues  (lines 464-473) — myDim_nod2D */
        for (int row = 0; row < mesh->myDim_nod2D; ++row) {
            if (mesh->ulevels_nod2D[row] > 1) continue;
            dm [row] = m_l [row];
            da [row] = a_l [row];
            dms[row] = ms_l[row];
        }

        /* (c) exchange (line 481) */
        fesom_exchange_nod2D(dm,  partit);
        fesom_exchange_nod2D(da,  partit);
        fesom_exchange_nod2D(dms, partit);
    }
}

/* ============================================================ */
/* E4 — ice_fem_fct (ice_fct.F90:498) — Zalesak limiter         */
/* ============================================================ */

static void ice_fem_fct(int                           log_id,
                        fesom_ice                    *ice,
                        const struct fesom_ssh_stiff *stiff,
                        struct fesom_partit          *partit,
                        const struct fesom_mesh      *mesh)
{
    /* Resolve the active tracer triple (values, valuesl, dvalues).
     * Fortran uses tr_array_id ∈ {1,2,3} with the mapping
     *   1 → m_ice, 2 → a_ice, 3 → m_snow
     * (NB: this is *not* the data() array index in Fortran either; it's a
     * logical tracer id). We mirror with FCT_LOG_M / FCT_LOG_A / FCT_LOG_MS. */
    const int  d_active  = data_idx_for_logical(log_id);
    real_t    *vals      = ice->data[d_active].values;     /* m_ice / a_ice / m_snow */
    real_t    *vals_l    = ice->data[d_active].valuesl;    /* low-order */
    real_t    *dvals     = ice->data[d_active].dvalues;    /* high-order increment */

    real_t    *icefluxes = ice->work.fct_fluxes;           /* [elem*3] */
    real_t    *icepplus  = ice->work.fct_plus;
    real_t    *icepminus = ice->work.fct_minus;
    real_t    *tmax      = ice->work.fct_tmax;
    real_t    *tmin      = ice->work.fct_tmin;

    const real_t gamma = ice->ice_gamma_fct;
    const int N_full = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl     = mesh->nl;

    /* Init tmax/tmin (lines 562-565) — myDim+eDim. */
    for (int n = 0; n < N_full; ++n) { tmax[n] = 0.0; tmin[n] = 0.0; }

    /* icoef(3,3) Fortran: 1 everywhere except diag = -2 (lines 575-582). */
    static const real_t icoef[3][3] = {
        { -2.0,  1.0,  1.0 },
        {  1.0, -2.0,  1.0 },
        {  1.0,  1.0, -2.0 },
    };

    /* Antidiffusive fluxes per element (lines 594-633) — myDim_elem2D, cavity skip.
     * icefluxes is sized myDim+eDim+eXDim; we only write myDim_elem2D, halo
     * entries left untouched (matches Fortran). */
    for (int elem = 0; elem < mesh->myDim_elem2D; ++elem) {
        if (mesh->ulevels[elem] > 1) continue;
        const int en[3] = { mesh->elem_nodes[3*elem + 0],
                            mesh->elem_nodes[3*elem + 1],
                            mesh->elem_nodes[3*elem + 2] };
        const real_t vol = mesh->elem_area[elem];

        /* For each q (vertex): flux_q = -Σ_n icoef(n,q) * (γ·val + dval)(en[n])
         *                              * vol / area(1, en[q]) / 12 */
        for (int q = 0; q < 3; ++q) {
            real_t s = 0.0;
            for (int n = 0; n < 3; ++n) {
                real_t v = gamma * vals[en[n]] + dvals[en[n]];
                s += icoef[n][q] * v;
            }
            real_t inv_area_q = 1.0 / mesh->area[en[q] * nl + 0];
            icefluxes[elem * 3 + q] = -s * (vol * inv_area_q) / 12.0;
        }
    }

    /* Cluster min/max (lines 654-718) — myDim_nod2D, cavity skip.
     * Only the active tracer is touched (Fortran has 3 separate if-blocks;
     * our log_id selects). */
    for (int row = 0; row < mesh->myDim_nod2D; ++row) {
        if (mesh->ulevels_nod2D[row] > 1) continue;
        int rs = stiff->rowptr[row];
        int re = stiff->rowptr[row + 1];

        real_t lo_l = +INFINITY, hi_l = -INFINITY;
        real_t lo_v = +INFINITY, hi_v = -INFINITY;
        for (int k = rs; k < re; ++k) {
            int j = stiff->colind[k];
            real_t a = vals_l[j], b = vals[j];
            if (a > hi_l) hi_l = a;
            if (a < lo_l) lo_l = a;
            if (b > hi_v) hi_v = b;
            if (b < lo_v) lo_v = b;
        }
        real_t hi = hi_l > hi_v ? hi_l : hi_v;
        real_t lo = lo_l < lo_v ? lo_l : lo_v;
        tmax[row] = hi - vals_l[row];
        tmin[row] = lo - vals_l[row];
    }

    /* Init icepplus/icepminus (lines 751-754) — myDim+eDim */
    for (int n = 0; n < N_full; ++n) { icepplus[n] = 0.0; icepminus[n] = 0.0; }

    /* Sum positive/negative fluxes per node (lines 770-805) — myDim_elem2D,
     * cavity skip, scatter to elnodes (which can be halo). icepplus/icepminus
     * have halo coverage and the exchange below sweeps owner sums. */
    for (int elem = 0; elem < mesh->myDim_elem2D; ++elem) {
        if (mesh->ulevels[elem] > 1) continue;
        const int en[3] = { mesh->elem_nodes[3*elem + 0],
                            mesh->elem_nodes[3*elem + 1],
                            mesh->elem_nodes[3*elem + 2] };
        for (int q = 0; q < 3; ++q) {
            int n = en[q];
            real_t flux = icefluxes[elem * 3 + q];
            if (flux > 0.0) icepplus [n] += flux;
            else            icepminus[n] += flux;
        }
    }

    /* Correction factors (lines 824-841) — myDim_nod2D, cavity skip */
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) continue;
        real_t flux = icepplus[n];
        if (fabs(flux) > 0.0) {
            real_t denom = flux > 1e-12 ? flux : 1e-12;
            real_t v = tmax[n] / denom;
            icepplus[n] = v < 1.0 ? v : 1.0;
        } else {
            icepplus[n] = 0.0;
        }
        flux = icepminus[n];
        if (fabs(flux) > 0.0) {
            real_t denom = flux < -1e-12 ? flux : -1e-12;
            real_t v = tmin[n] / denom;
            icepminus[n] = v < 1.0 ? v : 1.0;
        } else {
            icepminus[n] = 0.0;
        }
    }
    /* Fortran line 853: exchange_nod(icepminus, icepplus) — halo filled. */
    fesom_exchange_nod2D(icepplus,  partit);
    fesom_exchange_nod2D(icepminus, partit);

    /* Limit element fluxes (lines 865-879) — myDim_elem2D, cavity skip.
     * Reads icepplus/icepminus at halo elnodes, which are valid post-exchange. */
    for (int elem = 0; elem < mesh->myDim_elem2D; ++elem) {
        if (mesh->ulevels[elem] > 1) continue;
        const int en[3] = { mesh->elem_nodes[3*elem + 0],
                            mesh->elem_nodes[3*elem + 1],
                            mesh->elem_nodes[3*elem + 2] };
        real_t ae = 1.0;
        for (int q = 0; q < 3; ++q) {
            int n = en[q];
            real_t flux = icefluxes[elem * 3 + q];
            real_t cand = (flux >= 0.0) ? icepplus[n] : icepminus[n];
            if (cand < ae) ae = cand;
        }
        icefluxes[elem * 3 + 0] *= ae;
        icefluxes[elem * 3 + 1] *= ae;
        icefluxes[elem * 3 + 2] *= ae;
    }

    /* Apply update (lines 893-947, plus the analogous blocks for tr_array_id 2,3).
     * Two-phase: first overwrite values <- valuesl on myDim, then accumulate
     * limited fluxes from elements (myDim_elem2D), then exchange owners. */
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) continue;
        vals[n] = vals_l[n];
    }
    for (int elem = 0; elem < mesh->myDim_elem2D; ++elem) {
        if (mesh->ulevels[elem] > 1) continue;
        const int en[3] = { mesh->elem_nodes[3*elem + 0],
                            mesh->elem_nodes[3*elem + 1],
                            mesh->elem_nodes[3*elem + 2] };
        for (int q = 0; q < 3; ++q) {
            vals[en[q]] += icefluxes[elem * 3 + q];
        }
    }
    /* Note: Fortran's outer exchange_nod(m_ice, a_ice, m_snow) at line 1130
     * happens after ALL three fem_fct calls. We emit the per-tracer exchange
     * here too because each tracer's `vals` is independent and the next
     * fem_fct call may read a different tracer's halo (it doesn't, but the
     * post-driver exchange in Fortran covers all three regardless). */
    fesom_exchange_nod2D(vals, partit);
}

/* ============================================================ */
/* E5 — ice_fct_solve driver (ice_fct.F90:210)                  */
/* ============================================================ */

void fesom_ice_fct_solve(fesom_ice                    *ice,
                         const struct fesom_ssh_stiff *stiff,
                         struct fesom_partit          *partit,
                         const struct fesom_mesh      *mesh)
{
    /* Fortran order:
     *   1. ice_solve_high_order  (uses valuesl as scratch — must precede low_order)
     *   2. ice_solve_low_order
     *   3. ice_fem_fct(1)  ! m_ice
     *   4. ice_fem_fct(2)  ! a_ice
     *   5. ice_fem_fct(3)  ! m_snow
     */
    ice_solve_high_order(ice, stiff, partit, mesh);
    ice_solve_low_order (ice, stiff, partit, mesh);

    ice_fem_fct(FCT_LOG_M,  ice, stiff, partit, mesh);
    ice_fem_fct(FCT_LOG_A,  ice, stiff, partit, mesh);
    ice_fem_fct(FCT_LOG_MS, ice, stiff, partit, mesh);
}
