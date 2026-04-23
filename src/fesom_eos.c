#include "fesom_eos.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_mesh.h"
#include "fesom_tracers.h"

#include <math.h>

/*--- JM-EOS components -----------------------------------------------------
 * Literal port of densityJM_components (oce_ale_pressure_bv.F90:2605-2669).
 * All constants reproduced verbatim — do NOT round, fold, or "simplify."
 */
void fesom_eos_jm_components(real_t t, real_t s,
                             real_t *bulk_0, real_t *bulk_pz, real_t *bulk_pz2,
                             real_t *rhopot)
{
    static const real_t a0    = 19092.56,    at   = 209.8925;
    static const real_t at2   = -3.041638,   at3  = -1.852732e-3;
    static const real_t at4   = -1.361629e-5;
    static const real_t as    = 104.4077,    ast  = -6.500517;
    static const real_t ast2  = 0.1553190,   ast3 = 2.326469e-4;
    static const real_t ass   = -5.587545,   asst = 0.7390729;
    static const real_t asst2 = -1.909078e-2;
    static const real_t ap    = -4.721788e-1, apt  = -1.028859e-2;
    static const real_t apt2  = 2.512549e-4,  apt3 = 5.939910e-7;
    static const real_t aps   = 1.571896e-2,  apst = 2.598241e-4;
    static const real_t apst2 = -7.267926e-6, apss = -2.042967e-3;
    static const real_t ap2   = 1.045941e-5,  ap2t = -5.782165e-10;
    static const real_t ap2t2 = 1.296821e-7;
    static const real_t ap2s  = -2.595994e-7, ap2st = -1.248266e-9;
    static const real_t ap2st2= -3.508914e-9;

    static const real_t b0  = 999.842594,   bt   = 6.793952e-2;
    static const real_t bt2 = -9.095290e-3, bt3  = 1.001685e-4;
    static const real_t bt4 = -1.120083e-6, bt5  = 6.536332e-9;
    static const real_t bs  = 0.824493,     bst  = -4.08990e-3;
    static const real_t bst2 = 7.64380e-5,  bst3 = -8.24670e-7;
    static const real_t bst4 = 5.38750e-9;
    static const real_t bss  = -5.72466e-3, bsst = 1.02270e-4;
    static const real_t bsst2 = -1.65460e-6, bss2 = 4.8314e-4;

    real_t s_sqrt = sqrt(s);

    *bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))
             + s* (as  + t*(ast  + t*(ast2 + t*ast3))
                  + s_sqrt*(ass  + t*(asst + t*asst2)));

    *bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))
                    + s*(aps + t*(apst + t*apst2) + s_sqrt*apss);

    *bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)
                   + s *(ap2s + t*(ap2st + t*ap2st2));

    *rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))
                  + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4))))
                  + s*s_sqrt*(bss + t*(bsst + t*bsst2))
                  + s*s* bss2;
    /* Note: Fortran wrote s*(... + s_sqrt*(...) + s*bss2) which is identical
       to s*(...) + s*s_sqrt*(...) + s*s*bss2 — expanded here for clarity. */
}

/*--- pressure_bv (Phase 1 subset) -------------------------------------------
 * Mirror of pressure_bv (oce_ale_pressure_bv.F90:194-501) restricted to:
 *   - state_equation = 1 (JM-EOS)
 *   - which_ALE = 'linfs' (uses the hpressure linfs branch)
 *   - no cavity (nzmin = 0 in C / 1 in Fortran, ulevels_nod2D = 1)
 *   - use_density_ref = .false.  → density_ref ≡ density_0
 *   - mix_scheme ≠ KPP (skip dbsfc)
 *   - no MLD computation, no N² smoothing (deferred)
 *
 * Per-node temporaries are stack-allocated VLAs (one column at a time).
 * For Phase 1 nl ≤ 48 — fine on the stack.
 */
void fesom_pressure_bv(const struct fesom_tracers *tracers,
                       const struct fesom_mesh    *mesh,
                       struct fesom_aux           *aux)
{
    const int nl = mesh->nl;
    const real_t g       = (real_t)FESOM_G;
    const real_t rho_ref = (real_t)FESOM_DENSITY_0;
    /* state_equation_int = 1 in Phase 1 (JM-EOS).
       Multiplied through `0.1 * z * state_eq_int` exactly as Fortran does. */
    const real_t state_eq_int = 1.0;

    const real_t *T = tracers->data[FESOM_TRACER_T].values;
    const real_t *S = tracers->data[FESOM_TRACER_S].values;

    for (int n = 0; n < mesh->nod2D; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;     /* 1-based → 0-based */
        int nzmax = mesh->nlevels_nod2D[n] - 1;     /* exclusive bound; layers 0..nzmax-1 */

        /* Skip nodes with no wet layers (shouldn't happen in Phase 1 but guard). */
        if (nzmax <= nzmin) continue;

        real_t bulk_0[64], bulk_pz[64], bulk_pz2[64], rhopot[64], rho[64];
        /* nl≤48 in Phase 1; static cap of 64 leaves headroom and keeps the
           compiler from generating dynamic-stack code. */

        /* Pass 1: JM-EOS components per layer. */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            size_t i = FESOM_NODE3D(n, nz, nl);
            fesom_eos_jm_components(T[i], S[i],
                                    &bulk_0[nz], &bulk_pz[nz], &bulk_pz2[nz],
                                    &rhopot[nz]);
        }

        /* Pass 2: in-situ density at mid-layer depth Z[nz].
           Phase 1 has no partial cells / no cavity → Z_3d_n[nz][n] = Z[nz]. */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t z = mesh->Z[nz];
            real_t bulk = bulk_0[nz] + z*(bulk_pz[nz] + z*bulk_pz2[nz]);
            real_t r = bulk * rhopot[nz] / (bulk + 0.1 * z * state_eq_int) - rho_ref;
            rho[nz] = r;
            aux->density_m_rho0[FESOM_NODE3D(n, nz, nl)] = r;
        }

        /* hpressure — linfs branch, no cavity (nzmin == 0 in C).
           Mirror of oce_ale_pressure_bv.F90 lines 369-403. Surface boundary:
             hpressure[nzmin] = -Z[nzmin] * rho[nzmin] * g
           Then accumulate downward by 0.5 g * (rho_above*h_above + rho*h). */
        aux->hpressure[FESOM_NODE3D(n, nzmin, nl)] =
            -mesh->Z[nzmin] * rho[nzmin] * g;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t h_up   = mesh->hnode[FESOM_NODE3D(n, nz - 1, nl)];
            real_t h_this = mesh->hnode[FESOM_NODE3D(n, nz,     nl)];
            real_t a = 0.5 * g * (rho[nz - 1] * h_up + rho[nz] * h_this);
            aux->hpressure[FESOM_NODE3D(n, nz, nl)] =
                aux->hpressure[FESOM_NODE3D(n, nz - 1, nl)] + a;
        }

        /* bvfreq — N² between layers nzmin+1 and nzmax-1, then padded.
           Mirror of lines 427-475. Both rho_up and rho_dn are evaluated at the
           *same depth* zmean to cancel compressibility. */
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t zmean   = 0.5 * (mesh->Z[nz - 1] + mesh->Z[nz]);
            real_t bulk_up = bulk_0[nz - 1] + zmean*(bulk_pz[nz - 1] + zmean*bulk_pz2[nz - 1]);
            real_t bulk_dn = bulk_0[nz    ] + zmean*(bulk_pz[nz    ] + zmean*bulk_pz2[nz    ]);
            real_t rho_up  = bulk_up * rhopot[nz - 1] / (bulk_up + 0.1 * zmean * state_eq_int);
            real_t rho_dn  = bulk_dn * rhopot[nz    ] / (bulk_dn + 0.1 * zmean * state_eq_int);
            real_t dz_inv  = 1.0 / (mesh->Z[nz - 1] - mesh->Z[nz]);
            aux->bvfreq[FESOM_NODE3D(n, nz, nl)] =
                -g * dz_inv * (rho_up - rho_dn) / rho_ref;
        }
        /* Pad surface and bottom — Fortran lines 474-475 */
        if (nzmin + 1 < nzmax) {
            aux->bvfreq[FESOM_NODE3D(n, nzmin, nl)] =
                aux->bvfreq[FESOM_NODE3D(n, nzmin + 1, nl)];
            aux->bvfreq[FESOM_NODE3D(n, nzmax - 1, nl)] =
                aux->bvfreq[FESOM_NODE3D(n, nzmax - 2, nl)];
        }
    }
}

/*--- pressure_force_4_linfs_fullcell ---------------------------------------
 * Mirror of oce_ale_pressure_bv.F90:575-614. Five-line element loop:
 *   pgf_x[e][nz] = Σ_i gradient_sca[1..3] * hpressure[nz][V_i(e)] / density_0
 *   pgf_y[e][nz] = Σ_i gradient_sca[4..6] * hpressure[nz][V_i(e)] / density_0
 * Vertex layout in our gradient_sca: [6e + i] for i = 0..5.
 */
void fesom_pressure_force_linfs_fullcell(const struct fesom_mesh *mesh,
                                         struct fesom_aux        *aux)
{
    const int nl       = mesh->nl;
    const real_t inv_r = 1.0 / (real_t)FESOM_DENSITY_0;

    for (int e = 0; e < mesh->elem2D; ++e) {
        int nzmin = mesh->ulevels[e]   - 1;     /* 1-based → 0-based */
        int nzmax = mesh->nlevels[e]   - 1;     /* exclusive bound  */
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];

        const real_t *g = &mesh->gradient_sca[6*e];     /* [dN1x..dN3x, dN1y..dN3y] */

        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t hp0 = aux->hpressure[FESOM_NODE3D(n0, nz, nl)];
            real_t hp1 = aux->hpressure[FESOM_NODE3D(n1, nz, nl)];
            real_t hp2 = aux->hpressure[FESOM_NODE3D(n2, nz, nl)];
            aux->pgf_x[FESOM_ELEM3D(e, nz, nl)] =
                (g[0]*hp0 + g[1]*hp1 + g[2]*hp2) * inv_r;
            aux->pgf_y[FESOM_ELEM3D(e, nz, nl)] =
                (g[3]*hp0 + g[4]*hp1 + g[5]*hp2) * inv_r;
        }
    }
}
