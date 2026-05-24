#include "fesom_eos.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_mesh.h"
#include "fesom_tracers.h"
#include "fesom_halo.h"

#include <math.h>
#include <stdlib.h>

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
 *   - MLD1 (Large et al. 1997 / FESOM 1.4) is computed (Phase G2a);
 *     MLD2 / MLD3 (Levitus, Griffies) still deferred.
 *   - no N² smoothing
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

    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;     /* 1-based → 0-based */
        int nzmax = mesh->nlevels_nod2D[n] - 1;     /* exclusive bound; layers 0..nzmax-1 */

        /* Initialise MLD1_ind (G2a). Fortran sets MLD1_ind=nzmin+1 (1-based),
         * which in C 0-based is just nzmin+1. Stored 0-based throughout —
         * G3 readers add +1 to mirror Fortran's `bvfreq(MLD1_ind+1, n)`. */
        if (aux->MLD1_ind) aux->MLD1_ind[n] = nzmin + 1;

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

        /* Pass 2: in-situ density at mid-layer depth Z[nz], dbsfc1, db_max.
           Phase 1 has no partial cells / no cavity → Z_3d_n[nz][n] = Z[nz].
           db_max accumulates the max buoyancy-gradient over the full column,
           used by MLD1 (Large et al. 1997). Fortran lines 314-335. */
        real_t db_max = 0.0;
        real_t z_nzmin = mesh->Z[nzmin];
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t z = mesh->Z[nz];
            real_t bulk = bulk_0[nz] + z*(bulk_pz[nz] + z*bulk_pz2[nz]);
            real_t r = bulk * rhopot[nz] / (bulk + 0.1 * z * state_eq_int) - rho_ref;
            rho[nz] = r;
            aux->density_m_rho0[FESOM_NODE3D(n, nz, nl)] = r;

            /* Surface-density adiabatically brought to depth z (Fortran 326-328) */
            real_t bulk_surf = bulk_0[nzmin] + z*(bulk_pz[nzmin] + z*bulk_pz2[nzmin]);
            real_t rho_surf = bulk_surf * rhopot[nzmin]
                            / (bulk_surf + 0.1 * z * state_eq_int);

            /* dbsfc1 (Fortran line 332): use_density_ref=false in our config
             * so density_ref(nz, n) = density_0 = rho_ref. So
             * rho(nz)+density_ref = (r) + rho_ref = rho_in_situ_full. */
            real_t r_full = r + rho_ref;
            real_t dbsfc1 = -g * (rho_surf - r_full) / r_full;

            /* Store dbsfc for KPP bldepth (Fortran oce_ale_pressure_bv.F90:332,339).
             * Written unconditionally: PP never reads aux->dbsfc, so the PP result is
             * unaffected (Fortran's `if (mixing_kpp)` gate is a write-skip with no
             * numerical consequence). Bottom interface filled after the loop (:337). */
            aux->dbsfc[FESOM_NODE3D(n, nz, nl)] = dbsfc1;

            /* db_max accumulator (Fortran line 334). At nz==nzmin the
             * divisor is the |Z[nzmin] - Z[nzmin+1]| placeholder — a
             * non-zero divisor; at that level dbsfc1 is 0 anyway since
             * rho_surf == r_full, so the term contributes nothing. */
            int nz_eff = (nz > nzmin) ? nz : (nzmin + 1);
            real_t denom = fabs(z_nzmin - mesh->Z[nz_eff]);
            real_t cand = dbsfc1 / denom;
            if (cand > db_max) db_max = cand;
        }
        /* dbsfc bottom fill: dbsfc(nzmax)=dbsfc(nzmax-1) (Fortran :337). nzmax is the
         * deepest interface (0-based nlevels-1); the loop filled nzmin..nzmax-1. */
        aux->dbsfc[FESOM_NODE3D(n, nzmax, nl)] =
            aux->dbsfc[FESOM_NODE3D(n, nzmax - 1, nl)];

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
           *same depth* zmean to cancel compressibility. Also locates MLD1
           (Phase G2a): the shallowest level where N² > db_max. */
        int mld1_done = 0;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t zmean   = 0.5 * (mesh->Z[nz - 1] + mesh->Z[nz]);
            real_t bulk_up = bulk_0[nz - 1] + zmean*(bulk_pz[nz - 1] + zmean*bulk_pz2[nz - 1]);
            real_t bulk_dn = bulk_0[nz    ] + zmean*(bulk_pz[nz    ] + zmean*bulk_pz2[nz    ]);
            real_t rho_up  = bulk_up * rhopot[nz - 1] / (bulk_up + 0.1 * zmean * state_eq_int);
            real_t rho_dn  = bulk_dn * rhopot[nz    ] / (bulk_dn + 0.1 * zmean * state_eq_int);
            real_t dz_inv  = 1.0 / (mesh->Z[nz - 1] - mesh->Z[nz]);
            real_t bv = -g * dz_inv * (rho_up - rho_dn) / rho_ref;
            aux->bvfreq[FESOM_NODE3D(n, nz, nl)] = bv;

            /* MLD1: Large et al. (1997). Fortran lines 447-451.
             * Stored 0-based (G3 readers add +1, mirroring Fortran's
             * `bvfreq(MLD1_ind+1, n)`). */
            if (!mld1_done && bv > db_max && aux->MLD1_ind) {
                aux->MLD1_ind[n] = nz;
                mld1_done = 1;
            }
        }
        /* Pad surface and bottom — Fortran lines 474-475:
         *   bvfreq(nzmin) = bvfreq(nzmin+1)        [surface]
         *   bvfreq(nzmax) = bvfreq(nzmax-1)        [bottom interface]
         * nzmax is the deepest interface (0-based nlevels-1); the loop filled the
         * interior nzmin+1..nzmax-1. The horizontal N² smoothing (N2smth_h=.true.,
         * Fortran pressure_bv:499) is applied by fesom_smooth_nod3D from fesom_step
         * after the bvfreq halo exchange. */
        if (nzmin + 1 < nzmax) {
            aux->bvfreq[FESOM_NODE3D(n, nzmin, nl)] =
                aux->bvfreq[FESOM_NODE3D(n, nzmin + 1, nl)];
            aux->bvfreq[FESOM_NODE3D(n, nzmax, nl)] =
                aux->bvfreq[FESOM_NODE3D(n, nzmax - 1, nl)];
        }
    }
}

/*--- smooth_nod3D — area-weighted node-patch horizontal smoother --------------
 * Literal port of smooth_nod3D (gen_support.F90:99-198). Each owned node's value
 * at level nz becomes the area-weighted mean of the three vertices of every
 * surrounding element that reaches level nz:
 *     arr(nz,n) = Σ_el area_el·(arr(nz,V1)+arr(nz,V2)+arr(nz,V3)) / (3·Σ_el area_el)
 * Elements that bottom out above nz drop from both sums (so the patch shrinks
 * with depth, matching Fortran's `nle = min(nln, nlevels(el))`). One full sweep
 * per smoothing cycle, each followed by a halo exchange. The inverse patch area
 * `vol` is built on the first sweep and reused by later sweeps, exactly as the
 * Fortran. The caller must have a valid halo `arr` on entry (the sweep reads halo
 * vertices); the bvfreq use provides that via the preceding fesom_exchange_nod3D
 * (Fortran fills bvfreq on myDim+eDim instead, then smooths). Used for N² with
 * N2smth_h=.true., N2smth_hidx=1 (oce_modules.F90:105). */
void fesom_smooth_nod3D(real_t *arr, int nl, int n_smooth,
                        const struct fesom_mesh *mesh, struct fesom_partit *p)
{
    if (n_smooth < 1) return;
    const int Nmy = mesh->myDim_nod2D;
    const size_t sz = (size_t)Nmy * (size_t)nl;
    real_t *vol  = malloc(sz * sizeof(real_t));
    real_t *work = malloc(sz * sizeof(real_t));
    FESOM_CHECK(vol && work, "fesom_smooth_nod3D: out of memory");

    for (int sweep = 0; sweep < n_smooth; ++sweep) {
        for (int n = 0; n < Nmy; ++n) {
            int uln  = mesh->ulevels_nod2D[n] - 1;     /* 0-based first level   */
            int nlnz = mesh->nlevels_nod2D[n] - 1;     /* 0-based deepest level */
            for (int nz = uln; nz <= nlnz; ++nz) {
                if (sweep == 0) vol[(size_t)n*nl + nz] = 0.0;
                work[(size_t)n*nl + nz] = 0.0;
            }
            int o0 = mesh->nod_in_elem2D_offsets[n];
            int o1 = mesh->nod_in_elem2D_offsets[n + 1];
            for (int k = o0; k < o1; ++k) {
                int el  = mesh->nod_in_elem2D[k];
                int ule = mesh->ulevels[el] - 1;  if (ule < uln)  ule = uln;
                int nle = mesh->nlevels[el] - 1;  if (nle > nlnz) nle = nlnz;
                real_t a = mesh->elem_area[el];
                int v0 = mesh->elem_nodes[3*el + 0];
                int v1 = mesh->elem_nodes[3*el + 1];
                int v2 = mesh->elem_nodes[3*el + 2];
                for (int nz = ule; nz <= nle; ++nz) {
                    if (sweep == 0) vol[(size_t)n*nl + nz] += a;
                    work[(size_t)n*nl + nz] += a * ( arr[(size_t)v0*nl + nz]
                                                   + arr[(size_t)v1*nl + nz]
                                                   + arr[(size_t)v2*nl + nz] );
                }
            }
            /* invert the patch area once (1/(3·Σarea)); reused on later sweeps */
            if (sweep == 0)
                for (int nz = uln; nz <= nlnz; ++nz)
                    vol[(size_t)n*nl + nz] = 1.0 / (3.0 * vol[(size_t)n*nl + nz]);
        }
        /* combined scale + copy back to owned nodes (halo overwritten by exchange) */
        for (int n = 0; n < Nmy; ++n) {
            int uln  = mesh->ulevels_nod2D[n] - 1;
            int nlnz = mesh->nlevels_nod2D[n] - 1;
            for (int nz = uln; nz <= nlnz; ++nz)
                arr[(size_t)n*nl + nz] = work[(size_t)n*nl + nz] * vol[(size_t)n*nl + nz];
        }
        fesom_exchange_nod3D(arr, nl, p);
    }
    free(vol);
    free(work);
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

    /* Interior elements only — gradient_sca and elem_nodes are sized
     * myDim_elem2D in the MPI port. The PGF on halo elements arrives via
     * the elem_area exchange path inside fesom_step's per-step exchanges. */
    for (int e = 0; e < mesh->myDim_elem2D; ++e) {
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

/*--- sw_alpha_beta (oce_ale_pressure_bv.F90:2751-2846) -------------------------
 * McDougall (1987) thermal expansion and saline contraction coefficients per
 * node per level. Inputs: potential T (°C), S (PSU); pressure proxy = |Z[nz]|.
 *
 * For our linfs / no-cavity / no-partial-cells config Z_3d_n[nz, n] = Z[nz].
 *
 * Halo: Fortran exchange_nod on sw_alpha and sw_beta at end. We mirror.
 * ----------------------------------------------------------------------------- */
void fesom_compute_sw_alpha_beta(const struct fesom_tracers *tracers,
                                 const struct fesom_mesh    *mesh,
                                 struct fesom_aux           *aux)
{
    const int nl     = mesh->nl;
    const int myDim  = mesh->myDim_nod2D;
    const real_t *T  = tracers->data[FESOM_TRACER_T].values;
    const real_t *S  = tracers->data[FESOM_TRACER_S].values;

    for (int n = 0; n < myDim; ++n) {
        int nzmax = mesh->nlevels_nod2D[n] - 1;     /* layer count exclusive */
        int nzmin = mesh->ulevels_nod2D[n] - 1;     /* 0-based */
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t t1 = T[FESOM_NODE3D(n, nz, nl)] * 1.00024;
            real_t s1 = S[FESOM_NODE3D(n, nz, nl)];
            real_t p1 = fabs(mesh->Z[nz]);

            real_t t1_2 = t1*t1;
            real_t t1_3 = t1_2*t1;
            real_t t1_4 = t1_3*t1;
            real_t p1_2 = p1*p1;
            real_t p1_3 = p1_2*p1;
            real_t s35  = s1 - 35.0;
            real_t s35_2 = s35*s35;

            real_t beta = 0.785567e-3
                        - 0.301985e-5*t1
                        + 0.555579e-7*t1_2
                        - 0.415613e-9*t1_3
                        + s35*(-0.356603e-6 + 0.788212e-8*t1
                              + 0.408195e-10*p1 - 0.602281e-15*p1_2)
                        + s35_2*(0.515032e-8)
                        + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2)
                        + p1_2*(0.176621e-12 - 0.175379e-14*t1)
                        + p1_3*(0.121551e-17);

            real_t a_over_b = 0.665157e-1
                            + 0.170907e-1*t1
                            - 0.203814e-3*t1_2
                            + 0.298357e-5*t1_3
                            - 0.255019e-7*t1_4
                            + s35*(0.378110e-2 - 0.846960e-4*t1
                                  - 0.164759e-6*p1 - 0.251520e-11*p1_2)
                            + s35_2*(-0.678662e-5)
                            + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2)
                            + p1_2*t1_2*(0.512857e-12)
                            - p1_3*(0.302285e-13);

            aux->sw_beta [FESOM_NODE3D(n, nz, nl)] = beta;
            aux->sw_alpha[FESOM_NODE3D(n, nz, nl)] = a_over_b * beta;
        }
    }
}
