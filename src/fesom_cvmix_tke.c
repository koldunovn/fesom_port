/*
 * cvmix_tke.F90 — classical-TKE column core. Literal port; every formula,
 * loop bound, branch and evaluation order mirrors the Fortran (line refs
 * inline). Dead arguments dropped per the verified list in the header.
 *
 * Bit-fidelity notes:
 *  - the Intel reference build compiles with -r8 (src/CMakeLists.txt:335):
 *    DEFAULT-REAL literals like `6.6` (:717) and `0.5` are DOUBLE — port them
 *    as plain C double literals (verified by controlled replay; the
 *    single-precision-promotion guess was measurably wrong).
 *  - `forc_tke_surf**(3./2.)` (:793,:886): 3./2. is single 1.5 = exactly 1.5;
 *    Intel under -fp-model precise emits a pow call → C pow(x, 1.5).
 *    If the dump-diff ever shows last-bit noise exactly here, the compiler
 *    used x*sqrt(x) — switch TKE_POW32 accordingly.
 *  - solve_tridiag (cvmix_utils_addon.F90:116-149) divides via reciprocal
 *    (fxa = 1D0/m; cp = c*fxa) — ported as-is, NOT as cp = c/m.
 *  - Fortran max/min map to compare-select ternaries (matching maxsd/minsd
 *    operand semantics for finite values), the codebase convention.
 */
#include "fesom_cvmix_tke.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* per-column stack-scratch cap (matches KPP_NL_MAX / fesom_gm NL_MAX). */
enum { TKE_NL_MAX = 128 };

/* the 6.6 literal at :717 — the Intel build compiles with -r8
 * (src/CMakeLists.txt:335), so default-real literals are DOUBLE; plain 6.6.
 * (Empirically pinned by the T2 controlled replay: (double)6.6f produced
 * prandtl diffs of exactly (6.6 − (double)6.6f)·Rinum at 27k entries.) */
#define TKE_C66  6.6

/* x**(3./2.) at :793/:886 */
#define TKE_POW32(x)  pow((x), 1.5)

#define TKE_MIN2(a, b)  ((a) < (b) ? (a) : (b))
#define TKE_MAX2(a, b)  ((a) > (b) ? (a) : (b))

/* tke_constants_saved (:95) — one static instance; the FESOM call site never
 * passes tke_userdef_constants. */
static fesom_cvmix_tke_params tke_constants_saved;

const fesom_cvmix_tke_params *fesom_cvmix_tke_constants(void)
{
    return &tke_constants_saved;
}

/* init_tke (:101-281). The FESOM call site (gen:258-272) passes every
 * parameter below, so each `present()` branch is taken; the range checks are
 * ported verbatim (print + stop 1 → fprintf + abort). handle_old_vals is NOT
 * ported (dead — see header); its default (1) had no executed reader. */
void fesom_cvmix_init_tke(real_t c_k, real_t c_eps, real_t cd,
                          real_t alpha_tke, real_t mxl_min,
                          real_t kappaM_min, real_t kappaM_max,
                          int    tke_mxl_choice,
                          int    use_ubound_dirichlet,
                          int    use_lbound_dirichlet,
                          int    only_tke, int l_lc, real_t clc,
                          real_t tke_min, real_t tke_surf_min)
{
    fesom_cvmix_tke_params *p = &tke_constants_saved;

    if (c_k < 0.0 || c_k > 1.5) {                                /* :135 */
        fprintf(stderr, "ERROR:c_k can only be allowed_range\n"); abort();
    }
    p->c_k = c_k;
    if (c_eps < 0.0 || c_eps > 10.0) {                           /* :145 */
        fprintf(stderr, "ERROR:c_eps can only be allowed_range\n"); abort();
    }
    p->c_eps = c_eps;
    if (cd < 0.1 || cd > 30.0) {                                 /* :155 */
        fprintf(stderr, "ERROR:cd can only be allowed_range\n"); abort();
    }
    p->cd = cd;
    if (alpha_tke < 1.0 || alpha_tke > 90.0) {                   /* :165 */
        fprintf(stderr, "ERROR:alpha_tke can only be allowed_range\n"); abort();
    }
    p->alpha_tke = alpha_tke;
    if (mxl_min < 1.0e-12 || mxl_min > 0.4) {                    /* :175 */
        fprintf(stderr, "ERROR:mxl_min can only be allowed_range\n"); abort();
    }
    p->mxl_min = mxl_min;
    if (kappaM_min < 0.0 || kappaM_min > 1.0) {                  /* :185 */
        fprintf(stderr, "ERROR:KappaM_min can only be allowed_range\n"); abort();
    }
    p->kappaM_min = kappaM_min;
    if (kappaM_max < 10.0 || kappaM_max > 1000.0) {              /* :195 */
        fprintf(stderr, "ERROR:kappaM_max can only be allowed_range\n"); abort();
    }
    p->kappaM_max = kappaM_max;
    if (tke_mxl_choice < 1 || tke_mxl_choice > 2) {              /* :205 */
        fprintf(stderr, "ERROR:tke_mxl_choice can only be 1 or 2\n"); abort();
    }
    p->tke_mxl_choice = tke_mxl_choice;
    if (clc < 0.0 || clc > 30.0) {                               /* :225 */
        fprintf(stderr, "ERROR:clc can only be allowed_range\n"); abort();
    }
    p->clc = clc;
    if (tke_min < 1.0e-9 || tke_min > 1.0e-2) {                  /* :236 */
        fprintf(stderr, "ERROR:tke_min can only be allowed_range\n"); abort();
    }
    p->tke_min = tke_min;
    if (tke_surf_min < 1.0e-7 || tke_surf_min > 1.0e-2) {        /* :246 */
        fprintf(stderr, "ERROR:tke_surf_min can only be allowed_range\n"); abort();
    }
    p->tke_surf_min = tke_surf_min;

    p->use_ubound_dirichlet = use_ubound_dirichlet;              /* :255 */
    p->use_lbound_dirichlet = use_lbound_dirichlet;              /* :261 */
    p->only_tke             = only_tke;                          /* :267 */
    p->l_lc                 = l_lc;                              /* :274 */
}

/* solve_tridiag (cvmix_utils_addon.F90:116-149) — Thomas algorithm with the
 * exact reciprocal-multiply form. n = nlev+1 equations, 0-based here. */
static void tke_solve_tridiag(const real_t *a, const real_t *b,
                              const real_t *c, const real_t *d,
                              real_t *x, int n)
{
    real_t cp[TKE_NL_MAX + 1], dp[TKE_NL_MAX + 1];
    cp[0] = c[0] / b[0];                                         /* :134 */
    dp[0] = d[0] / b[0];                                         /* :135 */
    for (int i = 1; i < n; ++i) {                                /* :137 */
        real_t m   = b[i] - cp[i - 1] * a[i];
        real_t fxa = 1.0 / m;                                    /* :139 */
        cp[i] = c[i] * fxa;
        dp[i] = (d[i] - dp[i - 1] * a[i]) * fxa;
    }
    x[n - 1] = dp[n - 1];                                        /* :144 */
    for (int i = n - 2; i >= 0; --i)                             /* :146 */
        x[i] = dp[i] - cp[i] * x[i + 1];
}

/*===========================================================================
 * integrate_tke (:415-987). 0-based k ↔ Fortran k+1; arrays span
 * k = 0..nlev (interfaces) except dzw (0..nlev-1, layers).
 *===========================================================================*/
void fesom_cvmix_integrate_tke(int           nlev,
                               real_t        dtime,
                               real_t        rho_ref,
                               real_t        grav,
                               const real_t *dzw,
                               const real_t *dzt,
                               const real_t *tke_old,
                               const real_t *Ssqr,
                               const real_t *Nsqr,
                               const real_t *alpha_c,
                               const real_t *E_iw,
                               const real_t *iw_diss,
                               const real_t *tke_plc,
                               real_t        forc_tke_surf,
                               real_t        forc_rho_surf,
                               real_t        bottom_fric,
                               real_t       *tke_new,
                               real_t       *KappaM_out,
                               real_t       *KappaH_out,
                               const fesom_cvmix_tke_diag *diag)
{
    (void)alpha_c;      /* read only under .not.only_tke (:711) — gate-only */
    (void)E_iw;         /* read only under .not.only_tke (:711) — gate-only */
    (void)bottom_fric;  /* never read in the executed body */

    /* local column scratch (:535-581). The 13 diag terms are ALWAYS computed
     * here in Fortran order and copied to the nullable targets at the END —
     * never write through `diag` mid-computation. */
    real_t tke_unrest[TKE_NL_MAX + 1], tke_upd[TKE_NL_MAX + 1];
    real_t mxl[TKE_NL_MAX + 1], sqrttke[TKE_NL_MAX + 1];
    real_t prandtl[TKE_NL_MAX + 1], Rinum[TKE_NL_MAX + 1];
    real_t K_diss_v[TKE_NL_MAX + 1], P_diss_v[TKE_NL_MAX + 1];
    real_t forc[TKE_NL_MAX + 1];
    real_t a_dif[TKE_NL_MAX + 1], b_dif[TKE_NL_MAX + 1], c_dif[TKE_NL_MAX + 1];
    real_t a_tri[TKE_NL_MAX + 1], b_tri[TKE_NL_MAX + 1], c_tri[TKE_NL_MAX + 1];
    real_t d_tri[TKE_NL_MAX + 1], ke[TKE_NL_MAX + 1];
    real_t d_Tbpr[TKE_NL_MAX + 1], d_Tspr[TKE_NL_MAX + 1];
    real_t d_Tdif[TKE_NL_MAX + 1], d_Tdis[TKE_NL_MAX + 1];
    real_t d_Twin[TKE_NL_MAX + 1], d_Tiwf[TKE_NL_MAX + 1];
    real_t d_Tbck[TKE_NL_MAX + 1], d_Ttot[TKE_NL_MAX + 1];
    real_t d_Lmix[TKE_NL_MAX + 1], d_Pr[TKE_NL_MAX + 1];
    real_t d_int1[TKE_NL_MAX + 1], d_int2[TKE_NL_MAX + 1], d_int3[TKE_NL_MAX + 1];
    real_t tke_surf, diff_surf_forc, diff_bott_forc;
    int    k;

    const fesom_cvmix_tke_params *cst = &tke_constants_saved;  /* :590 */
    const real_t alpha_tke    = cst->alpha_tke;                /* :629-642 */
    const real_t c_eps        = cst->c_eps;
    const real_t cd           = cst->cd;
    const real_t KappaM_max   = cst->kappaM_max;
    const real_t mxl_min      = cst->mxl_min;
    const real_t c_k          = cst->c_k;
    const real_t tke_min      = cst->tke_min;
    const int    only_tke     = cst->only_tke;
    const int    l_lc         = cst->l_lc;
    /* tke_surf_min loaded at :636 — used only in the (gated) Dirichlet
     * branch; kappaM_min = 0.0 at :662 with the :663 clamp commented out. */

    if (nlev + 1 > TKE_NL_MAX) {
        fprintf(stderr, "fesom_cvmix_integrate_tke: nlev=%d exceeds TKE_NL_MAX\n",
                nlev);
        abort();
    }

    /* initialize diagnostics + work arrays (:602-623). tke_Twin/Tiwf/Tbck and
     * friends are zeroed here; entries not later overwritten keep 0. */
    for (k = 0; k <= nlev; ++k) {
        d_Tbpr[k] = 0.0; d_Tspr[k] = 0.0; d_Tdif[k] = 0.0; d_Tdis[k] = 0.0;
        d_Twin[k] = 0.0; d_Tiwf[k] = 0.0; d_Tbck[k] = 0.0; d_Ttot[k] = 0.0;
        d_int1[k] = 0.0; d_int2[k] = 0.0; d_int3[k] = 0.0;
        tke_new[k] = 0.0;
        tke_upd[k] = 0.0;
        a_dif[k] = 0.0; b_dif[k] = 0.0; c_dif[k] = 0.0;
        a_tri[k] = 0.0; b_tri[k] = 0.0; c_tri[k] = 0.0;
    }
    tke_surf = 0.0;
    (void)tke_surf;   /* consumed only by the gated Dirichlet branch (:783) */

    /*-----------------------------------------------------------------------
     * Part 1: mixing length scale (:666-698)
     *---------------------------------------------------------------------*/
    for (k = 0; k <= nlev; ++k)                                  /* :668 */
        sqrttke[k] = sqrt(TKE_MAX2(0.0, tke_old[k]));
    for (k = 0; k <= nlev; ++k)                                  /* :671 */
        mxl[k] = sqrt(2.0) * sqrttke[k] / sqrt(TKE_MAX2(1.0e-12, Nsqr[k]));

    /* executed: tke_mxl_choice==2 (:674-685); choice 3 not ported (the
     * reference pins 2; init aborts on anything else, mirroring :695-697). */
    mxl[0]    = 0.0;                                             /* :676 */
    mxl[nlev] = 0.0;                                             /* :677 */
    for (k = 1; k <= nlev - 1; ++k)                              /* :678 */
        mxl[k] = TKE_MIN2(mxl[k], mxl[k - 1] + dzw[k - 1]);
    mxl[nlev - 1] = TKE_MIN2(mxl[nlev - 1], mxl_min + dzw[nlev - 1]); /* :681 */
    for (k = nlev - 2; k >= 1; --k)                              /* :682 */
        mxl[k] = TKE_MIN2(mxl[k], mxl[k + 1] + dzw[k]);
    for (k = 0; k <= nlev; ++k)                                  /* :685 */
        mxl[k] = TKE_MAX2(mxl[k], mxl_min);

    /*-----------------------------------------------------------------------
     * Part 2: diffusivities (:700-720)
     *---------------------------------------------------------------------*/
    for (k = 0; k <= nlev; ++k)                                  /* :704 */
        KappaM_out[k] = TKE_MIN2(KappaM_max, c_k * mxl[k] * sqrttke[k]);
    for (k = 0; k <= nlev; ++k)                                  /* :705 */
        Rinum[k] = Nsqr[k] / TKE_MAX2(Ssqr[k], 1.0e-12);

    /* :710-712 — IDEMIX Ri modification: gate-only (.not.only_tke never
     * executes; alpha_c/E_iw unread). */

    for (k = 0; k <= nlev; ++k)                                  /* :717 */
        prandtl[k] = TKE_MAX2(1.0, TKE_MIN2(10.0, TKE_C66 * Rinum[k]));
    for (k = 0; k <= nlev; ++k)                                  /* :720 */
        KappaH_out[k] = KappaM_out[k] / prandtl[k];

    /*-----------------------------------------------------------------------
     * Part 3: tke forcing (:723-744)
     *---------------------------------------------------------------------*/
    for (k = 0; k <= nlev; ++k)                                  /* :726 */
        forc[k] = 0.0;
    for (k = 0; k <= nlev; ++k)                                  /* :729-730 */
        K_diss_v[k] = Ssqr[k] * KappaM_out[k];
    for (k = 0; k <= nlev; ++k)
        P_diss_v[k] = Nsqr[k] * KappaH_out[k];
    P_diss_v[0] = -forc_rho_surf * grav / rho_ref;               /* :733 */
    for (k = 0; k <= nlev; ++k)                                  /* :734 */
        forc[k] = forc[k] + K_diss_v[k] - P_diss_v[k];

    if (l_lc) {                                                  /* :737-739 */
        for (k = 0; k <= nlev; ++k)
            forc[k] = forc[k] + tke_plc[k];
    }

    /* :742-744 — iw_diss forcing: gate-only (.not.only_tke). */

    /*-----------------------------------------------------------------------
     * Part 4: implicit vertical diffusion + dissipation (:746-827)
     *---------------------------------------------------------------------*/
    for (k = 0; k <= nlev; ++k)                                  /* :749 */
        ke[k] = 0.0;
    for (k = 0; k <= nlev - 1; ++k) {                            /* :750-754 */
        /* Fortran k_F = k+1; kp1 = min(k_F+1, nlev); kk = max(k_F, 2) —
         * 0-based: */
        int kp1 = TKE_MIN2(k + 1, nlev - 1);
        int kk  = TKE_MAX2(k, 1);
        ke[k] = alpha_tke * 0.5 * (KappaM_out[kp1] + KappaM_out[kk]);
    }

    for (k = 0; k <= nlev - 1; ++k)                              /* :757-760 */
        c_dif[k] = ke[k] / (dzt[k] * dzw[k]);
    c_dif[nlev] = 0.0;                                           /* :761 */

    for (k = 1; k <= nlev - 1; ++k)                              /* :764-767 */
        b_dif[k] = ke[k - 1] / (dzt[k] * dzw[k - 1]) + ke[k] / (dzt[k] * dzw[k]);

    for (k = 1; k <= nlev; ++k)                                  /* :770-773 */
        a_dif[k] = ke[k - 1] / (dzt[k] * dzw[k - 1]);
    a_dif[0] = 0.0;                                              /* :774 */

    for (k = 0; k <= nlev; ++k)                                  /* :777 */
        tke_upd[k] = tke_old[k];

    /* upper boundary: the executed NEUMANN branch (:791-796); the Dirichlet
     * body (:780-790) is NOT ported (init aborts if enabled). */
    forc[0] = forc[0] + (cd * TKE_POW32(forc_tke_surf)) / dzt[0];  /* :793 */
    b_dif[0] = ke[0] / (dzt[0] * dzw[0]);                          /* :794 */
    diff_surf_forc = 0.0;                                          /* :795 */

    /* lower boundary: the executed NEUMANN branch (:812-815); Dirichlet body
     * (:799-811) NOT ported. */
    b_dif[nlev] = ke[nlev - 1] / (dzt[nlev] * dzw[nlev - 1]);      /* :813 */
    diff_bott_forc = 0.0;                                          /* :814 */

    for (k = 0; k <= nlev; ++k) {                                /* :818-821 */
        a_tri[k] = -dtime * a_dif[k];
        b_tri[k] = 1.0 + dtime * b_dif[k];
        c_tri[k] = -dtime * c_dif[k];
    }
    for (k = 1; k <= nlev - 1; ++k)                              /* :820 */
        b_tri[k] = b_tri[k] + dtime * c_eps * sqrttke[k] / mxl[k];

    for (k = 0; k <= nlev; ++k)                                  /* :824 */
        d_tri[k] = tke_upd[k] + dtime * forc[k];

    tke_solve_tridiag(a_tri, b_tri, c_tri, d_tri, tke_new, nlev + 1); /* :827 */

    /* implicit-tendency diagnostics (:829-849) */
    for (k = 1; k <= nlev - 1; ++k)                              /* :831-833 */
        d_Tdif[k] = a_dif[k] * tke_new[k - 1] - b_dif[k] * tke_new[k]
                  + c_dif[k] * tke_new[k + 1];
    d_Tdif[0]    = -b_dif[0] * tke_new[0] + c_dif[0] * tke_new[1];    /* :834 */
    d_Tdif[nlev] = a_dif[nlev] * tke_new[nlev - 1] - b_dif[nlev] * tke_new[nlev]; /* :835 */
    d_Tdif[1]        = d_Tdif[1] + diff_surf_forc;                    /* :836 */
    d_Tdif[nlev - 1] = d_Tdif[nlev - 1] + diff_bott_forc;             /* :837 */
    /* :841-849 — Dirichlet Tdif overrides: gate-only. */

    /* dissipation of TKE (:852-853) */
    for (k = 0; k <= nlev; ++k)
        d_Tdis[k] = 0.0;
    for (k = 1; k <= nlev - 1; ++k)
        d_Tdis[k] = -c_eps / mxl[k] * sqrttke[k] * tke_new[k];

    /*-----------------------------------------------------------------------
     * Part 5: reset tke to bounding values (:856-869)
     *---------------------------------------------------------------------*/
    for (k = 0; k <= nlev; ++k)                                  /* :860 */
        tke_unrest[k] = tke_new[k];

    if (only_tke) {                                              /* :867-869 */
        for (k = 0; k <= nlev; ++k)
            tke_new[k] = TKE_MAX2(tke_new[k], tke_min);
    }

    /*-----------------------------------------------------------------------
     * Part 6: assign diagnostic variables (:871-904)
     *---------------------------------------------------------------------*/
    for (k = 0; k <= nlev; ++k)                                  /* :876-877 */
        d_Tbpr[k] = -P_diss_v[k];
    for (k = 0; k <= nlev; ++k)
        d_Tspr[k] = K_diss_v[k];
    for (k = 0; k <= nlev; ++k)                                  /* :880 */
        d_Tbck[k] = (tke_new[k] - tke_unrest[k]) / dtime;
    /* surface: executed Neumann else-branch (:885-886); Dirichlet (:882-883)
     * gate-only. */
    d_Twin[0] = (cd * TKE_POW32(forc_tke_surf)) / dzt[0];        /* :886 */
    /* bottom: executed Neumann else-branch (:893-895). */
    d_Twin[nlev] = 0.0;                                          /* :895 */

    for (k = 0; k <= nlev; ++k)                                  /* :898 */
        d_Tiwf[k] = iw_diss[k];   /* unconditional read — zero col under only_tke */
    for (k = 0; k <= nlev; ++k)                                  /* :899 */
        d_Ttot[k] = (tke_new[k] - tke_old[k]) / dtime;
    /* :901-904 — tke_Lmix(nlev+1:)=0 / tke_Pr(nlev+1:)=0 touch only element
     * nlev of the slice, then the full-range assignment overwrites it: */
    for (k = 0; k <= nlev; ++k)
        d_Lmix[k] = mxl[k];
    for (k = 0; k <= nlev; ++k)
        d_Pr[k] = prandtl[k];

    /* debugging copies (:909-911) */
    for (k = 0; k <= nlev; ++k) {
        d_int1[k] = KappaH_out[k];
        d_int2[k] = KappaM_out[k];
        d_int3[k] = Nsqr[k];
    }
    /* :916-986 — if(.false.) debug print block: dead. */

    /* end-of-column diag copy-out (the ONLY writes through `diag`) */
    if (diag) {
        size_t nb = (size_t)(nlev + 1) * sizeof(real_t);
        if (diag->Tbpr) memcpy(diag->Tbpr, d_Tbpr, nb);
        if (diag->Tspr) memcpy(diag->Tspr, d_Tspr, nb);
        if (diag->Tdif) memcpy(diag->Tdif, d_Tdif, nb);
        if (diag->Tdis) memcpy(diag->Tdis, d_Tdis, nb);
        if (diag->Twin) memcpy(diag->Twin, d_Twin, nb);
        if (diag->Tiwf) memcpy(diag->Tiwf, d_Tiwf, nb);
        if (diag->Tbck) memcpy(diag->Tbck, d_Tbck, nb);
        if (diag->Ttot) memcpy(diag->Ttot, d_Ttot, nb);
        if (diag->Lmix) memcpy(diag->Lmix, d_Lmix, nb);
        if (diag->Pr)   memcpy(diag->Pr,   d_Pr,   nb);
        if (diag->int1) memcpy(diag->int1, d_int1, nb);
        if (diag->int2) memcpy(diag->int2, d_int2, nb);
        if (diag->int3) memcpy(diag->int3, d_int3, nb);
    }
}
