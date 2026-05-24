/*
 * KPP vertical mixing — faithful port of oce_ale_mixing_kpp.F90.
 * See fesom_kpp.h for the struct layout and docs/plans/20260524-kpp-vertical-mixing.md
 * for the phased port plan. Phase K0 lands the scaffolding (struct alloc/free,
 * abort-stub driver, dump harness); the routine bodies arrive K1..K8.
 */
#include "fesom_kpp.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================
 * KPP constants — verified against oce_ale_mixing_kpp.F90 + the CORE2
 * reference namelists (docs/kpp_reference_namelists/). Module parameters and
 * namelist values that the CORE2/KPP run uses. A_ver/K_ver reuse the shared
 * FESOM_PHASE1_{A,K}_VER (== the CORE2 1e-4 / 1e-5; same constant as PP).
 *===========================================================================*/
/* module parameters (oce_ale_mixing_kpp.F90:50-59) */
#define KPP_EPSLN        1.0e-40    /* epsln  (:50)  */
#define KPP_EPSILON      0.1        /* epsilon_kpp (:52) */
#define KPP_VONK         0.4        /* vonk   (:53)  */
#define KPP_CONC1        5.0        /* conc1  (:54)  */
#define KPP_ZMIN        (-4.0e-7)   /* zmin   (:56)  lookup-table zehat min */
#define KPP_ZMAX         0.0        /* zmax   (:57)  */
#define KPP_UMIN         0.0        /* umin   (:58)  lookup-table ustar min */
#define KPP_UMAX         0.04       /* umax   (:59)  */
/* namelist values (work_core; module defaults coincide except A_ver) */
#define KPP_RICR         0.3        /* Ricr   (namelist.oce:70) */
#define KPP_CONCV        1.6        /* concv  (namelist.oce:71) */
#define KPP_VISC_SH_LIMIT 5.0e-3    /* visc_sh_limit (namelist.oce:66) */
#define KPP_DIFF_SH_LIMIT 5.0e-3    /* diff_sh_limit (namelist.tra:97) */
/* double diffusion (ddmix) — CORE2 has double_diffusion=.false. (namelist.tra:105),
 * so the IF(double_diffusion) CALL ddmix gate (driver :420-422) is never taken.
 * Port-what-CORE2-uses: the ddmix body (:1012-1085) is DEFERRED. If a config ever
 * sets this true, port ddmix first — the #error below makes that explicit. */
#define KPP_DOUBLE_DIFFUSION 0
/* in-routine constants used in later phases (K3/K5/K7) */
#define KPP_RIINFTY      0.8        /* ri_iwmix shear-Ri shape limit (:737) */
#define KPP_MINMIX       3.0e-3     /* surface viscAE floor (driver :408)   */
#define KPP_CEKMAN       0.7        /* bldepth Ekman limit (:478)           */
#define KPP_CMONOB       1.0        /* bldepth Monin-Obukhov limit (:479)   */

/* wm/ws lookup-table index: Fortran wmt(0:nni+1, 0:nnj+1), i in [0,nni+1],
 * j in [0,nnj+1]. Row-major with i as the row (nnj+2 columns). */
#define KPP_TBL(i, j)  ((size_t)(i) * (size_t)(FESOM_KPP_NNJ + 2) + (size_t)(j))

/* fesom_kpp_mixing call index (= ocean step); mirrors the Fortran module's
 * kpp_call_count. Used for dump-step gating. */
static int s_kpp_call = 0;

/*===========================================================================
 * Allocation — mirrors oce_mixing_kpp_init's allocate block
 * (oce_ale_mixing_kpp.F90:128-156), sizing every array exactly as Fortran
 * (myDim+eDim nodes; ghats on nl-1, the rest on nl). All zeroed (calloc),
 * matching the Fortran init loop that sets every element to 0.
 *===========================================================================*/
void fesom_kpp_alloc(fesom_kpp *k, const struct fesom_mesh *mesh, int num_tracers)
{
    memset(k, 0, sizeof(*k));
    int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int nl = mesh->nl;
    k->n_nod       = N;
    k->nl          = nl;
    k->num_tracers = num_tracers;

    size_t n_nl   = (size_t)N * (size_t)nl;
    size_t n_nl1  = (size_t)N * (size_t)(nl - 1);
    size_t tbl    = (size_t)(FESOM_KPP_NNI + 2) * (size_t)(FESOM_KPP_NNJ + 2);

    k->diffK  = calloc((size_t)num_tracers * n_nl, sizeof(real_t));
    k->viscA  = calloc(n_nl,                       sizeof(real_t));
    k->blmc   = calloc((size_t)3 * n_nl,           sizeof(real_t));
    k->ghats  = calloc(n_nl1,                      sizeof(real_t));
    k->dVsq   = calloc(n_nl,                       sizeof(real_t));

    k->dkm1   = calloc((size_t)3 * (size_t)N,      sizeof(real_t));
    k->hbl    = calloc((size_t)N,                  sizeof(real_t));
    k->bfsfc  = calloc((size_t)N,                  sizeof(real_t));
    k->caseA  = calloc((size_t)N,                  sizeof(real_t));
    k->stable = calloc((size_t)N,                  sizeof(real_t));
    k->ustar  = calloc((size_t)N,                  sizeof(real_t));
    k->Bo     = calloc((size_t)N,                  sizeof(real_t));
    k->kbl    = calloc((size_t)N,                  sizeof(int));

    k->wmt    = calloc(tbl,                        sizeof(real_t));
    k->wst    = calloc(tbl,                        sizeof(real_t));

    FESOM_CHECK(k->diffK && k->viscA && k->blmc && k->ghats && k->dVsq
             && k->dkm1 && k->hbl && k->bfsfc && k->caseA && k->stable
             && k->ustar && k->Bo && k->kbl && k->wmt && k->wst,
             "fesom_kpp alloc: out of memory");
}

void fesom_kpp_free(fesom_kpp *k)
{
    free(k->diffK);  free(k->viscA);  free(k->blmc);  free(k->ghats);
    free(k->dVsq);   free(k->dkm1);   free(k->hbl);   free(k->bfsfc);
    free(k->caseA);  free(k->stable); free(k->ustar); free(k->Bo);
    free(k->kbl);    free(k->wmt);    free(k->wst);
    memset(k, 0, sizeof(*k));
}

/*===========================================================================
 * One-time init — Vtc, cg, deltaz/deltau, wm/ws lookup tables.
 * Literal port of oce_mixing_kpp_init (:111-220), constant block + table build.
 * (The Fortran allocate+zero is handled by fesom_kpp_alloc's calloc.)
 *===========================================================================*/
void fesom_kpp_init(fesom_kpp *k)
{
    /* table-build parameters (oce_ale_mixing_kpp.F90:115-123) */
    const real_t cstar = 10.0,    conam = 1.257, concm = 8.380;
    const real_t conc2 = 16.0,    zetam = -0.2;
    const real_t conas = -28.86,  concs = 98.96, conc3 = 16.0, zetas = -1.0;

    /* Vtc (eqn 23, :177): concv * sqrt(0.2/concs/epsilon_kpp) / vonk^2 / Ricr */
    k->Vtc = KPP_CONCV * sqrt(0.2 / concs / KPP_EPSILON)
             / (KPP_VONK * KPP_VONK) / KPP_RICR;

    /* cg (eqn 20, :188): cstar * vonk * (concs * vonk * epsilon_kpp)^(1/3) */
    k->cg = cstar * KPP_VONK * pow(concs * KPP_VONK * KPP_EPSILON, 1.0 / 3.0);

    /* table steps (:194-195): real(nni+1)/real(nnj+1) divisor */
    k->deltaz = (KPP_ZMAX - KPP_ZMIN) / (real_t)(FESOM_KPP_NNI + 1);
    k->deltau = (KPP_UMAX - KPP_UMIN) / (real_t)(FESOM_KPP_NNJ + 1);

    /* wm/ws tables (eqn 13 & B1, :197-219) */
    for (int i = 0; i <= FESOM_KPP_NNI + 1; ++i) {
        real_t zehat = k->deltaz * (real_t)i + KPP_ZMIN;
        for (int j = 0; j <= FESOM_KPP_NNJ + 1; ++j) {
            real_t usta = k->deltau * (real_t)j + KPP_UMIN;
            real_t zeta = zehat / (usta * usta * usta + KPP_EPSLN);
            real_t wm, ws;
            if (zehat >= 0.0) {
                wm = KPP_VONK * usta / (1.0 + KPP_CONC1 * zeta);
                ws = wm;
            } else {
                if (zeta > zetam)
                    wm = KPP_VONK * usta * pow(1.0 - conc2 * zeta, 1.0 / 4.0);
                else
                    wm = KPP_VONK * pow(conam * usta*usta*usta - concm * zehat,
                                        1.0 / 3.0);
                if (zeta > zetas)
                    ws = KPP_VONK * usta * pow(1.0 - conc3 * zeta, 1.0 / 2.0);
                else
                    ws = KPP_VONK * pow(conas * usta*usta*usta - concs * zehat,
                                        1.0 / 3.0);
            }
            k->wmt[KPP_TBL(i, j)] = wm;
            k->wst[KPP_TBL(i, j)] = ws;
        }
    }
}

/*===========================================================================
 * wscale — turbulent velocity scales wm, ws from the 2D lookup table.
 * Literal port of wscale (oce_ale_mixing_kpp.F90:828-877). Lookup for
 * zehat <= zmax (unstable); stable formula otherwise. Private (static).
 *===========================================================================*/
static void kpp_wscale(const fesom_kpp *k, real_t zehat, real_t us,
                       real_t *wm, real_t *ws)
{
    if (zehat <= KPP_ZMAX) {
        real_t zdiff = zehat - KPP_ZMIN;
        int iz = (int)(zdiff / k->deltaz);          /* INT (trunc) */
        if (iz > FESOM_KPP_NNI) iz = FESOM_KPP_NNI;
        if (iz < 0)             iz = 0;
        int izp1 = iz + 1;

        real_t udiff = us - KPP_UMIN;
        real_t uq = udiff / k->deltau;               /* MIN(.,nnj) then INT */
        if (uq > (real_t)FESOM_KPP_NNJ) uq = (real_t)FESOM_KPP_NNJ;
        int ju = (int)uq;
        if (ju < 0) ju = 0;
        int jup1 = ju + 1;

        real_t zfrac = zdiff / k->deltaz - (real_t)iz;
        real_t ufrac = udiff / k->deltau - (real_t)ju;
        real_t fzfrac = 1.0 - zfrac;

        real_t wam = fzfrac * k->wmt[KPP_TBL(iz, jup1)]
                   + zfrac  * k->wmt[KPP_TBL(izp1, jup1)];
        real_t wbm = fzfrac * k->wmt[KPP_TBL(iz, ju)]
                   + zfrac  * k->wmt[KPP_TBL(izp1, ju)];
        *wm = (1.0 - ufrac) * wbm + ufrac * wam;

        real_t was = fzfrac * k->wst[KPP_TBL(iz, jup1)]
                   + zfrac  * k->wst[KPP_TBL(izp1, jup1)];
        real_t wbs = fzfrac * k->wst[KPP_TBL(iz, ju)]
                   + zfrac  * k->wst[KPP_TBL(izp1, ju)];
        *ws = (1.0 - ufrac) * wbs + ufrac * was;
    } else {
        real_t u3 = us * us * us;
        *wm = KPP_VONK * us * u3 / (u3 + KPP_CONC1 * zehat + KPP_EPSLN);
        *ws = *wm;
    }
}

/*===========================================================================
 * ri_iwmix — interior viscosity/diffusivity: shear instability (local Ri) +
 * constant background (A_ver/K_ver) + static instability. Literal port of
 * ri_iwmix (oce_ale_mixing_kpp.F90:890-999). Loops myDim only (the driver
 * exchanges viscA/diffK afterward). Ri is stored in diffK ch0 as scratch.
 * smooth_Ri_ver / smooth_Ri_hor are .false. (CORE2) → omitted (no-op gates).
 *===========================================================================*/
static void kpp_ri_iwmix(fesom_kpp *k, const struct fesom_aux *aux,
                         const struct fesom_dyn *dyn, const struct fesom_mesh *mesh)
{
    const int nl   = k->nl;
    const int Nmy  = mesh->myDim_nod2D;
    const size_t slab = (size_t)k->n_nod * (size_t)nl;
    real_t *viscA  = k->viscA;
    real_t *diffKt = k->diffK + 0 * slab;   /* channel 0 = T (diffK(:,:,1)) */
    real_t *diffKs = k->diffK + 1 * slab;   /* channel 1 = S (diffK(:,:,2)) */
    const real_t A_bg = (real_t)FESOM_PHASE1_A_VER;
    const real_t K_bg = (real_t)FESOM_PHASE1_K_VER;

    /* Loop 1: Richardson number Ri = MAX(N²,0)/(shear²+epsln), stored in
     * diffKt as scratch (Fortran :917-946). shear at Z (mid-layer). */
    for (int n = 0; n < Nmy; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t dz_inv = 1.0 / (mesh->Z[nz - 1] - mesh->Z[nz]);   /* > 0 */
            real_t du = dyn->uvnode[FESOM_ELEMVEC(n, nz - 1, nl) + 0]
                      - dyn->uvnode[FESOM_ELEMVEC(n, nz,     nl) + 0];
            real_t dv = dyn->uvnode[FESOM_ELEMVEC(n, nz - 1, nl) + 1]
                      - dyn->uvnode[FESOM_ELEMVEC(n, nz,     nl) + 1];
            real_t shear = (du * du + dv * dv) * dz_inv * dz_inv;
            real_t Nsq = aux->bvfreq[FESOM_NODE3D(n, nz, nl)];
            real_t Nsq_pos = Nsq > 0.0 ? Nsq : 0.0;
            diffKt[FESOM_NODE3D(n, nz, nl)] = Nsq_pos / (shear + KPP_EPSLN);
        }
        /* edge copies of the Ri scratch (:932-933) */
        diffKt[FESOM_NODE3D(n, nzmin, nl)] = diffKt[FESOM_NODE3D(n, nzmin + 1, nl)];
        diffKt[FESOM_NODE3D(n, nzmax, nl)] = diffKt[FESOM_NODE3D(n, nzmax - 1, nl)];
    }

    /* Loop 2: viscA / diffK from the shear-Ri shape function (:954-997). */
    for (int n = 0; n < Nmy; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            size_t i = FESOM_NODE3D(n, nz, nl);
            real_t Rigg = diffKt[i] > 0.0 ? diffKt[i] : 0.0;        /* AMAX1(Ri,0) */
            real_t ratio = Rigg / KPP_RIINFTY;
            if (ratio > 1.0) ratio = 1.0;                          /* AMIN1(.,1) */
            real_t frit = 1.0 - ratio * ratio;
            frit = frit * frit * frit;
            viscA[i]  = KPP_VISC_SH_LIMIT * frit + A_bg;
            real_t dK = KPP_DIFF_SH_LIMIT * frit + K_bg;           /* Kv0_const branch */
            diffKt[i] = dK;
            diffKs[i] = dK;                                        /* diffK(2)=diffK(1) */
        }
        /* edge copies (:989-995) */
        size_t a0 = FESOM_NODE3D(n, nzmin, nl),     a1 = FESOM_NODE3D(n, nzmin + 1, nl);
        size_t b0 = FESOM_NODE3D(n, nzmax, nl),     b1 = FESOM_NODE3D(n, nzmax - 1, nl);
        viscA[a0]  = viscA[a1];  diffKt[a0] = diffKt[a1]; diffKs[a0] = diffKs[a1];
        viscA[b0]  = viscA[b1];  diffKt[b0] = diffKt[b1]; diffKs[b0] = diffKs[b1];
    }
}

/* generic per-node column getter for fesom_kpp_dump_nodes: arr is [N*nl]
 * node-major; comp == layer interface nz. */
typedef struct { const real_t *arr; int nl; } kpp_col_src;
static double kpp_get_col(int node, int comp, void *user)
{
    const kpp_col_src *s = (const kpp_col_src *)user;
    return (double)s->arr[(size_t)node * (size_t)s->nl + (size_t)comp];
}

/* K3 dump: raw ri_iwmix node outputs (viscA, diffK T/S) as full columns, plus
 * bvfreq (the input whose sign drives the step-1 background/static-instability
 * branch) so the validation can confirm any C-vs-Fortran viscA disagreement is
 * exactly a bvfreq sign flip (cross-code IC noise), not an algebra bug. */
static void kpp_dump_riiwmix(const fesom_kpp *k, const struct fesom_aux *aux,
                             const struct fesom_mesh *mesh, struct fesom_partit *partit)
{
    size_t slab = (size_t)k->n_nod * (size_t)k->nl;
    kpp_col_src sv = { k->viscA,           k->nl };
    kpp_col_src st = { k->diffK + 0 * slab, k->nl };
    kpp_col_src ss = { k->diffK + 1 * slab, k->nl };
    kpp_col_src sb = { aux->bvfreq,         k->nl };
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "ri_viscA",  k->nl, kpp_get_col, &sv);
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "ri_diffKt", k->nl, kpp_get_col, &st);
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "ri_diffKs", k->nl, kpp_get_col, &ss);
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "ri_bvfreq", k->nl, kpp_get_col, &sb);
}

/*===========================================================================
 * bldepth — oceanic boundary-layer depth hbl + kbl, bfsfc/stable/caseA.
 * Literal port of bldepth (oce_ale_mixing_kpp.F90:674-854). HIGHEST RISK
 * (historical OBL bug, FRESH_START:489). Loops myDim only (Fortran :697/708/813
 * are myDim with the +eDim commented out). Uses zbar (interface depths, NOT Z),
 * the static kpp_wscale (K2), and the driver pre-step ustar/Bo/dVsq + the
 * EOS-filled dbsfc/bvfreq/sw_alpha + forcing sw_3d.
 *
 * CORE2 use_sw_pene = .true. (Key invariant #6) — the sw-penetration branches
 * are always taken, so they are ported inline (port-what-CORE2-uses); the
 * else-branch (bfsfc=Bo constant) is not reachable in CORE2. smooth_hbl=.false.
 * (skip). kbl is stored 0-based (the dump emits kbl+1 to match the Fortran's
 * 1-based real(kbl,WP)).
 *===========================================================================*/
static void kpp_bldepth(fesom_kpp *k, const struct fesom_aux *aux,
                        const struct fesom_forcing *forcing,
                        const struct fesom_mesh *mesh)
{
    const int nl  = k->nl;
    const int Nmy = mesh->myDim_nod2D;
    const real_t Vtc = k->Vtc;

    /* init hbl/kbl to bottomed-out values (:696-702) */
    for (int n = 0; n < Nmy; ++n) {
        int nzmax = mesh->nlevels_nod2D[n] - 1;          /* 0-based deepest interface */
        k->kbl[n] = nzmax;
        k->hbl[n] = fabs(mesh->zbar_3d_n[FESOM_NODE3D(n, nzmax, nl)]);
    }

    /* bulk-Richardson search for hbl (:708-802) */
    for (int n = 0; n < Nmy; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        real_t Bo = k->Bo[n];
        real_t us = k->ustar[n];
        real_t coeff_sw = (real_t)FESOM_G * aux->sw_alpha[FESOM_NODE3D(n, nzmin, nl)]; /* :713 */
        real_t sw_surf  = forcing->sw_3d[FESOM_NODE3D(n, nzmin, nl)];
        real_t Rib_km1  = 0.0;
        k->bfsfc[n] = Bo;                                /* :717 */

        for (int nz = nzmin + 1; nz <= nzmax; ++nz) {    /* :719 (nz=nzmin+1..nzmax) */
            real_t zk   = fabs(mesh->zbar_3d_n[FESOM_NODE3D(n, nz,     nl)]);
            real_t zkm1 = fabs(mesh->zbar_3d_n[FESOM_NODE3D(n, nz - 1, nl)]);

            /* bfsfc = Bo + sw contribution (:725-727) */
            k->bfsfc[n] = Bo + coeff_sw * (sw_surf - forcing->sw_3d[FESOM_NODE3D(n, nz, nl)]);

            real_t bfs    = k->bfsfc[n];
            real_t stable = 0.5 + copysign(0.5, bfs);                  /* :729 */
            k->stable[n]  = stable;
            real_t sigma  = stable + (1.0 - stable) * KPP_EPSILON;     /* :730 */

            real_t zehat = KPP_VONK * sigma * zk * bfs;                /* :736 */
            real_t wm, ws;
            kpp_wscale(k, zehat, us, &wm, &ws);                        /* :737 */

            real_t bvsq = aux->bvfreq[FESOM_NODE3D(n, nz, nl)];        /* :747 */
            real_t Vtsq = zk * ws * sqrt(fabs(bvsq)) * Vtc;            /* :750 */

            real_t Ritop = zk * aux->dbsfc[FESOM_NODE3D(n, nz, nl)];   /* :757 */
            real_t Rib_k = Ritop / (k->dVsq[FESOM_NODE3D(n, nz, nl)] + Vtsq + KPP_EPSLN); /* :758 */
            real_t dzup  = zk - zkm1;                                  /* :759 */

            if (Rib_k > KPP_RICR) {                                    /* :761 */
                k->hbl[n] = zkm1 + dzup * (KPP_RICR - Rib_km1)
                                       / (Rib_k - Rib_km1 + KPP_EPSLN);/* :763 */
                k->kbl[n] = nz;                                        /* :764 */
                break;                                                 /* :765 (EXIT) */
            } else {
                Rib_km1 = Rib_k;                                       /* :767 */
            }

            /* linear interp of sw_3d to hbl, refine bfsfc (:774-785) */
            {
                real_t sw_km1 = forcing->sw_3d[FESOM_NODE3D(n, nz - 1, nl)];
                real_t sw_k   = forcing->sw_3d[FESOM_NODE3D(n, nz,     nl)];
                k->bfsfc[n] = Bo + coeff_sw *
                    ( sw_surf - ( sw_km1 + (sw_k - sw_km1) * (k->hbl[n] - zkm1) / dzup ) );
                real_t st = 0.5 + copysign(0.5, k->bfsfc[n]);          /* :783 */
                k->stable[n] = st;
                k->bfsfc[n]  = k->bfsfc[n] + st * KPP_EPSLN;           /* :784 */
            }
        }

        /* hekman / hmonob limits, gated bfsfc>0 .and. nzmin==1 (:792-801).
         * Fortran nzmin==1 (no cavity) → C nzmin==0. */
        if (k->bfsfc[n] > 0.0 && nzmin == 0) {
            real_t hekman = KPP_CEKMAN * us / fmax(fabs(mesh->coriolis_node[n]), KPP_EPSLN); /* :795 */
            real_t hmonob = KPP_CMONOB * us * us * us
                            / KPP_VONK / (k->bfsfc[n] + KPP_EPSLN);    /* :796-797 */
            real_t hlimit = k->stable[n] * fmin(hekman, hmonob);       /* :798 */
            k->hbl[n] = fmin(k->hbl[n], hlimit);                       /* :799 */
            k->hbl[n] = fmax(k->hbl[n],
                             fabs(mesh->zbar_3d_n[FESOM_NODE3D(n, 1, nl)])); /* :800 (zbar(2)) */
        }
    }

    /* smooth_hbl = .false. → skip (:806-809) */

    /* find new kbl from final hbl, refine bfsfc, set caseA (:812-852) */
    for (int n = 0; n < Nmy; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;

        k->kbl[n] = nzmax;                                            /* :820 */
        for (int nz = nzmin + 1; nz <= nzmax; ++nz) {                 /* :822 */
            if (fabs(mesh->zbar_3d_n[FESOM_NODE3D(n, nz, nl)]) > k->hbl[n]) {
                k->kbl[n] = nz;                                       /* :824 */
                break;                                               /* :825 (EXIT) */
            }
        }

        /* final bfsfc: linear interp of sw_3d to hbl using kbl (:832-844) */
        int kbl = k->kbl[n];
        real_t coeff_sw = (real_t)FESOM_G * aux->sw_alpha[FESOM_NODE3D(n, nzmin, nl)];
        real_t sw_surf  = forcing->sw_3d[FESOM_NODE3D(n, nzmin,  nl)];
        real_t sw_km1   = forcing->sw_3d[FESOM_NODE3D(n, kbl - 1, nl)];
        real_t sw_k     = forcing->sw_3d[FESOM_NODE3D(n, kbl,     nl)];
        real_t zbar_km1 = mesh->zbar_3d_n[FESOM_NODE3D(n, kbl - 1, nl)];
        real_t zbar_k   = mesh->zbar_3d_n[FESOM_NODE3D(n, kbl,     nl)];
        k->bfsfc[n] = k->Bo[n] + coeff_sw *
            ( sw_surf - ( sw_km1 + (sw_k - sw_km1)
                                   * (k->hbl[n] + zbar_km1) / (zbar_km1 - zbar_k) ) );
        real_t st = 0.5 + copysign(0.5, k->bfsfc[n]);                 /* :842 */
        k->stable[n] = st;
        k->bfsfc[n]  = k->bfsfc[n] + st * KPP_EPSLN;                  /* :843 */

        /* caseA: =1 if hbl above mid-point of level kbl, else 0 (:850-851) */
        real_t dzup = zbar_km1 - zbar_k;
        real_t arg  = fabs(zbar_k) - 0.5 * dzup - k->hbl[n];
        k->caseA[n] = 0.5 + copysign(0.5, arg);
    }
}

/* K5 dump: driver pre-step inputs to bldepth — prestep(ustar,Bo) + dVsq + dbsfc
 * columns. Mirrors the Fortran kpp_dump_prestep (oce_ale_mixing_kpp.F90:554-578). */
static double kpp_get_prestep(int node, int comp, void *user)
{
    const fesom_kpp *k = *(const fesom_kpp *const *)user;
    return comp == 0 ? (double)k->ustar[node] : (double)k->Bo[node];
}
static void kpp_dump_prestep(const fesom_kpp *k, const struct fesom_aux *aux,
                             const struct fesom_mesh *mesh, struct fesom_partit *partit)
{
    const fesom_kpp *kp = k;
    kpp_col_src dv = { k->dVsq,    k->nl };
    kpp_col_src db = { aux->dbsfc, k->nl };
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "prestep", 2, kpp_get_prestep, &kp);
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "dVsq",  k->nl, kpp_get_col, &dv);
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "dbsfc", k->nl, kpp_get_col, &db);
}

/* K5 dump: bldepth outputs — hbl, kbl(+1 → Fortran 1-based), bfsfc, stable, caseA.
 * Mirrors the Fortran kpp_dump_final's 'bldepth' tag (:596-604). */
static double kpp_get_bldepth(int node, int comp, void *user)
{
    const fesom_kpp *k = *(const fesom_kpp *const *)user;
    switch (comp) {
        case 0:  return (double)k->hbl[node];
        case 1:  return (double)(k->kbl[node] + 1);   /* 0-based → Fortran 1-based */
        case 2:  return (double)k->bfsfc[node];
        case 3:  return (double)k->stable[node];
        default: return (double)k->caseA[node];
    }
}
static void kpp_dump_bldepth(const fesom_kpp *k, const struct fesom_mesh *mesh,
                             struct fesom_partit *partit)
{
    const fesom_kpp *kp = k;
    fesom_kpp_dump_nodes(mesh, partit, s_kpp_call, "bldepth", 5, kpp_get_bldepth, &kp);
}

/*===========================================================================
 * KPP driver — mirror of oce_mixing_KPP (oce_ale_mixing_kpp.F90:249-494),
 * built incrementally. Through K5: zero viscA + pre-step (dVsq/ustar/Bo) +
 * ri_iwmix + bldepth (+ dumps). Downstream (blmix_kpp/enhance/combine/
 * single-Kv) land K6..K8 — guarded so only a FESOM_KPP_DUMP_DIR dump run
 * (1 step) may proceed until the driver is complete.
 *===========================================================================*/
void fesom_kpp_mixing(fesom_kpp                  *k,
                      struct fesom_aux           *aux,
                      const struct fesom_tracers *tracers,
                      const struct fesom_forcing *forcing,
                      const struct fesom_dyn     *dyn,
                      const struct fesom_mesh    *mesh,
                      struct fesom_partit        *partit)
{
    ++s_kpp_call;                            /* = ocean step (mirrors Fortran kpp_call_count) */
    const int nl  = k->nl;
    const int Nmy = mesh->myDim_nod2D;

    /* zero viscA over myDim+eDim (driver :340-343) */
    memset(k->viscA, 0, (size_t)k->n_nod * (size_t)nl * sizeof(real_t));

    /* ---- driver pre-step (oce_ale_mixing_kpp.F90:344-411) -------------------
     * dVsq (eqn 21: shear of UVnode re the surface), ustar (eqn 2), Bo (A2c/A2d).
     * dbsfc is filled by the EOS (fesom_eos.c); here only its surface level is
     * zeroed, matching the Fortran driver (:367). */
    const real_t *S = tracers->data[FESOM_TRACER_S].values;
    const real_t density_0_r = 1.0 / (real_t)FESOM_DENSITY_0;

    for (int n = 0; n < Nmy; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        k->dVsq[FESOM_NODE3D(n, nzmin, nl)]    = 0.0;     /* :366 */
        aux->dbsfc[FESOM_NODE3D(n, nzmin, nl)] = 0.0;     /* :367 */
        real_t usurf = dyn->uvnode[FESOM_ELEMVEC(n, nzmin, nl) + 0];   /* :370 */
        real_t vsurf = dyn->uvnode[FESOM_ELEMVEC(n, nzmin, nl) + 1];   /* :371 */
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {      /* :372 (nz=nzmin+1..nzmax-1) */
            real_t u_loc = 0.5 * ( dyn->uvnode[FESOM_ELEMVEC(n, nz - 1, nl) + 0]
                                 + dyn->uvnode[FESOM_ELEMVEC(n, nz,     nl) + 0] );
            real_t v_loc = 0.5 * ( dyn->uvnode[FESOM_ELEMVEC(n, nz - 1, nl) + 1]
                                 + dyn->uvnode[FESOM_ELEMVEC(n, nz,     nl) + 1] );
            real_t du = usurf - u_loc, dv = vsurf - v_loc;
            k->dVsq[FESOM_NODE3D(n, nz, nl)] = du * du + dv * dv;       /* :377 */
        }
        k->dVsq[FESOM_NODE3D(n, nzmax, nl)] =
            k->dVsq[FESOM_NODE3D(n, nzmax - 1, nl)];                    /* :380 */
    }

    for (int n = 0; n < Nmy; ++n) {                       /* :401-411 */
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        real_t sx = forcing->stress_node_surf[2*n + 0];
        real_t sy = forcing->stress_node_surf[2*n + 1];
        k->ustar[n] = sqrt( sqrt(sx*sx + sy*sy) * density_0_r );        /* eqn 2 (:404) */
        k->Bo[n] = -(real_t)FESOM_G
                   * ( aux->sw_alpha[FESOM_NODE3D(n, nzmin, nl)]
                         * forcing->heat_flux[n] / (real_t)FESOM_VCPW
                     + aux->sw_beta[FESOM_NODE3D(n, nzmin, nl)]
                         * forcing->water_flux[n] * S[FESOM_NODE3D(n, nzmin, nl)] ); /* :409-410 */
    }

    /* interior mixing (ri_iwmix, :419) */
    kpp_ri_iwmix(k, aux, dyn, mesh);

    /* K4: double diffusion gate (driver :423-425). CORE2 double_diffusion=.false.
     * → ddmix not called (body deferred, port-what-CORE2-uses). */
#if KPP_DOUBLE_DIFFUSION
#error "KPP double_diffusion=.true. unsupported: port ddmix (oce_ale_mixing_kpp.F90:1012-1085)"
    /* CALL ddmix(diffK, ...) */
#endif

    /* OBL depth (bldepth, :428) — K5 */
    kpp_bldepth(k, aux, forcing, mesh);

    /* dual-instrumentation dumps (step 1): pre-step + ri_iwmix + bldepth */
    if (fesom_kpp_dump_this_step(s_kpp_call)) {
        kpp_dump_prestep(k, aux, mesh, partit);
        kpp_dump_riiwmix(k, aux, mesh, partit);
        kpp_dump_bldepth(k, mesh, partit);
    }

    /* Driver incomplete through K5: blmix_kpp/enhance/combine/single-Kv land
     * K6..K8 (aux->Kv/Av untouched). Only a FESOM_KPP_DUMP_DIR dump run (1 step)
     * may proceed; the partial result is NOT valid for dynamics. */
    if (!fesom_kpp_dump_dir())
        FESOM_DIE("fesom_kpp_mixing: KPP driver incomplete (ported through K5 bldepth; "
                  "blmix_kpp/enhance/combine/single-Kv pending K6..K8). Only a "
                  "FESOM_KPP_DUMP_DIR dump run may proceed. Use FESOM_MIX_SCHEME=PP otherwise.");
}

/*===========================================================================
 * Dump harness — FESOM_KPP_DUMP_DIR-gated, mirrors the EVP dump diagnostic
 * (reference_evp_dump_diagnostic). Off by default (one NULL check / run).
 *===========================================================================*/
static const char *s_kpp_dump_dir    = NULL;
static int         s_kpp_dump_loaded  = 0;
static int         s_kpp_dump_step    = 1;

static void kpp_dump_load_env(void)
{
    if (s_kpp_dump_loaded) return;
    s_kpp_dump_dir = getenv("FESOM_KPP_DUMP_DIR");
    const char *e = getenv("FESOM_KPP_DUMP_STEP");
    if (e) { int v = atoi(e); if (v > 0) s_kpp_dump_step = v; }
    s_kpp_dump_loaded = 1;
}

const char *fesom_kpp_dump_dir(void)
{
    kpp_dump_load_env();
    return s_kpp_dump_dir;
}

int fesom_kpp_dump_this_step(int step_n)
{
    kpp_dump_load_env();
    return s_kpp_dump_dir && (step_n == s_kpp_dump_step);
}

void fesom_kpp_dump_init(const fesom_kpp *k, int rank)
{
    kpp_dump_load_env();
    if (!s_kpp_dump_dir || rank != 0) return;
    char path[1024];
    snprintf(path, sizeof(path), "%s/kpp_init_rank0.txt", s_kpp_dump_dir);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "kpp_dump_init: cannot open %s\n", path); return; }
    fprintf(f, "# Vtc %.17g\n# cg %.17g\n# deltaz %.17g\n# deltau %.17g\n",
            k->Vtc, k->cg, k->deltaz, k->deltau);
    fprintf(f, "# i j wmt wst  (nni=%d nnj=%d)\n", FESOM_KPP_NNI, FESOM_KPP_NNJ);
    for (int i = 0; i <= FESOM_KPP_NNI + 1; ++i)
        for (int j = 0; j <= FESOM_KPP_NNJ + 1; ++j)
            fprintf(f, "%d %d %.17g %.17g\n", i, j,
                    k->wmt[KPP_TBL(i, j)], k->wst[KPP_TBL(i, j)]);
    fclose(f);
}

/* K2 wscale sweep dump: rank 0 evaluates wscale over a fixed (zehat, ustar)
 * grid spanning unstable-table / stable-branch / clamp regions; writes
 * <dir>/kpp_wscale_rank0.txt (`i j zehat ustar wm ws`). The Fortran gate uses
 * the identical grid; both compared vs scripts/kpp_wscale_reference.py. */
#define KPP_SWEEP_NZ 201
#define KPP_SWEEP_NU 101
void fesom_kpp_dump_wscale_sweep(const fesom_kpp *k, int rank)
{
    kpp_dump_load_env();
    if (!s_kpp_dump_dir || rank != 0) return;
    char path[1024];
    snprintf(path, sizeof(path), "%s/kpp_wscale_rank0.txt", s_kpp_dump_dir);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "kpp_dump_wscale: cannot open %s\n", path); return; }
    fprintf(f, "# i j zehat ustar wm ws  (sweep %dx%d)\n",
            KPP_SWEEP_NZ, KPP_SWEEP_NU);
    for (int i = 0; i < KPP_SWEEP_NZ; ++i) {
        real_t zehat = -1.0e-6 + (real_t)i * (2.0e-6 / (real_t)(KPP_SWEEP_NZ - 1));
        for (int j = 0; j < KPP_SWEEP_NU; ++j) {
            real_t ustar = (real_t)j * (0.05 / (real_t)(KPP_SWEEP_NU - 1));
            real_t wm, ws;
            kpp_wscale(k, zehat, ustar, &wm, &ws);
            fprintf(f, "%d %d %.17g %.17g %.17g %.17g\n", i, j,
                    zehat, ustar, wm, ws);
        }
    }
    fclose(f);
}

void fesom_kpp_dump_nodes(const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          int step_n, const char *tag, int ncomp,
                          double (*get)(int node, int comp, void *user),
                          void *user)
{
    kpp_dump_load_env();
    if (!s_kpp_dump_dir) return;
    char path[1024];
    snprintf(path, sizeof(path), "%s/kpp_dump_s%d_%s_rank%d.txt",
             s_kpp_dump_dir, step_n, tag, partit->mype);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "kpp_dump: cannot open %s\n", path); return; }
    fprintf(f, "# step=%d tag=%s rank=%d N=%d ncomp=%d\n",
            step_n, tag, partit->mype, mesh->myDim_nod2D, ncomp);
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        fprintf(f, "%d", partit->myList_nod2D[n]);   /* 1-based gid */
        for (int c = 0; c < ncomp; ++c) fprintf(f, " %.17g", get(n, c, user));
        fputc('\n', f);
    }
    fclose(f);
}
