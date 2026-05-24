/*
 * KPP vertical mixing — faithful port of oce_ale_mixing_kpp.F90.
 * See fesom_kpp.h for the struct layout and docs/plans/20260524-kpp-vertical-mixing.md
 * for the phased port plan. Phase K0 lands the scaffolding (struct alloc/free,
 * abort-stub driver, dump harness); the routine bodies arrive K1..K8.
 */
#include "fesom_kpp.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"

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
/* in-routine constants used in later phases (K3/K5/K7) */
#define KPP_RIINFTY      0.8        /* ri_iwmix shear-Ri shape limit (:737) */
#define KPP_MINMIX       3.0e-3     /* surface viscAE floor (driver :408)   */
#define KPP_CEKMAN       0.7        /* bldepth Ekman limit (:478)           */
#define KPP_CMONOB       1.0        /* bldepth Monin-Obukhov limit (:479)   */

/* wm/ws lookup-table index: Fortran wmt(0:nni+1, 0:nnj+1), i in [0,nni+1],
 * j in [0,nnj+1]. Row-major with i as the row (nnj+2 columns). */
#define KPP_TBL(i, j)  ((size_t)(i) * (size_t)(FESOM_KPP_NNJ + 2) + (size_t)(j))

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
 * KPP driver — abort stub until the routine bodies (K1..K8) are in place.
 * Selected only when FESOM_MIX_SCHEME=KPP (fesom_step.c); the default PP path
 * never reaches here, so K0 stays byte-identical to HEAD.
 *===========================================================================*/
void fesom_kpp_mixing(fesom_kpp               *k,
                      struct fesom_aux        *aux,
                      const struct fesom_tracers *tracers,
                      const struct fesom_dyn  *dyn,
                      const struct fesom_mesh *mesh,
                      struct fesom_partit     *partit)
{
    (void)k; (void)aux; (void)tracers; (void)dyn; (void)mesh; (void)partit;
    FESOM_DIE("fesom_kpp_mixing: KPP not yet implemented (Phase K0 scaffolding). "
              "Use FESOM_MIX_SCHEME=PP. KPP routines land in phases K1..K8.");
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
