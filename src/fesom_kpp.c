/*
 * KPP vertical mixing — faithful port of oce_ale_mixing_kpp.F90.
 * See fesom_kpp.h for the struct layout and docs/plans/20260524-kpp-vertical-mixing.md
 * for the phased port plan. Phase K0 lands the scaffolding (struct alloc/free,
 * abort-stub driver, dump harness); the routine bodies arrive K1..K8.
 */
#include "fesom_kpp.h"
#include "fesom_aux.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
 * One-time init — Vtc, cg, deltaz/deltau, wm/ws lookup tables
 * (oce_mixing_kpp_init, :157-209). PHASE K1 — body lands next.
 *===========================================================================*/
void fesom_kpp_init(fesom_kpp *k)
{
    (void)k;
    /* Phase K1: build wmt/wst (nni=890, nnj=480) + derive Vtc (:167), cg (:178).
     * Stub for K0 (KPP path aborts before this matters). */
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
