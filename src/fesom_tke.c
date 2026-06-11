/*
 * gen_modules_cvmix_tke.F90 (module g_cvmix_tke) — TKE state + the
 * calc_cvmix_tke driver. Literal port; Fortran line refs inline.
 *
 * Loop-bound invariants (the project's signature bug class —
 * feedback_write_loops_halo cuts BOTH ways here):
 *   - the column loop runs OWNED nodes only (gen:296,311 node_size =
 *     myDim_nod2D) — do NOT extend to the halo;
 *   - the element-Av loop runs OWNED elements only (gen:500);
 *   - tke_Kv/tke_Av halos come from exchange_nod BEFORE the Kv copy /
 *     node→elem average; `tke` itself is NEVER exchanged;
 *   - NO element exchange in the driver — the element-Av halo comes from the
 *     shared post-mo_convect fesom_exchange_elem3D(aux->Av) in fesom_step.c.
 */
#include "fesom_tke.h"
#include "fesom_cvmix_tke.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* per-column stack scratch cap (= fesom_cvmix_tke.c TKE_NL_MAX) */
enum { TKE_NL_MAX = 128 };

#define TKE_CHECK(p, what) do {                                      \
    if (!(p)) {                                                      \
        fprintf(stderr, "fesom_tke_alloc: out of memory (%s)\n",     \
                (what));                                             \
        abort();                                                     \
    }                                                                \
} while (0)

/* fesom_tke_mixing call index (= ocean step, one call per step); mirrors the
 * Fortran tke_dump_mod counter. Used for dump-step gating. */
static int s_tke_call = 0;

/* shared per-column zero array for the gate-only inputs (tke_plc, E_iw,
 * alpha_c) AND iw_diss — the latter is READ unconditionally (:898) into
 * tke_Tiwf, so it must be real zeros, never garbage. */
static const real_t s_zero_col[TKE_NL_MAX + 1]; /* zero-initialized */

/*===========================================================================
 * dump harness (FESOM_TKE_DUMP_DIR) — the KPP/ALE dump convention
 *===========================================================================*/
static const char *s_tke_dump_dir   = NULL;
static int         s_tke_dump_steps = 3;
static int         s_tke_dump_env_loaded = 0;

static void tke_dump_load_env(void)
{
    if (s_tke_dump_env_loaded) return;
    s_tke_dump_dir = getenv("FESOM_TKE_DUMP_DIR");
    if (s_tke_dump_dir && !s_tke_dump_dir[0]) s_tke_dump_dir = NULL;
    const char *e = getenv("FESOM_TKE_DUMP_STEPS");
    if (e && atoi(e) > 0) s_tke_dump_steps = atoi(e);
    s_tke_dump_env_loaded = 1;
}

const char *fesom_tke_dump_dir(void)
{
    tke_dump_load_env();
    return s_tke_dump_dir;
}

int fesom_tke_dump_steps(void)
{
    tke_dump_load_env();
    return s_tke_dump_steps;
}

/* one [N*nl]-strided node field, all nl comps, owned nodes, gid-keyed */
static void tke_dump_nod_col(const struct fesom_mesh *mesh,
                             struct fesom_partit *partit,
                             const char *tag, const real_t *arr)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/tke_dump_s%d_%s_rank%d.txt",
             s_tke_dump_dir, s_tke_call, tag, partit->mype);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "tke_dump: cannot open %s\n", path); return; }
    int nl = mesh->nl;
    fprintf(f, "# step=%d tag=%s rank=%d N=%d ncomp=%d\n",
            s_tke_call, tag, partit->mype, mesh->myDim_nod2D, nl);
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        fprintf(f, "%d", partit->myList_nod2D[n]);   /* 1-based gid */
        for (int c = 0; c < nl; ++c)
            fprintf(f, " %.17g", (double)arr[FESOM_NODE3D(n, c, nl)]);
        fputc('\n', f);
    }
    fclose(f);
}

/* one [N] node scalar field */
static void tke_dump_nod_scalar(const struct fesom_mesh *mesh,
                                struct fesom_partit *partit,
                                const char *tag, const real_t *arr)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/tke_dump_s%d_%s_rank%d.txt",
             s_tke_dump_dir, s_tke_call, tag, partit->mype);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "tke_dump: cannot open %s\n", path); return; }
    fprintf(f, "# step=%d tag=%s rank=%d N=%d ncomp=%d\n",
            s_tke_call, tag, partit->mype, mesh->myDim_nod2D, 1);
    for (int n = 0; n < mesh->myDim_nod2D; ++n)
        fprintf(f, "%d %.17g\n", partit->myList_nod2D[n], (double)arr[n]);
    fclose(f);
}

/* one [E*nl]-strided element field, owned elements, gid-keyed */
static void tke_dump_elem_col(const struct fesom_mesh *mesh,
                              struct fesom_partit *partit,
                              const char *tag, const real_t *arr)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/tke_dump_s%d_%s_rank%d.txt",
             s_tke_dump_dir, s_tke_call, tag, partit->mype);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "tke_dump: cannot open %s\n", path); return; }
    int nl = mesh->nl;
    fprintf(f, "# step=%d tag=%s rank=%d N=%d ncomp=%d\n",
            s_tke_call, tag, partit->mype, mesh->myDim_elem2D, nl);
    for (int e = 0; e < mesh->myDim_elem2D; ++e) {
        fprintf(f, "%d", partit->myList_elem2D[e]);  /* 1-based gid */
        for (int c = 0; c < nl; ++c)
            fprintf(f, " %.17g", (double)arr[FESOM_ELEM3D(e, c, nl)]);
        fputc('\n', f);
    }
    fclose(f);
}

/*===========================================================================
 * controlled-input replay (FESOM_TKE_REPLAY_DIR) — the K6/K7 technique
 * (feedback_controlled_replay_validation). When set, the driver OVERRIDES
 * its per-column inputs (normstress, vshear2, bvfreq2, dz_trr, tke_old) with
 * the FORTRAN-dumped values for this call index, so the output diff is PURE
 * column algebra — isolated from the live forcing/EOS noise that threshold
 * flips amplify from step 2 on. Injecting tke_old also replaces the
 * recurrent state, making every step's comparison independent. The Fortran
 * rank-R dump lists rank R's owned nodes in myList order (same 16r
 * partition as the C), so file line i == local node i (asserted via gid).
 *===========================================================================*/
static void tke_replay_file(const char *dir, int step, const char *tag,
                            struct fesom_partit *partit,
                            int Nmy, int ncomp, real_t *out)
{
    char path[1024], line[16384];
    snprintf(path, sizeof(path), "%s/tke_dump_s%d_%s_rank%d.txt",
             dir, step, tag, partit->mype);
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "tke_replay: cannot open %s\n", path); abort(); }
    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "tke_replay: empty %s\n", path); abort();
    }
    for (int i = 0; i < Nmy; ++i) {
        if (!fgets(line, sizeof(line), f)) {
            fprintf(stderr, "tke_replay: short %s at %d\n", path, i); abort();
        }
        char *p = line, *e;
        long gid = strtol(p, &e, 10); p = e;
        if (gid != partit->myList_nod2D[i]) {
            fprintf(stderr, "tke_replay %s: gid/order mismatch at line %d "
                    "(%ld vs %d)\n", tag, i, gid, partit->myList_nod2D[i]);
            abort();
        }
        for (int c = 0; c < ncomp; ++c) {
            out[(size_t)i * ncomp + c] = (real_t)strtod(p, &e); p = e;
        }
    }
    fclose(f);
}

/*===========================================================================
 * Allocation — mirror of init_cvmix_tke's allocate block (gen:151-254) with
 * the reference &param_tke values (docs/tke_reference_namelists/), plus the
 * gate-only aborts. All arrays calloc'd (Fortran inits every slab to 0).
 *===========================================================================*/
void fesom_tke_alloc(fesom_tke *t, const struct fesom_mesh *mesh, int diag_on)
{
    memset(t, 0, sizeof(*t));
    int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;   /* gen:153 node_size */
    int nl = mesh->nl;
    t->n_nod   = N;
    t->nl      = nl;
    t->diag_on = diag_on;

    if (nl + 1 > TKE_NL_MAX) {
        fprintf(stderr, "fesom_tke_alloc: nl=%d exceeds TKE_NL_MAX\n", nl);
        abort();
    }

    size_t n_nl = (size_t)N * (size_t)nl;
    t->tke    = calloc(n_nl, sizeof(real_t));  TKE_CHECK(t->tke,    "tke");
    t->tke_Av = calloc(n_nl, sizeof(real_t));  TKE_CHECK(t->tke_Av, "tke_Av");
    t->tke_Kv = calloc(n_nl, sizeof(real_t));  TKE_CHECK(t->tke_Kv, "tke_Kv");

    t->forc_normstress = calloc((size_t)N, sizeof(real_t));
    t->forc_botfrict   = calloc((size_t)N, sizeof(real_t));
    t->forc_rhosurf    = calloc((size_t)N, sizeof(real_t));
    TKE_CHECK(t->forc_normstress, "forc_normstress");
    TKE_CHECK(t->forc_botfrict,   "forc_botfrict");
    TKE_CHECK(t->forc_rhosurf,    "forc_rhosurf");

    if (diag_on) {
        real_t **slabs[13] = { &t->Tbpr, &t->Tspr, &t->Tdif, &t->Tdis,
                               &t->Twin, &t->Tiwf, &t->Tbck, &t->Ttot,
                               &t->Lmix, &t->Pr,
                               &t->dummy1, &t->dummy2, &t->dummy3 };
        for (int i = 0; i < 13; ++i) {
            *slabs[i] = calloc(n_nl, sizeof(real_t));
            TKE_CHECK(*slabs[i], "tke diag slab");
        }
    }

    /* init_tke param plumbing (gen:258-272) with the namelist.cvmix
     * reference values — tke_cd = 3.75 is the NAMELIST value (module default
     * 1.0 loses; feedback_namelist_over_codedefault). tke_clangmuir keeps
     * its module default 0.3 (not in the namelist; gate-only anyway). */
    fesom_cvmix_init_tke(/*c_k*/        0.1,
                         /*c_eps*/      0.7,
                         /*cd*/         3.75,
                         /*alpha_tke*/  30.0,
                         /*mxl_min*/    1.0e-8,
                         /*kappaM_min*/ 0.0,
                         /*kappaM_max*/ 100.0,
                         /*tke_mxl_choice*/ 2,
                         /*use_ubound_dirichlet*/ 0,
                         /*use_lbound_dirichlet*/ 0,
                         /*only_tke*/   1,
                         /*l_lc*/       0,
                         /*clc*/        0.3,
                         /*tke_min*/    1.0e-6,
                         /*tke_surf_min*/ 1.0e-4);

    /* gate-only aborts: the port covers exactly the reference configuration
     * (plan "Out of scope"). If a future edit flips one of these constants,
     * fail loudly rather than run un-ported physics. */
    const fesom_cvmix_tke_params *p = fesom_cvmix_tke_constants();
    if (!p->only_tke || p->l_lc || p->use_ubound_dirichlet ||
        p->use_lbound_dirichlet || p->tke_mxl_choice != 2) {
        fprintf(stderr, "fesom_tke_alloc: configuration outside the ported "
                "path (IDEMIX/Langmuir/Dirichlet/mxl_choice!=2 are "
                "gate-only, not ported)\n");
        abort();
    }
}

void fesom_tke_free(fesom_tke *t)
{
    free(t->tke); free(t->tke_Av); free(t->tke_Kv);
    free(t->forc_normstress); free(t->forc_botfrict); free(t->forc_rhosurf);
    free(t->Tbpr); free(t->Tspr); free(t->Tdif); free(t->Tdis);
    free(t->Twin); free(t->Tiwf); free(t->Tbck); free(t->Ttot);
    free(t->Lmix); free(t->Pr);
    free(t->dummy1); free(t->dummy2); free(t->dummy3);
    memset(t, 0, sizeof(*t));
}

/*===========================================================================
 * calc_cvmix_tke (gen:279-507)
 *===========================================================================*/
void fesom_tke_mixing(fesom_tke                  *t,
                      struct fesom_aux           *aux,
                      const struct fesom_forcing *forcing,
                      const struct fesom_dyn     *dyn,
                      const struct fesom_mesh    *mesh,
                      struct fesom_partit        *partit)
{
    const int N_own  = mesh->myDim_nod2D;            /* gen:296 node_size */
    const int N_full = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int E_own  = mesh->myDim_elem2D;
    const int E_full = mesh->myDim_elem2D + mesh->eDim_elem2D;
    const int nl     = mesh->nl;

    real_t dz_trr[TKE_NL_MAX], bvfreq2[TKE_NL_MAX], vshear2[TKE_NL_MAX];
    real_t tke_old[TKE_NL_MAX];                      /* gen:287 column copy */

    ++s_tke_call;
    tke_dump_load_env();
    const int dmp_on = (s_tke_dump_dir != NULL && s_tke_call <= s_tke_dump_steps);

    /* gen:299-301 — zero the 2D forcing arrays every call (full extent) */
    memset(t->forc_normstress, 0, (size_t)N_full * sizeof(real_t));
    memset(t->forc_botfrict,   0, (size_t)N_full * sizeof(real_t));
    memset(t->forc_rhosurf,    0, (size_t)N_full * sizeof(real_t));

    /* gen:304-308 — IDEMIX input load: gate-only (only_tke pinned true). */

    /* dump staging for the per-column INPUTS (allocated only when dumping) */
    real_t *dmp_vsh = NULL, *dmp_bvf = NULL, *dmp_dzt = NULL, *dmp_tko = NULL;
    if (dmp_on) {
        size_t n_nl = (size_t)N_own * (size_t)nl;
        dmp_vsh = calloc(n_nl, sizeof(real_t));
        dmp_bvf = calloc(n_nl, sizeof(real_t));
        dmp_dzt = calloc(n_nl, sizeof(real_t));
        dmp_tko = calloc(n_nl, sizeof(real_t));
    }

    /* controlled replay: load this call's Fortran-dumped inputs */
    const char *rep_dir = getenv("FESOM_TKE_REPLAY_DIR");
    real_t *rep_norm = NULL, *rep_vsh = NULL, *rep_bvf = NULL,
           *rep_dzt = NULL, *rep_tko = NULL;
    const int rep_on = (rep_dir && rep_dir[0] && s_tke_call <= s_tke_dump_steps);
    if (rep_on) {
        size_t n_nl = (size_t)N_own * (size_t)nl;
        rep_norm = malloc((size_t)N_own * sizeof(real_t));
        rep_vsh  = malloc(n_nl * sizeof(real_t));
        rep_bvf  = malloc(n_nl * sizeof(real_t));
        rep_dzt  = malloc(n_nl * sizeof(real_t));
        rep_tko  = malloc(n_nl * sizeof(real_t));
        tke_replay_file(rep_dir, s_tke_call, "normstress", partit, N_own, 1,  rep_norm);
        tke_replay_file(rep_dir, s_tke_call, "vshear2",    partit, N_own, nl, rep_vsh);
        tke_replay_file(rep_dir, s_tke_call, "bvfreq2",    partit, N_own, nl, rep_bvf);
        tke_replay_file(rep_dir, s_tke_call, "dztrr",      partit, N_own, nl, rep_dzt);
        tke_replay_file(rep_dir, s_tke_call, "tkeold",     partit, N_own, nl, rep_tko);
        if (partit->mype == 0)
            fprintf(stderr, "[tke_replay] s%d: injected Fortran inputs from %s "
                    "(output diff = pure column algebra)\n", s_tke_call, rep_dir);
    }

    /*------------------------------------------------------------------
     * the per-node column loop — OWNED nodes only (gen:311)
     *-----------------------------------------------------------------*/
    for (int node = 0; node < N_own; ++node) {
        /* gen:314-315 — nln = nlevels-1 (last wet layer), nun = ulevels.
         * 0-based: uln0 = first wet interface, nln0 = last wet layer. */
        const int nln0 = mesh->nlevels_nod2D[node] - 2;
        const int uln0 = mesh->ulevels_nod2D[node] - 1;
        const int nlev = nln0 - uln0 + 1;            /* layer count at the call */

        /* gen:330 — surface momentum forcing from the NODAL wind stress
         * (the element-interpolation variant is commented out). */
        {
            real_t sx = forcing->stress_node_surf[2 * node + 0];
            real_t sy = forcing->stress_node_surf[2 * node + 1];
            t->forc_normstress[node] =
                sqrt(sx * sx + sy * sy) / (real_t)FESOM_DENSITY_0;
        }

        /* gen:335-340 — vertical velocity shear² at interfaces nun+1..nln */
        for (int k = 0; k < nl; ++k) vshear2[k] = 0.0;
        for (int nz = uln0 + 1; nz <= nln0; ++nz) {
            real_t du = dyn->uvnode[FESOM_ELEMVEC(node, nz - 1, nl) + 0]
                      - dyn->uvnode[FESOM_ELEMVEC(node, nz,     nl) + 0];
            real_t dv = dyn->uvnode[FESOM_ELEMVEC(node, nz - 1, nl) + 1]
                      - dyn->uvnode[FESOM_ELEMVEC(node, nz,     nl) + 1];
            real_t dz = mesh->Z_3d_n[FESOM_NODE3D(node, nz - 1, nl)]
                      - mesh->Z_3d_n[FESOM_NODE3D(node, nz,     nl)];
            vshear2[nz] = (du * du + dv * dv) / (dz * dz);
        }

        /* gen:352-354 — N² (bvfreq already holds the squared frequency,
         * horizontally smoothed upstream) */
        for (int k = 0; k < nl; ++k) bvfreq2[k] = 0.0;
        for (int nz = uln0 + 1; nz <= nln0; ++nz)
            bvfreq2[nz] = aux->bvfreq[FESOM_NODE3D(node, nz, nl)];

        /* gen:359-362 — dz_trr: distance between tracer points; surface and
         * bottom entries are half the layer thickness */
        for (int k = 0; k < nl; ++k) dz_trr[k] = 0.0;
        for (int nz = uln0 + 1; nz <= nln0; ++nz)
            dz_trr[nz] = fabs(mesh->Z_3d_n[FESOM_NODE3D(node, nz - 1, nl)]
                            - mesh->Z_3d_n[FESOM_NODE3D(node, nz,     nl)]);
        dz_trr[uln0]     = mesh->hnode[FESOM_NODE3D(node, uln0, nl)] / 2.0;
        dz_trr[nln0 + 1] = mesh->hnode[FESOM_NODE3D(node, nln0, nl)] / 2.0;

        /* gen:367-429 — Langmuir block: gate-only (tke_dolangmuir=F). */

        /* gen:433-435 — per-column copy of tke BEFORE the call (the ONLY
         * persistent input; tke_Av_old/tke_Kv_old fed verified-dead args
         * and are NOT ported). */
        for (int k = 0; k < nl; ++k)
            tke_old[k] = t->tke[FESOM_NODE3D(node, k, nl)];

        /* controlled replay: override every column input with the Fortran's */
        if (rep_on) {
            t->forc_normstress[node] = rep_norm[node];
            for (int k = 0; k < nl; ++k) {
                vshear2[k] = rep_vsh[FESOM_NODE3D(node, k, nl)];
                bvfreq2[k] = rep_bvf[FESOM_NODE3D(node, k, nl)];
                dz_trr[k]  = rep_dzt[FESOM_NODE3D(node, k, nl)];
                tke_old[k] = rep_tko[FESOM_NODE3D(node, k, nl)];
            }
        }

        if (dmp_on) {
            for (int k = 0; k < nl; ++k) {
                dmp_vsh[FESOM_NODE3D(node, k, nl)] = vshear2[k];
                dmp_bvf[FESOM_NODE3D(node, k, nl)] = bvfreq2[k];
                dmp_dzt[FESOM_NODE3D(node, k, nl)] = dz_trr[k];
                dmp_tko[FESOM_NODE3D(node, k, nl)] = tke_old[k];
            }
        }

        /* gen:437-482 — the cvmix call with nun:nln(+1) slices,
         * nlev = nln-nun+1. Outputs write straight into the global slabs
         * (tke_new/KappaM_out/KappaH_out were slices of tke/tke_Av/tke_Kv). */
        fesom_cvmix_tke_diag  dcols;
        fesom_cvmix_tke_diag *dptr = NULL;
        if (t->diag_on) {
            size_t off = FESOM_NODE3D(node, uln0, nl);
            dcols.Tbpr = t->Tbpr + off;  dcols.Tspr = t->Tspr + off;
            dcols.Tdif = t->Tdif + off;  dcols.Tdis = t->Tdis + off;
            dcols.Twin = t->Twin + off;  dcols.Tiwf = t->Tiwf + off;
            dcols.Tbck = t->Tbck + off;  dcols.Ttot = t->Ttot + off;
            dcols.Lmix = t->Lmix + off;  dcols.Pr   = t->Pr   + off;
            dcols.int1 = t->dummy1 + off;
            dcols.int2 = t->dummy2 + off;
            dcols.int3 = t->dummy3 + off;
            dptr = &dcols;
        }

        fesom_cvmix_integrate_tke(
            nlev,
            (real_t)FESOM_PHASE1_DT,                   /* dtime   (gen:444) */
            (real_t)FESOM_DENSITY_0,                   /* rho_ref (gen:445) */
            (real_t)FESOM_G,                           /* grav    (gen:446) */
            &mesh->hnode[FESOM_NODE3D(node, uln0, nl)],/* dzw = hnode(nun:nln) */
            &dz_trr[uln0],                             /* dzt(nun:nln+1)   */
            &tke_old[uln0],                            /* tke_old          */
            &vshear2[uln0],                            /* Ssqr             */
            &bvfreq2[uln0],                            /* Nsqr             */
            s_zero_col,                                /* alpha_c (gated)  */
            s_zero_col,                                /* E_iw    (gated)  */
            s_zero_col,                                /* iw_diss — REAL zeros,
                                                          read unconditionally */
            s_zero_col,                                /* tke_plc (gated)  */
            t->forc_normstress[node],                  /* forc_tke_surf    */
            t->forc_rhosurf[node],                     /* forc_rho_surf    */
            t->forc_botfrict[node],                    /* bottom_fric      */
            &t->tke   [FESOM_NODE3D(node, uln0, nl)],  /* tke_new          */
            &t->tke_Av[FESOM_NODE3D(node, uln0, nl)],  /* KappaM_out       */
            &t->tke_Kv[FESOM_NODE3D(node, uln0, nl)],  /* KappaH_out       */
            dptr);

        /* gen:484-487 — zero visc/diff at the surface and below-bottom
         * interfaces of the column */
        t->tke_Av[FESOM_NODE3D(node, nln0 + 1, nl)] = 0.0;
        t->tke_Kv[FESOM_NODE3D(node, nln0 + 1, nl)] = 0.0;
        t->tke_Av[FESOM_NODE3D(node, uln0,     nl)] = 0.0;
        t->tke_Kv[FESOM_NODE3D(node, uln0,     nl)] = 0.0;
    }

    /* dump: inputs + post-loop outputs (pre-exchange; owned rows only, so
     * exchange state is irrelevant) + the 10 budget diagnostics (read from
     * the stored slabs — diag_on is forced on when the dump dir is set). */
    if (dmp_on) {
        tke_dump_nod_scalar(mesh, partit, "normstress", t->forc_normstress);
        tke_dump_nod_col(mesh, partit, "vshear2", dmp_vsh);
        tke_dump_nod_col(mesh, partit, "bvfreq2", dmp_bvf);
        tke_dump_nod_col(mesh, partit, "dztrr",   dmp_dzt);
        tke_dump_nod_col(mesh, partit, "tkeold",  dmp_tko);
        tke_dump_nod_col(mesh, partit, "tke",   t->tke);
        tke_dump_nod_col(mesh, partit, "tkeav", t->tke_Av);
        tke_dump_nod_col(mesh, partit, "tkekv", t->tke_Kv);
        if (t->diag_on) {
            tke_dump_nod_col(mesh, partit, "tbpr", t->Tbpr);
            tke_dump_nod_col(mesh, partit, "tspr", t->Tspr);
            tke_dump_nod_col(mesh, partit, "tdif", t->Tdif);
            tke_dump_nod_col(mesh, partit, "tdis", t->Tdis);
            tke_dump_nod_col(mesh, partit, "twin", t->Twin);
            tke_dump_nod_col(mesh, partit, "tiwf", t->Tiwf);
            tke_dump_nod_col(mesh, partit, "tbck", t->Tbck);
            tke_dump_nod_col(mesh, partit, "ttot", t->Ttot);
            tke_dump_nod_col(mesh, partit, "lmix", t->Lmix);
            tke_dump_nod_col(mesh, partit, "pr",   t->Pr);
        }
    }
    free(dmp_vsh); free(dmp_bvf); free(dmp_dzt); free(dmp_tko);
    free(rep_norm); free(rep_vsh); free(rep_bvf); free(rep_dzt); free(rep_tko);

    /*------------------------------------------------------------------
     * gen:491-505 — wire the outputs
     *-----------------------------------------------------------------*/
    /* diffusivity: exchange THEN full-array copy (halo included) */
    fesom_exchange_nod3D(t->tke_Kv, nl, partit);     /* gen:493 */
    memcpy(aux->Kv, t->tke_Kv,
           (size_t)N_full * (size_t)nl * sizeof(real_t));  /* gen:494 Kv = tke_Kv */

    /* viscosity: exchange nodes, then node→elem mean on OWNED elements over
     * interior levels; Av zeroed over its FULL extent first (gen:499 Av=0).
     * Halo Av is restored by the shared post-mo_convect element exchange. */
    fesom_exchange_nod3D(t->tke_Av, nl, partit);     /* gen:498 */
    memset(aux->Av, 0, (size_t)E_full * (size_t)nl * sizeof(real_t));
    for (int elem = 0; elem < E_own; ++elem) {       /* gen:500 */
        int n0 = mesh->elem_nodes[3 * elem + 0];
        int n1 = mesh->elem_nodes[3 * elem + 1];
        int n2 = mesh->elem_nodes[3 * elem + 2];
        int ule0 = mesh->ulevels[elem] - 1;
        int nle0 = mesh->nlevels[elem] - 1;
        for (int nz = ule0 + 1; nz <= nle0 - 1; ++nz) {  /* gen:502 */
            aux->Av[FESOM_ELEM3D(elem, nz, nl)] =
                (t->tke_Av[FESOM_NODE3D(n0, nz, nl)]
               + t->tke_Av[FESOM_NODE3D(n1, nz, nl)]
               + t->tke_Av[FESOM_NODE3D(n2, nz, nl)]) / 3.0;
        }
    }

    if (dmp_on) {
        tke_dump_nod_col (mesh, partit, "kv", aux->Kv);
        tke_dump_elem_col(mesh, partit, "av", aux->Av);
    }

    /* exchange-and-compare halo probe (FESOM_TKE_HALO_PROBE=1, T4 gate;
     * feedback_write_loops_halo): aux->Kv halo must already be owner-fresh
     * here (copied from the exchanged tke_Kv) — re-exchanging a copy must
     * change NOTHING. A nonzero probe = a stale-halo write loop upstream. */
    if (getenv("FESOM_TKE_HALO_PROBE")) {
        size_t n_nl = (size_t)N_full * (size_t)nl;
        real_t *probe = malloc(n_nl * sizeof(real_t));
        if (probe) {
            memcpy(probe, aux->Kv, n_nl * sizeof(real_t));
            fesom_exchange_nod3D(probe, nl, partit);
            real_t worst = 0.0;
            for (size_t i = 0; i < n_nl; ++i) {
                real_t d = probe[i] - aux->Kv[i];
                if (d < 0) d = -d;
                if (d > worst) worst = d;
            }
            fprintf(stderr, "[tke_halo_probe r%d] s%d Kv re-exchange worst |d| = %.3e%s\n",
                    partit->mype, s_tke_call, (double)worst,
                    worst > 0.0 ? "  <-- STALE HALO" : "");
            free(probe);
        }
    }
}
