#ifndef FESOM_KPP_H
#define FESOM_KPP_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_partit;
struct fesom_aux;
struct fesom_tracers;
struct fesom_forcing;
struct fesom_dyn;

/*
 * KPP (K-Profile Parameterization) vertical mixing — strictly faithful port
 * of Fortran o_mixing_KPP_mod (oce_ale_mixing_kpp.F90, 1185 lines). This is
 * the vertical mixing scheme production CORE2 actually uses; the C also keeps
 * PP (fesom_pp.c) and selects at runtime via FESOM_MIX_SCHEME (see fesom_step.c).
 *
 * Spec: docs/plans/20260524-kpp-vertical-mixing.md. Reference namelists:
 * docs/kpp_reference_namelists/ (work_core; mix_scheme='KPP', dt=1800).
 *
 * This struct mirrors the o_mixing_KPP_mod module-level allocatables and
 * scalars (the Fortran routines USE them; we thread the struct, mirroring the
 * fesom_gm precedent). One exception: dbsfc lives in fesom_aux (it is written
 * by the EOS — fesom_eos.c — which only carries aux; Fortran's dbsfc is the
 * o_mixing_KPP_mod member filled from oce_ale_pressure_bv.F90 via USE).
 *
 * Memory layout — flat, matching Fortran column-major exactly so that
 * per-component slabs (blmc(:,:,j), diffK(:,:,tr)) are contiguous [N*nl] arrays
 * directly passable to fesom_exchange_nod3D (the driver exchanges them so).
 *   N  = n_nod = myDim_nod2D + eDim_nod2D   (halo coverage required)
 *   nl = number of layer interfaces
 *   Fortran (nz, node)        -> C [node*nl + nz]               (= FESOM_NODE3D)
 *   Fortran (nz, node, j)     -> C [j*N*nl + node*nl + nz]      (j-major slabs)
 *   Fortran (node, j)         -> C [j*N + node]                 (j-major)
 *   Fortran ghats(nz, node)   -> C [node*(nl-1) + nz]   (nl-1 levels, NOT nl!)
 *
 *   diffK [ntr * N * nl]   per-tracer diffusivity (T=ch0, S=ch1) — Kv_double
 *   viscA [N * nl]         node viscosity scratch (averaged to elem -> Av)
 *   blmc  [3 * N * nl]     boundary-layer mixing coeffs (viscA/Kv_T/Kv_S)
 *   dkm1  [3 * N]          BL coeffs at kbl-1 (consumed by enhance)
 *   ghats [N * (nl-1)]     nonlocal transport (computed, gated OFF in CORE2)
 *   dVsq  [N * nl]         (velocity shear re surface)^2
 *   hbl,bfsfc,caseA,stable,ustar,Bo  [N]
 *   kbl   [N] int          first grid level below hbl
 *   wmt,wst [(nni+2)*(nnj+2)]  turbulent-velocity-scale lookup tables
 *
 * ghats stride is (nl-1) because Fortran allocates ghats(nl-1, N) (the only
 * KPP array on the mid-layer count, not nl) — see feedback_tracer_stride_nl
 * rule 3 (nl-1 stride is correct only for arrays explicitly sized [N*(nl-1)]).
 */

/* lookup-table dimensions (oce_ale_mixing_kpp.F90:69-70) */
#define FESOM_KPP_NNI 890   /* zehat axis  */
#define FESOM_KPP_NNJ 480   /* ustar axis  */

typedef struct fesom_kpp {
    /* dimensions (cached for slab strides) */
    int n_nod;          /* myDim+eDim nodes */
    int nl;             /* layer interfaces */
    int num_tracers;    /* diffK channel count */

    /* per-tracer / per-component 3D node fields */
    real_t *diffK;      /* [ntr * N * nl]  T=ch0, S=ch1 (Kv_double) */
    real_t *viscA;      /* [N * nl]        node viscosity */
    real_t *blmc;       /* [3 * N * nl]    BL mixing coeffs */
    real_t *ghats;      /* [N * (nl-1)]    nonlocal transport */
    real_t *dVsq;       /* [N * nl]        shear^2 re surface */

    /* per-node 2D fields */
    real_t *dkm1;       /* [3 * N]         BL coeffs at kbl-1 */
    real_t *hbl;        /* [N]             boundary-layer depth */
    real_t *bfsfc;      /* [N]             surface buoyancy forcing */
    real_t *caseA;      /* [N]             =1 case A, =0 case B */
    real_t *stable;     /* [N]             =1 stable forcing, =0 unstable */
    real_t *ustar;      /* [N]             surface friction velocity */
    real_t *Bo;         /* [N]             surface turb buoyancy forcing */
    int    *kbl;        /* [N]             first grid level below hbl */

    /* turbulent-velocity-scale lookup tables (built in fesom_kpp_init) */
    real_t *wmt;        /* [(nni+2)*(nnj+2)] momentum scale */
    real_t *wst;        /* [(nni+2)*(nnj+2)] scalar scale */

    /* init-time scalars */
    real_t Vtc;         /* eqn 23 velocity-scale coefficient */
    real_t cg;          /* eqn 20 nonlocal-transport coefficient */
    real_t deltaz;      /* zehat table step */
    real_t deltau;      /* ustar table step */
} fesom_kpp;

void fesom_kpp_alloc(fesom_kpp *k, const struct fesom_mesh *mesh, int num_tracers);
void fesom_kpp_free (fesom_kpp *k);

/* One-time init: lookup tables + Vtc/cg + constants (Fortran
 * oce_mixing_kpp_init). Implemented in Phase K1; a no-op-safe stub until then. */
void fesom_kpp_init(fesom_kpp *k);

/*
 * KPP driver — mirror of oce_mixing_KPP (oce_ale_mixing_kpp.F90:249-425).
 * Fills aux->Av (element viscosity = viscAE) and k->diffK (per-tracer node
 * diffusivity = Kv_double); the caller copies the single aux->Kv = diffK
 * T-channel (oce_ale.F90:3518-3522) and runs mo_convect afterwards.
 *
 * Implemented incrementally K1..K8; an abort stub until then so the dispatch
 * fails loudly if KPP is selected before the port is complete.
 */
void fesom_kpp_mixing(fesom_kpp                  *k,
                      struct fesom_aux           *aux,
                      const struct fesom_tracers *tracers,
                      const struct fesom_forcing *forcing,
                      const struct fesom_dyn     *dyn,
                      const struct fesom_mesh    *mesh,
                      struct fesom_partit        *partit);

/* ---- dump harness (FESOM_KPP_DUMP_DIR-gated dual-instrumentation) --------
 * Mirror of the EVP dump diagnostic (reference_evp_dump_diagnostic). Off by
 * default (zero overhead). When FESOM_KPP_DUMP_DIR is set, the KPP routines
 * dump named per-node fields keyed by 1-based global id at a chosen step, to
 * be diffed against the Fortran gate by scripts/kpp_dump_diff.py. */

/* Returns the dump dir (cached) or NULL if FESOM_KPP_DUMP_DIR unset. */
const char *fesom_kpp_dump_dir(void);

/* K1 init dump: rank 0 writes <dir>/kpp_init_rank0.txt — the 4 init scalars
 * (Vtc, cg, deltaz, deltau) + the full wm/ws lookup table (`i j wmt wst`).
 * No-op unless FESOM_KPP_DUMP_DIR is set. Compared vs the Fortran gate +
 * scripts/kpp_init_reference.py. */
void fesom_kpp_dump_init(const fesom_kpp *k, int rank);

/* K2 wscale sweep dump: rank 0 writes <dir>/kpp_wscale_rank0.txt — wm/ws over a
 * fixed (zehat, ustar) grid. No-op unless FESOM_KPP_DUMP_DIR is set. */
void fesom_kpp_dump_wscale_sweep(const fesom_kpp *k, int rank);

/* True on the step we dump (default step 1; override FESOM_KPP_DUMP_STEP). */
int fesom_kpp_dump_this_step(int step_n);

/*
 * Dump `ncomp` real fields per owned node to
 *   <dir>/kpp_dump_s<step>_<tag>_rank<R>.txt
 * one line per owned node: "<gid> <v0> <v1> ...". `get(node, comp, user)`
 * returns component `comp` of field at local `node`. gids are 1-based
 * (mesh->myList_nod2D). Used by all phase dumps; comp layout is per-tag and
 * mirrored in scripts/kpp_dump_diff.py.
 */
void fesom_kpp_dump_nodes(const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          int step_n, const char *tag, int ncomp,
                          double (*get)(int node, int comp, void *user),
                          void *user);

/* As fesom_kpp_dump_nodes but per owned ELEMENT (keyed by 1-based element gid,
 * mesh->myList_elem2D); used for the K7 viscAE (=aux->Av) element field. */
void fesom_kpp_dump_elems(const struct fesom_mesh *mesh,
                          struct fesom_partit     *partit,
                          int step_n, const char *tag, int ncomp,
                          double (*get)(int elem, int comp, void *user),
                          void *user);

#endif /* FESOM_KPP_H */
