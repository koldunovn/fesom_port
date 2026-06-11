#ifndef FESOM_TKE_H
#define FESOM_TKE_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_partit;
struct fesom_aux;
struct fesom_forcing;
struct fesom_dyn;

/*
 * CVMix classical-TKE vertical mixing (Gaspar et al. 1990; Blanke & Delecluse
 * 1993 mixing length) — strictly faithful port of the ICON/Brüggemann CVMix
 * fork FESOM2 vendors in src/cvmix_driver/:
 *
 *   fesom_tke.{c,h}       = gen_modules_cvmix_tke.F90 (module g_cvmix_tke):
 *                           state arrays, namelist params, the per-node driver
 *                           calc_cvmix_tke, Av/Kv wiring, MPI exchanges.
 *   fesom_cvmix_tke.{c,h} = cvmix_tke.F90 (module cvmix_tke): init_tke +
 *                           integrate_tke, the pure column math (no MPI/mesh).
 *
 * Spec: docs/plans/20260610-tke-vertical-mixing.md. Reference namelists:
 * docs/tke_reference_namelists/ (work_linfs_tke; mix_scheme='cvmix_TKE',
 * dt=1800, namelist.cvmix &param_tke with tke_cd=3.75).
 *
 * Selected at runtime via FESOM_MIX_SCHEME=TKE (alias cvmix_TKE); KPP stays
 * the default, PP stays available. The TKE branch mirrors oce_ale.F90:3749-52
 * (calc_cvmix_tke → mo_convect → shared Kv/Av exchanges in fesom_step.c).
 *
 * State (mirrors the g_cvmix_tke module allocatables, gen:130-273):
 *   Always allocated, [N*nl] stride nl, N = myDim+eDim:
 *     tke      — PROGNOSTIC TKE field (init 0 → floored to tke_min by the
 *                first call; per-column tke_old is read each step). The
 *                Fortran comment calls it a "diagnostic" — it is not.
 *     tke_Av/tke_Kv — node visc/diff; global only because exchange_nod + the
 *                full Kv copy need halo storage. Prior-step values are NEVER
 *                read (no old-value blending — `handle_old_vals` and the
 *                old_KappaM/old_KappaH args are verified DEAD, see
 *                fesom_cvmix_tke.h).
 *   2D [N]: forc_normstress / forc_botfrict / forc_rhosurf (botfrict/rhosurf
 *                zeroed every call and passed — as the Fortran does).
 *   Gated by diag_on (FESOM_TKE_DIAG=1, or FESOM_TKE_DUMP_DIR set — the dump
 *   reads the stored diag slabs): the 10 budget arrays Tbpr/Tspr/Tdif/Tdis/
 *   Twin/Tiwf/Tbck/Ttot/Lmix/Pr + 3 cvmix_dummy_*, each [N*nl]. The column
 *   core ALWAYS computes them into local scratch (internal read-after-write
 *   deps); the flag gates only global storage. Diag-on vs diag-off model
 *   state is byte-identical (T3 gate).
 *
 *   NOT allocated (zero/dead by construction under only_tke=T, dolangmuir=F):
 *   tke_in3d_iwe/iwdis/iwealphac, tke_langmuir, langmuir_* — the column call
 *   passes a shared per-column zero array instead (value-identical; iw_diss
 *   is read unconditionally into tke_Tiwf, so it must be a REAL zero array).
 *   Module-level 1D tke_Av_old/tke_Kv_old/tke_old allocatables + tke_diss:
 *   shadowed dead code, not ported.
 */
typedef struct fesom_tke {
    int n_nod;          /* myDim+eDim nodes */
    int nl;             /* layer interfaces */
    int diag_on;        /* 1 → the 13 diag slabs below are allocated+stored */

    /* always allocated [N*nl] */
    real_t *tke;        /* prognostic TKE at interfaces */
    real_t *tke_Av;     /* node viscosity  (KappaM) */
    real_t *tke_Kv;     /* node diffusivity (KappaH) */

    /* always allocated [N] */
    real_t *forc_normstress;  /* |stress_node_surf|/rho0 */
    real_t *forc_botfrict;    /* zeroed every call (no function yet — Fortran) */
    real_t *forc_rhosurf;     /* zeroed every call (no function yet — Fortran) */

    /* diag slabs [N*nl] — NULL unless diag_on (gen_modules allocatables) */
    real_t *Tbpr;       /* buoyancy production */
    real_t *Tspr;       /* shear production */
    real_t *Tdif;       /* vertical diffusion of TKE */
    real_t *Tdis;       /* dissipation */
    real_t *Twin;       /* wind forcing */
    real_t *Tiwf;       /* internal-wave forcing (≡0 under only_tke) */
    real_t *Tbck;       /* background restoring (tke_min reset) */
    real_t *Ttot;       /* total tendency */
    real_t *Lmix;       /* mixing length */
    real_t *Pr;         /* Prandtl number */
    real_t *dummy1;     /* cvmix_dummy_1 = KappaH_out */
    real_t *dummy2;     /* cvmix_dummy_2 = KappaM_out */
    real_t *dummy3;     /* cvmix_dummy_3 = Nsqr */
} fesom_tke;

/* Allocate + zero the state (mirror of init_cvmix_tke's allocate block,
 * gen:151-207). diag_on additionally allocates the 13 diag slabs. Also runs
 * the C equivalent of the init_tke param plumbing (fesom_cvmix_tke.c) with
 * the reference namelist.cvmix values — tke_cd=3.75 (NAMELIST value; the
 * module default 1.0 loses, feedback_namelist_over_codedefault) — and
 * aborts if a gate-only option would be required (Langmuir, Dirichlet BCs,
 * IDEMIX). */
void fesom_tke_alloc(fesom_tke *t, const struct fesom_mesh *mesh, int diag_on);
void fesom_tke_free (fesom_tke *t);

/*
 * TKE driver — mirror of calc_cvmix_tke (gen_modules_cvmix_tke.F90:279-507).
 * Per OWNED node (gen:296,311 node_size=myDim_nod2D): build normstress,
 * vshear2 (UVnode/Z_3d_n), bvfreq2, dz_trr, per-column tke_old copy; call
 * integrate_tke with nun:nln+1 slices; zero tke_Av/tke_Kv at nun and nln+1.
 * Then exchange_nod(tke_Kv) → aux->Kv full copy; exchange_nod(tke_Av) →
 * node→elem mean into aux->Av over OWNED elements (gen:500), interior levels
 * only. NO element exchange here — the element-Av halo comes from the shared
 * post-mo_convect fesom_exchange_elem3D(aux->Av) in fesom_step.c, exactly as
 * for KPP/PP. tke itself is NEVER exchanged (no reader needs its halo).
 *
 * Abort stub until Phase T2 wires the body.
 */
void fesom_tke_mixing(fesom_tke                  *t,
                      struct fesom_aux           *aux,
                      const struct fesom_forcing *forcing,
                      const struct fesom_dyn     *dyn,
                      const struct fesom_mesh    *mesh,
                      struct fesom_partit        *partit);

/* ---- dump harness (FESOM_TKE_DUMP_DIR-gated dual-instrumentation) --------
 * Mirror of the KPP/ALE dump convention. Off by default (zero overhead).
 * When FESOM_TKE_DUMP_DIR is set, fesom_tke_mixing dumps per-node columns at
 * TKE calls 1..FESOM_TKE_DUMP_STEPS (default 3), to be diffed against the
 * Fortran tke_dump_mod gate by scripts/tke_dump_diff.py. Tags (= Fortran):
 *   inputs : normstress(1) vshear2 bvfreq2 dztrr tkeold        (nl)
 *   outputs: tke tkeav tkekv                                    (nl)
 *   diag   : tbpr tspr tdif tdis twin tiwf tbck ttot lmix pr    (nl)
 *   wired  : kv (nodes, nl)  av (elements, nl)
 * Format: tke_dump_s<step>_<tag>_rank<R>.txt, header
 * "# step= tag= rank= N= ncomp=", gid-keyed rows. */
const char *fesom_tke_dump_dir(void);
int         fesom_tke_dump_steps(void);

#endif /* FESOM_TKE_H */
