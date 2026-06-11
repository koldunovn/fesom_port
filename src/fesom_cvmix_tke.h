#ifndef FESOM_CVMIX_TKE_H
#define FESOM_CVMIX_TKE_H

#include "fesom_types.h"

/*
 * cvmix_tke.F90 (module cvmix_tke) — the classical-TKE column core. Literal
 * port: init_tke (:101-281, namelist params → the saved constants struct) and
 * integrate_tke (:415-987, the ~570-line prognostic TKE column step). Pure
 * column math: no MPI, no mesh — the driver (fesom_tke.c) assembles columns.
 *
 * == Arguments DROPPED from the Fortran signature (verified dead; plan
 *    "Dropped as dead", plan-review 2026-06-10) ==
 *  - old_KappaM / old_KappaH: declared intent(in) (:474-475), never appear in
 *    the body. There is NO old-value blending: KappaM_out is a fresh
 *    overwrite (:704), KappaH_out likewise (:720). The driver locals
 *    tke_Av_old/tke_Kv_old (gen:433-434) fed only these dead args.
 *  - handle_old_vals: consumed only by cvmix_update_tke, which is called only
 *    from the commented-out block in the never-called tke_wrap (:399-409).
 *  - max_nlev: unused in the body (only the commented Kappa_GM decl).
 *  - i / j / tstep_count: read only inside the `if (.false.)` debug block
 *    (:916-986).
 *  - tke_userdef_constants: absent at the FESOM call site → the module saved
 *    struct is used; here one static struct set by fesom_cvmix_init_tke.
 *  - tke_wrap (:285-411), put_tke machinery (:991-1093): never called /
 *    collapsed into struct assignment.
 * `bottom_fric` is KEPT for call-site parity although the executed body never
 * reads it (the Fortran driver passes tke_forc2d_botfrict(node) ≡ 0).
 *
 * == Gate-only branches (reference config: only_tke=T, l_lc=F, Dirichlet
 *    BCs F/F — fesom_tke_alloc aborts if a port-excluded option is set) ==
 *  - IDEMIX coupling (.not.only_tke): Rinum modification (:710-712) and
 *    iw_diss forcing (:742-744) — bodies not ported; E_iw/alpha_c are read
 *    only under the gate (pass the shared zero column).
 *  - Langmuir (l_lc): forc += tke_plc (:737-739) — gated, zero column.
 *  - Dirichlet upper/lower BC branches (:780-790, :799-811, :841-849,
 *    :881-883, :890-892): NOT ported (abort at init if enabled); the executed
 *    Neumann branches are.
 *  - iw_diss is READ UNCONDITIONALLY at :898 (tke_Tiwf = iw_diss) → the
 *    caller must pass a REAL zero array of length nlev+1, never NULL.
 */

/* tke_type (:67-93) — set once by fesom_cvmix_init_tke, read by every
 * integrate call (the `tke_constants_saved` module variable). */
typedef struct fesom_cvmix_tke_params {
    real_t c_k;            /* {0.1}   KappaM = c_k·L·√e                  */
    real_t c_eps;          /* {0.7}   dissipation coefficient            */
    real_t cd;             /* {3.75}  surface-TKE coefficient — NAMELIST
                              value; the executed NEUMANN BC uses it
                              (:793), the "3.75=Dirichlet" comment in the
                              Fortran is advisory only                    */
    real_t alpha_tke;      /* {30.0}  TKE diffusivity multiplier + L bound */
    real_t clc;            /* {0.3}   Langmuir factor (gate-only)         */
    real_t mxl_min;        /* {1e-8}  mixing-length floor                 */
    real_t kappaM_min;     /* {0.0}   (read into a local, never applied —
                              the commented :663 clamp; kept for parity)  */
    real_t kappaM_max;     /* {100.0} viscosity clamp                     */
    real_t tke_surf_min;   /* {1e-4}  min surface TKE                     */
    real_t tke_min;        /* {1e-6}  min interior TKE (only_tke floor)   */
    int    tke_mxl_choice; /* {2}     Blanke–Delecluse (only option)      */
    int    only_tke;       /* {1}                                          */
    int    l_lc;           /* {0}     Langmuir gate                       */
    int    use_ubound_dirichlet;  /* {0} */
    int    use_lbound_dirichlet;  /* {0} */
} fesom_cvmix_tke_params;

/* init_tke (:101-281): range-check every parameter exactly as the Fortran
 * (stop 1 on violation) and store into the saved struct. The FESOM call site
 * passes every optional except handle_old_vals/tke_userdef_constants, so the
 * "absent → default" branches are dead — values come from namelist.cvmix
 * &param_tke via gen:258-272. */
void fesom_cvmix_init_tke(real_t c_k, real_t c_eps, real_t cd,
                          real_t alpha_tke, real_t mxl_min,
                          real_t kappaM_min, real_t kappaM_max,
                          int    tke_mxl_choice,
                          int    use_ubound_dirichlet,
                          int    use_lbound_dirichlet,
                          int    only_tke, int l_lc, real_t clc,
                          real_t tke_min, real_t tke_surf_min);

/* Read-only access to the saved constants (for gate checks / dumps). */
const fesom_cvmix_tke_params *fesom_cvmix_tke_constants(void);

/* Budget-diagnostic copy-out targets — each pointer is a [nlev+1] column
 * slice or NULL. integrate_tke ALWAYS computes the full set into local
 * column scratch in Fortran order (internal read-after-write deps:
 * P_diss_v/K_diss_v seed Tbpr/Tspr; diag writes NEVER go through these
 * pointers mid-computation) and copies each non-NULL target at the END of
 * the column only. Pass NULL (or a NULL struct pointer) to discard. */
typedef struct fesom_cvmix_tke_diag {
    real_t *Tbpr;   /* -P_diss_v          buoyancy production  */
    real_t *Tspr;   /*  K_diss_v          shear production     */
    real_t *Tdif;   /* tridiag residual   vertical diffusion   */
    real_t *Tdis;   /* -c_eps/L·√e·e_new  dissipation          */
    real_t *Twin;   /* surface wind forcing (Neumann: cd·f^1.5/dzt1) */
    real_t *Tiwf;   /* = iw_diss (≡0 under only_tke)           */
    real_t *Tbck;   /* (tke_new-tke_unrest)/dt  tke_min reset  */
    real_t *Ttot;   /* (tke_new-tke_old)/dt                    */
    real_t *Lmix;   /* mixing length mxl                       */
    real_t *Pr;     /* Prandtl number                          */
    real_t *int1;   /* cvmix_int_1 = KappaH_out                */
    real_t *int2;   /* cvmix_int_2 = KappaM_out                */
    real_t *int3;   /* cvmix_int_3 = Nsqr                      */
} fesom_cvmix_tke_diag;

/*
 * integrate_tke (:415-987). All column arrays are interface-indexed
 * [nlev+1] except dzw [nlev] (layer thicknesses, = hnode slice); dzt is the
 * tracer-point distance with half-thickness end values. 0-based k here ↔
 * Fortran k-1. Executed path: tke_mxl_choice=2, Neumann BCs both ends,
 * only_tke (tke_min floor applied), l_lc off.
 *
 * Outputs: tke_new (the solved+floored TKE), KappaM_out (min(KappaM_max,
 * c_k·mxl·sqrttke) — FRESH overwrite), KappaH_out (= KappaM_out/prandtl).
 */
void fesom_cvmix_integrate_tke(int           nlev,
                               real_t        dtime,
                               real_t        rho_ref,
                               real_t        grav,
                               const real_t *dzw,        /* [nlev]   */
                               const real_t *dzt,        /* [nlev+1] */
                               const real_t *tke_old,    /* [nlev+1] */
                               const real_t *Ssqr,       /* [nlev+1] */
                               const real_t *Nsqr,       /* [nlev+1] */
                               const real_t *alpha_c,    /* [nlev+1] gated (.not.only_tke) */
                               const real_t *E_iw,       /* [nlev+1] gated (.not.only_tke) */
                               const real_t *iw_diss,    /* [nlev+1] READ UNCONDITIONALLY  */
                               const real_t *tke_plc,    /* [nlev+1] gated (l_lc)          */
                               real_t        forc_tke_surf,
                               real_t        forc_rho_surf,
                               real_t        bottom_fric, /* unread in executed body */
                               real_t       *tke_new,    /* out [nlev+1] */
                               real_t       *KappaM_out, /* out [nlev+1] */
                               real_t       *KappaH_out, /* out [nlev+1] */
                               const fesom_cvmix_tke_diag *diag); /* nullable */

#endif /* FESOM_CVMIX_TKE_H */
