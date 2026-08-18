#ifndef FESOM_STEP_H
#define FESOM_STEP_H

#include "fesom_types.h"
#include "fesom_ssh.h"
#include "fesom_tracer_adv.h"

struct fesom_mesh;
struct fesom_aux;
struct fesom_dyn;
struct fesom_tracers;
struct fesom_forcing;

/*
 * Phase 1 timestep driver. Mirrors the sequence in oce_timestep_ale
 * (oce_ale.F90:3339-3943) restricted to the Phase 1 subset:
 *
 *  1. compute_pressure_bv     (EOS, density, hpressure, bvfreq)
 *  2. compute_pressure_force  (PGF — linfs + full cells)
 *  3. mixing                  (compute_vel_nodes → pp_mixing → mo_convect)
 *  4. compute_vel_rhs         (Coriolis AB2 + SSH grad + PGF)
 *  5. impl_vert_visc          (TDMA: implicit visc + wind stress + bottom drag)
 *  6. compute_ssh_rhs         (linfs branch)
 *  7. ssh_solve (CG)          (preconditioned conjugate gradient)
 *  8. update_vel              (uv += uv_rhs - g·θ·dt·∇d_eta)
 *  9. compute_hbar            (transport divergence → ssh_rhs_old, hbar)
 * 10. eta_n update inline     (eta_n = α·hbar + (1-α)·hbar_old)
 * 11. ALE step (linfs)        (hnode_new = hnode; vert_vel from divergence)
 * 12. solve_tracers           (upwind T, then upwind S; ALE reconstruction)
 * 13. commit thickness        (hnode = hnode_new; helem = mean of vertex hnode)
 *
 * Deferred (off in Phase 1, gated by namelist in Fortran):
 *   - sw_alpha_beta (only needed by KPP/GM)
 *   - sigma_xy / neutral_slope (GM/Redi)
 *   - momentum_advection (momadv_opt=2)
 *   - viscosity_filter (no horizontal viscosity)
 *   - GM/Redi bolus velocity update
 *   - update_stiff_mat_ale (linfs branch skips it)
 */
struct fesom_partit;
struct fesom_ice;
struct fesom_jra55;
struct fesom_sss_runoff;
struct fesom_gm;
struct fesom_kpp;
struct fesom_tke;

typedef struct fesom_step_ctx {
    fesom_ssh_stiff           *stiff;
    fesom_solverinfo          *solver;
    fesom_tracer_adv_scratch  *tra_sc;
    struct fesom_partit       *partit;     /* needed for halo exchanges */
    struct fesom_ice          *ice;        /* sea-ice state; NULL = no ice */
    struct fesom_gm           *gm;         /* GM/Redi state; NULL = stub-only path */
    struct fesom_kpp          *kpp;        /* KPP mixing state; used when FESOM_MIX_SCHEME=KPP */
    struct fesom_tke          *tke;        /* CVMix-TKE mixing state; non-NULL only when
                                              FESOM_MIX_SCHEME=TKE|cvmix_TKE (allocated in main) */
    /* Optional pointers for the sea-ice thermodynamics path. NULL is fine
     * if the FESOM_NO_ICE_THERMO env knob is set. */
    const struct fesom_jra55      *jra;
    const struct fesom_sss_runoff *sr;

    /* which_ALE switch (Phase Z0) — copies of the fesom_ale.h globals
     * (fesom_ale_zstar / fesom_use_virt_salt / fesom_is_nonlinfs), set by
     * main after fesom_ale_mode_init. linfs defaults keep every consumer on
     * the v1.0 code path. */
    /* Set when the run resumed from a restart. Mirrors Fortran's
     * `if (lfirst .and. (.not. r_restart))` guards: step 1 of a resumed run is
     * not the first step of the integration, so the Adams-Bashforth
     * cold-start treatment must not fire again. Without this the resumed leg
     * uses ff_step = 1 where the uninterrupted run uses ab2. */
    int    restarted;

    int    ale_zstar;       /* 0 = linfs (default), 1 = zstar */
    int    use_virt_salt;   /* 1 under linfs, 0 under zstar   */
    real_t is_nonlinfs;     /* 0.0 under linfs, 1.0 under zstar */
} fesom_step_ctx;

/* Bitwise state dump used by the restart round-trip diagnosis; see
 * restart_probe in fesom_step.c. Gated by FESOM_RESTART_PROBE=<step> for the
 * in-step call and FESOM_RESTART_PROBE_TOP=<step> for the call the driver makes
 * at the top of the iteration, before the forcing and the ice step. */
struct fesom_aux;
struct fesom_forcing;
void fesom_restart_probe_state(const char *gate_env, int step_n,
                               const struct fesom_step_ctx *ctx,
                               const struct fesom_mesh *mesh,
                               struct fesom_dyn *dyn,
                               struct fesom_tracers *tracers,
                               const struct fesom_forcing *forcing,
                               const struct fesom_aux *aux);


/*
 * step_n: the current step index, 1-based to mirror Fortran. The very first
 * call must use step_n = 1 — compute_vel_rhs uses lfirst-equivalent gating.
 *
 * Returns the CG iteration count for diagnostics.
 */
int fesom_timestep(int                          step_n,
                   const fesom_step_ctx        *ctx,
                   struct fesom_mesh           *mesh,
                   struct fesom_aux            *aux,
                   struct fesom_dyn            *dyn,
                   struct fesom_tracers        *tracers,
                   const struct fesom_forcing  *forcing);

#endif /* FESOM_STEP_H */
