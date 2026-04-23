/*
 * Phase 1 timestep driver — wires the 11 substeps in FRESH_START.md §5 order.
 * Each call performs one full timestep.
 */
#include "fesom_step.h"
#include "fesom_ale.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_eos.h"
#include "fesom_forcing.h"
#include "fesom_mesh.h"
#include "fesom_momentum.h"
#include "fesom_pp.h"
#include "fesom_ssh.h"
#include "fesom_tracer_adv.h"
#include "fesom_tracers.h"

int fesom_timestep(int                          step_n,
                   const fesom_step_ctx        *ctx,
                   struct fesom_mesh           *mesh,
                   struct fesom_aux            *aux,
                   struct fesom_dyn            *dyn,
                   struct fesom_tracers        *tracers,
                   const struct fesom_forcing  *forcing)
{
    /*  1. EOS + hydrostatic pressure + N²  */
    fesom_pressure_bv(tracers, mesh, aux);

    /*  2. PGF (linfs + full cells)  */
    fesom_pressure_force_linfs_fullcell(mesh, aux);

    /*  3. mixing: UVnode → PP → convective adjustment  */
    fesom_compute_vel_nodes(mesh, dyn);
    fesom_pp_mixing(mesh, dyn, aux);
    fesom_mo_convect(mesh, aux);

    /*  4. momentum RHS (Coriolis AB2 + SSH gradient + PGF)  */
    fesom_compute_vel_rhs(mesh, aux, dyn, /*is_first_step=*/(step_n == 1));

    /*  5. horizontal viscosity (opt_visc=5 with backscatter) — FRESH_START.md §5  */
    fesom_visc_filt_bcksct(mesh, dyn);

    /*  6. implicit vertical viscosity TDMA  */
    fesom_impl_vert_visc(mesh, aux, forcing, dyn);

    /*  6. SSH RHS (linfs)  */
    fesom_compute_ssh_rhs_linfs(mesh, dyn);

    /*  7. CG SSH solve  */
    int cg_iters = fesom_ssh_solve_cg(ctx->stiff, ctx->solver, mesh, dyn);

    /*  8. velocity update with SSH-gradient correction  */
    fesom_update_vel(mesh, dyn);

    /*  9. transport-divergence → ssh_rhs_old, then hbar update  */
    fesom_compute_hbar(mesh, dyn);

    /* 10. eta_n inline (oce_ale.F90:3771-3775).
          eta_n = α·hbar + (1-α)·hbar_old, but ONLY where ulevels_nod2D == 1.  */
    {
        const real_t alpha = (real_t)FESOM_PHASE1_ALPHA;
        for (int n = 0; n < mesh->nod2D; ++n) {
            if (mesh->ulevels_nod2D[n] == 1) {
                dyn->eta_n[n] = alpha * mesh->hbar[n]
                              + (1.0 - alpha) * mesh->hbar_old[n];
            }
        }
    }

    /* 11. ALE step (linfs):
         hnode_new := hnode (no thickness motion)
         w        := vert. divergence integration
         CFLz + wsplit: split w → w_e (explicit, for tracers) and w_i (implicit) */
    fesom_ale_thickness_linfs(mesh);
    fesom_ale_vert_vel_linfs(mesh, dyn);
    fesom_ale_compute_cflz(mesh, dyn);
    fesom_ale_compute_wvel_split(mesh, dyn);

    /* 12. tracer advection: T then S, FCT (central HO + QR4C vert + Zalesak). */
    fesom_tracer_advect_one_fct(ctx->tra_sc, FESOM_TRACER_T, mesh, dyn, tracers);
    fesom_tracer_advect_one_fct(ctx->tra_sc, FESOM_TRACER_S, mesh, dyn, tracers);

    /* 13. commit thickness: hnode := hnode_new, helem from vertex mean.  */
    fesom_ale_commit_thickness(mesh);

    return cg_iters;
}
