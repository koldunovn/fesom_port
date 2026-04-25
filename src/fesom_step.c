/*
 * Phase 1 timestep driver — wires the substeps in FRESH_START.md §5 order.
 * For Phase 4 (MPI) every kernel that writes a field other ranks read is
 * followed by a halo exchange. Cheat sheet documented inline.
 */
#include "fesom_step.h"
#include "fesom_ale.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_eos.h"
#include "fesom_forcing.h"
#include "fesom_gm.h"
#include "fesom_halo.h"
#include "fesom_ice.h"
#include "fesom_mesh.h"
#include "fesom_momentum.h"
#include "fesom_partit.h"
#include "fesom_pp.h"
#include "fesom_ssh.h"
#include "fesom_tracer_adv.h"
#include "fesom_tracer_diff.h"
#include "fesom_tracers.h"

#include <stdlib.h>

int fesom_timestep(int                          step_n,
                   const fesom_step_ctx        *ctx,
                   struct fesom_mesh           *mesh,
                   struct fesom_aux            *aux,
                   struct fesom_dyn            *dyn,
                   struct fesom_tracers        *tracers,
                   const struct fesom_forcing  *forcing)
{
    fesom_partit *p = ctx->partit;
    int nl = mesh->nl;

    /*  1. EOS + hydrostatic pressure + N²  */
    fesom_pressure_bv(tracers, mesh, aux);
    /* exchange the per-node 3D outputs (Fortran oce_ale_pressure_bv.F90:2844-) */
    fesom_exchange_nod3D(aux->density_m_rho0, nl, p);
    fesom_exchange_nod3D(aux->hpressure,      nl, p);
    fesom_exchange_nod3D(aux->bvfreq,         nl, p);
    fesom_exchange_nod3D(aux->sw_alpha,       nl, p);
    fesom_exchange_nod3D(aux->sw_beta,        nl, p);

    /*  1b. GM/Redi prerequisites + per-step coefficient builder
     *      (Phases G2b + G3). Outputs fer_K, fer_C, Ki, fer_tapfac,
     *      neutral_slope, slope_tapered. No readers yet (G4 / G5 / G7
     *      will plumb them). G8 will gate this block under gmredi_on. */
    if (ctx->gm) {
        fesom_compute_sigma_xy     (aux, tracers, mesh, ctx->gm, p);
        fesom_compute_neutral_slope(aux,           mesh, ctx->gm, p);
        fesom_init_redi_gm         (aux,           mesh, ctx->gm, p);
    }

    /*  2. PGF (linfs + full cells)  */
    fesom_pressure_force_linfs_fullcell(mesh, aux);
    fesom_exchange_elem3D(aux->pgf_x, nl, p);
    fesom_exchange_elem3D(aux->pgf_y, nl, p);

    /*  3. mixing: UVnode → PP → convective adjustment  */
    fesom_compute_vel_nodes(mesh, dyn);
    fesom_halo_exchange(dyn->uvnode, FESOM_HALO_NOD3D, nl, 2, p);

    fesom_pp_mixing(mesh, dyn, aux);
    fesom_exchange_nod3D (aux->Kv, nl, p);
    fesom_exchange_elem3D(aux->Av, nl, p);

    fesom_mo_convect(mesh, aux);
    fesom_exchange_nod3D (aux->Kv, nl, p);
    fesom_exchange_elem3D(aux->Av, nl, p);

    /*  4. momentum RHS (Coriolis AB2 + SSH gradient + PGF)  */
    fesom_compute_vel_rhs(mesh, aux, dyn, /*is_first_step=*/(step_n == 1));
    /* uv_rhs needed by visc_filt_bcksct on halo elements */
    fesom_halo_exchange(dyn->uv_rhs,   FESOM_HALO_ELEM3D, nl, 2, p);
    fesom_halo_exchange(dyn->uv_rhsAB, FESOM_HALO_ELEM3D, nl, 2, p);

    /*  5. horizontal viscosity (opt_visc=5 with backscatter)  */
    fesom_visc_filt_bcksct(mesh, dyn, p);
    /* uv_rhs is the final output; needed at halo for impl_vert_visc neighbour
     * reads through the TDMA SpMV. */
    fesom_halo_exchange(dyn->uv_rhs, FESOM_HALO_ELEM3D, nl, 2, p);

    /*  6. implicit vertical viscosity TDMA  */
    fesom_impl_vert_visc(mesh, aux, forcing, dyn);
    fesom_halo_exchange(dyn->uv_rhs, FESOM_HALO_ELEM3D, nl, 2, p);

    /*  7. SSH RHS (linfs)  */
    fesom_compute_ssh_rhs_linfs(mesh, dyn);
    fesom_exchange_nod2D(dyn->ssh_rhs, p);   /* Fortran oce_ale.F90:1954 */

    /*  8. CG SSH solve  */
    int cg_iters = fesom_ssh_solve_cg(ctx->stiff, ctx->solver, mesh, dyn);
    /* fesom_ssh_solve_cg already exchanges its own X (= d_eta) at exit
     * (Fortran solver.F90:279). For now we keep the explicit exchange here
     * because the parallel CG modifications are slice 30e, not 30d. */
    fesom_exchange_nod2D(dyn->d_eta, p);

    /*  9. velocity update with SSH-gradient correction  */
    fesom_update_vel(mesh, dyn);
    fesom_halo_exchange(dyn->uv, FESOM_HALO_ELEM3D, nl, 2, p);

    /* 10. transport-divergence → ssh_rhs_old, then hbar update  */
    fesom_compute_hbar(mesh, dyn);
    fesom_exchange_nod2D(dyn->ssh_rhs_old, p);   /* Fortran oce_ale.F90:2078 */
    fesom_exchange_nod2D(mesh->hbar,       p);   /* Fortran oce_ale.F90:2102 */

    /* 11. eta_n inline (oce_ale.F90:3771-3775).
          eta_n = α·hbar + (1-α)·hbar_old, but ONLY where ulevels_nod2D == 1.  */
    {
        const real_t alpha = (real_t)FESOM_PHASE1_ALPHA;
        int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
        for (int n = 0; n < N; ++n) {
            if (mesh->ulevels_nod2D[n] == 1) {
                dyn->eta_n[n] = alpha * mesh->hbar[n]
                              + (1.0 - alpha) * mesh->hbar_old[n];
            }
        }
    }
    /* eta_n already covers myDim+eDim because hbar/hbar_old are exchanged. */

    /* 12. ALE step (linfs). */
    fesom_ale_thickness_linfs(mesh);
    /* hnode_new = hnode (no exchange needed; both already cover halo). */

    fesom_ale_vert_vel_linfs(mesh, dyn);
    fesom_exchange_nod3D(dyn->w, nl, p);     /* Fortran oce_ale.F90:2679 */

    fesom_ale_compute_cflz(mesh, dyn);
    fesom_exchange_nod3D(dyn->cfl_z, nl, p);

    fesom_ale_compute_wvel_split(mesh, dyn);
    fesom_exchange_nod3D(dyn->w_e, nl, p);
    fesom_exchange_nod3D(dyn->w_i, nl, p);

    /* 13. tracer advection: T then S.
     * FESOM_NO_TRADV=1 → skip advection (diagnostic for MPI drift hunt). */
    static int nt_adv_checked = 0, nt_adv_skip = 0;
    if (!nt_adv_checked) {
        const char *e = getenv("FESOM_NO_TRADV");
        nt_adv_skip = (e && atoi(e));
        nt_adv_checked = 1;
    }
    if (!nt_adv_skip) {
        fesom_tracer_advect_one_fct(ctx->tra_sc, FESOM_TRACER_T, mesh, dyn, tracers);
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_T].values, nl, p);

        fesom_tracer_advect_one_fct(ctx->tra_sc, FESOM_TRACER_S, mesh, dyn, tracers);
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_S].values, nl, p);
    }

    /* 13b. implicit vertical diffusion of tracers + surface heat/water flux BC.
     * FESOM_NO_TRDIFF=1 → skip diffusion. */
    static int nt_dif_checked = 0, nt_dif_skip = 0;
    if (!nt_dif_checked) {
        const char *e = getenv("FESOM_NO_TRDIFF");
        nt_dif_skip = (e && atoi(e));
        nt_dif_checked = 1;
    }
    if (!nt_dif_skip) {
        fesom_impl_vert_diff_tracers(mesh, aux, forcing, tracers);
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_T].values, nl, p);
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_S].values, nl, p);
    }

    /* 14. commit thickness: hnode := hnode_new, helem from vertex mean. */
    fesom_ale_commit_thickness(mesh);
    fesom_exchange_nod3D (mesh->hnode, nl, p);   /* both already same — but be explicit */
    fesom_exchange_elem3D(mesh->helem, nl, p);   /* Fortran oce_ale.F90:1027,1249 */

    /* Sea-ice step is now called from fesom_main BEFORE the ocean step
     * (ice writes heat_flux/water_flux that the ocean step consumes). */

    return cg_iters;
}
