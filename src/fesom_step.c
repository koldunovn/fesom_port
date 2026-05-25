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
#include "fesom_kpp.h"
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

    /* Phase G8 — env-gated master switch for GM/Redi.
     *   FESOM_NO_GMREDI=1 → treat ctx->gm as NULL for the rest of this step.
     *                      Skips compute_sigma_xy, compute_neutral_slope,
     *                      init_Redi_GM, fer_solve_Gamma, fer_gamma2vel,
     *                      Redi K33 in tracer diff, bolus add/sub, and
     *                      G7a/G7b. Used by G9 bit-identity off-switch test.
     * Read once and cache (matches FESOM_NO_TRADV style elsewhere). */
    static int s_gmredi_env_loaded = 0;
    static int s_no_gmredi         = 0;
    if (!s_gmredi_env_loaded) {
        const char *e = getenv("FESOM_NO_GMREDI");
        s_no_gmredi = (e && atoi(e));
        s_gmredi_env_loaded = 1;
    }
    /* Local alias so the original ctx->gm isn't mutated. */
    struct fesom_gm *gm = s_no_gmredi ? NULL : ctx->gm;

    /* Vertical-mixing scheme dispatch (mirror oce_ale.F90:3515 mix_scheme_nmb).
     *   FESOM_MIX_SCHEME=KPP (DEFAULT) → fesom_kpp_mixing (K-Profile) — the CORE2
     *                                    production scheme (mix_scheme='KPP'),
     *                                    validated end-to-end (K0-K10).
     *   FESOM_MIX_SCHEME=PP            → fesom_pp_mixing (Pacanowski-Philander), opt-out
     * Env knob mirroring the FESOM_NO_GMREDI pattern. mo_convect runs after either
     * scheme (Fortran calls it in both branches, :3524 / :3531). */
    static int s_mix_env_loaded = 0;
    static int s_use_kpp        = 1;
    if (!s_mix_env_loaded) {
        const char *e = getenv("FESOM_MIX_SCHEME");
        s_use_kpp = !(e && (e[0] == 'P' || e[0] == 'p'));   /* default KPP; =PP to opt out */
        s_mix_env_loaded = 1;
    }

    /*  1. EOS + hydrostatic pressure + N²  */
    fesom_pressure_bv(tracers, mesh, aux);
    /* sw_alpha / sw_beta — McDougall (1987). Needed by GM/Redi (and KPP).
     * Mirror of Fortran oce_ale.F90:3475 sw_alpha_beta. */
    fesom_compute_sw_alpha_beta(tracers, mesh, aux);
    /* exchange the per-node 3D outputs (Fortran oce_ale_pressure_bv.F90:2844-) */
    fesom_exchange_nod3D(aux->density_m_rho0, nl, p);
    fesom_exchange_nod3D(aux->hpressure,      nl, p);
    fesom_exchange_nod3D(aux->bvfreq,         nl, p);
    fesom_exchange_nod3D(aux->sw_alpha,       nl, p);
    fesom_exchange_nod3D(aux->sw_beta,        nl, p);
    /* horizontal N² smoothing — N2smth_h=.true., N2smth_hidx=1 (Fortran
     * oce_ale_pressure_bv.F90:499 smooth_nod3D(bvfreq,1)). The exchange above
     * populated the halo bvfreq the sweep reads (Fortran fills bvfreq on
     * myDim+eDim then smooths); fesom_smooth_nod3D re-exchanges after. Must precede
     * every bvfreq consumer (GM, PP/KPP, mo_convect). */
    fesom_smooth_nod3D(aux->bvfreq, nl, 1, mesh, p);

    /*  1b. GM/Redi prerequisites + per-step coefficient builder +
     *      streamfunction solve + bolus velocity reconstruction
     *      (Phases G2b + G3 + G4). Outputs sigma_xy / neutral_slope /
     *      slope_tapered / fer_tapfac / fer_K / fer_C / Ki / fer_gamma /
     *      dyn->fer_uv. Still no readers (G5 / G6 / G7 will plumb the
     *      tracer-side use). G8 will gate this block under gmredi_on. */
    if (gm) {
        fesom_compute_sigma_xy     (aux, tracers, mesh, gm, p);
        fesom_compute_neutral_slope(aux,           mesh, gm, p);
        fesom_init_redi_gm         (aux,           mesh, gm, p);
        fesom_fer_solve_gamma      (aux,           mesh, gm, p);
        fesom_fer_gamma2vel        (dyn,           mesh, gm, p);
    }

    /*  2. PGF (linfs + full cells)  */
    fesom_pressure_force_linfs_fullcell(mesh, aux);
    fesom_exchange_elem3D(aux->pgf_x, nl, p);
    fesom_exchange_elem3D(aux->pgf_y, nl, p);

    /*  3. mixing: UVnode → (PP or KPP) → convective adjustment.
     *     Dispatch mirrors oce_ale.F90:3515-3531 (mix_scheme_nmb 1=KPP / 2=PP);
     *     mo_convect is shared (called after either scheme).  */
    fesom_compute_vel_nodes(mesh, dyn);
    fesom_halo_exchange(dyn->uvnode, FESOM_HALO_NOD3D, nl, 2, p);

    if (s_use_kpp) {
        /* KPP writes aux->Av (elements) + the single aux->Kv (T-channel,
         * oce_ale.F90:3518-3522) and does its own internal halo exchanges. */
        fesom_kpp_mixing(ctx->kpp, aux, tracers, forcing, dyn, mesh, p);
    } else {
        fesom_pp_mixing(mesh, dyn, aux);
        fesom_exchange_nod3D (aux->Kv, nl, p);
        fesom_exchange_elem3D(aux->Av, nl, p);
    }

    fesom_mo_convect(mesh, aux);
    fesom_exchange_nod3D (aux->Kv, nl, p);
    fesom_exchange_elem3D(aux->Av, nl, p);

    /*  4. momentum RHS (Coriolis AB2 + SSH gradient + PGF)  */
    fesom_compute_vel_rhs(mesh, aux, dyn, /*is_first_step=*/(step_n == 1), p);
    /* uv_rhs needed by visc_filt_bcksct on halo elements */
    fesom_halo_exchange(dyn->uv_rhs,   FESOM_HALO_ELEM3D, nl, 2, p);
    fesom_halo_exchange(dyn->uv_rhsAB, FESOM_HALO_ELEM3D, nl, 2, p);

    /*  5. horizontal viscosity. CORE2 runs opt_visc=7 (biharmonic, flow-aware);
     *     the biharmonic ∇⁴ damps grid-scale 2Δx modes that opt_visc=5 lets grow
     *     — required for dt=1800 stability (work_core/namelist.dyn opt_visc=7).  */
    fesom_visc_filt_bidiff(mesh, dyn, p);
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

    /* DEBUG (FESOM_DIAG_SSHSLV=<gid>): dump ssh_rhs / d_eta / hbar at one node,
       matching the Fortran [FSSH] format, to compare the implicit free-surface
       DAMPING STRENGTH (|d_eta|/|ssh_rhs|) C-vs-Fortran. A systematically
       larger C d_eta = the stiffness operator under-damps the 2dx mode. */
    if (getenv("FESOM_DIAG_SSHSLV")) {
        int tgid = atoi(getenv("FESOM_DIAG_SSHSLV"));
        for (int row = 0; row < mesh->myDim_nod2D; ++row) {
            if (p->myList_nod2D[row] != tgid) continue;
            fprintf(stderr, "[sshslv r%d] gid %d step %d ssh_rhs=%.7e d_eta=%.7e hbar=%.7e\n",
                    p->mype, tgid, step_n, (double)dyn->ssh_rhs[row],
                    (double)dyn->d_eta[row], (double)mesh->hbar[row]);
            break;
        }
        fflush(stderr);
    }

    /* DEBUG (FESOM_DIAG_SPREAD=<gid>): per-step 1-ring spread of the fields that
       could seed the dt=1800 central-Arctic 2dx blow-up, so we can see WHICH one
       erupts first (SSH-solve side vs element-velocity side) and WHEN. For the
       target node: hbar/d_eta/ssh_rhs spread over the 1-ring (node + edge
       neighbours via surrounding elements), plus the surface |uv| min/max over
       the surrounding ELEMENTS (a checkerboard in element velocity = the 2dx). */
    if (getenv("FESOM_DIAG_SPREAD")) {
        int tgid = atoi(getenv("FESOM_DIAG_SPREAD"));
        for (int row = 0; row < mesh->myDim_nod2D; ++row) {
            if (p->myList_nod2D[row] != tgid) continue;
            double hmin = mesh->hbar[row],  hmax = hmin;
            double emin = dyn->d_eta[row],  emax = emin;
            double rmin = dyn->ssh_rhs[row], rmax = rmin;
            double vmin = 1e30, vmax = -1e30;   /* surface |uv| over surrounding elems */
            double pmin = 1e30, pmax = -1e30;   /* surface |PGF| over surrounding elems */
            double qmin = 1e30, qmax = -1e30;   /* surface |uv_rhs| over surrounding elems */
            int o0 = mesh->nod_in_elem2D_offsets[row];
            int o1 = mesh->nod_in_elem2D_offsets[row + 1];
            for (int k = o0; k < o1; ++k) {
                int e = mesh->nod_in_elem2D[k];
                double u = dyn->uv[FESOM_ELEMVEC(e, 0, nl) + 0];
                double v = dyn->uv[FESOM_ELEMVEC(e, 0, nl) + 1];
                double spd = sqrt(u*u + v*v);
                if (spd < vmin) vmin = spd;
                if (spd > vmax) vmax = spd;
                double px = aux->pgf_x[FESOM_ELEM3D(e, 0, nl)];
                double py = aux->pgf_y[FESOM_ELEM3D(e, 0, nl)];
                double pmag = sqrt(px*px + py*py);
                if (pmag < pmin) pmin = pmag;
                if (pmag > pmax) pmax = pmag;
                double qx = dyn->uv_rhs[FESOM_ELEMVEC(e, 0, nl) + 0];
                double qy = dyn->uv_rhs[FESOM_ELEMVEC(e, 0, nl) + 1];
                double qmag = sqrt(qx*qx + qy*qy);
                if (qmag < qmin) qmin = qmag;
                if (qmag > qmax) qmax = qmag;
                for (int j = 0; j < 3; ++j) {
                    int nd = mesh->elem_nodes[3*e + j];
                    if (mesh->hbar[nd]  < hmin) hmin = mesh->hbar[nd];
                    if (mesh->hbar[nd]  > hmax) hmax = mesh->hbar[nd];
                    if (dyn->d_eta[nd]  < emin) emin = dyn->d_eta[nd];
                    if (dyn->d_eta[nd]  > emax) emax = dyn->d_eta[nd];
                    if (dyn->ssh_rhs[nd] < rmin) rmin = dyn->ssh_rhs[nd];
                    if (dyn->ssh_rhs[nd] > rmax) rmax = dyn->ssh_rhs[nd];
                }
            }
            fprintf(stderr, "[spread r%d] gid %d step %d hbar_sprd=%.4e deta_sprd=%.4e "
                    "sshrhs_sprd=%.4e velsurf[%.3e,%.3e] velsprd=%.3e pgf_sprd=%.4e "
                    "uvrhs_sprd=%.4e\n",
                    p->mype, tgid, step_n, hmax - hmin, emax - emin,
                    rmax - rmin, vmin, vmax, vmax - vmin, pmax - pmin, qmax - qmin);
            break;
        }
        fflush(stderr);
    }

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

    fesom_ale_vert_vel_linfs(mesh, dyn, gm ? 1 : 0);
    fesom_exchange_nod3D(dyn->w, nl, p);     /* Fortran oce_ale.F90:2679 */
    if (gm) {
        /* Mirror Fortran oce_ale.F90:2681 — exchange_nod(fer_Wvel). */
        fesom_exchange_nod3D(dyn->fer_w, nl, p);
    }

    fesom_ale_compute_cflz(mesh, dyn);
    fesom_exchange_nod3D(dyn->cfl_z, nl, p);

    fesom_ale_compute_wvel_split(mesh, dyn);
    fesom_exchange_nod3D(dyn->w_e, nl, p);
    fesom_exchange_nod3D(dyn->w_i, nl, p);

    /* 13a. Phase G6b — bolus velocity add (Fortran oce_ale_tracer.F90:199-211).
     * Wraps BOTH advection (13) and diffusion (13b) so each tracer sees the
     * bolus-augmented velocity field. Subtracted back after diffusion. */
    if (gm) {
        const int E_loop = mesh->myDim_elem2D + mesh->eDim_elem2D;
        const int N_loop = mesh->myDim_nod2D + mesh->eDim_nod2D;
        for (int e = 0; e < E_loop; ++e) {
            for (int nz = 0; nz < nl; ++nz) {
                size_t k = (size_t)e * nl * 2 + nz * 2;
                dyn->uv[k + 0] += dyn->fer_uv[k + 0];
                dyn->uv[k + 1] += dyn->fer_uv[k + 1];
            }
        }
        for (int n = 0; n < N_loop; ++n) {
            for (int nz = 0; nz < nl; ++nz) {
                size_t k = (size_t)n * nl + nz;
                dyn->w  [k] += dyn->fer_w[k];
                dyn->w_e[k] += dyn->fer_w[k];
            }
        }
    }

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
        if (gm) {
            /* G7a builds gm->tr_xy + adds vertical-explicit Redi flux to T. */
            fesom_diff_ver_part_redi_expl(FESOM_TRACER_T, gm, aux, mesh, tracers, p);
            /* G7b reuses gm->tr_xy + builds gm->tr_z + horizontal Redi flux. */
            fesom_diff_part_hor_redi    (FESOM_TRACER_T, gm, aux, mesh, tracers, p);
        }
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_T].values, nl, p);

        fesom_tracer_advect_one_fct(ctx->tra_sc, FESOM_TRACER_S, mesh, dyn, tracers);
        if (gm) {
            fesom_diff_ver_part_redi_expl(FESOM_TRACER_S, gm, aux, mesh, tracers, p);
            fesom_diff_part_hor_redi    (FESOM_TRACER_S, gm, aux, mesh, tracers, p);
        }
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
        fesom_impl_vert_diff_tracers(mesh, aux, forcing, tracers, gm);
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_T].values, nl, p);
        fesom_exchange_nod3D(tracers->data[FESOM_TRACER_S].values, nl, p);
    }

    /* Salinity floor — port of the *intent* behind Fortran's tracer_cutoff
     * namelist (declared in oce_modules.F90:172-176 but never actually
     * enforced in the Fortran source either; declared dead-code).
     *
     * Why we need it: confined brackish seas (Baltic, Hudson Bay, etc.) can
     * receive enough spring meltwater + runoff + precipitation to drive
     * surface salinity to near-zero locally. Once S < ~0.2 PSU, the JM
     * equation of state produces NaN density, NaN propagates into the SSH
     * CG solver (pp·App = NaN), and the model aborts. The 2-yr CORE2 run
     * hit this in spring 1959 at Gulf-of-Bothnia / south-Baltic nodes.
     *
     * The clamp is deterministic over myDim+eDim and idempotent, so it's
     * safe to apply without an extra halo exchange. Conservative floor of
     * 0.5 PSU — well below natural Baltic dynamics (typical surface SSS
     * 3-8 PSU there) so it only fires when the model is clearly broken.
     *
     * Net non-conservation introduced is <0.1% of global salt mass over
     * a multi-year run (verified post-hoc on the 2-yr CORE2 case).
     * Suppress with FESOM_NO_SFLOOR=1 to bisect. */
    {
        static int sf_checked = 0, sf_skip = 0;
        if (!sf_checked) {
            const char *e = getenv("FESOM_NO_SFLOOR");
            sf_skip = (e && atoi(e));
            sf_checked = 1;
        }
        if (!sf_skip) {
            const real_t S_FLOOR = (real_t)0.5;
            real_t *S = tracers->data[FESOM_TRACER_S].values;
            const int N_full = mesh->myDim_nod2D + mesh->eDim_nod2D;
            for (int n = 0; n < N_full; ++n) {
                const int nzmax = mesh->nlevels_nod2D[n] - 1;
                for (int nz = 0; nz < nzmax; ++nz) {
                    const size_t i = FESOM_NODE3D(n, nz, mesh->nl);
                    if (S[i] < S_FLOOR) S[i] = S_FLOOR;
                }
            }
        }
    }

    /* 13c. Phase G6b — bolus velocity sub (Fortran oce_ale_tracer.F90:284-295).
     * Restore dyn->uv / dyn->w / dyn->w_e to their pre-add values so the
     * remainder of the timestep (and the next step) sees the original
     * velocity field. */
    if (gm) {
        const int E_loop = mesh->myDim_elem2D + mesh->eDim_elem2D;
        const int N_loop = mesh->myDim_nod2D + mesh->eDim_nod2D;
        for (int e = 0; e < E_loop; ++e) {
            for (int nz = 0; nz < nl; ++nz) {
                size_t k = (size_t)e * nl * 2 + nz * 2;
                dyn->uv[k + 0] -= dyn->fer_uv[k + 0];
                dyn->uv[k + 1] -= dyn->fer_uv[k + 1];
            }
        }
        for (int n = 0; n < N_loop; ++n) {
            for (int nz = 0; nz < nl; ++nz) {
                size_t k = (size_t)n * nl + nz;
                dyn->w  [k] -= dyn->fer_w[k];
                dyn->w_e[k] -= dyn->fer_w[k];
            }
        }
    }

    /* 14. commit thickness: hnode := hnode_new, helem from vertex mean. */
    fesom_ale_commit_thickness(mesh);
    fesom_exchange_nod3D (mesh->hnode, nl, p);   /* both already same — but be explicit */
    fesom_exchange_elem3D(mesh->helem, nl, p);   /* Fortran oce_ale.F90:1027,1249 */

    /* Sea-ice step is now called from fesom_main BEFORE the ocean step
     * (ice writes heat_flux/water_flux that the ocean step consumes). */

    return cg_iters;
}
