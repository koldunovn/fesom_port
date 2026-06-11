#include "fesom_ice_maevp.h"
#include "fesom_constants.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * mEVP (whichEVP=1) — literal port of EVPdynamics_m, ice_maEVP.F90:429-882
 * (the NR-optimized routine with ssh2rhs / stress_tensor_m / stress2rhs_m
 * inlined). Branches compiled away by the reference config and excluded per
 * the plan: use_floatice ssh2rhs branch (we port the levitating else-branch,
 * which covers linfs AND zstar), use_cavity edge-BC block, __icepack blocks.
 * The ulevels/ulevels_nod2D cavity cycles ARE kept as written (no-op here).
 *
 * mEVP-vs-stdEVP asymmetries ported deliberately (plan fidelity-traps
 * checklist — do NOT "fix" them against the fesom_ice_evp.c template):
 *   rdt = FULL ice_dt (not /evp_rheol_steps) and the drag term carries rdt;
 *   pressure_fac has NO 0.5 — the 0.5 sits in the sigma11/22 updates only;
 *   NO theta_io rotation; element mask is mean-msum>0.01, node mask
 *   a_ice>=0.01; non-ice nodes are SKIPPED (velocity retained), not zeroed;
 *   uice_old/vice_old never touched; sigma11/12/22 NOT zeroed on entry
 *   (persist across calls); bc_index_nod2D multiplies det in the node solve
 *   (redundant with the edge-BC loop but both are ported).
 */

/* --- env-gated reference dumps (FESOM_EVP_DUMP_DIR) ---
 * Same file format + naming as the Fortran fesom_evp_dump module
 * (ice_EVP.F90) and the std-EVP C dumps; the mEVP Fortran working-copy
 * patch writes points Q/U0/F, P, it<substep>, UF with the same value
 * order. Own static call counter — only one rheology runs per execution. */
static int         s_maevp_call_step  = 0;
static const char *s_maevp_dump_dir   = NULL;
static int         s_maevp_dump_loaded = 0;

static void maevp_dump_load_env(void)
{
    if (s_maevp_dump_loaded) return;
    s_maevp_dump_dir = getenv("FESOM_EVP_DUMP_DIR");
    s_maevp_dump_loaded = 1;
}

static void maevp_dump(int step, const char *point, const char *cls,
                       const real_t * const *vals, int nvals,
                       int N, const int *gids, int rank)
{
    if (!s_maevp_dump_dir) return;
    char path[1024];
    snprintf(path, sizeof(path),
             "%s/evp_dump_s%d_%s_%s_rank%d.txt",
             s_maevp_dump_dir, step, point, cls, rank);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "maevp_dump: cannot open %s\n", path); return; }
    fprintf(f, "# step=%d point=%s array=%s rank=%d N=%d\n", step, point, cls, rank, N);
    for (int n = 0; n < N; ++n) {
        fprintf(f, "%d", gids[n]);
        for (int k = 0; k < nvals; ++k) fprintf(f, " %.17g", vals[k][n]);
        fputc('\n', f);
    }
    fclose(f);
}

void fesom_ice_evp_dynamics_m(fesom_ice            *ice,
                              struct fesom_partit  *partit,
                              struct fesom_mesh    *mesh)
{
    int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;   /* node_size */
    int nl = mesh->nl;

    /* pointer associations (ice_maEVP.F90:481-510) */
    real_t *u_ice  = ice->uice;
    real_t *v_ice  = ice->vice;
    real_t *a_ice  = ice->data[FESOM_ICE_AICE].values;
    real_t *m_ice  = ice->data[FESOM_ICE_MICE].values;
    real_t *m_snow = ice->data[FESOM_ICE_MSNOW].values;
    real_t *eps11  = ice->work.eps11;
    real_t *eps12  = ice->work.eps12;
    real_t *eps22  = ice->work.eps22;
    real_t *sigma11 = ice->work.sigma11;
    real_t *sigma12 = ice->work.sigma12;
    real_t *sigma22 = ice->work.sigma22;
    real_t *u_rhs_ice = ice->uice_rhs;
    real_t *v_rhs_ice = ice->vice_rhs;
    real_t *rhs_a  = ice->data[FESOM_ICE_AICE].values_rhs;  /* ssh-gradient scratch — */
    real_t *rhs_m  = ice->data[FESOM_ICE_MICE].values_rhs;  /* Fortran reuses exactly these */
    real_t *u_w    = ice->srfoce_u;
    real_t *v_w    = ice->srfoce_v;
    real_t *elevation       = ice->srfoce_ssh;
    real_t *stress_atmice_x = ice->stress_atmice_x;
    real_t *stress_atmice_y = ice->stress_atmice_y;
    real_t *u_ice_aux = ice->uice_aux;
    real_t *v_ice_aux = ice->vice_aux;
    real_t  rhoice = ice->thermo.rhoice;
    real_t  rhosno = ice->thermo.rhosno;

    /* per-call scratch (Fortran automatic arrays, heap-persistent in C) */
    real_t *inv_thickness = ice->work.mevp_inv_thickness;  /* [myDim_nod2D]  */
    real_t *mass          = ice->work.mevp_mass;           /* [myDim_nod2D]  */
    int    *ice_nod       = ice->work.mevp_ice_nod;        /* [myDim_nod2D]  */
    real_t *pressure_fac  = ice->work.mevp_pressure_fac;   /* [myDim_elem2D] */
    int    *ice_el        = ice->work.mevp_ice_el;         /* [myDim_elem2D] */

    /* setup (:513-518): rdt is the FULL ice step — std EVP divides by
     * evp_rheol_steps, mEVP does NOT (trap 1) */
    const real_t val3 = 1.0 / 3.0;
    const real_t vale = 1.0 / (ice->ellipse * ice->ellipse);
    const real_t det2 = 1.0 / (1.0 + ice->alpha_evp);
    const real_t det1 = ice->alpha_evp * det2;
    const real_t rdt  = ice->ice_dt;
    const int    steps = ice->evp_rheol_steps;

    /* u_ice_aux = u_ice over the FULL allocation extent (:521-522) — whole-
     * array Fortran assignment, node_size = myDim+eDim (trap 8) */
    for (int row = 0; row < N; ++row) {
        u_ice_aux[row] = u_ice[row];
        v_ice_aux[row] = v_ice[row];
    }

    /* --- reference dumps: entry inputs (mirror of the Fortran patch) --- */
    maevp_dump_load_env();
    ++s_maevp_call_step;
    int dump_now = (s_maevp_dump_dir != NULL) && (s_maevp_call_step <= 4);
    if (dump_now) {
        const real_t *qn[4] = { a_ice, m_ice, m_snow, elevation };
        maevp_dump(s_maevp_call_step, "Q", "node", qn, 4,
                   mesh->myDim_nod2D, partit->myList_nod2D, partit->mype);
        const real_t *uv0[2] = { u_ice, v_ice };
        maevp_dump(s_maevp_call_step, "U0", "node", uv0, 2,
                   mesh->myDim_nod2D, partit->myList_nod2D, partit->mype);
        const real_t *fv[4] = { stress_atmice_x, stress_atmice_y, u_w, v_w };
        maevp_dump(s_maevp_call_step, "F", "node", fv, 4,
                   mesh->myDim_nod2D, partit->myList_nod2D, partit->mype);
    }

    /* --- inlined ssh2rhs, levitating branch (:535-543, 589-621) ---
     * rhs_a/rhs_m zeroed over myDim ONLY (:539-542, trap 7); the element
     * assembly writes ALL 3 nodes UNGUARDED (:604-617, trap 6 — halo entries
     * written-never-read; arrays are node_size so in-bounds). */
    for (int row = 0; row < mesh->myDim_nod2D; ++row) {
        rhs_a[row] = 0.0;
        rhs_m[row] = 0.0;
    }
    /* (use_floatice && !linfs branch :547-587 not ported — use_floatice=false;
     * the else-branch below is the executed path for linfs AND zstar) */
    for (int el = 0; el < mesh->myDim_elem2D; ++el) {
        /* if element has any cavity node skip it (:596) */
        if (mesh->ulevels[el] > 1) continue;
        int n0 = mesh->elem_nodes[3*el + 0];
        int n1 = mesh->elem_nodes[3*el + 1];
        int n2 = mesh->elem_nodes[3*el + 2];
        real_t vol = mesh->elem_area[el];
        real_t *gs = &mesh->gradient_sca[6*el];   /* dx=gs[0..2], dy=gs[3..5] */
        real_t bb = FESOM_G * val3 * vol;
        real_t aa = bb * (gs[0]*elevation[n0] + gs[1]*elevation[n1] + gs[2]*elevation[n2]);
        bb        = bb * (gs[3]*elevation[n0] + gs[4]*elevation[n1] + gs[5]*elevation[n2]);
        rhs_a[n0] -= aa;  rhs_m[n0] -= bb;
        rhs_a[n1] -= aa;  rhs_m[n1] -= bb;
        rhs_a[n2] -= aa;  rhs_m[n2] -= bb;
    }

    /* --- node precompute (:624-647, myDim): masks + mass regularization.
     * Zero-then-cycle order as in Fortran; rhs scaling happens INSIDE the
     * ice branch only (:641-642, trap 7); mass = M/((1+M²)·area) verbatim
     * (trap 5); max(·, 9.0) limiter (:635). */
    for (int i = 0; i < mesh->myDim_nod2D; ++i) {
        inv_thickness[i] = 0.0;
        mass[i] = 0.0;
        ice_nod[i] = 0;
        /* if cavity node skip it (:631) */
        if (mesh->ulevels_nod2D[i] > 1) continue;
        if (a_ice[i] >= 0.01) {
            real_t it = (rhoice * m_ice[i] + rhosno * m_snow[i]) / a_ice[i];
            inv_thickness[i] = 1.0 / ((it > 9.0) ? it : 9.0);   /* limit the mass */

            real_t ms = m_ice[i] * rhoice + m_snow[i] * rhosno;
            mass[i] = ms / ((1.0 + ms * ms) * mesh->area[i*nl + 0]);

            /* scale rhs_a, rhs_m, too */
            rhs_a[i] = rhs_a[i] / mesh->area[i*nl + 0];
            rhs_m[i] = rhs_m[i] / mesh->area[i*nl + 0];

            ice_nod[i] = 1;
        }
    }

    /* --- element precompute (:650-667): mean-msum element mask (trap 4) and
     * pressure_fac WITH det2 and NO 0.5 (trap 2 — the 0.5 of std-EVP
     * ice_strength does NOT belong here; it sits in the sigma11/22 update). */
    for (int el = 0; el < mesh->myDim_elem2D; ++el) {
        pressure_fac[el] = 0.0;
        ice_el[el] = 0;
        /* if element has any cavity node skip it (:658) */
        if (mesh->ulevels[el] > 1) continue;
        int n0 = mesh->elem_nodes[3*el + 0];
        int n1 = mesh->elem_nodes[3*el + 1];
        int n2 = mesh->elem_nodes[3*el + 2];
        real_t msum = (m_ice[n0] + m_ice[n1] + m_ice[n2]) * val3;
        if (msum > 0.01) {
            ice_el[el] = 1;
            real_t asum = (a_ice[n0] + a_ice[n1] + a_ice[n2]) * val3;
            pressure_fac[el] = det2 * ice->pstar * msum
                               * exp(-ice->c_pressure * (1.0 - asum));
        }
    }
    /* zero u_rhs/v_rhs over myDim (:668-673) */
    for (int row = 0; row < mesh->myDim_nod2D; ++row) {
        u_rhs_ice[row] = 0.0;
        v_rhs_ice[row] = 0.0;
    }

    if (dump_now) {
        const real_t *pn[4] = { inv_thickness, mass, rhs_a, rhs_m };
        maevp_dump(s_maevp_call_step, "P", "node", pn, 4,
                   mesh->myDim_nod2D, partit->myList_nod2D, partit->mype);
        const real_t *pe[1] = { pressure_fac };
        maevp_dump(s_maevp_call_step, "P", "elem", pe, 1,
                   mesh->myDim_elem2D, partit->myList_elem2D, partit->mype);
    }

    /* --- mEVP pseudo-time iteration main loop (:680-875) --- */
    for (int shortstep = 1; shortstep <= steps; ++shortstep) {

        /* inlined stress_tensor_m + stress2rhs_m: element loop (:689-792) */
        for (int el = 0; el < mesh->myDim_elem2D; ++el) {
            if (mesh->ulevels[el] > 1) continue;     /* (:690) */
            if (!ice_el[el]) continue;               /* (:693) */

            int eln[3] = { mesh->elem_nodes[3*el + 0],
                           mesh->elem_nodes[3*el + 1],
                           mesh->elem_nodes[3*el + 2] };
            real_t *gs = &mesh->gradient_sca[6*el];  /* dx=gs[0..2], dy=gs[3..5] */
            /* METRICS (:699) */
            real_t meancos = val3 * mesh->metric_factor[el];
            real_t U[3] = { u_ice_aux[eln[0]], u_ice_aux[eln[1]], u_ice_aux[eln[2]] };
            real_t V[3] = { v_ice_aux[eln[0]], v_ice_aux[eln[1]], v_ice_aux[eln[2]] };

            /* deformation rate tensor (:702-705) */
            eps11[el] = gs[0]*U[0] + gs[1]*U[1] + gs[2]*U[2]
                        - (V[0] + V[1] + V[2]) * meancos;             /* metrics */
            eps22[el] = gs[3]*V[0] + gs[4]*V[1] + gs[5]*V[2];
            eps12[el] = 0.5 * ( gs[3]*U[0] + gs[4]*U[1] + gs[5]*U[2]
                              + gs[0]*V[0] + gs[1]*V[1] + gs[2]*V[2]
                              + (U[0] + U[1] + U[2]) * meancos );     /* metrics */

            /* switch to eps1, eps2 (:708-709) */
            real_t eps1 = eps11[el] + eps22[el];
            real_t eps2 = eps11[el] - eps22[el];

            /* moduli (:712); delta_min is ADDITIVE (:714, trap 10) */
            real_t delta = sqrt(eps1*eps1 + vale*(eps2*eps2 + 4.0*eps12[el]*eps12[el]));
            real_t pressure = pressure_fac[el] / (delta + ice->delta_min);

            /* implicit sigma update (:721-723): the 0.5 lives in the
             * sigma11/22 lines, NOT in sigma12, NOT in pressure_fac (trap 2) */
            sigma12[el] = det1*sigma12[el] +     pressure*eps12[el]*vale;
            sigma11[el] = det1*sigma11[el] + 0.5*pressure*(eps1 - delta + eps2*vale);
            sigma22[el] = det1*sigma22[el] + 0.5*pressure*(eps1 - delta - eps2*vale);

            /* inlined stress2rhs_m: stress-divergence assembly GUARDED to
             * owned nodes (:740-790, trap 6 asymmetry vs the unguarded
             * ssh2rhs assembly above) */
            for (int k = 0; k < 3; ++k) {
                if (eln[k] < mesh->myDim_nod2D) {
                    u_rhs_ice[eln[k]] -= mesh->elem_area[el] *
                        (sigma11[el]*gs[k]   + sigma12[el]*gs[k+3] + sigma12[el]*meancos);  /* metrics */
                    v_rhs_ice[eln[k]] -= mesh->elem_area[el] *
                        (sigma12[el]*gs[k]   + sigma22[el]*gs[k+3] - sigma11[el]*meancos);  /* metrics */
                }
            }
        }

        /* node solve (:795-819): SKIP non-ice nodes — NO else branch; a
         * non-ice node keeps its u_ice_aux entry value (trap 13 — std EVP
         * zeroes here; mEVP must NOT). bc_index_nod2D multiplies det
         * (:814, trap 9). Drag CARRIES rdt; drag·u_w sits OUTSIDE the
         * rdt·(...) group; +beta·u_aux closes the rhs (:806-811, trap 1).
         * NO theta_io rotation anywhere (trap 3). */
        for (int i = 0; i < mesh->myDim_nod2D; ++i) {
            if (mesh->ulevels_nod2D[i] > 1) continue;    /* (:797) */
            if (ice_nod[i]) {                            /* skip if ice is absent */
                u_rhs_ice[i] = u_rhs_ice[i]*mass[i] + rhs_a[i];
                v_rhs_ice[i] = v_rhs_ice[i]*mass[i] + rhs_m[i];

                real_t du = u_ice_aux[i] - u_w[i];
                real_t dv = v_ice_aux[i] - v_w[i];
                real_t umod = sqrt(du*du + dv*dv);
                real_t drag = rdt * ice->cd_oce_ice * umod * FESOM_DENSITY_0
                              * inv_thickness[i];

                /* rhs for water stress, air stress, u_rhs_ice (internal stress + ssh) */
                real_t rhsu = u_ice[i] + drag*u_w[i]
                              + rdt*(inv_thickness[i]*stress_atmice_x[i] + u_rhs_ice[i])
                              + ice->beta_evp*u_ice_aux[i];
                real_t rhsv = v_ice[i] + drag*v_w[i]
                              + rdt*(inv_thickness[i]*stress_atmice_y[i] + v_rhs_ice[i])
                              + ice->beta_evp*v_ice_aux[i];

                /* solve (Coriolis and water stress treated implicitly) */
                real_t obd = 1.0 + ice->beta_evp + drag;
                real_t rf  = rdt * mesh->coriolis_node[i];
                real_t det = mesh->bc_index_nod2D[i] / (obd*obd + rf*rf);

                u_ice_aux[i] = det * (obd*rhsu + rf*rhsv);
                v_ice_aux[i] = det * (obd*rhsv - rf*rhsu);
            }
        }

        /* edge BC (:823-857): coastal boundary edges -> zero aux velocity at
         * both endpoints. Same global-id convention as the std-EVP BC loop
         * and the M1 bc_index_nod2D build. (use_cavity block :842-856 not
         * ported — use_cavity=false.) */
        for (int ed = 0; ed < mesh->myDim_edge2D; ++ed) {
            if (partit->myList_edge2D[ed] > mesh->edge2D_in) {
                int e0 = mesh->edges[ed*2 + 0];
                int e1 = mesh->edges[ed*2 + 1];
                u_ice_aux[e0] = 0.0;  v_ice_aux[e0] = 0.0;
                u_ice_aux[e1] = 0.0;  v_ice_aux[e1] = 0.0;
            }
        }

        /* reference dump: per-substep state AFTER node solve + edge BC,
         * BEFORE the halo exchange — owned entries only (mirrors the
         * Fortran patch placement) */
        if (dump_now && (shortstep == 1 || shortstep == 2 ||
                         shortstep == 60 || shortstep == steps)) {
            char itlab[16];
            snprintf(itlab, sizeof(itlab), "it%d", shortstep);
            const real_t *uvit[2] = { u_ice_aux, v_ice_aux };
            maevp_dump(s_maevp_call_step, itlab, "node", uvit, 2,
                       mesh->myDim_nod2D, partit->myList_nod2D, partit->mype);
            const real_t *sgit[3] = { sigma11, sigma12, sigma22 };
            maevp_dump(s_maevp_call_step, itlab, "elem", sgit, 3,
                       mesh->myDim_elem2D, partit->myList_elem2D, partit->mype);
        }

        /* halo exchange of (u_ice_aux, v_ice_aux): Fortran overlaps
         * exchange_nod_begin/end around the rhs zeroing (:861-872); two
         * blocking exchanges + the zero loop are result-identical (the zero
         * loop touches only u_rhs/v_rhs, not the exchanged fields). */
        fesom_exchange_nod2D(u_ice_aux, partit);
        fesom_exchange_nod2D(v_ice_aux, partit);
        for (int row = 0; row < mesh->myDim_nod2D; ++row) {
            u_rhs_ice[row] = 0.0;
            v_rhs_ice[row] = 0.0;
        }
    }

    /* final copy over myDim+eDim (:876-881, trap 8) — no extra exchange:
     * the aux halo is current from the last substep's exchange. */
    for (int row = 0; row < N; ++row) {
        u_ice[row] = u_ice_aux[row];
        v_ice[row] = v_ice_aux[row];
    }

    if (dump_now) {
        const real_t *uvf[2] = { u_ice, v_ice };
        maevp_dump(s_maevp_call_step, "UF", "node", uvf, 2,
                   mesh->myDim_nod2D, partit->myList_nod2D, partit->mype);
    }
}
