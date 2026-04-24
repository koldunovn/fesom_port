#include "fesom_ale.h"
#include "fesom_aux.h"
#include "fesom_bulk.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_eos.h"
#include "fesom_forcing.h"
#include "fesom_forcing_analytical.h"
#include "fesom_ic.h"
#include "fesom_io.h"
#include "fesom_halo.h"
#include "fesom_jra55.h"
#include "fesom_mesh.h"
#include "fesom_sss_runoff.h"
#include "fesom_momentum.h"
#include "fesom_mpi.h"
#include "fesom_phc.h"
#include "fesom_pp.h"
#include "fesom_ssh.h"
#include "fesom_step.h"
#include "fesom_tracer_adv.h"
#include "fesom_tracers.h"
#include "fesom_types.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_sanity(const fesom_mesh *m)
{
    /*
     * FRESH_START.md §20 expectations:
     *   - pi:     nod2D=3140  elem2D=5839  nl=~23
     *   - CORE2:  nod2D=126858 elem2D=244659 nl=48
     *   - first 3 nodes' coords are lat/lon in degrees
     *   - first 3 elements' nodes are 1-based on disk → here 0-based
     *   - all elem_nodes in [0, nod2D)
     */
    printf("[fesom_port] mesh counts: nod2D=%d  elem2D=%d  edge2D=%d  nl=%d\n",
           m->nod2D, m->elem2D, m->edge2D, m->nl);
    printf("[fesom_port] partition  : myDim_nod2D=%d  eDim_nod2D=%d  myDim_elem2D=%d\n",
           m->myDim_nod2D, m->eDim_nod2D, m->myDim_elem2D);

    /* zbar profile (first/last few) */
    printf("[fesom_port] zbar:");
    int zshow = m->nl < 6 ? m->nl : 6;
    for (int k = 0; k < zshow; ++k) printf("  %g", (double)m->zbar[k]);
    if (m->nl > 6) printf("  ...  %g", (double)m->zbar[m->nl - 1]);
    printf("\n");

    /* Show first 3 nodes: rotated-radian + geographic-degree */
    int show = m->nod2D < 3 ? m->nod2D : 3;
    for (int i = 0; i < show; ++i) {
        double rlon = (double)m->coord_nod2D[2*i + 0] / FESOM_RAD;
        double rlat = (double)m->coord_nod2D[2*i + 1] / FESOM_RAD;
        double glon = (double)m->geo_coord_nod2D[2*i + 0] / FESOM_RAD;
        double glat = (double)m->geo_coord_nod2D[2*i + 1] / FESOM_RAD;
        printf("  node %d: geo=(%9.4f,%9.4f)  rot=(%9.4f,%9.4f)  coast=%d  depth=%g  nlev=%d\n",
               i, glon, glat, rlon, rlat,
               m->coast_flag[i], (double)m->depth[i], m->nlevels_nod2D[i]);
    }

    show = m->elem2D < 3 ? m->elem2D : 3;
    for (int e = 0; e < show; ++e) {
        printf("  elem %d: nodes (%d, %d, %d)  nlev=%d\n", e,
               m->elem_nodes[3*e + 0],
               m->elem_nodes[3*e + 1],
               m->elem_nodes[3*e + 2],
               m->nlevels[e]);
    }

    /* Show first 3 edges and adjacency */
    show = m->edge2D < 3 ? m->edge2D : 3;
    for (int eg = 0; eg < show; ++eg) {
        printf("  edge %d: nodes (%d, %d)  tri (%d, %d)\n", eg,
               m->edges[2*eg + 0],   m->edges[2*eg + 1],
               m->edge_tri[2*eg + 0], m->edge_tri[2*eg + 1]);
    }

    /* §20 sanity — local extents (myDim+eDim for nodes, myDim only for elem
     * connectivity since halo elements don't carry elem_nodes). */
    int N_loc  = m->myDim_nod2D  + m->eDim_nod2D;
    int Ei_loc = m->myDim_elem2D;                                /* elem_nodes range */
    int EG_loc = m->myDim_edge2D + m->eDim_edge2D;
    int oob_e = 0, oob_t = 0, edge_boundaries = 0;
    int nlmin = m->nlevels_nod2D[0], nlmax = nlmin;
    for (int e = 0; e < Ei_loc; ++e) {
        for (int k = 0; k < 3; ++k) {
            int v = m->elem_nodes[3*e + k];
            if (v < 0 || v >= N_loc) ++oob_e;
        }
    }
    int E_loc = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    for (int i = 0; i < EG_loc; ++i) {
        for (int k = 0; k < 2; ++k) {
            int v = m->edge_tri[2*i + k];
            if (v == -1) ++edge_boundaries;
            else if (v < -1 || v >= E_loc) ++oob_t;
        }
    }
    for (int i = 1; i < N_loc; ++i) {
        if (m->nlevels_nod2D[i] < nlmin) nlmin = m->nlevels_nod2D[i];
        if (m->nlevels_nod2D[i] > nlmax) nlmax = m->nlevels_nod2D[i];
    }
    printf("[fesom_port] sanity: elem_nodes_oob=%d (must be 0); "
           "edge_tri_oob=%d (must be 0); edge_boundary_endpoints=%d\n",
           oob_e, oob_t, edge_boundaries);
    printf("[fesom_port] sanity: nlevels_nod2D min=%d  max=%d  (nl=%d)\n",
           nlmin, nlmax, m->nl);

    if (m->nlevels_nod2D_min) {
        int kvmin_min = m->nlevels_nod2D_min[0], kvmin_max = kvmin_min;
        for (int i = 1; i < N_loc; ++i) {
            int v = m->nlevels_nod2D_min[i];
            if (v < kvmin_min) kvmin_min = v;
            if (v > kvmin_max) kvmin_max = v;
        }
        printf("[fesom_port] sanity: nlevels_nod2D_min (K_v⁻) min=%d max=%d\n",
               kvmin_min, kvmin_max);
    }

    /* elem_area: range over INTERIOR elements only (halo got 0 then
     * exchange — but only if multi-rank; for 1-rank, halo doesn't exist). */
    if (m->elem_area && Ei_loc > 0) {
        real_t amin = m->elem_area[0], amax = amin, asum = 0.0;
        int eminidx = 0, emaxidx = 0;
        for (int e = 0; e < Ei_loc; ++e) {
            real_t a = m->elem_area[e];
            asum += a;
            if (a < amin) { amin = a; eminidx = e; }
            if (a > amax) { amax = a; emaxidx = e; }
        }
        real_t amean = asum / (real_t)Ei_loc;
        printf("[fesom_port] elem_area (m²): min=%.4e (elem %d)  max=%.4e (elem %d)  mean=%.4e\n",
               (double)amin, eminidx, (double)amax, emaxidx, (double)amean);
        printf("[fesom_port] equiv resolution: min=%.2f km  max=%.2f km  mean=%.2f km\n",
               sqrt(2.0 * (double)amin) / 1000.0,
               sqrt(2.0 * (double)amax) / 1000.0,
               sqrt(2.0 * (double)amean) / 1000.0);
    }
    if (m->areasvol && m->nl > 0) {
        printf("[fesom_port] areasvol[node0, nz=0]=%.4e m²  area[node0, nz=0]=%.4e m²\n",
               (double)m->areasvol[0], (double)m->area[0]);
    }

    if (m->coriolis_node) {
        real_t fmin = m->coriolis_node[0], fmax = fmin;
        for (int n = 1; n < N_loc; ++n) {
            if (m->coriolis_node[n] < fmin) fmin = m->coriolis_node[n];
            if (m->coriolis_node[n] > fmax) fmax = m->coriolis_node[n];
        }
        printf("[fesom_port] coriolis_node (1/s): min=%.4e  max=%.4e\n",
               (double)fmin, (double)fmax);
    }

    if (m->elem_cos && Ei_loc > 0) {
        real_t cmin = m->elem_cos[0], cmax = cmin;
        for (int e = 1; e < Ei_loc; ++e) {
            if (m->elem_cos[e] < cmin) cmin = m->elem_cos[e];
            if (m->elem_cos[e] > cmax) cmax = m->elem_cos[e];
        }
        printf("[fesom_port] elem_cos: min=%.4f  max=%.4f  (must be in (0,1])\n",
               (double)cmin, (double)cmax);
    }

    if (m->edge_dxdy && EG_loc > 0) {
        real_t lmax = 0.0;
        for (int e = 0; e < EG_loc; ++e) {
            real_t dx = m->edge_dxdy[2*e + 0];
            real_t dy = m->edge_dxdy[2*e + 1];
            real_t l  = sqrt(dx*dx + dy*dy);
            if (l > lmax) lmax = l;
        }
        printf("[fesom_port] max |edge_dxdy| (radians) = %.4e  (~%.1f km)\n",
               (double)lmax, (double)lmax * FESOM_R_EARTH / 1000.0);
    }

    if (m->edge_cross_dxdy && EG_loc > 0) {
        printf("[fesom_port] edge_cross_dxdy[0] (m): "
               "left=(%.2e, %.2e)  right=(%.2e, %.2e)\n",
               (double)m->edge_cross_dxdy[0],
               (double)m->edge_cross_dxdy[1],
               (double)m->edge_cross_dxdy[2],
               (double)m->edge_cross_dxdy[3]);
    }

    if (m->gradient_sca && Ei_loc > 0) {
        real_t residual_max = 0.0;
        for (int e = 0; e < Ei_loc; ++e) {
            real_t sx = m->gradient_sca[6*e + 0]
                      + m->gradient_sca[6*e + 1]
                      + m->gradient_sca[6*e + 2];
            real_t sy = m->gradient_sca[6*e + 3]
                      + m->gradient_sca[6*e + 4]
                      + m->gradient_sca[6*e + 5];
            real_t r = (fabs(sx) > fabs(sy)) ? fabs(sx) : fabs(sy);
            if (r > residual_max) residual_max = r;
        }
        printf("[fesom_port] gradient_sca: max |Σ ∂N_i/∂{x,y}| = %.3e (must be ≈ 0)\n",
               (double)residual_max);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <mesh_dir> [output_dir] [dt_seconds] [nsteps] [snap_every] [phc_nc_path] [jra55_year]\n",
                argv[0]);
        fprintf(stderr, "  jra55_year > 0 → use JRA55-do daily forcing for that year (replaces analytical wind)\n");
        return 2;
    }
    const char *mesh_dir = argv[1];
    const char *out_dir  = (argc >= 3 && argv[2][0] != '\0') ? argv[2] : NULL;
    if (argc >= 4) fesom_phase1_dt = (real_t)atof(argv[3]);
    int nsteps_cli     = (argc >= 5) ? atoi(argv[4]) : 0;
    int snap_every_cli = (argc >= 6) ? atoi(argv[5]) : 0;
    const char *phc_path = (argc >= 7 && argv[6][0] != '\0') ? argv[6] : NULL;
    int jra55_year = (argc >= 8) ? atoi(argv[7]) : 0;  /* 0 disables JRA55 */
    if (out_dir) printf("[fesom_port] snapshots → %s/snap_NNNNNN.nc\n", out_dir);
    printf("[fesom_port] dt = %.1f s\n", (double)fesom_phase1_dt);
    if (phc_path) printf("[fesom_port] PHC IC source: %s\n", phc_path);
    if (jra55_year > 0) printf("[fesom_port] JRA55-do forcing year: %d\n", jra55_year);

    fesom_mpi mpi;
    fesom_mpi_init(&mpi, mesh_dir, argc, argv);

    fesom_mesh mesh;
    fesom_mesh_init(&mesh);
    /* Rank 0 reads global mesh files; others fill via Bcast then extract
     * their slice using partit->myList_*. For npes==1 this is identity. */
    /* But we need partit->myDim_* etc set first. For npes>1 it's populated
     * by fesom_partit_init from my_list*.out; for npes==1 we synthesise them
     * AFTER reading mesh dims (so the synthesise call comes between read and
     * scatter). To handle both: do read+scatter inline. */
    if (mpi.npes == 1) {
        /* Synthesised partit needs global counts up front (before any halo
         * code can run). Read global on rank 0 directly. */
        fesom_mesh_read(&mesh, mesh_dir, &mpi);
        fesom_partit_set_global_counts_serial(&mpi, mesh.nod2D, mesh.elem2D, mesh.edge2D);
    } else {
        fesom_mesh_read(&mesh, mesh_dir, &mpi);
    }
    fesom_mesh_compute_metrics(&mesh, &mpi);
    fesom_mesh_alloc_state(&mesh);

    /* Phase 4 step 30b: identity test of halo exchange. No-op on 1 rank.
     * On 2+ ranks: positive test (halo == owning rank), negative test
     * (corruption recovery). Aborts on mismatch. */
    fesom_halo_identity_test(&mpi);

    print_sanity(&mesh);

    /* Phase 1 step 2: allocate the rest of the model state. All zeros until
       step 3 sets initial conditions. */
    fesom_dyn     dyn;     fesom_dyn_alloc    (&dyn,    &mesh);
    fesom_tracers tracers; fesom_tracers_alloc(&tracers, &mesh);
    fesom_aux     aux;     fesom_aux_alloc    (&aux,    &mesh);
    fesom_forcing forcing; fesom_forcing_alloc(&forcing, &mesh);

    /* Cheap accounting: total bytes allocated for the 3D state arrays (rough). */
    size_t bytes = 0;
    bytes += (size_t)mesh.nod2D  * mesh.nl  * sizeof(real_t) * 12; /* aux+tracer+mesh */
    bytes += (size_t)mesh.elem2D * mesh.nl  * sizeof(real_t) * 6;  /* dyn+aux at elem */
    bytes += (size_t)mesh.elem2D * mesh.nl  * 2 * sizeof(real_t) * 3;
    printf("[fesom_port] state alloc: dyn AB_order=%d, tracers=%d, "
           "approx 3D bytes=%.1f MiB\n",
           dyn.AB_order, tracers.num_tracers, (double)bytes / (1024.0 * 1024.0));

    /* Phase 1 step 3: initial conditions. */
    fesom_ic_thickness(&mesh, &mesh, &dyn);
    fesom_ic_tracers_constant(&mesh, &tracers, 10.0, 35.0);
    /* (Phase 2 T blob is added below, AFTER the IC/rest-state sanity tests so
       those still see the unperturbed constant field.) */

    /* IC sanity: T/S at probe node 0 should be exactly 10 / 35 in valid layers,
       and 0 below. Column sum of hnode at probe node 0 should equal
       -zbar[nlevels_nod2D[0] - 1] (since zbar[0]=0). */
    {
        int n = 0;
        int nlev_n = mesh.nlevels_nod2D[n];
        double Tsurf = (double)tracers.data[FESOM_TRACER_T].values[FESOM_NODE3D(n, 0, mesh.nl)];
        double Ssurf = (double)tracers.data[FESOM_TRACER_S].values[FESOM_NODE3D(n, 0, mesh.nl)];
        double Tbelow = (nlev_n - 1 < mesh.nl)
            ? (double)tracers.data[FESOM_TRACER_T].values[FESOM_NODE3D(n, nlev_n - 1, mesh.nl)]
            : 0.0;
        double hsum = 0.0;
        for (int nz = 0; nz < nlev_n - 1; ++nz) {
            hsum += (double)mesh.hnode[FESOM_NODE3D(n, nz, mesh.nl)];
        }
        double depth_from_zbar = -(double)mesh.zbar[nlev_n - 1];
        printf("[fesom_port] IC at node 0 (nlev=%d): T_surf=%.3f  S_surf=%.3f  "
               "T_below=%.3g  Σhnode=%.3f m  -zbar[nlev-1]=%.3f m\n",
               nlev_n, Tsurf, Ssurf, Tbelow, hsum, depth_from_zbar);
    }
    /* Global checks: T should be exactly 10.0 in every valid layer, hbar should be 0. */
    {
        int t_bad = 0, s_bad = 0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 0; nz < nzmax; ++nz) {
                if (tracers.data[FESOM_TRACER_T].values[FESOM_NODE3D(n, nz, mesh.nl)] != 10.0) ++t_bad;
                if (tracers.data[FESOM_TRACER_S].values[FESOM_NODE3D(n, nz, mesh.nl)] != 35.0) ++s_bad;
            }
        }
        real_t hbar_max = 0.0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            real_t a = fabs(mesh.hbar[n]);
            if (a > hbar_max) hbar_max = a;
        }
        printf("[fesom_port] IC global: T≠10 count=%d  S≠35 count=%d  max|hbar|=%g (must all be 0)\n",
               t_bad, s_bad, (double)hbar_max);
    }
    /* Cell-vs-vertex consistency: helem[e][0] should equal mean of hnode at the 3 vertices
       (full-cell linfs → trivially equal; verify it nonetheless). */
    {
        real_t maxdiff = 0.0;
        for (int e = 0; e < mesh.elem2D; ++e) {
            int nzmax = mesh.nlevels[e] - 1;
            for (int nz = 0; nz < nzmax; ++nz) {
                real_t hsum = 0.0;
                for (int k = 0; k < 3; ++k) {
                    int v = mesh.elem_nodes[3*e + k];
                    hsum += mesh.hnode[FESOM_NODE3D(v, nz, mesh.nl)];
                }
                real_t hmean = hsum / 3.0;
                real_t d = fabs(mesh.helem[FESOM_ELEM3D(e, nz, mesh.nl)] - hmean);
                if (d > maxdiff) maxdiff = d;
            }
        }
        printf("[fesom_port] helem vs (1/3) Σ hnode: max diff = %.3e (should be ~0 for full cells)\n",
               (double)maxdiff);
    }

    /* Phase 1 step 4: JM-EOS + pressure_bv. */
    fesom_pressure_bv(&tracers, &mesh, &aux);
    /* Phase 1 step 5: PGF (linfs + full cells branch). */
    fesom_pressure_force_linfs_fullcell(&mesh, &aux);

    /* PGF sanity: at constant T/S, hpressure has no horizontal variation, so
       pgf should be exactly 0 (modulo machine epsilon × the partition-of-unity
       residual we measured ~1e-20). Anything materially larger means a bug. */
    {
        real_t pgmax = 0.0;
        size_t total_e = (size_t)mesh.elem2D * (size_t)mesh.nl;
        for (size_t i = 0; i < total_e; ++i) {
            real_t a = fabs(aux.pgf_x[i]);
            if (a > pgmax) pgmax = a;
            a = fabs(aux.pgf_y[i]);
            if (a > pgmax) pgmax = a;
        }
        printf("[fesom_port] PGF sanity (constant T/S): max|pgf_{x,y}| = %.3e m/s² "
               "(must be ≈ 0; ~1e-18 is normal)\n", (double)pgmax);
    }

    /* Phase 1 step 6: linfs ALE thickness + vertical velocity. */
    fesom_ale_thickness_linfs(&mesh);
    fesom_ale_vert_vel_linfs(&mesh, &dyn);

    /* ALE sanity: with UV ≡ 0, w must be exactly 0. hnode_new must equal hnode. */
    {
        real_t wmax = 0.0;
        size_t total = (size_t)mesh.nod2D * (size_t)mesh.nl;
        for (size_t i = 0; i < total; ++i) {
            real_t a = fabs(dyn.w[i]);
            if (a > wmax) wmax = a;
        }
        int diff = memcmp(mesh.hnode_new, mesh.hnode, total * sizeof(real_t));
        printf("[fesom_port] ALE sanity (UV=0): max|w| = %.3e m/s (must be 0)  "
               "hnode_new ≡ hnode: %s\n",
               (double)wmax, diff == 0 ? "yes" : "NO");
    }

    /* Phase 1 step 7: SSH stiffness matrix + preconditioner + CG solver. */
    fesom_ssh_stiff   stiff;
    fesom_solverinfo  solver;
    fesom_ssh_stiff_alloc_and_build(&stiff, &mesh);
    fesom_ssh_preconditioner(&stiff, &mesh);
    fesom_solverinfo_alloc(&solver, &mesh);

    /* Stiffness diagonal sanity — should mirror Fortran's "STIFF diag" debug
       block in init_stiff_mat_ale (lines 1584-1599). All diagonals must be
       positive (mass term + sum of -gradient² style off-diag couplings; the
       MITgcm preconditioner only converges if so). */
    {
        int neg = 0;
        real_t dmin = stiff.values[stiff.rowptr[0]];
        real_t dmax = dmin;
        real_t mass_min = 1e30, mass_max = -1e30;
        for (int row = 0; row < mesh.nod2D; ++row) {
            real_t d = stiff.values[stiff.rowptr[row]];
            if (d < 0.0)  ++neg;
            if (d < dmin) dmin = d;
            if (d > dmax) dmax = d;
            real_t mass = mesh.areasvol[FESOM_NODE3D(row, 0, mesh.nl)] / FESOM_PHASE1_DT;
            if (mass < mass_min) mass_min = mass;
            if (mass > mass_max) mass_max = mass;
        }
        printf("[fesom_port] SSH stiff: nnz=%d  diag neg=%d  range=[%.4e, %.4e]  "
               "mass=[%.4e, %.4e]\n",
               stiff.nnz, neg, (double)dmin, (double)dmax,
               (double)mass_min, (double)mass_max);
    }

    /* Symmetry check: matrix must be symmetric within machine epsilon for
       CG to converge. For each (row, col) pair with row<col, find matching
       (col, row) entry and compare. */
    {
        real_t asym_max = 0.0;
        int    asym_pairs = 0;
        for (int row = 0; row < mesh.nod2D; ++row) {
            for (int n = stiff.rowptr[row]; n < stiff.rowptr[row+1]; ++n) {
                int col = stiff.colind[n];
                if (col <= row) continue;
                real_t a_rc = stiff.values[n];
                real_t a_cr = 0.0;
                int found = 0;
                for (int m = stiff.rowptr[col]; m < stiff.rowptr[col+1]; ++m) {
                    if (stiff.colind[m] == row) {
                        a_cr = stiff.values[m];
                        found = 1;
                        break;
                    }
                }
                if (!found) { ++asym_pairs; continue; }
                real_t d = fabs(a_rc - a_cr);
                if (d > asym_max) asym_max = d;
            }
        }
        printf("[fesom_port] SSH stiff symmetry: max|A_rc - A_cr| = %.3e  "
               "missing-pairs = %d\n", (double)asym_max, asym_pairs);
    }

    /* RHS=0 test: solve Ax=0, expect d_eta=0 in 0 iterations.
       Currently UV=0 and ssh_rhs_old=0, so compute_ssh_rhs gives ssh_rhs=0. */
    fesom_compute_ssh_rhs_linfs(&mesh, &dyn);
    {
        real_t rhs_max = 0.0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            real_t a = fabs(dyn.ssh_rhs[n]);
            if (a > rhs_max) rhs_max = a;
        }
        int iters = fesom_ssh_solve_cg(&stiff, &solver, &mesh, &dyn);
        real_t deta_max = 0.0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            real_t a = fabs(dyn.d_eta[n]);
            if (a > deta_max) deta_max = a;
        }
        printf("[fesom_port] CG zero-RHS test: max|ssh_rhs|=%.3e  iters=%d  "
               "max|d_eta|=%.3e (must be 0)\n",
               (double)rhs_max, iters, (double)deta_max);
    }

    /* Bump-RHS test: inject a single non-zero on the RHS and solve.
       Pick an interior node (rough proxy: max node id minus a quarter) and
       set ssh_rhs to area*1.0 there. Should converge in O(10-50) iters and
       give d_eta with the magnitude near 1/diagonal. */
    {
        int n_probe = mesh.nod2D / 2;
        memset(dyn.ssh_rhs, 0, (size_t)mesh.nod2D * sizeof(real_t));
        memset(dyn.d_eta,   0, (size_t)mesh.nod2D * sizeof(real_t));
        real_t area_n = mesh.areasvol[FESOM_NODE3D(n_probe, 0, mesh.nl)];
        dyn.ssh_rhs[n_probe] = area_n;          /* unit-amplitude perturbation */

        int iters = fesom_ssh_solve_cg(&stiff, &solver, &mesh, &dyn);
        real_t deta_max = 0.0, deta_at_probe = dyn.d_eta[n_probe];
        for (int n = 0; n < mesh.nod2D; ++n) {
            real_t a = fabs(dyn.d_eta[n]);
            if (a > deta_max) deta_max = a;
        }
        /* Verify residual: ‖rhs - A x‖ / ‖rhs‖ should be ≤ soltol */
        real_t *Ax = malloc((size_t)mesh.nod2D * sizeof(real_t));
        for (int row = 0; row < mesh.nod2D; ++row) {
            real_t s = 0;
            for (int j = stiff.rowptr[row]; j < stiff.rowptr[row+1]; ++j) {
                s += stiff.values[j] * dyn.d_eta[stiff.colind[j]];
            }
            Ax[row] = s;
        }
        real_t res2 = 0, rhs2 = 0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            real_t r = dyn.ssh_rhs[n] - Ax[n];
            res2 += r * r;
            rhs2 += dyn.ssh_rhs[n] * dyn.ssh_rhs[n];
        }
        free(Ax);
        printf("[fesom_port] CG bump-RHS test (node %d): iters=%d  "
               "d_eta@probe=%.4e  max|d_eta|=%.4e  ‖res‖/‖rhs‖=%.3e\n",
               n_probe, iters, (double)deta_at_probe, (double)deta_max,
               (double)sqrt(res2 / rhs2));
    }

    /* Phase 1 step 8: momentum + rest-state preservation test (FRESH_START.md §20).
       With T=const, S=const, η=0, UV=0, no wind, no bottom drag (UV=0 → C_d*|U|*U=0):
         - compute_vel_rhs: AB history zero, no Coriolis (UV=0), no SSH grad, PGF≈0 (PoU)
                            → uv_rhs ≈ 0
         - impl_vert_visc:  no wind stress (Phase 1 forcing is zero), no bottom drag
                            (UV=0), no advection (w_i=0), TDMA on zero RHS → uv_rhs ≈ 0
         - compute_ssh_rhs: zero UV → zero divergence → ssh_rhs ≈ 0
         - solve_ssh:       zero RHS → d_eta = 0 (exact short-circuit)
         - update_vel:      no SSH grad, uv_rhs ≈ 0 → uv stays ≈ 0
         - compute_hbar:    zero divergence → hbar stays 0
       Everything must stay below machine epsilon for the literal port to be correct. */
    fesom_compute_vel_rhs(&mesh, &aux, &dyn, /* is_first_step = */ 1);
    fesom_impl_vert_visc(&mesh, &aux, &forcing, &dyn);
    {
        real_t uvr_max = 0.0;
        size_t total = (size_t)mesh.elem2D * (size_t)mesh.nl * 2;
        for (size_t i = 0; i < total; ++i) {
            real_t a = fabs(dyn.uv_rhs[i]);
            if (a > uvr_max) uvr_max = a;
        }
        printf("[fesom_port] momentum rest-state: max|uv_rhs| after vel_rhs+impl_visc = %.3e m/s "
               "(must be ≈ 0)\n", (double)uvr_max);
    }
    fesom_compute_ssh_rhs_linfs(&mesh, &dyn);
    int iters_rest = fesom_ssh_solve_cg(&stiff, &solver, &mesh, &dyn);
    {
        real_t deta_max = 0.0, ssh_max = 0.0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            if (fabs(dyn.d_eta[n])   > deta_max) deta_max = fabs(dyn.d_eta[n]);
            if (fabs(dyn.ssh_rhs[n]) > ssh_max)  ssh_max  = fabs(dyn.ssh_rhs[n]);
        }
        printf("[fesom_port] momentum rest-state: max|ssh_rhs|=%.3e  CG iters=%d  "
               "max|d_eta|=%.3e (must be ≈ 0)\n",
               (double)ssh_max, iters_rest, (double)deta_max);
    }
    fesom_update_vel(&mesh, &dyn);
    fesom_compute_hbar(&mesh, &dyn);
    {
        real_t uv_max = 0.0;
        size_t te = (size_t)mesh.elem2D * (size_t)mesh.nl * 2;
        for (size_t i = 0; i < te; ++i) {
            real_t a = fabs(dyn.uv[i]);
            if (a > uv_max) uv_max = a;
        }
        real_t hbar_max = 0.0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            real_t a = fabs(mesh.hbar[n]);
            if (a > hbar_max) hbar_max = a;
        }
        printf("[fesom_port] momentum rest-state: max|uv|=%.3e m/s  max|hbar|=%.3e m "
               "(both must be ≈ 0)\n", (double)uv_max, (double)hbar_max);
    }

    /* Phase 1 step 11: PP mixing + convective adjustment.
       Sequence: compute_vel_nodes (UV → UVnode), pp_mixing (Kv at nodes,
       Av at elements), mo_convect (raise Kv/Av where N²<0).
       With UV=0 and constant T/S: shear=0 → factor=0 → Kv = K_ver, Av = A_ver;
       no convection (bvfreq=0 → no instabmix). */
    fesom_compute_vel_nodes(&mesh, &dyn);
    fesom_pp_mixing(&mesh, &dyn, &aux);
    fesom_mo_convect(&mesh, &aux);
    {
        real_t kv_min = 1e30, kv_max = -1e30;
        real_t av_min = 1e30, av_max = -1e30;
        int kv_off_bg = 0, av_off_bg = 0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 1; nz < nzmax; ++nz) {
                real_t k = aux.Kv[FESOM_NODE3D(n, nz, mesh.nl)];
                if (k < kv_min) kv_min = k;
                if (k > kv_max) kv_max = k;
                if (fabs(k - FESOM_PHASE1_K_VER) > 1e-15) ++kv_off_bg;
            }
        }
        for (int e = 0; e < mesh.elem2D; ++e) {
            int nzmax = mesh.nlevels[e] - 1;
            for (int nz = 1; nz < nzmax; ++nz) {
                real_t a = aux.Av[FESOM_ELEM3D(e, nz, mesh.nl)];
                if (a < av_min) av_min = a;
                if (a > av_max) av_max = a;
                if (fabs(a - FESOM_PHASE1_A_VER) > 1e-15) ++av_off_bg;
            }
        }
        printf("[fesom_port] PP rest-state: Kv ∈ [%.4e, %.4e]  Av ∈ [%.4e, %.4e]  "
               "Kv≠K_ver:%d Av≠A_ver:%d\n",
               (double)kv_min, (double)kv_max, (double)av_min, (double)av_max,
               kv_off_bg, av_off_bg);
    }

    /* Phase 1 step 10: simple upwind tracer advection.
       With UV=0, W=0, constant T=10, S=35 IC: every flux is zero, so the
       tracer field must be unchanged bit-for-bit. */
    fesom_tracer_adv_scratch tra_sc;
    fesom_tracer_adv_init(&tra_sc, &mesh);
    fesom_tracer_advect_one(&tra_sc, FESOM_TRACER_T, &mesh, &dyn, &tracers);
    fesom_tracer_advect_one(&tra_sc, FESOM_TRACER_S, &mesh, &dyn, &tracers);
    {
        int t_off = 0, s_off = 0;
        real_t t_dev = 0.0, s_dev = 0.0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 0; nz < nzmax; ++nz) {
                size_t k = FESOM_NODE3D(n, nz, mesh.nl);
                real_t t_dev_local = fabs(tracers.data[FESOM_TRACER_T].values[k] - 10.0);
                real_t s_dev_local = fabs(tracers.data[FESOM_TRACER_S].values[k] - 35.0);
                if (t_dev_local > 0)        ++t_off;
                if (s_dev_local > 0)        ++s_off;
                if (t_dev_local > t_dev) t_dev = t_dev_local;
                if (s_dev_local > s_dev) s_dev = s_dev_local;
            }
        }
        printf("[fesom_port] tracer rest-state (UV=W=0, const IC): "
               "T deviations = %d (max %.3e),  S deviations = %d (max %.3e)\n",
               t_off, (double)t_dev, s_off, (double)s_dev);
    }
    /* Slice 14a sanity: compute fct_LO from the upwind fluxes and verify it
       equals the input T (=10) at machine ε. The upwind fluxes from the call
       above are still in tra_sc.adv_flux_*; re-run upwind on the current T and
       check the LO solution. */
    {
        /* Reset T to exact 10/35 in case the previous advect_one introduced
           round-off. (The fct_LO test is about the LO formula itself.) */
        for (int n = 0; n < mesh.nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 0; nz < nzmax; ++nz) {
                size_t k = FESOM_NODE3D(n, nz, mesh.nl);
                tracers.data[FESOM_TRACER_T].values[k] = 10.0;
            }
        }
        /* Re-run upwind fluxes (UV=W=0 → all fluxes 0) and compute fct_LO. */
        fesom_tracer_advect_one(&tra_sc, FESOM_TRACER_T, &mesh, &dyn, &tracers);
        fesom_tracer_compute_fct_LO(&tra_sc, FESOM_TRACER_T, &mesh, &tracers);
        real_t lo_dev = 0.0;
        int    lo_off = 0;
        for (int n = 0; n < mesh.nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 0; nz < nzmax; ++nz) {
                real_t d = fabs(tra_sc.fct_LO[FESOM_NODE3D(n, nz, mesh.nl)] - 10.0);
                if (d > 0) ++lo_off;
                if (d > lo_dev) lo_dev = d;
            }
        }
        printf("[fesom_port] FCT LO sanity (UV=W=0, T=10): "
               "fct_LO≠10 count=%d  max|fct_LO-10|=%.3e\n",
               lo_off, (double)lo_dev);
    }

    /* Phase 3 step 21: load PHC IC if a path was given. Replaces the
       constant T=10 / S=35 IC. PHC files are in-situ T (CF
       standard_name=sea_water_temperature) — t_insitu=1 triggers ptheta
       conversion. */
    if (phc_path) {
        fesom_phc_load_ic(phc_path, &mesh, &tracers, /*t_insitu=*/ 1);
    } else {
        /* Phase 2 visualisation aid: add a Gaussian +5°C T blob in the North
           Atlantic-ish region. The wind-driven flow will advect and deform it. */
        fesom_ic_tracer_T_blob(&mesh, &tracers,
                               /*lon0_deg=*/ -45.0,
                               /*lat0_deg=*/  40.0,
                               /*sigma_deg=*/ 10.0,
                               /*sigma_z=*/  300.0,
                               /*amp_C=*/      5.0);
    }

    /* Phase 3 step 24: JRA55-do daily forcing.
       If jra55_year > 0, init JRA55 reader + bulk formulae + SSS restoring +
       CORE2 runoff; the timestep loop below calls fesom_jra55_step,
       fesom_bulk_compute, and fesom_sss_runoff_step every step.
       Otherwise fall back to Phase 2 analytical wind. */
    fesom_jra55 jra; int use_jra = 0;
    fesom_sss_runoff sr; int use_sr = 0;
    if (jra55_year > 0) {
        fesom_jra55_init(&jra, &mesh);
        fesom_jra55_open_year(&jra, &mesh, jra55_year);
        use_jra = 1;
        /* Phase 3 step 25 paths from work_core/namelist.forcing. */
        const char *sss_path    = "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/PHC2_salx.nc";
        const char *runoff_path = "/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/CORE2_runoff.nc";
        fesom_sss_runoff_init(&sr, &mesh, &forcing, sss_path, runoff_path);
        use_sr = 1;
        printf("[fesom_port] SSS restoring: %s\n", sss_path);
        {
            real_t rmn = forcing.runoff[0], rmx = rmn;
            for (int k = 1; k < mesh.nod2D; ++k) {
                if (forcing.runoff[k] < rmn) rmn = forcing.runoff[k];
                if (forcing.runoff[k] > rmx) rmx = forcing.runoff[k];
            }
            printf("[fesom_port] CORE2 runoff:  %s  (range %.4e..%.4e m/s)\n",
                   runoff_path, (double)rmn, (double)rmx);
        }
        printf("[fesom_port] mesh.ocean_area = %.4e m²\n", (double)mesh.ocean_area);
        printf("[fesom_port] JRA55 init: 8 fields opened for year %d, "
               "Nlon=%d Nlat=%d Ntime=%d  cal=%s\n",
               jra55_year,
               jra.fld[0].Nlon, jra.fld[0].Nlat, jra.fld[0].Ntime,
               jra.fld[0].calendar);
        /* Compute the very-first surface state for this set of bulk inputs
           so the first timestep sees forcing immediately. uvnode is zero at
           IC (UV=0); T_oc from initial T field. */
        fesom_jra55_step(&jra, &mesh, jra55_year, /*daynew=*/1, /*timenew=*/0.0);
        fesom_compute_vel_nodes(&mesh, &dyn);
        fesom_bulk_compute(&jra, &mesh, &dyn, &tracers, &forcing);
        real_t tx_min = forcing.stress_surf[0], tx_max = tx_min;
        for (int e = 1; e < mesh.elem2D; ++e) {
            real_t t = forcing.stress_surf[2*e + 0];
            if (t < tx_min) tx_min = t;
            if (t > tx_max) tx_max = t;
        }
        real_t qmin = forcing.heat_flux[0], qmax = qmin;
        for (int n = 1; n < mesh.nod2D; ++n) {
            real_t q = forcing.heat_flux[n];
            if (q < qmin) qmin = q;
            if (q > qmax) qmax = q;
        }
        printf("[fesom_port] JRA55 first state: stress_surf_x ∈ [%+.3e, %+.3e] N/m²  "
               "heat_flux ∈ [%+.2e, %+.2e] W/m²\n",
               (double)tx_min, (double)tx_max, (double)qmin, (double)qmax);
    } else {
        /* Phase 2 step 17: analytical wind forcing (Slice A — constant zonal
           cosine pattern, no thermal). Should drive an Ekman drift on top of
           the rest state from IC. */
        fesom_forcing_set_analytical(&mesh, &forcing,
                                     /*tau0     = */ 0.05,    /* N/m² */
                                     /*Ly_factor= */ 2.0);    /* one cosine over 90° */
        real_t tx_min = forcing.stress_surf[0], tx_max = tx_min;
        for (int e = 1; e < mesh.elem2D; ++e) {
            real_t t = forcing.stress_surf[2*e + 0];
            if (t < tx_min) tx_min = t;
            if (t > tx_max) tx_max = t;
        }
        printf("[fesom_port] analytical forcing: stress_surf_x ∈ [%+.3e, %+.3e] N/m²\n",
               (double)tx_min, (double)tx_max);
    }

    /* Phase 1 step 12 / Phase 2 driver loop. */
    {
        fesom_step_ctx ctx = { .stiff = &stiff,
                               .solver = &solver,
                               .tra_sc = &tra_sc };
        const int nsteps      = (nsteps_cli > 0)     ? nsteps_cli     : 500;
        const int snap_every  = (snap_every_cli > 0) ? snap_every_cli : 25;
        const int print_every = snap_every;
        printf("[fesom_port] timestep loop: %d steps, dt=%.0f s, print every %d, snapshot every %d\n",
               nsteps, FESOM_PHASE1_DT, print_every, snap_every);
        if (out_dir) {
            char path[1024];
            snprintf(path, sizeof(path), "%s/snap_%06d.nc", out_dir, 0);
            fesom_io_write_snapshot(path, 0, FESOM_PHASE1_DT,
                                    &mesh, &dyn, &tracers, &aux);
        }
        printf("  step    iters    max|uv|       max|eta|      max|w|        "
               "min(T)    max(T)    min(S)    max(S)\n");
        /* Cumulative days at the start of each month (non-leap year, since
           the JRA55 calendar conversion above is gregorian — but for monthly
           SSS climatology we just need a 1-12 month index). */
        static const int days_in_month[12] = {
            31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
        };
        int month_prev = -1;
        for (int n = 1; n <= nsteps; ++n) {
            int    daynew  = 1;
            double timenew = 0.0;
            if (use_jra) {
                double t_sec = (double)(n - 1) * (double)FESOM_PHASE1_DT;
                daynew  = (int)floor(t_sec / 86400.0) + 1;
                timenew = t_sec - (daynew - 1) * 86400.0;
                fesom_jra55_step(&jra, &mesh, jra55_year, daynew, (real_t)timenew);
                fesom_bulk_compute(&jra, &mesh, &dyn, &tracers, &forcing);
            }
            if (use_sr) {
                /* Convert daynew → month_now (1..12), non-leap year. */
                int month_now = 1;
                int doy_remain = daynew;
                for (int mi = 0; mi < 12; ++mi) {
                    if (doy_remain <= days_in_month[mi]) { month_now = mi + 1; break; }
                    doy_remain -= days_in_month[mi];
                }
                /* Fortran update_monthly_flag fires at start (mstep==1) and at
                   midnight at end of each month. We approximate with: first
                   step OR month change. */
                int update_monthly_flag = (n == 1) || (month_now != month_prev);
                fesom_sss_runoff_step(&sr, &mesh, &tracers, &forcing,
                                      jra55_year, month_now, update_monthly_flag);
                month_prev = month_now;
            }
            int iters = fesom_timestep(n, &ctx, &mesh, &aux, &dyn,
                                       &tracers, &forcing);
            if (n == 1 || n % print_every == 0 || n == nsteps) {
                real_t uv_max = 0.0, eta_max = 0.0, w_max = 0.0;
                size_t te = (size_t)mesh.elem2D * (size_t)mesh.nl * 2;
                for (size_t i = 0; i < te; ++i) {
                    real_t a = fabs(dyn.uv[i]);
                    if (a > uv_max) uv_max = a;
                }
                for (int k = 0; k < mesh.nod2D; ++k) {
                    real_t a = fabs(dyn.eta_n[k]);
                    if (a > eta_max) eta_max = a;
                }
                size_t tn = (size_t)mesh.nod2D * (size_t)mesh.nl;
                for (size_t i = 0; i < tn; ++i) {
                    real_t a = fabs(dyn.w[i]);
                    if (a > w_max) w_max = a;
                }
                real_t T_min = 1e30, T_max = -1e30;
                real_t S_min = 1e30, S_max = -1e30;
                for (int k = 0; k < mesh.nod2D; ++k) {
                    int nzmax = mesh.nlevels_nod2D[k] - 1;
                    for (int nz = 0; nz < nzmax; ++nz) {
                        real_t T = tracers.data[FESOM_TRACER_T].values[FESOM_NODE3D(k, nz, mesh.nl)];
                        real_t S = tracers.data[FESOM_TRACER_S].values[FESOM_NODE3D(k, nz, mesh.nl)];
                        if (T < T_min) T_min = T;
                        if (T > T_max) T_max = T;
                        if (S < S_min) S_min = S;
                        if (S > S_max) S_max = S;
                    }
                }
                printf("  %4d   %5d   %.3e   %.3e   %.3e   %.5f  %.5f  %.5f  %.5f\n",
                       n, iters, (double)uv_max, (double)eta_max, (double)w_max,
                       (double)T_min, (double)T_max,
                       (double)S_min, (double)S_max);
                if (uv_max > 5.0 || !(uv_max == uv_max)) {
                    printf("[fesom_port] BLOWUP detected at step %d — aborting\n", n);
                    break;
                }
            }
            if (out_dir && (n % snap_every == 0)) {
                char path[1024];
                snprintf(path, sizeof(path), "%s/snap_%06d.nc", out_dir, n);
                fesom_io_write_snapshot(path, n, FESOM_PHASE1_DT,
                                        &mesh, &dyn, &tracers, &aux);
            }
        }
    }

    fesom_tracer_adv_free(&tra_sc);

    fesom_solverinfo_free(&solver);
    fesom_ssh_stiff_free(&stiff);

    /* Probe nodes 1001 and 1500 (Fortran 1-based) — same probes the dump shim
       uses. C 0-based ids: 1000 and 1499. With T=10°C, S=35 PSU constant,
       density at the surface should be ≈ 1027 kg/m³ → density_m_rho0 ≈ -3,
       hpressure[surface] ≈ -Z[0]*rho*g (≈ -75 Pa for Z[0]=-2.5 m), and
       bvfreq should be ~0 except for the cancel-compressibility residual. */
    {
        const int probes[] = { 1000, 1499 };  /* C 0-based */
        for (size_t k = 0; k < sizeof(probes)/sizeof(probes[0]); ++k) {
            int n = probes[k];
            if (n >= mesh.nod2D) continue;
            int nlev = mesh.nlevels_nod2D[n];
            printf("[fesom_port] EOS probe node %d (Fortran %d, nlev=%d):\n",
                   n, n + 1, nlev);
            for (int nz = 0; nz < (nlev - 1 < 5 ? nlev - 1 : 5); ++nz) {
                size_t i = FESOM_NODE3D(n, nz, mesh.nl);
                printf("  nz=%d  Z=%7.2f  ρ-ρ0=%9.5f  P=%12.4e  N²=%12.4e\n",
                       nz, (double)mesh.Z[nz],
                       (double)aux.density_m_rho0[i],
                       (double)aux.hpressure[i],
                       (double)aux.bvfreq[i]);
            }
        }
        /* Global bounds */
        real_t rmin = aux.density_m_rho0[0], rmax = rmin;
        real_t pmin = aux.hpressure[0],      pmax = pmin;
        real_t bmin = aux.bvfreq[0],         bmax = bmin;
        size_t total = (size_t)mesh.nod2D * (size_t)mesh.nl;
        for (size_t i = 1; i < total; ++i) {
            if (aux.density_m_rho0[i] < rmin) rmin = aux.density_m_rho0[i];
            if (aux.density_m_rho0[i] > rmax) rmax = aux.density_m_rho0[i];
            if (aux.hpressure[i]      < pmin) pmin = aux.hpressure[i];
            if (aux.hpressure[i]      > pmax) pmax = aux.hpressure[i];
            if (aux.bvfreq[i]         < bmin) bmin = aux.bvfreq[i];
            if (aux.bvfreq[i]         > bmax) bmax = aux.bvfreq[i];
        }
        printf("[fesom_port] EOS global: ρ-ρ0 ∈ [%g, %g]  P ∈ [%g, %g]  N² ∈ [%g, %g]\n",
               (double)rmin, (double)rmax, (double)pmin, (double)pmax,
               (double)bmin, (double)bmax);
    }

    if (use_sr)  fesom_sss_runoff_free(&sr);
    if (use_jra) fesom_jra55_free(&jra);
    fesom_forcing_free(&forcing);
    fesom_aux_free    (&aux);
    fesom_tracers_free(&tracers);
    fesom_dyn_free    (&dyn);
    fesom_mesh_free(&mesh);
    fesom_mpi_finalize(&mpi);
    return 0;
}
