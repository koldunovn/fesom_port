#include "fesom_ale.h"
#include "fesom_aux.h"
#include "fesom_bulk.h"
#include "fesom_calendar.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_eos.h"
#include "fesom_forcing.h"
#include "fesom_forcing_analytical.h"
#include "fesom_gm.h"
#include "fesom_kpp.h"
#include "fesom_ic.h"
#include "fesom_ice.h"
#include "fesom_ice_coupling.h"
#include "fesom_ice_fct.h"
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

    /* GM/Redi prereq sanity (Phase G0). */
    if (m->ulevels_nod2D_max) {
        int u_min = m->ulevels_nod2D_max[0], u_max = u_min;
        for (int i = 1; i < N_loc; ++i) {
            int v = m->ulevels_nod2D_max[i];
            if (v < u_min) u_min = v;
            if (v > u_max) u_max = v;
        }
        printf("[fesom_port] sanity: ulevels_nod2D_max min=%d max=%d\n",
               u_min, u_max);
    }
    if (m->mesh_resolution) {
        real_t r_min = m->mesh_resolution[0], r_max = r_min;
        for (int i = 1; i < N_loc; ++i) {
            real_t v = m->mesh_resolution[i];
            if (v < r_min) r_min = v;
            if (v > r_max) r_max = v;
        }
        printf("[fesom_port] sanity: mesh_resolution (Voronoi diam) min=%.1f km  max=%.1f km\n",
               (double)r_min / 1000.0, (double)r_max / 1000.0);
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
    /* Force line buffering so SLURM-redirected stdout/stderr show progress
     * immediately rather than waiting for a 4KB block. Critical for diagnosing
     * hangs. */
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

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

    /* MESH-METADATA CONSISTENCY DIAGNOSTIC.
     * For each rank that holds GID_PROBE locally (myDim or eDim), print
     * area[L,0], areasvol[L,0], and the surrounding-element count so the
     * OWNER's value can be compared to the NEIGHBOR's halo computation.
     * If they differ → mesh metadata is inconsistent across ranks (the
     * suspected source of remaining partition-boundary drift). */
    if (mpi.npes > 1) {
        /* Pick 5 GIDs from rank 0's halo (eDim) — these are guaranteed to be
         * interior on some neighbour ranks and halo on rank 0. Broadcast the
         * list, then every rank prints OWN/HALO entries for those GIDs. */
        int probe_gids[5] = { -1, -1, -1, -1, -1 };
        if (mpi.mype == 0 && mesh.eDim_nod2D >= 5) {
            int stride = mesh.eDim_nod2D / 5;
            if (stride < 1) stride = 1;
            for (int p = 0; p < 5; ++p) {
                int slot = mesh.myDim_nod2D + p * stride;
                if (slot < mesh.myDim_nod2D + mesh.eDim_nod2D)
                    probe_gids[p] = mpi.myList_nod2D[slot] - 1;
            }
        }
        MPI_Bcast(probe_gids, 5, MPI_INT, 0, mpi.MPI_COMM_FESOM);

        int N_alloc = mesh.myDim_nod2D + mesh.eDim_nod2D;
        for (int p = 0; p < 5; ++p) {
            int gid = probe_gids[p];
            if (gid < 0 || gid >= mesh.nod2D) continue;
            for (int i = 0; i < N_alloc; ++i) {
                if (mpi.myList_nod2D[i] - 1 != gid) continue;
                int is_owner = (i < mesh.myDim_nod2D);
                int o0 = mesh.nod_in_elem2D_offsets[i];
                int o1 = mesh.nod_in_elem2D_offsets[i + 1];
                fprintf(stderr,
                    "[probe rank %d] gid=%d local=%d %-4s nlev=%d "
                    "area[0]=%.10e areasvol[0]=%.10e nsurround=%d\n",
                    mpi.mype, gid, i,
                    is_owner ? "OWN" : "HALO",
                    mesh.nlevels_nod2D[i],
                    (double)mesh.area[i * mesh.nl + 0],
                    (double)mesh.areasvol[i * mesh.nl + 0],
                    o1 - o0);
                fflush(stderr);
                break;
            }
        }
    }

    print_sanity(&mesh);

    /* Phase 1 step 2: allocate the rest of the model state. All zeros until
       step 3 sets initial conditions. */
    fesom_dyn     dyn;     fesom_dyn_alloc    (&dyn,    &mesh);
    fesom_tracers tracers; fesom_tracers_alloc(&tracers, &mesh);
    fesom_aux     aux;     fesom_aux_alloc    (&aux,    &mesh);
    fesom_forcing forcing; fesom_forcing_alloc(&forcing, &mesh);
    fesom_gm      gm;      fesom_gm_alloc     (&gm,     &mesh);
    /* KPP mixing state (mirrors o_mixing_KPP_mod). Allocated unconditionally
     * like gm; used only when FESOM_MIX_SCHEME=KPP. fesom_kpp_init builds the
     * wm/ws lookup tables + Vtc/cg once (Fortran oce_mixing_kpp_init in setup);
     * harmless under PP (touches only the kpp struct). */
    fesom_kpp     kpp;     fesom_kpp_alloc    (&kpp,    &mesh, tracers.num_tracers);
    fesom_kpp_init(&kpp);
    fesom_kpp_dump_init(&kpp, mpi.mype);          /* K1 validation dump (gated) */
    fesom_kpp_dump_wscale_sweep(&kpp, mpi.mype);  /* K2 validation dump (gated) */

    /* Sea-ice state (Phase A: allocator + no-op driver only).
     * fesom_ice_setup needs the ocean dt to compute Tevp_inv, so it runs
     * after the dt CLI override has been applied. */
    fesom_ice ice;
    fesom_ice_init (&ice, &mpi, &mesh);
    fesom_ice_setup(&ice, &mpi, &mesh);

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
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 0; nz < nzmax; ++nz) {
                if (tracers.data[FESOM_TRACER_T].values[FESOM_NODE3D(n, nz, mesh.nl)] != 10.0) ++t_bad;
                if (tracers.data[FESOM_TRACER_S].values[FESOM_NODE3D(n, nz, mesh.nl)] != 35.0) ++s_bad;
            }
        }
        real_t hbar_max = 0.0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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
        for (int e = 0; e < mesh.myDim_elem2D; ++e) {
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

    /* Phase 4 (MPI): standalone sanity blocks below were designed for 1-rank
     * serial validation. On multi-rank they would require halo exchanges
     * between every kernel call (only fesom_step does that) and the parallel
     * CG (slice 30e). Each test is now gated on `mpi.npes == 1`; kernel
     * allocations and calls outside the tests still run. */
    int do_sanity = (mpi.npes == 1);
    /* Phase 1 step 4: JM-EOS + pressure_bv. */
    fesom_pressure_bv(&tracers, &mesh, &aux);
    /* Phase 1 step 5: PGF (linfs + full cells branch). */
    fesom_pressure_force_linfs_fullcell(&mesh, &aux);

    /* PGF sanity: at constant T/S, hpressure has no horizontal variation, so
       pgf should be exactly 0 (modulo machine epsilon × the partition-of-unity
       residual we measured ~1e-20). Anything materially larger means a bug. */
    if (do_sanity) {
        real_t pgmax = 0.0;
        size_t total_e = (size_t)(mesh.myDim_elem2D + mesh.eDim_elem2D + mesh.eXDim_elem2D) * (size_t)mesh.nl;
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
    fesom_ale_vert_vel_linfs(&mesh, &dyn, /* gm_on */ 0);

    /* ALE sanity: with UV ≡ 0, w must be exactly 0. hnode_new must equal hnode. */
    if (do_sanity) {
        real_t wmax = 0.0;
        size_t total = (size_t)(mesh.myDim_nod2D + mesh.eDim_nod2D) * (size_t)mesh.nl;
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
    fesom_ssh_preconditioner(&stiff, &mesh, &mpi);
    fesom_solverinfo_alloc(&solver, &mesh);
    solver.partit = &mpi;

    /* Phase E1: build the FCT mass matrix once, sharing ssh_stiff's CSR
     * sparsity. Must run AFTER fesom_ssh_stiff_alloc_and_build. */
    fesom_ice_mass_matrix_fill(&ice, &stiff, &mpi, &mesh);

    /* Stiffness diagonal sanity. */
    if (do_sanity) {
        int neg = 0;
        real_t dmin = stiff.values[stiff.rowptr[0]];
        real_t dmax = dmin;
        real_t mass_min = 1e30, mass_max = -1e30;
        for (int row = 0; row < mesh.myDim_nod2D; ++row) {
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

    /* Symmetry check. */
    if (do_sanity) {
        real_t asym_max = 0.0;
        int    asym_pairs = 0;
        for (int row = 0; row < mesh.myDim_nod2D; ++row) {
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

    /* RHS=0 test (and the bump-RHS test below): need parallel CG on multi-rank
     * which is slice 30e. Skip on multi-rank. */
    if (do_sanity) {
        fesom_compute_ssh_rhs_linfs(&mesh, &dyn);
        real_t rhs_max = 0.0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
            real_t a = fabs(dyn.ssh_rhs[n]);
            if (a > rhs_max) rhs_max = a;
        }
        int iters = fesom_ssh_solve_cg(&stiff, &solver, &mesh, &dyn);
        real_t deta_max = 0.0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
            real_t a = fabs(dyn.d_eta[n]);
            if (a > deta_max) deta_max = a;
        }
        printf("[fesom_port] CG zero-RHS test: max|ssh_rhs|=%.3e  iters=%d  "
               "max|d_eta|=%.3e (must be 0)\n",
               (double)rhs_max, iters, (double)deta_max);
    }

    /* Bump-RHS test. */
    if (do_sanity) {
        int n_probe = mesh.myDim_nod2D / 2;
        memset(dyn.ssh_rhs, 0, (size_t)(mesh.myDim_nod2D + mesh.eDim_nod2D) * sizeof(real_t));
        memset(dyn.d_eta,   0, (size_t)(mesh.myDim_nod2D + mesh.eDim_nod2D) * sizeof(real_t));
        real_t area_n = mesh.areasvol[FESOM_NODE3D(n_probe, 0, mesh.nl)];
        dyn.ssh_rhs[n_probe] = area_n;          /* unit-amplitude perturbation */

        int iters = fesom_ssh_solve_cg(&stiff, &solver, &mesh, &dyn);
        real_t deta_max = 0.0, deta_at_probe = dyn.d_eta[n_probe];
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
            real_t a = fabs(dyn.d_eta[n]);
            if (a > deta_max) deta_max = a;
        }
        /* Verify residual: ‖rhs - A x‖ / ‖rhs‖ should be ≤ soltol */
        real_t *Ax = malloc((size_t)mesh.myDim_nod2D * sizeof(real_t));
        for (int row = 0; row < mesh.myDim_nod2D; ++row) {
            real_t s = 0;
            for (int j = stiff.rowptr[row]; j < stiff.rowptr[row+1]; ++j) {
                s += stiff.values[j] * dyn.d_eta[stiff.colind[j]];
            }
            Ax[row] = s;
        }
        real_t res2 = 0, rhs2 = 0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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

    /* Phase 1 step 8: momentum + rest-state preservation test.
     * Skip on multi-rank — the kernels modify state; without the rest of
     * the timestep machinery + halo exchanges, the resulting state is
     * inconsistent and the next sanity test sees garbage. */
    if (!do_sanity) goto skip_rest_state;

    fesom_compute_vel_rhs(&mesh, &aux, &dyn, /* is_first_step = */ 1, &mpi);
    fesom_impl_vert_visc(&mesh, &aux, &forcing, &dyn);
    {
        real_t uvr_max = 0.0;
        size_t total = (size_t)(mesh.myDim_elem2D + mesh.eDim_elem2D + mesh.eXDim_elem2D) * (size_t)mesh.nl * 2;
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
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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
        size_t te = (size_t)(mesh.myDim_elem2D + mesh.eDim_elem2D + mesh.eXDim_elem2D) * (size_t)mesh.nl * 2;
        for (size_t i = 0; i < te; ++i) {
            real_t a = fabs(dyn.uv[i]);
            if (a > uv_max) uv_max = a;
        }
        real_t hbar_max = 0.0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
            real_t a = fabs(mesh.hbar[n]);
            if (a > hbar_max) hbar_max = a;
        }
        printf("[fesom_port] momentum rest-state: max|uv|=%.3e m/s  max|hbar|=%.3e m "
               "(both must be ≈ 0)\n", (double)uv_max, (double)hbar_max);
    }

    /* Phase 1 step 11: PP mixing rest-state test. Same caveat — skipped
     * for multi-rank since state needs halo-exchanges. */
    fesom_compute_vel_nodes(&mesh, &dyn);
    fesom_pp_mixing(&mesh, &dyn, &aux);
    fesom_mo_convect(&mesh, &aux);
    {
        real_t kv_min = 1e30, kv_max = -1e30;
        real_t av_min = 1e30, av_max = -1e30;
        int kv_off_bg = 0, av_off_bg = 0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
            int nzmax = mesh.nlevels_nod2D[n] - 1;
            for (int nz = 1; nz < nzmax; ++nz) {
                real_t k = aux.Kv[FESOM_NODE3D(n, nz, mesh.nl)];
                if (k < kv_min) kv_min = k;
                if (k > kv_max) kv_max = k;
                if (fabs(k - FESOM_PHASE1_K_VER) > 1e-15) ++kv_off_bg;
            }
        }
        for (int e = 0; e < mesh.myDim_elem2D; ++e) {
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

skip_rest_state:
    /* tra_sc allocation runs unconditionally — needed by timestep loop. */
    fesom_tracer_adv_scratch tra_sc;
    fesom_tracer_adv_init(&tra_sc, &mesh);
    tra_sc.partit = &mpi;       /* used by FCT mid-pipeline exchange */
    if (do_sanity) {
        fesom_tracer_advect_one(&tra_sc, FESOM_TRACER_T, &mesh, &dyn, &tracers);
        fesom_tracer_advect_one(&tra_sc, FESOM_TRACER_S, &mesh, &dyn, &tracers);
    }
    if (do_sanity) {
        int t_off = 0, s_off = 0;
        real_t t_dev = 0.0, s_dev = 0.0;
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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
    /* Slice 14a sanity: fct_LO test. */
    if (do_sanity) {
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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
        for (int n = 0; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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
        fesom_phc_load_ic(phc_path, &mesh, &tracers, /*t_insitu=*/ 1, &mpi);
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

    /* Sea-ice cold-start IC: must run AFTER tracer IC so SST is set.
     * Mirrors Fortran ice_initial_state at ice_setup_step.F90:500-521. */
    fesom_ice_initial_state(&ice, &tracers, &mpi, &mesh);

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
            for (int k = 1; k < mesh.myDim_nod2D + mesh.eDim_nod2D; ++k) {
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
        fesom_jra55_step(&jra, &mesh, &mpi, jra55_year, /*daynew=*/1, /*timenew=*/0.0);
        fesom_compute_vel_nodes(&mesh, &dyn);
        fesom_bulk_compute(&jra, &mesh, &dyn, &tracers, &forcing, &ice, &mpi);
        real_t tx_min = forcing.stress_surf[0], tx_max = tx_min;
        for (int e = 1; e < mesh.myDim_elem2D; ++e) {
            real_t t = forcing.stress_surf[2*e + 0];
            if (t < tx_min) tx_min = t;
            if (t > tx_max) tx_max = t;
        }
        real_t qmin = forcing.heat_flux[0], qmax = qmin;
        for (int n = 1; n < mesh.myDim_nod2D + mesh.eDim_nod2D; ++n) {
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
        /* The node interpolation in set_analytical area-averages over
           nod_in_elem2D, whose neighbour list is INCOMPLETE for eDim halo
           nodes -> divergent eDim values. Exchange so eDim matches the owner,
           mirroring the bulk path (fesom_bulk.c) and Fortran's
           exchange_nod(stress_atmoce). Required for SSH-RHS conservation:
           oce_fluxes_mom re-averages these node values onto (replicated)
           elements every step. */
        fesom_halo_exchange(forcing.stress_node_surf, FESOM_HALO_NOD2D, 1, 2, &mpi);
        real_t tx_min = forcing.stress_surf[0], tx_max = tx_min;
        for (int e = 1; e < mesh.myDim_elem2D; ++e) {
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
                               .tra_sc = &tra_sc,
                               .partit = &mpi,
                               .ice    = &ice,
                               .gm     = &gm,
                               .kpp    = &kpp,
                               .jra    = use_jra ? &jra : NULL,
                               .sr     = use_sr  ? &sr  : NULL };
        const int nsteps      = (nsteps_cli > 0)     ? nsteps_cli     : 500;
        /* snap_every semantics:
         *   > 0  → snapshot every N steps
         *   = 0  → use default (25)              [legacy default]
         *   < 0  → DISABLE snapshots entirely    [for production runs] */
        const int snap_every  = (snap_every_cli > 0) ? snap_every_cli
                              : (snap_every_cli == 0 ? 25 : 0);
        /* FESOM_PRINT_EVERY env override — useful for tracking instability
         * onset without changing snapshot cadence. Defaults to snap_every. */
        const char *pe_env = getenv("FESOM_PRINT_EVERY");
        /* Default print_every to snap_every if positive, else a safe
         * non-zero fallback (1000) — n % print_every must never divide by 0. */
        const int print_every = (pe_env && atoi(pe_env) > 0) ? atoi(pe_env)
                              : (snap_every > 0 ? snap_every : 1000);

        /* Physics-bisect knobs (for MPI drift hunt — developer hint: one of
         * wind / fluxes / tracer update is amplifying partition-dependent
         * noise). Each knob is applied every step after bulk+runoff and
         * before fesom_timestep. */
        const int freeze_ts = (getenv("FESOM_FREEZE_TS") && atoi(getenv("FESOM_FREEZE_TS")));
        const int no_wind   = (getenv("FESOM_NO_WIND")   && atoi(getenv("FESOM_NO_WIND")));
        const int no_hflux  = (getenv("FESOM_NO_HFLUX")  && atoi(getenv("FESOM_NO_HFLUX")));
        /* Snapshot T,S at IC for FESOM_FREEZE_TS restore. Covers myDim+eDim
         * so halo copies are also reset — the step writes both. */
        real_t *T_ic = NULL, *S_ic = NULL;
        size_t ts_nbytes = 0;
        if (freeze_ts) {
            ts_nbytes = (size_t)(mesh.myDim_nod2D + mesh.eDim_nod2D)
                      * (size_t)mesh.nl * sizeof(real_t);
            T_ic = malloc(ts_nbytes); S_ic = malloc(ts_nbytes);
            FESOM_CHECK(T_ic && S_ic, "FESOM_FREEZE_TS: alloc snapshot");
            memcpy(T_ic, tracers.data[FESOM_TRACER_T].values, ts_nbytes);
            memcpy(S_ic, tracers.data[FESOM_TRACER_S].values, ts_nbytes);
        }
        if (mpi.mype == 0 && (freeze_ts || no_wind || no_hflux)) {
            printf("[fesom_port] physics-bisect knobs: FREEZE_TS=%d NO_WIND=%d NO_HFLUX=%d\n",
                   freeze_ts, no_wind, no_hflux);
        }

        printf("[fesom_port] timestep loop: %d steps, dt=%.0f s, print every %d, snapshot every %d\n",
               nsteps, FESOM_PHASE1_DT, print_every, snap_every);

        /* Task 6: orchestrator owns the calendar + the monthly streams. Init
         * BEFORE the step-0 snapshot so that snapshot's `time` attribute
         * can be anchored to the model calendar (Task 9). */
        fesom_io_t io;
        {
            /* When JRA55 is on (typical), the model calendar must start
             * at jra55_year so that fesom_jra55_step_cal indexes the
             * forcing files correctly. With JRA55 disabled (analytical-
             * forcing smoke runs only) the calendar still needs a start
             * date for CF time stamps; pick 1958 to match the project's
             * de-facto default and avoid confusion when files written by
             * an analytical run sit alongside JRA55-driven output. */
            int start_year = (jra55_year > 0) ? jra55_year : 1958;
            const char *io_out = (out_dir != NULL) ? out_dir : ".";
            /* Optional override file: FESOM_IO_CONFIG=<path>. Empty / unset
             * → use compiled-in monthly-only default. */
            const char *cfg_path = getenv("FESOM_IO_CONFIG");
            fesom_io_init(&io, io_out, FESOM_CAL_GREGORIAN,
                          start_year, 1, 1,
                          cfg_path,
                          &mesh, &mpi);
        }

        if (out_dir && snap_every > 0) {
            char path[1024];
            snprintf(path, sizeof(path), "%s/snap_%06d.nc", out_dir, 0);
            fesom_io_write_snapshot(path, 0, FESOM_PHASE1_DT, &io.calendar,
                                    &mesh, &dyn, &tracers, &aux, &ice, &mpi);
        }
        if (mpi.mype == 0) {
            printf("# step iters | uv eta w | T-range S-range | stress hf wf vs rs | hp pgf rho bv-range Kv Av\n");
            fflush(stdout);
        }

        for (int n = 1; n <= nsteps; ++n) {
            if (use_jra) {
                /* Calendar crosses year boundaries cleanly now that
                 * fesom_jra55_step_cal calls open_year on rollover; the
                 * old "cal.year == jra55_year" safety overlap from Task 2
                 * was retired with multi-year support. */
                fesom_jra55_step_cal(&jra, &mesh, &mpi, &io.calendar);
                fesom_bulk_compute(&jra, &mesh, &dyn, &tracers, &forcing, &ice, &mpi);
            }
            if (use_sr) {
                fesom_sss_runoff_step_cal(&sr, &mesh, &tracers, &forcing, &mpi,
                                          n, &io.prev_calendar, &io.calendar);
            }
            /* Physics-bisect toggles (env-gated). Applied AFTER bulk+runoff
             * writes forcing, BEFORE fesom_timestep reads it. */
            if (no_wind) {
                int E = mesh.myDim_elem2D + mesh.eDim_elem2D + mesh.eXDim_elem2D;
                int N = mesh.myDim_nod2D + mesh.eDim_nod2D;
                memset(forcing.stress_surf,      0, (size_t)E * 2 * sizeof(real_t));
                memset(forcing.stress_node_surf, 0, (size_t)N * 2 * sizeof(real_t));
                memset(ice.stress_atmice_x,      0, (size_t)N * sizeof(real_t));
                memset(ice.stress_atmice_y,      0, (size_t)N * sizeof(real_t));
            }
            if (no_hflux) {
                int N = mesh.myDim_nod2D + mesh.eDim_nod2D;
                memset(forcing.heat_flux,    0, (size_t)N * sizeof(real_t));
                memset(forcing.water_flux,   0, (size_t)N * sizeof(real_t));
                memset(forcing.virtual_salt, 0, (size_t)N * sizeof(real_t));
                memset(forcing.relax_salt,   0, (size_t)N * sizeof(real_t));
            }
            if (freeze_ts) {
                memcpy(tracers.data[FESOM_TRACER_T].values, T_ic, ts_nbytes);
                memcpy(tracers.data[FESOM_TRACER_S].values, S_ic, ts_nbytes);
            }

            /* Sea-ice step runs BEFORE the ocean step within each iteration,
             * mirroring the Fortran flow in oce_timestep_ale. The ice step:
             *   1. ocean2ice — copies ocean surface state into ice->srfoce_*
             *   2. thermodynamics — writes ice tracers + flx_h + flx_fw
             *   3. cut_off — clamps ice tracers
             *   4. oce_fluxes — overwrites heat_flux/water_flux/virtual_salt/relax_salt
             *      from the ice-mediated fluxes (Phase C2).
             * The ocean step that follows then sees the updated forcing. */
            fesom_ice_step(n, &ice, &mpi, &mesh,
                           &dyn, &tracers, &forcing,
                           use_jra ? &jra : NULL,
                           use_sr  ? &sr  : NULL,
                           &stiff);

            /* Phase D5 — ice-mediated surface momentum flux. Must be after
             * fesom_ice_step (uses post-EVP ice->uice/vice halo) and before
             * fesom_timestep (compute_vel_rhs / impl_vert_visc consume
             * forcing->stress_surf). Mirrors Fortran oce_fluxes_mom called
             * from oce_timestep_ale before the ocean step.
             *
             * Skipped under FESOM_NO_WIND: that toggle's intent is "no
             * surface stress on the ocean", which D5 would otherwise
             * re-populate from ice-ocean drag. */
            if (!no_wind) {
                fesom_ice_oce_fluxes_mom(&ice, &mpi, &mesh, &forcing);
            }

            /* Shortwave penetration (oce_shortwave_pene.F90 cal_shortwave_rad,
             * called from ice_oce_coupling.F90:833): AFTER fesom_ice_step finalizes
             * heat_flux, BEFORE fesom_timestep's tracer diffusion consumes
             * heat_flux + sw_3d. chl: constant or Sweeney monthly climatology read
             * via fesom_read_other_NetCDF (gen_surface_forcing.F90:1585). Needs
             * jra->shortwave, so gated on use_jra. */
            if (use_jra && FESOM_PHASE1_USE_SW_PENE) {
                static int chl_sweeney = -1, chl_const_done = 0;
                static const char *chl_file = NULL;
                if (chl_sweeney < 0) {
                    const char *src = getenv("FESOM_CHL_SOURCE");   /* default Sweeney */
                    chl_sweeney = (src && strcmp(src, "None") == 0) ? 0 : 1;
                    chl_file = getenv("FESOM_CHL_FILE");
                    if (!chl_file)
                        chl_file = "/pool/data/AWICM/FESOM2/FORCING/Sweeney/Sweeney_2005.nc";
                }
                int Nn = mesh.myDim_nod2D + mesh.eDim_nod2D;
                if (!chl_sweeney) {
                    if (!chl_const_done) {
                        for (int i = 0; i < Nn; ++i)
                            forcing.chl[i] = (real_t)FESOM_PHASE1_CHL_CONST;
                        chl_const_done = 1;
                    }
                } else if (n == 1 ||
                           fesom_calendar_crossed(&io.prev_calendar, &io.calendar,
                                                  FESOM_PERIOD_MONTHLY)) {
                    /* month_now is already advanced (calendar_crossed fires on the
                     * first step of the new month), so NO +1 — unlike the Fortran
                     * which fires at end-of-month and adds +1. Adding +1 here skips a
                     * month (was loading 1,3,4,..,12, never February). See the SSS
                     * read in fesom_sss_runoff.c for the full explanation. */
                    int mi = io.calendar.month;
                    if (mpi.mype == 0)
                        fprintf(stderr, "[chl] Updating chlorophyll climatology for month %d\n", mi);
                    fesom_read_other_NetCDF(chl_file, "chl", mi, forcing.chl, 1, 1, &mesh);
                }
                fesom_cal_shortwave_rad(&mesh, &jra, &ice, &forcing);
            }

            /* Per-step heartbeat on rank 0 — independent of print_every so
             * we always see SOMETHING happening even if the model is hung
             * inside a single step. */
            if (mpi.mype == 0) {
                fprintf(stderr, "[step %d] entering fesom_timestep\n", n);
            }
            int iters = fesom_timestep(n, &ctx, &mesh, &aux, &dyn,
                                       &tracers, &forcing);
            if (mpi.mype == 0) {
                fprintf(stderr, "[step %d] done — %d CG iters\n", n, iters);
            }

            /* Drive output streams + advance calendar. Per-step cost is
             * a local-only sum into accumulators (no MPI on hot path);
             * gather + write happens only on calendar boundaries. */
            {
                fesom_state st = { &mesh, &dyn, &tracers, &aux, &ice, &forcing };
                fesom_io_step(&io, (double)FESOM_PHASE1_DT, &st, &mpi);
            }
            if (n == 1 || n % print_every == 0 || n == nsteps) {
                /* Iterate myDim only for stat collection — halo entries are
                 * just owner copies, so global max via MPI_Allreduce is
                 * exact. Reduces noise from divergent halos before exchange. */
                real_t uv_max = 0.0, eta_max = 0.0, w_max = 0.0;
                size_t te = (size_t)mesh.myDim_elem2D * (size_t)mesh.nl * 2;
                for (size_t i = 0; i < te; ++i) {
                    real_t a = fabs(dyn.uv[i]);
                    if (a > uv_max) uv_max = a;
                }
                for (int k = 0; k < mesh.myDim_nod2D; ++k) {
                    real_t a = fabs(dyn.eta_n[k]);
                    if (a > eta_max) eta_max = a;
                }
                size_t tn = (size_t)mesh.myDim_nod2D * (size_t)mesh.nl;
                for (size_t i = 0; i < tn; ++i) {
                    real_t a = fabs(dyn.w[i]);
                    if (a > w_max) w_max = a;
                }
                real_t T_min = 1e30, T_max = -1e30;
                real_t S_min = 1e30, S_max = -1e30;
                for (int k = 0; k < mesh.myDim_nod2D; ++k) {
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
                /* Forcing fields (myDim only). */
                real_t stress_max = 0.0, hflux_max = 0.0, wflux_max = 0.0;
                real_t vsalt_max = 0.0, rsalt_max = 0.0;
                for (int k = 0; k < mesh.myDim_elem2D; ++k) {
                    real_t s = fabs(forcing.stress_surf[2*k+0]);
                    real_t t = fabs(forcing.stress_surf[2*k+1]);
                    if (s > stress_max) stress_max = s;
                    if (t > stress_max) stress_max = t;
                }
                for (int k = 0; k < mesh.myDim_nod2D; ++k) {
                    real_t a = fabs(forcing.heat_flux[k]);
                    if (a > hflux_max) hflux_max = a;
                    a = fabs(forcing.water_flux[k]);
                    if (a > wflux_max) wflux_max = a;
                    if (forcing.virtual_salt) {
                        a = fabs(forcing.virtual_salt[k]);
                        if (a > vsalt_max) vsalt_max = a;
                    }
                    if (forcing.relax_salt) {
                        a = fabs(forcing.relax_salt[k]);
                        if (a > rsalt_max) rsalt_max = a;
                    }
                }
                /* Aux fields (myDim only). */
                real_t hpres_max = 0.0, pgf_max = 0.0, dens_max = 0.0;
                real_t bv_min = 1e30, bv_max = -1e30;
                real_t Kv_max = 0.0, Av_max = 0.0;
                for (int k = 0; k < mesh.myDim_nod2D; ++k) {
                    int nzmax = mesh.nlevels_nod2D[k] - 1;
                    for (int nz = 0; nz < nzmax; ++nz) {
                        size_t i = FESOM_NODE3D(k, nz, mesh.nl);
                        real_t a = fabs(aux.hpressure[i]);
                        if (a > hpres_max) hpres_max = a;
                        a = fabs(aux.density_m_rho0[i]);
                        if (a > dens_max) dens_max = a;
                        a = aux.bvfreq[i];
                        if (a < bv_min) bv_min = a;
                        if (a > bv_max) bv_max = a;
                        a = fabs(aux.Kv[i]);
                        if (a > Kv_max) Kv_max = a;
                    }
                }
                for (int e = 0; e < mesh.myDim_elem2D; ++e) {
                    int nzmax = mesh.nlevels[e] - 1;
                    for (int nz = 0; nz < nzmax; ++nz) {
                        size_t i = FESOM_ELEM3D(e, nz, mesh.nl);
                        real_t a = fabs(aux.pgf_x[i]);
                        if (a > pgf_max) pgf_max = a;
                        a = fabs(aux.pgf_y[i]);
                        if (a > pgf_max) pgf_max = a;
                        a = fabs(aux.Av[i]);
                        if (a > Av_max) Av_max = a;
                    }
                }
                /* Global reductions across ranks. */
                if (mpi.npes > 1) {
                    real_t buf_max[16] = {uv_max, eta_max, w_max, T_max, S_max,
                        stress_max, hflux_max, wflux_max, vsalt_max, rsalt_max,
                        hpres_max, pgf_max, dens_max, bv_max, Kv_max, Av_max};
                    real_t buf_min[3]  = {T_min, S_min, bv_min};
                    real_t out_max[16], out_min[3];
                    MPI_Allreduce(buf_max, out_max, 16, MPI_DOUBLE, MPI_MAX, mpi.MPI_COMM_FESOM);
                    MPI_Allreduce(buf_min, out_min, 3,  MPI_DOUBLE, MPI_MIN, mpi.MPI_COMM_FESOM);
                    uv_max=out_max[0]; eta_max=out_max[1]; w_max=out_max[2];
                    T_max=out_max[3];  S_max=out_max[4];
                    stress_max=out_max[5]; hflux_max=out_max[6]; wflux_max=out_max[7];
                    vsalt_max=out_max[8];  rsalt_max=out_max[9];
                    hpres_max=out_max[10]; pgf_max=out_max[11]; dens_max=out_max[12];
                    bv_max=out_max[13]; Kv_max=out_max[14]; Av_max=out_max[15];
                    T_min=out_min[0]; S_min=out_min[1]; bv_min=out_min[2];
                }
                if (mpi.mype == 0) {
                    printf("  %4d  it=%3d  uv=%.2e eta=%.2e w=%.2e | T[%.2f,%.2f] S[%.2f,%.2f] | "
                           "stress=%.2e hf=%.2e wf=%.2e vs=%.2e rs=%.2e | "
                           "hp=%.2e pgf=%.2e rho=%.2e bv[%.1e,%.1e] Kv=%.2e Av=%.2e\n",
                           n, iters, (double)uv_max, (double)eta_max, (double)w_max,
                           (double)T_min, (double)T_max, (double)S_min, (double)S_max,
                           (double)stress_max, (double)hflux_max, (double)wflux_max,
                           (double)vsalt_max, (double)rsalt_max,
                           (double)hpres_max, (double)pgf_max, (double)dens_max,
                           (double)bv_min, (double)bv_max, (double)Kv_max, (double)Av_max);
                    fflush(stdout);
                }
                if (uv_max > 5.0 || !(uv_max == uv_max)) {
                    if (mpi.mype == 0) {
                        fprintf(stderr, "[fesom_port] BLOWUP at step %d "
                                "(uv=%.3e eta=%.3e w=%.3e) — aborting all ranks\n",
                                n, (double)uv_max, (double)eta_max, (double)w_max);
                        fflush(stderr);
                        fflush(stdout);
                    }
                    MPI_Abort(mpi.MPI_COMM_FESOM, 99);
                }
            }
            if (out_dir && snap_every > 0 && (n % snap_every == 0)) {
                char path[1024];
                snprintf(path, sizeof(path), "%s/snap_%06d.nc", out_dir, n);
                fesom_io_write_snapshot(path, n, FESOM_PHASE1_DT, &io.calendar,
                                        &mesh, &dyn, &tracers, &aux, &ice, &mpi);
            }
        }
        free(T_ic); free(S_ic);

        /* Mid-window flush + close monthly stream files. */
        fesom_io_finalize(&io, &mpi);
    }

    fesom_tracer_adv_free(&tra_sc);

    fesom_solverinfo_free(&solver);
    fesom_ssh_stiff_free(&stiff);

    if (use_sr)  fesom_sss_runoff_free(&sr);
    if (use_jra) fesom_jra55_free(&jra);
    fesom_forcing_free(&forcing);
    fesom_ice_free    (&ice);
    fesom_gm_free     (&gm);
    fesom_kpp_free    (&kpp);
    fesom_aux_free    (&aux);
    fesom_tracers_free(&tracers);
    fesom_dyn_free    (&dyn);
    fesom_mesh_free(&mesh);
    fesom_mpi_finalize(&mpi);
    return 0;
}
