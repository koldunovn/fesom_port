#ifndef FESOM_MESH_H
#define FESOM_MESH_H

#include "fesom_types.h"

/*
 * Mirror of the Fortran t_mesh subset used in Phase 1. Naming follows
 * FRESH_START.md §7. All counts and indices are 0-based in the C
 * representation; on-disk Fortran 1-based IDs are converted on read.
 *
 * Halo-aware allocation from day 1: every node-sized array is sized
 * (myDim_nod2D + eDim_nod2D); same for elements. In serial mode
 * eDim_* = 0 and myDim_* equals the global count.
 */
typedef struct fesom_mesh {
    /* counts */
    int nod2D;        /* total nodes on this rank: myDim + eDim */
    int elem2D;       /* total elements on this rank             */
    int edge2D;       /* total edges on this rank                */
    int nl;           /* max number of vertical levels           */

    int myDim_nod2D,  eDim_nod2D;
    int myDim_elem2D, eDim_elem2D, eXDim_elem2D;
    int myDim_edge2D, eDim_edge2D;

    /* horizontal topology — coordinates in ROTATED RADIANS after read */
    real_t *coord_nod2D;     /* [nod2D * 2]   (lon, lat) rotated radians   */
    real_t *geo_coord_nod2D; /* [nod2D * 2]   (lon, lat) geographic radians;
                                              kept for IC/forcing interpolation */
    int    *coast_flag;      /* [nod2D]       0=interior, 1=coast          */
    int    *elem_nodes;      /* [elem2D * 3]  node IDs (0-based, CW after orient fix) */

    /* edges */
    int    *edges;         /* [edge2D * 2] node IDs (0-based)                */
    int    *edge_tri;      /* [edge2D * 2] adjacent element IDs (-1 = either
                              true boundary OR halo elem not in local map).
                              For boundary checks use partit->myList_edge2D
                              compared to mesh->edge2D_in (Fortran convention)
                              — relying on edge_tri[2nd] < 0 misclassifies any
                              interior edge whose halo element is missing. */
    int     edge2D_in;     /* global count of interior edges. Mesh files put
                              interior edges at IDs [1, edge2D_in] and boundary
                              edges at IDs [edge2D_in+1, edge2D]. Mirror of
                              Fortran mesh%edge2D_in. */

    /* vertical */
    int    *nlevels_nod2D;     /* [nod2D]   K_v⁺ (MAX over surrounding cells, from nlvls.out) */
    int    *nlevels_nod2D_min; /* [nod2D]   K_v⁻ (MIN over surrounding cells) — ALE limit */
    int    *ulevels_nod2D;     /* [nod2D]   upper level (=1 without cavities) */
    int    *ulevels_nod2D_max; /* [nod2D]   MAX of ulevels over surrounding cells.
                                  Used by GM fer_solve_Gamma. Mirror of
                                  Fortran mesh%ulevels_nod2D_max, oce_mesh.F90:1657. */
    int    *nlevels;           /* [elem2D]  K_c, per-cell level count, from elvls.out */
    int    *ulevels;           /* [elem2D]  upper level per cell (=1 without cavities) */
    real_t *zbar;              /* [nl]      interface depths, negative downward */
    real_t *Z;                 /* [nl-1]    mid-layer depths = 0.5*(zbar[nz]+zbar[nz+1]) */
    real_t *zbar_3d_n;         /* [nod2D * nl]  per-node interface depths.
                                  In our linfs/no-cavity/no-partial-cell config this
                                  collapses to zbar[nz] for valid levels and 0 elsewhere
                                  (constant in time). Mirror of Fortran mesh%zbar_3d_n,
                                  oce_ale.F90:266+340. Used by GM scaling_GMzexp. */
    real_t *depth;             /* [nod2D]   bathymetry per node — INPUT METADATA only;
                                  the operative cellwise depth lives in nlevels(elem) */
    real_t *mesh_resolution;   /* [nod2D]   Voronoi-cell diameter, smoothed 3 passes.
                                  Mirror of Fortran mesh%mesh_resolution,
                                  oce_mesh.F90:2199+2353-2378. Used by GM/Redi. */

    /* node→element inverse (CSR). For node n, its surrounding cells are
       nod_in_elem2D[nod_in_elem2D_offsets[n] .. nod_in_elem2D_offsets[n+1]-1]. */
    int    *nod_in_elem2D_offsets; /* [nod2D + 1] */
    int    *nod_in_elem2D;         /* [Σ counts ] */

    /* derived geometry */
    real_t *elem_area;        /* [elem2D]            cell area, m² */
    real_t *area;             /* [nod2D * nl]        upper-edge scalar CV area, m² */
    real_t *areasvol;         /* [nod2D * nl]        "mid" CV area (= area without cavity), m² */
    real_t  ocean_area;       /* total open-ocean surface area Σ areasvol(ulevels_nod2D(n), n)
                                  for ulevels_nod2D(n)==1, m² (Fortran oce_mesh.F90:2380-2393) */

    /* Per-cell metric */
    real_t *elem_cos;         /* [elem2D]   cos(rotated lat at centroid) */
    real_t *metric_factor;    /* [elem2D]   tan(rotated lat)/R_earth, m⁻¹ */
    real_t *coriolis;         /* [elem2D]   2Ω sin(geographic lat at centroid), s⁻¹ */
    real_t *coriolis_node;    /* [nod2D]    2Ω sin(geographic lat), s⁻¹ */
    real_t *elem_center_x;    /* [elem2D]   centroid x in rotated radians (cyclic-aware) */
    real_t *elem_center_y;    /* [elem2D]   centroid y in rotated radians */

    /* Edge geometry */
    real_t *edge_dxdy;        /* [edge2D * 2]  node2 − node1 in rotated radians */
    /* MFCT up/down-wind triangle per edge (oce_muscl_adv.F90:185-333). For each
       owned edge: [2*edge+0] = upwind triangle (around node1, contains −x),
       [2*edge+1] = downwind triangle (around node2, contains +x); -1 = absent
       (boundary / unresolved → fill_up_dn_grad falls back). 0-based elem ids. */
    int    *edge_up_dn_tri;   /* [myDim_edge2D * 2] */
    real_t *edge_cross_dxdy;  /* [edge2D * 4]  (cell1_center − edge_mid, cell2_center − edge_mid)
                                                in physical METERS, lon scaled by elem_cos.
                                                Components 2,3 are zero on boundary edges. */

    /* Per-cell linear-shape-function gradient (∂N_i/∂x and ∂N_i/∂y for the
       3 vertices, in 1/m). Layout:
         gradient_sca[0..2] = dN1/dx, dN2/dx, dN3/dx
         gradient_sca[3..5] = dN1/dy, dN2/dy, dN3/dy
       Mirror of mesh%gradient_sca(6, elem) in oce_mesh.F90. */
    real_t *gradient_sca;     /* [elem2D * 6] */

    /* Time-evolving layer thicknesses & SSH (mesh-owned in Fortran t_mesh).
       Allocated zero by fesom_mesh_alloc_state; initialised by Phase 1 step 3. */
    real_t *hnode;            /* [nod2D * nl]   layer thickness at vertices  */
    real_t *hnode_new;        /* [nod2D * nl]   thickness after SSH update   */
    real_t *helem;            /* [elem2D * nl]  layer thickness at cells     */
    real_t *hbar;             /* [nod2D]                                     */
    real_t *hbar_old;         /* [nod2D]                                     */

    /* Boundary mask: 1 on interior nodes, 0 on open-boundary edge endpoints.
       Allocated and populated by fesom_ice_init (mirrors MOD_ICE.F90:889-895
       which moved the alloc into ice_init for namelist-timing reasons).
       Read by ice EVP. NULL until fesom_ice_init runs. */
    real_t *bc_index_nod2D;   /* [myDim+eDim] */
} fesom_mesh;

void fesom_mesh_init(fesom_mesh *m);
void fesom_mesh_free(fesom_mesh *m);

/* Parallel mesh reader. Rank 0 reads global files, broadcasts to all ranks;
 * each rank then extracts its slice using partit->myList_* indices. After
 * this call, all per-rank arrays are sized to local extent (see size table
 * in the .c file). For npes==1 the reader behaves bit-identically to the
 * pre-MPI serial path. */
struct fesom_partit;
void fesom_mesh_read(fesom_mesh *m, const char *mesh_dir,
                     struct fesom_partit *partit);

void fesom_mesh_compute_metrics(fesom_mesh *m, struct fesom_partit *partit);

/* MFCT: build edge_up_dn_tri (per-edge up/down-wind triangle). Called at the end
 * of fesom_mesh_compute_metrics (needs edges, nod_in_elem2D, coord_nod2D + the
 * wide element halo). oce_muscl_adv.F90:185-333. */
void fesom_mesh_build_edge_up_dn_tri(fesom_mesh *m, struct fesom_partit *partit);

/* Allocate time-evolving mesh state (hnode, hnode_new, helem, hbar, hbar_old).
   Counts must already be set (call after fesom_mesh_read). All zeros. */
void fesom_mesh_alloc_state(fesom_mesh *m);

/* Rotate a geographic (east,north) vector into the model's rotated frame at a
   node, given its geographic and rotated lon/lat [rad]. Mirror of Fortran
   vector_g2r; no-op if FESOM_FORCE_ROTATION is 0. Used to rotate JRA winds. */
void fesom_vector_g2r(real_t *u, real_t *v,
                      real_t glon, real_t glat, real_t rlon, real_t rlat);

#endif /* FESOM_MESH_H */
