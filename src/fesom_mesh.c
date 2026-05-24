#include "fesom_mesh.h"
#include "fesom_constants.h"
#include "fesom_halo.h"
#include "fesom_partit.h"

#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MPI_CHECK(call) do {                                              \
    int _ec = (call);                                                     \
    if (_ec != MPI_SUCCESS) {                                             \
        char _s[MPI_MAX_ERROR_STRING]; int _n=0;                          \
        MPI_Error_string(_ec, _s, &_n);                                   \
        FESOM_DIE("MPI error: %s", _s);                                   \
    }                                                                     \
} while (0)

void fesom_mesh_init(fesom_mesh *m)
{
    memset(m, 0, sizeof(*m));
}

void fesom_mesh_free(fesom_mesh *m)
{
    free(m->coord_nod2D);
    free(m->geo_coord_nod2D);
    free(m->coast_flag);
    free(m->elem_nodes);
    free(m->edges);
    free(m->edge_tri);
    free(m->nlevels_nod2D);
    free(m->nlevels_nod2D_min);
    free(m->ulevels_nod2D);
    free(m->ulevels_nod2D_max);
    free(m->nlevels);
    free(m->ulevels);
    free(m->zbar);
    free(m->Z);
    free(m->zbar_3d_n);
    free(m->depth);
    free(m->mesh_resolution);
    free(m->nod_in_elem2D_offsets);
    free(m->nod_in_elem2D);
    free(m->edge_up_dn_tri);
    free(m->elem_area);
    free(m->area);
    free(m->areasvol);
    free(m->elem_cos);
    free(m->metric_factor);
    free(m->coriolis);
    free(m->coriolis_node);
    free(m->elem_center_x);
    free(m->elem_center_y);
    free(m->edge_dxdy);
    free(m->edge_cross_dxdy);
    free(m->gradient_sca);
    free(m->hnode);
    free(m->hnode_new);
    free(m->helem);
    free(m->hbar);
    free(m->hbar_old);
    free(m->bc_index_nod2D);
    memset(m, 0, sizeof(*m));
}

void fesom_mesh_alloc_state(fesom_mesh *m)
{
    int N = m->myDim_nod2D + m->eDim_nod2D;
    int E = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    size_t n_nl = (size_t)N * (size_t)m->nl;
    size_t e_nl = (size_t)E * (size_t)m->nl;
    m->hnode     = calloc(n_nl, sizeof(real_t));
    m->hnode_new = calloc(n_nl, sizeof(real_t));
    m->helem     = calloc(e_nl, sizeof(real_t));
    m->hbar      = calloc((size_t)N, sizeof(real_t));
    m->hbar_old  = calloc((size_t)N, sizeof(real_t));
    FESOM_CHECK(m->hnode && m->hnode_new && m->helem && m->hbar && m->hbar_old,
                "mesh state alloc: out of memory");
}

/*--- I/O helpers ------------------------------------------------------------*/

static FILE *open_mesh_file(const char *mesh_dir, const char *name)
{
    char path[1024];
    int n = snprintf(path, sizeof(path), "%s/%s", mesh_dir, name);
    FESOM_CHECK(n > 0 && (size_t)n < sizeof(path),
                "mesh path too long: %s/%s", mesh_dir, name);
    FILE *fh = fopen(path, "r");
    FESOM_CHECK(fh != NULL, "cannot open %s", path);
    return fh;
}

/*--- Rotation ---------------------------------------------------------------
 * Mirror of g_rotate_grid::set_mesh_transform_matrix and g2r in Fortran.
 * Convention: rotate around z by alpha, around new x by beta, around new z
 * by gamma. Default Euler angles set in fesom_constants.h come from
 * work_pi and work_core namelist.config.
 */

static void build_rotation_matrix(real_t M[9])
{
    const real_t al = FESOM_ALPHA_EULER_DEG * FESOM_RAD;
    const real_t be = FESOM_BETA_EULER_DEG  * FESOM_RAD;
    const real_t ga = FESOM_GAMMA_EULER_DEG * FESOM_RAD;
    /* row-major flat 3x3 — same numeric layout as Fortran r2g_matrix(i,j) */
    M[0] = cos(ga)*cos(al) - sin(ga)*cos(be)*sin(al);
    M[1] = cos(ga)*sin(al) + sin(ga)*cos(be)*cos(al);
    M[2] = sin(ga)*sin(be);
    M[3] = -sin(ga)*cos(al) - cos(ga)*cos(be)*sin(al);
    M[4] = -sin(ga)*sin(al) + cos(ga)*cos(be)*cos(al);
    M[5] = cos(ga)*sin(be);
    M[6] = sin(be)*sin(al);
    M[7] = -sin(be)*cos(al);
    M[8] = cos(be);
}

/* Geographic (rad) -> rotated (rad). Mirrors Fortran g2r. */
static void g2r(const real_t M[9], real_t glon, real_t glat,
                real_t *rlon, real_t *rlat)
{
    real_t xg = cos(glat) * cos(glon);
    real_t yg = cos(glat) * sin(glon);
    real_t zg = sin(glat);

    real_t xr = M[0]*xg + M[1]*yg + M[2]*zg;
    real_t yr = M[3]*xg + M[4]*yg + M[5]*zg;
    real_t zr = M[6]*xg + M[7]*yg + M[8]*zg;

    *rlat = asin(zr);
    if (yr == 0.0 && xr == 0.0) {
        *rlon = 0.0;
    } else {
        *rlon = atan2(yr, xr);
    }
}

/* Rotated (rad) -> geographic (rad). Inverse of g2r; for an orthogonal
 * rotation matrix this is M^T applied to (xr, yr, zr). */
static void r2g(const real_t M[9], real_t rlon, real_t rlat,
                real_t *glon, real_t *glat)
{
    real_t xr = cos(rlat) * cos(rlon);
    real_t yr = cos(rlat) * sin(rlon);
    real_t zr = sin(rlat);

    real_t xg = M[0]*xr + M[3]*yr + M[6]*zr;
    real_t yg = M[1]*xr + M[4]*yr + M[7]*zr;
    real_t zg = M[2]*xr + M[5]*yr + M[8]*zr;

    *glat = asin(zg);
    if (yg == 0.0 && xg == 0.0) {
        *glon = 0.0;
    } else {
        *glon = atan2(yg, xg);
    }
}

/* Rotate a 2D vector with geographic (east, north) components into the model's
 * rotated frame, at a node whose geographic (glon,glat) and rotated (rlon,rlat)
 * coordinates [radians] are known. Mirror of Fortran vector_g2r
 * (gen_modules_rotate_grid.F90:120). M maps geo Cartesian -> rot Cartesian
 * (== Fortran r2g_matrix). No-op when rotation is disabled. Magnitude-preserving,
 * so wind-speed-based thermodynamics (|wind|) is unaffected; only direction. */
void fesom_vector_g2r(real_t *u, real_t *v,
                      real_t glon, real_t glat, real_t rlon, real_t rlat)
{
#if FESOM_FORCE_ROTATION
    static real_t M[9];
    static int built = 0;
    if (!built) { build_rotation_matrix(M); built = 1; }

    real_t tlon = *u, tlat = *v;              /* geographic east/north comps */
    /* geographic vector -> Cartesian geographic */
    real_t txg = -tlat*sin(glat)*cos(glon) - tlon*sin(glon);
    real_t tyg = -tlat*sin(glat)*sin(glon) + tlon*cos(glon);
    real_t tzg =  tlat*cos(glat);
    /* Cartesian geographic -> Cartesian rotated */
    real_t txr = M[0]*txg + M[1]*tyg + M[2]*tzg;
    real_t tyr = M[3]*txg + M[4]*tyg + M[5]*tzg;
    real_t tzr = M[6]*txg + M[7]*tyg + M[8]*tzg;
    /* Cartesian rotated -> rotated east/north components */
    *v = -sin(rlat)*cos(rlon)*txr - sin(rlat)*sin(rlon)*tyr + cos(rlat)*tzr;
    *u = -sin(rlon)*txr + cos(rlon)*tyr;
#else
    (void)u; (void)v; (void)glon; (void)glat; (void)rlon; (void)rlat;
#endif
}

/*--- nod2d.out --------------------------------------------------------------
 * Format:
 *   <count>
 *   <id> <lon_deg> <lat_deg> <coast_flag>      ... count lines
 *
 * Stored: geo_coord_nod2D in geographic radians; coord_nod2D in rotated
 * radians (same convention as Fortran mesh%coord_nod2D).
 */
static void read_nod2d(fesom_mesh *m, const char *mesh_dir)
{
    FILE *fh = open_mesh_file(mesh_dir, "nod2d.out");
    int n = 0;
    FESOM_CHECK(fscanf(fh, "%d", &n) == 1, "nod2d.out: missing header count");
    FESOM_CHECK(n > 0, "nod2d.out: invalid count %d", n);

    m->nod2D = n;
    m->myDim_nod2D = n;
    m->eDim_nod2D  = 0;
    m->coord_nod2D     = malloc((size_t)n * 2 * sizeof(real_t));
    m->geo_coord_nod2D = malloc((size_t)n * 2 * sizeof(real_t));
    m->coast_flag      = malloc((size_t)n * sizeof(int));
    FESOM_CHECK(m->coord_nod2D && m->geo_coord_nod2D && m->coast_flag,
                "nod2d.out: out of memory");

    real_t M[9];
    build_rotation_matrix(M);

    for (int i = 0; i < n; ++i) {
        int id, cf;
        double lon_deg, lat_deg;
        int got = fscanf(fh, "%d %lf %lf %d", &id, &lon_deg, &lat_deg, &cf);
        FESOM_CHECK(got == 4, "nod2d.out: parse error at row %d (got %d)", i + 1, got);
        FESOM_CHECK(id == i + 1, "nod2d.out: expected id=%d at row %d, got %d",
                    i + 1, i + 1, id);

        real_t glon = (real_t)lon_deg * FESOM_RAD;
        real_t glat = (real_t)lat_deg * FESOM_RAD;
        m->geo_coord_nod2D[2*i + 0] = glon;
        m->geo_coord_nod2D[2*i + 1] = glat;

        if (FESOM_FORCE_ROTATION) {
            real_t rlon, rlat;
            g2r(M, glon, glat, &rlon, &rlat);
            m->coord_nod2D[2*i + 0] = rlon;
            m->coord_nod2D[2*i + 1] = rlat;
        } else {
            m->coord_nod2D[2*i + 0] = glon;
            m->coord_nod2D[2*i + 1] = glat;
        }
        m->coast_flag[i] = cf;
    }
    fclose(fh);
}

/*--- elem2d.out -------------------------------------------------------------
 * Format:
 *   <count>
 *   <n1> <n2> <n3>          ... count lines, 1-based node IDs
 */
static void read_elem2d(fesom_mesh *m, const char *mesh_dir)
{
    FILE *fh = open_mesh_file(mesh_dir, "elem2d.out");
    int n = 0;
    FESOM_CHECK(fscanf(fh, "%d", &n) == 1, "elem2d.out: missing header");
    FESOM_CHECK(n > 0, "elem2d.out: invalid count %d", n);

    m->elem2D       = n;
    m->myDim_elem2D = n;
    m->eDim_elem2D  = 0;
    m->elem_nodes   = malloc((size_t)n * 3 * sizeof(int));
    FESOM_CHECK(m->elem_nodes, "elem2d.out: out of memory");

    for (int i = 0; i < n; ++i) {
        int a, b, c;
        int got = fscanf(fh, "%d %d %d", &a, &b, &c);
        FESOM_CHECK(got == 3, "elem2d.out: parse error at row %d (got %d)", i + 1, got);
        a -= 1; b -= 1; c -= 1;
        FESOM_CHECK(a >= 0 && a < m->nod2D
                 && b >= 0 && b < m->nod2D
                 && c >= 0 && c < m->nod2D,
                 "elem2d.out row %d: node id out of range (%d %d %d), nod2D=%d",
                 i + 1, a, b, c, m->nod2D);
        m->elem_nodes[3*i + 0] = a;
        m->elem_nodes[3*i + 1] = b;
        m->elem_nodes[3*i + 2] = c;
    }
    fclose(fh);
}

/*--- aux3d.out --------------------------------------------------------------
 * Format (FESOM2.0, depth-at-node):
 *   <nl>
 *   <zbar(1)> <zbar(2)> ... <zbar(nl)>            (free-format, may span lines)
 *   <depth(1)> <depth(2)> ... <depth(nod2D)>      (one per node)
 *
 * Mirror of oce_mesh.F90:567-707. Fortran-side: zbar is forced negative
 * (line 578 if zbar(2)>0, negate); depth same (line 700: if x>0, x=-x).
 */
static void read_aux3d(fesom_mesh *m, const char *mesh_dir)
{
    FILE *fh = open_mesh_file(mesh_dir, "aux3d.out");
    int nl = 0;
    FESOM_CHECK(fscanf(fh, "%d", &nl) == 1, "aux3d.out: missing nl");
    FESOM_CHECK(nl >= 3, "aux3d.out: nl=%d too small", nl);
    m->nl   = nl;
    m->zbar = malloc((size_t)nl * sizeof(real_t));
    FESOM_CHECK(m->zbar, "aux3d.out: out of memory (zbar)");

    for (int k = 0; k < nl; ++k) {
        double v;
        FESOM_CHECK(fscanf(fh, "%lf", &v) == 1, "aux3d.out: zbar truncated at k=%d", k);
        m->zbar[k] = (real_t)v;
    }
    /* Force zbar to be negative downward (matches Fortran line 578) */
    if (m->zbar[1] > 0.0) {
        for (int k = 0; k < nl; ++k) m->zbar[k] = -m->zbar[k];
    }

    /* Mid-layer depths Z[nz] = 0.5*(zbar[nz]+zbar[nz+1]) for nz=0..nl-2.
       Mirror of mesh%Z = 0.5*(zbar(1:nl-1)+zbar(2:nl)) at oce_mesh.F90:580. */
    m->Z = malloc((size_t)(nl - 1) * sizeof(real_t));
    FESOM_CHECK(m->Z, "aux3d.out: out of memory (Z)");
    for (int k = 0; k < nl - 1; ++k) {
        m->Z[k] = 0.5 * (m->zbar[k] + m->zbar[k + 1]);
    }

    m->depth = malloc((size_t)m->nod2D * sizeof(real_t));
    FESOM_CHECK(m->depth, "aux3d.out: out of memory (depth)");
    for (int i = 0; i < m->nod2D; ++i) {
        double v;
        FESOM_CHECK(fscanf(fh, "%lf", &v) == 1,
                    "aux3d.out: depth truncated at node %d", i);
        if (v > 0.0) v = -v;             /* Fortran line 700 */
        m->depth[i] = (real_t)v;
    }
    fclose(fh);
}

/*--- nlvls.out / elvls.out --------------------------------------------------*/
static void read_int_per_line(const char *path, int *out, int n)
{
    FILE *fh = fopen(path, "r");
    FESOM_CHECK(fh != NULL, "cannot open %s", path);
    for (int i = 0; i < n; ++i) {
        FESOM_CHECK(fscanf(fh, "%d", &out[i]) == 1,
                    "%s: parse error at row %d", path, i + 1);
    }
    fclose(fh);
}

static void read_nlvls(fesom_mesh *m, const char *mesh_dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/nlvls.out", mesh_dir);
    m->nlevels_nod2D = malloc((size_t)m->nod2D * sizeof(int));
    FESOM_CHECK(m->nlevels_nod2D, "nlvls.out: out of memory");
    read_int_per_line(path, m->nlevels_nod2D, m->nod2D);
}

static void read_elvls(fesom_mesh *m, const char *mesh_dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/elvls.out", mesh_dir);
    m->nlevels = malloc((size_t)m->elem2D * sizeof(int));
    FESOM_CHECK(m->nlevels, "elvls.out: out of memory");
    read_int_per_line(path, m->nlevels, m->elem2D);
}

/*--- edges.out / edge_tri.out -----------------------------------------------
 * One edge per line; both files must agree on edge count.
 */
static int count_lines(const char *path)
{
    FILE *fh = fopen(path, "r");
    FESOM_CHECK(fh != NULL, "cannot open %s", path);
    int lines = 0, c, prev = '\n';
    while ((c = fgetc(fh)) != EOF) {
        if (c == '\n') ++lines;
        prev = c;
    }
    if (prev != '\n' && prev != EOF) ++lines;
    fclose(fh);
    return lines;
}

static void read_edges(fesom_mesh *m, const char *mesh_dir)
{
    char path_e[1024], path_t[1024];
    snprintf(path_e, sizeof(path_e), "%s/edges.out",    mesh_dir);
    snprintf(path_t, sizeof(path_t), "%s/edge_tri.out", mesh_dir);

    int ne = count_lines(path_e);
    int nt = count_lines(path_t);
    FESOM_CHECK(ne == nt, "edges.out (%d) and edge_tri.out (%d) row mismatch", ne, nt);
    m->edge2D   = ne;
    m->edges    = malloc((size_t)ne * 2 * sizeof(int));
    m->edge_tri = malloc((size_t)ne * 2 * sizeof(int));
    FESOM_CHECK(m->edges && m->edge_tri, "edges: out of memory");

    FILE *fe = fopen(path_e, "r");
    FILE *ft = fopen(path_t, "r");
    FESOM_CHECK(fe && ft, "edges/edge_tri: cannot reopen");
    for (int i = 0; i < ne; ++i) {
        int a, b;
        FESOM_CHECK(fscanf(fe, "%d %d", &a, &b) == 2,
                    "edges.out: parse error at row %d", i + 1);
        m->edges[2*i + 0] = a - 1;
        m->edges[2*i + 1] = b - 1;

        int e1, e2;
        FESOM_CHECK(fscanf(ft, "%d %d", &e1, &e2) == 2,
                    "edge_tri.out: parse error at row %d", i + 1);
        /* Fortran stores -1 (or sometimes 0) for boundary edges; keep -1 in C */
        m->edge_tri[2*i + 0] = (e1 > 0) ? e1 - 1 : -1;
        m->edge_tri[2*i + 1] = (e2 > 0) ? e2 - 1 : -1;
    }
    fclose(fe);
    fclose(ft);

    /* Count interior edges = those with both adjacent elements valid in the
     * GLOBAL mesh. Mirrors Fortran mesh%edge2D_in (mesh-prep convention is
     * to list interior edges first, boundary edges last). */
    int interior = 0;
    for (int i = 0; i < ne; ++i) {
        if (m->edge_tri[2*i + 0] >= 0 && m->edge_tri[2*i + 1] >= 0) ++interior;
    }
    m->edge2D_in = interior;
}

/*--- Element orientation (test_tri port) ------------------------------------
 * Mirrors oce_mesh.F90:1679-1729. For each triangle, compute the cross
 * product of edges (n1->n2) x (n1->n3) in the (lon, lat) plane after a
 * cyclic-length wrap. If the result is positive (CCW), swap nodes 2 and 3
 * to enforce CW orientation. This is critical for correct stiffness-matrix
 * sign in the SSH solver — FRESH_START.md §4 + §11.
 */
static int orient_cw(fesom_mesh *m)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    int swapped = 0;
    /* elem_nodes is sized myDim_elem2D only (Fortran convention) */
    int Ei = m->myDim_elem2D;
    for (int e = 0; e < Ei; ++e) {
        int n0 = m->elem_nodes[3*e + 0];
        int n1 = m->elem_nodes[3*e + 1];
        int n2 = m->elem_nodes[3*e + 2];
        real_t ax = m->coord_nod2D[2*n0 + 0];
        real_t ay = m->coord_nod2D[2*n0 + 1];
        real_t bx = m->coord_nod2D[2*n1 + 0] - ax;
        real_t by = m->coord_nod2D[2*n1 + 1] - ay;
        real_t cx = m->coord_nod2D[2*n2 + 0] - ax;
        real_t cy = m->coord_nod2D[2*n2 + 1] - ay;
        if (bx >  half_cyc) bx -= cyc;
        if (bx < -half_cyc) bx += cyc;
        if (cx >  half_cyc) cx -= cyc;
        if (cx < -half_cyc) cx += cyc;
        real_t r = bx * cy - by * cx;
        if (r > 0.0) {
            m->elem_nodes[3*e + 1] = n2;
            m->elem_nodes[3*e + 2] = n1;
            ++swapped;
        }
    }
    return swapped;
}

/*--- Build node→element CSR adjacency ---------------------------------------
 * Mirror of how Fortran's nod_in_elem2D is constructed (oce_mesh.F90 fills
 * mesh%nod_in_elem2D_num then mesh%nod_in_elem2D). Two-pass: count, then
 * scatter using a running offset.
 */
static void build_nod_in_elem2D(fesom_mesh *m)
{
    /* CSR rows = local nodes (myDim+eDim). Source elements are ALL local
     * elements (myDim+eDim+eXDim). Mirror of Fortran oce_mesh.F90:2010-2074:
     * the Fortran does a multi-pass scheme with exchange_nod on global elem
     * IDs to fill nod_in_elem2D for halo nodes too. We do it directly: since
     * scatter_mesh now fills elem_nodes for myDim+eDim+eXDim with local node
     * IDs (or -1 for outer-halo unmappable), we can just iterate every local
     * element and build the adjacency.
     *
     * This is critical for the compute_areas pipeline: a partition-boundary
     * node's area is the sum of (1/3 * elem_area) over ALL surrounding
     * elements — including halo elements on this rank. Without those, the
     * boundary area is underestimated → SSH stiffness mass term too small →
     * CG converges to oversized eta → blow-up. (oce_mesh.F90:2266) */
    const int N = m->myDim_nod2D + m->eDim_nod2D;
    const int E = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    m->nod_in_elem2D_offsets = calloc((size_t)N + 1, sizeof(int));
    FESOM_CHECK(m->nod_in_elem2D_offsets, "nod_in_elem2D_offsets: out of memory");

    /* Pass 1: count per node. Skip elem entries with v=-1 (outer-halo unmappable). */
    for (int e = 0; e < E; ++e) {
        for (int k = 0; k < 3; ++k) {
            int v = m->elem_nodes[3*e + k];
            if (v < 0 || v >= N) continue;
            ++m->nod_in_elem2D_offsets[v + 1];
        }
    }
    for (int i = 1; i <= N; ++i) {
        m->nod_in_elem2D_offsets[i] += m->nod_in_elem2D_offsets[i - 1];
    }
    int total = m->nod_in_elem2D_offsets[N];

    /* Pass 2: scatter. */
    int *cursor = calloc((size_t)N, sizeof(int));
    FESOM_CHECK(cursor, "nod_in_elem2D: cursor OOM");
    m->nod_in_elem2D = malloc((size_t)total * sizeof(int));
    FESOM_CHECK(m->nod_in_elem2D, "nod_in_elem2D: out of memory");

    for (int e = 0; e < E; ++e) {
        for (int k = 0; k < 3; ++k) {
            int v = m->elem_nodes[3*e + k];
            if (v < 0 || v >= N) continue;
            int slot = m->nod_in_elem2D_offsets[v] + cursor[v]++;
            m->nod_in_elem2D[slot] = e;
        }
    }
    free(cursor);
}

/*--- K_v⁻ and ulevels (Phase 1: no cavity, ulevels = 1 everywhere) ----------
 * K_v⁻ = MIN_{c ∈ C(v)} K_c — used as the ALE deformation limit
 * (FRESH_START.md §2 GOTCHA, geometry.rst Partial cells section).
 */
static void compute_vertical_levels_aux(fesom_mesh *m)
{
    int N = m->myDim_nod2D + m->eDim_nod2D;
    int E = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    m->ulevels           = malloc((size_t)E * sizeof(int));
    m->ulevels_nod2D     = malloc((size_t)N * sizeof(int));
    m->nlevels_nod2D_min = malloc((size_t)N * sizeof(int));
    m->ulevels_nod2D_max = malloc((size_t)N * sizeof(int));
    FESOM_CHECK(m->ulevels && m->ulevels_nod2D && m->nlevels_nod2D_min
             && m->ulevels_nod2D_max,
                "vertical aux: out of memory");

    for (int e = 0; e < E; ++e) m->ulevels[e]       = 1;
    for (int n = 0; n < N; ++n) m->ulevels_nod2D[n] = 1;

    for (int n = 0; n < N; ++n) {
        int o0 = m->nod_in_elem2D_offsets[n];
        int o1 = m->nod_in_elem2D_offsets[n + 1];
        if (o1 == o0) {
            /* Halo node may have no local-element neighbours if it sits at the
             * extreme outer boundary of the partition; mirror nlevels_nod2D
             * so K_v⁻ is at least sensible. */
            m->nlevels_nod2D_min[n] = m->nlevels_nod2D[n];
            m->ulevels_nod2D_max[n] = m->ulevels_nod2D[n];
            continue;
        }
        int mn = m->nlevels[m->nod_in_elem2D[o0]];
        int mx_u = m->ulevels[m->nod_in_elem2D[o0]];
        for (int k = o0 + 1; k < o1; ++k) {
            int e  = m->nod_in_elem2D[k];
            int v  = m->nlevels[e];
            int vu = m->ulevels[e];
            if (v  < mn  ) mn   = v;
            if (vu > mx_u) mx_u = vu;
        }
        m->nlevels_nod2D_min[n] = mn;
        m->ulevels_nod2D_max[n] = mx_u;
    }
}

/*--- GM/Redi mesh prerequisites -------------------------------------------- *
 *  mesh_resolution[n] — Voronoi diameter, smoothed 3 passes (Fortran
 *      oce_mesh.F90:2353-2378).  All loops match Fortran bounds exactly.
 *  zbar_3d_n[n*nl + nz] — per-node interface depths.  In our
 *      linfs / no-cavity / no-partial-cells config this collapses to
 *      zbar[nz] for nz < nlevels_nod2D[n], 0 elsewhere.  Built once at
 *      setup; constant in time for linfs.
 * ------------------------------------------------------------------------ */
static void compute_mesh_resolution(fesom_mesh *m, fesom_partit *partit)
{
    int N      = m->myDim_nod2D + m->eDim_nod2D;
    int N_own  = m->myDim_nod2D;
    int nl     = m->nl;

    m->mesh_resolution = malloc((size_t)N * sizeof(real_t));
    FESOM_CHECK(m->mesh_resolution, "mesh_resolution: out of memory");

    /* Voronoi diameter: 2 * sqrt(areasvol(ulevels_nod2D(n), n) / pi).
     * Fortran loop bound is myDim+eDim (oce_mesh.F90:2356). */
    const real_t inv_pi = 1.0 / 3.14159265358979323846;
    for (int n = 0; n < N; ++n) {
        int ul = m->ulevels_nod2D[n] - 1;          /* 0-based */
        real_t a = m->areasvol[(size_t)n * nl + ul];
        m->mesh_resolution[n] = 2.0 * sqrt(a * inv_pi);
    }

    /* 3 passes of mass-matrix smoothing: per-node weighted mean of
     * surrounding-element vertex resolutions (oce_mesh.F90:2361-2377).
     * Fortran writes only myDim then exchange_nod, three times. */
    real_t *work = calloc((size_t)N_own, sizeof(real_t));
    FESOM_CHECK(work, "mesh_resolution smoothing: out of memory");
    for (int pass = 0; pass < 3; ++pass) {
        for (int n = 0; n < N_own; ++n) {
            int o0 = m->nod_in_elem2D_offsets[n];
            int o1 = m->nod_in_elem2D_offsets[n + 1];
            real_t vol = 0.0, acc = 0.0;
            for (int k = o0; k < o1; ++k) {
                int elem = m->nod_in_elem2D[k];
                int v0 = m->elem_nodes[3*elem + 0];
                int v1 = m->elem_nodes[3*elem + 1];
                int v2 = m->elem_nodes[3*elem + 2];
                real_t mean3 = (m->mesh_resolution[v0]
                              + m->mesh_resolution[v1]
                              + m->mesh_resolution[v2]) / 3.0;
                acc += mean3 * m->elem_area[elem];
                vol += m->elem_area[elem];
            }
            work[n] = (vol > 0.0) ? acc / vol : 0.0;
        }
        for (int n = 0; n < N_own; ++n) m->mesh_resolution[n] = work[n];
        if (partit && partit->npes > 1) {
            fesom_exchange_nod2D(m->mesh_resolution, partit);
        }
    }
    free(work);
}

static void compute_zbar_3d_n(fesom_mesh *m)
{
    int N  = m->myDim_nod2D + m->eDim_nod2D;
    int nl = m->nl;
    m->zbar_3d_n = calloc((size_t)N * nl, sizeof(real_t));
    FESOM_CHECK(m->zbar_3d_n, "zbar_3d_n: out of memory");

    /* Fortran oce_ale.F90:342-375 with our config (no cavity, no partial
     * cells): zbar_3d_n(nz, n) = zbar(nz) for nz in [ulevels..nlevels].
     * ulevels_nod2D = 1 always for us, so for nz from 0 to nlevels-1
     * we copy zbar[nz]. Beyond nlevels we leave 0. Same for halo nodes. */
    for (int n = 0; n < N; ++n) {
        int nzmax = m->nlevels_nod2D[n];   /* 1-based count */
        for (int nz = 0; nz < nzmax; ++nz) {
            m->zbar_3d_n[(size_t)n * nl + nz] = m->zbar[nz];
        }
    }
}

/*--- elem_area, area, areasvol ---------------------------------------------
 * Mirror of oce_mesh.F90:2202-2315. Sequence:
 *   ay = cos(mean lat over the 3 vertices)             [cellwise constant]
 *   a, b = vertex2-vertex1, vertex3-vertex1             (rotated radians)
 *   trim_cyclic on the lon (x) components
 *   scale lon components by ay (local-Cartesian projection)
 *   elem_area = 0.5 * |a.x * b.y - b.x * a.y|           (radian² for now)
 *   area(nz, n) += elem_area / 3 for each layer where the cell exists
 *   non-cavity: areasvol(nz, n) = area(nz, n)
 *   * R_earth² to convert all areas to m²
 */
/* Stage 1: compute elem_area for myDim only and allocate node arrays.
 * Stage 2 must be called after elem_area has been halo-exchanged. */
static void compute_elem_area_only(fesom_mesh *m)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    const real_t r2       = (real_t)FESOM_R_EARTH * (real_t)FESOM_R_EARTH;

    int N  = m->myDim_nod2D + m->eDim_nod2D;
    int E  = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    int Ei = m->myDim_elem2D;
    m->elem_area = calloc((size_t)E, sizeof(real_t));
    m->area      = calloc((size_t)N * (size_t)m->nl, sizeof(real_t));
    m->areasvol  = calloc((size_t)N * (size_t)m->nl, sizeof(real_t));
    FESOM_CHECK(m->elem_area && m->area && m->areasvol,
                "areas: out of memory");

    for (int e = 0; e < Ei; ++e) {
        int n0 = m->elem_nodes[3*e + 0];
        int n1 = m->elem_nodes[3*e + 1];
        int n2 = m->elem_nodes[3*e + 2];
        real_t ay = (m->coord_nod2D[2*n0 + 1] +
                     m->coord_nod2D[2*n1 + 1] +
                     m->coord_nod2D[2*n2 + 1]) / 3.0;
        ay = cos(ay);
        real_t ax = m->coord_nod2D[2*n1 + 0] - m->coord_nod2D[2*n0 + 0];
        real_t aly = m->coord_nod2D[2*n1 + 1] - m->coord_nod2D[2*n0 + 1];
        real_t bx = m->coord_nod2D[2*n2 + 0] - m->coord_nod2D[2*n0 + 0];
        real_t bly = m->coord_nod2D[2*n2 + 1] - m->coord_nod2D[2*n0 + 1];
        if (ax >  half_cyc) ax -= cyc;
        if (ax < -half_cyc) ax += cyc;
        if (bx >  half_cyc) bx -= cyc;
        if (bx < -half_cyc) bx += cyc;
        ax *= ay;
        bx *= ay;
        real_t cross = ax * bly - bx * aly;
        if (cross < 0) cross = -cross;
        m->elem_area[e] = 0.5 * cross;
    }
    /* Convert elem_area to m² before exchange (so halo receives m² values). */
    for (int e = 0; e < Ei; ++e) m->elem_area[e] *= r2;
}

/* Stage 2: build per-node area / areasvol from elem_area (halo-exchanged).
 * Mirror of Fortran oce_mesh.F90:2266-2293. */
static void compute_node_areas(fesom_mesh *m)
{
    int N  = m->myDim_nod2D + m->eDim_nod2D;

    /* distribute cell area to surrounding nodes per layer (no cavity).
     * Iterate over all local nodes (interior + halo). elem_area now valid
     * on halo elements after exchange, so boundary nodes get correct area. */
    for (int n = 0; n < N; ++n) {
        int o0 = m->nod_in_elem2D_offsets[n];
        int o1 = m->nod_in_elem2D_offsets[n + 1];
        for (int k = o0; k < o1; ++k) {
            int e = m->nod_in_elem2D[k];
            int nzmin = m->ulevels[e] - 1;
            int nzmax = m->nlevels[e] - 1;
            real_t third = m->elem_area[e] / 3.0;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                m->area[FESOM_NODE3D(n, nz, m->nl)] += third;
            }
        }
        int nzmin_n = m->ulevels_nod2D[n] - 1;
        int nzmax_n = m->nlevels_nod2D[n] - 1;
        for (int nz = nzmin_n; nz < nzmax_n; ++nz) {
            m->areasvol[FESOM_NODE3D(n, nz, m->nl)] = m->area[FESOM_NODE3D(n, nz, m->nl)];
        }
    }

    /* ocean_area — local sum over INTERIOR nodes (Allreduce in caller). */
    m->ocean_area = 0.0;
    for (int n = 0; n < m->myDim_nod2D; ++n) {
        if (m->ulevels_nod2D[n] > 1) continue;       /* cavity */
        m->ocean_area += m->areasvol[FESOM_NODE3D(n, 0, m->nl)];
    }
}

/*--- elem/edge centers, cyclic-aware ----------------------------------------
 * Mirror of edge_center / elem_center in oce_mesh.F90:2115-2155.
 */
static void elem_center_xy(const fesom_mesh *m, int e, real_t *x, real_t *y)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    int n0 = m->elem_nodes[3*e + 0];
    int n1 = m->elem_nodes[3*e + 1];
    int n2 = m->elem_nodes[3*e + 2];
    real_t ax[3];
    ax[0] = m->coord_nod2D[2*n0 + 0];
    ax[1] = m->coord_nod2D[2*n1 + 0];
    ax[2] = m->coord_nod2D[2*n2 + 0];
    real_t amin = ax[0];
    if (ax[1] < amin) amin = ax[1];
    if (ax[2] < amin) amin = ax[2];
    for (int k = 0; k < 3; ++k) {
        if (ax[k] - amin >=  half_cyc) ax[k] -= cyc;
        if (ax[k] - amin <  -half_cyc) ax[k] += cyc;
    }
    *x = (ax[0] + ax[1] + ax[2]) / 3.0;
    *y = (m->coord_nod2D[2*n0 + 1] +
          m->coord_nod2D[2*n1 + 1] +
          m->coord_nod2D[2*n2 + 1]) / 3.0;
}

static void edge_center_xy(const fesom_mesh *m, int n1, int n2, real_t *x, real_t *y)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    real_t ax = m->coord_nod2D[2*n1 + 0];
    real_t bx = m->coord_nod2D[2*n2 + 0];
    if (ax - bx >  half_cyc) ax -= cyc;
    if (ax - bx < -half_cyc) bx -= cyc;
    *x = 0.5 * (ax + bx);
    *y = 0.5 * (m->coord_nod2D[2*n1 + 1] + m->coord_nod2D[2*n2 + 1]);
}

/*--- elem_cos, metric_factor, coriolis(_node) -------------------------------
 * Mirror of oce_mesh.F90:2475-2528.
 *   elem_cos[e]      = cos(rotated lat at centroid)
 *   metric_factor[e] = tan(rotated lat at centroid) / R_earth
 *   coriolis[e]      = 2Ω sin(geographic lat at centroid)
 *   coriolis_node[n] = 2Ω sin(geographic lat at node)
 */
static void compute_metric_and_coriolis(fesom_mesh *m)
{
    real_t M[9];
    build_rotation_matrix(M);

    int N  = m->myDim_nod2D + m->eDim_nod2D;
    int E  = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    int Ei = m->myDim_elem2D;
    m->elem_cos      = calloc((size_t)E, sizeof(real_t));
    m->metric_factor = calloc((size_t)E, sizeof(real_t));
    m->coriolis      = calloc((size_t)E, sizeof(real_t));
    m->coriolis_node = malloc((size_t)N * sizeof(real_t));
    m->elem_center_x = calloc((size_t)E, sizeof(real_t));
    m->elem_center_y = calloc((size_t)E, sizeof(real_t));
    FESOM_CHECK(m->elem_cos && m->metric_factor && m->coriolis && m->coriolis_node
             && m->elem_center_x && m->elem_center_y,
                "metric/coriolis: out of memory");

    for (int n = 0; n < N; ++n) {
        m->coriolis_node[n] =
            2.0 * (real_t)FESOM_OMEGA * sin(m->geo_coord_nod2D[2*n + 1]);
    }

    for (int e = 0; e < Ei; ++e) {
        real_t cx, cy;
        elem_center_xy(m, e, &cx, &cy);
        m->elem_center_x[e] = cx;
        m->elem_center_y[e] = cy;
        m->elem_cos[e]      = cos(cy);
        m->metric_factor[e] = tan(cy) / (real_t)FESOM_R_EARTH;
        real_t glon, glat;
        if (FESOM_FORCE_ROTATION) {
            r2g(M, cx, cy, &glon, &glat);
        } else {
            glon = cx; glat = cy;
        }
        m->coriolis[e] = 2.0 * (real_t)FESOM_OMEGA * sin(glat);
    }
}

/*--- edge_dxdy / edge_cross_dxdy --------------------------------------------
 * Mirror of oce_mesh.F90:2534-2573. Note unit asymmetry per geometry.rst:
 *   edge_dxdy is in RADIANS (not multiplied by R_earth)
 *   edge_cross_dxdy is in METERS (lon scaled by elem_cos, then * R_earth)
 */
static void compute_edges_geometry(fesom_mesh *m)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    const real_t r_earth  = (real_t)FESOM_R_EARTH;

    int E  = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    int EG = m->myDim_edge2D + m->eDim_edge2D;
    m->edge_dxdy       = malloc((size_t)EG * 2 * sizeof(real_t));
    m->edge_cross_dxdy = malloc((size_t)EG * 4 * sizeof(real_t));
    FESOM_CHECK(m->edge_dxdy && m->edge_cross_dxdy, "edges geometry: out of memory");

    /* Cell centroids — use the persistent mesh fields. For interior elements
     * these are filled by compute_metric_and_coriolis above; for halo elements
     * (eDim/eXDim) they came in via the elem_center_x/y halo exchange in
     * fesom_mesh_compute_metrics. We use them directly (no local cache).
     *
     * Note: elem_nodes is only sized myDim_elem2D — calling elem_center_xy
     * for halo element indices would segfault. Always use mesh->elem_center_*
     * for any code path that needs halo-element centers. */
    real_t *cx = m->elem_center_x;
    real_t *cy = m->elem_center_y;

    for (int eg = 0; eg < EG; ++eg) {
        int n1 = m->edges[2*eg + 0];
        int n2 = m->edges[2*eg + 1];

        /* edge_dxdy = node2 − node1 (radians, lon-cyclic-aware) */
        real_t dx = m->coord_nod2D[2*n2 + 0] - m->coord_nod2D[2*n1 + 0];
        real_t dy = m->coord_nod2D[2*n2 + 1] - m->coord_nod2D[2*n1 + 1];
        if (dx >  half_cyc) dx -= cyc;
        if (dx < -half_cyc) dx += cyc;
        m->edge_dxdy[2*eg + 0] = dx;
        m->edge_dxdy[2*eg + 1] = dy;

        /* edge_cross_dxdy: from edge midpoint to each adjacent cell center,
           in meters with lon scaled by that cell's elem_cos. */
        real_t emx, emy;
        edge_center_xy(m, n1, n2, &emx, &emy);
        int el1 = m->edge_tri[2*eg + 0];
        int el2 = m->edge_tri[2*eg + 1];

        if (el1 >= 0) {
            real_t bx = cx[el1] - emx;
            real_t by = cy[el1] - emy;
            if (bx >  half_cyc) bx -= cyc;
            if (bx < -half_cyc) bx += cyc;
            bx *= m->elem_cos[el1];
            m->edge_cross_dxdy[4*eg + 0] = bx * r_earth;
            m->edge_cross_dxdy[4*eg + 1] = by * r_earth;
        } else {
            /* Boundary on the "left" side — Fortran would not normally hit this
               (edge_tri[0] is always defined); guard anyway. */
            m->edge_cross_dxdy[4*eg + 0] = 0.0;
            m->edge_cross_dxdy[4*eg + 1] = 0.0;
        }

        if (el2 >= 0) {
            real_t bx = cx[el2] - emx;
            real_t by = cy[el2] - emy;
            if (bx >  half_cyc) bx -= cyc;
            if (bx < -half_cyc) bx += cyc;
            bx *= m->elem_cos[el2];
            m->edge_cross_dxdy[4*eg + 2] = bx * r_earth;
            m->edge_cross_dxdy[4*eg + 3] = by * r_earth;
        } else {
            m->edge_cross_dxdy[4*eg + 2] = 0.0;
            m->edge_cross_dxdy[4*eg + 3] = 0.0;
        }
    }

    /* No free — cx/cy point into mesh-owned arrays. */
}

/*--- gradient_sca -----------------------------------------------------------
 * Linear-shape-function gradient at each cell. Mirror of oce_mesh.F90:2619-2641.
 * Layout: gradient_sca[6e + i] for i = 0..5 = (dN1/dx, dN2/dx, dN3/dx, dN1/dy, dN2/dy, dN3/dy).
 *
 * Derivation (from Fortran comments lines 2589-2617):
 *   dN2/dx = dy31 / (2 A_tri)   dN3/dx = -dy21 / (2 A_tri)
 *   dN1/dx = -(dN2/dx + dN3/dx) = (-dy31 + dy21) / (2 A_tri)
 *   dN2/dy = dx31 / (2 A_tri)   dN3/dy = -dx21 / (2 A_tri)
 *   dN1/dy = -(dN2/dy + dN3/dy) = ( dx31 - dx21) / (2 A_tri)
 *
 * Fortran factors out: dfactor = -0.5 * R_earth / elem_area, with deltaX
 * components scaled by elem_cos. Sign convention: this dfactor is what the
 * physics expects. We mirror it literally — do NOT "fix" the sign.
 */
static void compute_gradient_sca(fesom_mesh *m)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    const real_t r_earth  = (real_t)FESOM_R_EARTH;

    /* Fortran sizes gradient_sca for myDim_elem2D ONLY (oce_mesh.F90:2461) */
    int Ei = m->myDim_elem2D;
    m->gradient_sca = malloc((size_t)Ei * 6 * sizeof(real_t));
    FESOM_CHECK(m->gradient_sca, "gradient_sca: out of memory");

    for (int e = 0; e < Ei; ++e) {
        int n0 = m->elem_nodes[3*e + 0];
        int n1 = m->elem_nodes[3*e + 1];
        int n2 = m->elem_nodes[3*e + 2];

        real_t dX31 = m->coord_nod2D[2*n2 + 0] - m->coord_nod2D[2*n0 + 0];
        if (dX31 >  half_cyc) dX31 -= cyc;
        if (dX31 < -half_cyc) dX31 += cyc;
        dX31 *= m->elem_cos[e];

        real_t dX21 = m->coord_nod2D[2*n1 + 0] - m->coord_nod2D[2*n0 + 0];
        if (dX21 >  half_cyc) dX21 -= cyc;
        if (dX21 < -half_cyc) dX21 += cyc;
        dX21 *= m->elem_cos[e];

        real_t dY31 = m->coord_nod2D[2*n2 + 1] - m->coord_nod2D[2*n0 + 1];
        real_t dY21 = m->coord_nod2D[2*n1 + 1] - m->coord_nod2D[2*n0 + 1];

        real_t dfactor = -0.5 * r_earth / m->elem_area[e];

        m->gradient_sca[6*e + 0] = (-dY31 + dY21) * dfactor;
        m->gradient_sca[6*e + 1] = ( dY31)        * dfactor;
        m->gradient_sca[6*e + 2] = (-dY21)        * dfactor;
        m->gradient_sca[6*e + 3] = ( dX31 - dX21) * dfactor;
        m->gradient_sca[6*e + 4] = (-dX31)        * dfactor;
        m->gradient_sca[6*e + 5] = ( dX21)        * dfactor;
    }
}

/*--- Bcast + scatter to per-rank slices -------------------------------------
 * After this call, every per-rank array on every rank is sized to local
 * extent (myDim+eDim for nodes, myDim+eDim+eXDim for elements,
 * myDim+eDim for edges). Global node/element/edge IDs in connectivity
 * arrays are translated to LOCAL indices via global→local hash tables.
 * mesh->nod2D / elem2D / edge2D continue to hold GLOBAL counts after
 * scatter (used for snapshot writes that gather to rank 0).
 *
 * For npes==1 with synthesised partit (myList_*[i]=i+1), the extraction is
 * an identity copy and the result is bit-identical to the pre-MPI serial
 * mesh. */
static void bcast_real_array(real_t **buf, int n, fesom_partit *p)
{
    if (p->mype != 0) {
        free(*buf);
        *buf = malloc((size_t)n * sizeof(real_t));
        FESOM_CHECK(*buf, "bcast_real_array: alloc on rank %d", p->mype);
    }
    MPI_CHECK(MPI_Bcast(*buf, n, MPI_DOUBLE, 0, p->MPI_COMM_FESOM));
}

/* Build global-id → local-index hash for nodes. lookup[gid_0_based] = local
 * index in [0, myDim+eDim), or -1 if the global node is not in this rank's
 * local set. We allocate a flat lookup table sized to global nod2D — at
 * CORE2 (126858) this is ~500 KB per rank, trivial. */
static void build_node_g2l(int *lookup, fesom_partit *p, int global_nod2D)
{
    for (int i = 0; i < global_nod2D; ++i) lookup[i] = -1;
    int n = p->myDim_nod2D + p->eDim_nod2D;
    for (int i = 0; i < n; ++i) {
        int gid = p->myList_nod2D[i] - 1;          /* 1-based → 0-based */
        FESOM_CHECK(gid >= 0 && gid < global_nod2D,
                    "build_node_g2l: rank %d local %d → gid %d out of [0,%d)",
                    p->mype, i, gid, global_nod2D);
        lookup[gid] = i;
    }
}

static void scatter_mesh(fesom_mesh *m, fesom_partit *p)
{
    /* Broadcast scalar dims (rank 0 has them after read_*). */
    int dims[5];
    if (p->mype == 0) {
        dims[0] = m->nod2D;  dims[1] = m->elem2D;
        dims[2] = m->edge2D; dims[3] = m->nl;
        dims[4] = m->edge2D_in;
    }
    MPI_CHECK(MPI_Bcast(dims, 5, MPI_INT, 0, p->MPI_COMM_FESOM));
    int g_nod2D  = dims[0];
    int g_elem2D = dims[1];
    int g_edge2D = dims[2];
    int g_nl     = dims[3];
    if (p->mype != 0) {
        m->nod2D = g_nod2D; m->elem2D = g_elem2D;
        m->edge2D = g_edge2D; m->nl = g_nl;
        m->edge2D_in = dims[4];
    }

    /* zbar / Z are global (depth column shared across ranks) — small. */
    bcast_real_array(&m->zbar, g_nl, p);
    if (p->mype != 0) { free(m->Z); m->Z = malloc((size_t)(g_nl - 1) * sizeof(real_t));
                        FESOM_CHECK(m->Z, "Z alloc"); }
    MPI_CHECK(MPI_Bcast(m->Z, g_nl - 1, MPI_DOUBLE, 0, p->MPI_COMM_FESOM));

    /* Broadcast the global per-node and per-elem arrays from rank 0. */
    /* coord_nod2D + geo_coord_nod2D are 2D (lon, lat). */
    if (p->mype != 0) {
        free(m->coord_nod2D);     m->coord_nod2D     = malloc((size_t)g_nod2D * 2 * sizeof(real_t));
        free(m->geo_coord_nod2D); m->geo_coord_nod2D = malloc((size_t)g_nod2D * 2 * sizeof(real_t));
        free(m->coast_flag);      m->coast_flag      = malloc((size_t)g_nod2D * sizeof(int));
        free(m->depth);           m->depth           = malloc((size_t)g_nod2D * sizeof(real_t));
        free(m->nlevels_nod2D);   m->nlevels_nod2D   = malloc((size_t)g_nod2D * sizeof(int));
        free(m->elem_nodes);      m->elem_nodes      = malloc((size_t)g_elem2D * 3 * sizeof(int));
        free(m->nlevels);         m->nlevels         = malloc((size_t)g_elem2D * sizeof(int));
        free(m->edges);           m->edges           = malloc((size_t)g_edge2D * 2 * sizeof(int));
        free(m->edge_tri);        m->edge_tri        = malloc((size_t)g_edge2D * 2 * sizeof(int));
        FESOM_CHECK(m->coord_nod2D && m->geo_coord_nod2D && m->coast_flag
                 && m->depth && m->nlevels_nod2D && m->elem_nodes
                 && m->nlevels && m->edges && m->edge_tri,
                    "scatter_mesh: alloc on rank %d", p->mype);
    }
    MPI_CHECK(MPI_Bcast(m->coord_nod2D,     g_nod2D * 2, MPI_DOUBLE, 0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->geo_coord_nod2D, g_nod2D * 2, MPI_DOUBLE, 0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->coast_flag,      g_nod2D,     MPI_INT,    0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->depth,           g_nod2D,     MPI_DOUBLE, 0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->nlevels_nod2D,   g_nod2D,     MPI_INT,    0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->elem_nodes,      g_elem2D * 3, MPI_INT,   0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->nlevels,         g_elem2D,    MPI_INT,    0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->edges,           g_edge2D * 2, MPI_INT,   0, p->MPI_COMM_FESOM));
    MPI_CHECK(MPI_Bcast(m->edge_tri,        g_edge2D * 2, MPI_INT,   0, p->MPI_COMM_FESOM));

    /* Now extract per-rank slices, replacing the global arrays. */
    int my_n  = p->myDim_nod2D + p->eDim_nod2D;
    int my_e  = p->myDim_elem2D + p->eDim_elem2D + p->eXDim_elem2D;
    int my_eg = p->myDim_edge2D + p->eDim_edge2D;

    /* Build global→local node lookup (used to translate elem_nodes). */
    int *node_g2l = malloc((size_t)g_nod2D * sizeof(int));
    FESOM_CHECK(node_g2l, "scatter_mesh: node_g2l alloc");
    build_node_g2l(node_g2l, p, g_nod2D);

    /* Per-node arrays — extract by myList_nod2D[i]-1. */
    real_t *new_coord     = malloc((size_t)my_n * 2 * sizeof(real_t));
    real_t *new_geo_coord = malloc((size_t)my_n * 2 * sizeof(real_t));
    int    *new_coast     = malloc((size_t)my_n * sizeof(int));
    real_t *new_depth     = malloc((size_t)my_n * sizeof(real_t));
    int    *new_nlevels_n = malloc((size_t)my_n * sizeof(int));
    FESOM_CHECK(new_coord && new_geo_coord && new_coast && new_depth && new_nlevels_n,
                "scatter_mesh: per-node alloc");
    for (int i = 0; i < my_n; ++i) {
        int gid = p->myList_nod2D[i] - 1;
        new_coord    [2*i + 0] = m->coord_nod2D    [2*gid + 0];
        new_coord    [2*i + 1] = m->coord_nod2D    [2*gid + 1];
        new_geo_coord[2*i + 0] = m->geo_coord_nod2D[2*gid + 0];
        new_geo_coord[2*i + 1] = m->geo_coord_nod2D[2*gid + 1];
        new_coast    [i]       = m->coast_flag     [gid];
        new_depth    [i]       = m->depth          [gid];
        new_nlevels_n[i]       = m->nlevels_nod2D  [gid];
    }

    /* Per-elem arrays — Fortran convention (oce_mesh.F90):
     *   elem_nodes : sized myDim+eDim+eXDim — Fortran's elem2D_nodes is
     *                allocated for the full local extent, and ssh_stiff
     *                build at oce_ale.F90:1535 reads elem2D_nodes(:, el(i))
     *                for el(i) potentially in halo. Halo entries are filled
     *                directly from the broadcast global elem_nodes; their
     *                local node IDs are valid because (eDim_elem2D + eXDim
     *                touches at most 2 rings out, which is contained in
     *                myDim+eDim_nod2D for first ring; for outer entries
     *                whose vertex isn't in our local list we mark -1 and
     *                trust callers to skip).
     *   nlevels    : sized myDim+eDim+eXDim (line 943).
     */
    int    *new_elem_nodes = malloc((size_t)my_e * 3 * sizeof(int));
    int    *new_nlevels_e  = malloc((size_t)my_e * sizeof(int));
    FESOM_CHECK(new_elem_nodes && new_nlevels_e, "scatter_mesh: per-elem alloc");
    for (int i = 0; i < my_e; ++i) {
        int gid = p->myList_elem2D[i] - 1;
        new_nlevels_e[i] = m->nlevels[gid];
    }
    int halo_elem_unmappable = 0;
    for (int i = 0; i < my_e; ++i) {
        int gid = p->myList_elem2D[i] - 1;
        int interior = (i < p->myDim_elem2D);
        for (int k = 0; k < 3; ++k) {
            int g_node = m->elem_nodes[3*gid + k];   /* 0-based global node */
            int l_node = node_g2l[g_node];
            if (interior) {
                FESOM_CHECK(l_node >= 0,
                            "scatter_mesh: rank %d INTERIOR elem %d (gid=%d) references "
                            "global node %d which is NOT in local node list — "
                            "partition file inconsistency",
                            p->mype, i, gid + 1, g_node + 1);
            } else if (l_node < 0) {
                /* Halo element vertex not in our node list — only happens for
                 * eXDim elements 2+ rings out. SSH stiffness only reads first-
                 * ring halo (eDim) elements; mark and skip. */
                ++halo_elem_unmappable;
            }
            new_elem_nodes[3*i + k] = l_node;        /* may be -1 for outer halo */
        }
    }
    if (halo_elem_unmappable > 0 && p->mype == 0) {
        printf("[fesom_mesh] scatter: %d halo-element vertex refs unmappable "
               "(outer eXDim ring) — left as -1\n", halo_elem_unmappable);
    }

    /* Per-edge arrays — extract by myList_edge2D[i]-1; translate node ids to
     * local; translate edge_tri (global elem id) to local elem id. For the
     * latter we'd need an elem global→local hash too. */
    int *elem_g2l = malloc((size_t)g_elem2D * sizeof(int));
    FESOM_CHECK(elem_g2l, "scatter_mesh: elem_g2l alloc");
    for (int i = 0; i < g_elem2D; ++i) elem_g2l[i] = -1;
    for (int i = 0; i < my_e; ++i) {
        int gid = p->myList_elem2D[i] - 1;
        elem_g2l[gid] = i;
    }

    int *new_edges    = malloc((size_t)my_eg * 2 * sizeof(int));
    int *new_edge_tri = malloc((size_t)my_eg * 2 * sizeof(int));
    FESOM_CHECK(new_edges && new_edge_tri, "scatter_mesh: per-edge alloc");
    for (int i = 0; i < my_eg; ++i) {
        int gid = p->myList_edge2D[i] - 1;
        for (int k = 0; k < 2; ++k) {
            int g_node = m->edges[2*gid + k];      /* 0-based global node */
            int l_node = (g_node >= 0) ? node_g2l[g_node] : -1;
            FESOM_CHECK(l_node >= 0,
                        "scatter_mesh: rank %d local edge %d (gid=%d) end-node %d "
                        "(global %d) not in local node list",
                        p->mype, i, gid + 1, k, g_node + 1);
            new_edges[2*i + k] = l_node;

            int g_elem = m->edge_tri[2*gid + k];   /* may be -1 for boundary */
            new_edge_tri[2*i + k] = (g_elem >= 0) ? elem_g2l[g_elem] : -1;
            /* Halo elem may be missing on this rank; that's fine for boundary
             * checks but not for stencils — not a fatal error here, viscosity
             * code already handles -1 as "no neighbour". */
        }
    }

    free(node_g2l);
    free(elem_g2l);

    /* Swap in new local arrays, free the old global ones. */
    free(m->coord_nod2D);     m->coord_nod2D     = new_coord;
    free(m->geo_coord_nod2D); m->geo_coord_nod2D = new_geo_coord;
    free(m->coast_flag);      m->coast_flag      = new_coast;
    free(m->depth);           m->depth           = new_depth;
    free(m->nlevels_nod2D);   m->nlevels_nod2D   = new_nlevels_n;
    free(m->elem_nodes);      m->elem_nodes      = new_elem_nodes;
    free(m->nlevels);         m->nlevels         = new_nlevels_e;
    free(m->edges);           m->edges           = new_edges;
    free(m->edge_tri);        m->edge_tri        = new_edge_tri;

    /* Sync mesh's myDim_ / eDim_ fields from partit. */
    m->myDim_nod2D   = p->myDim_nod2D;
    m->eDim_nod2D    = p->eDim_nod2D;
    m->myDim_elem2D  = p->myDim_elem2D;
    m->eDim_elem2D   = p->eDim_elem2D;
    m->eXDim_elem2D  = p->eXDim_elem2D;
    m->myDim_edge2D  = p->myDim_edge2D;
    m->eDim_edge2D   = p->eDim_edge2D;
}

/*--- Public entry points ----------------------------------------------------*/

void fesom_mesh_read(fesom_mesh *m, const char *mesh_dir,
                     fesom_partit *partit)
{
    if (partit->mype == 0) {
        read_nod2d(m, mesh_dir);
        read_elem2d(m, mesh_dir);
        read_aux3d(m, mesh_dir);
        read_nlvls(m, mesh_dir);
        read_elvls(m, mesh_dir);
        read_edges(m, mesh_dir);
    }

    if (partit->npes > 1) {
        scatter_mesh(m, partit);
    } else {
        /* npes==1, synthesised partit: m already has global=local arrays from
         * the read_* functions; just sync the dim fields. */
        m->myDim_nod2D  = m->nod2D;   m->eDim_nod2D  = 0;
        m->myDim_elem2D = m->elem2D;  m->eDim_elem2D = 0;  m->eXDim_elem2D = 0;
        m->myDim_edge2D = m->edge2D;  m->eDim_edge2D = 0;
    }

    int swapped = orient_cw(m);
    int total_swapped = swapped;
    if (partit->npes > 1) {
        MPI_CHECK(MPI_Reduce(&swapped, &total_swapped, 1, MPI_INT, MPI_SUM, 0,
                             partit->MPI_COMM_FESOM));
    }
    int total_elem = m->elem2D;     /* global, all ranks have it */
    if (partit->mype == 0) {
        printf("[fesom_mesh] orient_cw: swapped %d / %d elements (CW after fix)\n",
               total_swapped, total_elem);
    }
}

void fesom_mesh_compute_metrics(fesom_mesh *m, fesom_partit *partit)
{
    build_nod_in_elem2D(m);
    compute_vertical_levels_aux(m);

    /* Fortran order (oce_mesh.F90:2213-2266):
     *   1. compute elem_area for myDim
     *   2. exchange_elem(elem_area)               <-- critical for boundary nodes
     *   3. build node area / areasvol from elem_area
     * Doing the exchange BEFORE step 3 ensures partition-boundary nodes
     * receive correct contributions from their off-rank surrounding elements.
     * Without this, areasvol on boundary nodes is underestimated → SSH mass
     * matrix term too small → CG converges to oversized eta → blow-up. */
    compute_elem_area_only(m);
    if (partit->npes > 1) {
        fesom_halo_exchange(m->elem_area, 2 /* ELEM2D */,      1, 1, partit);
        fesom_halo_exchange(m->elem_area, 4 /* ELEM2D_FULL */, 1, 1, partit);
    }
    compute_node_areas(m);

    compute_metric_and_coriolis(m);   /* fills elem_cos / metric_factor / coriolis */

    /* Halo exchange of remaining element-scalar fields (elem_area already done). */
    if (partit->npes > 1) {
        fesom_halo_exchange(m->elem_cos,       2, 1, 1, partit);
        fesom_halo_exchange(m->metric_factor,  2, 1, 1, partit);
        fesom_halo_exchange(m->coriolis,       2, 1, 1, partit);
        fesom_halo_exchange(m->elem_center_x,  2, 1, 1, partit);
        fesom_halo_exchange(m->elem_center_y,  2, 1, 1, partit);
        fesom_halo_exchange(m->elem_cos,       4, 1, 1, partit);
        fesom_halo_exchange(m->metric_factor,  4, 1, 1, partit);
        fesom_halo_exchange(m->coriolis,       4, 1, 1, partit);
        fesom_halo_exchange(m->elem_center_x,  4, 1, 1, partit);
        fesom_halo_exchange(m->elem_center_y,  4, 1, 1, partit);
    }

    compute_edges_geometry(m);        /* needs elem_cos */
    compute_gradient_sca(m);          /* needs elem_cos + elem_area */

    /* GM/Redi mesh prerequisites — needs areasvol + elem_area + nlevels. */
    compute_mesh_resolution(m, partit);
    compute_zbar_3d_n(m);

    /* ocean_area: local sum (computed in compute_node_areas) → MPI_Allreduce. */
    if (partit->npes > 1) {
        real_t local = m->ocean_area;
        MPI_CHECK(MPI_Allreduce(&local, &m->ocean_area, 1, MPI_DOUBLE,
                                MPI_SUM, partit->MPI_COMM_FESOM));
    }

    /* MFCT: per-edge up/down-wind triangle (oce_muscl_adv.F90:185-333). */
    fesom_mesh_build_edge_up_dn_tri(m, partit);
}

/*===========================================================================
 * edge_up_dn_tri — port of oce_muscl_adv.F90:185-333. For each owned edge,
 * find the upwind triangle (around node1, containing −x) and the downwind
 * triangle (around node2, containing +x), where x = coord(n2)−coord(n1).
 * "Containing" = the edge vector lies between the triangle's two edge vectors
 * (b,c) emanating from the shared node (atan2 angle test). Uses coord_elem +
 * e_nodes built on the WIDE element halo so 2-ring triangles are covered.
 *===========================================================================*/

/* Is the edge vector (xx,xy) between the triangle edge vectors b=(bx,by),
 * c=(cx,cy)? Decompose along c and 90°-CCW, compare atan2 angles (F:268-288). */
static int edge_x_between_bc(real_t bx, real_t by, real_t cx, real_t cy,
                             real_t xx, real_t xy)
{
    real_t cr = cx*cx + cy*cy;
    if (cr == 0.0) return 0;
    real_t bxp = ( bx*cx + by*cy)/cr;
    real_t byp = (-bx*cy + by*cx)/cr;
    real_t xxp = ( xx*cx + xy*cy)/cr;
    real_t xyp = (-xx*cy + xy*cx)/cr;
    real_t ab = atan2(byp, bxp);
    real_t ax = atan2(xyp, xxp);
    if (ab > 0.0 && ax > 0.0 && ax < ab) return 1;
    if (ab < 0.0 && ax < 0.0 && ax > ab) return 1;
    if (ab == ax || ax == 0.0) return 1;
    return 0;
}

/* Scan the triangles around `node` for the one whose two edge vectors from the
 * shared node bracket (xx,xy). ce = element node coords [E*6], en = element node
 * global ids [E*3]. Returns the 0-based elem id or -1. */
static int edge_up_dn_scan(const fesom_mesh *m, const fesom_partit *partit,
                           const real_t *ce, const real_t *en,
                           int node, real_t xx, real_t xy)
{
    const real_t cyc = (real_t)FESOM_CYCLIC_LENGTH_RAD, hcyc = 0.5*cyc;
    long long gnode = (long long)partit->myList_nod2D[node];   /* shared node global id */
    int b0 = m->nod_in_elem2D_offsets[node], b1 = m->nod_in_elem2D_offsets[node + 1];
    for (int k = b0; k < b1; ++k) {
        int el = m->nod_in_elem2D[k];
        int js = -1;
        for (int j = 0; j < 3; ++j)
            if ((long long)(en[(size_t)el*3 + j] + 0.5) == gnode) { js = j; break; }
        if (js < 0) continue;
        int jb = (js == 0) ? 1 : 0;     /* the two non-shared nodes, Fortran order */
        int jc = (js == 2) ? 1 : 2;
        real_t bx = ce[(size_t)el*6 + 2*jb + 0] - ce[(size_t)el*6 + 2*js + 0];
        real_t by = ce[(size_t)el*6 + 2*jb + 1] - ce[(size_t)el*6 + 2*js + 1];
        real_t cx = ce[(size_t)el*6 + 2*jc + 0] - ce[(size_t)el*6 + 2*js + 0];
        real_t cy = ce[(size_t)el*6 + 2*jc + 1] - ce[(size_t)el*6 + 2*js + 1];
        if (bx >  hcyc) bx -= cyc;  if (bx < -hcyc) bx += cyc;
        if (cx >  hcyc) cx -= cyc;  if (cx < -hcyc) cx += cyc;
        if (edge_x_between_bc(bx, by, cx, cy, xx, xy)) return el;
    }
    return -1;
}

void fesom_mesh_build_edge_up_dn_tri(fesom_mesh *m, fesom_partit *partit)
{
    const int Eme = m->myDim_edge2D;
    const int E   = m->myDim_elem2D + m->eDim_elem2D + m->eXDim_elem2D;
    const real_t cyc = (real_t)FESOM_CYCLIC_LENGTH_RAD, hcyc = 0.5*cyc;

    m->edge_up_dn_tri = malloc((size_t)Eme * 2 * sizeof(int));
    FESOM_CHECK(m->edge_up_dn_tri, "edge_up_dn_tri alloc");
    for (int i = 0; i < Eme * 2; ++i) m->edge_up_dn_tri[i] = -1;

    /* coord_elem [E*6] + e_nodes (global ids) [E*3], on the wide element halo
     * (F:192-235: exchange_elem of the element node coords + ids). */
    real_t *ce = calloc((size_t)E * 6, sizeof(real_t));
    real_t *en = calloc((size_t)E * 3, sizeof(real_t));
    FESOM_CHECK(ce && en, "edge_up_dn_tri scratch alloc");
    for (int el = 0; el < m->myDim_elem2D; ++el)
        for (int j = 0; j < 3; ++j) {
            int nd = m->elem_nodes[3*el + j];
            ce[(size_t)el*6 + 2*j + 0] = m->coord_nod2D[2*nd + 0];
            ce[(size_t)el*6 + 2*j + 1] = m->coord_nod2D[2*nd + 1];
            en[(size_t)el*3 + j]       = (real_t)partit->myList_nod2D[nd];
        }
    if (partit && partit->npes > 1) {
        fesom_halo_exchange(ce, FESOM_HALO_ELEM2D_FULL, 6, 1, partit);
        fesom_halo_exchange(en, FESOM_HALO_ELEM2D_FULL, 3, 1, partit);
    }

    for (int ed = 0; ed < Eme; ++ed) {
        int n1 = m->edges[2*ed + 0];
        int n2 = m->edges[2*ed + 1];
        real_t x0 = m->coord_nod2D[2*n2 + 0] - m->coord_nod2D[2*n1 + 0];
        real_t x1 = m->coord_nod2D[2*n2 + 1] - m->coord_nod2D[2*n1 + 1];
        if (x0 >  hcyc) x0 -= cyc;  if (x0 < -hcyc) x0 += cyc;
        m->edge_up_dn_tri[2*ed + 0] = edge_up_dn_scan(m, partit, ce, en, n1, -x0, -x1);
        m->edge_up_dn_tri[2*ed + 1] = edge_up_dn_scan(m, partit, ce, en, n2,  x0,  x1);
    }
    free(ce); free(en);
}
