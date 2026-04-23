#include "fesom_mesh.h"
#include "fesom_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    free(m->nlevels);
    free(m->ulevels);
    free(m->zbar);
    free(m->Z);
    free(m->depth);
    free(m->nod_in_elem2D_offsets);
    free(m->nod_in_elem2D);
    free(m->elem_area);
    free(m->area);
    free(m->areasvol);
    free(m->elem_cos);
    free(m->metric_factor);
    free(m->coriolis);
    free(m->coriolis_node);
    free(m->edge_dxdy);
    free(m->edge_cross_dxdy);
    free(m->gradient_sca);
    free(m->hnode);
    free(m->hnode_new);
    free(m->helem);
    free(m->hbar);
    free(m->hbar_old);
    memset(m, 0, sizeof(*m));
}

void fesom_mesh_alloc_state(fesom_mesh *m)
{
    size_t n_nl = (size_t)m->nod2D  * (size_t)m->nl;
    size_t e_nl = (size_t)m->elem2D * (size_t)m->nl;
    m->hnode     = calloc(n_nl, sizeof(real_t));
    m->hnode_new = calloc(n_nl, sizeof(real_t));
    m->helem     = calloc(e_nl, sizeof(real_t));
    m->hbar      = calloc((size_t)m->nod2D, sizeof(real_t));
    m->hbar_old  = calloc((size_t)m->nod2D, sizeof(real_t));
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
    for (int e = 0; e < m->elem2D; ++e) {
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
    const int N = m->nod2D;
    m->nod_in_elem2D_offsets = calloc((size_t)N + 1, sizeof(int));
    FESOM_CHECK(m->nod_in_elem2D_offsets, "nod_in_elem2D_offsets: out of memory");

    /* Pass 1: count per node into offsets[1..N], prefix-sum after. */
    for (int e = 0; e < m->elem2D; ++e) {
        for (int k = 0; k < 3; ++k) {
            int v = m->elem_nodes[3*e + k];
            ++m->nod_in_elem2D_offsets[v + 1];
        }
    }
    for (int i = 1; i <= N; ++i) {
        m->nod_in_elem2D_offsets[i] += m->nod_in_elem2D_offsets[i - 1];
    }
    int total = m->nod_in_elem2D_offsets[N];
    FESOM_CHECK(total == 3 * m->elem2D, "nod_in_elem2D: total mismatch %d vs %d",
                total, 3 * m->elem2D);

    /* Pass 2: scatter using a temporary cursor per node. */
    int *cursor = calloc((size_t)N, sizeof(int));
    FESOM_CHECK(cursor, "nod_in_elem2D: cursor OOM");
    m->nod_in_elem2D = malloc((size_t)total * sizeof(int));
    FESOM_CHECK(m->nod_in_elem2D, "nod_in_elem2D: out of memory");

    for (int e = 0; e < m->elem2D; ++e) {
        for (int k = 0; k < 3; ++k) {
            int v = m->elem_nodes[3*e + k];
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
    m->ulevels      = malloc((size_t)m->elem2D * sizeof(int));
    m->ulevels_nod2D    = malloc((size_t)m->nod2D  * sizeof(int));
    m->nlevels_nod2D_min = malloc((size_t)m->nod2D  * sizeof(int));
    FESOM_CHECK(m->ulevels && m->ulevels_nod2D && m->nlevels_nod2D_min,
                "vertical aux: out of memory");

    for (int e = 0; e < m->elem2D; ++e) m->ulevels[e]      = 1;
    for (int n = 0; n < m->nod2D;  ++n) m->ulevels_nod2D[n] = 1;

    for (int n = 0; n < m->nod2D; ++n) {
        int o0 = m->nod_in_elem2D_offsets[n];
        int o1 = m->nod_in_elem2D_offsets[n + 1];
        FESOM_CHECK(o1 > o0, "node %d has no surrounding elements", n);
        int mn = m->nlevels[m->nod_in_elem2D[o0]];
        for (int k = o0 + 1; k < o1; ++k) {
            int v = m->nlevels[m->nod_in_elem2D[k]];
            if (v < mn) mn = v;
        }
        m->nlevels_nod2D_min[n] = mn;
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
static void compute_areas(fesom_mesh *m)
{
    const real_t cyc      = FESOM_CYCLIC_LENGTH_RAD;
    const real_t half_cyc = 0.5 * cyc;
    const real_t r2       = (real_t)FESOM_R_EARTH * (real_t)FESOM_R_EARTH;

    m->elem_area = malloc((size_t)m->elem2D * sizeof(real_t));
    m->area      = calloc((size_t)m->nod2D * (size_t)m->nl, sizeof(real_t));
    m->areasvol  = calloc((size_t)m->nod2D * (size_t)m->nl, sizeof(real_t));
    FESOM_CHECK(m->elem_area && m->area && m->areasvol,
                "areas: out of memory");

    /* per-cell area in radian² */
    for (int e = 0; e < m->elem2D; ++e) {
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

    /* distribute cell area to surrounding nodes per layer (no cavity) */
    for (int n = 0; n < m->nod2D; ++n) {
        int o0 = m->nod_in_elem2D_offsets[n];
        int o1 = m->nod_in_elem2D_offsets[n + 1];
        for (int k = o0; k < o1; ++k) {
            int e = m->nod_in_elem2D[k];
            int nzmin = m->ulevels[e] - 1;          /* 1-based → 0-based */
            int nzmax = m->nlevels[e] - 1;          /* exclusive bound = nlevels-1 levels */
            real_t third = m->elem_area[e] / 3.0;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                m->area[FESOM_NODE3D(n, nz, m->nl)] += third;
            }
        }
        /* non-cavity: areasvol == area */
        int nzmin_n = m->ulevels_nod2D[n] - 1;
        int nzmax_n = m->nlevels_nod2D[n] - 1;
        for (int nz = nzmin_n; nz < nzmax_n; ++nz) {
            m->areasvol[FESOM_NODE3D(n, nz, m->nl)] = m->area[FESOM_NODE3D(n, nz, m->nl)];
        }
    }

    /* radian² → m² */
    for (int e = 0; e < m->elem2D; ++e) m->elem_area[e] *= r2;
    size_t tot = (size_t)m->nod2D * (size_t)m->nl;
    for (size_t i = 0; i < tot; ++i) {
        m->area[i]     *= r2;
        m->areasvol[i] *= r2;
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

    m->elem_cos      = malloc((size_t)m->elem2D * sizeof(real_t));
    m->metric_factor = malloc((size_t)m->elem2D * sizeof(real_t));
    m->coriolis      = malloc((size_t)m->elem2D * sizeof(real_t));
    m->coriolis_node = malloc((size_t)m->nod2D  * sizeof(real_t));
    FESOM_CHECK(m->elem_cos && m->metric_factor && m->coriolis && m->coriolis_node,
                "metric/coriolis: out of memory");

    for (int n = 0; n < m->nod2D; ++n) {
        m->coriolis_node[n] =
            2.0 * (real_t)FESOM_OMEGA * sin(m->geo_coord_nod2D[2*n + 1]);
    }

    for (int e = 0; e < m->elem2D; ++e) {
        real_t cx, cy;
        elem_center_xy(m, e, &cx, &cy);
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

    m->edge_dxdy       = malloc((size_t)m->edge2D * 2 * sizeof(real_t));
    m->edge_cross_dxdy = malloc((size_t)m->edge2D * 4 * sizeof(real_t));
    FESOM_CHECK(m->edge_dxdy && m->edge_cross_dxdy, "edges geometry: out of memory");

    /* Cell centroids cached so we don't recompute per-edge. */
    real_t *cx = malloc((size_t)m->elem2D * sizeof(real_t));
    real_t *cy = malloc((size_t)m->elem2D * sizeof(real_t));
    FESOM_CHECK(cx && cy, "elem-center cache: out of memory");
    for (int e = 0; e < m->elem2D; ++e) elem_center_xy(m, e, &cx[e], &cy[e]);

    for (int eg = 0; eg < m->edge2D; ++eg) {
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

    free(cx);
    free(cy);
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

    m->gradient_sca = malloc((size_t)m->elem2D * 6 * sizeof(real_t));
    FESOM_CHECK(m->gradient_sca, "gradient_sca: out of memory");

    for (int e = 0; e < m->elem2D; ++e) {
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

/*--- Public entry points ----------------------------------------------------*/

void fesom_mesh_read(fesom_mesh *m, const char *mesh_dir)
{
    read_nod2d(m, mesh_dir);
    read_elem2d(m, mesh_dir);
    read_aux3d(m, mesh_dir);
    read_nlvls(m, mesh_dir);
    read_elvls(m, mesh_dir);
    read_edges(m, mesh_dir);

    int swapped = orient_cw(m);
    printf("[fesom_mesh] orient_cw: swapped %d / %d elements (CW after fix)\n",
           swapped, m->elem2D);
}

void fesom_mesh_compute_metrics(fesom_mesh *m)
{
    build_nod_in_elem2D(m);
    compute_vertical_levels_aux(m);
    compute_areas(m);                 /* needs nod_in_elem2D */
    compute_metric_and_coriolis(m);   /* needs nothing extra */
    compute_edges_geometry(m);        /* needs elem_cos */
    compute_gradient_sca(m);          /* needs elem_cos + elem_area */
}
