/*
 * Unit tests for the restart file layer (src/fesom_io_restart.c).
 *
 * Builds a six-node, four-element, three-level mesh by hand — enough for the
 * gather/scatter path and the netCDF layout, small enough to reason about —
 * writes a restart, and reads it back into a second set of arrays.
 *
 * The guard cases (wrong format version, wrong mesh, missing variable,
 * truncated file) end in MPI_Abort by design, so each is run in a child
 * process and judged by its exit status. `test_restart_format <case>` runs one
 * such case directly; with no argument the binary runs the whole suite.
 *
 * Build target `test_restart_format`, run via `ctest -V`.
 */
#include "fesom_dyn.h"
#include "fesom_ice_types.h"
#include "fesom_io_restart.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_ssh.h"
#include "fesom_tracers.h"

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NOD  6
#define ELE  4
#define NL   3

static int n_passed = 0;
static int n_failed = 0;

/* The build is -DNDEBUG, so assert() would compile the call inside it away.
 * MUST always evaluates its argument. */
#define MUST(expr)                                                           \
    do {                                                                     \
        if (!(expr)) {                                                       \
            fprintf(stderr, "fatal %s:%d: %s\n", __FILE__, __LINE__, #expr); \
            exit(2);                                                         \
        }                                                                    \
    } while (0)

#define CHECK(expr)                                                          \
    do {                                                                     \
        if (expr) { ++n_passed; }                                            \
        else      {                                                          \
            ++n_failed;                                                      \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr);  \
        }                                                                    \
    } while (0)

/* nlevels_nod2D for the toy mesh; its sum is the mesh signature the reader
 * checks, so it must differ between the two meshes built below. */
static const int NLEV_A[NOD] = { 3, 3, 2, 3, 2, 3 };   /* sum 16 */
static const int NLEV_B[NOD] = { 3, 3, 3, 3, 2, 3 };   /* sum 17 */

typedef struct {
    fesom_mesh    mesh;
    fesom_partit  partit;
    fesom_dyn     dyn;
    fesom_tracers tracers;
    fesom_ice     ice;
    fesom_ssh_stiff stiff;
} model;

static real_t *rbuf(size_t n) { real_t *p = calloc(n, sizeof(real_t)); MUST(p); return p; }

/* A model whose every restart-carried array is filled with a distinct,
 * exactly representable pattern, so a swapped or transposed field shows up as
 * a mismatch rather than as plausible numbers. */
static void model_init(model *m, const int *nlev, int fill)
{
    memset(m, 0, sizeof(*m));

    m->partit.npes = 1;
    m->partit.mype = 0;
    m->partit.MPI_COMM_FESOM = MPI_COMM_WORLD;
    m->partit.myDim_nod2D  = NOD;
    m->partit.myDim_elem2D = ELE;
    m->partit.myList_nod2D  = malloc(NOD * sizeof(int));
    m->partit.myList_elem2D = malloc(ELE * sizeof(int));
    MUST(m->partit.myList_nod2D && m->partit.myList_elem2D);
    for (int i = 0; i < NOD; ++i) m->partit.myList_nod2D[i]  = i + 1;  /* 1-based */
    for (int i = 0; i < ELE; ++i) m->partit.myList_elem2D[i] = i + 1;

    m->mesh.nod2D  = NOD;
    m->mesh.elem2D = ELE;
    m->mesh.nl     = NL;
    m->mesh.myDim_nod2D  = NOD; m->mesh.eDim_nod2D  = 0;
    m->mesh.myDim_elem2D = ELE; m->mesh.eDim_elem2D = 0; m->mesh.eXDim_elem2D = 0;
    m->mesh.nlevels_nod2D = malloc(NOD * sizeof(int));
    MUST(m->mesh.nlevels_nod2D);
    memcpy(m->mesh.nlevels_nod2D, nlev, NOD * sizeof(int));

    m->mesh.hbar      = rbuf(NOD);
    m->mesh.hbar_old  = rbuf(NOD);
    m->mesh.hnode     = rbuf((size_t)NOD * NL);
    m->mesh.hnode_new = rbuf((size_t)NOD * NL);
    m->mesh.zbar_3d_n = rbuf((size_t)NOD * NL);
    m->mesh.Z_3d_n    = rbuf((size_t)NOD * NL);
    m->mesh.helem     = rbuf((size_t)ELE * NL);
    m->mesh.dhe       = rbuf(ELE);

    m->dyn.eta_n       = rbuf(NOD);
    m->dyn.d_eta       = rbuf(NOD);
    m->dyn.ssh_rhs_old = rbuf(NOD);
    m->dyn.w           = rbuf((size_t)NOD * NL);
    m->dyn.w_e         = rbuf((size_t)NOD * NL);
    m->dyn.w_i         = rbuf((size_t)NOD * NL);
    m->dyn.uvnode      = rbuf((size_t)NOD * NL * 2);
    m->dyn.uv          = rbuf((size_t)ELE * NL * 2);
    m->dyn.uv_rhsAB    = rbuf((size_t)ELE * NL * 2);

    m->tracers.num_tracers = 2;
    for (int t = 0; t < 2; ++t) {
        m->tracers.data[t].values    = rbuf((size_t)NOD * NL);
        m->tracers.data[t].valuesAB  = rbuf((size_t)NOD * NL);
        m->tracers.data[t].valuesold = rbuf((size_t)NOD * NL);
    }

    for (int i = 0; i < FESOM_NUM_ICE_TRACERS; ++i)
        m->ice.data[i].values = rbuf(NOD);
    m->ice.uice = rbuf(NOD);
    m->ice.vice = rbuf(NOD);
    m->ice.thermo.t_skin = rbuf(NOD);
    m->ice.work.sigma11  = rbuf(ELE);
    m->ice.work.sigma12  = rbuf(ELE);
    m->ice.work.sigma22  = rbuf(ELE);
    m->ice.ice_steps_since_upd = 7;
    m->ice.ice_update          = 1;

    /* A stiffness matrix with rows of unequal length and columns deliberately
     * NOT in global order, so the sort in stiff_row_order is exercised. */
    m->stiff.dim = NOD;
    m->stiff.rowptr = malloc((NOD + 1) * sizeof(int));
    MUST(m->stiff.rowptr);
    static const int rowlen[NOD] = { 3, 2, 4, 1, 3, 2 };
    m->stiff.rowptr[0] = 0;
    for (int i = 0; i < NOD; ++i) m->stiff.rowptr[i + 1] = m->stiff.rowptr[i] + rowlen[i];
    m->stiff.nnz = m->stiff.rowptr[NOD];
    m->stiff.colind = malloc((size_t)m->stiff.nnz * sizeof(int));
    m->stiff.values = rbuf((size_t)m->stiff.nnz);
    MUST(m->stiff.colind);
    for (int r = 0, k = 0; r < NOD; ++r)
        for (int i = 0; i < rowlen[r]; ++i, ++k)
            m->stiff.colind[k] = (r + 3 * (rowlen[r] - i)) % NOD;   /* descending-ish */

    if (!fill) return;

    for (int k = 0; k < m->stiff.nnz; ++k) m->stiff.values[k] = (2000 + k) / 64.0;

    /* Every value is k/64 for an integer k, so equality is exact and a
     * mis-shuffle is visible by eye in ncdump. */
    for (int n = 0; n < NOD; ++n) {
        m->mesh.hbar[n]         = (n + 1) / 64.0;
        m->mesh.hbar_old[n]     = (n + 11) / 64.0;
        m->dyn.eta_n[n]         = (n + 21) / 64.0;
        m->dyn.d_eta[n]         = (n + 81) / 64.0;
        m->dyn.ssh_rhs_old[n]   = (n + 31) / 64.0;
        m->ice.uice[n]          = (n + 41) / 64.0;
        m->ice.vice[n]          = (n + 51) / 64.0;
        m->ice.thermo.t_skin[n] = (n + 61) / 64.0;
        for (int i = 0; i < FESOM_NUM_ICE_TRACERS; ++i)
            m->ice.data[i].values[n] = (n + 71 + 10 * i) / 64.0;
        for (int z = 0; z < NL; ++z) {
            size_t k = (size_t)n * NL + z;
            m->mesh.hnode[k]     = (100 + 10 * n + z) / 64.0;
            m->mesh.zbar_3d_n[k] = (200 + 10 * n + z) / 64.0;
            m->mesh.Z_3d_n[k]    = (300 + 10 * n + z) / 64.0;
            m->dyn.w[k]          = (400 + 10 * n + z) / 64.0;
            m->dyn.w_e[k]        = (500 + 10 * n + z) / 64.0;
            m->dyn.w_i[k]        = (600 + 10 * n + z) / 64.0;
            for (int c = 0; c < 2; ++c)
                m->dyn.uvnode[(size_t)n * NL * 2 + (size_t)z * 2 + c] =
                    (1800 + 100 * c + 10 * n + z) / 64.0;
            for (int t = 0; t < 2; ++t) {
                m->tracers.data[t].values[k]    = (700 + 100 * t + 10 * n + z) / 64.0;
                m->tracers.data[t].valuesAB[k]  = (900 + 100 * t + 10 * n + z) / 64.0;
                m->tracers.data[t].valuesold[k] = (1100 + 100 * t + 10 * n + z) / 64.0;
            }
        }
    }
    for (int e = 0; e < ELE; ++e) {
        m->mesh.dhe[e]          = (e + 3) / 64.0;
        m->ice.work.sigma11[e]  = (e + 13) / 64.0;
        m->ice.work.sigma12[e]  = (e + 23) / 64.0;
        m->ice.work.sigma22[e]  = (e + 33) / 64.0;
        for (int z = 0; z < NL; ++z) {
            m->mesh.helem[(size_t)e * NL + z] = (1300 + 10 * e + z) / 64.0;
            for (int c = 0; c < 2; ++c) {
                size_t k = (size_t)e * NL * 2 + (size_t)z * 2 + c;
                m->dyn.uv[k]       = (1400 + 100 * c + 10 * e + z) / 64.0;
                m->dyn.uv_rhsAB[k] = (1600 + 100 * c + 10 * e + z) / 64.0;
            }
        }
    }
}

static void model_free(model *m)
{
    free(m->partit.myList_nod2D); free(m->partit.myList_elem2D);
    free(m->mesh.nlevels_nod2D);
    free(m->mesh.hbar); free(m->mesh.hbar_old); free(m->mesh.hnode);
    free(m->mesh.hnode_new); free(m->mesh.zbar_3d_n); free(m->mesh.Z_3d_n);
    free(m->mesh.helem); free(m->mesh.dhe);
    free(m->dyn.eta_n); free(m->dyn.d_eta); free(m->dyn.ssh_rhs_old);
    free(m->dyn.w); free(m->dyn.w_e); free(m->dyn.w_i);
    free(m->dyn.uv); free(m->dyn.uv_rhsAB); free(m->dyn.uvnode);
    for (int t = 0; t < 2; ++t) {
        free(m->tracers.data[t].values);
        free(m->tracers.data[t].valuesAB);
        free(m->tracers.data[t].valuesold);
    }
    for (int i = 0; i < FESOM_NUM_ICE_TRACERS; ++i) free(m->ice.data[i].values);
    free(m->ice.uice); free(m->ice.vice); free(m->ice.thermo.t_skin);
    free(m->ice.work.sigma11); free(m->ice.work.sigma12); free(m->ice.work.sigma22);
    free(m->stiff.rowptr); free(m->stiff.colind); free(m->stiff.values);
}

static fesom_calendar_t a_calendar(void)
{
    fesom_calendar_t c;
    memset(&c, 0, sizeof(c));
    c.kind = FESOM_CAL_GREGORIAN;
    c.year = 1963; c.month = 7; c.day = 19;
    c.hour = 13; c.minute = 30; c.second = 0.0;
    c.origin_year = 1958; c.origin_month = 1; c.origin_day = 1;
    return c;
}

static void write_reference(const char *path, model *m)
{
    fesom_calendar_t cal = a_calendar();
    fesom_restart_write(path, 4242, 1800.0, &cal, &m->mesh, &m->dyn,
                        &m->tracers, &m->ice, NULL, &m->stiff, &m->partit);
}

static int arrays_equal(const real_t *a, const real_t *b, size_t n)
{
    for (size_t i = 0; i < n; ++i) if (a[i] != b[i]) return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
static void test_round_trip(const char *dir)
{
    char path[512];
    snprintf(path, sizeof path, "%s/roundtrip.restart.nc", dir);

    model src, dst;
    model_init(&src, NLEV_A, /*fill=*/1);
    model_init(&dst, NLEV_A, /*fill=*/0);
    write_reference(path, &src);

    int step = 0; double dt = 0.0;
    fesom_calendar_t cal;
    memset(&cal, 0, sizeof(cal));
    fesom_restart_read(path, &step, &dt, &cal, &dst.mesh, &dst.dyn,
                       &dst.tracers, &dst.ice, NULL, &dst.stiff, &dst.partit);

    CHECK(step == 4242);
    CHECK(dt == 1800.0);
    CHECK(cal.year == 1963 && cal.month == 7 && cal.day == 19);
    CHECK(cal.hour == 13 && cal.minute == 30 && cal.second == 0.0);
    CHECK(cal.origin_year == 1958 && cal.origin_month == 1 && cal.origin_day == 1);
    CHECK(cal.kind == FESOM_CAL_GREGORIAN);
    CHECK(dst.ice.ice_steps_since_upd == 7 && dst.ice.ice_update == 1);

    /* Every carried array, bit for bit. */
    CHECK(arrays_equal(src.mesh.hbar,      dst.mesh.hbar,      NOD));
    CHECK(arrays_equal(src.mesh.hbar_old,  dst.mesh.hbar_old,  NOD));
    CHECK(arrays_equal(src.dyn.eta_n,      dst.dyn.eta_n,      NOD));
    CHECK(arrays_equal(src.dyn.d_eta,      dst.dyn.d_eta,      NOD));
    CHECK(arrays_equal(src.dyn.ssh_rhs_old,dst.dyn.ssh_rhs_old,NOD));
    CHECK(arrays_equal(src.ice.uice,       dst.ice.uice,       NOD));
    CHECK(arrays_equal(src.ice.vice,       dst.ice.vice,       NOD));
    CHECK(arrays_equal(src.ice.thermo.t_skin, dst.ice.thermo.t_skin, NOD));
    for (int i = 0; i < FESOM_NUM_ICE_TRACERS; ++i)
        CHECK(arrays_equal(src.ice.data[i].values, dst.ice.data[i].values, NOD));

    CHECK(arrays_equal(src.mesh.hnode,     dst.mesh.hnode,     (size_t)NOD * NL));
    CHECK(arrays_equal(src.mesh.zbar_3d_n, dst.mesh.zbar_3d_n, (size_t)NOD * NL));
    CHECK(arrays_equal(src.mesh.Z_3d_n,    dst.mesh.Z_3d_n,    (size_t)NOD * NL));
    CHECK(arrays_equal(src.dyn.w,          dst.dyn.w,          (size_t)NOD * NL));
    CHECK(arrays_equal(src.dyn.w_e,        dst.dyn.w_e,        (size_t)NOD * NL));
    CHECK(arrays_equal(src.dyn.w_i,        dst.dyn.w_i,        (size_t)NOD * NL));
    CHECK(arrays_equal(src.dyn.uvnode,     dst.dyn.uvnode,     (size_t)NOD * NL * 2));
    for (int t = 0; t < 2; ++t) {
        CHECK(arrays_equal(src.tracers.data[t].values,    dst.tracers.data[t].values,    (size_t)NOD * NL));
        CHECK(arrays_equal(src.tracers.data[t].valuesAB,  dst.tracers.data[t].valuesAB,  (size_t)NOD * NL));
        CHECK(arrays_equal(src.tracers.data[t].valuesold, dst.tracers.data[t].valuesold, (size_t)NOD * NL));
    }

    CHECK(arrays_equal(src.mesh.dhe,         dst.mesh.dhe,         ELE));
    CHECK(arrays_equal(src.ice.work.sigma11, dst.ice.work.sigma11, ELE));
    CHECK(arrays_equal(src.ice.work.sigma12, dst.ice.work.sigma12, ELE));
    CHECK(arrays_equal(src.ice.work.sigma22, dst.ice.work.sigma22, ELE));
    CHECK(arrays_equal(src.mesh.helem,    dst.mesh.helem,    (size_t)ELE * NL));
    CHECK(arrays_equal(src.dyn.uv,        dst.dyn.uv,        (size_t)ELE * NL * 2));
    CHECK(arrays_equal(src.dyn.uv_rhsAB,  dst.dyn.uv_rhsAB,  (size_t)ELE * NL * 2));

    /* The accumulated stiffness matrix travels with its global column ids, so
     * the values land back on the same entries however the rows are ordered. */
    CHECK(arrays_equal(src.stiff.values, dst.stiff.values, (size_t)src.stiff.nnz));

    /* hnode_new is not in the file; the reader mirrors hnode into it. */
    CHECK(arrays_equal(dst.mesh.hnode, dst.mesh.hnode_new, (size_t)NOD * NL));

    model_free(&src); model_free(&dst);
}

/* The two velocity components must not be swapped: u and v come from one
 * interleaved array and are written as two variables. */
static void test_components_are_not_swapped(const char *dir)
{
    char path[512];
    snprintf(path, sizeof path, "%s/comp.restart.nc", dir);
    model src;
    model_init(&src, NLEV_A, 1);
    write_reference(path, &src);

    int ncid, vu, vv;
    double u[ELE * NL], v[ELE * NL];
    MUST(nc_open(path, NC_NOWRITE, &ncid) == NC_NOERR);
    MUST(nc_inq_varid(ncid, "u", &vu) == NC_NOERR);
    MUST(nc_inq_varid(ncid, "v", &vv) == NC_NOERR);
    MUST(nc_get_var_double(ncid, vu, u) == NC_NOERR);
    MUST(nc_get_var_double(ncid, vv, v) == NC_NOERR);
    nc_close(ncid);

    /* File layout is (nz, elem); the source is [elem][nz][comp]. */
    for (int e = 0; e < ELE; ++e)
        for (int z = 0; z < NL; ++z) {
            CHECK(u[z * ELE + e] == src.dyn.uv[(size_t)e * NL * 2 + (size_t)z * 2 + 0]);
            CHECK(v[z * ELE + e] == src.dyn.uv[(size_t)e * NL * 2 + (size_t)z * 2 + 1]);
        }
    model_free(&src);
}

static void test_tke_is_optional(const char *dir)
{
    char path[512];
    snprintf(path, sizeof path, "%s/notke.restart.nc", dir);
    model src;
    model_init(&src, NLEV_A, 1);
    write_reference(path, &src);            /* tke = NULL */

    int ncid, vid;
    MUST(nc_open(path, NC_NOWRITE, &ncid) == NC_NOERR);
    CHECK(nc_inq_varid(ncid, "tke", &vid) == NC_ENOTVAR);   /* absent, as in Fortran */
    CHECK(nc_inq_varid(ncid, "temp", &vid) == NC_NOERR);
    nc_close(ncid);
    model_free(&src);
}

static void test_config_and_cadence(void)
{
    fesom_restart_cfg cfg;

    unsetenv("FESOM_RESTART_OUT"); unsetenv("FESOM_RESTART_IN");
    unsetenv("FESOM_RESTART_EVERY");
    fesom_restart_config(&cfg);
    CHECK(cfg.out_dir == NULL && cfg.in_path == NULL && cfg.every == 0);
    CHECK(fesom_restart_due(&cfg, 10, 10) == 0);     /* writing off entirely */

    setenv("FESOM_RESTART_OUT", "/tmp/x", 1);
    setenv("FESOM_RESTART_EVERY", "5", 1);
    fesom_restart_config(&cfg);
    CHECK(cfg.every == 5);
    CHECK(fesom_restart_due(&cfg, 5, 20) == 1);
    CHECK(fesom_restart_due(&cfg, 6, 20) == 0);
    CHECK(fesom_restart_due(&cfg, 20, 20) == 1);     /* last step always */

    setenv("FESOM_RESTART_EVERY", "0", 1);
    fesom_restart_config(&cfg);
    CHECK(fesom_restart_due(&cfg, 5, 20) == 0);
    CHECK(fesom_restart_due(&cfg, 20, 20) == 1);
    unsetenv("FESOM_RESTART_OUT"); unsetenv("FESOM_RESTART_EVERY");

    fesom_calendar_t cal = a_calendar();
    char path[256];
    fesom_restart_path(path, sizeof path, "/work/run", &cal);
    CHECK(strcmp(path, "/work/run/fesom.19630719_48600.restart.nc") == 0);
}

/* ------------------------------------------------------------------ *
 * Guard cases — each must abort. Run as a child process.              *
 * ------------------------------------------------------------------ */
static void run_guard_case(const char *which)
{
    char dir[] = "/tmp/fesom_rst_guard_XXXXXX";
    MUST(mkdtemp(dir));
    char path[512];
    snprintf(path, sizeof path, "%s/g.restart.nc", dir);

    model src;
    model_init(&src, NLEV_A, 1);
    write_reference(path, &src);

    if (strcmp(which, "version") == 0) {
        int ncid, bogus = FESOM_RESTART_FORMAT_VERSION + 99;
        MUST(nc_open(path, NC_WRITE, &ncid) == NC_NOERR);
        MUST(nc_put_att_int(ncid, NC_GLOBAL, "restart_format_version",
                              NC_INT, 1, &bogus) == NC_NOERR);
        nc_close(ncid);
    } else if (strcmp(which, "missing_var") == 0) {
        /* Rewrite the file without `salt`: netCDF cannot delete a variable, so
         * the fixture is a fresh file that simply never defines it. */
        int ncid, dn, de, dz;
        MUST(nc_create(path, NC_CLOBBER | NC_64BIT_DATA, &ncid) == NC_NOERR);
        nc_def_dim(ncid, "node", NOD, &dn);
        nc_def_dim(ncid, "elem", ELE, &de);
        nc_def_dim(ncid, "nz",   NL,  &dz);
        int v = FESOM_RESTART_FORMAT_VERSION;
        nc_put_att_int(ncid, NC_GLOBAL, "restart_format_version", NC_INT, 1, &v);
        int step = 1; nc_put_att_int(ncid, NC_GLOBAL, "step", NC_INT, 1, &step);
        double dt = 1800.0; nc_put_att_double(ncid, NC_GLOBAL, "dt", NC_DOUBLE, 1, &dt);
        int cali[6] = { 0, 1963, 7, 19, 13, 30 };
        nc_put_att_int(ncid, NC_GLOBAL, "calendar", NC_INT, 6, cali);
        double sec = 0.0; nc_put_att_double(ncid, NC_GLOBAL, "second", NC_DOUBLE, 1, &sec);
        int org[3] = { 1958, 1, 1 };
        nc_put_att_int(ncid, NC_GLOBAL, "calendar_origin", NC_INT, 3, org);
        long long sig = 16;
        nc_put_att_longlong(ncid, NC_GLOBAL, "sum_nlevels_nod2D", NC_INT64, 1, &sig);
        nc_close(ncid);
    } else if (strcmp(which, "truncated") == 0) {
        MUST(truncate(path, 512) == 0);
    }
    /* "wrong_mesh" needs no doctoring of the file — the reader below uses a
     * mesh with a different bathymetry. */

    model dst;
    model_init(&dst, strcmp(which, "wrong_mesh") == 0 ? NLEV_B : NLEV_A, 0);
    int step = 0; double dt = 0.0;
    fesom_calendar_t cal;
    memset(&cal, 0, sizeof(cal));
    fesom_restart_read(path, &step, &dt, &cal, &dst.mesh, &dst.dyn,
                       &dst.tracers, &dst.ice, NULL, &dst.stiff, &dst.partit);

    fprintf(stderr, "guard case '%s' returned instead of aborting\n", which);
    exit(0);   /* the parent treats a zero exit as a failure */
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc > 1) {          /* child: one guard case, expected to abort */
        run_guard_case(argv[1]);
        MPI_Finalize();
        return 0;
    }

    char dir[] = "/tmp/fesom_rst_XXXXXX";
    MUST(mkdtemp(dir));

    test_round_trip(dir);
    test_components_are_not_swapped(dir);
    test_tke_is_optional(dir);
    test_config_and_cadence();

    static const char *guards[] = { "version", "missing_var", "truncated", "wrong_mesh" };
    for (size_t i = 0; i < sizeof guards / sizeof guards[0]; ++i) {
        char cmd[1024];
        snprintf(cmd, sizeof cmd, "%s %s >/dev/null 2>&1", argv[0], guards[i]);
        int rc = system(cmd);
        if (rc == 0) {
            ++n_failed;
            fprintf(stderr, "FAIL: guard case '%s' did not abort\n", guards[i]);
        } else {
            ++n_passed;
        }
    }

    printf("\n===== %d passed, %d failed =====\n", n_passed, n_failed);
    MPI_Finalize();
    return n_failed == 0 ? 0 : 1;
}
