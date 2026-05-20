/*
 * Unit test for fesom_io_stream dispatcher functions. Catches stride
 * bugs at unit-test time — particularly the FESOM_VAR_3D_NODE_MID
 * accumulator-size invariant of `local_n * nl` (NOT `local_n * (nl-1)`),
 * see feedback_tracer_stride_nl.md.
 *
 * Builds against fesom_io_stream_dispatch.c only (no netCDF, no fesom_io.c).
 * MPI is linked because FESOM_DIE references MPI_Abort, but never invoked
 * in the passing test path.
 */
#include "fesom_io_stream.h"

#include <stdio.h>
#include <stdlib.h>

static int n_passed = 0, n_failed = 0;

#define CHECK(expr)                                                         \
    do {                                                                    \
        if (expr) ++n_passed;                                               \
        else      {                                                         \
            ++n_failed;                                                     \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr); \
        }                                                                   \
    } while (0)

#define CHECK_EQ(actual, expected)                                          \
    do {                                                                    \
        size_t _a = (size_t)(actual);                                       \
        size_t _e = (size_t)(expected);                                     \
        if (_a == _e) ++n_passed;                                           \
        else          {                                                     \
            ++n_failed;                                                     \
            fprintf(stderr, "FAIL %s:%d: %s = %zu, expected %zu\n",         \
                    __FILE__, __LINE__, #actual, _a, _e);                   \
        }                                                                   \
    } while (0)

int main(void)
{
    /* Standard small fixture: 4 owned nodes, 5 owned elements, 6 layers. */
    const int N  = 4;
    const int E  = 5;
    const int nl = 6;

    /* ---- accum_size: per-rank local accumulator buffer size ------- */
    /* The CRITICAL invariant: 3D_NODE_MID buffer is sized N*nl, not
     * N*(nl-1). Tracers store columns at stride nl; resolver returns
     * `tracers->data[..].values + node*nl`, and the accumulator must
     * have the same stride or sums shift across layer boundaries. */
    CHECK_EQ(fesom_io_stream_accum_size(FESOM_VAR_2D_NODE,       N, E, nl), N);
    CHECK_EQ(fesom_io_stream_accum_size(FESOM_VAR_3D_NODE_MID,   N, E, nl), N * nl);
    CHECK_EQ(fesom_io_stream_accum_size(FESOM_VAR_3D_NODE_IFACE, N, E, nl), N * nl);
    CHECK_EQ(fesom_io_stream_accum_size(FESOM_VAR_2D_ELEM,       N, E, nl), E);
    CHECK_EQ(fesom_io_stream_accum_size(FESOM_VAR_3D_ELEM_MID,   N, E, nl), E * nl);
    CHECK_EQ(fesom_io_stream_accum_size(FESOM_VAR_2D_ELEM_VEC,   N, E, nl), E * 2);

    /* Stride bug regression: accidentally using nl-1 yields E*(nl-1) = 25
     * instead of E*nl = 30. Make sure these are not equal. */
    CHECK(fesom_io_stream_accum_size(FESOM_VAR_3D_NODE_MID, N, E, nl)
          != (size_t)(N * (nl - 1)));
    CHECK(fesom_io_stream_accum_size(FESOM_VAR_3D_ELEM_MID, N, E, nl)
          != (size_t)(E * (nl - 1)));

    /* ---- global_size: rank-0 global buffer size for one record ---- */
    const int gN = 100, gE = 220;
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_2D_NODE,       gN, gE, nl), gN);
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_3D_NODE_MID,   gN, gE, nl), gN * nl);
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_3D_NODE_IFACE, gN, gE, nl), gN * nl);
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_2D_ELEM,       gN, gE, nl), gE);
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_3D_ELEM_MID,   gN, gE, nl), gE * nl);
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_2D_ELEM_VEC,   gN, gE, nl), gE * 2);

    /* ---- n_layers_written: how many layers go into the netCDF file */
    /* This differs from accum stride: file only writes (nl-1) for mid,
     * nl for interface, 0 for 2D. The mismatch with accum_size on the
     * mid kinds is intentional (write loops index k=0..(nl-2) with the
     * accum stride still being nl). */
    CHECK_EQ(fesom_io_stream_n_layers_written(FESOM_VAR_2D_NODE,       nl), 0);
    CHECK_EQ(fesom_io_stream_n_layers_written(FESOM_VAR_3D_NODE_MID,   nl), nl - 1);
    CHECK_EQ(fesom_io_stream_n_layers_written(FESOM_VAR_3D_NODE_IFACE, nl), nl);
    CHECK_EQ(fesom_io_stream_n_layers_written(FESOM_VAR_2D_ELEM,       nl), 0);
    CHECK_EQ(fesom_io_stream_n_layers_written(FESOM_VAR_3D_ELEM_MID,   nl), nl - 1);
    CHECK_EQ(fesom_io_stream_n_layers_written(FESOM_VAR_2D_ELEM_VEC,   nl), 0);

    /* ---- Realistic CORE2 sizes (smoke check the arithmetic) ------- */
    /* CORE2: nod2D=126858, elem2D=244659, nl=48.                     */
    const int CORE2_N = 126858, CORE2_E = 244659, CORE2_NL = 48;
    /* T per record: 126858 * 48 = 6_089_184 elements (≈ 48 MB at f64). */
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_3D_NODE_MID,
             CORE2_N, CORE2_E, CORE2_NL), 6089184);
    /* uv per record: 244659 * 48 = 11_743_632 elements per component,
     *  × 2 components → 23_487_264 elements total. */
    CHECK_EQ(fesom_io_stream_global_size(FESOM_VAR_3D_ELEM_MID,
             CORE2_N, CORE2_E, CORE2_NL), 11743632);

    fprintf(stderr, "test_io_stream_unit: %d passed, %d failed\n",
            n_passed, n_failed);
    return n_failed == 0 ? 0 : 1;
}
