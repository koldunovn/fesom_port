/*
 * IO subsystem - DOCUMENTED EXCEPTION to the literal Fortran->C port rule.
 * See fesom_io_stream.h banner for context.
 *
 * fesom_io_stream_dispatch.c - small set of pure dispatcher functions
 * carved out of fesom_io_stream.c so test_io_stream_unit can link them
 * without dragging in netCDF + the full IO TU. Keeps the stride
 * invariant for FESOM_VAR_3D_NODE_MID under unit-test coverage.
 */
#include "fesom_io_stream.h"
#include "fesom_types.h"

size_t fesom_io_stream_accum_size(fesom_var_kind_t kind, int local_n, int local_e, int nl)
{
    /* Stride for 3D node fields is `nl`, NOT `nl-1` — see
     * feedback_tracer_stride_nl.md (cost a session of debug time). */
    switch (kind) {
        case FESOM_VAR_2D_NODE:        return (size_t)local_n;
        case FESOM_VAR_3D_NODE_MID:    return (size_t)local_n * (size_t)nl;
        case FESOM_VAR_3D_NODE_IFACE:  return (size_t)local_n * (size_t)nl;
        case FESOM_VAR_2D_ELEM:        return (size_t)local_e;
        case FESOM_VAR_3D_ELEM_MID:    return (size_t)local_e * (size_t)nl;
        case FESOM_VAR_3D_ELEM_IFACE:  return (size_t)local_e * (size_t)nl;
        case FESOM_VAR_2D_ELEM_VEC:    return (size_t)local_e * 2;
    }
    FESOM_DIE("io_stream: unknown var kind %d", (int)kind);
    return 0;
}

size_t fesom_io_stream_global_size(fesom_var_kind_t kind, int nod2D, int elem2D, int nl)
{
    switch (kind) {
        case FESOM_VAR_2D_NODE:        return (size_t)nod2D;
        case FESOM_VAR_3D_NODE_MID:    return (size_t)nod2D * (size_t)nl;
        case FESOM_VAR_3D_NODE_IFACE:  return (size_t)nod2D * (size_t)nl;
        case FESOM_VAR_2D_ELEM:        return (size_t)elem2D;
        case FESOM_VAR_3D_ELEM_MID:    return (size_t)elem2D * (size_t)nl;
        case FESOM_VAR_3D_ELEM_IFACE:  return (size_t)elem2D * (size_t)nl;
        case FESOM_VAR_2D_ELEM_VEC:    return (size_t)elem2D * 2;
    }
    FESOM_DIE("io_stream: unknown var kind %d", (int)kind);
    return 0;
}

int fesom_io_stream_n_layers_written(fesom_var_kind_t kind, int nl)
{
    switch (kind) {
        case FESOM_VAR_2D_NODE:        return 0;
        case FESOM_VAR_3D_NODE_MID:    return nl - 1;
        case FESOM_VAR_3D_NODE_IFACE:  return nl;
        case FESOM_VAR_2D_ELEM:        return 0;
        case FESOM_VAR_3D_ELEM_MID:    return nl - 1;
        case FESOM_VAR_3D_ELEM_IFACE:  return nl;
        case FESOM_VAR_2D_ELEM_VEC:    return 0;
    }
    return 0;
}
