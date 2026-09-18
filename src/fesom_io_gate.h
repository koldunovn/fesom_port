#ifndef FESOM_IO_GATE_H
#define FESOM_IO_GATE_H
#include "fesom_types.h"
struct fesom_mesh; struct fesom_dyn; struct fesom_tracers; struct fesom_aux;
struct fesom_ice; struct fesom_forcing; struct fesom_jra55; struct fesom_partit;
/* NULL unless FESOM_GATE_DUMP_DIR is set. */
const char *fesom_gate_dir(void);
void fesom_io_write_gate(int step_n,
                         const struct fesom_mesh     *mesh,
                         const struct fesom_dyn      *dyn,
                         const struct fesom_tracers  *tracers,
                         const struct fesom_aux      *aux,
                         const struct fesom_ice      *ice,
                         const struct fesom_forcing  *forcing,
                         const struct fesom_jra55    *jra,
                         struct fesom_partit         *partit);
void fesom_io_write_gate_2d(int step_n, int sub, int nv, const char **names, const real_t **ptrs,
                            const struct fesom_mesh *mesh, struct fesom_partit *partit);
int fesom_gate_inject_node(const char *path, const char *var, real_t *local,
                           const struct fesom_mesh *mesh, struct fesom_partit *partit);
#endif
