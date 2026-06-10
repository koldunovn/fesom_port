#include "fesom_ale_dump.h"
#include "fesom_aux.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_types.h"

#include <stdio.h>
#include <stdlib.h>

static const char *s_dir    = NULL;
static int         s_steps  = 3;
static int         s_loaded = 0;

static void ale_dump_load_env(void)
{
    if (s_loaded) return;
    s_dir = getenv("FESOM_ALE_DUMP_DIR");
    if (s_dir && !s_dir[0]) s_dir = NULL;
    const char *e = getenv("FESOM_ALE_DUMP_STEPS");
    if (e && e[0]) {
        int v = atoi(e);
        if (v > 0) s_steps = v;
    }
    s_loaded = 1;
}

int fesom_ale_dump_active(int step)
{
    ale_dump_load_env();
    return s_dir != NULL && step <= s_steps;
}

/* getter signature: value of component c (0-based) at local row i */
typedef double (*ale_get_fn)(int i, int c, const void *user);

static void ale_dump_write(int step, const char *tag, int ncomp, int nrows,
                           const int *gids, ale_get_fn get, const void *user,
                           int mype)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/ale_dump_s%d_%s_rank%d.txt",
             s_dir, step, tag, mype);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "ale_dump: cannot open %s\n", path); return; }
    fprintf(f, "# step=%d tag=%s rank=%d N=%d ncomp=%d\n",
            step, tag, mype, nrows, ncomp);
    for (int i = 0; i < nrows; ++i) {
        fprintf(f, "%d", gids[i]);                       /* 1-based gid */
        for (int c = 0; c < ncomp; ++c) fprintf(f, " %.17g", get(i, c, user));
        fputc('\n', f);
    }
    fclose(f);
}

/* --- getters -------------------------------------------------------------- */

/* up to 4 independent 2D arrays as components */
typedef struct { const real_t *a[4]; } ale_src_2d;
static double get_2d(int i, int c, const void *user)
{
    const ale_src_2d *s = (const ale_src_2d *)user;
    return (double)s->a[c][i];
}

/* one [N*nl] node/elem-major 3D array; comp = level */
typedef struct { const real_t *arr; int nl; } ale_src_col;
static double get_col(int i, int c, const void *user)
{
    const ale_src_col *s = (const ale_src_col *)user;
    return (double)s->arr[(size_t)i * (size_t)s->nl + (size_t)c];
}

/* --- bookends (mirror the Fortran ale_dump_* subroutines) ----------------- */

void fesom_ale_dump_forcing(int step, const struct fesom_forcing *forcing,
                            const struct fesom_mesh *mesh,
                            struct fesom_partit *partit)
{
    if (!fesom_ale_dump_active(step)) return;
    ale_src_2d s = { { forcing->water_flux, forcing->virtual_salt,
                       forcing->relax_salt, forcing->real_salt_flux } };
    ale_dump_write(step, "forcing", 4, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_2d, &s, partit->mype);
}

void fesom_ale_dump_pgf(int step, const struct fesom_aux *aux,
                        const struct fesom_mesh *mesh,
                        struct fesom_partit *partit)
{
    if (!fesom_ale_dump_active(step)) return;
    ale_src_col sx = { aux->pgf_x, mesh->nl };
    ale_src_col sy = { aux->pgf_y, mesh->nl };
    ale_dump_write(step, "pgf_x", mesh->nl - 1, mesh->myDim_elem2D,
                   partit->myList_elem2D, get_col, &sx, partit->mype);
    ale_dump_write(step, "pgf_y", mesh->nl - 1, mesh->myDim_elem2D,
                   partit->myList_elem2D, get_col, &sy, partit->mype);
}

void fesom_ale_dump_sshsolve(int step, const struct fesom_dyn *dyn,
                             const struct fesom_mesh *mesh,
                             struct fesom_partit *partit)
{
    if (!fesom_ale_dump_active(step)) return;
    ale_src_2d s = { { dyn->ssh_rhs, dyn->d_eta, NULL, NULL } };
    ale_dump_write(step, "sshsolve", 2, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_2d, &s, partit->mype);
}

void fesom_ale_dump_hbar(int step, const struct fesom_dyn *dyn,
                         const struct fesom_mesh *mesh,
                         struct fesom_partit *partit)
{
    if (!fesom_ale_dump_active(step)) return;
    ale_src_2d s = { { mesh->hbar, mesh->hbar_old, dyn->ssh_rhs_old,
                       dyn->eta_n } };
    ale_dump_write(step, "hbar", 4, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_2d, &s, partit->mype);
    ale_src_2d sd = { { mesh->dhe, NULL, NULL, NULL } };
    ale_dump_write(step, "dhe", 1, mesh->myDim_elem2D,
                   partit->myList_elem2D, get_2d, &sd, partit->mype);
}

void fesom_ale_dump_vertvel(int step, const struct fesom_dyn *dyn,
                            const struct fesom_mesh *mesh,
                            struct fesom_partit *partit)
{
    if (!fesom_ale_dump_active(step)) return;
    ale_src_col sw = { dyn->w,          mesh->nl };
    ale_src_col sh = { mesh->hnode_new, mesh->nl };
    ale_dump_write(step, "Wvel", mesh->nl, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_col, &sw, partit->mype);
    ale_dump_write(step, "hnode_new", mesh->nl - 1, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_col, &sh, partit->mype);
}

void fesom_ale_dump_thickness(int step, const struct fesom_mesh *mesh,
                              struct fesom_partit *partit)
{
    if (!fesom_ale_dump_active(step)) return;
    ale_src_col sh = { mesh->hnode,     mesh->nl };
    ale_src_col sz = { mesh->zbar_3d_n, mesh->nl };
    ale_src_col sZ = { mesh->Z_3d_n,    mesh->nl };
    ale_src_col se = { mesh->helem,     mesh->nl };
    ale_dump_write(step, "hnode", mesh->nl - 1, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_col, &sh, partit->mype);
    ale_dump_write(step, "zbar_3d_n", mesh->nl, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_col, &sz, partit->mype);
    ale_dump_write(step, "Z_3d_n", mesh->nl - 1, mesh->myDim_nod2D,
                   partit->myList_nod2D, get_col, &sZ, partit->mype);
    ale_dump_write(step, "helem", mesh->nl - 1, mesh->myDim_elem2D,
                   partit->myList_elem2D, get_col, &se, partit->mype);
}
