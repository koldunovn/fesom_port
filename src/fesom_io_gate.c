/*
 * fesom_io_gate.c — per-step full-field dump for the matched-state C-vs-Fortran
 * gate (drift investigation 2026-09-18, docs/plans/20260918-longrun-drift-strategy.md
 * Track B). Env FESOM_GATE_DUMP_DIR=<dir>: after every step, write
 * <dir>/gate_NNNNNN.nc with the surface fluxes, the forcing arrays, the ice
 * state and the 3-D T/S/w/Kv and u/v, all gathered to global numbering (double).
 * Table-driven; nothing here touches model state. Rank 0 writes.
 */
#include "fesom_io_gate.h"
#include "fesom_aux.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_ice_types.h"
#include "fesom_io_gather.h"
#include "fesom_jra55.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NCG(call) do { int _e = (call); if (_e != NC_NOERR) \
    FESOM_DIE("gate dump: NetCDF error: %s", nc_strerror(_e)); } while (0)

typedef struct { const char *name; const real_t *p; int is_elem; int levels; int ncomp; int comp; } gate_var;

const char *fesom_gate_dir(void)
{
    static int checked = 0; static const char *dir = NULL;
    if (!checked) { dir = getenv("FESOM_GATE_DUMP_DIR"); if (dir && !dir[0]) dir = NULL; checked = 1; }
    return dir;
}

void fesom_io_write_gate(int step_n,
                         const struct fesom_mesh     *mesh,
                         const struct fesom_dyn      *dyn,
                         const struct fesom_tracers  *tracers,
                         const struct fesom_aux      *aux,
                         const struct fesom_ice      *ice,
                         const struct fesom_forcing  *forcing,
                         const struct fesom_jra55    *jra,
                         struct fesom_partit         *partit)
{
    const char *dir = fesom_gate_dir();
    if (!dir) return;
    const int nl = mesh->nl;
    gate_var v[64]; int nv = 0;
#define N2(nm, ptr)        v[nv++] = (gate_var){ nm, ptr, 0, 1,  1, 0 }
#define E2(nm, ptr)        v[nv++] = (gate_var){ nm, ptr, 1, 1,  1, 0 }
#define N3(nm, ptr)        v[nv++] = (gate_var){ nm, ptr, 0, nl, 1, 0 }
#define E3(nm, ptr)        v[nv++] = (gate_var){ nm, ptr, 1, nl, 1, 0 }
#define E3V(nm, ptr, c)    v[nv++] = (gate_var){ nm, ptr, 1, nl, 2, c }
#define N2V(nm, ptr, c)    v[nv++] = (gate_var){ nm, ptr, 0, 1,  2, c }
    if (forcing) {
        N2("heat_flux",       forcing->heat_flux);
        N2("water_flux",      forcing->water_flux);
        N2("real_salt_flux",  forcing->real_salt_flux);
        N2("relax_salt",      forcing->relax_salt);
        N2("virtual_salt",    forcing->virtual_salt);
        N2("evaporation",     forcing->evaporation);
        N2("ice_sublimation", forcing->ice_sublimation);
        N2("runoff",          forcing->runoff);
        N2("Ssurf",           forcing->Ssurf);
        N2("Ch_atm_oce",      forcing->Ch_atm_oce);
        N2("Ce_atm_oce",      forcing->Ce_atm_oce);
        N2V("stress_node_x",  forcing->stress_node_surf, 0);
        N2V("stress_node_y",  forcing->stress_node_surf, 1);
        v[nv++] = (gate_var){ "tx_sur", forcing->stress_surf, 1, 1, 2, 0 };
        v[nv++] = (gate_var){ "ty_sur", forcing->stress_surf, 1, 1, 2, 1 };
        N3("sw_3d",           forcing->sw_3d);
    }
    if (jra) {
        N2("u_wind", jra->u_wind); N2("v_wind", jra->v_wind); N2("Tair", jra->Tair);
        N2("shum", jra->shum); N2("shortwave", jra->shortwave); N2("longwave", jra->longwave);
        N2("prec_rain", jra->prec_rain); N2("prec_snow", jra->prec_snow);
    }
    if (ice) {
        N2("a_ice", ice->data[FESOM_ICE_AICE].values);
        N2("m_ice", ice->data[FESOM_ICE_MICE].values);
        N2("m_snow", ice->data[FESOM_ICE_MSNOW].values);
        N2("uice", ice->uice); N2("vice", ice->vice);
        N2("t_skin", ice->thermo.t_skin); N2("ustar_io", ice->thermo.ustar);
        N2("thdgr", ice->thermo.thdgr); N2("thdgrsn", ice->thermo.thdgrsn); N2("thdgra", ice->thermo.thdgra);
        N2("flx_h", ice->flx_h); N2("flx_fw", ice->flx_fw);
        N2("srfoce_u", ice->srfoce_u); N2("srfoce_v", ice->srfoce_v);
        N2("srfoce_temp", ice->srfoce_temp); N2("srfoce_salt", ice->srfoce_salt);
        N2("stress_atmice_x", ice->stress_atmice_x); N2("stress_atmice_y", ice->stress_atmice_y);
        N2("stress_iceoce_x", ice->stress_iceoce_x); N2("stress_iceoce_y", ice->stress_iceoce_y);
    }
    N2("eta_n", dyn->eta_n); N2("hbar", mesh->hbar); N2("ssh_rhs_old", dyn->ssh_rhs_old);
    N3("temp", tracers->data[FESOM_TRACER_T].values);
    N3("salt", tracers->data[FESOM_TRACER_S].values);
    N3("w", dyn->w); N3("Kv", aux->Kv); N3("bvfreq", aux->bvfreq); N3("hnode", mesh->hnode);
    N3("uvnode_u", dyn->uvnode); v[nv-1].ncomp = 2; v[nv-1].comp = 0;
    N3("uvnode_v", dyn->uvnode); v[nv-1].ncomp = 2; v[nv-1].comp = 1;
    E3V("u", dyn->uv, 0); E3V("v", dyn->uv, 1);
    E3("Av", aux->Av);

    gather_plan gp; gather_plan_init(&gp, mesh, partit);
    MPI_Comm comm = partit->MPI_COMM_FESOM;
    int ncid = -1, d_node, d_elem, d_nz, vid[64];
    char path[1024];
    if (partit->mype == 0) {
        snprintf(path, sizeof path, "%s/gate_%06d.nc", dir, step_n);
        NCG(nc_create(path, NC_CLOBBER | NC_NETCDF4, &ncid));
        NCG(nc_def_dim(ncid, "node", (size_t)mesh->nod2D, &d_node));
        NCG(nc_def_dim(ncid, "elem", (size_t)mesh->elem2D, &d_elem));
        NCG(nc_def_dim(ncid, "nz",   (size_t)nl, &d_nz));
        NCG(nc_put_att_int(ncid, NC_GLOBAL, "step", NC_INT, 1, &step_n));
        for (int i = 0; i < nv; ++i) {
            int dims[2] = { d_nz, v[i].is_elem ? d_elem : d_node };
            NCG(nc_def_var(ncid, v[i].name, NC_DOUBLE, v[i].levels > 1 ? 2 : 1,
                           v[i].levels > 1 ? dims : dims + 1, &vid[i]));
        }
        NCG(nc_enddef(ncid));
    }
    for (int i = 0; i < nv; ++i) {
        const size_t cnt = v[i].is_elem ? (size_t)mesh->elem2D : (size_t)mesh->nod2D;
        const int stride = v[i].levels * v[i].ncomp;
        real_t *g = NULL, *plane = NULL;
        if (partit->mype == 0) {
            g = malloc(cnt * (size_t)stride * sizeof(real_t));
            plane = malloc(cnt * (size_t)v[i].levels * sizeof(real_t));
            FESOM_CHECK(g && plane, "gate dump: alloc");
        }
        if (v[i].is_elem) gather_elem(v[i].p, stride, &gp, g, comm);
        else              gather_node(v[i].p, stride, &gp, g, comm);
        if (partit->mype == 0) {
            for (size_t e = 0; e < cnt; ++e)
                for (int lev = 0; lev < v[i].levels; ++lev)
                    plane[(size_t)lev * cnt + e] = g[e * (size_t)stride + (size_t)lev * v[i].ncomp + v[i].comp];
            NCG(nc_put_var_double(ncid, vid[i], plane));
            free(g); free(plane);
        }
    }
    if (partit->mype == 0) NCG(nc_close(ncid));
    gather_plan_free(&gp);
    MPI_Barrier(comm);
}

/* Per-substep 2-D node fields: gate_<step>_s<sub>.nc. Same env gate. */
void fesom_io_write_gate_2d(int step_n, int sub, int nv, const char **names, const real_t **ptrs,
                            const struct fesom_mesh *mesh, struct fesom_partit *partit)
{
    const char *dir = fesom_gate_dir();
    if (!dir) return;
    gather_plan gp; gather_plan_init(&gp, mesh, partit);
    MPI_Comm comm = partit->MPI_COMM_FESOM;
    int ncid = -1, d_node, vid[32];
    if (partit->mype == 0) {
        char path[1024];
        snprintf(path, sizeof path, "%s/gate_%06d_s%02d.nc", dir, step_n, sub);
        NCG(nc_create(path, NC_CLOBBER | NC_NETCDF4, &ncid));
        NCG(nc_def_dim(ncid, "node", (size_t)mesh->nod2D, &d_node));
        for (int i = 0; i < nv; ++i) NCG(nc_def_var(ncid, names[i], NC_DOUBLE, 1, &d_node, &vid[i]));
        NCG(nc_enddef(ncid));
    }
    for (int i = 0; i < nv; ++i) {
        real_t *g = NULL;
        if (partit->mype == 0) { g = malloc((size_t)mesh->nod2D * sizeof(real_t)); FESOM_CHECK(g, "gate2d alloc"); }
        gather_node(ptrs[i], 1, &gp, g, comm);
        if (partit->mype == 0) { NCG(nc_put_var_double(ncid, vid[i], g)); free(g); }
    }
    if (partit->mype == 0) NCG(nc_close(ncid));
    gather_plan_free(&gp);
}

/* Controlled replay: overwrite a node field from a global netCDF variable
 * (rank 0 reads, scatter, halo exchange). Returns 1 if done. */
#include "fesom_halo.h"
int fesom_gate_inject_node(const char *path, const char *var, real_t *local,
                           const struct fesom_mesh *mesh, struct fesom_partit *partit)
{
    gather_plan gp; gather_plan_init(&gp, mesh, partit);
    real_t *g = NULL;
    if (partit->mype == 0) {
        int ncid, vid;
        NCG(nc_open(path, NC_NOWRITE, &ncid));
        NCG(nc_inq_varid(ncid, var, &vid));
        g = malloc((size_t)mesh->nod2D * sizeof(real_t)); FESOM_CHECK(g, "inject alloc");
        NCG(nc_get_var_double(ncid, vid, g));
        NCG(nc_close(ncid));
    }
    scatter_node(g, 1, &gp, local, partit->MPI_COMM_FESOM);
    if (partit->mype == 0) free(g);
    fesom_halo_exchange(local, FESOM_HALO_NOD2D, 1, 1, partit);
    gather_plan_free(&gp);
    if (partit->mype == 0) printf("[gate] injected %s from %s\n", var, path);
    return 1;
}
