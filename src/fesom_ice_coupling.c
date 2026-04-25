#include "fesom_ice_coupling.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_sss_runoff.h"
#include "fesom_tracers.h"
#include "fesom_types.h"

#include <math.h>
#include <mpi.h>

/*
 * integrate_nod_2D — gen_support.F90:318-351.
 * Local sum over interior owned nodes weighted by surface areasvol, then
 * global Allreduce. (Inlined here; identical to the static function in
 * fesom_sss_runoff.c — a small duplication kept to avoid widening any API
 * for a 10-line helper.)
 */
static real_t integrate_nod_2D(const real_t            *data,
                               const struct fesom_mesh *mesh,
                               struct fesom_partit     *partit)
{
    real_t lval = 0.0;
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        int ul = mesh->ulevels_nod2D[n] - 1;   /* 0-based */
        lval += data[n] * mesh->areasvol[FESOM_NODE3D(n, ul, mesh->nl)];
    }
    if (partit && partit->npes > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &lval, 1, MPI_DOUBLE, MPI_SUM,
                      partit->MPI_COMM_FESOM);
    }
    return lval;
}

/*
 * Mirror of Fortran ocean2ice (ice_oce_coupling.F90:154).
 * Cavity skip uses ulevels_nod2D > 1 (matches Fortran).
 *
 * Uses the area-weighted-mean recipe directly from Fortran (lines 211-242)
 * rather than reading dyn->uvnode (which would also be valid, since uvnode
 * is computed by the same recipe earlier in the step). Literal port.
 */
void fesom_ocean2ice(fesom_ice                  *ice,
                     const struct fesom_dyn     *dyn,
                     const struct fesom_tracers *tracers,
                     struct fesom_partit        *partit,
                     struct fesom_mesh          *mesh)
{
    int N      = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int nl     = mesh->nl;
    real_t *u_w = ice->srfoce_u;
    real_t *v_w = ice->srfoce_v;

    /* T_oc, S_oc, elevation = surface tracer + hbar (Fortran lines 191-198) */
    if (ice->ice_update) {
        for (int n = 0; n < N; ++n) {
            if (mesh->ulevels_nod2D[n] > 1) continue;
            ice->srfoce_temp[n] = tracers->data[FESOM_TRACER_T].values[n * nl + 0];
            ice->srfoce_salt[n] = tracers->data[FESOM_TRACER_S].values[n * nl + 0];
            ice->srfoce_ssh [n] = mesh->hbar[n];
        }
    } else {
        /* Time-averaging branch — inactive when ice_ave_steps=1; mirror anyway. */
        real_t sf = (real_t)ice->ice_steps_since_upd;
        real_t sn = sf + 1.0;
        for (int n = 0; n < N; ++n) {
            if (mesh->ulevels_nod2D[n] > 1) continue;
            ice->srfoce_temp[n] = (ice->srfoce_temp[n] * sf
                                   + tracers->data[FESOM_TRACER_T].values[n * nl + 0]) / sn;
            ice->srfoce_salt[n] = (ice->srfoce_salt[n] * sf
                                   + tracers->data[FESOM_TRACER_S].values[n * nl + 0]) / sn;
            ice->srfoce_ssh [n] = (ice->srfoce_ssh [n] * sf + mesh->hbar[n]) / sn;
        }
    }

    /* Zero u_w/v_w then accumulate area-weighted surface UV per node.
       Fortran lines 211-216 zero over myDim+eDim, lines 218-242 compute over myDim. */
    for (int n = 0; n < N; ++n) { u_w[n] = 0.0; v_w[n] = 0.0; }

    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) continue;
        real_t uw = 0.0, vw = 0.0, vol = 0.0;
        int beg = mesh->nod_in_elem2D_offsets[n];
        int end = mesh->nod_in_elem2D_offsets[n + 1];
        for (int k = beg; k < end; ++k) {
            int elem = mesh->nod_in_elem2D[k];
            if (mesh->ulevels[elem] > 1) continue;
            real_t a = mesh->elem_area[elem];
            vol += a;
            uw  += dyn->uv[FESOM_ELEMVEC(elem, 0, nl) + 0] * a;
            vw  += dyn->uv[FESOM_ELEMVEC(elem, 0, nl) + 1] * a;
        }
        if (vol > 0.0) {
            uw /= vol;
            vw /= vol;
        }
        if (ice->ice_update) {
            u_w[n] = uw;
            v_w[n] = vw;
        } else {
            real_t sf = (real_t)ice->ice_steps_since_upd;
            real_t sn = sf + 1.0;
            u_w[n] = (u_w[n] * sf + uw) / sn;
            v_w[n] = (v_w[n] * sf + vw) / sn;
        }
    }

    /* Halo exchange — Fortran line 244 */
    fesom_exchange_nod2D(u_w, partit);
    fesom_exchange_nod2D(v_w, partit);
}

/*
 * Mirror of Fortran oce_fluxes (ice_oce_coupling.F90:249), non-ICEPACK
 * branch (lines 393-398). Also includes the virtual_salt block (436-457)
 * and relax_salt block (498-528) which used to live in fesom_sss_runoff.c
 * but properly belong here per the Fortran source.
 *
 * Cavity, wiso, iceberg branches all skipped (out of scope).
 */
void fesom_ice_oce_fluxes(fesom_ice                     *ice,
                          struct fesom_partit           *partit,
                          struct fesom_mesh             *mesh,
                          const struct fesom_tracers    *tracers,
                          struct fesom_forcing          *forcing,
                          const struct fesom_sss_runoff *sr)
{
    int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    int nl = mesh->nl;

    /* Non-ICEPACK branch — Fortran lines 393-398.
     * NOTE: NO `-runoff` term. Runoff is already inside ice->flx_fw via therm_ice. */
    for (int n = 0; n < N; ++n) {
        forcing->heat_flux [n] = -ice->flx_h [n];
        forcing->water_flux[n] = -ice->flx_fw[n];
    }

    /* heat_flux_in (used later by sw_pene if shortwave separation lands).
     * Not stored in our forcing struct yet; deferred. Fortran lines 416-420. */

    /* Virtual salt flux — Fortran lines 436-457. Linfs uses_virt_salt=true. */
    if (sr->use_virt_salt) {
        real_t rsss_default = sr->ref_sss;
        for (int n = 0; n < N; ++n) {
            real_t rsss = sr->ref_sss_local
                            ? tracers->data[FESOM_TRACER_S].values[n * nl + 0]
                            : rsss_default;
            forcing->virtual_salt[n] = rsss * forcing->water_flux[n];
        }
        real_t net = integrate_nod_2D(forcing->virtual_salt, mesh, partit)
                     / mesh->ocean_area;
        for (int n = 0; n < mesh->myDim_nod2D; ++n) {
            if (mesh->ulevels_nod2D[n] > 1) continue;  /* leave cavity nodes alone */
            forcing->virtual_salt[n] -= net;
        }
        fesom_exchange_nod2D(forcing->virtual_salt, partit);
    }

    /* SSS restoring — Fortran lines 498-528. */
    for (int n = 0; n < N; ++n) {
        forcing->relax_salt[n] = sr->surf_relax_S
                                 * (forcing->Ssurf[n]
                                    - tracers->data[FESOM_TRACER_S].values[n * nl + 0]);
    }
    real_t net_relax = integrate_nod_2D(forcing->relax_salt, mesh, partit)
                       / mesh->ocean_area;
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        forcing->relax_salt[n] -= net_relax;
    }
    fesom_exchange_nod2D(forcing->relax_salt, partit);

    /* Halo exchange of heat_flux/water_flux — used by tracer diffusion BC. */
    fesom_exchange_nod2D(forcing->heat_flux,  partit);
    fesom_exchange_nod2D(forcing->water_flux, partit);
}

