#include "fesom_ice_coupling.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_forcing.h"
#include "fesom_halo.h"
#include "fesom_kpp.h"     /* [forcing-diff dig] reuse the KPP dump harness for the iceforce dump */
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

/* [forcing-diff dig] iceforce dump getter: 10 comps/node — a_ice, m_ice, u_ice,
 * v_ice, srfoce_u, srfoce_v, stress_iceoce_x/y, stress_node_surf_x/y. */
typedef struct {
    const real_t *aice, *mice, *uice, *vice, *sru, *srv, *siox, *sioy, *sns;
} fesom_icf_src;
static double fesom_get_icf(int n, int c, void *u)
{
    const fesom_icf_src *s = (const fesom_icf_src *)u;
    switch (c) {
        case 0:  return (double)s->aice[n];      case 1:  return (double)s->mice[n];
        case 2:  return (double)s->uice[n];      case 3:  return (double)s->vice[n];
        case 4:  return (double)s->sru[n];       case 5:  return (double)s->srv[n];
        case 6:  return (double)s->siox[n];      case 7:  return (double)s->sioy[n];
        case 8:  return (double)s->sns[2*n + 0]; default: return (double)s->sns[2*n + 1];
    }
}

/*
 * Mirror of Fortran oce_fluxes_mom (ice_oce_coupling.F90:53-149).
 *
 * Two loops, both literal:
 *   per-node loop: 0..myDim+eDim — compute stress_iceoce_x/y, blend with
 *     atmoce (snapshot from forcing->stress_node_surf as written by
 *     fesom_bulk_compute), write back the blended stress_node_surf.
 *   per-element loop: 0..myDim_elem2D — stress_surf = mean of the 3 vertex
 *     stress_node_surf. Matches Fortran's myDim-only bound; impl_vert_visc
 *     also reads only myDim.
 *
 * No halo exchanges are emitted: the per-node loop covers myDim+eDim, so
 * the per-element loop reads stress_node_surf at owned-element vertices
 * (myDim or eDim) against fresh values.
 */
void fesom_ice_oce_fluxes_mom(fesom_ice              *ice,
                              struct fesom_partit    *partit,
                              struct fesom_mesh      *mesh,
                              struct fesom_forcing   *forcing)
{
    (void)partit;  /* No halo exchange needed; see header. */

    const int N      = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int Eown   = mesh->myDim_elem2D;
    const real_t rho_cd = (real_t)FESOM_DENSITY_0 * ice->cd_oce_ice;
    real_t *u_ice = ice->uice;
    real_t *v_ice = ice->vice;
    real_t *u_w   = ice->srfoce_u;
    real_t *v_w   = ice->srfoce_v;
    real_t *a_ice = ice->data[0].values;            /* ice concentration */
    real_t *sicx  = ice->stress_iceoce_x;
    real_t *sicy  = ice->stress_iceoce_y;
    real_t *sns   = forcing->stress_node_surf;      /* [2*n+0,1] */
    real_t *ss    = forcing->stress_surf;           /* [2*el+0,1] */

    /* Per-node — Fortran lines 105-122. */
    for (int n = 0; n < N; ++n) {
        if (mesh->ulevels_nod2D[n] > 1) continue;     /* cavity */

        real_t a = a_ice[n];
        real_t atmx = sns[2*n + 0];                   /* atmoce snapshot */
        real_t atmy = sns[2*n + 1];

        if (a > 0.001) {
            real_t du  = u_ice[n] - u_w[n];
            real_t dv  = v_ice[n] - v_w[n];
            real_t aux = sqrt(du*du + dv*dv) * rho_cd;
            sicx[n] = aux * du;
            sicy[n] = aux * dv;
        } else {
            sicx[n] = 0.0;
            sicy[n] = 0.0;
        }
        sns[2*n + 0] = sicx[n] * a + atmx * (1.0 - a);
        sns[2*n + 1] = sicy[n] * a + atmy * (1.0 - a);
    }

    /* Per-element — Fortran lines 127-143. Bound is myDim_elem2D only. */
    for (int e = 0; e < Eown; ++e) {
        if (mesh->ulevels[e] > 1) continue;          /* cavity */

        int v0 = mesh->elem_nodes[3*e + 0];
        int v1 = mesh->elem_nodes[3*e + 1];
        int v2 = mesh->elem_nodes[3*e + 2];
        ss[2*e + 0] = (sns[2*v0 + 0] + sns[2*v1 + 0] + sns[2*v2 + 0]) / 3.0;
        ss[2*e + 1] = (sns[2*v0 + 1] + sns[2*v1 + 1] + sns[2*v2 + 1]) / 3.0;
    }

    /* [forcing-diff dig] step-1 ice/stress dump (FESOM_KPP_DUMP_DIR-gated, tag
     * 'iceforce') — pinpoint the sea-ice source of the C-vs-Fortran step-1 ustar
     * diff. Mirrors the Fortran oce_fluxes_mom dump; called once per step. */
    static int s_icf_call = 0;
    ++s_icf_call;
    if (fesom_kpp_dump_this_step(s_icf_call)) {
        fesom_icf_src src = { a_ice, ice->data[1].values, u_ice, v_ice,
                              u_w, v_w, sicx, sicy, sns };
        fesom_kpp_dump_nodes(mesh, partit, s_icf_call, "iceforce", 10, fesom_get_icf, &src);
    }
}
