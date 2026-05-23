/*
 * Implicit vertical diffusion of tracers + surface heat/water flux BC.
 * Literal port of diff_ver_part_impl_ale (oce_ale_tracer.F90:562-1082)
 * and bc_surface (oce_ale_tracer.F90:1475-1524).
 *
 * Indexing convention (matches existing port, e.g. fesom_momentum.c):
 *   Fortran 1-based index k  ↔  C 0-based index k-1.
 *   nzmin / nzmax are converted to 0-based at the top of each node loop,
 *   so the loop bounds and array accesses can mirror the Fortran code with
 *   the same `nz` variable name without any -1 sprinkled at use sites.
 *
 * Fields not yet ported are taken as literal zero with a comment naming the
 * Fortran symbol:
 *   virtual_salt(n)     → 0   (zlevel/zstar only; we are linfs)
 *   relax_salt(n)       → 0   (Phase 3 step 25 will fill from SSS restoring)
 *   real_salt_flux(n)   → 0   (sea-ice; not ported)
 *   is_nonlinfs         → 0   (we are linfs)
 *   do_wimpl            → 0   (Fortran sets false when tra_adv_lim=='FCT')
 *   KPP nonlocal flux   → skipped (use_kpp_nonlclflx=false)
 *
 * Phase G5: Ki(nz,n) and slope_tapered(...) are now read when `gm` is
 * non-NULL. The Ty / Ty1 projection of horizontal-Redi onto the vertical
 * axis augments the diagonal diffusivity. With gm=NULL the additions
 * collapse to zero — preserving the pre-G5 numerical path.
 */
#include "fesom_tracer_diff.h"

#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_forcing.h"
#include "fesom_gm.h"
#include "fesom_mesh.h"
#include "fesom_tracers.h"
#include "fesom_types.h"

enum { NL_MAX = 64 };

/* ------------------------------------------------------------------------ *
 * bc_surface — oce_ale_tracer.F90:1475-1524                                *
 * Returns the surface forcing increment added to tr(nzmin).                *
 * id == 1 → temperature (Fortran tracers%data%ID==1)                       *
 * id == 2 → salinity                                                       *
 * ------------------------------------------------------------------------ */
static real_t bc_surface(int n,
                         int id,
                         real_t sval,
                         const struct fesom_forcing *forcing)
{
    const real_t dt          = (real_t)FESOM_PHASE1_DT;
    const real_t vcpw        = (real_t)FESOM_VCPW;
    const real_t is_nonlinfs = 0.0;     /* linfs */

    real_t bc = 0.0;
    if (id == 1) {
        /*___temperature_____________________________________________________*/
        bc = -dt * (forcing->heat_flux[n] / vcpw
                   + sval * forcing->water_flux[n] * is_nonlinfs);
    } else if (id == 2) {
        /*___salinity________________________________________________________
         * Fortran line 1523:
         *   bc_surface = dt*(virtual_salt(n) + relax_salt(n)
         *                    + real_salt_flux(n)*is_nonlinfs)
         * For linfs (us): virtual_salt is rsss·water_flux balanced;
         * relax_salt comes from SSS climatology restoring; real_salt_flux is
         * a sea-ice contribution (still 0 — sea ice not yet ported). */
        const real_t real_salt_flux = 0.0;  /* sea ice */
        bc = dt * (forcing->virtual_salt[n]
                  + forcing->relax_salt[n]
                  + real_salt_flux * is_nonlinfs);
    }
    return bc;
}

/* ------------------------------------------------------------------------ *
 * diff_ver_part_impl_ale — oce_ale_tracer.F90:562-1082                     *
 * (FCT branch: do_wimpl == .false.; Redi off: isredi == 0.)                *
 *                                                                          *
 * In Fortran the surface layer index is nzmin (=1) and the bottom layer    *
 * index is nzmax-1 (=nlevels_nod2D(n)-1). After the 1-based→0-based shift  *
 * the same identifiers refer to:                                           *
 *   nzmin = ulevels_nod2D[n] - 1   (typically 0)                           *
 *   nzmax = nlevels_nod2D[n] - 1   (the bottom interface index)            *
 * Layers run nz ∈ [nzmin, nzmax-1]; interfaces nz ∈ [nzmin, nzmax].         *
 * ------------------------------------------------------------------------ */
static void diff_ver_part_impl_ale(int                          tr_num,
                                   const struct fesom_mesh     *mesh,
                                   const struct fesom_aux      *aux,
                                   const struct fesom_forcing  *forcing,
                                   struct fesom_tracers        *tracers,
                                   const struct fesom_gm       *gm)
{
    const int    myDim_nod2D = mesh->myDim_nod2D;
    const int    nl          = mesh->nl;
    const int    nl1         = nl - 1;          /* slope_tapered second axis size */
    const real_t dt          = (real_t)FESOM_PHASE1_DT;

    /* Fortran:
     *   if ((trim(tracers%data(tr_num)%tra_adv_lim)=='FCT')
     *       .OR. (.not. dynamics%use_wsplit)) do_wimpl=.false.
     *   if (Redi) isredi=1._WP
     * For us: tra_adv_lim is FCT for both T and S → do_wimpl=.false. */
    const int    do_wimpl    = 0;
    /* Phase G5: isredi switches on iff a GM/Redi state struct is provided.
     * gm=NULL → falls back to pre-G5 numerics (Ty/Ty1 ≡ 0). */
    const real_t isredi      = (gm != NULL) ? 1.0 : 0.0;
    const real_t *st         = gm ? gm->slope_tapered : NULL;  /* [N * nl1 * 3] */
    const real_t *Ki         = gm ? gm->Ki            : NULL;  /* [N * nl] */

    real_t *trarr = tracers->data[tr_num].values;
    /* Fortran tracers%data(tr_num)%ID is 1 for T, 2 for S. Our enum:
     *   FESOM_TRACER_T = 0, FESOM_TRACER_S = 1, so id = tr_num + 1. */
    const int id = tr_num + 1;

    real_t a[NL_MAX], b[NL_MAX], c[NL_MAX], tr[NL_MAX];
    real_t cp[NL_MAX], tp[NL_MAX];
    real_t zbar_n[NL_MAX], Z_n[NL_MAX];

    for (int n = 0; n < myDim_nod2D; ++n) {
        /* initialise (Fortran lines 727-733) */
        for (int k = 0; k < NL_MAX; ++k) {
            a[k] = 0.0; b[k] = 0.0; c[k] = 0.0;
            tr[k] = 0.0; tp[k] = 0.0; cp[k] = 0.0;
            zbar_n[k] = 0.0; Z_n[k] = 0.0;
        }
        /* Fortran lines 735-736 — note the 1-based→0-based shift. */
        const int nzmax = mesh->nlevels_nod2D[n] - 1;   /* bottom interface (0-based) */
        const int nzmin = mesh->ulevels_nod2D[n] - 1;   /* top    interface (0-based) */
        if (nzmax - nzmin < 1) continue;
        FESOM_CHECK(nzmax < NL_MAX, "tracer impl_diff: nzmax %d >= NL_MAX", nzmax);

        /* zbar_n & Z_n on the NEW vertical mesh (Fortran lines 743-751).
         * In Fortran:
         *     zbar_n(nzmax)   = zbar_n_bot(n)
         *     Z_n(nzmax-1)    = zbar_n(nzmax) + hnode_new(nzmax-1,n)/2
         *     do nz=nzmax-1,nzmin+1,-1
         *         zbar_n(nz)  = zbar_n(nz+1) + hnode_new(nz,n)
         *         Z_n(nz-1)   = zbar_n(nz)   + hnode_new(nz-1,n)/2
         *     end do
         *     zbar_n(nzmin)   = zbar_n(nzmin+1) + hnode_new(nzmin,n)
         *
         * After the index shift each occurrence of the Fortran index k means
         * the 0-based slot k. zbar_n_bot(n) is the bottom depth at this node;
         * for full cells it equals mesh%zbar(nlevels_nod2D(n)) which in C is
         * mesh->zbar[nzmax]. */
        zbar_n[nzmax]   = mesh->zbar[nzmax];
        Z_n[nzmax - 1]  = zbar_n[nzmax]
                        + mesh->hnode_new[FESOM_NODE3D(n, nzmax - 1, nl)] / 2.0;
        for (int nz = nzmax - 1; nz >= nzmin + 1; --nz) {
            zbar_n[nz]    = zbar_n[nz + 1]
                          + mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];
            Z_n[nz - 1]   = zbar_n[nz]
                          + mesh->hnode_new[FESOM_NODE3D(n, nz - 1, nl)] / 2.0;
        }
        zbar_n[nzmin]   = zbar_n[nzmin + 1]
                        + mesh->hnode_new[FESOM_NODE3D(n, nzmin, nl)];

        real_t zinv1 = 0.0, zinv2 = 0.0, zinv = 0.0;

        /*___________________________________________________________________
         * Regular part of coefficients: --> surface layer  (Fortran 752-786) */
        {
            const int nz = nzmin;
            /* 1/dz(nz) (Fortran line 757) */
            zinv2 = 1.0 / (Z_n[nz] - Z_n[nz + 1]);
            zinv  = 1.0 * dt;       /* Fortran line 758 */

            /* Isoneutral K33 augmentation (Fortran 761-763, Phase G5). */
            real_t Ty1 = 0.0;
            if (gm) {
                real_t st_nz   = st[(size_t)n * nl1 * 3 + nz       * 3 + 2];
                real_t st_nz1  = st[(size_t)n * nl1 * 3 + (nz + 1) * 3 + 2];
                real_t Ki_nz   = Ki[(size_t)n * nl + nz];
                real_t Ki_nz1  = Ki[(size_t)n * nl + (nz + 1)];
                Ty1 = (Z_n[nz] - zbar_n[nz + 1]) * zinv2 * st_nz  * st_nz  * Ki_nz
                    + (zbar_n[nz + 1] - Z_n[nz + 1]) * zinv2 * st_nz1 * st_nz1 * Ki_nz1;
            }
            Ty1 = Ty1 * isredi;

            /* a/b/c at the surface (Fortran 766-769) */
            a[nz] = 0.0;
            c[nz] = -(aux->Kv[FESOM_NODE3D(n, nz + 1, nl)] + Ty1) * zinv2 * zinv
                   * mesh->area[FESOM_NODE3D(n, nz + 1, nl)]
                   / mesh->areasvol[FESOM_NODE3D(n, nz, nl)];
            b[nz] = -c[nz] + mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];

            if (do_wimpl) {
                /* Fortran 772-784 (intentionally inactive: do_wimpl == false). */
            }
            zinv1 = zinv2;          /* Fortran line 786 */
        }

        /*___________________________________________________________________
         * Regular part of coefficients: --> 2nd...nl-2 layer (Fortran 788-828) */
        for (int nz = nzmin + 1; nz <= nzmax - 2; ++nz) {
            zinv2 = 1.0 / (Z_n[nz] - Z_n[nz + 1]);

            /* Isoneutral K33 augmentation (Fortran 794-799, Phase G5). */
            real_t Ty  = 0.0;
            real_t Ty1 = 0.0;
            if (gm) {
                real_t st_m1  = st[(size_t)n * nl1 * 3 + (nz - 1) * 3 + 2];
                real_t st_nz  = st[(size_t)n * nl1 * 3 + nz       * 3 + 2];
                real_t st_nz1 = st[(size_t)n * nl1 * 3 + (nz + 1) * 3 + 2];
                real_t Ki_m1  = Ki[(size_t)n * nl + (nz - 1)];
                real_t Ki_nz  = Ki[(size_t)n * nl + nz];
                real_t Ki_nz1 = Ki[(size_t)n * nl + (nz + 1)];
                Ty  = (Z_n[nz - 1] - zbar_n[nz]) * zinv1 * st_m1 * st_m1 * Ki_m1
                    + (zbar_n[nz]  - Z_n[nz])    * zinv1 * st_nz * st_nz * Ki_nz;
                Ty1 = (Z_n[nz]     - zbar_n[nz + 1]) * zinv2 * st_nz  * st_nz  * Ki_nz
                    + (zbar_n[nz + 1] - Z_n[nz + 1]) * zinv2 * st_nz1 * st_nz1 * Ki_nz1;
            }
            Ty  = Ty  * isredi;
            Ty1 = Ty1 * isredi;

            a[nz] = -(aux->Kv[FESOM_NODE3D(n, nz, nl)] + Ty) * zinv1 * zinv
                   * (mesh->area[FESOM_NODE3D(n, nz, nl)]
                      / mesh->areasvol[FESOM_NODE3D(n, nz, nl)]);
            c[nz] = -(aux->Kv[FESOM_NODE3D(n, nz + 1, nl)] + Ty1) * zinv2 * zinv
                   * mesh->area[FESOM_NODE3D(n, nz + 1, nl)]
                   / mesh->areasvol[FESOM_NODE3D(n, nz, nl)];
            b[nz] = -a[nz] - c[nz] + mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];

            zinv1 = zinv2;

            if (do_wimpl) {
                /* Fortran 814-827 (intentionally inactive). */
            }
        }

        /*___________________________________________________________________
         * Regular part of coefficients: --> nl-1 layer (Fortran 830-859) */
        {
            const int nz = nzmax - 1;
            zinv = 1.0 * dt;        /* Fortran line 834 */

            /* Isoneutral K33 augmentation (Fortran 837-839, Phase G5). */
            real_t Ty = 0.0;
            if (gm) {
                real_t st_m1 = st[(size_t)n * nl1 * 3 + (nz - 1) * 3 + 2];
                real_t st_nz = st[(size_t)n * nl1 * 3 + nz       * 3 + 2];
                real_t Ki_m1 = Ki[(size_t)n * nl + (nz - 1)];
                real_t Ki_nz = Ki[(size_t)n * nl + nz];
                Ty = (Z_n[nz - 1] - zbar_n[nz]) * zinv1 * st_m1 * st_m1 * Ki_m1
                   + (zbar_n[nz]  - Z_n[nz])    * zinv1 * st_nz * st_nz * Ki_nz;
            }
            Ty = Ty * isredi;

            a[nz] = -(aux->Kv[FESOM_NODE3D(n, nz, nl)] + Ty) * zinv1 * zinv
                   * (mesh->area[FESOM_NODE3D(n, nz, nl)]
                      / mesh->areasvol[FESOM_NODE3D(n, nz, nl)]);
            c[nz] = 0.0;
            b[nz] = -a[nz] + mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];

            if (do_wimpl) {
                /* Fortran 851-859 (intentionally inactive). */
            }
        }

        /*___________________________________________________________________
         * RHS construction (Fortran 862-884) */
        {
            const int nz = nzmin;
            const real_t dz = mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];
            tr[nz] = -(b[nz] - dz) * trarr[FESOM_NODE3D(n, nz,     nl)]
                    - c[nz]        * trarr[FESOM_NODE3D(n, nz + 1, nl)];
        }
        for (int nz = nzmin + 1; nz <= nzmax - 2; ++nz) {
            const real_t dz = mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];
            tr[nz] = -a[nz]        * trarr[FESOM_NODE3D(n, nz - 1, nl)]
                    -(b[nz] - dz)  * trarr[FESOM_NODE3D(n, nz,     nl)]
                    -c[nz]         * trarr[FESOM_NODE3D(n, nz + 1, nl)];
        }
        {
            const int nz = nzmax - 1;
            const real_t dz = mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];
            tr[nz] = -a[nz]       * trarr[FESOM_NODE3D(n, nz - 1, nl)]
                    -(b[nz] - dz) * trarr[FESOM_NODE3D(n, nz,     nl)];
        }

        /*___________________________________________________________________
         * KPP nonlocal flux contributions (Fortran 892-1015) — not used:
         * use_kpp_nonlclflx == .false. */

        /*___________________________________________________________________
         * Surface BC (Fortran 1018-1030):
         *   tr(nzmin) = tr(nzmin) + bc_surface(...) */
        {
            const int    nz   = nzmin;
            const real_t sval = trarr[FESOM_NODE3D(n, nz, nl)];
            tr[nz] += bc_surface(n, id, sval, forcing);
        }

        /*___________________________________________________________________
         * Shortwave penetration (Fortran oce_ale_tracer.F90:990): for the
         * temperature tracer (id==1), add the sw_3d flux divergence per layer.
         * zinv = dt for linfs (the /(zbar(nz)-zbar(nz+1)) is commented out in
         * the ALE Fortran). sw_3d is built by fesom_cal_shortwave_rad. */
        if (id == 1 && FESOM_PHASE1_USE_SW_PENE) {
            const real_t *sw = forcing->sw_3d;
            const real_t  dtl = (real_t)FESOM_PHASE1_DT;
            for (int nz = nzmin; nz <= nzmax - 1; ++nz) {
                const real_t top = sw[FESOM_NODE3D(n, nz,     nl)];
                const real_t bot = sw[FESOM_NODE3D(n, nz + 1, nl)];
                const real_t ar  = mesh->area[FESOM_NODE3D(n, nz + 1, nl)]
                                 / mesh->areasvol[FESOM_NODE3D(n, nz, nl)];
                tr[nz] += (top - bot * ar) * dtl;
            }
        }

        /*___________________________________________________________________
         * Forward sweep (Fortran 1052-1061) */
        cp[nzmin] = c[nzmin] / b[nzmin];
        tp[nzmin] = tr[nzmin] / b[nzmin];
        for (int nz = nzmin + 1; nz <= nzmax - 1; ++nz) {
            const real_t m = b[nz] - cp[nz - 1] * a[nz];
            cp[nz] = c[nz] / m;
            tp[nz] = (tr[nz] - tp[nz - 1] * a[nz]) / m;
        }

        /*___________________________________________________________________
         * Back substitution (Fortran 1063-1069) */
        tr[nzmax - 1] = tp[nzmax - 1];
        for (int nz = nzmax - 2; nz >= nzmin; --nz) {
            tr[nz] = tp[nz] - cp[nz] * tr[nz + 1];
        }

        /*___________________________________________________________________
         * Update tracer (Fortran 1073-1077) */
        for (int nz = nzmin; nz <= nzmax - 1; ++nz) {
            trarr[FESOM_NODE3D(n, nz, nl)] += tr[nz];
        }
    }
}

/* ------------------------------------------------------------------------ *
 * Public driver — call once per tracer per timestep.                       *
 * ------------------------------------------------------------------------ */
void fesom_impl_vert_diff_tracers(const struct fesom_mesh    *mesh,
                                  const struct fesom_aux     *aux,
                                  const struct fesom_forcing *forcing,
                                  struct fesom_tracers       *tracers,
                                  const struct fesom_gm      *gm)
{
    diff_ver_part_impl_ale(FESOM_TRACER_T, mesh, aux, forcing, tracers, gm);
    diff_ver_part_impl_ale(FESOM_TRACER_S, mesh, aux, forcing, tracers, gm);
}
