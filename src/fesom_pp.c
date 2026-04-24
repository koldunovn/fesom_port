/*
 * PP mixing + convective adjustment. Literal ports — see header for
 * Fortran source line references.
 */
#include "fesom_pp.h"
#include "fesom_aux.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_mesh.h"

#include <math.h>
#include <stdlib.h>

/*===========================================================================
 * compute_vel_nodes (oce_dyn.F90:177-226)
 *
 * For each node n and each layer nz:
 *   tx, ty, tvol = 0
 *   for each surrounding element el (via nod_in_elem2D):
 *       if nz outside (ulevels(el) .. nlevels(el)-1): skip
 *       tvol += elem_area(el)
 *       tx   += UV(1, nz, el) * elem_area(el)
 *       ty   += UV(2, nz, el) * elem_area(el)
 *   UVnode(1, nz, n) = tx / tvol
 *   UVnode(2, nz, n) = ty / tvol
 *===========================================================================*/

void fesom_compute_vel_nodes(const struct fesom_mesh *mesh,
                             struct fesom_dyn        *dyn)
{
    const int N    = mesh->myDim_nod2D;       /* interior only — halo gets exchange */
    const int nl   = mesh->nl;

    for (int n = 0; n < N; ++n) {
        int uln = mesh->ulevels_nod2D[n] - 1;     /* 1-based → 0-based */
        int nln = mesh->nlevels_nod2D[n] - 1;     /* exclusive bound  */
        for (int nz = uln; nz < nln; ++nz) {
            real_t tvol = 0.0, tx = 0.0, ty = 0.0;
            int o0 = mesh->nod_in_elem2D_offsets[n];
            int o1 = mesh->nod_in_elem2D_offsets[n + 1];
            for (int k = o0; k < o1; ++k) {
                int e   = mesh->nod_in_elem2D[k];
                int ule = mesh->ulevels[e] - 1;
                int nle = mesh->nlevels[e] - 1;
                if (nz < ule || nz >= nle) continue;
                real_t a = mesh->elem_area[e];
                tvol += a;
                tx   += dyn->uv[FESOM_ELEMVEC(e, nz, nl) + 0] * a;
                ty   += dyn->uv[FESOM_ELEMVEC(e, nz, nl) + 1] * a;
            }
            if (tvol > 0.0) {
                dyn->uvnode[FESOM_ELEMVEC(n, nz, nl) + 0] = tx / tvol;
                dyn->uvnode[FESOM_ELEMVEC(n, nz, nl) + 1] = ty / tvol;
            } else {
                dyn->uvnode[FESOM_ELEMVEC(n, nz, nl) + 0] = 0.0;
                dyn->uvnode[FESOM_ELEMVEC(n, nz, nl) + 1] = 0.0;
            }
        }
    }
}

/*===========================================================================
 * oce_mixing_pp (oce_ale_mixing_pp.F90:2-93)
 *
 * Three loops:
 *   1. Per node, per interior interface nz ∈ [nzmin+1, nzmax-1):
 *        dz_inv = 1 / (Z(nz-1) - Z(nz))             ; > 0 because Z negative-down
 *        shear  = (Δu)² + (Δv)²                     ; node values, not element
 *        shear *= dz_inv²
 *        Kv(nz, node) = shear / (shear + 5*max(N², 0) + 1e-14)
 *      ← THIS IS THE "factor" — overwrites Kv with the dimensionless ratio.
 *
 *   2. Per element, per interior interface nz:
 *        Av(nz, elem) = mix_coeff_PP * Σ Kv(nz, elnodes)² / 3 + A_ver
 *      ← Reads `Kv` BEFORE it has been raised to power 3 and offset by K_ver.
 *
 *   3. Per node again (Kv0_const branch):
 *        Kv(nz, node) = mix_coeff_PP * Kv(nz, node)³ + K_ver
 *      ← Order of (2) and (3) matters: Av is built from the (factor)² of Kv,
 *        before Kv itself is raised to (factor)³ and offset.
 *
 * For Phase 1 with no partial cells, Z_3d_n[nz][n] = Z[nz] (1D mid-level depth).
 *===========================================================================*/

void fesom_pp_mixing(const struct fesom_mesh *mesh,
                     const struct fesom_dyn  *dyn,
                     struct fesom_aux        *aux)
{
    /* Fortran oce_ale_mixing_pp.F90:41,72 iterates `myDim+eDim_nod2D` for the
     * Kv loops (compute Kv at all local nodes, including halo). The element
     * Av loop reads Kv at all 3 vertices — if halo Kv is stale (only myDim
     * computed locally), Av on boundary elements is wrong → wrong vertical
     * mixing → tracer divergence at partition boundaries → S explosion. */
    const int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int E  = mesh->myDim_elem2D;
    const int nl = mesh->nl;
    const real_t mix_coeff = (real_t)FESOM_PHASE1_MIX_COEFF_PP;
    const real_t K_bg      = (real_t)FESOM_PHASE1_K_VER;
    const real_t A_bg      = (real_t)FESOM_PHASE1_A_VER;

    /* Loop 1: factor = shear² / (shear² + 5N² + ε) → write into Kv */
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t dz = mesh->Z[nz - 1] - mesh->Z[nz];      /* > 0 */
            real_t dz_inv = 1.0 / dz;
            real_t du = dyn->uvnode[FESOM_ELEMVEC(n, nz - 1, nl) + 0]
                      - dyn->uvnode[FESOM_ELEMVEC(n, nz,     nl) + 0];
            real_t dv = dyn->uvnode[FESOM_ELEMVEC(n, nz - 1, nl) + 1]
                      - dyn->uvnode[FESOM_ELEMVEC(n, nz,     nl) + 1];
            real_t shear = (du * du + dv * dv) * dz_inv * dz_inv;
            real_t Nsq   = aux->bvfreq[FESOM_NODE3D(n, nz, nl)];
            real_t Nsq_pos = Nsq > 0.0 ? Nsq : 0.0;
            aux->Kv[FESOM_NODE3D(n, nz, nl)]
                = shear / (shear + 5.0 * Nsq_pos + 1.0e-14);
        }
    }

    /* Loop 2: Av at elements from (Kv-as-factor)² averaged over 3 vertices */
    for (int e = 0; e < E; ++e) {
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];
        int nzmin = mesh->ulevels[e] - 1;
        int nzmax = mesh->nlevels[e] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t k0 = aux->Kv[FESOM_NODE3D(n0, nz, nl)];
            real_t k1 = aux->Kv[FESOM_NODE3D(n1, nz, nl)];
            real_t k2 = aux->Kv[FESOM_NODE3D(n2, nz, nl)];
            aux->Av[FESOM_ELEM3D(e, nz, nl)]
                = mix_coeff * (k0*k0 + k1*k1 + k2*k2) / 3.0 + A_bg;
        }
    }

    /* Loop 3: Kv = mix_coeff * factor³ + K_ver  (Kv0_const branch) */
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t f = aux->Kv[FESOM_NODE3D(n, nz, nl)];
            aux->Kv[FESOM_NODE3D(n, nz, nl)] = mix_coeff * f * f * f + K_bg;
        }
    }
}

/*===========================================================================
 * mo_convect — use_instabmix branch only (oce_mo_conv.F90:79-94 for Kv,
 * 99-120 for Av). Phase 1 skips momix and windmix entirely.
 *===========================================================================*/

void fesom_mo_convect(const struct fesom_mesh *mesh,
                      struct fesom_aux        *aux)
{
    if (!FESOM_PHASE1_USE_INSTABMIX) return;

    /* Fortran oce_mo_conv.F90:36,79 iterates `myDim+eDim_nod2D` for the Kv
     * loop. The Av loop reads Kv at element vertices, so halo Kv must be
     * locally fresh — same reasoning as pp_mixing. */
    const int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int E  = mesh->myDim_elem2D;
    const int nl = mesh->nl;
    const real_t imix_kv = (real_t)FESOM_PHASE1_INSTABMIX_KV;

    /* Kv: where bvfreq < 0, set Kv = max(Kv, instabmix_kv).
       Mirror of oce_mo_conv.F90:85. */
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            if (aux->bvfreq[k] < 0.0) {
                if (aux->Kv[k] < imix_kv) aux->Kv[k] = imix_kv;
            }
        }
    }

    /* Av: where ANY of the three vertex bvfreq < 0, set Av = max(Av, imix_kv).
       Mirror of oce_mo_conv.F90:108. */
    for (int e = 0; e < E; ++e) {
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];
        int nzmin = mesh->ulevels[e] - 1;
        int nzmax = mesh->nlevels[e] - 1;
        for (int nz = nzmin + 1; nz < nzmax; ++nz) {
            real_t b0 = aux->bvfreq[FESOM_NODE3D(n0, nz, nl)];
            real_t b1 = aux->bvfreq[FESOM_NODE3D(n1, nz, nl)];
            real_t b2 = aux->bvfreq[FESOM_NODE3D(n2, nz, nl)];
            if (b0 < 0.0 || b1 < 0.0 || b2 < 0.0) {
                size_t ke = FESOM_ELEM3D(e, nz, nl);
                if (aux->Av[ke] < imix_kv) aux->Av[ke] = imix_kv;
            }
        }
    }
}
