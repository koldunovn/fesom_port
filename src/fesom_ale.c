#include "fesom_ale.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_mesh.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void fesom_ale_thickness_linfs(struct fesom_mesh *mesh)
{
    /* hnode_new = hnode in linfs. Memcpy is the literal port. */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    size_t bytes = (size_t)N * (size_t)mesh->nl * sizeof(real_t);
    memcpy(mesh->hnode_new, mesh->hnode, bytes);
}

void fesom_ale_commit_thickness(struct fesom_mesh *mesh)
{
    const int nl = mesh->nl;
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    /* helem on interior elements; halo elements get values via the
     * exchange in fesom_step. */
    int Ei = mesh->myDim_elem2D;
    /* hnode := hnode_new */
    size_t bytes = (size_t)N * (size_t)nl * sizeof(real_t);
    memcpy(mesh->hnode, mesh->hnode_new, bytes);

    for (int e = 0; e < Ei; ++e) {
        int nzmin = mesh->ulevels[e]   - 1;
        int nzmax = mesh->nlevels[e]   - 1;
        int n0 = mesh->elem_nodes[3*e + 0];
        int n1 = mesh->elem_nodes[3*e + 1];
        int n2 = mesh->elem_nodes[3*e + 2];
        for (int nz = nzmin; nz < nzmax; ++nz) {
            mesh->helem[FESOM_ELEM3D(e, nz, nl)] =
                (mesh->hnode[FESOM_NODE3D(n0, nz, nl)] +
                 mesh->hnode[FESOM_NODE3D(n1, nz, nl)] +
                 mesh->hnode[FESOM_NODE3D(n2, nz, nl)]) / 3.0;
        }
    }
}

/*
 * Vertical velocity from horizontal divergence — linfs branch only.
 *
 * Algorithm (vert_vel_ale lines 2189-2346, linfs path):
 *
 * 1. Zero w at all nodes.
 *
 * 2. For each edge ed:
 *      let n1, n2     = edges[ed]
 *      let el1, el2   = edge_tri[ed]   (el2 may be -1 for boundary)
 *      let dxL, dyL   = edge_cross_dxdy[ed][0..1]   (centroid(el1) − edge_mid)
 *      let dxR, dyR   = edge_cross_dxdy[ed][2..3]
 *
 *      For el1 (always present):
 *        c1[nz] = (uv[el1].v * dxL − uv[el1].u * dyL) * helem[el1][nz]
 *        w[n1][nz] += c1[nz]
 *        w[n2][nz] −= c1[nz]
 *
 *      For el2 if present (note negation of dxR,dyR sign):
 *        c1[nz] = -(uv[el2].v * dxR − uv[el2].u * dyR) * helem[el2][nz]
 *        w[n1][nz] += c1[nz]
 *        w[n2][nz] −= c1[nz]
 *
 * 3. Cumulative sum upward (bottom → top):
 *      for nz from nzmax-1 down to nzmin:
 *          w[n][nz] += w[n][nz + 1]
 *
 * 4. Divide by area(nz, n) → m/s.
 *
 * For Phase 1 the upper layer surface-flux correction (zlevel/zstar branches)
 * is omitted — linfs has dh/dt = 0, so w is the physical vertical velocity
 * directly.
 */
void fesom_ale_vert_vel_linfs(const struct fesom_mesh *mesh,
                              struct fesom_dyn        *dyn,
                              int                      gm_on)
{
    const int nl = mesh->nl;

    /* Step 1: zero w everywhere (all local nodes); fer_w too if GM on. */
    int N = mesh->myDim_nod2D + mesh->eDim_nod2D;
    /* Fortran vert_vel_ale (oce_ale.F90:2200) iterates `do ed=1, myDim_edge2D`
     * only. The earlier "myDim+eDim" was wrong — iterating halo edges adds
     * contributions to interior boundary nodes that the owner's loop also
     * adds → small per-step error → CG NaN after ~85 steps. */
    int EG = mesh->myDim_edge2D;
    memset(dyn->w, 0, (size_t)N * (size_t)nl * sizeof(real_t));
    if (gm_on) {
        memset(dyn->fer_w, 0, (size_t)N * (size_t)nl * sizeof(real_t));
    }

    /* Step 2: edge-flux accumulation. Each interior boundary node receives
     * contributions from all its incident edges via this rank's myDim_edge2D
     * (the partition replicates cross-rank edges); exchange_nod(w) afterward
     * refreshes halo entries with their owners' values.
     *
     * Phase G6a: when gm_on, the same edge iteration also accumulates fer_w
     * from dyn->fer_uv. The c1 (regular) and c2 (GM) contributions write to
     * dyn->w and dyn->fer_w respectively — disjoint targets, so the dyn->w
     * accumulator is bit-identical to the gm_on=0 path. */
    for (int ed = 0; ed < EG; ++ed) {
        int n1 = mesh->edges[2*ed + 0];
        int n2 = mesh->edges[2*ed + 1];
        int el1 = mesh->edge_tri[2*ed + 0];
        int el2 = mesh->edge_tri[2*ed + 1];

        /* Left cell — always present in a well-formed mesh */
        if (el1 >= 0) {
            real_t dxL = mesh->edge_cross_dxdy[4*ed + 0];
            real_t dyL = mesh->edge_cross_dxdy[4*ed + 1];
            int nzmin = mesh->ulevels[el1] - 1;
            int nzmax = mesh->nlevels[el1] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                /* Mirror of c1=(UV(2)*dx - UV(1)*dy)*helem; UV index 1=u, 2=v. */
                real_t u = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
                real_t v = dyn->uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
                real_t h = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
                real_t c1 = (v * dxL - u * dyL) * h;
                dyn->w[FESOM_NODE3D(n1, nz, nl)] += c1;
                dyn->w[FESOM_NODE3D(n2, nz, nl)] -= c1;
                if (gm_on) {
                    real_t fu = dyn->fer_uv[FESOM_ELEMVEC(el1, nz, nl) + 0];
                    real_t fv = dyn->fer_uv[FESOM_ELEMVEC(el1, nz, nl) + 1];
                    real_t c2 = (fv * dxL - fu * dyL) * h;
                    dyn->fer_w[FESOM_NODE3D(n1, nz, nl)] += c2;
                    dyn->fer_w[FESOM_NODE3D(n2, nz, nl)] -= c2;
                }
            }
        }
        /* Right cell — only on interior edges */
        if (el2 >= 0) {
            real_t dxR = mesh->edge_cross_dxdy[4*ed + 2];
            real_t dyR = mesh->edge_cross_dxdy[4*ed + 3];
            int nzmin = mesh->ulevels[el2] - 1;
            int nzmax = mesh->nlevels[el2] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                real_t u = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v = dyn->uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                real_t c1 = -(v * dxR - u * dyR) * h;
                dyn->w[FESOM_NODE3D(n1, nz, nl)] += c1;
                dyn->w[FESOM_NODE3D(n2, nz, nl)] -= c1;
                if (gm_on) {
                    real_t fu = dyn->fer_uv[FESOM_ELEMVEC(el2, nz, nl) + 0];
                    real_t fv = dyn->fer_uv[FESOM_ELEMVEC(el2, nz, nl) + 1];
                    real_t c2 = -(fv * dxR - fu * dyR) * h;
                    dyn->fer_w[FESOM_NODE3D(n1, nz, nl)] += c2;
                    dyn->fer_w[FESOM_NODE3D(n2, nz, nl)] -= c2;
                }
            }
        }
    }

    /* Step 3: cumulative sum bottom → top (all local nodes; halo nodes will
     * be re-synced by exchange_nod(w) in the timestep wrapper). Mirror for
     * fer_w when gm_on. */
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmax - 1; nz >= nzmin; --nz) {
            dyn->w[FESOM_NODE3D(n, nz, nl)] +=
                dyn->w[FESOM_NODE3D(n, nz + 1, nl)];
            if (gm_on) {
                dyn->fer_w[FESOM_NODE3D(n, nz, nl)] +=
                    dyn->fer_w[FESOM_NODE3D(n, nz + 1, nl)];
            }
        }
    }

    /* Step 4: divide by area(nz, n) → m/s */
    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t a = mesh->area[FESOM_NODE3D(n, nz, nl)];
            if (a > 0.0) {
                dyn->w[FESOM_NODE3D(n, nz, nl)] /= a;
                if (gm_on) {
                    dyn->fer_w[FESOM_NODE3D(n, nz, nl)] /= a;
                }
            }
        }
    }
}

/*===========================================================================
 * compute_CFLz (oce_ale.F90:2935-3022)
 *
 * At each interface nz the vertical CFL accumulates contributions from
 * both adjacent layers:
 *   iter nz=k: CFL[k]   += |W[k]|·dt / h[k]
 *              CFL[k+1]  = |W[k+1]|·dt / h[k]       (overwrite — first touch)
 *   iter nz=k+1: CFL[k+1] += |W[k+1]|·dt / h[k+1]
 * Net: CFL[interior k] = |W[k]|·dt · (1/h[k-1] + 1/h[k])
 *      CFL[surface]   = |W[nzmin]|·dt / h[nzmin]
 *      CFL[bottom]    = |W[nzmax]|·dt / h[nzmax-1]
 * We zero cfl_z first and use += throughout (mathematically equivalent).
 *===========================================================================*/

void fesom_ale_compute_cflz(const struct fesom_mesh *mesh,
                            struct fesom_dyn        *dyn)
{
    const int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl = mesh->nl;
    const real_t dt = (real_t)FESOM_PHASE1_DT;

    memset(dyn->cfl_z, 0, (size_t)N * (size_t)nl * sizeof(real_t));

    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = nzmin; nz < nzmax; ++nz) {
            real_t h = mesh->hnode_new[FESOM_NODE3D(n, nz, nl)];
            if (h <= 0.0) continue;
            real_t w_top = dyn->w[FESOM_NODE3D(n, nz,     nl)];
            real_t w_bot = dyn->w[FESOM_NODE3D(n, nz + 1, nl)];
            real_t c1 = fabs(w_top) * dt / h;
            real_t c2 = fabs(w_bot) * dt / h;
            dyn->cfl_z[FESOM_NODE3D(n, nz,     nl)] += c1;
            dyn->cfl_z[FESOM_NODE3D(n, nz + 1, nl)] += c2;
        }
    }
}

/*===========================================================================
 * compute_Wvel_split (oce_ale.F90:3026-3074)
 *
 *   w_e, w_i : explicit and implicit parts of w
 *   For each (n, nz):
 *     w_e[n, nz] = w[n, nz]
 *     w_i[n, nz] = 0
 *     if (use_wsplit && cfl_z > maxcfl):
 *       dd     = max(cfl_z − maxcfl, 0) / max(maxcfl, 1e-12)
 *       w_e   *= 1 / (1 + dd)
 *       w_i    = dd / (1 + dd) · w
 *===========================================================================*/

void fesom_ale_compute_wvel_split(const struct fesom_mesh *mesh,
                                  struct fesom_dyn        *dyn)
{
    const int N  = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl = mesh->nl;
    const real_t maxcfl     = (real_t)FESOM_PHASE1_WSPLIT_MAXCFL;
    const real_t use_wsplit = (real_t)FESOM_PHASE1_USE_WSPLIT;
    const real_t inv_maxcfl = 1.0 / ((maxcfl > 1e-12) ? maxcfl : 1e-12);

    for (int n = 0; n < N; ++n) {
        int nzmin = mesh->ulevels_nod2D[n] - 1;
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        /* Loop over interfaces: nzmin..nzmax inclusive */
        for (int nz = nzmin; nz <= nzmax; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, nl);
            real_t w   = dyn->w[k];
            real_t cfl = dyn->cfl_z[k];
            if (use_wsplit > 0.5 && cfl > maxcfl) {
                real_t excess = cfl - maxcfl;
                real_t dd = (excess > 0.0 ? excess : 0.0) * inv_maxcfl;
                real_t inv_1_dd = 1.0 / (1.0 + dd);
                dyn->w_e[k] = w * inv_1_dd;
                dyn->w_i[k] = w * dd * inv_1_dd;
            } else {
                dyn->w_e[k] = w;
                dyn->w_i[k] = 0.0;
            }
        }
    }
}
