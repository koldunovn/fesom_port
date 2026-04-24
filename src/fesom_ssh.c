/*
 * SSH stiffness matrix + diagonal-symmetric preconditioner + CG iteration.
 * Literal port of Fortran:
 *   - init_stiff_mat_ale       (oce_ale.F90:1393-1684)
 *   - ssh_solve_preconditioner (solver.F90:31-95)
 *   - compute_ssh_rhs_ale      (oce_ale.F90:1821-1956, linfs branch)
 *   - ssh_solve_cg             (solver.F90:98-281)
 *
 * Phase 1 simplifications (literal port — these branches are explicitly the
 * Fortran behaviour when the corresponding flags are off):
 *   - linfs ALE: update_stiff_mat_ale is NOT called (gated at oce_ale.F90:3722).
 *     The stiffness matrix is built once and reused every step.
 *   - No partial cells: zbar_e_srf[e] = zbar[0] = 0, zbar_e_bot[e] = zbar[nlevels[e]-1].
 *   - No cavity: ulevels_nod2D = ulevels = 1 → 0-based 0 → diagonal mass term
 *     is added at every node.
 *   - Serial (1 MPI rank): exchange_nod and MPI_Allreduce are no-ops.
 */
#include "fesom_ssh.h"
#include "fesom_constants.h"
#include "fesom_dyn.h"
#include "fesom_mesh.h"

#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "fesom_halo.h"
#include "fesom_partit.h"

/*===========================================================================
 * Stiffness matrix: build CSR + fill values
 *===========================================================================*/

void fesom_ssh_stiff_alloc_and_build(fesom_ssh_stiff       *S,
                                     const struct fesom_mesh *mesh)
{
    /* Stiffness matrix has myDim rows (interior nodes only); column indices
     * reach into the halo (myDim..myDim+eDim) for SpMV. Mirrors Fortran
     * where ssh_stiff is sized myDim_nod2D. */
    const int N         = mesh->myDim_nod2D;
    const int N_alloc   = mesh->myDim_nod2D + mesh->eDim_nod2D;
    /* Loop over INTERIOR edges only — Fortran oce_ale.F90:1517 has
     * `do ed=1,myDim_edge2D` with the comment "!! Attention". Edges in eDim
     * belong to a neighbour rank and are accumulated by the owner. The owner's
     * interior loop adds to its myDim row; the same row is in our eDim_nod2D
     * (halo), where we don't have a matrix row. So skipping halo edges here
     * exactly matches Fortran. */
    const int E         = mesh->myDim_edge2D;
    memset(S, 0, sizeof(*S));
    S->dim = N;
    (void)N_alloc;

    /* ---- 1. Per-row neighbour count and local-id list ---------------------
     * Fortran lines 1435-1463. n_pos(1, n) = n is set first so the diagonal
     * is the first entry in each row's neighbour list. */
    /* n_num: counter per row in pass 1, and node-id → sparse-position lookup
     * in pass 2 (used as a sparse hash, so it must be sized for ALL local
     * nodes — interior + halo — since column indices reach into halo). */
    int *n_num = calloc((size_t)N_alloc, sizeof(int));
    int  max_neighbors = 12;     /* Fortran uses 12 — same hard cap */
    int *n_pos = calloc((size_t)N * (size_t)max_neighbors, sizeof(int));
    FESOM_CHECK(n_num && n_pos, "ssh_stiff: out of memory (neighbour build)");

    for (int n = 0; n < N; ++n) {
        n_num[n] = 1;
        n_pos[n * max_neighbors + 0] = n;
    }
    for (int ed = 0; ed < E; ++ed) {
        int n1 = mesh->edges[2*ed + 0];
        int n2 = mesh->edges[2*ed + 1];
        /* Only add edges where the row endpoint is INTERIOR (the matrix has
         * no halo-row entries); the column endpoint can be interior or halo. */
        if (n1 < N) {
            FESOM_CHECK(n_num[n1] < max_neighbors,
                        "ssh_stiff: node %d has > %d neighbours", n1, max_neighbors);
            n_pos[n1 * max_neighbors + n_num[n1]++] = n2;
        }
        if (n2 < N) {
            FESOM_CHECK(n_num[n2] < max_neighbors,
                        "ssh_stiff: node %d has > %d neighbours", n2, max_neighbors);
            n_pos[n2 * max_neighbors + n_num[n2]++] = n1;
        }
    }

    /* ---- 2. Build CSR rowptr (lines 1467-1477) ---------------------------- */
    S->rowptr = malloc((size_t)(N + 1) * sizeof(int));
    FESOM_CHECK(S->rowptr, "ssh_stiff: out of memory (rowptr)");
    S->rowptr[0] = 0;
    for (int n = 0; n < N; ++n) {
        S->rowptr[n + 1] = S->rowptr[n] + n_num[n];
    }
    S->nnz = S->rowptr[N];

    /* ---- 3. Fill colind + zero values + alloc pr_values (lines 1481-1494) - */
    S->colind    = malloc((size_t)S->nnz * sizeof(int));
    S->values    = calloc((size_t)S->nnz, sizeof(real_t));
    S->pr_values = calloc((size_t)S->nnz, sizeof(real_t));
    FESOM_CHECK(S->colind && S->values && S->pr_values,
                "ssh_stiff: out of memory (CSR arrays)");
    for (int n = 0; n < N; ++n) {
        int start = S->rowptr[n];
        for (int k = 0; k < n_num[n]; ++k) {
            S->colind[start + k] = n_pos[n * max_neighbors + k];
        }
    }
    free(n_pos);

    /* ---- 4. Stiffness term per edge × adjacent cell (lines 1503-1571) -----
     * factor = g * dt * alpha * theta
     * For each edge, for each adjacent cell el(i):
     *   fy[k] = (zbar_e_bot[el(i)] - zbar_e_srf[el(i)]) *
     *           ( gradient_sca[1..3](el(i)) * edge_cross_dxdy[2*i  ](ed)
     *           - gradient_sca[4..6](el(i)) * edge_cross_dxdy[2*i-1](ed) )
     *   if i==2: fy = -fy
     *   row = edges(1, ed): values[npos] += fy * factor
     *   row = edges(2, ed): values[npos] -= fy * factor
     * with npos = positions of elnodes(1..3) in row's neighbour list. */
    const real_t factor = (real_t)FESOM_G * (real_t)FESOM_PHASE1_DT
                        * (real_t)FESOM_PHASE1_ALPHA * (real_t)FESOM_PHASE1_THETA;

    /* n_num is reused as a "node-id → sparse-position-in-current-row" lookup
     * (Fortran line 1507 zeros it). Must reset across the FULL local extent. */
    for (int n = 0; n < N_alloc; ++n) n_num[n] = 0;

    for (int ed = 0; ed < E; ++ed) {
        int el[2] = { mesh->edge_tri[2*ed + 0], mesh->edge_tri[2*ed + 1] };
        int e_n[2] = { mesh->edges[2*ed + 0], mesh->edges[2*ed + 1] };
        for (int i = 0; i < 2; ++i) {
            /* For ed in myDim_edge2D, the FESOM partition guarantees both
             * adjacent triangles are in myDim_elem2D — so gradient_sca and
             * elem_nodes (both sized myDim_elem2D in Fortran) are always
             * safe to read. Skip true boundary (-1) and any element id past
             * myDim_elem2D as defensive (would indicate partition oddity). */
            if (el[i] < 0 || el[i] >= mesh->myDim_elem2D) continue;
            int en[3] = { mesh->elem_nodes[3*el[i] + 0],
                          mesh->elem_nodes[3*el[i] + 1],
                          mesh->elem_nodes[3*el[i] + 2] };
            const real_t *g = &mesh->gradient_sca[6 * el[i]];

            /* zbar_e_srf and zbar_e_bot (Phase 1: no partial cells, no cavity) */
            real_t zbar_e_srf = mesh->zbar[0];                            /* = 0 */
            real_t zbar_e_bot = mesh->zbar[mesh->nlevels[el[i]] - 1];     /* < 0 */
            real_t depth_diff = zbar_e_bot - zbar_e_srf;                  /* < 0 */

            /* edge_cross_dxdy[4*ed + 2*i + 0] = dx_i (Fortran 2*i-1)
               edge_cross_dxdy[4*ed + 2*i + 1] = dy_i (Fortran 2*i)
               Recall Fortran 1-based i ∈ {1,2}: 2*i-1 ∈ {1,3} → C 0-based offset {0,2};
                                                  2*i   ∈ {2,4} → C 0-based offset {1,3}. */
            real_t dx_i = mesh->edge_cross_dxdy[4*ed + 2*i + 0];
            real_t dy_i = mesh->edge_cross_dxdy[4*ed + 2*i + 1];

            real_t fy[3];
            for (int k = 0; k < 3; ++k) {
                fy[k] = depth_diff * (g[k] * dy_i - g[3 + k] * dx_i);
            }
            if (i == 1) {  /* Fortran i==2 → C i==1 */
                fy[0] = -fy[0]; fy[1] = -fy[1]; fy[2] = -fy[2];
            }

            /* Row = edges(1, ed) → C e_n[0]. Only build matrix entries for
             * interior rows; halo-row contributions are summed on the
             * owning rank instead. */
            if (e_n[0] < N) {
                int row = e_n[0];
                int rstart = S->rowptr[row];
                int rend   = S->rowptr[row + 1];
                for (int n = rstart; n < rend; ++n) {
                    n_num[S->colind[n]] = n;
                }
                int npos[3] = { n_num[en[0]], n_num[en[1]], n_num[en[2]] };
                for (int k = 0; k < 3; ++k) {
                    S->values[npos[k]] += fy[k] * factor;
                }
            }
            if (e_n[1] < N) {
                int row = e_n[1];
                int rstart = S->rowptr[row];
                int rend   = S->rowptr[row + 1];
                for (int n = rstart; n < rend; ++n) {
                    n_num[S->colind[n]] = n;
                }
                int npos[3] = { n_num[en[0]], n_num[en[1]], n_num[en[2]] };
                for (int k = 0; k < 3; ++k) {
                    S->values[npos[k]] -= fy[k] * factor;
                }
            }
        }
    }

    /* ---- 5. Mass matrix term on diagonal (lines 1575-1582) ---------------
     * For non-cavity nodes (ulevels_nod2D[row] == 1, all of them in Phase 1):
     *   values[diagonal] += areasvol[ulevels_nod2D[row], row] / dt
     * In C 0-based with ulevels_nod2D[row]=1, the surface-layer index is 0. */
    const real_t inv_dt = 1.0 / (real_t)FESOM_PHASE1_DT;
    for (int row = 0; row < N; ++row) {
        if (mesh->ulevels_nod2D[row] > 1) continue;        /* cavity guard */
        int diag = S->rowptr[row];
        int top_layer = mesh->ulevels_nod2D[row] - 1;       /* = 0 in Phase 1 */
        S->values[diag] += mesh->areasvol[FESOM_NODE3D(row, top_layer, mesh->nl)] * inv_dt;
    }

    free(n_num);
}

void fesom_ssh_stiff_free(fesom_ssh_stiff *S)
{
    free(S->rowptr);
    free(S->colind);
    free(S->values);
    free(S->pr_values);
    memset(S, 0, sizeof(*S));
}

/*===========================================================================
 * Preconditioner: ssh_solve_preconditioner (solver.F90:31-95)
 *===========================================================================*/

void fesom_ssh_preconditioner(fesom_ssh_stiff *S, const struct fesom_mesh *mesh,
                              struct fesom_partit *partit)
{
    const int N       = mesh->myDim_nod2D;
    const int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;

    /* Collect diagonal values. Sized for full local extent because column
     * indices in `S->colind` can refer to halo nodes (n >= N), and the
     * preconditioner formula reads diag_values[node] for every off-diagonal.
     * Halo entries are filled by halo exchange below (solver.F90 mirrors). */
    real_t *diag_values = calloc((size_t)N_alloc, sizeof(real_t));
    FESOM_CHECK(diag_values, "preconditioner: out of memory");
    for (int row = 0; row < N; ++row) {
        diag_values[row] = S->values[S->rowptr[row]];   /* diag at offset 0 */
    }
    if (partit && partit->npes > 1) {
        fesom_halo_exchange(diag_values, FESOM_HALO_NOD2D, 1, 1, partit);
    }

    /* MITgcm-style symmetric preconditioner (solver.F90:77-86):
     *   pr[diag]                 = 1 / diag(row)
     *   pr[row, col=node, n>0]   = -0.5 * (off / diag(row)) / (diag(row) + diag(node))
     */
    for (int row = 0; row < N; ++row) {
        int  start = S->rowptr[row];
        int  nend  = S->rowptr[row + 1] - start;
        real_t diag_row = S->values[start];
        S->pr_values[start] = 1.0 / diag_row;
        for (int n = 1; n < nend; ++n) {
            int node = S->colind[start + n];
            S->pr_values[start + n] = -0.5 * (S->values[start + n] / diag_row)
                                          / (diag_row + diag_values[node]);
        }
    }
    free(diag_values);
}

/*===========================================================================
 * compute_ssh_rhs_ale — linfs branch (oce_ale.F90:1821-1956)
 *===========================================================================*/

void fesom_compute_ssh_rhs_linfs(const struct fesom_mesh *mesh,
                                 struct fesom_dyn        *dyn)
{
    const int N  = mesh->myDim_nod2D;
    const int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;
    /* Fortran oce_ale.F90:1862 uses `do ed=1, myDim_edge2D` only. The edge
     * partition replicates each cross-rank edge in BOTH neighbour ranks'
     * myDim_edge2D, so each interior endpoint receives full contributions
     * from its rank alone. exchange_nod at the end of the routine refreshes
     * halo entries (which were only side-effects of the loop) with owner
     * values — same pattern as Fortran. */
    const int E  = mesh->myDim_edge2D;
    const int nl = mesh->nl;
    const real_t alpha = (real_t)FESOM_PHASE1_ALPHA;
    const real_t one_minus_alpha = 1.0 - alpha;
    (void)N;

    /* Zero ssh_rhs over full local extent. */
    memset(dyn->ssh_rhs, 0, (size_t)N_alloc * sizeof(real_t));

    /* Edge loop — accumulate fluxes into ssh_rhs (lines 1862-1919). */
    for (int ed = 0; ed < E; ++ed) {
        int n1 = mesh->edges[2*ed + 0];
        int n2 = mesh->edges[2*ed + 1];
        int el1 = mesh->edge_tri[2*ed + 0];
        int el2 = mesh->edge_tri[2*ed + 1];

        real_t c1 = 0.0;
        if (el1 >= 0) {
            real_t dx1 = mesh->edge_cross_dxdy[4*ed + 0];
            real_t dy1 = mesh->edge_cross_dxdy[4*ed + 1];
            int nzmin = mesh->ulevels[el1] - 1;
            int nzmax = mesh->nlevels[el1] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                real_t u  = dyn->uv    [FESOM_ELEMVEC(el1, nz, nl) + 0];
                real_t v  = dyn->uv    [FESOM_ELEMVEC(el1, nz, nl) + 1];
                real_t ur = dyn->uv_rhs[FESOM_ELEMVEC(el1, nz, nl) + 0];
                real_t vr = dyn->uv_rhs[FESOM_ELEMVEC(el1, nz, nl) + 1];
                real_t h  = mesh->helem[FESOM_ELEM3D(el1, nz, nl)];
                c1 += alpha * ((v + vr) * dx1 - (u + ur) * dy1) * h;
            }
        }
        real_t c2 = 0.0;
        if (el2 >= 0) {
            real_t dx2 = mesh->edge_cross_dxdy[4*ed + 2];
            real_t dy2 = mesh->edge_cross_dxdy[4*ed + 3];
            int nzmin = mesh->ulevels[el2] - 1;
            int nzmax = mesh->nlevels[el2] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                real_t u  = dyn->uv    [FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t v  = dyn->uv    [FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t ur = dyn->uv_rhs[FESOM_ELEMVEC(el2, nz, nl) + 0];
                real_t vr = dyn->uv_rhs[FESOM_ELEMVEC(el2, nz, nl) + 1];
                real_t h  = mesh->helem[FESOM_ELEM3D(el2, nz, nl)];
                c2 -= alpha * ((v + vr) * dx2 - (u + ur) * dy2) * h;
            }
        }
        dyn->ssh_rhs[n1] += (c1 + c2);
        dyn->ssh_rhs[n2] -= (c1 + c2);
    }

    /* linfs branch: ssh_rhs += (1 - alpha) * ssh_rhs_old (lines 1947-1951).
       For Phase 1 with no surface forcing yet, water_flux is zero so we skip
       the corresponding term in the non-linfs branch (which we do not exercise). */
    for (int n = 0; n < N; ++n) {
        dyn->ssh_rhs[n] += one_minus_alpha * dyn->ssh_rhs_old[n];
    }
}

/*===========================================================================
 * Solverinfo lifecycle
 *===========================================================================*/

void fesom_solverinfo_alloc(fesom_solverinfo       *si,
                            const struct fesom_mesh *mesh)
{
    memset(si, 0, sizeof(*si));
    si->maxiter = FESOM_PHASE1_MAXITER;
    si->soltol  = FESOM_PHASE1_SOLTOL;
    /* Allocate for full local extent — matrix-vector reads pp at column
     * indices that may reach into halo. */
    size_t n = (size_t)(mesh->myDim_nod2D + mesh->eDim_nod2D);
    si->rr  = calloc(n, sizeof(real_t));
    si->zz  = calloc(n, sizeof(real_t));
    si->pp  = calloc(n, sizeof(real_t));
    si->App = calloc(n, sizeof(real_t));
    FESOM_CHECK(si->rr && si->zz && si->pp && si->App,
                "solverinfo: out of memory");
}

void fesom_solverinfo_free(fesom_solverinfo *si)
{
    free(si->rr);
    free(si->zz);
    free(si->pp);
    free(si->App);
    memset(si, 0, sizeof(*si));
}

/*===========================================================================
 * CG iteration — ssh_solve_cg (solver.F90:98-281)
 *
 * Solves S * x = rhs where x = dyn->d_eta (in/out), rhs = dyn->ssh_rhs.
 *===========================================================================*/

/* CSR mat-vec: y = S * v. Same loop body as the Fortran sum() expression
   on line 160-161 / 203 / 239 (uses pr_values for the precond step). */
static void csr_matvec(const fesom_ssh_stiff *S,
                       const real_t          *Avals,
                       const real_t          *v,
                       real_t                *y, int dim)
{
    for (int row = 0; row < dim; ++row) {
        real_t s = 0.0;
        int rstart = S->rowptr[row];
        int rend   = S->rowptr[row + 1];
        for (int n = rstart; n < rend; ++n) {
            s += Avals[n] * v[S->colind[n]];
        }
        y[row] = s;
    }
}

int fesom_ssh_solve_cg(const fesom_ssh_stiff *S,
                       fesom_solverinfo      *si,
                       const struct fesom_mesh *mesh,
                       struct fesom_dyn        *dyn)
{
    /* Slice 30e — full parallel CG. Each iteration exchanges pp before SpMV,
     * rr after residual update; reductions use MPI_Allreduce. */
    /* Helper macro: exchange a nod2D field. */
    #define EXCH(field) fesom_halo_exchange((field), FESOM_HALO_NOD2D, 1, 1, si->partit)
    /* All-reduce sum of a single double in-place. */
    #define ALLREDUCE_SUM(var) MPI_Allreduce(MPI_IN_PLACE, &(var), 1, MPI_DOUBLE, MPI_SUM, si->partit->MPI_COMM_FESOM)

    const int    N      = mesh->myDim_nod2D;
    const real_t soltol = si->soltol;
    real_t      *X      = dyn->d_eta;
    const real_t *rhs   = dyn->ssh_rhs;
    real_t       *rr    = si->rr;
    real_t       *zz    = si->zz;
    real_t       *pp    = si->pp;
    real_t       *App   = si->App;

    /* Initial ‖rhs‖² and tolerance (Fortran solver.F90:142-154). */
    real_t s_old = 0.0;
    for (int row = 0; row < N; ++row) s_old += rhs[row] * rhs[row];
    if (si->partit && si->partit->npes > 1) ALLREDUCE_SUM(s_old);
    /* The global problem size for the rtol normalisation must be the GLOBAL
     * row count, not local (Fortran uses nod2D). */
    int N_global = (si->partit && si->partit->npes > 1) ? mesh->nod2D : N;
    real_t rtol = soltol * sqrt(s_old / (real_t)N_global);

    if (s_old == 0.0) {
        memset(X, 0, (size_t)N * sizeof(real_t));
        si->last_iters = 0;
        return 0;
    }

    /* r0 = rhs - A * X. Need X halo exchanged before SpMV. */
    if (si->partit && si->partit->npes > 1) EXCH(X);
    csr_matvec(S, S->values, X, rr, N);
    for (int row = 0; row < N; ++row) rr[row] = rhs[row] - rr[row];
    if (si->partit && si->partit->npes > 1) EXCH(rr);

    /* z0 = M^{-1} r0; pp = z0 */
    csr_matvec(S, S->pr_values, rr, zz, N);
    memcpy(pp, zz, (size_t)N * sizeof(real_t));

    /* s_old = r0·z0 */
    s_old = 0.0;
    for (int row = 0; row < N; ++row) s_old += rr[row] * zz[row];
    if (si->partit && si->partit->npes > 1) ALLREDUCE_SUM(s_old);

    int iter = 0;
    int   verbose = (getenv("FESOM_VERBOSE_CG") != NULL);
    /* Heartbeat from rank 0 every 100 iters regardless of env var, so a
     * hung CG always shows something in the log. */
    int   heartbeat_every = 100;
    for (iter = 1; iter <= si->maxiter; ++iter) {
        /* App = A * pp; need pp halo exchanged. */
        if (si->partit && si->partit->npes > 1) EXCH(pp);
        csr_matvec(S, S->values, pp, App, N);

        /* α = s_old / (pp·App) */
        real_t s_aux = 0.0;
        for (int row = 0; row < N; ++row) s_aux += pp[row] * App[row];
        if (si->partit && si->partit->npes > 1) ALLREDUCE_SUM(s_aux);
        if (s_aux == 0.0 || s_aux != s_aux) {       /* NaN/zero check */
            if (si->partit == NULL || si->partit->mype == 0) {
                fprintf(stderr, "[fesom_ssh] CG abort at iter %d: pp·App = %g (s_old=%g)\n",
                        iter, (double)s_aux, (double)s_old); fflush(stderr);
            }
            FESOM_DIE("CG: pp·App is %g — matrix singular or NaN propagated", s_aux);
        }
        real_t al = s_old / s_aux;

        for (int row = 0; row < N; ++row) {
            X [row] += al * pp [row];
            rr[row] -= al * App[row];
        }
        if (si->partit && si->partit->npes > 1) EXCH(rr);

        csr_matvec(S, S->pr_values, rr, zz, N);

        real_t sp0 = 0.0, sp1 = 0.0;
        for (int row = 0; row < N; ++row) {
            sp0 += rr[row] * zz[row];
            sp1 += rr[row] * rr[row];
        }
        if (si->partit && si->partit->npes > 1) {
            ALLREDUCE_SUM(sp0);
            ALLREDUCE_SUM(sp1);
        }

        real_t residual = sqrt(sp1 / (real_t)N_global);
        if ((si->partit == NULL || si->partit->mype == 0)
            && (verbose ? (iter <= 5 || iter % 50 == 0)
                        : (iter % heartbeat_every == 0))) {
            fprintf(stderr, "[fesom_ssh] CG iter %4d: res=%.4e rtol=%.4e\n",
                    iter, (double)residual, (double)rtol);
            fflush(stderr);
        }
        if (residual < rtol) break;
        if (residual != residual || residual > 1e30) {       /* NaN/divergence */
            if (si->partit == NULL || si->partit->mype == 0) {
                fprintf(stderr,
                    "[fesom_ssh] CG abort at iter %d: residual=%g (NaN or divergence)\n",
                    iter, (double)residual); fflush(stderr);
            }
            FESOM_DIE("CG residual diverged");
        }

        real_t be = sp0 / s_old;
        s_old = sp0;
        for (int row = 0; row < N; ++row) pp[row] = zz[row] + be * pp[row];
    }
    if (iter > si->maxiter) {
        if (si->partit == NULL || si->partit->mype == 0) {
            fprintf(stderr, "[fesom_ssh] CG hit maxiter=%d without converging "
                    "(last residual ~%.4e, rtol=%.4e)\n",
                    si->maxiter, (double)sqrt(s_old/(real_t)N_global), (double)rtol);
            fflush(stderr);
        }
        FESOM_DIE("CG did not converge");
    }
    if (si->partit && si->partit->npes > 1) EXCH(X);
    si->last_iters = iter;
    #undef EXCH
    #undef ALLREDUCE_SUM
    return iter;
}
