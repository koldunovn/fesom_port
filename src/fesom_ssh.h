#ifndef FESOM_SSH_H
#define FESOM_SSH_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;

/*
 * SSH stiffness matrix in CSR. Mirrors mesh%ssh_stiff (a t_sparse_matrix).
 *
 * Layout (literal port of init_stiff_mat_ale, oce_ale.F90:1393-1684):
 *   - dim    = nod2D (one row per local-owned node)
 *   - rowptr [dim + 1]  — 0-based start index of each row in colind/values
 *   - colind [nnz]      — local 0-based node index for each non-zero
 *   - values [nnz]      — matrix coefficient
 *   - pr_values [nnz]   — preconditioner coefficients (same sparsity pattern)
 *
 * The diagonal entry of each row is at offset 0:
 *   diag(row) = values[rowptr[row]],
 * because n_pos(1, n) = n is set first when building neighbours
 * (oce_ale.F90:1444-1445).
 *
 * Phase 1 (linfs) builds the matrix once at startup; it is not updated each
 * step (oce_ale.F90:3722 gates update by `.not. trim(which_ale)=='linfs'`).
 */
typedef struct fesom_ssh_stiff {
    int     dim;
    int     nnz;
    int    *rowptr;     /* [dim + 1] */
    int    *colind;     /* [nnz]     */
    real_t *values;     /* [nnz]     */
    real_t *pr_values;  /* [nnz]     — preconditioner */
} fesom_ssh_stiff;

/*
 * CG solver workspace. Mirror of t_solverinfo (MOD_DYN.F90:13-25).
 */
struct fesom_partit;

typedef struct fesom_solverinfo {
    int     maxiter;
    real_t  soltol;
    real_t *rr;
    real_t *zz;
    real_t *pp;
    real_t *App;
    int     last_iters;     /* iters used by the last solve (diagnostic) */
    struct fesom_partit *partit;  /* set after fesom_solverinfo_alloc() — used by parallel CG */
} fesom_solverinfo;

/*--- Lifecycle -------------------------------------------------------------*/
void fesom_ssh_stiff_alloc_and_build(fesom_ssh_stiff       *S,
                                     const struct fesom_mesh *mesh);
void fesom_ssh_stiff_free(fesom_ssh_stiff *S);

void fesom_solverinfo_alloc(fesom_solverinfo       *si,
                            const struct fesom_mesh *mesh);
void fesom_solverinfo_free(fesom_solverinfo *si);

/*
 * Compute the diagonal-Jacobi-with-symmetrisation preconditioner from the
 * current stiffness values. Mirror of ssh_solve_preconditioner
 * (solver.F90:31-95).
 */
void fesom_ssh_preconditioner(fesom_ssh_stiff *S, const struct fesom_mesh *mesh,
                              struct fesom_partit *partit);

/*
 * Build the SSH RHS into dyn->ssh_rhs (linfs branch). Mirror of
 * compute_ssh_rhs_ale (oce_ale.F90:1821-1956).
 *
 * For Phase 1: ssh_rhs(n) = -alpha * div_int((U+U_rhs)*helem) + (1-alpha)*ssh_rhs_old
 * (water_flux is zero in Phase 1 — surface forcing comes later.)
 */
void fesom_compute_ssh_rhs_linfs(const struct fesom_mesh *mesh,
                                 struct fesom_dyn        *dyn);

/*
 * Conjugate-gradient SSH solve. Mirror of ssh_solve_cg (solver.F90:98-281).
 * Result lands in dyn->d_eta. Returns the iteration count.
 */
int fesom_ssh_solve_cg(const fesom_ssh_stiff *S,
                       fesom_solverinfo      *si,
                       const struct fesom_mesh *mesh,
                       struct fesom_dyn        *dyn);

#endif /* FESOM_SSH_H */
