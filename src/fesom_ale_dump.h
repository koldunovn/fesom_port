#ifndef FESOM_ALE_DUMP_H
#define FESOM_ALE_DUMP_H

/*
 * ALE dual-instrumentation dump (Phase Z0 of the zstar port).
 *
 * Mirror of the Fortran `ale_dump_mod` (uncommitted in the Intel tree,
 * oce_ale.F90): same filename/format as the KPP dumps —
 *   <dir>/ale_dump_s<step>_<tag>_rank<R>.txt
 *   header "# step= tag= rank= N= ncomp="; gid-keyed rows, OWNED entities only
 * so scripts/ale_dump_diff.py reads both codes. Gated on FESOM_ALE_DUMP_DIR;
 * dumps ocean steps 1..FESOM_ALE_DUMP_STEPS (default 3). Every entry point is
 * a no-op when the env var is unset.
 *
 * Call sites in fesom_timestep mirror the Fortran oce_timestep_ale sites:
 *   forcing   — step entry (water_flux, virtual_salt, relax_salt, real_salt_flux)
 *   pgf       — after the pressure force (pgf_x, pgf_y; nl-1 comps, elements)
 *   sshsolve  — after the CG solve (ssh_rhs, d_eta)
 *   hbar      — after the eta_n update (hbar, hbar_old, ssh_rhs_old, eta_n + dhe)
 *   vertvel   — after vert_vel (Wvel nl comps, hnode_new nl-1 comps) — BEFORE
 *               the GM bolus add, which modifies dyn->w in place
 *   thickness — after the commit (hnode, zbar_3d_n, Z_3d_n, helem)
 */

struct fesom_mesh;
struct fesom_dyn;
struct fesom_aux;
struct fesom_forcing;
struct fesom_partit;

/* 1 if FESOM_ALE_DUMP_DIR is set and step <= FESOM_ALE_DUMP_STEPS. */
int fesom_ale_dump_active(int step);

void fesom_ale_dump_forcing(int step, const struct fesom_forcing *forcing,
                            const struct fesom_mesh *mesh,
                            struct fesom_partit *partit);
void fesom_ale_dump_pgf(int step, const struct fesom_aux *aux,
                        const struct fesom_mesh *mesh,
                        struct fesom_partit *partit);
void fesom_ale_dump_sshsolve(int step, const struct fesom_dyn *dyn,
                             const struct fesom_mesh *mesh,
                             struct fesom_partit *partit);
void fesom_ale_dump_hbar(int step, const struct fesom_dyn *dyn,
                         const struct fesom_mesh *mesh,
                         struct fesom_partit *partit);
void fesom_ale_dump_vertvel(int step, const struct fesom_dyn *dyn,
                            const struct fesom_mesh *mesh,
                            struct fesom_partit *partit);
void fesom_ale_dump_thickness(int step, const struct fesom_mesh *mesh,
                              struct fesom_partit *partit);

#endif /* FESOM_ALE_DUMP_H */
