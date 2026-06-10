#ifndef FESOM_ALE_H
#define FESOM_ALE_H

#include "fesom_types.h"

struct fesom_mesh;
struct fesom_dyn;

/*
 * which_ALE runtime switch (Phase Z0 of the zstar port).
 *
 * Mirrors the Fortran module globals: g_config which_ALE,
 * gen_modules_forcing use_virt_salt and o_PARAM is_nonlinfs — the latter
 * two are DERIVED from which_ALE in oce_setup_step.F90:115-125:
 *   linfs  → use_virt_salt = .true.,  is_nonlinfs = 0.0
 *   zstar  → use_virt_salt = .false., is_nonlinfs = 1.0
 *
 * Selected by FESOM_ALE=linfs|zstar (default linfs — byte-identical to v1.0).
 * Any other value (zlevel etc.) aborts: not ported (no reference run).
 * fesom_ale_mode_init must be called once at startup before any consumer.
 */
extern int    fesom_ale_zstar;      /* 0 = linfs (default), 1 = zstar      */
extern int    fesom_use_virt_salt;  /* 1 under linfs, 0 under zstar        */
extern real_t fesom_is_nonlinfs;    /* 0.0 under linfs, 1.0 under zstar    */
void fesom_ale_mode_init(void);

/*
 * ALE (Arbitrary Lagrangian-Eulerian) step. Phase 1 implements only the
 * `which_ALE = 'linfs'` branch:
 *
 *   - hnode_new = hnode (no thickness motion)
 *   - w computed by accumulating horizontal edge fluxes over the two
 *     adjacent cells (Gauss-law div integration), then cumsum vertically
 *     from bottom up, then division by area(nz, n) to get m/s.
 *
 * Mirror of the linfs path of vert_vel_ale (oce_ale.F90:2132+) plus
 * update_thickness_ale (oce_ale.F90:1035+).
 *
 * Phase G6a: when `gm_on != 0`, the same edge loop also accumulates
 *   dyn->fer_w from dyn->fer_uv. The dyn->w computation path is
 *   bit-identical to the gm_on=0 case, by construction (the fer_w
 *   block adds reads/writes only on dyn->fer_w / dyn->fer_uv — never
 *   touches dyn->w or its inputs). This is the bit-identity off-switch
 *   guarantee tested at G9.
 *
 * Future: zlevel surface-flux correction, zstar layer scaling,
 * partial cells.
 */
void fesom_ale_vert_vel_linfs(const struct fesom_mesh *mesh,
                              struct fesom_dyn        *dyn,
                              int                      gm_on);

/* hnode_new := hnode (linfs invariant). Called at the start of the ALE step. */
void fesom_ale_thickness_linfs(struct fesom_mesh *mesh);

/* Commit at end of timestep: hnode := hnode_new and helem := mean of 3 vertex hnode. */
void fesom_ale_commit_thickness(struct fesom_mesh *mesh);

struct fesom_partit;
struct fesom_dyn;

/*
 * Phase Z1 — zstar thickness machinery.
 *
 * fesom_ale_init_thickness_zstar: literal port of init_thickness_ale
 * (oce_ale.F90:962-1221), zstar case, including the mode-shared pre-block
 * (ssh_rhs_old = (hbar-hbar_old)·areasvol/dt; eta_n = α·hbar_old+(1-α)·hbar —
 * note the REVERSED weights vs the per-step update; all zeros at cold start).
 * Node loop (myDim+eDim): hnode = (zbar(nz)-zbar(nz+1))·(1+hbar/dd) for
 * nz=1..nlevels_nod2D_min-2 with dd = zbar(1)-zbar(nlevels_nod2D_min-1);
 * nominal spacing for the bottom-intersecting levels; bottom layer =
 * bottom_node_thickness. Element loop (myDim): dhe = mean(hbar), helem =
 * mean of vertex hnode, bottom = bottom_elem_thickness. Then
 * hnode_new := hnode and exchange_elem(helem). Replaces fesom_ic_thickness
 * when FESOM_ALE=zstar (cold start hbar=0 → values equal the linfs init by
 * construction; helem may differ from the linfs zbar-diff form by ≤1 ulp
 * because the mean-of-3 formula rounds differently — exactly as in Fortran).
 *
 * fesom_ale_update_thickness_zstar: literal port of update_thickness_ale
 * (oce_ale.F90:1226-1444), zstar case. Per-step COMMIT over myDim+eDim,
 * bottom→top for nz = nlevels_nod2D_min-2 .. 1 (1-based):
 *   hnode = hnode_new; zbar_3d_n(nz) = zbar_3d_n(nz+1) + hnode_new(nz);
 *   Z_3d_n(nz) = zbar_3d_n(nz+1) + hnode_new(nz)/2
 * (bottom-protected levels keep their init geometry); helem = mean of vertex
 * hnode for nz = 1..nlevels(elem)-2 (bottom helem stays bottom_elem_thickness);
 * exchange_elem(helem). ldiag_DVD rescue NOT ported (diagnostic off).
 * NOTE: until Z5 produces hnode_new, this commit is a no-op by construction.
 */
void fesom_ale_init_thickness_zstar(struct fesom_mesh   *mesh,
                                    struct fesom_dyn    *dyn,
                                    struct fesom_partit *partit);

void fesom_ale_update_thickness_zstar(struct fesom_mesh   *mesh,
                                      struct fesom_partit *partit);

/*
 * Phase Z5 — zstar branch of vert_vel_ale (oce_ale.F90:2755-2827): on top of
 * the shared divergence Wvel (fesom_ale_vert_vel_linfs), distribute
 * (hbar − hbar_old) over the stretched column: the vertically-integrated Wvel
 * correction −(zbar_3d_n(nz)−dd1)·(dd/dt) and hnode_new = hnode + Δz·dd for
 * nz = 1..nlevels_nod2D_min−2 (owned nodes), plus the surface water-flux BC
 * Wvel(1) −= water_flux. Caller must exchange Wvel AND hnode_new after
 * (oce_ale.F90:2870-2871). water_flux NULL → zero flux.
 */
void fesom_ale_vert_vel_zstar_distribute(const struct fesom_mesh *mesh,
                                         struct fesom_dyn        *dyn,
                                         const real_t            *water_flux);

/*
 * Vertical CFL (compute_CFLz, oce_ale.F90:2935-3022):
 *   CFL_z[nz, n] gets contributions from W at the two interfaces above and
 *   below layer nz, scaled by the appropriate layer thickness.
 *
 * wsplit (compute_Wvel_split, oce_ale.F90:3026-3074):
 *   For each interface where CFL_z > wsplit_maxcfl, split W into explicit
 *   (w_e, used by tracer advection) and implicit (w_i, used by impl_vert_visc)
 *   parts so the tracer scheme stays CFL-bounded. Default behaviour
 *   (CFL ≤ maxcfl) is w_e = w, w_i = 0.
 */
void fesom_ale_compute_cflz(const struct fesom_mesh *mesh,
                            struct fesom_dyn        *dyn);

void fesom_ale_compute_wvel_split(const struct fesom_mesh *mesh,
                                  struct fesom_dyn        *dyn);

#endif /* FESOM_ALE_H */
