# Port the mEVP sea-ice rheology solver (whichEVP=1)

> From the 2026-06-11 brainstorm. All design decisions are user-approved (2026-06-11):
> env knob `FESOM_WHICH_EVP=1` (numeric, mirrors the namelist integer); run matrix =
> **single-knob Fortran reference + 5yr C stability** (no combined zstar/TKE/mEVP leg);
> new file pair `fesom_ice_maevp.{c,h}`; dev on branch `mevp` off main@2714071 in the
> `fesom2_port_zstar` worktree. Default stays whichEVP=0 with off-switch byte-identity.
>
> **Revised 2026-06-11 after plan-review (2 agents). Folded in:**
> (1) **`bc_index_nod2D` already EXISTS in the C** (`fesom_mesh.h:141`, `real_t`, built
> unconditionally at `fesom_ice.c:199-209`) — but with the **multi-rank-UNSAFE
> `edge_tri[2nd]<0` boundary test** that `fesom_mesh.h:35-40` explicitly warns
> misclassifies interior edges whose halo element is missing. It has ZERO readers today
> (grep-verified); mEVP is the FIRST consumer → M1 REBUILDS it with the
> `myList_edge2D > edge2D_in` convention (byte-neutral: nothing reads it on the default
> path). Without this, mEVP would hard-zero velocity at inter-rank seam nodes — a
> multi-rank-only divergence.
> (2) Three new fidelity traps: non-ice nodes SKIPPED not zeroed (trap 13); drag carries
> rdt (trap 1 amended); uice_old never touched (trap 14).
> (3) Byte-gate baseline made executable: M0 captures a pinned 2714071 baseline run
> (16r/200 steps/dt=500, snapshots + run.log — the zstar-Z0/TKE-T0 config) BEFORE any
> src change; all later byte-gates compare against that capture.
> (4) Branch state VERIFIED with git: `main` == `zstar` == 2714071 (the TKE ff-merge
> brought zstar+TKE into main; `git merge-base --is-ancestor 8260dea main` passes) —
> creating `mevp` here is a same-commit branch, no working-tree change. A reviewer
> concern to the contrary was based on a stale memory, now corrected.
> (5) Env-knob parse cite fixed (numeric `atoi` at init, `fesom_ice.c:392-393` precedent —
> not the per-step presence pattern at :297).
> (6) Ice-velocity VECTOR gate added to M5 (`feedback_magnitude_invariant_masks_rotation`:
> trap 3 is a rotation-class difference; speed/extent metrics are rotation-blind).

## Overview

Port the modified-EVP sea-ice rheology (Bouillon et al., Ocean Modelling 2013; the
pseudo-time iteration with implicit α/β stabilization that replaces standard EVP
subcycling) from FESOM2 Fortran to C as a **strictly faithful line-by-line port**,
runtime-selected by `FESOM_WHICH_EVP=1` (unset/`0` = standard EVP, the default; `2`/garbage
aborts at init with "not ported").

The ENTIRE mEVP path is **one Fortran subroutine**: `EVPdynamics_m`
(`ice_maEVP.F90:429-882`). The NR optimization pass inlined `ssh2rhs`, `stress_tensor_m`,
and `stress2rhs_m` into it; the standalone copies in the same file are called only by aEVP
(whichEVP=2, out of scope) and `stress_tensor_m` has no call sites at all. So the port is
one new C function plus wiring: two aux velocity arrays, two scalar params, one mesh BC
array, and a dispatcher branch.

**GUIDING PRINCIPLE (non-negotiable):** THIS IS A PORT. Reproduce every formula, constant,
loop bound, mask, branch asymmetry, and exchange placement exactly. The fidelity-traps
checklist below names the places where mEVP *differs* from the already-ported standard EVP —
the highest-risk lines are the ones where copying the std-EVP template would be wrong.

**Standard EVP stays byte-identical to main@2714071** (the zstar+TKE-validated HEAD) —
verified by a byte-gate in every phase that touches shared code.

## Reference mEVP configuration (VERIFIED 2026-06-11 against work_core/namelist.ice)

Fortran reference = `work_core` config with ONLY `whichEVP` flipped 0→1. Everything mEVP
reads is already in the reference `namelist.ice` (no namelist-over-default trap: namelist
and module default are both 250 for α/β):

| Symbol | Value | Source | Role in EVPdynamics_m |
|---|---|---|---|
| `whichEVP` | **1** (the flip) | namelist.ice | dispatch `ice_setup_step.F90:216-218` |
| `alpha_evp` | `250` | namelist.ice == MOD_ICE.F90:198 default | `det2=1/(1+α)`, `det1=α·det2` |
| `beta_evp` | `250` | namelist.ice == module default | `+β·u_aux` in rhs; `(1+β+drag)` in det |
| `evp_rheol_steps` | `120` | namelist.ice | **mEVP iteration count** (`steps`) |
| `Pstar` | `30000.0` | namelist.ice | `pressure_fac` |
| `ellipse` | `2.0` | namelist.ice | `vale = 1/e² = 0.25` |
| `c_pressure` | `20.0` | namelist.ice | `exp(-c·(1-asum))` |
| `delta_min` | `1.0e-11` | namelist.ice | ADDITIVE: `pressure_fac/(delta+delta_min)` |
| `Cd_oce_ice` | `0.0055` | namelist.ice | drag |
| `theta_io` | `0.0` | namelist.ice | **NOT read by mEVP** (no rotation — unlike std EVP) |
| `ice_ave_steps` | `1` | namelist.ice | `ice_dt = dt = 1800`; `rdt = ice_dt` FULL step |
| `use_floatice` | `.false.` | g_config default | levitating ssh2rhs branch (covers linfs AND zstar) |
| `use_cavity` | `.false.` | config | cavity edge-BC block skipped (std-EVP precedent) |
| `density_0`, `g` | C constants | o_param | drag, elevation gradient |

⚠️ A stale `work_core_evp1/` exists in `fesom27/fesom2/` (January experiment: its log shows
`EVP scheme option= 0`, namelist.ice identical to work_core, OLD forcing settings
nm_nc_iyear=1948/tmid=1). **Do NOT reuse it** — create a fresh `work_mevp` from today's
`work_core`.

## Context (from discovery, 2026-06-11)

### Fortran source (port reference) — `port2/fesom2/src/`
| File / lines | Item | Status |
|---|---|---|
| `ice_maEVP.F90:429-882` | `EVPdynamics_m` — THE port (everything inlined) | port |
| `ice_maEVP.F90:86-193` | `stress_tensor_m` — **zero call sites** (dead) | not ported |
| `ice_maEVP.F90:200-319, 325-422` | `ssh2rhs`, `stress2rhs_m` — called only by `EVPdynamics_a` | not ported (aEVP) |
| `ice_maEVP.F90:891-1276` | aEVP family (`find_alpha/beta_field_a`, `stress_tensor_a`, `EVPdynamics_a`) | OUT OF SCOPE; knob=2 aborts |
| `MOD_ICE.F90:743-748` | `uice_aux`/`vice_aux` alloc (node_size, zero-init) when whichEVP≠0 | mirror |
| `MOD_ICE.F90:198` | `alpha_evp=250, beta_evp=250` defaults | mirror as C constants |
| `MOD_ICE.F90:889-895` | `bc_index_nod2D` build (in ice_init, unconditional) | mirror |
| `ice_setup_step.F90:209-226` | whichEVP dispatch SELECT CASE | mirror at fesom_ice.c:327 |
| `io_meandata.F90:1593-1594` | `whichEVP==1` IO block — **EMPTY** | nothing to port |
| `MOD_ICE.F90:455-459` | restart write of aux fields | N/A — C is cold-start only |
| `fesom_module.F90:323` | log echo `'EVP scheme option='` | M0 reference-run gate |
| `ice_EVP.F90:20-84` | `module fesom_evp_dump` (env-gated dump helpers, already in the working copy) | reuse from ice_maEVP.F90 |

### C integration points — `fesom2_port_zstar/src/`
| C file | Current state | Change |
|---|---|---|
| `fesom_ice_maevp.{c,h}` | NEW | `fesom_ice_evp_dynamics_m()` — literal port + dump points |
| `fesom_ice_types.h` | whichEVP locked 0, aux arrays "NOT mirrored" comment | + `alpha_evp`, `beta_evp`; + `uice_aux`/`vice_aux` pointers; + mEVP work arrays (see Technical Details); update header doc |
| `fesom_ice.c` | hard `whichEVP=0` at :68; numeric-env precedent `atoi(getenv("FESOM_DIAG_GID"))` at :392-393; dispatch abort :327; `bc_index_nod2D` build at :199-209 (**edge_tri-based — multi-rank-unsafe, must be rebuilt**) | parse `FESOM_WHICH_EVP` (atoi) in `fesom_ice_init`, replacing the hard 0 at :68; conditional alloc/free of aux + work arrays; REBUILD `bc_index_nod2D` with the `edge2D_in` convention; dispatch `==1` branch |
| `fesom_mesh.h` | has `edges` [2·edge2D], `edge2D_in`, `coriolis_node`, `gradient_sca` [6·elem], `metric_factor`, `area` [n·nl+0]; **`bc_index_nod2D` already declared at :141 (`real_t*`, [myDim+eDim], freed by fesom_mesh.c:70)**; :35-40 documents why `edge_tri[2nd]<0` misclassifies at multi-rank | NO new member needed — keep `real_t` (Fortran declares INTEGER but assigns `1._WP/0._WP` and uses it multiplicatively); fix the build, not the type |
| `fesom_ice_evp.c` | std EVP + `evp_dump` infra | UNTOUCHED (byte-gate); dump infra pattern reused in the new file |
| `fesom_halo.h` | `fesom_exchange_nod2D(field, partit)` | used ×2 per substep (u_aux, v_aux) — the Fortran `exchange_nod_begin/end(u,v)` pair collapses to two blocking calls |
| `CMakeLists.txt` | — | + `fesom_ice_maevp.c` |
| `jobs/` | zstar/TKE templates (`job_tke_t0_byteident`, `job_zstar_z7_crossrank`, `validation_2yr.sh`, `production_864_5yr.sh`) | + `job_mevp_*` clones with `FESOM_WHICH_EVP=1`; **unique OUT_DIR per job** |
| Fortran working copy (uncommitted, like the KPP/TKE gates) | `fesom_evp_dump` module exists | + dump calls inside `EVPdynamics_m` |

### Existing C arrays mEVP reuses (no new allocation)
`sigma11/12/22`, `eps11/12/22` (element, shared with std EVP — only one scheme runs per
execution), `uice_rhs`/`vice_rhs`, `srfoce_u/v/ssh`, `stress_atmice_x/y`,
`data[AICE/MICE].values_rhs` (**reused as ssh-gradient scratch** `rhs_a`/`rhs_m` — Fortran
does exactly this; safe because ice-advection `TG_rhs` overwrites them after EVP each step).

## Out of scope (gate-only / not ported)
- **aEVP (whichEVP=2)** and the four `_a` routines + `alpha_evp_array`/`beta_evp_array` —
  knob value 2 aborts at init.
- **`stress_tensor_m`/`ssh2rhs`/`stress2rhs_m` standalones** — dead or aEVP-only (inlined
  copies are what we port).
- **`use_floatice` branch** (ice_maEVP.F90:547-587) — levitating else-branch only
  (`:589-620`, covers linfs AND zstar per the Fortran comment); same exclusion as the
  std-EVP port.
- **Cavity edge-BC block** (`:842-856`, `if (use_cavity)`) — skipped like std EVP; the
  `ulevels>1`/`ulevels_nod2D>1` cycles ARE kept as written (no-op in our config).
- **ICEPACK blocks** (`#if defined (__icepack)`) — not compiled.
- **Restart I/O** of aux fields — no restart system in C.
- **IO additions** — the Fortran whichEVP==1 output block is empty; standard ice outputs
  (a_ice, m_ice, uice, vice) are unchanged.

## Fidelity-traps checklist (review EVERY item — all 14 — against the diff in M2)

The places where mEVP differs from the std-EVP template — copying `fesom_ice_evp.c`
idioms here would be a silent wrong-physics port:

1. **`rdt = ice_dt` — the FULL ice step** (`:517`). Std EVP divides by `evp_rheol_steps`.
   The #1 copy-from-template trap. AND the drag term **carries rdt structurally**:
   `drag = rdt·Cd·umod·density_0·inv_thickness` (`:807`) — std EVP's drag has NO rdt
   (ice_EVP.F90:812 / fesom_ice_evp.c:424). `drag·u_w` sits OUTSIDE the `rdt·(…)` group
   in rhsu (`:810`). Fixing the rdt value but copying the std drag line is a ×1800 error.
2. **0.5 factor placement**: `pressure_fac = det2·pstar·msum·exp(-c_pressure(1-asum))` has
   **NO 0.5** (`:664`; std EVP's ice_strength HAS one — `feedback_ice_strength_halffactor`).
   The 0.5 lives inside the sigma11/22 updates (`:722-723`); the sigma12 update (`:721`)
   has NO 0.5. Port the factor *placement*, not the std-EVP shape.
3. **NO theta_io rotation** anywhere in mEVP — drag is `drag·u_w` directly (`:810-811`);
   std EVP rotates with `ax/ay`.
4. **Masks**: element `ice_el = (msum > 0.01)` with msum = 3-vertex MEAN of m_ice (`:660-661`;
   std EVP uses per-vertex `m_ice>0 && a_ice>0`); node `ice_nod = (a_ice >= 0.01)` (`:633`).
5. **Mass regularization verbatim**: `mass = M/((1+M²)·area)` (`:637-638`) — do NOT
   "simplify" to `1/(M·area)`; the form is the small-mass regularizer.
6. **Asymmetric assembly guards**: ssh-gradient assembly writes ALL 3 element nodes
   UNGUARDED (`:604-617`, halo entries written-never-read; arrays are node_size so safe);
   stress-divergence assembly is GUARDED `elnodes ≤ myDim_nod2D` (`:740-790`). Preserve
   the asymmetry.
7. **`rhs_a`/`rhs_m` zeroed over myDim ONLY** (`:539-543`); `/area` scaling happens INSIDE
   the `ice_nod` branch only (`:641-642`). Halo garbage in rhs is faithful.
8. **`u_ice_aux = u_ice` init (`:521-522`) AND the final copy back (`:877-880`) both cover
   myDim+eDim** (`feedback_write_loops_halo`). No exchange after the final copy — the aux
   halo is current from the last substep's exchange.
9. **`bc_index_nod2D` multiplies `det`** (`:814`) — zeroes coastal velocity; redundant with
   the edge-BC loop but port BOTH. The build loop runs **owned edges only**
   (`do n=1, myDim_edge2D`, MOD_ICE.F90:891) — some eDim boundary nodes stay =1, which is
   harmless (mEVP reads bc_index at owned nodes only, `:795/:814`). Do NOT extend the
   build loop to halo edges. ⚠️ The EXISTING C build (fesom_ice.c:199-209) uses the
   `edge_tri[2nd]<0` test — at multi-rank that marks inter-rank seam nodes as boundary
   (`bc_index=0` → mEVP velocity hard-zeroed where Fortran has 1). REBUILD with
   `myList_edge2D[ed] > edge2D_in` (the std-EVP BC-loop convention,
   fesom_ice_evp.c:453-455, which documents exactly this misclassification).
10. **`delta_min` is ADDITIVE**: `pressure = pressure_fac/(delta + delta_min)` (`:714`).
11. **sigma NOT zeroed on entry** — `sigma11/12/22` persist across time steps (and across
    no-ice intervals for masked elements). Initial values are the alloc-time zeros.
12. **-r8 literals** (`feedback_r8_literals_double`): `0.01`, `9.0`, `1/3`, `0.5` are all
    double in C; never float-suffix. `area(1,i)` = `mesh->area[n*nl + 0]`;
    `gradient_sca(1:3/4:6, el)` = `&mesh->gradient_sca[6*el]` dx=[0..2], dy=[3..5];
    `meancos = metric_factor[el]/3`.
13. **Non-ice nodes are SKIPPED, not zeroed** (plan-review catch): the node solve is
    `if (ice_nod(i)) … end if` with **NO else** (`:800-818`) — a non-ice interior node
    keeps `u_ice_aux = u_ice` (its entry value) through to the final copy, so its
    velocity is RETAINED across steps. Std EVP does the opposite (`else: U_ice=0`,
    ice_EVP.F90:824-826 — the C template has it literally at fesom_ice_evp.c:436-439).
    Do NOT port the else; it would zero ice-edge velocities every step.
14. **mEVP never touches `uice_old`/`vice_old`** — std EVP saves them each substep
    (ice_EVP.F90:790-793). No pointer, no write in EVPdynamics_m. Do not insert the save.

## Algorithm blocks (C function structure ↔ Fortran lines)

1. **Setup** (`:513-522`): `val3=1/3`, `vale=1/ellipse²`, `det2=1/(1+α)`, `det1=α·det2`,
   `rdt=ice_dt`, `steps=evp_rheol_steps`; `u_ice_aux=u_ice`, `v_ice_aux=v_ice` over
   myDim+eDim.
2. **Inlined ssh2rhs, levitating branch** (`:538-543, 590-621`): zero rhs_a/rhs_m (myDim);
   per element (skip `ulevels>1`): `bb=g·vol/3`, `aa=bb·Σdx·η`, `bb=bb·Σdy·η`; subtract
   into rhs_a/rhs_m at all 3 elnodes (unguarded).
3. **Node precompute** (`:624-647`, myDim, skip cavity): zero inv_thickness/mass +
   `ice_nod=false`; if `a_ice≥0.01`: `inv_thickness=1/max((ρᵢm+ρₛmₛ)/a, 9.0)`,
   `mass=M/((1+M²)·area)`, `rhs_a/=area`, `rhs_m/=area`, `ice_nod=true`.
4. **Element precompute** (`:650-673`): `pressure_fac=0`, `ice_el=false`; if `ulevels==1`
   and `msum>0.01`: `ice_el=true`, `pressure_fac=det2·pstar·msum·exp(-c_pressure(1-asum))`.
   Zero `u_rhs/v_rhs` (myDim).
5. **Main loop** `shortstep=1..steps` (`:680-875`):
   - per element (skip cavity; if `ice_el`): `meancos=metric_factor/3`;
     `eps11=Σdx·u_aux − Σv_aux·meancos`; `eps22=Σdy·v_aux`;
     `eps12=0.5·(Σ(dy·u_aux+dx·v_aux) + Σu_aux·meancos)`;
     `eps1=eps11+eps22`, `eps2=eps11−eps22`;
     `delta=√(eps1²+vale·(eps2²+4·eps12²))`; `pressure=pressure_fac/(delta+delta_min)`;
     `sigma12=det1·sigma12 + pressure·eps12·vale`;
     `sigma11=det1·sigma11 + 0.5·pressure·(eps1−delta+eps2·vale)`;
     `sigma22=det1·sigma22 + 0.5·pressure·(eps1−delta−eps2·vale)`;
     then stress-divergence assembly per vertex k with `elnodes(k) ≤ myDim` guard:
     `u_rhs −= elem_area·(sigma11·dx(k)+sigma12·dy(k)+sigma12·meancos)`;
     `v_rhs −= elem_area·(sigma12·dx(k)+sigma22·dy(k)−sigma11·meancos)`.
   - per node (skip cavity; if `ice_nod`): `u_rhs=u_rhs·mass+rhs_a`, `v_rhs=v_rhs·mass+rhs_m`;
     `umod=|u_aux−u_w|`; `drag=rdt·Cd·umod·density_0·inv_thickness`;
     `rhsu=u_ice+drag·u_w+rdt·(inv_thickness·τx+u_rhs)+β·u_aux` (and v analog);
     `det=bc_index/((1+β+drag)²+(rdt·f)²)`;
     `u_aux=det·((1+β+drag)·rhsu+rdt·f·rhsv)`, `v_aux=det·((1+β+drag)·rhsv−rdt·f·rhsu)`.
   - edge BC (`:824-838`): boundary edges (`myList_edge2D > edge2D_in`) → zero u_aux/v_aux
     at both endpoints. (Cavity block skipped.)
   - exchange u_aux, v_aux (two `fesom_exchange_nod2D` calls); zero u_rhs/v_rhs (myDim).
6. **Final** (`:876-880`): `u_ice=u_ice_aux`, `v_ice=v_ice_aux` over myDim+eDim.

## Development Approach
- One new C function, literal port, single phase of code (M2) — small enough that
  splitting the routine across phases would hurt reviewability.
- **Validation = dual instrumentation vs Fortran** (the KPP/TKE playbook): mirrored
  env-gated dump points in both codes, substep-level diff, controlled replay
  (`feedback_controlled_replay_validation`) in reserve.
- **Off-switch byte-gate in every phase**: default run (knob unset), pinned config
  16r/200 steps/dt=500, byte-identical (snapshots + run.log) to the M0
  `baseline_2714071` capture.
- **Size every new array to its reader's loop bound**
  (`feedback_array_size_vs_reader_loop`): aux arrays node_size (read at halo by the
  final-copy/exchange path); `pressure_fac`/`ice_el` sized myDim_elem2D (read over
  myDim_elem2D); `mass`/`inv_thickness`/`ice_nod` myDim_nod2D (read owned-only);
  rhs scratch is node_size already (halo written unguarded, trap 6).
- Unique OUT_DIR per SLURM job (`feedback_unique_outdir_per_job`).
- Build via `bash -l configure.sh` (never mambaforge). Fortran build needs NO change —
  `ice_maEVP.F90` is always compiled; the reference flip is namelist-only.

## Validation Strategy
- **Dump gate** (`FESOM_EVP_DUMP_DIR`, existing infra both sides): inputs (a_ice, m_ice,
  m_snow, elevation, u_w/v_w, stress_atmice), precomputes (pressure_fac, mass,
  inv_thickness, scaled rhs_a/rhs_m), per-substep u_aux/v_aux + sigma11/12/22 at substeps
  1/2/60/120, final uice/vice. Same dump points, same filenames, both codes; diff with the
  established scripts.
- **Cold-start caveat** (KPP lesson): step 1 has uice=0 and spun-up-from-IC ocean state —
  inputs should be bit-identical (PHC rank-noise floor at 16r); 0 unexplained
  disagreements is the bar, input-cross-check for mask-threshold flips (a_ice≈0.01 nodes).
- **Cross-rank**: 1/8/16r short mEVP runs at PHC-IC-accepted scale.
- **Climate**: 2yr C-mEVP vs NEW Fortran `work_mevp` (864r, dt=1800). Bar: TKE-class
  agreement (SST/SSS RMS ~0.005/0.003, biases ~0, non-drifting) + ice metrics (extent,
  volume, drift speed NH+SH seasonal cycles; `nan_to_num→0` before any time-averaging —
  `feedback_nan_mask_temporal_average`).
- **Diff-of-diffs**: (C-mEVP − C-EVP) spatial pattern ≈ (Fortran-mEVP − Fortran-EVP).
  Catches "option silently does nothing" — the two rheologies SHOULD differ visibly
  (smoother sigma, different ridging), and the C must reproduce that contrast, not just
  the absolute fields.
- **Ice-velocity VECTOR gate** (`feedback_magnitude_invariant_masks_rotation`): trap 3
  is a rotation-class difference, and extent/volume/drift-SPEED metrics are all
  rotation-blind — a wrong-angle drag would pass every one. M5 compares monthly-mean
  (uice, vice) per node: magnitude ratio ≈1 AND inter-vector angle ≈0 vs Fortran-mEVP.
- **Stability**: 5yr C-mEVP dt=1800; peak |uv| margin vs the 5.0 trip (autumn transient
  watch per `project_5yr_milestone`).

## Implementation Steps

### Phase M0: Branch + baseline capture + Fortran reference runs (launch queue FIRST — wall-clock dominates)
**Files:** branch `mevp`; Fortran side `work_mevp/` (new), `ice_maEVP.F90` dump calls
(uncommitted); `jobs/job_mevp_fdump`, `jobs/job_mevp_baseline`; `docs/mevp_reference_namelists/`
- [ ] Create branch `mevp` in `fesom2_port_zstar` (VERIFIED: `main`==`zstar`==2714071, so
      this is a same-commit branch creation — no working-tree change, zstar/TKE code is
      merged content); commit this plan.
- [ ] **Capture the byte-gate baseline BEFORE any src change**: `bash -l configure.sh` at
      2714071, run the pinned gate config — **16r dist_16, 200 steps, dt=500, knob unset**
      — archive snapshots + run.log to `/work/ab0995/a270088/port/mevp/baseline_2714071/`
      (zstar-Z0/TKE-T0 precedent). ALL later byte-gates re-run this exact config and
      byte-compare against this capture.
- [ ] Copy `work_core` → `work_mevp` (fresh — NOT the stale `work_core_evp1`); flip
      `whichEVP=1` (ONLY change); verify `alpha_evp=250`, `beta_evp=250`,
      `evp_rheol_steps=120` present; unique ResultPath; submit the 2yr dt=1800 run.
- [ ] Gate: run log echoes `EVP scheme option= 1` (fesom_module.F90:323) and runs clean
      past the first ice step.
- [ ] Archive namelists → `docs/mevp_reference_namelists/` + PROVENANCE.md (note the
      stale work_core_evp1 explicitly).
- [ ] Fortran dump instrumentation: `use fesom_evp_dump` in `EVPdynamics_m`; add the
      env-gated dump calls at the points listed in Validation Strategy (substeps
      1/2/60/120 via the shortstep counter). Uncommitted working-copy patch, KPP/TKE-style.
- [ ] Short instrumented Fortran run (16r, dist_16, 2 steps, `FESOM_EVP_DUMP_DIR` set).
- [ ] **Validate:** dump files present, finite, plausible (pressure_fac scale ~1e4·m_ice;
      u_aux → cm/s-dm/s scale by substep 120).

### Phase M1: C skeleton — knob, state, bc_index_nod2D, dispatch (NO default-path change)
**Files:** Modify `src/fesom_ice_types.h`, `src/fesom_ice.c`, `src/fesom_mesh.h`,
`src/fesom_mesh.c` (free), `CMakeLists.txt`; Create `src/fesom_ice_maevp.{c,h}` (stub)
- [ ] Parse `FESOM_WHICH_EVP` in `fesom_ice_init` — numeric `atoi`, the
      `FESOM_DIAG_GID` precedent (fesom_ice.c:392-393), NOT the per-step presence
      pattern at :297 — replacing the hard `whichEVP=0` at :68: unset/`0` → 0; `1` → 1;
      anything else → `MPI_Abort` with "whichEVP=N not ported (0=EVP, 1=mEVP)". Init-time
      because the value gates the conditional allocations below.
- [ ] Struct fields: `alpha_evp=250.0`, `beta_evp=250.0` set unconditionally at init
      (namelist mirror); `uice_aux`/`vice_aux` alloc node_size **zero-init ONLY when
      whichEVP≠0** (MOD_ICE.F90:743-748 mirror) + free; mEVP work arrays (Technical
      Details) alloc'd under the same condition.
- [ ] `bc_index_nod2D` — **REBUILD the existing array** (it already exists:
      fesom_mesh.h:141 `real_t*`, built at fesom_ice.c:199-209; keep `real_t`, do NOT
      add an `int*` duplicate): replace the `edge_tri[ed*2+1] >= 0` interior test with
      `myList_edge2D[ed] <= edge2D_in` (MOD_ICE.F90:891-894 mirror; std-EVP BC-loop
      convention fesom_ice_evp.c:453-455 — the edge_tri test misclassifies inter-rank
      seam edges as boundary and would hard-zero mEVP velocity there). Loop bound stays
      `0..myDim_edge2D-1` (owned edges only, trap 9). Fix the now-wrong "Both schemes
      produce the same set" comment at fesom_ice.c:195-198. No exchange (Fortran has
      none; mEVP reads owned nodes only).
- [ ] Grep-verify `bc_index_nod2D` still has NO reader other than the new mEVP code
      (today: build + free only) — the rebuild is then provably byte-neutral for the
      default path.
- [ ] Dispatch at fesom_ice.c:327: `whichEVP==1` → `fesom_ice_evp_dynamics_m()` stub that
      aborts "mEVP port pending (M2)"; `==0` path UNTOUCHED.
- [ ] Update the header-comment contracts in fesom_ice_types.h/fesom_ice.c (whichEVP no
      longer "locked at 0").
- [ ] **Validate:** byte-gate — pinned config (16r/200/dt=500, knob unset) byte-identical
      to the M0 `baseline_2714071` capture (snapshots + run.log); `FESOM_WHICH_EVP=2` and
      `=1` both abort with the right messages; build clean.

### Phase M2: Port `EVPdynamics_m` → `fesom_ice_evp_dynamics_m()`
**Files:** `src/fesom_ice_maevp.c` (the function + C dump points)
- [ ] Block 1 — setup constants + aux init over myDim+eDim (traps 1, 8).
- [ ] Block 2 — inlined ssh2rhs, levitating branch only; myDim-only zeroing; unguarded
      3-node assembly (traps 6, 7).
- [ ] Block 3 — node precompute: masks, `max(·,9.0)` limiter, `M/((1+M²)·area)`,
      in-branch rhs scaling (traps 4, 5, 7).
- [ ] Block 4 — element precompute: mean-msum mask, pressure_fac with det2 and NO 0.5
      (traps 2, 4).
- [ ] Block 5 — main loop: strain rates with metric terms, additive delta_min, sigma
      updates with the 0.5 in 11/22 only, guarded stress assembly, full-rdt drag,
      `+β·u_aux`, `bc_index×det` implicit solve, NO theta_io, edge BC, two
      `fesom_exchange_nod2D` calls, end-of-substep rhs zeroing (traps 1, 2, 3, 9, 10).
- [ ] Block 6 — final copy over myDim+eDim, no extra exchange (trap 8).
- [ ] C dump points mirroring M0's Fortran points exactly (same names/substeps).
- [ ] **Fidelity review: walk the 14-item traps checklist against the final diff, item by
      item, citing the C line for each.**
- [ ] **Validate:** build clean; 1r smoke (200 steps dt=1800, `FESOM_WHICH_EVP=1`): ice
      drifts (uice nonzero, evolving), no NaN, no HUGE-ice_strength stderr fires, a_ice/m_ice
      ranges sane vs std-EVP run; byte-gate default re-run.

### Phase M3: Substep-level dump-diff C vs Fortran
**Files:** `jobs/job_mevp_cdump`; reuse the EVP diff scripts (`reference_evp_dump_diagnostic`)
- [ ] C dump run: 16r dist_16, 2 steps, dt=1800, `FESOM_WHICH_EVP=1` +
      `FESOM_EVP_DUMP_DIR` — same partition as M0's Fortran dump run.
- [ ] Diff sequence (stop at first unexplained): inputs (expect bit-identical) →
      precomputes (pressure_fac/mass/inv_thickness/rhs) → substep 1 u_aux/sigma →
      substeps 2/60/120 → final uice/vice. Bar: ~1e-12 rel at substep 1, small monotone
      growth to 120; 0 unexplained diffs; mask-threshold flips resolved by input
      cross-check.
- [ ] If a live diff is amplified/uninterpretable: controlled replay — inject Fortran
      inputs (a/m_ice, elevation, u_w, stress) into the C before the routine
      (`feedback_controlled_replay_validation`).
- [ ] If 16r shows velocity diffs CONCENTRATED AT INTER-RANK SEAMS that 1r lacks:
      suspect `bc_index_nod2D` first (the M1 rebuild; dump it from both codes and diff
      directly — it is step-constant).
- [ ] **Validate:** documented diff table in `docs/validation_mevp/M3_dumpdiff.md`;
      byte-gate.

### Phase M4: Cross-rank consistency + 30-day stability
**Files:** `jobs/job_mevp_crossrank`
- [ ] 1/8/16r 30-day dt=1800 mEVP runs (unique OUT_DIRs); cross-rank a_ice/m_ice/uice
      means at the PHC-IC-accepted scale (sea-ice Phase E precedent: ~1e-4 rel).
- [ ] max|uv| (ocean) and ice-speed sanity vs the std-EVP C run (ice/ocean speed ratio
      diagnostic from `feedback_ice_strength_halffactor`).
- [ ] **Validate:** clean exits, no NaN/uv>5, consistency table; byte-gate.

### Phase M5: 2yr climate validation vs work_mevp
**Files:** `jobs/job_mevp_2yr`; reuse `validation_2yr.sh` + climate-compare scripts
- [ ] C mEVP 2yr dt=1800 864r (unique OUT_DIR).
- [ ] Ocean: SST/SSS RMS + bias maps + drift check vs Fortran work_mevp. Bar: TKE-class
      (~0.005/0.003, biases ~0, non-drifting).
- [ ] Ice: NH+SH extent/volume/drift-speed seasonal cycles C vs Fortran-mEVP
      (`nan_to_num→0` before time means).
- [ ] Ice-velocity VECTOR check (Validation Strategy): per-node monthly-mean (uice,vice)
      magnitude ratio + inter-vector angle vs Fortran-mEVP — not just speed
      (`feedback_magnitude_invariant_masks_rotation`).
- [ ] **Diff-of-diffs**: (C-mEVP − C-EVP) vs (F-mEVP − F-EVP) spatial patterns for ice
      thickness/drift and SST — the contrasts must match, and must be visibly nonzero
      (proves the option is live and faithful, not inert).
- [ ] **Validate:** `docs/validation_mevp_2yr/RESULTS.md` with acceptance numbers.

### Phase M6: 5yr stability + acceptance + handoff
**Files:** `jobs/job_mevp_5yr` (clone production_864_5yr.sh)
- [ ] 5yr C mEVP dt=1800: clean run, peak |uv| margin vs the 5.0 trip (compare the 4.65
      autumn transient of the EVP baseline), ice fields bounded.
- [ ] Acceptance review against Overview: knob works (0/1/abort), off-switch
      byte-identical, dump-diff ≤ bars, cross-rank consistent, 2yr climate + ice metrics
      pass, diff-of-diffs matches, 5yr stable.
- [ ] Tag `mevp-validated-<date>`; update auto-memory (`project_mevp_port_state`, link
      lessons); move this plan to `docs/plans/completed/`; write next-session prompt.
- [ ] Merge `mevp` → `main` decision deferred to the user (zstar/TKE precedent).

## Technical Details
- **Function signature**: `void fesom_ice_evp_dynamics_m(fesom_ice *ice, struct
  fesom_partit *partit, struct fesom_mesh *mesh)` — same shape as the std-EVP entry.
- **Per-call scratch** (Fortran automatic arrays `inv_thickness/mass/ice_nod` [myDim_nod2D]
  and `pressure_fac/ice_el` [myDim_elem2D]): allocated as persistent `ice->work` arrays
  when whichEVP=1 (the `inv_areamass`/`inv_mass` precedent — heap, not stack: 1r CORE2 is
  ~127k nodes). All five are fully (re)written every call before use, so persistence is
  semantics-neutral.
- **Index translation**: Fortran `elnodes(k) <= myDim_nod2D` → C `elem_nodes[3*el+k] <
  myDim_nod2D`; 1-based `area(1,i)` → `area[n*nl+0]`; `gradient_sca(1:3,el)/(4:6,el)` →
  `gradient_sca[6*el+0..2] / [6*el+3..5]`.
- **Exchange semantics**: Fortran `exchange_nod_begin(u_aux,v_aux)+exchange_nod_end` with
  the rhs-zero loop in between (comm/compute overlap) → C: two blocking
  `fesom_exchange_nod2D` calls then the zero loop. Result-identical; the zero loop touches
  only u_rhs/v_rhs, not the exchanged fields.
- **Memory delta**: default path +1 int node array (`bc_index_nod2D`, built always like
  post-2023 Fortran); mEVP path +2 real node arrays (aux) + the 5 work arrays. No change
  to any existing array.
- **whichEVP=1 + FESOM_ALE=zstar**: composes by construction — the levitating ssh2rhs
  branch is the executed Fortran path for linfs AND zstar (use_floatice=false), and the
  zstar code IS in this tree (verified: zstar branch tip == main tip == 2714071, the TKE
  ff-merge carried it). Not a validated configuration in this plan's matrix
  (user-approved single-knob scope); noted in the results doc as untested-combination.

## Post-Completion
*Informational — external/manual:*

**Run matrix (approved 2026-06-11):** single-knob only — `work_mevp` (Fortran, 2yr) vs
C-mEVP 2yr + 5yr C-mEVP stability. NO combined zstar/TKE/mEVP Fortran reference (explicitly
descoped by the user).

**Manual verification:**
- mEVP-vs-EVP science contrast (drift speed distributions, deformation/ridging patterns
  vs the Bouillon/Danilov literature) — optional science check, not a port gate.

**Decisions deferred to the user:**
- Merge `mevp` → `main` timing.
- Whether aEVP (whichEVP=2) gets ported later (the `_a` routines + alpha/beta arrays are
  the prepared seam; the dispatcher abort names it).
- Whether mEVP ever becomes default (standard EVP remains default per current decision).
