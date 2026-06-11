# Port the CVMix TKE vertical mixing scheme (mix_scheme='cvmix_TKE', classical TKE only)

> Plan 2 of 2 from the 2026-06-10 brainstorm (companion: `20260610-zstar-vertical-coordinate.md`,
> which is implemented FIRST per the approved order).
> **Revised 2026-06-10 after plan-review.** Corrections folded in: `handle_old_vals` /
> `old_KappaM` / `old_KappaH` are DEAD in the executed path (no old-value blending —
> `integrate_tke` never reads them; the driver's `tke_Av_old`/`tke_Kv_old` locals feed unused
> args) → dropped from the port; explicit owned-only horizontal loop bounds pinned;
> `iw_diss` is read unconditionally into the `tke_Tiwf` diagnostic → pass a REAL zero column
> (not skip); diag scratch has internal read-after-write dependencies (`tke_Tdif`→`tke_Twin`;
> `P_diss_v`/`K_diss_v` seed `Tbpr`/`Tspr`) → compute the FULL diag set in local column
> scratch in Fortran order, copy out at the end only; element-Av halo comes from the shared
> post-`mo_convect` exchange (no in-driver element exchange); T0 must verify the reference
> run LOG echoes `tke_cd=3.75` (namelist actually consumed; `__cvmix` actually compiled).
> All design decisions are user-approved (2026-06-10): two C files mirroring the Fortran
> split; diagnostics always computed per-column but globally stored only under
> `FESOM_TKE_DIAG=1`; Langmuir/IDEMIX/Dirichlet gate-only; validated vs a NEW Fortran
> linfs+TKE reference.

## Overview

Port the classical TKE vertical mixing scheme (Gaspar et al. 1990, Blanke & Delecluse
mixing length — the ICON/Brüggemann CVMix-fork implementation FESOM2 vendors in
`src/cvmix_driver/`) to C as a **strictly faithful line-by-line port**, runtime-selected by
`FESOM_MIX_SCHEME=TKE` (KPP stays the default; PP stays available). TKE is a prognostic
column scheme: it advances a turbulent-kinetic-energy field `tke[nl×N]` each step (shear +
buoyancy production, Kolmogorov dissipation `c_eps·e^{3/2}/L`, implicit vertical diffusion of
TKE) and yields `KappaM`→`Av` and `KappaH`→`Kv`.

**GUIDING PRINCIPLE (non-negotiable):** THIS IS A PORT. Reproduce every formula, constant,
loop bound, branch, halo exchange, and ordering exactly. Namelist values beat module defaults
(`feedback_namelist_over_codedefault`) — the reference `namelist.cvmix` sets **`tke_cd=3.75`**
while the module default is `1.0`; port 3.75. "Port what the reference uses"
(`feedback_port_used_terms`): `only_tke=.true.`, Langmuir OFF, Dirichlet BCs OFF, IDEMIX
absent → those bodies are gate-only. Verified-dead arguments (never read in the executed
body) are dropped from the C signature with the verification cited — see "Dropped as dead".

**KPP and PP stay byte-identical to the pre-TKE baseline** (v1.0, or zstar-validated HEAD if
zstar landed first) — verified by a byte-gate in every phase that touches shared code.

## Reference TKE configuration (VERIFIED 2026-06-10 — values to port)

Fortran reference = `work_linfs_2yr` config with ONLY `mix_scheme` flipped to `'cvmix_TKE'`
(+ `namelist.cvmix` copied into the work dir). The Intel build already has `CVMIX:BOOL=ON` —
**no rebuild needed** (but T0 verifies the binary actually has `__cvmix`, see below).

From `work_core/namelist.cvmix` `&param_tke` (re-confirm against the actual reference work
dir in T0) + driver/module defaults:

| Symbol | Value | Source | Consequence |
|---|---|---|---|
| `mix_scheme` | `'cvmix_TKE'` → `mix_scheme_nmb=5` | oce_setup_step.F90:151 | dispatch branch ==5 |
| `tke_c_k` | `0.1` | namelist.cvmix | KappaM = c_k·L·√e |
| `tke_c_eps` | `0.7` | namelist.cvmix | dissipation coefficient |
| `tke_alpha` | `30.0` | namelist.cvmix | max L ≤ √(alpha·e)/N stability bound |
| `tke_mxl_min` | `1.0e-8` | namelist.cvmix | mixing-length floor |
| `tke_kappaM_min` / `max` | `0.0` / `100.0` | namelist.cvmix | viscosity clamps |
| `tke_cd` | **`3.75`** | namelist.cvmix (**module default 1.0 — namelist wins**) | surface-TKE coefficient in the **Neumann** BC `forc(1) += cd·forc_tke_surf^{3/2}/dzt(1)` (cvmix_tke.F90:791-796). The "3.75=Dirichlet" comment is advisory only — the Neumann branch executes with 3.75 (plan-review confirmed) |
| `tke_surf_min` | `1.0e-4` | namelist.cvmix | min surface TKE |
| `tke_min` | `1.0e-6` | namelist.cvmix | min interior TKE |
| `tke_mxl_choice` | `2` | namelist.cvmix | Blanke–Delecluse mixing length (only implemented option) |
| `tke_dolangmuir` | `.false.` | namelist.cvmix | Langmuir body GATE-ONLY |
| `tke_only` | `.true.` (default; flipped only by nmb==56) | gen_modules_cvmix_tke.F90:53,226 | IDEMIX inputs stay zero |
| `use_ubound_dirichlet` / `use_lbound_dirichlet` | `.false.` / `.false.` (defaults, not in namelist) | gen_modules_cvmix_tke.F90:54-55 | Neumann BC branches are the executed path |
| `timerelax_tke`, `relne`, `relax` | `.false.`, 0.4, 0.6 (defaults) | driver module | driver namelist params that never reach cvmix — confirm unused in T2, gate-only note |
| `i_vert_diff` | `.true.` | namelist.tra:74 | implicit TDMA consumes Kv (same as KPP) |
| `mo_convect` | called AFTER calc_cvmix_tke | oce_ale.F90 dispatch tail | keep convective adjustment |
| `density_0`, `g`, `dt` | as in C constants | passed as rho_ref/grav/dtime | identical values |

### Dropped as dead (plan-review verified — cite these, don't re-litigate)
- `handle_old_vals`: consumed only by `cvmix_update_tke`, which is called only from the
  commented-out block in `tke_wrap` (cvmix_tke.F90:399-409). `integrate_tke` never reads it.
  **There is NO old-value blending**: `KappaM_out = min(KappaM_max, c_k·mxl·sqrttke)` (:704),
  `KappaH_out = KappaM_out/prandtl` (:720) — fresh overwrites.
- `old_KappaM`/`old_KappaH` args: declared `intent(in)` (:474-475), never appear in the body.
  The driver's per-column locals `tke_Av_old`/`tke_Kv_old` (gen:433-434) feed only these dead
  args → NOT ported. The ONLY persistent input is `tke` itself (`tke_old`, read at :668).
- `max_nlev` arg: unused in the body (only dead `tke_wrap`/commented Kappa_GM).
- `i`/`j`/`tstep_count` args: used only inside an `if (.false.)` debug block (:916-986).
- `tke_userdef_constants`: absent at the call → module saved struct; C uses one static
  params struct.
- Module-level 1D `tke_Av_old/tke_Kv_old/tke_old` allocatables: shadowed dead code.
- `tke_wrap` (:285-411), put_tke machinery (:991-1093), `tke_diss` (declared, never
  allocated/used by the driver).

➕ **T0 sub-task:** archive the exact `namelist.cvmix` + `namelist.oce/tra` of the reference
run → `docs/tke_reference_namelists/` + PROVENANCE.md; re-confirm every value above; **grep
the reference run LOG for the `init_cvmix_tke` echo block** (gen:228-241 prints `tke_only`,
`tke_cd`, …) and confirm `tke_cd = 3.75` and `tke_only = T` — this proves both that
`namelist.cvmix` was actually read (silent fallback to cd=1.0 otherwise, gen:216-223) and
that the binary has `__cvmix` (otherwise `mix_scheme_nmb=5` silently falls through the
`#if`-guarded dispatch and NO mixing is computed — a plausible-looking broken reference).

## Context (from discovery, 2026-06-10; spot-checked by plan-review)

### Fortran source (port reference) — `port2/fesom2/src/cvmix_driver/`
| File / lines | Routine | Role | Status |
|---|---|---|---|
| `cvmix_tke.F90:101-281` | `init_tke` | stuff namelist params into the module `tke_type` via put_tke; defaults for absent optionals | port (collapse to static struct assignment) |
| `cvmix_tke.F90:415-987` | `integrate_tke` | THE CORE (~570 lines): mixing length (mxl_choice=2 with mxl_min, sqrt(2e)/N bound via alpha, bounded-by-distance walls), Prandtl number Pr(Ri), shear/buoyancy production, dissipation, implicit tridiagonal solve for tke_new (tke_old read at :668), tke_min/tke_surf_min floors, **Neumann** surface BC with cd (:791-796), kappaM clamps, KappaM/KappaH fresh outputs (:704,:720), only_tke/IDEMIX/Langmuir/Dirichlet branches, ALL budget diagnostics (`intent(out)`, unconditional) | port (gate-only the non-executed branches; dead args dropped) |
| `gen_modules_cvmix_tke.F90:130-273` | `init_cvmix_tke` | allocate state + diag arrays, read &param_tke, call init_tke | port (allocation split by diag flag) |
| `gen_modules_cvmix_tke.F90:279-507` | `calc_cvmix_tke` | per-node column loop (**`node_size=myDim_nod2D`**, :296,311 — OWNED ONLY): normstress=√(stress_node_surf²)/ρ₀; vshear2 from UVnode/(ΔZ_3d_n)²; bvfreq2=bvfreq(nun+1:nln); dz_trr (half-thickness surface/bottom); per-column tke_old copy (LOCAL); cvmix call with nun:nln(+1) slices, nlev=nln-nun+1; zero `tke_Av/Kv` at `nun` and `nln+1` after the call; `exchange_nod(tke_Kv)` (:493); `Kv = tke_Kv` (full copy); `exchange_nod(tke_Av)` (:498); `Av(nz,elem)=Σ tke_Av(nz,elnodes)/3` for **elem=1..myDim_elem2D** (:500), nz=ulevels+1..nlevels-1; NO element exchange in the driver | port |
| `oce_setup_step.F90:151,185-187` | dispatch init | nmb=5 → init_cvmix_tke (inside `#if __cvmix`) | port (C: init in fesom_main when TKE selected) |
| `oce_ale.F90` dispatch tail | `calc_cvmix_tke` + `mo_convect` after | same slot as KPP/PP in C `fesom_step.c` | port |

### Argument semantics pinned (cvmix_tke.F90:461-520; plan-review verified)
- Budget diagnostics (`tke_Tbpr/Tspr/Tdif/Tdis/Twin/Tiwf/Tbck/Ttot/Lmix/Pr`, `cvmix_int_1/2/3`)
  are plain `dimension(nlev+1), intent(out)` — **NOT optional, always computed**. None of them
  is read to produce `tke_new`/`KappaM_out`/`KappaH_out` (**plan-review confirmed: the
  diag-gating design is physics-airtight**). BUT there are internal read-after-write
  dependencies AMONG the scratch quantities: `tke_Tdif` is read when building `tke_Twin`
  (:836-837, :882, :891), and `P_diss_v`/`K_diss_v` are PHYSICS scratch (feed `forc` :734)
  that also seed `tke_Tbpr`/`tke_Tspr` (:876-877). ⇒ C computes the FULL diag set into LOCAL
  column scratch in Fortran order, every call; the `FESOM_TKE_DIAG` flag gates only the global
  `[nl×N]` arrays and the end-of-column copy-out. **Never write through nullable pointers
  mid-computation** — NULL would corrupt `tke_Twin` or segfault.
- `tke_plc` (Langmuir): read only under `if (l_lc)` (:737-739) — purely gated.
  `E_iw`/`alpha_c`: read only under `if (.not. only_tke)` (:710-712) — purely gated.
  **`iw_diss`: read UNCONDITIONALLY at :898** (`tke_Tiwf = iw_diss`) ⇒ the C must pass a REAL
  zero array of length nlev+1 (not skip it); result `tke_Tiwf=0`, matching the Fortran whose
  `tke_in3d_iwdis` stays 0 under only_tke.
- `tke_userdef_constants` optional → absent → module saved constants (C: one static struct).

### C integration points (target)
| C file | Current state | Change |
|---|---|---|
| `src/fesom_cvmix_tke.{c,h}` | NEW | literal port of init_tke + integrate_tke (pure column math, no MPI/mesh; dead args dropped per the verified list) |
| `src/fesom_tke.{c,h}` | NEW | driver port: `struct fesom_tke` state + per-node loop + Av/Kv wiring (the `fesom_kpp` precedent: state struct threaded via step ctx, alloc/free in main) |
| `src/fesom_step.c` | `FESOM_MIX_SCHEME=PP|KPP` dispatch (KPP default), `mo_convect` + the shared unconditional Kv/Av exchanges after (`fesom_exchange_elem3D(aux->Av)` at :266) | + third value `TKE` (accept alias `cvmix_TKE`); TKE branch calls `fesom_tke_mixing`; `mo_convect` after (slot exists); **element-Av halo comes from the existing shared post-mo_convect exchange — do NOT add an element exchange in the TKE driver** |
| `src/fesom_main.c` | env options block; kpp alloc/free precedent | + FESOM_TKE_DIAG parse; tke alloc/init/free when selected |
| `src/fesom_aux.c/.h` | single `aux->Kv[N·nl]` (TDMA consumer), `aux->Av[E·nl]` | NO change — TKE writes them like KPP does (Kv full copy, Av node→elem) |
| `src/fesom_tracer_diff.c` / `fesom_momentum.c` | implicit TDMA reads aux->Kv / aux->Av | NO change |
| inputs | `stress_node_surf`, `UVnode`, `bvfreq` (N²-smoothed, `feedback_bvfreq_smoothing_gap`), `zbar_3d_n`/`Z_3d_n`, `hnode` | all exist (KPP era); under zstar they are live — TKE reads them per step (no caching) |
| `src/fesom_io_stream*` | 17-var catalog + FESOM_IO_CONFIG | + register `tke` (+ diag fields when FESOM_TKE_DIAG=1), off by default |

### State & memory (the user's explicit concern)
Always allocated (`[nl×N]`, N=myDim+eDim, stride nl):
- `tke` — **PROGNOSTIC state** (init 0 → floored to tke_min by the first call; `tke_old` is
  read at :668 to drive mxl/KappaM). The Fortran comment calls it a "diagnostic" — it is not;
  never fold it into the diag set.
- `tke_Av`, `tke_Kv` — needed globally ONLY because `exchange_nod` + the full-array
  `Kv=tke_Kv` copy require halo storage; **their prior-step values are never read** (no
  old-value blending — see "Dropped as dead").
- 2D: `tke_forc2d_normstress/botfrict/rhosurf` (botfrict/rhosurf zeroed every call and
  passed — port as the Fortran does).

**Gated by `FESOM_TKE_DIAG=1`** (default OFF, pointers NULL): the 10 budget arrays
`tke_Tbpr/Tspr/Tdif/Tdis/Twin/Tiwf/Tbck/Ttot/Lmix/Pr` + 3 `cvmix_dummy_*` `[nl×N]` each.

NOT allocated at all (zero/dead by construction): `tke_in3d_iwe/iwdis/iwealphac`,
`tke_langmuir`, `langmuir_*` (per-column zero scratch of length nlev+1 passed instead —
value-identical under only_tke/dolangmuir=off, incl. the unconditional `iw_diss→tke_Tiwf`
copy), module-level `tke_Av_old/tke_Kv_old/tke_old` 1D allocatables (shadowed dead code),
`tke_diss` (declared, never allocated/used by the driver).

### Out of scope (gate-only / not ported)
- **IDEMIX coupling** (`tke_only=.false.`, nmb=56, `g_cvmix_idemix`) — abort-if-selected gate.
- **Langmuir** (`tke_dolangmuir=.true.` body, gen:367-429 + tke_plc use) — abort-if-set gate.
- **Dirichlet BC branches** (`use_*bound_dirichlet=.true.`) inside integrate_tke — gate +
  `#if`-excluded body (KPP `ddmix` precedent: build fails loudly if enabled).
- **`tke_wrap`**, put_tke machinery, `timerelax_tke`/`relne`/`relax` (driver namelist params
  that never reach cvmix — confirm in T2, then gate-only note).
- **cvmix library proper** (CVMix-src) — nothing else from it is needed; `cvmix_kinds` →
  `double`, `cvmix_put` → struct init.
- **Restart** of tke state — no restart system in C.

### Key Fortran invariants (must hold in C)
1. **Same dispatch slot as KPP/PP**: bvfreq (smoothed) → scheme → `mo_convect` → shared
   Kv/Av exchanges in `fesom_step.c`. TKE branch must leave PP/KPP paths untouched.
2. **Horizontal loop bounds (plan-review finding — the project's signature bug class):**
   main column loop runs **owned nodes only** (`node = 1..myDim_nod2D`, gen:296,311; C:
   `0..myDim_nod2D-1`); element-Av loop runs **owned elements only** (`myDim_elem2D`, :500).
   `tke_Kv`/`tke_Av` halos come from `exchange_nod` BEFORE the `Kv=tke_Kv` copy / node→elem
   average. Do NOT extend either loop to the halo, and do NOT exchange `tke`
   (`feedback_write_loops_halo` cuts both ways — no missing AND no extra exchanges).
3. **Column slicing**: nlev=nln-nun+1 layers (nln=nlevels_nod2D-1, nun=ulevels=1), interface
   arrays nlev+1; `dzw=hnode(nun:nln)` (layers), `dzt=dz_trr(nun:nln+1)` with half-thickness
   end values; vshear2/bvfreq2 zero outside nun+1..nln.
4. **Old-state semantics**: per-column copy of `tke` (only) taken BEFORE the call; outputs
   written back to the global slabs. KappaM/KappaH are fresh overwrites — no blending.
5. **Post-call zeroing**: `tke_Av(nln+1)=tke_Av(nun)=0`, same for tke_Kv — BEFORE the
   exchanges; element Av averaged over nz=ulevels(elem)+1..nlevels(elem)-1 only; element-Av
   halo filled by the SHARED post-`mo_convect` `fesom_exchange_elem3D(aux->Av)`
   (fesom_step.c:266), exactly as for KPP/PP — no in-driver element exchange.
6. **Surface forcing**: `forc_tke_surf = |stress_node_surf|/density_0` (per node, NOT the
   commented element-interpolation variant); forc_rho_surf=0, bottom_fric=0 every call.
7. **vshear2 denominator** uses `Z_3d_n` spacing (live under zstar), squared; computed at
   interfaces nun+1..nln from UVnode differences.
8. **i_vert_diff=.true.**: Kv feeds the implicit TDMA exactly as KPP's single Kv does;
   `fesom_tracer_diff.c` unchanged.
9. **Diag math is part of integrate_tke** — compute the FULL set into local column scratch
   in Fortran order every call (internal deps: Tdif→Twin; P_diss_v/K_diss_v→Tbpr/Tspr);
   copy out to nullable global pointers at the END of the column only.

## Development Approach
- Faithful port discipline: one Fortran routine per phase, dependency order, literal diff.
- **Validation = dual-instrumentation vs Fortran** (dump-diff; the KPP playbook applies
  nearly verbatim, incl. controlled replay if input noise amplifies).
- **KPP/PP byte-identical with TKE off** — discrete Validate item in every phase touching
  shared code.
- New node arrays `[N·nl]` stride nl, sized myDim+eDim (`feedback_tracer_stride_nl`,
  `feedback_array_size_vs_reader_loop`). Unique OUT_DIR per SLURM job.
- If zstar has landed: baseline byte-gate = the zstar-validated HEAD (linfs default), and the
  TKE validation runs stay **linfs**+TKE to isolate the mixing knob (approved matrix).

## Validation Strategy
- **Dual dump gate `FESOM_TKE_DUMP_DIR`** (both codes; Fortran gate uncommitted in
  `gen_modules_cvmix_tke.F90`, Intel tree, like the KPP gates): step-1 16r dt=1800 per-node
  columns of INPUTS (normstress, vshear2, bvfreq2, dz_trr, tke_old) and OUTPUTS
  (tke, tke_Av, tke_Kv, final aux Kv + elem Av), plus the 10 budget diagnostics (both codes
  always compute them — dump from the C's column scratch).
- **Step-1 caveats (KPP lessons)**: cold start uv=0 → vshear2≡0 → shear production untested at
  step 1 (covered by the multi-step dump + climate, K3 precedent); inputs differ only at the
  PHC rank-noise floor; threshold/floor outputs (tke_min flooring, kappaM clamps) validated by
  input-cross-check — 0 unexplained disagreements is the bar. Dump steps 2–3 as well (shear on).
- **Controlled replay** (`feedback_controlled_replay_validation`): if live diffs are
  amplified/uninterpretable, inject Fortran-dumped column inputs into the C before the call
  (FESOM_TKE_REPLAY_DIR, same 16r partition) — expect ≤1e-12 like KPP K6/K7.
- **Diag-flag invariance**: FESOM_TKE_DIAG=0 vs =1 runs byte-identical in all model state.
- **Byte-gates**: KPP (default) and PP runs bit-identical to baseline with TKE off.
- **Cross-rank**: 1/8/16r TKE short runs.
- **Climate**: 2yr C(linfs+TKE) 864r vs NEW Fortran `work_linfs_tke`. Bar: KPP-class
  agreement (SST/SSS RMS ~0.005/0.002, biases ~0). Scheme contrast: C+TKE vs Fortran-KPP must
  show a clearly larger signal than C+TKE vs Fortran-TKE (the comparison resolves schemes).

## Implementation Steps

### Phase T0: Scaffolding — dispatch, state struct, dump harness, reference runs (NO behavior change)
**Files:** Create `src/fesom_tke.{c,h}`, `scripts/tke_dump_diff.py`, `jobs/job_tke_*`; Modify `src/fesom_step.c`, `src/fesom_main.c`, `CMakeLists.txt`; Fortran `gen_modules_cvmix_tke.F90` (dump gate, uncommitted)
- [x] `struct fesom_tke` (tke/tke_Av/tke_Kv slabs + 2D forc arrays + nullable diag pointers +
      param struct), threaded via step ctx like `fesom_kpp`; `FESOM_MIX_SCHEME=TKE` third
      branch (abort stub), `FESOM_TKE_DIAG` parse.
- [x] **Launch the Fortran reference runs NOW** (queue time): (a) `work_tke_dump` 16r dt=1800
      3 steps with the dump gate; (b) `work_linfs_tke` 2yr dt=1800 (copy work_linfs_2yr; flip
      `mix_scheme='cvmix_TKE'`; copy namelist.cvmix in). Same Intel binary (CVMIX=ON ✓).
      Unique OUT_DIRs. Archive namelists → `docs/tke_reference_namelists/` + PROVENANCE.md.
- [x] **Provenance hardening (plan-review finding 6):** grep both reference run logs for the
      `init_cvmix_tke` echo (`tke_cd = 3.75`, `tke_only = T`, gen:228-241) — proves the
      namelist was consumed (else silent cd=1.0 fallback) AND `__cvmix` is compiled in (else
      nmb=5 silently computes NO mixing — a plausible-looking broken reference).
- [x] Add `FESOM_TKE_DUMP_DIR` C writer skeleton + Fortran dump gate + `scripts/tke_dump_diff.py`.
- [x] **Validate:** TKE off → byte-identical to baseline HEAD (16r/200/dt=500). Fortran TKE
      dump run sane (tke fields finite, Kv profiles plausible, echo block confirmed).

### Phase T1: `integrate_tke` — the column core (literal port)
**Files:** Create `src/fesom_cvmix_tke.{c,h}`; `scripts/tke_column_reference.py` (optional independent check)
- [x] Port `init_tke` as struct initialization (all params incl. cd=3.75 namelist value).
      `handle_old_vals` is NOT ported (dead — see "Dropped as dead").
- [x] Port `integrate_tke` (:415-987) literally: mixing length (mxl_choice=2: sqrt(2e)/N
      bound with alpha=30, mxl_min floor, wall-bounded growth), Pr(Ri), production terms,
      dissipation, the implicit tridiagonal TKE solve (tke_old :668), floors
      (tke_min/tke_surf_min), **Neumann** surface BC with cd (:791-796), kappaM clamps,
      fresh KappaM/KappaH outputs (:704,:720). C signature DROPS the verified-dead args
      (old_KappaM/old_KappaH, max_nlev, i/j/tstep_count, tke_userdef_constants) — each
      documented with the deadness citation. Diag terms computed into LOCAL column scratch
      in Fortran order (Tdif before Twin; P_diss_v/K_diss_v seed Tbpr/Tspr), copied to
      caller-provided NULLABLE arrays at the end only. Gate-only: Dirichlet branches
      (`#if`+abort); IDEMIX/Langmuir inputs as zero columns — `E_iw`/`alpha_c`/`tke_plc` are
      branch-gated (verified), `iw_diss` is a REAL zero array (read unconditionally at :898
      into tke_Tiwf).
- [x] **Validate:** offline column harness: feed BOTH codes identical synthetic columns
      (stable/unstable/sheared/shallow/deep, incl. nlev edge cases) via the dump/replay gates;
      C-vs-Fortran outputs (tke_new, KappaM, KappaH + all 10 budget terms) ≤1e-12 (KPP K1/K2
      precedent). Byte-gate: TKE still not callable live — baseline unchanged.

### Phase T2: `calc_cvmix_tke` driver — column assembly + Av/Kv wiring
**Files:** Modify `src/fesom_tke.{c,h}`, `src/fesom_step.c`
- [x] Port the per-node loop **over owned nodes only** (invariant 2): normstress, vshear2
      (UVnode/Z_3d_n), bvfreq2, dz_trr, per-column `tke_old` copy (local), the cvmix call with
      exact slicing (nun:nln+1, nlev), post-call zeroing (invariant 5), forc arrays zeroed
      each call. Confirm timerelax_tke/relne/relax never reach the cvmix call — gate-only note.
- [x] Wire outputs: `exchange_nod(tke_Kv)` → `aux->Kv` full copy; `exchange_nod(tke_Av)` →
      node→elem mean over interior levels into `aux->Av` (**owned elements only**; halo Av
      comes from the shared post-`mo_convect` exchange at fesom_step.c:266 — no new exchange);
      `mo_convect` after (dispatch slot).
- [x] Remove the T0 abort stub → TKE is LIVE.
- [x] **Validate (dump-diff step 1, 16r dt=1800 vs `work_tke_dump`):** inputs first
      (normstress/vshear2/bvfreq2/dz_trr — expect bit-identical or PHC-floor), then outputs
      (tke/Av/Kv + budget terms) — 0 unexplained, input-cross-check for floor/clamp flips;
      controlled replay if amplified. Steps 2–3 dumps exercise shear production. Byte-gates:
      KPP default run AND PP run bit-identical to baseline.

### Phase T3: Diagnostics storage + IO
**Files:** Modify `src/fesom_tke.c`, `src/fesom_io_stream_dispatch.c` (+ io config)
- [x] FESOM_TKE_DIAG=1: allocate the 13 arrays, end-of-column copy-out; register `tke` +
      budget fields in the stream catalog (selectable via FESOM_IO_CONFIG; default streams
      unchanged).
- [x] **Validate:** diag-on vs diag-off runs byte-identical in ALL model state (snapshots);
      diag outputs finite + budget closes (`tke_Ttot ≈ Σ terms` — the scheme's own
      consistency); memory delta measured and recorded (the point of the flag).

### Phase T4: Multi-rank + short stability
- [ ] 1/8/16r TKE 30-day runs — cross-rank consistency at PHC-noise scale; exchange-and-
      compare probe on tke_Kv/tke_Av halos (and confirm NO probe signal on `tke` — it must
      not need an exchange).
- [ ] 30-day dt=1800 stability smoke (max|uv|, MLD sanity in convection regions vs Fortran TKE).
- [ ] **Validate:** clean runs, no NaN/uv>5; byte-gates re-checked.

### Phase T5: 2yr climate validation vs work_linfs_tke
**Files:** `scripts/tke_climate_compare.py`, `jobs/job_tke_2yr`
- [x] Run C linfs+TKE 2yr dt=1800 864r (unique OUT_DIR).
- [x] Compare vs Fortran `work_linfs_tke`: SST/SSS RMS + bias maps + seasonal cycle + MLD.
      Bar: KPP-class (RMS ~0.005/0.002, biases ~0, non-drifting).
- [x] Scheme contrast: C+TKE vs Fortran-KPP ≫ C+TKE vs Fortran-TKE (resolves schemes).
- [x] **Validate:** acceptance numbers in `docs/validation_tke_2yr/`.

### Phase T6: Acceptance + handoff
- [ ] Overview requirements met: TKE selectable; KPP/PP byte-identical when off; column core
      ≤1e-12 vs Fortran; climate matches the TKE reference; diag flag inert + functional.
- [ ] Update docs + auto-memory (`project_tke_port_state`); commit per-phase; tag suggestion
      `tke-validated-<date>`; move this plan to `docs/plans/completed/`.

## Technical Details
- **Files**: `fesom_cvmix_tke.c` = cvmix_tke.F90 (column math only; params struct + column
  pointers; diag pointers nullable, written only at end-of-column; dead args dropped);
  `fesom_tke.c` = gen_modules_cvmix_tke.F90 (state, MPI, mesh, dispatch glue). 1:1
  traceability for dump-diff review.
- **Dispatch**: `FESOM_MIX_SCHEME` ∈ {KPP(default), PP, TKE|cvmix_TKE}. Init in main when
  selected (oce_setup_step.F90:185-187 mirror).
- **Memory**: always-on TKE state = 3×[nl×N]+3×[N]; diag flag adds 13×[nl×N]. At CORE2 16r
  (~8k nodes/rank, nl=48) diag ≈ 40 MB/rank — the flag matters at NG-mesh scale (nl=70).
- **Zero-input scratch**: per-column `[nlev+1]` zeros for tke_plc/E_iw/alpha_c/iw_diss — NOT
  global arrays. iw_diss MUST be a real array (unconditional read :898); the others are
  branch-gated (verified).
- **Combined-coordinate note**: under FESOM_ALE=zstar, TKE reads live hnode/Z_3d_n per step —
  no TKE code change (the driver re-derives dz_trr/vshear2 each call by construction).

## Post-Completion
*Informational — external/manual:*

**Run matrix (approved 2026-06-10, shared with the zstar plan):**
- Fortran `work_zstar_tke` 2yr (both knobs) + C zstar+TKE 2yr compare; then the 5yr C
  zstar+TKE stability run (autumn uv-transient margin re-measure). These close the matrix:
  zstar+KPP ✓ (zstar plan), linfs+TKE ✓ (this plan), zstar+TKE ✓ (combined).

**Manual verification:**
- MLD maps vs observations (Labrador/Weddell/GIN) — TKE-vs-KPP physics payoff, science check.
- TKE budget output (diag flag on) qualitative check vs ICON/FESOM literature profiles.

**Decisions deferred to the user:**
- Whether TKE ever becomes default (KPP remains default per current decision).
- Whether to port IDEMIX (`cvmix_TKE+cvmix_IDEMIX`) later — the iw-input plumbing is staged
  by the gate-only seams.
- Upstream PR for any Fortran bugs found during the port (step-1-forcing precedent:
  default remains "the Fortran is right — search the C first").
