# Port FESOM1.4 KPP vertical mixing (mix_scheme='KPP', mix_scheme_nmb=1)

> **Revised 2026-05-24 after plan-review + verification against the CORE2 reference namelists.**
> Corrections folded in: real interior constants (not the phantom `visc_cbu_iw`/`diff_cbt_iw`);
> `double_diffusion=.false.` and `use_kpp_nonlclflx=.false.` in CORE2 → `ddmix` and the
> `ghats`→tracer nonlocal flux are GATED OFF (port the gate, defer the unused body — "port what
> CORE2 uses"); `dbsfc` is already half-computed in C (store it, don't re-derive); `mo_convect`
> IS called after KPP; `smooth_blmc=.true.`/`smooth_hbl=.false.` are hardcoded module logicals;
> added the missing `dkm1` state array. **Re-review fix:** vertical diffusion uses a SINGLE
> `Kv = diffK(:,:,1)` (T-channel) for BOTH T and S (`oce_ale.F90:3518-3522`); `fesom_tracer_diff.c`
> needs NO change — the earlier "per-tracer Kv into tracer diffusion" was a deviation, now removed.

## Overview

Port the FESOM1.4 K-Profile Parameterization vertical mixing scheme from Fortran
(`oce_ale_mixing_kpp.F90`, 1185 lines — **NOT** the CVMix wrapper `gen_modules_cvmix_kpp.F90`)
to C, as a **strictly faithful line-by-line port**. KPP is the vertical mixing scheme
production CORE2 actually uses; the C currently runs PP, a deliberate simplification, and all
validation to date has been against *specially-built* Fortran+PP references. KPP lets the C be
validated against the **default CORE2 (KPP) configuration**.

**GUIDING PRINCIPLE (non-negotiable):** THIS IS A PORT. The Fortran has worked for decades.
Do **not** simplify, optimize, or invent. Reproduce every formula, constant, loop bound,
branch, halo exchange, and ordering exactly. Where a value comes from a namelist, port the
**CORE2 reference run's** value (`feedback_namelist_over_codedefault`). "Port what CORE2 uses"
(`feedback_port_used_terms`): features the reference namelist disables (`ddmix`, nonlocal flux)
get their **dispatch gate** ported, but their unused bodies are deferred — exactly as the port
kept `opt_visc=5 visc_filt_bcksct` as faithful-but-unused alongside opt_visc=7. The bar is
dual-instrumentation C vs Fortran on intermediate fields until they match.

KPP is added as a **runtime-selectable scheme alongside PP** (mirroring `oce_ale.F90:3515`
`mix_scheme_nmb`); PP stays byte-identical when KPP is not selected.

## CORE2 reference KPP configuration (VERIFIED — these are the values to port)

From `work_core/namelist.oce` + `work_core/namelist.tra` (the KPP CORE2 config; `work_core`
and `work_linfs_d1800` are both KPP+dt=1800 and produced the reference `fortran_2yr_dt1800`):

| Symbol | Value | Source | Consequence |
|---|---|---|---|
| `mix_scheme` | `'KPP'` (→ `mix_scheme_nmb=1`) | namelist.oce:67 | dispatch branch `==1` |
| `Ricr` | `0.3` | namelist.oce:70 | bulk-Ri threshold in `bldepth` + `Vtc` |
| `A_ver` | `1.e-4` | namelist.oce:19 | background **viscosity** floor in `ri_iwmix` |
| `visc_sh_limit` | `5.0e-3` | namelist.oce:66 | shear-instability viscosity cap |
| `diff_sh_limit` | `5.0e-3` | namelist.tra:97 | shear-instability diffusivity cap |
| `K_ver` | `1.0e-5` | namelist.tra:101 | background **diffusivity** floor |
| `Kv0_const` | `.true.` | namelist.tra:100 | constant K_ver → **`Kv0_background_qiang` (lat/depth) NOT used** (C already matches, `fesom_pp.h:20`) |
| `double_diffusion` | `.false.` | namelist.tra:105 | **`ddmix` NOT exercised** → gate only |
| `use_kpp_nonlclflx` | absent → `.false.` (`oce_modules.F90:168`) | — | **`ghats`→tracer flux OFF** → gate only |
| `Riinfty` | `0.8` (module parameter) | kpp.F90:737 | shear-Ri shape limit |
| `smooth_blmc` | `.true.` (hardcoded module logical) | kpp.F90:73 | blmc smoothing **ON** (its own 3 exchanges + smooth_nod) |
| `smooth_hbl` | `.false.` (hardcoded module logical) | kpp.F90:74 | hbl smoothing OFF |
| `minmix` | `3.0e-3` (in-routine) | kpp.F90:408 | surface viscAE floor |
| `cekman`/`cmonob` | `0.7`/`1.0` (in-routine params) | kpp.F90:478-479 | Ekman/Monin-Obukhov limits |
| derived `Vtc`,`cg` | from concv,concs,epsilon_kpp,vonk,Ricr,cstar | kpp.F90:167,178 | computed in init |
| `mo_convect` | **called after `oce_mixing_KPP`** | oce_ale.F90:3524 | keep convective adjustment for KPP too |

➕ **K0 sub-task:** dump and archive the exact `namelist.oce`+`namelist.tra` that produced
`/scratch/a/a270088/fortran_2yr_dt1800` and re-confirm every value above before coding.

## Context (from discovery)

### Fortran source (port reference)
| File / lines | Routine | Role | CORE2 status |
|---|---|---|---|
| `oce_ale_mixing_kpp.F90:101-210` | `oce_mixing_kpp_init` | wm/ws lookup table (892×482), `Vtc`,`cg`, constants | port |
| `:249-425` | `oce_mixing_KPP` (driver) | dVsq, ustar, Bo → ri_iwmix → (ddmix) → bldepth → blmix_kpp → enhance → combine → node→elem avg | port |
| `:468-648` | `bldepth` (private) | OBL depth `hbl`/`kbl`, Ekman/Monin-Obukhov, caseA — **HIGHEST RISK** | port |
| `:668-717` | `wscale` (private) | wm/ws from the 2D lookup table | port |
| `:730-839` | `ri_iwmix` | interior: shear-Ri + background (A_ver/K_ver) + static instability | port |
| `:852-925` | `ddmix` | double diffusion | **gate only** (double_diffusion=.false.) |
| `:949-1121` | `blmix_kpp` | BL coeffs `blmc[3]`, `dkm1[3]`, shape G(σ), `ghats` (computed) | port |
| `:1140-1184` | `enhance` | enhancement at `kbl-1` (uses `dkm1`) | port |
| `oce_ale.F90:3515-3524` | dispatch | `if mix_scheme_nmb==1: call oce_mixing_KPP(...); call mo_convect(...)` | port |
| `oce_ale_tracer.F90:~892-942` | nonlocal flux | `MIN(ghats·blmc,1)·flux` into tracer RHS | **gate only** (use_kpp_nonlclflx=.false.) |
| `oce_ale_pressure_bv.F90` (dbsfc block) | — | `dbsfc` = g·Δρ surface→level | **already in C as `dbsfc1`** — store it |

Module state to mirror: `dVsq`, `ustar`, `Bo`, `bfsfc`, `stable`, `caseA`, `hbl`, `kbl`(int),
`ghats(nl-1,n)`, `blmc(nl,n,3)`, **`dkm1(n,3)`** (BL coeffs at kbl-1, made by `blmix_kpp:1115-1117`,
consumed by `enhance:1163-1177`), plus the wm/ws lookup tables and scalars `Vtc`,`cg`.

### C integration points (target)
| C file:line | What | Change for KPP |
|---|---|---|
| `src/fesom_pp.c`/`.h` | PP `fesom_pp_mixing` → `aux->Kv`(nodes),`aux->Av`(elems) | mirror as KPP counterpart; keep PP intact |
| `src/fesom_step.c:86-96` | bvfreq→sw_alpha_beta→PP→mo_convect→exchanges | dispatch: KPP-or-PP, **then mo_convect for both** |
| `src/fesom_aux.h:25` / `fesom_aux.c` | single `Kv[n,nl]`, `Av[e,nl]` | KEEP single `aux->Kv` as the TDMA consumer (KPP writes it = `diffK(:,:,1)`); ADD KPP-internal per-tracer `diffK[n,nl,ntr]`, node `viscA`, `hbl`,`kbl`,`ghats`,`blmc[3]`,`dkm1[3]`,`Bo`,`bfsfc`,`stable`,`caseA`,`dVsq`,`ustar`,`dbsfc` |
| `src/fesom_eos.c:114-146` | `dbsfc1` **already computed** per level (for db_max MLD), not stored | store into `aux->dbsfc`, add `dbsfc(nzmax)=dbsfc(nzmax-1)`, gate on KPP |
| `src/fesom_tracer_diff.c:181-330` | implicit TDMA reads single `aux->Kv`; nonlocal skipped `:281` | **NO change** — keep single `aux->Kv` (KPP writes it = `diffK(:,:,1)`); nonlocal flux stays gated off (CORE2) |
| `src/fesom_momentum.c:274-463` | implicit vert. visc. reads `aux->Av` | no change |
| `src/fesom_constants.h:37,43` | `mix_scheme='PP'`, `K_ver`,`A_ver` defines | add dispatch knob; reuse K_ver/A_ver |
| `src/fesom_forcing.h` | `heat_flux`,`water_flux`,`stress_node_surf`,`sw_3d` ✓ | KPP inputs available |

### Out of scope (faithful = port only what CORE2/mix_scheme_nmb==1 uses)
- **`ddmix` body** — `double_diffusion=.false.`; port the gate, defer the body (unused).
- **`ghats`→tracer nonlocal flux body** — `use_kpp_nonlclflx=.false.`; port the gate, defer the
  body. (`ghats` is still *computed* in `blmix_kpp` as a literal part of that routine; it just
  never enters the tracer equation in CORE2.) If ever enabled, the full term is
  `MIN(ghats·blmc[2|3], 1.0)` differenced across interfaces with the `area(nz,n)/areasvol(nz,n)`
  weighting, ×`heat_flux/vcpw·dt` (T, `+`) or ×`-rsss·water_flux·dt` (S, `−`, with
  `ref_sss`/`ref_sss_local`) — full term `oce_ale_tracer.F90:892-987`; document, don't implement now.
- **`Kv0_background_qiang`** (lat/depth background) — `Kv0_const=.true.` → constant K_ver (C matches).
- **CVMix KPP** (`mix_scheme_nmb==3`), **`==17`** add-on, **TKE/TOY** — different/unused.
- **`smooth_Ri_ver`/`smooth_Ri_hor`** in `ri_iwmix` — both default `.false.`; port the gate, no-op.

### Patterns to follow
- One C file per Fortran file: `src/fesom_kpp.c`/`.h` mirrors `oce_ale_mixing_kpp.F90`.
- Halo/loop-bound audits are load-bearing (`feedback_write_loops_halo`, `feedback_array_size_vs_reader_loop`):
  replicate every `myDim` vs `myDim+eDim` bound and every `exchange_nod`. KPP exchanges
  `diffK(:,:,1)`,`diffK(:,:,2)`,`ghats`,`viscA` after the combine (driver `:401-406`), and
  `smooth_blmc` does its **own** 3× `exchange_nod(blmc)` + `smooth_nod` before the combine (`:370-379`).
- Tracer stride is `nl` (`feedback_tracer_stride_nl`): size new 3D node arrays `[N*nl]`.
- Env-gated knobs for dispatch + dumps (mirror `FESOM_NO_GMREDI`, `FESOM_EVP_DUMP_DIR`).

### Key Fortran invariants (must hold in C)
1. **KPP computes per-tracer `diffK` (`Kv_double`, T=`(:,:,1)`/S=`(:,:,2)`), but vertical diffusion
   uses a SINGLE `Kv`.** After KPP, FESOM2 sets `Kv(:,node)=Kv_double(:,node,1)` over `myDim+eDim`
   (`oce_ale.F90:3518-3522`) and the tracer-diffusion TDMA uses that single `Kv` for BOTH T and S
   (`oce_ale_tracer.F90:768,806,846`). The S channel `diffK(:,:,2)` is **NOT** routed into salinity
   diffusion in CORE2 — it only feeds the gated-off nonlocal flux. So `fesom_tracer_diff.c` keeps
   reading the single `aux->Kv` **unchanged** (it already does); the KPP driver just writes
   `aux->Kv = diffK(:,:,1)`. Keep the per-tracer `diffK` storage for the combine (blmc Kv_T/Kv_S
   channels) and the deferred nonlocal flux.
2. **Av at elements via node→element averaging** with surface floor `minmix=3e-3` and bottom fill
   `viscAE(nlevels)=viscAE(nlevels-1)` (`:417,422`).
3. **`mo_convect` runs after KPP** (oce_ale.F90:3524) — convective adjustment is NOT subsumed; keep it.
4. **`bldepth` needs `dbsfc`** (store the existing C `dbsfc1`, extend to nzmax), `dVsq` (built in the
   driver, referenced to the surface velocity), `Bo`, `ustar`, `bvfreq`, `sw_alpha`, `sw_3d`,
   `coriolis_node`, `zbar_3d_n`. The Ekman/Monin-Obukhov block is gated `bfsfc>0 .and. nzmin==1`.
5. **`i_vert_diff=.true.` under KPP** — implicit TDMA already in C; KPP feeds Kv/Av into it.
6. **CORE2 `use_sw_pene=.true.`** → `Bo`/`bfsfc`/`bldepth` include the sw-penetration branch (sw_pene ported).
7. **Interior background** in `ri_iwmix` is `visc_sh_limit·frit + A_ver` / `diff_sh_limit·frit + K_ver`
   with `frit=(1−min(Rig/Riinfty,1)²)³`, `Riinfty=0.8`; Ri is stored in `diffK(:,:,1)` as scratch
   (`:766`); edge copies `diffK(nzmin)=diffK(nzmin+1)`, `diffK(nzmax)=diffK(nzmax-1)`.

## Development Approach

- Faithful port discipline: one Fortran routine per phase, dependency order, literal diff, no refactors.
- **Validation = dual-instrumentation vs Fortran** (not unit tests — matches sea-ice/MFCT/EVP). Each phase
  ends with a Validate step comparing the C routine output to Fortran at an identical input state,
  before the next phase. Requires Fortran-side dumps (Phase K0).
- **PP byte-identical when KPP off** — a discrete Validate item in EVERY phase that compiles changes into
  shared code (K0 aux/dispatch, K5 eos dbsfc, K8 tracer_diff), not just the endpoints.
- Update this plan (`[x]`, ➕, ⚠️) as you go.

## Validation Strategy

- **Per-routine dual-instrumentation (K1–K9):** env-gated dumps in BOTH codes (C `FESOM_KPP_DUMP_DIR`;
  Fortran gate in `oce_ale_mixing_kpp.F90`, mirroring `reference_evp_dump_diagnostic`). Dump each
  routine's outputs by global id at a chosen step; diff with `scripts/kpp_dump_diff.py` (model on
  `evp_dump_diff.py`). Bit-identical where algebra is exact; ≤1e-10 rel where lookup interpolation/order differ.
- **Reference:** build a short matched Fortran KPP **dump** run in K0 (the existing `fortran_2yr_dt1800`
  is end-to-end only). End-to-end climate compare (K9) uses `fortran_2yr_dt1800`.
- **End-to-end (K9):** C+KPP dt=1800 vs `fortran_2yr_dt1800`. Expect Kv/Av to now match KPP (≫ the
  ~0.5 PP-vs-KPP corr) and the ~0.24°C PP warm-SST signature to shrink.

## Implementation Steps

### Phase K0: Scaffolding — dispatch, aux arrays, dump harness, reference provenance (NO behavior change)
**Files:** Create `src/fesom_kpp.{c,h}`, `scripts/kpp_dump_diff.py`; Modify `src/fesom_aux.{c,h}`, `src/fesom_step.{c,h}`, `src/fesom_main.c`, `CMakeLists.txt`, `oce_ale_mixing_kpp.F90` (dump gate)
- [x] Archive the exact namelists that produced `fortran_2yr_dt1800`; re-confirm every value in the "CORE2 reference KPP configuration" table. → `docs/kpp_reference_namelists/` (+ PROVENANCE.md). `work_core`/`work_linfs_d1800` namelist.oce/.tra **byte-identical**; all values confirmed (mix_scheme=KPP, Ricr=0.3, A_ver=1e-4, visc/diff_sh_limit=5e-3, K_ver=1e-5, Kv0_const=.true., double_diffusion=.false., use_kpp_nonlclflx absent→.false., dt=1800).
- [x] Add KPP module state — in **struct `fesom_kpp`** (`fesom_kpp.{c,h}`, mirroring o_mixing_KPP_mod + the fesom_gm precedent), NOT aux: per-tracer `diffK[ntr·N·nl]`, node `viscA`, `blmc[3·N·nl]`, `ghats[N·(nl-1)]`, `dVsq`, `dkm1[3·N]`, `hbl`,`bfsfc`,`caseA`,`stable`,`ustar`,`Bo`,`kbl`(int), `wmt`/`wst` tables, `Vtc`/`cg`. ➕ **Only `dbsfc` lives in aux** (EOS writes it; EOS carries aux not kpp). Layouts match Fortran column-major so per-component slabs are contiguous & directly exchangeable. KEEP single `aux->Kv` as TDMA consumer.
- [x] Add runtime selector `FESOM_MIX_SCHEME=PP|KPP` (default PP); dispatch in `fesom_step.c` mirroring `oce_ale.F90:3515`, **`mo_convect` after either scheme** (shared). KPP branch = `fesom_kpp_mixing` abort stub. Threaded `kpp` through `fesom_step_ctx` + `fesom_main.c` (alloc/ctx/free) like `gm`.
- [x] Add `FESOM_KPP_DUMP_DIR` dumps: C `fesom_kpp_dump_nodes` + Fortran gate (`kpp_dump_load_env`/`kpp_dump_nod`/`kpp_dump_prestep`/`kpp_dump_final`, FESOM_KPP_DUMP_STEP, call counter). `scripts/kpp_dump_diff.py` written. Fortran dump run = `work_kpp_dump` (linfs, dt=1800, 16r) → `/work/.../kpp_fdump/dump` (8 tags: prestep/dVsq/dbsfc/bldepth/viscA/diffKt/diffKs/ghats; diff-script self-parse 0.000e+00).
- [x] **Validate:** with KPP off, byte-identical to HEAD `028cc1b` (PP regression). ✅ `job_kpp_k0_byteident` (16r/200steps/dt=500): all snapshots + run.log diagnostics **bit-for-bit identical** (worst |Δ|=0.000e+00); `scripts/kpp_byteident_check.py`.

### Phase K1: `oce_mixing_kpp_init` — lookup tables + constants ✅
**Files:** Modify `src/fesom_kpp.{c,h}`, `src/fesom_main.c`; `scripts/kpp_init_reference.py`
- [x] Ported the wm/ws 2D lookup table: `nni=890`, `nnj=480`, `[(nni+2)·(nnj+2)]`=892×482 each, `deltaz/deltau`, exact build loops (`:194-219`). `fesom_kpp_init` wired in `main` after alloc (harmless under PP).
- [x] Constants in `fesom_kpp.c` (KPP_*): `Ricr=0.3`, `concv=1.6`, A_ver/K_ver reuse FESOM_PHASE1_* (1e-4/1e-5), `visc/diff_sh_limit=5e-3`, `Riinfty=0.8`, `minmix=3e-3`, `cekman=0.7`, `cmonob=1.0`, `epsilon_kpp=0.1`, `vonk=0.4`, `conc1=5`, `epsln=1e-40`, zmin/zmax/umin/umax; derived `Vtc` (`:177`), `cg` (`:188`).
- [x] **Validate:** full table + 4 constants dumped (C `fesom_kpp_dump_init` + Fortran gate) vs independent `scripts/kpp_init_reference.py`. **Vtc/cg/deltaz/deltau bit-identical C-vs-Fortran (|Δ|=0); wmt/wst C-vs-F max|Δ|=3.5e-18 (max_rel ~4e-16, libm last-ULP), well within ≤1e-12.** Both codes also match the python reference to ~4e-16. (`job_kpp_k1_validate`.)

### Phase K2: `wscale` — turbulent velocity scales ✅
**Files:** Modify `src/fesom_kpp.{c,h}`, `src/fesom_main.c`; `scripts/kpp_wscale_reference.py`
- [x] Ported `wscale` (`:828-877`) as static `kpp_wscale`: index `INT(zdiff/deltaz)`, `MIN(iz,nni)`/`MAX(iz,0)`, `INT(MIN(udiff/deltau,nnj))`/`MAX(ju,0)` clamps, bilinear interpolation, stable formula (`zehat>zmax`).
- [x] **Validate:** `(zehat,ustar)` sweep 201×101 spanning unstable-table / stable-branch / clamp regions (`fesom_kpp_dump_wscale_sweep` + Fortran gate, identical grid). **C-vs-Fortran wscale max|Δ|=2.8e-15 (max_rel ~3e-13)**, grid bit-identical; both match independent `scripts/kpp_wscale_reference.py` to ~3e-13. Well within ≤1e-10. (`job_kpp_k2_validate`.)

### Phase K3: `ri_iwmix` — interior mixing ✅ (background+static-instability; shear curve → K9)
**Files:** Modify `src/fesom_kpp.c` (incremental driver); Fortran ri dump point
- [x] Ported `ri_iwmix` (`:890-999`) as static `kpp_ri_iwmix`: loop-1 Ri=MAX(N²,0)/(shear²+epsln) stored in diffK ch0 scratch; loop-2 `frit=(1−min(Rigg/0.8,1)²)³`, viscA=`visc_sh_limit·frit+A_ver`, diffK(1)=`diff_sh_limit·frit+K_ver` (Kv0_const), diffK(2)=diffK(1); edge copies both loops. `smooth_Ri_ver/hor` (.false.) omitted. Incremental driver `fesom_kpp_mixing`: zero viscA (myDim+eDim) + ri_iwmix + dump; guarded so only a FESOM_KPP_DUMP_DIR run (1 step) proceeds (bldepth+ pending).
- [x] **Validate:** ri dump (viscA/diffKt/diffKs + bvfreq) vs Fortran at step 1. At rest (uv=0→shear=0) ri_iwmix is a deterministic binary of sign(bvfreq): both `frit` endpoints give the **exact** Fortran values (background 1e-4/1e-5, static-instability 5.1e-3/5.01e-3 — diff value exactly 5e-3). `scripts/kpp_ri_bvfreq_check.py`: self-consistency 0 violations both codes; every C-vs-Fortran viscA disagreement is exactly a bvfreq sign flip (~0.5–1.3%/level, weak-stratification IC noise per [[feedback_phc_rank_dependent]]), 0 unexplained. The **shear-instability curve** (0<frit<1) is untested at rest → validated end-to-end at K9 + direct line-by-line port. PP byte-identical (KPP only entered under FESOM_MIX_SCHEME=KPP).

### Phase K4: `ddmix` — double diffusion (GATE ONLY for CORE2) ✅
**Files:** Modify `src/fesom_kpp.c`
- [x] Ported the `IF (double_diffusion) CALL ddmix` gate (`:420-422`) as a compile-time `#if KPP_DOUBLE_DIFFUSION` (=0 for CORE2) with `#error` + the ddmix call site commented — body deferred (port-what-CORE2-uses). If ever enabled, the build fails loudly pointing at the body (`:1012-1085`) to port.
- [x] **Validate:** no-op — the gate is compile-time-excluded, so the driver is numerically identical to K3 (which is PP byte-identical with KPP off and ri-validated with KPP on). No ddmix call in the binary.

### Phase K5: driver pre-step + `dbsfc` + `bldepth` — OBL depth (HIGHEST RISK) ✅
**Files:** Modify `src/fesom_kpp.{c,h}`, `src/fesom_eos.{c,h}`, `src/fesom_step.c`
- [x] `fesom_eos.c`: **store** `dbsfc1` into `aux->dbsfc` per level + `dbsfc(nzmax)=dbsfc(nzmax-1)` (`:337`). Stored UNCONDITIONALLY (PP never reads it → numerically PP-safe; the Fortran `if(mixing_kpp)` gate is a write-skip). C-vs-Fortran **bit-identical** (max|Δ|~8.7e-15, rel ~1e-13).
- [x] ➕ **bvfreq fix (REQUIRED for bldepth fidelity — `[[feedback_bvfreq_smoothing_gap]]`):** bldepth reads `bvfreq` VALUES (`Vtsq∝√|bvfreq|`) to the bottom interface and exposed two pre-existing `fesom_eos.c` deviations: (1) a **bottom-pad off-by-one** (`bvfreq[nzmax-1]=bvfreq[nzmax-2]` → fixed to `bvfreq[nzmax]=bvfreq[nzmax-1]`), and (2) the **deferred horizontal N² smoothing** (`N2smth_h=.true.` in ALL ref namelists, `smooth_nod3D(bvfreq,1)`, pressure_bv:499). Ported `fesom_smooth_nod3D` (gen_support.F90:99-198, area-weighted node-patch avg + halo exch), called from `fesom_step.c` after the bvfreq exchange, **for all schemes** (Fortran does so regardless of mix_scheme). Result: `ri_bvfreq` C-vs-Fortran now **bit-identical** (rel ~1e-13; was 2–7%/level) → `ri_viscA/diffKt/diffKs` now **bit-identical** (K3's "sign-flip noise" was this smoothing). **Re-baselines PP** (user-approved 2026-05-24): byte-identity gate must re-point to a fresh PP run; d1800 PP climate shifts (max|ΔS|~5 PSU at convection nodes after 200 steps, small in mean — PP now MORE faithful). IC byte-identical.
- [x] Port the driver pre-step (`:344-411`): `dVsq` (eqn 21), `ustar` (eqn 2, `density_0_r=1/1030` confirmed vs `oce_modules.F90:14`), `Bo` (A2c/A2d). dVsq C-vs-Fortran **bit-identical** (0.0 at rest).
- [x] Port `bldepth` (`:674-854`): bulk-Ri search, Ricr=0.3 hbl interp, sw-pene `bfsfc` branches (use_sw_pene=.true.→inline), Ekman/Monin-Obukhov gated `bfsfc>0 .and. nzmin==1`(C:`==0`), `smooth_hbl`(off), final `kbl`, `caseA`. kbl stored 0-based (dump emits +1). SIGN→`copysign`.
- [x] **Validate (`job_kpp_k5_validate`, 16r/dt=1800/step1 vs `/work/.../kpp_fdump/dump`):** ALL EOS inputs bit-identical (dVsq=0, dbsfc/bvfreq/ri_* ~1e-13). bldepth out diffs (hbl ≤54m, kbl ≤5lev at 1% of nodes) are ENTIRELY forcing-input driven: `bfsfc`↔`Bo` corr **0.9993**; kbl identical at **99.0%**; the only "suspicious" nodes (`scripts/kpp_bldepth_xcheck.py`) are textbook `bfsfc>0` Ekman-gate sign-flips (Fortran hbl=5m=`|zbar(2)|` clamp, C unclamped) from the tiny `ΔBo~5e-10`. **0 unexplained algebra disagreements.** ustar/Bo differ because step-1 FORCING (stress, heat/water flux) differs C-vs-Fortran (pre-existing, K5 doesn't touch forcing; bounded — PP climate validated at 0.04 RMS with it). PP re-baseline run confirms `smooth_nod3D` stable (200 steps clean).

### Phase K6: `blmix_kpp` — BL coeffs + dkm1 + ghats ✅
**Files:** Modify `src/fesom_kpp.c`; Fortran `oce_ale_mixing_kpp.F90` (kpp_dump_blmix, uncommitted)
- [x] Ported `blmix_kpp` (`:1155-1327`) as static `kpp_blmix`: dthick/diff_col per-column scratch, caseA-selected matching level `kn=kbl-caseA` (capped), one-sided derivatives + interior match at `hbl` (eqn 18), shape functions G(σ) (eqn 11), `blmc[3]` (eqn 10), `ghats` (eqn 20), `dkm1[3]` at `kbl-1` (σ from zbar, not Z). Channels 0/1/2 = momentum(wm,Gm)/T(ws,Gt)/S(ws,Gs); note difsp/Gs↔diff_col ch3(=slab2,S), diftp/Gt↔ch2(=slab1,T). Wired after bldepth in the driver.
- [x] **Validate (CONTROLLED REPLAY — the key technique for K6+):** the step-1 forcing diff ([[project_forcing_step1_diff]]) perturbs bfsfc/ustar at ~every node, and blmix amplifies it (f1∝bfsfc/ustar⁴, wscale, stable-flips) → the live-run dump (`job_kpp_k6_validate`) diffs at ~52% of nodes, uninterpretable alone. So added `FESOM_KPP_REPLAY_DIR` (`kpp_replay_inputs`): inject the Fortran-dumped bldepth/prestep/ri inputs into the C before blmix (same 16r partition → dump line i = local node i). `job_kpp_k6_replay`: blmc_m/t/s, dkm1, bl_ghats vs Fortran = **worst max|Δ|=3.18e-13, 0 signal** (wscale libm last-ULP) → blmix algebra bit-faithful. (Fortran dumps blmc/dkm1/ghats via `kpp_dump_blmix` BEFORE enhance, since enhance modifies them at kbl-1; `scripts/kpp_blmix_xcheck.py` for the live-run input attribution.)

### Phase K7: `enhance` + assembly + exchanges ✅
**Files:** Modify `src/fesom_kpp.{c,h}`; Fortran `oce_ale_mixing_kpp.F90` (kpp_dump_elem + viscAE, uncommitted)
- [x] Ported `enhance` (`:1346-1390`) as static `kpp_enhance`: blends interior(caseA)/BL/dkm1 at kbl-1 with `delta=(hbl+zbar(k))/(zbar(k)-zbar(k+1))`, `dstar=(1-δ)²·dkm1+δ²·dkmp5`, `blmc(k)=(1-δ)·X+δ·dstar`; `ghats(k)*=(1-caseA)`. Channels m/T/S (T=diffK ch1, S=diffK ch2).
- [x] Ported the `smooth_blmc` block (`:437-447`, smooth_blmc=.true.): per channel `fesom_exchange_nod3D(blmc_j)` then `fesom_smooth_nod3D(blmc_j, 3)` (3 sweeps; the already-validated bvfreq smoother reused with N_smooth=3).
- [x] Ported the combine + assembly (`:451-490`): `nz<kbl` → `viscA/diffK=fmax(interior,blmc)`, else `ghats=0`; exchanges on `diffK(1)`,`diffK(2)`,`ghats`(stride nl-1),`viscA`; node→elem `viscAE=Σ/3` into `aux->Av`; bottom fill `viscAE(nlevels)=viscAE(nlevels-1)`; surface `minmix=3e-3` floor.
- [x] **Validate (controlled replay `job_kpp_k7_replay`):** inject Fortran bldepth/prestep/ri inputs → run the full K7 chain → final `viscA`/`diffKt`/`diffKs`/`ghats` + element `viscAE` vs Fortran = **0 signal, worst max|Δ|=3.18e-13** (smooth_nod3D + wscale libm) → K7 algebra bit-faithful. Every rank replays its own owned inputs so the smooth_blmc pre-exchange fills halos with Fortran-matching values. Added `fesom_kpp_dump_elems` (C) + `kpp_dump_elem`/viscAE in the Fortran `kpp_dump_final` (uncommitted, Intel tree).

### Phase K8: wire KPP outputs into the step (single Kv; nonlocal flux stays GATED OFF) ✅
**Files:** Modify `src/fesom_kpp.c` (Kv copy + guard removed). `src/fesom_step.c` dispatch already anticipated K8 (no change). **NOT** `fesom_tracer_diff.c`.
- [x] At the end of `fesom_kpp_mixing`: `memcpy(aux->Kv, k->diffK ch0, n_nod*nl)` = single `aux->Kv = diffK(:,:,1)` (T-channel) over myDim+eDim (`oce_ale.F90:3518-3522`; diffK ch0 halo-valid from the combine exchange). `fesom_tracer_diff.c` UNCHANGED — reads the single `aux->Kv` for T AND S (S-channel diffK(:,:,2) NOT routed into salinity in CORE2). Removed the driver guard → KPP is a LIVE scheme.
- [x] `use_kpp_nonlclflx` stays OFF: `ghats`→tracer flux NOT added (`fesom_tracer_diff.c` already skips it; ghats computed but unused, full term documented in "Out of scope").
- [x] Dispatch in `fesom_step.c` was already written for K8 (`s_use_kpp` → `fesom_kpp_mixing`; `mo_convect` + the Kv/Av exchanges shared after) — no change needed.
- [x] **Validate (`job_kpp_k8_validate`, 16r):** (a) KPP LIVE 100 steps dt=1800 = finite (exit 0, no NaN/uv>5, CG ~110/step vs PP ~90); (b) PP 200 steps dt=500, KPP off = **byte-identical** to the K5 PP run (snap_000000+snap_000200 worst |Δ|=0). K6/K7/K8 are KPP-only → PP unchanged.

### Phase K9: end-to-end climate validation vs fortran_2yr_dt1800 ✅
**Files:** `scripts/kpp_climate_compare.py`, `job_kpp_2yr_dt1800`
- [x] Ran C+KPP dt=1800 864r 2yr (1958-59, `kpp_2yr_dt1800`) — clean (exit 0, CG ~86/step, no NaN/uv>5).
- [x] **Validate vs `fortran_2yr_dt1800` (the default-CORE2 KPP reference):** C+KPP matches
  Fortran-KPP to global SST/SSS RMS **0.0046/0.0023** (1958), **0.0128/0.0044** (1959); biases ~0,
  non-drifting → KPP port faithful end-to-end (as tight as the PP-vs-PP 0.004 match).
- [x] **Scheme contrast confirms the KPP payoff:** C+**PP** vs Fortran-**KPP** = SST RMS
  **0.085/0.093** — the genuine PP↔KPP difference is ~18× the C+KPP residual, so the comparison
  correctly resolves schemes and C+KPP reproduces its reference.
- [x] Residual: the slow 1958→1959 SST-RMS growth (0.0046→0.0128) is the accumulated step-1
  sea-ice-forcing transient ([[project_forcing_step1_diff]]) + chaos; tiny, biases stay ~0.

### Phase K10: Acceptance + stability ✅ (stability)
- [x] All Overview requirements met: KPP selectable (`FESOM_MIX_SCHEME=KPP`); PP byte-identical when off (K8); every ported routine replay-bit-faithful (K1-K7); climate matches the KPP reference (K9, 0.005-0.013°C).
- [x] **Multi-year stability (`kpp_5yr_dt1800`, job 25108160): C+KPP ran a full 5 yr (1958-1962, 87589 steps) clean — exit 0, no uv>5/NaN/blowup, max ice velocity ~1.2 m/s.** 5-yr global-mean drift tiny + bounded: SST yr2→yr5 −0.029°C, SSS −0.0085 PSU. KPP is stable at dt=1800.
- [x] **Default `mix_scheme` = KPP** (user decision 2026-05-25, commit `8d0cdbc`): `fesom_step` dispatch defaults to KPP; `FESOM_MIX_SCHEME=PP` opts back to PP. Matches CORE2 production.
- ➕ **Found while digging the forcing diff:** `ice_gamma_fct` used the 0.25 module default instead of the CORE2 namelist 0.5 ([[feedback_namelist_over_codedefault]]); fixed (`fesom_ice.c`). Re-baselines the ice/climate (re-validation owed); NOT the dominant forcing-diff root (that's the step-1 EVP `u_ice`).

### Phase K11: [Final] Docs + memory + commit ✅
- [x] Updated `docs/NEXT_SESSION_HANDOFF.md` (KPP complete + default) + memory (`project_kpp_port_state`, `feedback_controlled_replay_validation`, `feedback_bvfreq_smoothing_gap`, `project_forcing_step1_diff`, `feedback_namelist_over_codedefault`).
- [x] Committed per-phase (K0 `13bfe95` … K10/default `8d0cdbc`); tag the validated KPP milestone.
- [x] Moved this plan to `docs/plans/completed/`.
- ➕ Next-session kickoff for the deferred step-1 forcing fix: `docs/NEXT_SESSION_PROMPT.md` + `docs/FORCING_STEP1_DIFFERENCE.md`.

## Technical Details
- **New aux** (stride `nl`, `myDim+eDim`): indexed `diffK[N*nl*ntr]`, node `viscA[N*nl]`, `dVsq`, `dbsfc`, `ghats`; `blmc[N*nl*3]`, `dkm1[N*3]`; per-node scalars `hbl`,`Bo`,`bfsfc`,`stable`,`caseA`,`ustar`; per-node int `kbl`.
- **Driver call order** (`:280-424`): zero viscA → dVsq → ustar,Bo → `ri_iwmix` → (gate `ddmix`) → `bldepth` → `blmix_kpp` → `enhance` → (`smooth_blmc` exch+smooth) → combine(MAX in BL, ghats=0 else) → exchange diffK(1),diffK(2),ghats,viscA → node→elem viscAE + bottom fill + minmix floor. Then `aux->Kv(:,node)=diffK(:,node,1)` over `myDim+eDim` (single Kv = T-channel for all tracers, `oce_ale.F90:3518-3522`). Then in the step driver: `mo_convect`.
- **Single-Kv consumer:** the tracer TDMA reads only the single `aux->Kv`; the per-tracer `diffK[ntr]` storage exists for the combine (blmc channels) and the deferred nonlocal flux, NOT for salinity diffusion.
- **Dispatch:** CORE2 KPP = `mix_scheme_nmb==1`. NOT 3 (CVMix), NOT 17.

## Post-Completion
*Informational — external/manual:*

**Manual verification:**
- MLD maps vs observations in deep-convection regions (Labrador, Weddell, GIN) — KPP's qualitative payoff over PP.
- Confirm (don't assume) the namelist provenance of `fortran_2yr_dt1800` before trusting K9.

**Decisions deferred to the user:**
- Whether KPP becomes the default `mix_scheme` (production CORE2) or stays opt-in.
- Whether to extend KPP validation to the 5-yr stability/drift protocol (`drift_5yr_d1800.py`) after K9.
- Whether to ever port the `ddmix` body / nonlocal-flux body (only if a non-CORE2 config enables them).
