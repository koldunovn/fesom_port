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

### Phase K4: `ddmix` — double diffusion (GATE ONLY for CORE2)
**Files:** Modify `src/fesom_kpp.c`
- [ ] Port the `IF (double_diffusion) CALL ddmix` gate (`:359`). Since CORE2 has it `.false.`, the body is **deferred (unused)** — add a clear stub/abort if the flag is ever set, mirroring the kept-but-unused pattern. (Port the body later only if a config needs it.)
- [ ] **Validate:** confirm no-op (byte-identical with `double_diffusion=.false.`).

### Phase K5: driver pre-step + `dbsfc` + `bldepth` — OBL depth (HIGHEST RISK)
**Files:** Modify `src/fesom_kpp.c`, `src/fesom_eos.c`
- [ ] In `fesom_eos.c:114-146`, **store** the already-computed `dbsfc1` into `aux->dbsfc[N*nl]`, add `dbsfc(nzmax)=dbsfc(nzmax-1)` (`:337`), gate on KPP so the PP path is byte-identical.
- [ ] Port the driver pre-step (`:280-352`): `dVsq` (shear vs surface velocity, eqn 21), `ustar` (from `stress_node_surf`, eqn 2), `Bo` (eqns A2c/A2d from `sw_alpha/beta`, `heat_flux`, `water_flux`, surface S).
- [ ] Port `bldepth` (`:468-648`): bulk-Ri search (`Ritop=zk·dbsfc`, `Vtsq=zk·ws·√|bvsq|·Vtc`, `Rib=Ritop/(dVsq+Vtsq+epsln)`), Ricr=0.3 linear interpolation of `hbl`, sw-pene `bfsfc` branch, Ekman/Monin-Obukhov limits gated `bfsfc>0 .and. nzmin==1` (`AMIN1/AMAX1`→`fmin/fmax`), `smooth_hbl` (off), final `kbl`, `caseA`.
- [ ] **Validate (most effort here):** diff `dbsfc`,`dVsq`,`ustar`,`Bo`,`bfsfc`,`hbl`,`kbl`,`caseA` node-by-node; bin by region; chase any `hbl`/`kbl` disagreement (historical bug: wrong OBL → surface flux concentrated, `FRESH_START:489`). Re-run the PP regression (eos touched).

### Phase K6: `blmix_kpp` — BL coeffs + dkm1 + ghats
**Files:** Modify `src/fesom_kpp.c`
- [ ] Port `blmix_kpp` (`:949-1121`): shape functions G(σ), `blmc[3]` (viscA/Kv_T/Kv_S BL coeffs) matched to interior at `hbl` via `caseA`, `dkm1[3]` at `kbl-1`, and `ghats` (computed even though CORE2 won't consume it downstream).
- [ ] **Validate:** diff `blmc[1..3]`, `dkm1[1..3]`, `ghats` node-by-node.

### Phase K7: `enhance` + assembly + exchanges
**Files:** Modify `src/fesom_kpp.c`
- [ ] Port `enhance` (`:1140-1184`) using `dkm1` (enhancement at `kbl-1`).
- [ ] Port the `smooth_blmc` block (ON): 3× `exchange_nod(blmc(:,:,j))` + `smooth_nod(blmc,3)` (`:370-379`).
- [ ] Port the combine (`:382-396`): for `nz<kbl` `viscA/diffK=MAX(interior,blmc)`, else `ghats=0`; then `exchange_nod` on `diffK(1)`,`diffK(2)`,`ghats`,`viscA`; node→elem `viscAE=Σ/3`; bottom fill `viscAE(nlevels)=viscAE(nlevels-1)`; surface floor `minmix=3e-3`.
- [ ] **Validate:** diff final `viscAE`(elements), `diffK(T)`, `diffK(S)`, `ghats` full fields — the module-level gate.

### Phase K8: wire KPP outputs into the step (single Kv; nonlocal flux stays GATED OFF)
**Files:** Modify `src/fesom_step.c` (dispatch); `src/fesom_kpp.c` (Kv copy). **NOT** `fesom_tracer_diff.c`.
- [ ] After the KPP combine, set the single `aux->Kv(:,node) = diffK(:,node,1)` (T-channel) over `myDim+eDim`, mirroring `oce_ale.F90:3518-3522`. **Do NOT** change `fesom_tracer_diff.c` — it keeps reading the single `aux->Kv` for all tracers (faithful: S uses the T-channel Kv in CORE2). The per-tracer `diffK` storage remains for the combine + deferred nonlocal flux only.
- [ ] Keep `use_kpp_nonlclflx` gated **off** (CORE2): the `ghats`→tracer flux is NOT added. Leave a documented, flag-gated insertion point (full term in "Out of scope") for a future config.
- [ ] Wire the KPP branch into `fesom_step.c` dispatch; keep `mo_convect` after KPP (both schemes).
- [ ] **Validate:** PP byte-identical when KPP off; with KPP on, a few-step run finite; dump-compare the **post-copy single `aux->Kv`** (the array the TDMA actually reads) against the Fortran module `Kv` after `oce_ale.F90:3520` — NOT `Kv_double(:,:,2)` — plus `Av`/`viscAE`. (Routine-level K6/K7 dumps still compare both `diffK(:,:,1)` and `(:,:,2)`.)

### Phase K9: end-to-end climate validation vs fortran_2yr_dt1800
**Files:** Create `scripts/kpp_climate_compare.py` (model on `clim_validation_2yr.py`), `job_kpp_2yr_dt1800`
- [ ] Run C+KPP dt=1800 864r 2yr (1958-59); compare to `/scratch/a/a270088/fortran_2yr_dt1800`.
- [ ] Confirm Kv/Av now correlate with KPP, MLD/convection regions match, SST bias shrinks vs PP.
- [ ] Localize any divergence via the K1–K7 dumps.
- [ ] **Validate:** SST/SSS area-weighted RMS + maps vs the KPP reference; document residuals.

### Phase K10: Acceptance + stability
- [ ] All Overview requirements met; KPP selectable; PP byte-identical when off; every ported routine dual-instrument-matched; climate matches the KPP reference.
- [ ] Multi-year stability at dt=1800 with KPP (stronger mixing; watch the autumn transient, `project_5yr_milestone`).
- [ ] Decide default `mix_scheme` (PP vs KPP) — separate, after validation.

### Phase K11: [Final] Docs + memory + commit
- [ ] Update `docs/NEXT_SESSION_HANDOFF.md`, memory (`project_dt1800_state`, `feedback_port_used_terms`, new `project_kpp_port_state`).
- [ ] Commit per-phase; tag if it becomes a validated milestone.
- [ ] Move this plan to `docs/plans/completed/`.

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
