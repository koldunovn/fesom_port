# Step-1 forcing difference C-vs-Fortran (ustar / Bo) — ADDRESS IMMEDIATELY AFTER KPP

**Status:** OPEN. Deferred by user decision (2026-05-24) to *immediately after KPP is fully
ported* (after K9/K11). Recorded here so the investigation can start from full context.

## What

At ocean step 1 (cold start from PHC IC, JRA55 1958, dt=1800, 16r), the **surface forcing
fields differ C-vs-Fortran**, surfaced by the KPP K5 `bldepth` dual-instrumentation
(`job_kpp_k5_validate`, dumps in `/work/ab0995/a270088/port/kpp_k5_validate/kpp_dump` vs the
Fortran reference `/work/ab0995/a270088/port/kpp_fdump/dump`). The two KPP pre-step quantities
derived from forcing both diverge:

| field | formula | depends on | C-vs-Fortran @ step1 |
|---|---|---|---|
| `ustar` | `sqrt(sqrt(τx²+τy²)·ρ0⁻¹)` | `stress_node_surf` (wind stress) | max\|Δ\|=3.006e-3, max\|val\|=3.483e-2, **rel 8.6%**, #(\|Δ\|>1e-9)=**37004 / 126858 (29%)** |
| `Bo` | `−g·(α·Q/vcpw + β·W·S)` | `heat_flux` Q, `water_flux` W | max\|Δ\|=2.522e-6, max\|val\|=4.049e-6, **rel 62%** (tiny abs), #(\|Δ\|>1e-9)=**123899 / 126858 (98%)** |

The **non-forcing** ingredients of these formulas are all **bit-identical** C-vs-Fortran
(verified in the same dump run): `sw_alpha`/`sw_beta` (EOS, rel ~1e-13), surface `S`, `ρ0`
(=1030, `density_0_r=1/1030` confirmed vs `oce_modules.F90:14`), `vcpw`, `g`. So:

- **`ustar` diverges ⇔ `stress_node_surf` (wind stress) diverges** — at ~29% of nodes.
  It is **bit-identical at the other ~71%** (e.g. all 4 "clean" bldepth nodes had Δustar ~1e-17).
- **`Bo` diverges ⇔ `heat_flux` and/or `water_flux` diverge** — at ~98% of nodes, by small
  absolute amounts (Bo itself is O(1e-6) m²/s³, so even tiny flux diffs give large *relative* Bo
  diffs and can flip its sign).

## Consequence (why it matters)

1. **bldepth (K5):** the differing ustar/Bo *threshold-flip* the OBL diagnosis at ~1% of nodes —
   `bfsfc↔Bo` correlate at 0.9993, `kbl` is identical at 99.0%, and the residual `hbl` diffs
   (≤54 m) are all bulk-Ri-exit or `bfsfc>0` Ekman-gate flips (Fortran `hbl=|zbar(2)|` clamp fires,
   C's doesn't). The bldepth **algebra is faithful** (`scripts/kpp_bldepth_xcheck.py`: 0
   unexplained); the diffs are pure forcing-input noise.
2. **K6/K7 per-routine validation:** every downstream KPP routine validated at step 1 inherits this
   input noise → validate the **algebra via input-cross-check**, not bit-identity of threshold
   outputs.
3. **K9 end-to-end KPP climate:** forcing diffs propagate, so C-KPP will be climate-*close* to
   Fortran-KPP (like PP), not bit-identical. PP climate validated at 0.040/0.016 SST/SSS RMS WITH
   this same forcing noise ([[project_climate_validation_d1800]]), so it is bounded — but it caps
   how tight K9 can get and should be removed for a clean KPP climate match.

## Scope / provenance

**Pre-existing.** K5 touches no forcing code (only eos `dbsfc`/`bvfreq` + KPP). The forcing path
was ported in earlier phases and accepted at the *climate* level (RMS ~0.04), never verified
bit-identical at step 1. So this is a long-standing residual now made visible by the first
consumer that reads the forcing-derived `ustar`/`Bo` per node at step 1.

## Hypotheses (ranked) + investigation plan

The two fields fail differently, suggesting two causes:

- **`stress_node_surf` differs at ~29% (region/state-dependent), bit-identical elsewhere** →
  most likely the **sea-ice ocean-surface stress**: ice-covered nodes use ice–ocean stress (a
  function of the ice state + ice velocity), open-water nodes use atm–ocean stress. A small step-1
  ice-state or ice-velocity difference (EVP, a_ice/m_ice IC interp) would perturb stress exactly on
  the ice-covered subset. **Test:** map the 37004 differing gids → lat/ice-cover (`a_ice` at step1);
  expect them concentrated under ice. Compare `stress_surf`/`stress_node_surf` and `a_ice`/`uice`
  C-vs-Fortran at step 1. Files: `fesom_ice*.c` (EVP, ice–ocean stress), `fesom_bulk.c` (atm–ocean τ).
- **`heat_flux`/`water_flux` differ ~everywhere by small amounts** → a systematic difference in the
  **bulk-flux assembly** (L&Y09): a longwave/sensible/latent component, the SST/humidity/air-temp
  inputs, the transfer-coefficient iteration, or runoff/precip in `water_flux`. **Test:** dump
  `heat_flux`, `water_flux` (and their components: SW/LW/sensible/latent) from BOTH codes at step 1
  and diff component-by-component to localize. Files: `fesom_bulk.c` (`fesom_bulk_compute`,
  `Ch_atm_oce`/`Ce_atm_oce`), `fesom_forcing.c` (JRA interp), `fesom_ice_thermodynamics`
  (heat_flux assembly incl. the runoff fold, [[feedback_runoff_fold_when_ice]]).

**To reproduce / extend the evidence:** the K5 dump already has C `ustar`/`Bo`; add Fortran+C dumps
of `heat_flux`/`water_flux`/`stress_node_surf` (+ components) at step 1, 16r, and diff with
`scripts/kpp_dump_diff.py`. Bin the ustar-differing gids by latitude / ice cover.

## Pointers

- Evidence run: `job_kpp_k5_validate` → `/work/ab0995/a270088/port/kpp_k5_validate/`.
- Cross-check script: `scripts/kpp_bldepth_xcheck.py` (input↔output attribution).
- Related memory: [[feedback_bvfreq_smoothing_gap]] (the *other* gap K5 surfaced, now FIXED),
  [[feedback_phc_rank_dependent]], [[feedback_runoff_fold_when_ice]],
  [[feedback_array_size_vs_reader_loop]] (a prior forcing-halo OOB lesson).
