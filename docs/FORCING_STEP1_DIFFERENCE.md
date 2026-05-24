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

## Dig results (2026-05-25) — CONFIRMED: it is the sea-ice ocean-surface stress

**Geographic test (`/tmp/forcing_geo.py`, C kpp_k6_validate vs Fortran kpp_fdump, gid→lat via the
mesh diag) settles it.** ustar diffs are concentrated in the polar/ice regions and ZERO (bit-
identical, max|Δ|~1e-15) in ice-free latitudes:

| band | frac(\|Δustar\|>1e-9) | max\|Δustar\| |
|---|---|---|
| Antarctic <−60 | 77% | 3.0e-3 |
| Arctic >60 | 78% | 3.0e-3 |
| N mid 40–60 (ice edge) | 11% | 9.9e-4 |
| Tropics / subtropics / S-mid 40 | **~0%** | **~1e-15 (bit-identical)** |

Large ustar diffs (>1e-4): **100% at |lat|>50, 95% at |lat|>60**, median |lat|=74. Bo follows the
same pattern (polar max ~2.5e-6 vs tropical ~2.5e-9). So **ice-free forcing is bit-identical; the
whole difference is the sea-ice contribution** (over-ice stress + over-ice heat/water flux).

**Mechanism + narrowing (code-verified).** `stress_node_surf` over ice (`fesom_ice_coupling.c:223-234`,
matching Fortran `ice_oce_coupling.F90:112-120` line-for-line):
`sns = stress_iceoce·a_ice + stress_atmoce·(1−a_ice)`, with
`stress_iceoce = ρ·cd_oce_ice·|u_ice−u_w|·(u_ice−u_w)`. The atm-ocean term is bit-identical
(ice-free ustar proves it). At **step 1 the ocean is at rest → `u_w` (=srfoce) = 0**, so the
over-ice stress collapses to `ρ·cd·|u_ice|·u_ice·a_ice`. Hence the difference is the **step-1 EVP
ice velocity `u_ice`/`v_ice`** (and/or the `a_ice` IC blend weight). The EVP is iterative + sensitive
and the sea-ice port was validated at the *climate/advection* level ([[project_sea_ice_port_state]]),
NOT step-1 bit-identity — exactly the kind of small step-1 transient that perturbs `u_ice`. Bo's
polar part is the analogous over-ice heat-flux difference (ice thermodynamics, depends on ice state).

## Definitive next step (ready to run)

Cross-dump the step-1 ice state both codes, keyed by node gid, 16r (same partition):
`a_ice, m_ice, u_ice, v_ice, srfoce_u, srfoce_v, stress_iceoce_x/y, stress_node_surf_x/y`.
- **C:** add an env-gated dump in `fesom_step.c` (has `ctx->ice` + `forcing`) right at the KPP
  point, reusing `fesom_kpp_dump_nodes`.
- **Fortran:** add the matching dump in `ice_oce_coupling.F90` (the `stress2oce` routine has
  `ice%`/`u_ice`/`a_ice`/`u_w` + `stress_node_surf` in scope) at step 1.
Then diff with `scripts/kpp_dump_diff.py` + bin by lat. Expected: `srfoce≈0` and `a_ice` ~matching →
the signal is `u_ice`/`v_ice` (EVP). If so, the root is the step-1 EVP solve; chase the EVP inputs
(`stress_atmice`, ice strength P*, substep count, the EVP α/β) C-vs-Fortran — reuse
[[reference_evp_dump_diagnostic]]. A quick climate-level cross-check: add `a_ice`/`uice` to
`scripts/kpp_climate_compare.py` (K9 outputs them) to see whether the ice state stays close in the
mean (transient) or drifts (structural).

**Heat/water-flux (Bo) tiny non-polar residual (~1e-9, rel ~1e-3):** a secondary, much smaller
bulk-flux-assembly difference (SW/LW/sensible/latent or SST/humidity). Dump `heat_flux`/`water_flux`
components if the ice-stress fix doesn't also resolve it. Files: `fesom_bulk.c`, `fesom_forcing.c`,
`fesom_ice_thermodynamics` ([[feedback_runoff_fold_when_ice]]).

## Pointers

- Evidence run: `job_kpp_k5_validate` → `/work/ab0995/a270088/port/kpp_k5_validate/`.
- Cross-check script: `scripts/kpp_bldepth_xcheck.py` (input↔output attribution).
- Related memory: [[feedback_bvfreq_smoothing_gap]] (the *other* gap K5 surfaced, now FIXED),
  [[feedback_phc_rank_dependent]], [[feedback_runoff_fold_when_ice]],
  [[feedback_array_size_vs_reader_loop]] (a prior forcing-halo OOB lesson).
