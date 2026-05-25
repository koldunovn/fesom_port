# Step-1 forcing difference C-vs-Fortran (ustar / Bo) — RESOLVED (Fortran cold-start bug)

**Status: RESOLVED 2026-05-25.** Root cause = a **Fortran (FESOM2) cold-start bug**, not a C
porting bug and not cold-start sensitivity. The C port is the *correct* one. See the "RESOLUTION"
section at the end for the mechanism, proof, and fix. The investigation notes below are kept as the
historical record of how it was narrowed (the earlier "EVP u_ice" framing was the symptom; the cause
is one rotation upstream of the EVP).

## RESOLUTION (2026-05-25): Fortran cold-start wind-rotation bug

**The bug** (`gen_surface_forcing.F90`): on a rotated grid the geo→rotated (g2r) wind/stress rotation
is applied **only in `sbc_do`**, gated on `do_rotation_wind`, which is set **only inside the
`getcoeffld` (coefficient-refresh) branch** (~line 1542; rotation block ~1547). The cold-start init
`nc_sbc_ini` (~603-700) loads the first coefficients (`getcoeffld`) and builds the wind
(`data_timeinterp`) but **never rotates them**. At step 1 `force_newcoeff=.false.` (same year) and
`rdate` is inside the first forcing window, so `sbc_do` does not refresh → the rotation block is
skipped → the first forcing window's wind stress on **both ocean and ice** stays in the **geographic
frame**, mis-directed by the local g2r angle.

**Proof** (per-substep EVP dump, both codes, 16r, step 1; `scripts/evp_dump_diff.py` +
`evp_dump_geo.py`, new `F` dump point = `stress_atmice`): every EVP input bit-identical (a_ice/m_ice
IC, ice_strength, inv_mass, sigma, u_rhs) **except `stress_atmice`**, which differed at all 126858
nodes. The diff was a **pure rotation** — per-node `|stress|` ratio C/F = 1.0000, direction off by
**exactly the local g2r vector angle** (`|delta − theta| = 0.000` at every node). C rotates
`jra->u_wind` from step 1 (`fesom_jra55.c:692`); Fortran's first-window wind is unrotated, so
`F_stress = r2g(C_stress)`.

**Why it hid the whole port:** rotation is magnitude-preserving, so `ustar`, fluxes, SST/SSS/SSH,
ice area/volume and ice/ocean speed *ratios* are all blind to it; only the EVP velocity *vector*
exposes it. And it is transient (self-corrects after the first forcing-window crossing, ~6 steps for
3-hourly JRA at dt=1800), so multi-year climate matched regardless (re-baseline 2-yr RMS ~0.004).

**Fix:** rotate the initial coefficients in `nc_sbc_ini` (mirror the `sbc_do` rotation, guarded by
`l_xwind`/`l_xstre` `.and. rotated_grid`). Applied to `port2/fesom2`; the C port needs **no change**.
PR to upstream `FESOM/fesom2` (bug confirmed in the latest main).

---
*Historical investigation notes (how it was narrowed) follow.*

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

## PINPOINTED (2026-05-25) — it's the step-1 sea-ice STATE (dominantly EVP u_ice)

Cross-dumped the step-1 ice state from BOTH codes (tag `iceforce`, 16r same partition, C in
`fesom_ice_coupling.c` `fesom_ice_oce_fluxes_mom` + Fortran `ice_oce_coupling.F90` `oce_fluxes_mom`;
`/tmp/icf_cmp.py`). Polar (|lat|>50) C-vs-Fortran:

| field | max\|Δ\| | typical val | comment |
|---|---|---|---|
| `srfoce_u/v` | **0.0 exactly** | 0 | ocean at rest → ice-ocean stress = ρ·cd·\|u_ice\|·u_ice·a_ice |
| `u_ice`/`v_ice` | 0.24 / 0.27 | ~0.3 | **rel ~80% — the dominant difference (EVP velocity)** |
| `a_ice` | 1.49e-2 | ~0.9 | ice concentration differs ~1.5% (36803 nodes) |
| `m_ice` | 3.23e-2 | ~2.0 | ice thickness differs ~3% |

Ice-free latitudes: every field bit-identical. So the **root is the step-1 sea-ice state** —
dominantly the **EVP ice velocity `u_ice`** (rel ~80%, from a `u_ice`=0 cold-start IC the EVP is
maximally sensitive), with a smaller `a_ice`/`m_ice` (concentration/thickness) contribution. This
propagates through `stress_iceoce → stress_node_surf → ustar` (the polar 3e-3 ustar diff that KPP
sees). The sea-ice was validated at the climate/advection level ([[project_sea_ice_port_state]]),
NOT step-1 bit-identity — exactly this.

**⚠️ dump caveat:** the `iceforce` dump's `stress_node_surf` column is captured at the ice-coupling
point and is NOT the stress KPP reads — it disagrees with the `prestep` ustar by ~0.12 even in the
tropics in BOTH codes (the atm-ocean stress is re-finalized between the ice coupling and the ocean
step). Use `prestep` ustar (bit-identical tropics 8e-16, polar 3e-3) as the KPP-input ground truth;
the `iceforce` a_ice/m_ice/u_ice/v_ice/srfoce columns are the reliable ice-state.

## Remaining next level (the actual fix, post-KPP)

Two sub-questions, both sea-ice-subsystem step-1 transients:
1. **Why `u_ice` differs (dominant):** chase the EVP inputs/iteration C-vs-Fortran at step 1 —
   `stress_atmice` (wind-on-ice, [[feedback_wind_rotation_g2r]]), ice strength P* (a_ice/m_ice,
   [[feedback_ice_strength_halffactor]]), the EVP substep count + α/β. Reuse
   [[reference_evp_dump_diagnostic]] (per-substep EVP dumps) to localize where the C and Fortran
   ice velocity diverge within the step.
2. **Why `a_ice`/`m_ice` differ (~1.5-3%):** ice IC interpolation vs one thermo/advection step —
   compare the step-0 IC ice (C snap_000000) to the Fortran IC, then one step.

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
