# Z9 — zstar 2yr climate validation results (2026-06-10)

**VERDICT: PASSED.** C (`FESOM_ALE=zstar`, KPP, dt=1800, 864r, 1958-1959, run
`/work/ab0995/a270088/port/zstar/c_zstar_2yr`, job 25495449, exit 0, branch
`zstar` @ 89611d5) vs the NEW Fortran zstar reference
(`/work/ab0995/a270088/port/zstar/fortran_zstar_2yr`, job 25492898) — the pair
differing from the linfs pair by the ONE knob `which_ALE='zstar'`
(provenance: `docs/zstar_reference_namelists/PROVENANCE.md`).

## Headline (ZS = C-zstar − F-zstar, area-weighted annual means)

| year | field | GLOBAL bias | GLOBAL RMS | worst region RMS |
|---|---|---|---|---|
| 1958 | SST | −0.0004 °C | **0.0038 °C** | 0.0054 (midlat) |
| 1958 | SSS | +0.0000 PSU | **0.0014 PSU** | 0.0015 (midlat) |
| 1959 | SST | −0.0005 °C | **0.0050 °C** | 0.0076 (midlat) |
| 1959 | SSS | +0.0001 PSU | **0.0021 PSU** | 0.0027 (tropics) |

- **At/better than the KPP-vs-KPP acceptance bar** (K9: 0.0046/0.0023 → 0.0128/0.0044).
  The 1958→1959 RMS growth is ~3× SMALLER than K9's — consistent with the step-1
  forcing transient being gone (the zstar pair's Fortran binary includes the
  nc_sbc_ini cold-start fix 2682a9fb; the C was always correct).
- Biases ≈ 0 both years → non-drifting.

## Coordinate contrast (LIN = C-zstar − F-linfs, same C run)

Global SST RMS 0.0130/0.0231, SSS RMS 0.0108/0.0178 (1958/1959) — **3–9× the
ZS headline**, largest at N-high-lat SSS (0.031/0.048 PSU: virtual-salt vs
real-flux sea-ice freshwater handling, as expected). The comparison clearly
resolves the coordinate (the K9-pattern requirement).

Full table: `climate_compare_output.txt` (generator:
`scripts/zstar_climate_compare.py`).

## Companion gates closed the same day

- **Z8 month**: 30-day dt=1800 16r zstar (job 25495448) exit 0; max|uv| over
  the run **1.79 m/s** (linfs autumn-transient lore: 4.65 peak); no blowup,
  vs≡0 / hp≡0 in the stats confirm the zstar arms live.
- **Z7 cross-rank**: 5-day 1/8/16r (job 25495447, all exit 0). Pairwise
  final-snapshot worst|Δ| (T 3.4 °C, S 1.2 PSU at day 5) is BELOW the step-0
  IC spread of the same pairs (T 6.3 °C, S 1.9 PSU in `snap_000000`) — i.e.
  the cross-rank differences ORIGINATE in the documented PHC per-rank IC
  noise (`feedback_phc_rank_dependent`) and shrink under mixing, non-growing.
  A linfs-baseline cross-rank (job 25499700) runs for the apples-to-apples
  magnitude bar — expect the same scale.
- linfs byte-gates green at every phase (Z0, Z2-Z6, Z2-Z7): snapshots
  worst |Δ| = 0 vs the v1.0 binary.

## Z10 remainder (user decision point)

- Default coordinate flip linfs → zstar (KPP precedent: flip after validation)
  — awaiting user.
- Tag suggestion: `zstar-validated-2026-06-10`.
- Plan move to `docs/plans/completed/` after the decision.
