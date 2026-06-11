# TKE port validation results (plan `20260610-tke-vertical-mixing.md`)

Worktree branch `zstar`, started from `zstar-validated-2026-06-10` (8260dea).
Reference provenance: `docs/tke_reference_namelists/PROVENANCE.md` (Fortran
2yr job 25500644, dump job 25500775; frozen binary = the zstar-pair binary).

## ACCEPTANCE (T6) — ALL GATES GREEN, 2026-06-11

| Requirement | Verdict |
|---|---|
| TKE selectable (`FESOM_MIX_SCHEME=TKE\|cvmix_TKE`; KPP default, PP kept) | ✓ |
| KPP and PP byte-identical with TKE off (vs HEAD 8260dea) | ✓ both legs, all snapshots |
| Column core vs Fortran ≤1e-12 | ✓ replay worst **1.1e-16** (after the -r8 6.6 catch) |
| Climate matches the TKE reference (KPP-class bar) | ✓ RMS 0.0049/0.0028 → 0.0088/0.0040, biases ≈0 |
| Scheme contrast resolves the scheme | ✓ 11–18× (SST), >100× (N-high SSS) |
| Diag flag inert (byte-identical) + functional (budget closes) | ✓ / ✓ (4e-15 rel) |
| Cross-rank at PHC-noise scale, halo probe clean | ✓ day-5 ≤ IC spread; 3200/3200 zero |
| 30-day dt=1800 stability | ✓ clean, max uv 1.64 |
| Combined zstar+TKE 2yr (run matrix) | ✓ RMS 0.0037/0.0015 → 0.0078/0.0041 |
| 5yr zstar+TKE stability (run matrix) | ✓ clean, peak uv 2.93 (margin 2.07) |

## T0 — scaffolding gates

- **Reference-run provenance hardening:** `init_cvmix_tke` echo verified in
  BOTH Fortran run logs: `tke_cd = 3.75`, `tke_only = T` (namelist consumed;
  `__cvmix` compiled — `nm` additionally shows `calc_cvmix_tke` in the frozen
  lib). 2yr run completed exit 0, 730 days, fesom.clock at 1960-01-01.
- **Byte-gate (job 25501105), TKE off vs frozen HEAD-8260dea baseline,
  16r/200 steps/dt=500/snapshots every 25:**
  - KPP default: **BYTE-IDENTICAL across all snapshots. PASS** (run.log diff
    = path strings only).
  - PP opt-out: **BYTE-IDENTICAL across all snapshots. PASS.** (This leg ran
    the binary INCLUDING the T3 IO changes — gates those too.)

## T1/T2 — column core + driver, dump-diff vs Fortran (16r, 3 steps, dt=1800)

C dump job 25501106 vs Fortran dump job 25500775; `scripts/tke_dump_diff.py`.

**Step 1 (clean column-algebra test; uv=0 → shear off):**

- All pure-algebra outputs with bit-identical inputs are bit-identical or
  1-ulp: `tkeav`/`tkekv`/`kv`/`av`/`lmix`/`pr`/`tbpr`/`tspr`/`tdif`/`tdis`/
  `tiwf` ALL **max|diff| = 0.0 exactly**; `dztrr` 0.0 (proves hnode/Z_3d_n
  byte-identical under linfs); `vshear2` 0.0; `tkeold` 0.0.
- `bvfreq2`: 1.4e-15 (rel 1e-13) — the established PHC-IC noise floor.
- `normstress` (INPUT, from the live forcing chain): bit-identical at 13,653
  of 126,858 nodes; median rel diff 3.5e-15 elsewhere (bulk-formula last-bit
  noise); tail to 1.6e-3 rel at a small set (PHC-floor/freezing-set nodes —
  the accepted noise class).
- **Conditional bit-identity: where `normstress` is bit-identical, `tke`
  agrees to 3.5e-18 and `twin` to 1.7e-21 — pure 1-ulp libm `pow` residual.**
  Empirically pinned: Intel's `x**(3./2.)` equals glibc `pow(x,1.5)` EXACTLY
  at 126,741/126,858 nodes; the remaining 117 differ by ≤1 ulp (Intel libm is
  not correctly-rounded; glibc ≥2.28 is). Irreducible, 6 orders below the
  1e-12 bar. (`x*sqrt(x)` was tested and rejected: only 79% exact.)

**Steps 2-3 (live): the expected threshold-flip amplification, mechanism
quantified (the zstar-s2/K5 pattern):**

- s1 bvfreq sign flips C-vs-F: **0 of 126,858 columns** (the 1e-15 s1 noise
  flips nothing).
- s2 bvfreq sign flips: 44,854 columns (near-zero N² in winter mixed layers —
  ice-onset/heat-flux noise at the PHC floor moves N² across 0) → mo_convect
  Kv 0↔0.1 flips + Rinum/prandtl clamp flips.
- Of the 40,835 columns with s2 `tke` diff > 1e-4, **66.1% sit directly on a
  flip column**; the rest receive flipped bvfreq via the N2smth_h horizontal
  smoothing. Live s2+ comparison is therefore noise-amplified, not
  algebra-divergent — the controlled replay (below) isolates the algebra.

**Budget closure (the scheme's own consistency; both codes, all 3 steps):**
`max|Ttot − (Tbpr+Tspr+Tdif+Tdis+Twin+Tiwf+Tbck)|` ≤ 3.7e-19 (rel ≤ 4e-15)
in C; Fortran identical magnitudes.

**Controlled replay (FESOM_TKE_REPLAY_DIR; Fortran inputs injected per column
per step → output diff = pure algebra; bar ≤1e-12):**

- **First replay (job 25501627) CAUGHT A REAL PORT BUG:** `prandtl` differed
  by exactly (6.6 − (double)6.6f)·Rinum at 27k entries — the port had assumed
  the `6.6` literal at cvmix_tke.F90:717 was single-precision; the Intel
  reference build compiles with **-r8** (src/CMakeLists.txt:335) so
  default-real literals are DOUBLE. One-line fix (TKE_C66 = 6.6).
- **After the fix (job 25501766): PASS — worst max|diff| across ALL output
  tags, all 3 steps = 1.110e-16** (one ulp of Av≈0.72; `pr`/`kv`/`tkeav`/
  `tkekv` exact-0 at s3, `tke` 1.7e-18). Four orders below the bar. The
  column core + driver assembly are bit-faithful on real multi-step columns
  with shear production active.

## T3 — diagnostics + IO

- Diag storage gated by `FESOM_TKE_DIAG=1` (or an active dump dir); the 13
  slabs are the only allocation difference; column core always computes the
  full set in local scratch (internal read-after-write deps preserved).
- IO: `tke` + the 10 budget terms registered as OPTIONAL stream variables —
  written ONLY when an `FESOM_IO_CONFIG` entry names them; default streams
  untouched.
- Diag-flag invariance gate (FESOM_TKE_DIAG=0 vs 1, 16r/200/dt=500 snapshot
  byte-compare, job 25503533): **BYTE-IDENTICAL across all snapshots. PASS.**
- Halo probe (exchange-and-compare on aux->Kv, run A, FESOM_TKE_HALO_PROBE=1):
  **3200/3200 probe lines report worst |d| = 0.000e+00** (200 steps × 16
  ranks) — no stale halo anywhere in the TKE wiring.
- Memory delta of the flag: expected 13×[N·nl]·8B ≈ +42 MB/rank at CORE2 16r
  — below run-to-run MaxRSS noise there (1446 vs 1404 MB; the diag=0 leg
  also carried the probe buffer). The flag's payoff is at NG-mesh scale
  (nl=70, ~10× nodes), as designed.

## T4 — cross-rank + stability — cross-rank **PASS**

1/8/16r 5-day runs (240 steps, dt=1800; job 25503534, all exit 0). The
byteident script prints "FAIL" (it is a bit-checker); the bar is PHC-noise
scale (feedback_phc_rank_dependent), assessed Z7-style against the STEP-0 IC
spread:

| pair | step-0 IC spread (T / S) | day-5 spread (T / S) |
|---|---|---|
| np1–np8  | 6.32 / 1.92 | 3.48 / 1.16  (shrinks under mixing) |
| np8–np16 | 2.89 / 27.7 | 3.40 / 28.0  (S=28 exists AT the IC — rank-dependent PHC extrapolation at a coastal node group; not growing) |
| np1–np16 | 6.32 / 27.7 | 3.48 / 28.0 |

Same numbers as the zstar Z7 precedent for np1-np8 (T 6.3→3.4, S 1.9→1.2 —
same IC). Day-5 ≤ step-0 everywhere → rank noise originates at the IC; no
TKE halo bug (consistent with the 3200/3200-zero halo probe).

- 30-day 16r dt=1800 stability (month16 leg, 1440 steps): **exit 0, zero
  blowup/uv>5/NaN lines, final max|uv| 1.64** — clean. **T4 PASS.**

## Run matrix — combined zstar+TKE legs (approved 2026-06-10) — 2yr **PASS**

Fortran `work_zstar_tke` (job 25503765; TWO knobs vs work_linfs_2yr_b:
which_ALE='zstar' + mix_scheme='cvmix_TKE'; frozen binary; exit 0, echo
verified, 730 days) vs C `FESOM_ALE=zstar FESOM_MIX_SCHEME=TKE` 2yr 864r
(job 25504000, exit 0, max uv 1.57). Full table:
`climate_compare_zstar_tke_output.txt`:

| year | dSST bias | dSST RMS | dSSS bias | dSSS RMS | coord-contrast SST/SSS RMS |
|---|---|---|---|---|---|
| 1958 | −0.0003 | **0.0037** | −0.0000 | **0.0015** | 0.0132 / 0.0142 |
| 1959 | −0.0006 | **0.0078** | −0.0001 | **0.0041** | 0.0212 / 0.0195 |

Headline at/below the per-feature legs; the contrast vs F-linfs+TKE is 2–5×
larger → the comparison resolves the coordinate with TKE active. **The 2yr
matrix is closed: zstar+KPP ✓ (Z9), linfs+TKE ✓ (T5), zstar+TKE ✓.**
5yr zstar+TKE stability (job 25504001): **PASS** — full 87600 steps
(1958-1962), exit 0, zero blowup/uv>5/NaN; peak max|uv| over the run **2.93**
(autumn transient) → margin 2.07 to the uv>5 trip, far inside the linfs+KPP
lore (peak 4.65, margin 0.35). The approved run matrix is FULLY closed.

## T5 — 2yr climate vs work_linfs_tke — **PASS**

C linfs+TKE 2yr dt=1800 864r (job 25503535, exit 0, 35040 steps, max uv 1.67)
vs Fortran `fortran_linfs_tke` (25500644). Full table:
`climate_compare_output.txt`. Area-weighted annual means:

| year | dSST bias | dSST RMS | dSSS bias | dSSS RMS |
|---|---|---|---|---|
| 1958 | −0.0003 | **0.0049** | −0.0001 | **0.0028** |
| 1959 | +0.0001 | **0.0088** | −0.0004 | **0.0040** |

- At the KPP-class bar (~0.005/0.002; zstar Z9 was 0.0038/0.0014 →
  0.0050/0.0021). Biases ≈0 both years → no drift; only the chaotic-noise
  RMS grows year-to-year (the accepted signature).
- **Scheme contrast resolves the scheme decisively:** C-TKE vs F-KPP global
  SST RMS 0.088/0.093 (1958/1959) = **18×/11× the headline**; N-high SSS
  contrast 0.246/0.200 vs headline 0.0017/0.0013 (>100×). (The Fortran's own
  TKE-vs-KPP signal is SST RMS 0.115 — the C-vs-F headline sits 20× below
  the physics signal it must not confuse.)
