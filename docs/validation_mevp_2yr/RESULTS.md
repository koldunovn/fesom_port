# M5 — 2yr climate validation: C mEVP vs Fortran mEVP (2026-06-11)

Runs (all dt=1800, 864r, cold start 1958, KPP+linfs; C binary = mevp branch
@7ec565d-era build, byte-gated default; F binary = the frozen pair shared by
BOTH Fortran legs):

| leg | run | job |
|---|---|---|
| Cm (C mEVP) | /work/ab0995/a270088/port/mevp/c_mevp_2yr | 25516446 |
| Fm (F mEVP, reference) | /work/ab0995/a270088/port/mevp/fortran_mevp_2yr | 25514027 |
| Ce (C std-EVP, same binary) | /work/ab0995/a270088/port/mevp/c_evp_2yr | 25516447 |
| Fe (F std-EVP, same frozen lib) | /work/ab0995/a270088/port/zstar/fortran_linfs_2yr_b | (zstar matrix) |

Tool: `scripts/mevp_climate_compare.py` (incl. the Gate-3 r2g frame fix).

## Gate 1 — ocean climate (bar: TKE-class, biases ~0, non-drifting): PASS

dSST/dSSS (Cm−Fm), area-weighted:

| year | SST bias | SST RMS | SSS bias | SSS RMS |
|---|---|---|---|---|
| 1958 | −0.0002 | **0.0049** | −0.0000 | **0.0024** |
| 1959 | −0.0007 | 0.0125 | −0.0003 | 0.0042 |

Year-1 numbers are exactly TKE-class (TKE: 0.0049/0.0028; zstar: 0.004/0.001);
year-2 growth to 0.0125/0.0042 is chaotic-divergence accumulation, within the
KPP-5yr accepted envelope (0.005–0.017) and dominated by the TROPICS (0.0176),
where the rheology is irrelevant and the cross-rheology column shows the same
number (0.0182) — i.e. the comparison saturates at the C-vs-F two-run noise
floor there. Where the rheology acts, the mEVP match clearly beats the
cross-rheology distance: N>60 SST RMS 0.0052 (Cm−Fm) vs 0.0198 (Cm−Fe), ~4×;
S<−60 0.0092 vs 0.0190. Biases ≈0 both years ⇒ non-drifting.

## Gate 2 — ice seasonal cycles (NH+SH extent / volume / drift): PASS

1959 monthly, `nan_to_num→0` before all time means
(feedback_nan_mask_temporal_average): extent matches to ≤0.04×10⁶ km²
(≲0.3%) in every month, both hemispheres; volume to ≤0.02×10³ km³;
drift-speed means to ≤0.07 cm/s (worst: Feb SH 10.11 vs 10.04). Full tables
in /tmp/m5_results.txt (regenerate with the script).

## Gate 3 — ice-velocity VECTOR check: PASS (after frame alignment) + FINDING

**Finding (diagnostic-layer, not a model bug):** the C IO stream writes
uice/vice (and u/v) in the NATIVE ROTATED frame (`resolve_uice` accumulates
`ice->uice` raw); Fortran io_meandata sets `do_rotation=.TRUE.` for the
(uice,vice) pair and writes GEOGRAPHIC components. Raw comparison shows the
exact signature `feedback_magnitude_invariant_masks_rotation` predicts:
|u| ratio ≡ 1.000 (isometry) with position-dependent angles to ±40°.
Model internals are proven in the rotated frame at 1e-14 (M3), so this is
purely an output-convention difference. M5 compares with the C vectors
rotated r2g (vector_r2g port, now in the script).

With frames aligned (Mar+Sep 1959, NH+SH, a>0.15, |u|>0.005 m/s):

| | ratio median | ratio IQR | angle median | angle IQR | angle wmean |
|---|---|---|---|---|---|
| Mar NH | 1.0006 | ±0.1% | −0.014° | ±0.05° | −0.04° |
| Mar SH | 1.0008 | ±0.3% | −0.002° | ±0.16° | −0.42° |
| Sep NH | 1.0004 | ±0.1% | +0.004° | ±0.07° | +0.05° |
| Sep SH | 1.0001 | ±0.1% | −0.000° | ±0.07° | +0.03° |

No systematic rotation: the port's drag/Coriolis/stress algebra is
angle-faithful. (Residual vecRMS/spdRMS ≈ 0.10–0.15 is the same year-2
two-run chaotic divergence as Gate 1, angle-unstructured.)

**Open decision (user):** make the C IO rotate vector pairs r2g at output
(matching Fortran convention + standard analysis expectations)? Touching it
changes every C vector output ⇒ breaks byte-gate lineage; needs its own
re-baseline. Until then: C vector outputs are rotated-frame; analyses must
rotate (the compare script does).

## Gate 4 — diff-of-diffs (option live & faithful): PASS

(Cm−Ce) vs (Fm−Fe), 1959 annual means, area-weighted pattern correlation
over the ice mask (SST: global wet):

| field | pattern corr | wstd C | wstd F | max\|C\| | max\|F\| |
|---|---|---|---|---|---|
| m_ice | **+0.962** | 0.0152 | 0.0165 | 0.43 | 0.51 |
| drift speed | **+0.904** | 0.0729 | 0.0721 | 0.67 | 0.69 |
| SST | **+0.881** | 0.0072 | 0.0074 | 0.36 | 0.36 |

The C reproduces the Fortran's mEVP−EVP rheology CONTRAST (pattern +
amplitude), and the contrast is visibly nonzero ⇒ the knob is live and
faithful, not inert.

## Verdict

**M5 PASS on all four gates.**
