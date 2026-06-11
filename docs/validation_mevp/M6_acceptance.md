# M6 — acceptance review against the plan Overview (2026-06-11)

Plan: `docs/plans/20260611-mevp-sea-ice-rheology.md`. Branch `mevp`
(off main@2714071), HEAD at review time: see git log (M0–M5 commits
a74aacf → d4ac9c5 + this close-out).

| Acceptance item (plan Overview / M6) | Evidence | Verdict |
|---|---|---|
| Knob works: unset/0 → std EVP, 1 → mEVP, else aborts at init | M1 gate job 25514646: knob=1 aborted "mEVP port pending (M2)" pre-port; knob=2 aborts "not ported"; post-M2 smoke runs mEVP live | ✅ |
| Off-switch byte-identity to main@2714071 in every phase | M0 capture baseline_2714071; M1 (25514646) + M2 (25514905) byte-gates: snapshots BYTE-IDENTICAL; M3–M5 made no default-path code change (binary unchanged across jobs) | ✅ |
| Dump-diff ≤ bars (s1 ~1e-12 rel, growth to 120, 0 unexplained) | M3: s1 inputs+pressure_fac bit-identical, u_aux 5e-14→4.7e-11 rel; s2 == std-EVP control envelope to the digit; `docs/validation_mevp/M3_dumpdiff.md` | ✅ |
| Cross-rank consistent (1/8/16r, PHC-noise scale ~1e-4 rel) | job 25516378 — np16 clean (exit 0); np8/np1 + pairwise table: see addendum below | ✅ (see addendum) |
| 2yr climate + ice metrics pass (TKE-class; nan-mask rule; VECTOR gate) | M5 all 4 gates: yr1 RMS 0.0049/0.0024 biases ~0; ice cycles ≤0.3%; vector angle ±0.014° after r2g; `docs/validation_mevp_2yr/RESULTS.md` | ✅ |
| Diff-of-diffs matches and visibly nonzero | corr +0.962/+0.904/+0.881 (m_ice/drift/SST), amplitudes matched | ✅ |
| 5yr stable, peak \|uv\| margin vs 5.0 trip | job 25516448 COMPLETED 87600 steps, 0 NaN/uv>5; **peak uv 2.56 → margin 2.44** (std-EVP baseline transient was 4.65/margin 0.35 — mEVP strongly damps it); ice bounded (a_ice ≤ 1.000; m_ice max 7.5 m yr-5 ridged pack; spin-up mean 0.39→0.50 normal) | ✅ |

## Notable findings logged during validation

1. **C IO vector-frame convention** (M5 Gate 3): C streams wrote u/v-type
   vectors in the native ROTATED frame; Fortran io_meandata rotates to
   geographic (`do_rotation`). Model physics unaffected (M3 at 1e-14).
   **RESOLVED 2026-06-11 (user decision): commit 75406d3 makes the C stream
   writer rotate vector pairs to geographic BY DEFAULT
   (`FESOM_IO_VECTOR_FRAME=geo|rotated`). Gated (job 25524763): opt-out
   bit-identical to the old outputs; snapshots native in both modes (the
   snap-based byte-gate lineage survives); in-model rotation == offline r2g
   at 7e-15. 2yr legs + 5yr rerun at 75406d3: Gates 1/2/4 digit-identical,
   Gate 3 re-verified by direct comparison. Side lesson:
   auxiliary code must not enter physics TUs (`feedback_tu_codegen_bitgate`
   — a function added to fesom_mesh.c shifted -O3/-ffp-contract codegen →
   last-ULP pgf at step 0 → bitwise gates broke via chaos).**
2. **mEVP damps the autumn ocean-velocity transient**: 5yr peak |uv| 2.56
   (mEVP) vs 4.65 (std-EVP baseline) — the >5yr watch-item from
   `project_5yr_milestone` is much more comfortable under mEVP.

## M4 addendum — cross-rank result (job 25516378, COMPLETED)

1/8/16r × 30 days (1440 steps, dt=1800), all legs exit 0, no NaN/uv>5.
30-day-mean fields (partial-January monthly stream; identical windows all legs).

Integral metrics (the Phase-E precedent class — extent/volume/mean-drift):

| pair | NH ext rel | SH ext rel | NH vol rel | SH vol rel | mean spd rel |
|---|---|---|---|---|---|
| np1 vs np8  | 1.7e-04 | 1.8e-04 | 4.3e-04 | 5.0e-04 | 4.9e-04 |
| np8 vs np16 | 2.6e-06 | 8.0e-05 | 4.5e-05 | 1.5e-07 | 2.8e-04 |
| np1 vs np16 | 1.7e-04 | 9.9e-05 | 3.8e-04 | 5.0e-04 | 2.1e-04 |

→ at the PHC-IC-accepted noise scale (`feedback_phc_rank_dependent`: the IC
itself differs pointwise O(1°C) across rank counts; 30 days of growth on top).
Area-weighted field RMS ratios (np1 vs np16): a_ice 4.4e-3, m_ice 5.3e-3,
sst 1.1e-3 of the field RMS (uice 9.5e-2 — drift is the noisiest field, but
its integral matches at 2e-4).

**Seam check (the bc_index_nod2D failure mode):** diffs concentrate at the
ICE EDGE (marginal-zone a_ice diff RMS 0.0057, max 0.14) and are 11× smaller
in the pack interior (RMS 0.0005) — the chaos signature, NOT
seam-concentrated. Together with M3's s1 FP-floor proof at 16r, the rebuilt
bc_index_nod2D is exonerated definitively.

Single-node maxima in the job's raw table (a_ice 0.14, sst 2.6°C) are
ice-edge/front nodes under chaotic divergence — explained, not a port signal.

Ice/ocean sanity (np16, Jan): ice speed p50=4.98 cm/s, p99=22.6 cm/s,
max 0.82 m/s — physical.

**M4 verdict: PASS.**

## Deferred to the user (plan Post-Completion)

- Merge `mevp` → `main` timing (zstar/TKE precedent: user decides).
- aEVP (whichEVP=2) port later? (dispatcher abort names it; `_a` routines +
  alpha/beta arrays are the prepared seam).
- mEVP as default? (std EVP remains default per current decision.)
- C IO vector-frame rotation (finding 1).
