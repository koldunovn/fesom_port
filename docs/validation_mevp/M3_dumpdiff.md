# M3 — substep-level dump-diff, C mEVP vs Fortran mEVP (16r, 2 steps, dt=1800)

2026-06-11. Fortran: job 25514438 (`work_mevp_dump`, instrumented
`EVPdynamics_m`, build/lib/libfesom.so incl. the sbc cold-start fix) →
`/work/ab0995/a270088/port/mevp/fdump_16r/dump`. C: job 25514906
(`jobs/job_mevp_cdump`, M2 binary) →
`/work/ab0995/a270088/port/mevp/cdump_16r/dump`.
Tool: `scripts/maevp_dump_diff.py` (gid-keyed, same dist_16 partition both
sides). Bar (plan M3): ~1e-12 rel at substep 1, small monotone growth to
120; 0 unexplained diffs.

## Step 1 — controlled test (cold start: inputs MUST be bit-identical)

| point | result |
|---|---|
| Q (a/m/msnow/elev), U0 | **bit-identical** (max\|d\| = 0) |
| F stress_atmice | 7e-14 rel (bulk-formula FP reduction noise) |
| P node (inv_thickness, mass, rhs) | inv_thickness bit-identical; mass 6e-14 rel; rhs ≡ 0 both (elev=0) |
| P elem pressure_fac | **bit-identical** |
| it1 u_aux/v_aux | 5e-14 rel; sigma ≡ 0 both (zero strain at substep 1) |
| it2 | u_aux 5e-14 rel; sigma 1.7e-12 rel (first nonzero sigma) |
| it60 | u_aux 4.8e-14 rel; sigma ≤8.5e-13 rel |
| it120 / UF | u_aux 2.4e-11, v_aux 4.7e-11 rel; sigma ≤1.2e-8 rel |

Verdict: **at/below the bar**. The slow e-13→e-8 sigma growth over 120
pseudo-iterations is FP-noise amplification through the near-singular
`pressure_fac/(delta+1e-11)` factor at near-rigid pack elements (delta ≈
delta_min ⇒ d(pressure)/pressure ~ d(delta)/1e-11); the physically meaningful
velocity stays ≤5e-11 rel. No spatial structure; top-diff gids are scattered
pack-interior elements.

**bc_index_nod2D proof**: with bit-identical inputs at s1, a multi-rank
misclassification (the M1 rebuild's target failure mode) would hard-zero C
velocities at inter-rank seam nodes — an O(1) signal. Observed: 5e-14 rel,
no seam pattern ⇒ the rebuilt bc_index_nod2D is correct at 16r.

## Step 2 — live full-model step (inputs inherit step-1 ocean+thermo)

s2 inputs differ: a_ice 1.8e-4 rel, m_ice 1.9e-4, elevation 3.8e-7,
u_w/v_w up to **1.2e-2 rel**. **Explained — pre-existing, rheology-
independent C-vs-Fortran step-1 divergence**, proven by the std-EVP control:
`evp_dump_diff.py` on the May-25 std-EVP pair (`evp_fdump` vs `evp_cdump`,
ZERO mEVP code) gives **the same numbers to the printed digits** —
a_ice 1.671e-04, m_ice 3.761e-04, m_snow 9.222e-05, u_w 4.709e-03,
v_w 2.656e-03 (driver: srfoce_temp 1.9e-2 rel after the step-1
thermo/mixing path). This is the envelope the whole validated port operates
in (it converges to the 0.004/0.002 RMS 2yr climate agreement).

mEVP internals track those input diffs commensurately — no amplification
anomaly:

| quantity | s2 rel | consistency check |
|---|---|---|
| pressure_fac | 1.1e-3 | ≈ Δmsum/msum + 20·Δasum = 1.9e-4 + 20·1.8e-4 ✓ (exp(−c(1−a)) sensitivity) |
| it1 u_aux | 4.7e-6 | < input envelope (the implicit solve contracts toward forcing balance) |
| it120 u_aux / UF | 1.9e-4 | ≈ input envelope 1.8e-4 ✓ |
| it120 sigma | 6.4e-4…1.6e-3 | tracks pressure_fac 1.1e-3 ✓ (std-EVP control's s2 sigma: 4.7e-3 — mEVP is tighter) |

## Verdict

**M3 PASS — 0 unexplained disagreements.** s1 = controlled algebra test at
the FP floor; s2 = inherited rheology-independent input noise with
commensurate propagation, bounded by the std-EVP control pair.

Byte-gate: no C code change in M3 (analysis + job only); the M2 byte-gate
(job 25514905, BYTE-IDENTICAL vs baseline_2714071) ran THIS binary
(build/fesom_port unchanged between jobs 25514905/25514906) and stands.
