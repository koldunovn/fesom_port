# Next-session prompt — mEVP port COMPLETE; open user decisions

(2026-06-11: the mEVP port was executed M0→M6 in one session and PASSED every
gate. Branch **`mevp`** in THIS worktree (`fesom2_port_zstar`), tag
**`mevp-validated-2026-06-11`**, history a74aacf→close-out. NOT merged to main.)

## State

- `FESOM_WHICH_EVP=1` selects mEVP (Bouillon 2013, `src/fesom_ice_maevp.c` =
  literal port of `EVPdynamics_m`); unset/0 = std EVP (default, byte-identical
  to main@2714071 — gated M0/M1/M2); 2/garbage abort at init.
- Validation: M3 dump-diff at the FP floor (s1 u_aux ≤5e-14→4.7e-11 rel);
  M4 cross-rank 1/8/16r integrals ≤5e-4 (edge-concentrated chaos, no seams);
  M5 2yr climate yr1 RMS 0.0049/0.0024 + ice cycles ≤0.3% + vector gate
  ±0.014° + diff-of-diffs corr 0.96/0.90/0.88; M6 5yr clean, peak |uv| 2.56
  (margin 2.44 vs the 5.0 trip — std-EVP baseline was 4.65/0.35).
- Docs: `docs/validation_mevp/{M2_fidelity_review,M3_dumpdiff,M6_acceptance}.md`,
  `docs/validation_mevp_2yr/RESULTS.md`, `docs/mevp_reference_namelists/`.
- Memory: `project_mevp_port_state`.
- Fortran side: `work_mevp` (2yr reference, frozen binary) + `work_mevp_dump`;
  the EVPdynamics_m dump patch stays UNCOMMITTED in `port2/fesom2` (KPP/TKE
  convention).

## Decisions taken 2026-06-11 (same day, user-approved)

1. **MERGED**: `main` fast-forwarded to the mEVP close-out, then to the
   IO-frame commit (`git fetch . mevp:main`; main == mevp).
2. **C IO vector-frame: DONE** (75406d3): stream vector pairs ((u,v),
   (uice,vice) + the io_meandata pair list) are written in GEOGRAPHIC
   components by default; `FESOM_IO_VECTOR_FRAME=rotated` opts out
   bit-exactly (gate job 25524763). Snapshots/debug dumps stay
   native-frame — the snap-based byte-gate lineage vs baseline_2714071
   REMAINS VALID. 2yr legs + 5yr rerun at 75406d3: all M5/M6 numbers
   reproduced digit-for-digit; Gate 3 re-verified by DIRECT comparison.
   ⚠️ Implementation lesson: never add code to physics TUs —
   `feedback_tu_codegen_bitgate` (fesom_mesh.c edit shifted
   -O3/-ffp-contract codegen → last-ULP pgf at step 0 → bitwise gates
   broke via chaos; the transform is a private copy in fesom_io_stream.c).

## Still open (USER)

- aEVP (whichEVP=2) later? mEVP as default ever? (std EVP stays default.)

## Quick reference

- Outputs root `/work/ab0995/a270088/port/mevp/`: baseline_2714071,
  fdump_16r + cdump_16r (M3), m1_byteident, m2_smoke, m4_xrank,
  c_mevp_2yr + c_evp_2yr (the matched C 2×2 legs), c_mevp_5yr,
  fortran_mevp_2yr.
- Switches: `FESOM_MIX_SCHEME=KPP(default)|PP|TKE`, `FESOM_ALE=linfs(default)|zstar`,
  `FESOM_WHICH_EVP=0(default)|1`, `FESOM_IO_VECTOR_FRAME=geo(default)|rotated`, `FESOM_EVP_DUMP_DIR` (mEVP dumps: Q/U0/F, P,
  it1/2/60/120, UF — diff with `scripts/maevp_dump_diff.py`).
- Build C: `source /sw/etc/profile.levante && source env.sh && make -C build
  fesom_port` (fresh: `bash -l configure.sh`).
  Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`.
