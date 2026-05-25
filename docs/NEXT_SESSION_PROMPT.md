# Next-session prompt — fix the step-1 forcing difference (sea-ice EVP `u_ice`)

Copy the block below as the next session's opening prompt.

---

> Continue the FESOM2→C port (CORE2, dt=1800, linfs). **The KPP vertical-mixing port is COMPLETE
> and validated** (K0–K10, `fesom2_port/main`; KPP is now the DEFAULT `mix_scheme`, `FESOM_MIX_SCHEME=PP`
> opts back to PP). Per-routine bit-faithful vs Fortran (controlled-replay, ~3e-13); end-to-end
> climate C+KPP vs Fortran+KPP = SST/SSS RMS 0.005–0.017 / 0.002–0.007 over 5 yr, biases ~0,
> stable. **This session's task: fix the one remaining open issue — the step-1 forcing difference.**
>
> **GUIDING PRINCIPLE: strict faithful port — reproduce Fortran exactly, no simplification; port the
> reference run's NAMELIST values, not module defaults (`feedback_namelist_over_codedefault` — bitten
> 3×: ref_sss_local, ice_gamma_fct, …).**
>
> 1. **Read `docs/FORCING_STEP1_DIFFERENCE.md` (the authoritative spec + diagnosis)** and the
>    auto-memory: `project_forcing_step1_diff` (the target), `feedback_controlled_replay_validation`
>    (the validation technique), `reference_evp_dump_diagnostic` (the EVP per-substep dump tool),
>    `feedback_wind_rotation_g2r`, `feedback_ice_strength_halffactor`, `project_sea_ice_port_state`,
>    `feedback_namelist_over_codedefault`.
> 2. **The diagnosis so far (don't re-derive):** at step 1, `ustar`/`Bo` (KPP inputs) differ C-vs-Fortran
>    ONLY in the polar/sea-ice regions (ice-free latitudes bit-identical). Mechanism is the ice-ocean
>    stress `ρ·cd·|u_ice−u_w|·(u_ice−u_w)·a_ice + atmoce·(1−a_ice)` (C `fesom_ice_coupling.c:223` ==
>    Fortran `ice_oce_coupling.F90:120`, line-for-line). `srfoce`(u_w)=0 (ocean at rest), atm-ocean is
>    bit-identical → the root is the **step-1 EVP ice velocity `u_ice`** (polar rel ~80%; `a_ice`/`m_ice`
>    differ only ~1.5–3%, transported by the differing `u_ice`). **All EVP parameters already match**
>    C-vs-Fortran (whichEVP=0, Pstar=30000, c_pressure=20, delta_min=1e-11, evp_rheol_steps=120,
>    Cd_oce_ice=0.0055, theta_io=0) → it is NOT a parameter mismatch. So `u_ice` differs either from
>    **(a) genuine cold-start sensitivity** (faithful EVP amplifying tiny IC-interp noise from u_ice=0;
>    the climate validating at ~0.01°C + 5-yr stability argues for this) **or (b) a subtle EVP code
>    difference** (rheology/iteration ordering the climate mean hides).
> 3. **Decisive next step — per-substep EVP dumps.** Use `reference_evp_dump_diagnostic`
>    (`FESOM_EVP_DUMP_DIR`, env-gated per-substep `u_ice`/`v_ice`/stress dumps, BOTH codes) to compare
>    where the C and Fortran ice velocity diverge across the 120 EVP subcycles at step 1, 16r (same
>    partition → dump line i = local node i). Diverge from substep 1 ⇒ an INPUT differs — re-dump
>    `stress_atmice` (wind-on-ice) + ice strength `P*` (built from a_ice/m_ice) and diff. Diverge
>    gradually ⇒ the iteration/rheology (compare the C `fesom_ice_evp.c` rheology + stress-update
>    ordering line-by-line vs `ice_EVP.F90`). The `iceforce` dump (C `fesom_ice_coupling.c` +
>    Fortran `ice_oce_coupling.F90`, tag `iceforce`; `scripts/kpp_iceforce_cmp.py`,
>    `scripts/kpp_forcing_geo.py`) already gives a_ice/m_ice/u_ice/v_ice/srfoce/stress at step 1.
> 4. If it's (a) sensitivity: document it as accepted (bounded, climate-validated) and close the issue.
>    If it's (b) a code diff: fix it faithfully, re-run the iceforce dump to confirm `u_ice` matches,
>    then a 2-yr KPP climate re-check vs `fortran_2yr_dt1800`.
>
> **Also owed (re-baseline bookkeeping from the KPP session), do early:** (i) a fresh PP byte-identity
> baseline — the bvfreq N² fix (`feedback_bvfreq_smoothing_gap`) superseded `validated-d1800`, so the
> K0/K8 byte-identity gate must re-point to a post-fix PP run; (ii) a KPP 2-yr climate re-validation
> with `ice_gamma_fct=0.5` (the K10/figure runs used the pre-fix 0.25; effect is small but unconfirmed).
> Both fold into one fresh PP+KPP 2-yr run (864r, dt=1800) vs `fortran_pp_2yr` / `fortran_2yr_dt1800`.
>
> Build: C `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`. Fortran
> (port2/fesom2, Intel): `source /sw/etc/profile.levante && source env/levante.dkrz.de/shell &&
> make -C build fesom.x -j8 && make -C build install` (install mandatory — `bin/fesom.x` loads
> `lib64/libfesom.so`). Python: `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python`
> (PYTHONPATH=/home/a/a270088/PYTHON).

---

## Quick reference (paths + commands)

- **Forcing-diff spec:** `docs/FORCING_STEP1_DIFFERENCE.md` (mechanism, geographic data, the ready
  per-substep EVP plan). **General handoff:** `docs/NEXT_SESSION_HANDOFF.md`.
- **KPP plan (DONE):** `docs/plans/completed/20260524-kpp-vertical-mixing.md`.
- **Forcing diagnostics:** `scripts/kpp_iceforce_cmp.py` (ice-state C-vs-F by field, polar),
  `scripts/kpp_forcing_geo.py` (ustar/Bo diff by latitude band), `job_kpp_iceforce` (C iceforce dump run).
- **EVP code:** C `src/fesom_ice_evp.c` (whichEVP=0); Fortran `src/ice_EVP.F90`. Ice-ocean stress:
  C `src/fesom_ice_coupling.c` / Fortran `src/ice_oce_coupling.F90` (both carry the `iceforce` dump).
- **References:** Fortran KPP `/scratch/a/a270088/fortran_2yr_dt1800` (2yr) + `fortran_kpp_5yr_d1800`
  (5yr, `work_kpp_5yr_d1800`); C KPP 5yr `/work/ab0995/a270088/port/kpp_5yr_dt1800`.
- **5-yr figures:** `docs/kpp_5yr_figures/` + `scripts/kpp_5yr_figures.py`.
