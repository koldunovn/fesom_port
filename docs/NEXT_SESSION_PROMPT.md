# Next-session prompt — KPP port, continue at Phase K5 (`bldepth`)

Copy the block below as the next session's opening prompt.

---

> Continue the FESOM2→C port (CORE2, dt=1800, linfs): the **KPP vertical-mixing port**,
> strictly faithfully. **GUIDING PRINCIPLE: THIS IS A PORT. The Fortran
> (`/home/a/a270088/port2/fesom2/src/oce_ale_mixing_kpp.F90`) has worked for decades — do
> NOT simplify, optimize, or invent. Reproduce every formula/constant/loop-bound/branch/
> halo-exchange/ordering exactly; port the CORE2 reference namelist values, not code defaults.**
>
> **State: Phases K0–K4 are DONE, validated, and committed** on `fesom2_port/main`
> (HEAD `c63474a`; commits K0 `13bfe95`, K1 `88812a0`, K2 `8114a8b`, K3 `9c80362`,
> K4 `c63474a`). Scaffolding + `oce_mixing_kpp_init` + `wscale` + `ri_iwmix` + the `ddmix`
> gate are in; PP stays byte-identical with KPP off. Build is current.
>
> 1. Read **`docs/plans/20260524-kpp-vertical-mixing.md`** (authoritative spec; K0–K4 boxes
>    are `[x]`, K5 is next). Then the auto-memory — especially **`project_kpp_port_state`**
>    (has the K5 setup + the validation methodology learned in K3), `feedback_port_used_terms`,
>    `feedback_namelist_over_codedefault`, `feedback_phc_rank_dependent` (the IC-noise floor),
>    `feedback_write_loops_halo`, `feedback_tracer_stride_nl`, `feedback_node_vs_cell_geometry`,
>    `reference_evp_dump_diagnostic`.
> 2. Build. C: `source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`.
>    Fortran (port2/fesom2, **Intel** build): `source /sw/etc/profile.levante && source
>    env/levante.dkrz.de/shell && make -C build fesom.x -j8 && make -C build install`
>    — the **install step is mandatory**: `bin/fesom.x` loads `lib64/libfesom.so` (RPATH
>    `$ORIGIN/../lib64`), NOT the fresh `build/lib/`, so without `install` you run stale KPP.
> 3. **Implement Phase K5** (`oce_ale_mixing_kpp.F90`; highest risk — OBL depth, historically
>    buggy per `FRESH_START.md:489`). Three pieces, extend the incremental driver
>    `fesom_kpp_mixing` (currently zero-viscA → ri_iwmix → dump → guard):
>    - **(a) store `dbsfc`** in `src/fesom_eos.c` (`fesom_pressure_bv`): the local `dbsfc1`
>      is already computed at `:136` for db_max — write it into `aux->dbsfc[FESOM_NODE3D(n,nz,nl)]`
>      and add the bottom fill `dbsfc(nzmax)=dbsfc(nzmax-1)` (Fortran oce_ale_pressure_bv.F90:337).
>      **Store it UNCONDITIONALLY** — PP never reads `aux->dbsfc`, so the PP path stays
>      byte-identical without a gate (re-confirm with the byte-identity job).
>    - **(b) driver pre-step** (`:307-352`, before ri_iwmix in the Fortran order): `dVsq`
>      (eqn 21, shear of UVnode re the surface, `dVsq(nzmax)=dVsq(nzmax-1)`), `ustar`
>      (`sqrt(sqrt(stress_node_surf_x²+y²)·density_0_r)`, eqn 2), `Bo`
>      (`-g·(sw_alpha(nzmin)·heat_flux/vcpw + sw_beta(nzmin)·water_flux·S(nzmin))`, eqns A2c/d).
>      Inputs all present: `forcing->{heat_flux,water_flux,stress_node_surf}`, `aux->{sw_alpha,
>      sw_beta}`, `tracers` S, `density_0_r = 1/FESOM_DENSITY_0`, `FESOM_VCPW`.
>    - **(c) `bldepth`** (`:628-808`): bulk-Ri search (`Ritop=zk·dbsfc`, `Vtsq=zk·ws·√|bvsq|·Vtc`,
>      `Rib=Ritop/(dVsq+Vtsq+epsln)`, `Ricr=0.3` linear interp of `hbl`), the `use_sw_pene`
>      `bfsfc` branch (sw_3d interpolated to hbl), Ekman/Monin-Obukhov limits gated
>      `bfsfc>0 .and. nzmin==1` (`cekman=0.7`/`cmonob=1.0`, `KPP_CEKMAN`/`KPP_CMONOB` already
>      defined), `smooth_hbl` (.false., skip), final `kbl` search, `caseA`. Uses `wscale`
>      (K2, static `kpp_wscale`), `aux->{dbsfc,bvfreq,sw_alpha}`, `forcing->sw_3d`, `kpp->
>      {ustar,Bo,dVsq}`, `mesh->{zbar (|zbar_3d_n|, interface depths), coriolis_node}`.
>      Note bldepth loops `myDim` and uses `zbar` (interfaces), not `Z` (mid-layer).
> 4. **Validate K5** by dual-instrumentation at step 1 vs the EXISTING Fortran reference dumps
>    in `/work/ab0995/a270088/port/kpp_fdump/dump/` (already current: `kpp_dump_s1_prestep`
>    [ustar,Bo], `_dVsq`, `_dbsfc`, `_bldepth` [hbl,kbl,bfsfc,stable,caseA]). Add matching C
>    dumps in `fesom_kpp_mixing` (mirror `kpp_dump_riiwmix`); run the C KPP dump-run
>    (`FESOM_MIX_SCHEME=KPP FESOM_KPP_DUMP_DIR=... nsteps=1`, model on `job_kpp_k3_validate`),
>    diff with `scripts/kpp_dump_diff.py`. **Caveat (from K3):** at step 1 from rest `uv=0 →
>    shear=0 → dVsq=0` (trivial); `ustar/Bo/dbsfc/hbl/kbl/bfsfc/caseA` are meaningful. Inputs
>    differ at the PHC-rank IC-noise floor → threshold fields (`kbl`, `hbl` near weak
>    stratification) show O(1) marginal-node diffs that are INPUT noise, not bugs. **Bin
>    `hbl`/`kbl` diffs by region**; to isolate the algebra, dump the driving inputs (dbsfc,
>    bvfreq) too and confirm disagreements track input differences (mirror
>    `scripts/kpp_ri_bvfreq_check.py`, **interior-only — exclude edge-copy levels**). The
>    velocity-driven part of bldepth (dVsq term) is untested at rest → covered by K9.
>    **Re-run the byte-identity job** (`job_kpp_k3_validate` PP part, or a 200-step PP run vs
>    `/work/ab0995/a270088/port/kpp_k0_byteident/head`) — K5 touches `fesom_eos.c` (shared).
> 5. Commit K5. Then **K6** (`blmix_kpp` :1109-1281 — BL coeffs `blmc[3]`, `dkm1[3]`, shape
>    G(σ), `ghats`), **K7** (`enhance` :1300-1344 + the `smooth_blmc` 3×exchange+smooth_nod +
>    the combine MAX-in-BL + the 4 exchanges + node→elem `viscAE`=Σ/3 + bottom fill + `minmix`
>    floor), **K8** (wire: `aux->Kv(:,node)=diffK(:,node,1)` over myDim+eDim, oce_ale.F90:3518-22;
>    `fesom_tracer_diff.c` UNCHANGED — single Kv for T and S; nonlocal flux stays gated off),
>    **K9** (C+KPP dt=1800 2yr 864r vs `/scratch/a/a270088/fortran_2yr_dt1800`). One routine
>    per phase, dual-instrument before the next, commit per phase.
>
> Reminders: gate-only (don't implement) `ddmix` (done, K4) and the `ghats`→tracer nonlocal
> flux (`use_kpp_nonlclflx=.false.`); vertical diffusion uses a SINGLE `Kv=diffK(:,:,1)` for T
> and S — do NOT route `diffK(:,:,2)` into salinity. The incremental driver guards itself so
> only a `FESOM_KPP_DUMP_DIR` run (1 step) proceeds until K8 completes it.

---

## Quick reference (paths + commands)

- **KPP plan (spec):** `docs/plans/20260524-kpp-vertical-mixing.md` (K5 section has every detail).
- **C source:** `src/fesom_kpp.{c,h}` (driver `fesom_kpp_mixing`, static `kpp_wscale`/`kpp_ri_iwmix`,
  dump helpers `fesom_kpp_dump_*`, `kpp_dump_riiwmix`, constants `KPP_*`). `src/fesom_eos.c:136`
  (`dbsfc1` to store). `src/fesom_step.c` (dispatch, `mo_convect` shared after either scheme).
- **Fortran ref:** `/home/a/a270088/port2/fesom2/src/oce_ale_mixing_kpp.F90` — instrumented with
  `kpp_dump_*` gates (UNCOMMITTED in that tree, like the `[FSSH]` debug; Intel build).
  Add a `bldepth`/prestep dump point in the driver if you want fresh Fortran dumps; the current
  ones in `/work/.../kpp_fdump/dump/` already cover prestep+bldepth+dbsfc at step 1.
- **Fortran KPP dump run:** `cd /home/a/a270088/port2/fesom2/work_kpp_dump && sbatch job_kpp_dump`
  (KPP/linfs/dt=1800/16r, dumps at step 1 to `/work/ab0995/a270088/port/kpp_fdump/dump/`).
- **C validation jobs (models):** `job_kpp_k3_validate` (KPP dump-run nsteps=1 + PP byte-identity),
  `job_kpp_k2_validate`, `job_kpp_k1_validate`. HEAD byte-identity baseline:
  `/work/ab0995/a270088/port/kpp_k0_byteident/head`.
- **Diff/reference scripts:** `scripts/kpp_dump_diff.py` (gid-keyed, auto-discovers tags),
  `scripts/kpp_ri_bvfreq_check.py` (threshold/input-noise isolation, **interior-only**),
  `scripts/kpp_byteident_check.py`, `scripts/kpp_{init,wscale}_reference.py`.
- **Mesh/PHC/forcing:** see `docs/NEXT_SESSION_HANDOFF.md` §7.
- **Python:** `/work/ab0995/a270088/mambaforge/envs/nereus/bin/python` (PYTHONPATH=/home/a/a270088/PYTHON).
