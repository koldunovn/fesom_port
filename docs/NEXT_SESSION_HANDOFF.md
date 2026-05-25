# FESOM2 C-port — session handoff (2026-05-25)

Literal one-to-one C port of FESOM2 (CORE2 config, dt=1800, linfs). **KPP is now the default
vertical-mixing scheme** (production CORE2); `FESOM_MIX_SCHEME=PP` opts back to PP. Goal: C
climate ≈ Fortran climate.

---

## 0. TL;DR status

> **KPP SESSION 2026-05-25 — DONE. The KPP vertical-mixing port is COMPLETE, validated, and is now
> the DEFAULT scheme.** Phases K0–K10 on `fesom2_port/main` (HEAD ~`8d0cdbc`):
> - **Per-routine bit-faithful vs Fortran** (bldepth/blmix/enhance/combine/viscAE) via the
>   **controlled-replay** technique (inject Fortran inputs → isolate algebra from forcing noise),
>   worst ~3.18e-13.
> - **End-to-end climate** C+KPP vs Fortran+KPP (`fortran_2yr_dt1800` + `fortran_kpp_5yr_d1800`):
>   global SST/SSS RMS 0.005–0.017 / 0.002–0.007 over 5 yr, biases ~0, NON-growing; 5-yr stable
>   (no uv>5). Figures: `docs/kpp_5yr_figures/`.
> - **Two faithfulness fixes** found+committed: bvfreq N² horizontal smoothing + bottom-pad
>   (`feedback_bvfreq_smoothing_gap`, improved PP climate ~10×); `ice_gamma_fct` 0.25→0.5
>   (`feedback_namelist_over_codedefault`). Both re-baseline the climate (re-validation owed).
> - **One open issue, deferred to next session:** the step-1 forcing difference = the step-1 EVP
>   ice velocity `u_ice` (polar only; bounded). Spec + kickoff: `docs/FORCING_STEP1_DIFFERENCE.md`
>   + `docs/NEXT_SESSION_PROMPT.md`.
> - **Owed:** a fresh PP+KPP 2-yr re-validation run to bank the two re-baselines.

> **PRIOR — VALIDATION SESSION 2026-05-24(b) (PP, pre-KPP) — DONE. The port is validated at the production config.**
> Ran the definitive climate + stability runs at HEAD `6d65bd3` (first full runs with the
> cyclic-seam fix), all at dt=1800/864r/PP/linfs:
> - **2-yr climate** (`clim_2yr_dt1800` vs `fortran_pp_2yr`): global SST/SSS area-wtd RMS
>   **0.040 / 0.016** (on target), biases ~0 and non-drifting yr1→yr2. **W Med cyclic-seam
>   fix confirmed in the climate**: seam ΔSSS rms 0.055→0.006; the seam nodes moved off the
>   corrupt flat-Atlantic 37.351 to 37.384 == Fortran 37.388.
>   (`scripts/clim_validation_2yr.py` → `docs/validation_2yr_d1800/`.)
> - **5-yr stability** (`stab_5yr_dt1800`): **completed all 5 yr (1958-62) clean, exit 0,
>   max uv 4.65 m/s — never tripped the uv>5 guard.** Resolves the old dt=1200 year-4.9 trip.
>   CAVEAT: the 4.65 peak is the recurring autumn transient (dt=1200 hit 6.84 there); 0.35
>   m/s margin to the guard → watch / raise the guard for runs >5 yr.
> - **5-yr drift vs a NEW matched Fortran reference** (`fortran_pp_5yr_d1800`, built this
>   session = `work_linfs_pp_5yr_d1800`, dt=1800 PP linfs 5yr): EXCELLENT. Surface SST/SSS
>   C−F bounded ±0.05°C/±0.004PSU no trend; interior ΔT 0-100m −0.012°C flat, 100-700m
>   −0.004, ≥700m **LOCKED 0.000**; ice area/volume track; spatial maps near-zero yr1→yr5
>   (no regional runaway). (`scripts/drift_5yr_d1800.py` → `docs/drift_5yr_d1800/`.)

The port reproduces the Fortran CORE2 climate to **global SST/SSS RMS ≈ 0.04 / 0.016 PSU**
(annual 1958, vs `fortran_pp_2yr`). This session closed the major systematic surface
biases (river-mouth SSS, equatorial SST, Arctic SSS, Red Sea seasonal cycle) and fixed a
cyclic-longitude IC bug. The remaining surface residual is one **explained, non-bug**
feature: the Danish-straits/Kattegat front (~1 PSU, a resolution limit, dynamical — the
IC there is proven faithful).

`git HEAD = 6d65bd3`, working tree clean (only untracked docs/plots/jobs). The built
binary `build/fesom_port` matches HEAD (all fixes below).

---

## 1. What we did this session (commits, newest first)

| commit | what | validated effect (vs `fortran_pp_2yr`, annual 1958) |
|---|---|---|
| `6d65bd3` | **IC cyclic-lon fix**: PHC lon axis is periodic; the C mis-ported the `Nlon+2` wrap-halo (shifted real end cells ±360 + dropped the 0.5/359.5°E columns) → corrupt IC within 2° of lon 0/360 (W Med flat 37.35 vs PHC 37.7). Now proper halo padding. | IC nodes \|ΔS\|>0.1 vs Fortran: **262→21**; W Med seam gone. IC-only (no full-year climate run yet — see §2). |
| `575d812` | **IO**: added `FESOM_VAR_3D_ELEM_IFACE`; `Av` now written on `nz` (48 interfaces). | output-only (physics unchanged) |
| `2cb2bdb` | **IO**: `Kv`/`bvfreq` written on `nz` (48 interfaces), not `nz_1` (47). | output-only |
| `de4ced1` | **SSS+chl monthly read skipped February** (calendar-timing: Fortran updates end-of-month w/ `+1`→next; C updates start-of-month via `calendar_crossed`, so the ported `+1` double-counted → read months 1,3,4,…,12). Fix: drop the `+1`. | Red Sea RMS **0.357→0.029**, Arctic SSS 0.056→0.012, **global SSS 0.028→0.016** |
| `b8e1d36` | **`ref_sss_local`**: virtual-salt-flux reference salinity = LOCAL surface SSS (namelist.tra `ref_sss_local=.true.`); C had used global const 34.7 (the code default). | Arctic SSS RMS **−50%** (0.112→0.056) |
| `1703aab` | **MFCT** 3rd-order horizontal tracer advection (replaced interim 2nd-order central). | river-mouth SSS RMS **−70%**, equatorial SST **−61%**, global SST 0.110→0.041 |

(Pre-session, already committed: `86511f0` AB2-epsilon dt=1800 fix, `20349d9` sw_pene,
`7dc5150` albedo 0.1, `8dac997` wind rotation, `ad1d3b7` ice_strength 0.5.)

### Investigations that did NOT need a code fix
- **Baltic/Kattegat ~1 PSU residual** = **dynamical, NOT a port bug.** Proven: geometry
  bit-identical; currents/shear match; halocline Kv identical; PP formula correct (the
  surface Kv difference is a Ri feedback, not a bug); and crucially the **IC there matches
  the Fortran** (apples-to-apples at dt=1s, ΔS≈0.001). It develops over the year from the
  autumn PP-mixing feedback on a sub-grid front at the Danish straits. Memory:
  `feedback_phc_rank_dependent`, `feedback_periodic_read_timing`.
- Method that settled it: dump C step-0 IC vs run Fortran at **dt=1s** (1 step,
  ~1e-3 evolution) → clean IC compare on the same partition. Binning IC diffs by longitude
  is what exposed the cyclic-seam bug (86% within 2° of lon 0/360). Memory:
  `feedback_cyclic_seam_interp`.

---

## 2. What's left to do (open items, roughly in priority)

1. ~~Definitive climate run with the FINAL binary.~~ **DONE 2026-05-24(b)** — `clim_2yr_dt1800`
   vs `fortran_pp_2yr`: global SST/SSS RMS 0.040/0.016 ✓, W Med cyclic-seam fix confirmed
   in the climate ✓ (see §0 TL;DR + `docs/validation_2yr_d1800/`). Tagged `validated-d1800-5yr`.
2. ~~Long-term stability / 5-yr.~~ **DONE 2026-05-24(b)** — `stab_5yr_dt1800` ran a full 5 yr
   clean at dt=1800 (max uv 4.65, no `uv>5` trip; guard now at `fesom_main.c:1139`), and the
   drift vs the new matched Fortran ref `fortran_pp_5yr_d1800` is tiny+non-growing
   (`docs/drift_5yr_d1800/`). CAVEAT carried forward: the recurring autumn transient peaks
   4.65 m/s (0.35 margin) — for runs **>5 yr** watch it / raise the guard.
3. **Baltic** (optional): documented as a resolution limit, not a bug. Only revisit if a
   Danish-straits mesh/resolution change is on the table — it's not a port issue.
4. **Sea-ice Task F** (memory `project_sea_ice_port_state`): 256-rank month validation —
   the only remaining sea-ice task.
5. Housekeeping: scratch IC-dump runs can be deleted — `/work/ab0995/a270088/port/{ic_864,
   ic_8,fort_ic864}` (diagnostic only).

---

## 3. How to BUILD the C

```bash
source /sw/etc/profile.levante
source /home/a/a270088/port2/fesom2_port/env.sh
make -C /home/a/a270088/port2/fesom2_port/build fesom_port
# binary: /home/a/a270088/port2/fesom2_port/build/fesom_port
```
(NEVER use mambaforge for the build env; use serial netcdf — Fortran does gather+serial I/O.
On a login node, strip `OMPI_MCA_*` if set. See memory `feedback_levante_build`.)

## 4. How to RUN the C (SLURM, 864 ranks, 7 nodes)

Invocation: `fesom_port  MESH  OUT  dt  nsteps  snap_every  PHC  year`
- `dt=1800`, `nsteps=17600` ≈ 1 yr; `snap_every=-1` → monthly streams only; `>0` → also
  writes `snap_NNNNNN.nc` (step-0 snapshot = the interpolated IC, before any timestep).
- `MESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2`
- `PHC=/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc`

Copy an existing job (they're the template), edit `OUT`/`--job-name`, `sbatch`:
```bash
cd /home/a/a270088/port2/fesom2_port
cp job_monthfix_1yr job_<name>          # edit OUT=/work/ab0995/a270088/port/<name>, job-name, comment
mkdir -p /work/ab0995/a270088/port/<name>
sbatch job_<name>
```
Job header: `-p compute --nodes=7 --ntasks=864 -A ab0995`. Runtime ≈ 11 min/yr.
**Use a unique OUT_DIR per job** (concurrent jobs sharing OUT clobber logs — memory
`feedback_unique_outdir_per_job`).
Output per run: `{sst,sss,ssh,temp,salt,u,v,w,Kv,Av,bvfreq,density,a_ice,m_ice,m_snow}.fesom.<year>.monthly.nc`.

## 5. How to RUN the Fortran reference

- **Existing reference** (use this for comparison): `/scratch/a/a270088/fortran_pp_2yr`
  (864 ranks, dt=1800, 1958–1959, PP, Sweeney chl). Mesh diag: `fortran_pp_2yr/fesom.mesh.diag.nc`.
- **Executable**: `/home/a/a270088/fesom27/fesom2/bin/fesom.x`. Setup: `…/work_core/`
  (namelists + `job_levante`). The two Fortran source trees `port2/fesom2/src` (port
  reference) and `fesom27/fesom2/src` (built fesom.x) are byte-identical for the routines
  we checked.
- **To re-run** (e.g. a fresh IC dump or a matched experiment): copy work_core namelists to
  a fresh `ResultPath` dir, set `ResultPath`, `run_length`/`run_length_unit`, `step_per_day`
  (dt=86400/step_per_day; `=86400`→dt=1s for an IC dump), output via `&nml_list io_list`
  (e.g. `'salt',1,'s',8`; one entry per var, no duplicates), and a `fesom.clock`
  (`0.0 1 1958\n0.0 1 1958` for cold start). Submit at `--ntasks=864`. Working example:
  `/work/ab0995/a270088/port/fort_ic864/` (job_ic + namelists).
  - **CAVEAT — verify `which_ALE`**: `work_core/namelist.config` currently shows
    `which_ALE='zstar'`, but the C port is **linfs** and `fortran_pp_2yr` matches the C
    globally (so it was linfs). Set `which_ALE='linfs'` before any climate re-run, or you
    won't be comparing like-for-like. (IC interpolation is ALE-independent, so the dt=1s IC
    dumps were fine as-is.)

## 6. How to PLOT comparisons

Python env + path (always):
```bash
PYTHONPATH=/home/a/a270088/PYTHON \
  /work/ab0995/a270088/mambaforge/envs/nereus/bin/python  <script>
# pipe through:  2>&1 | grep -vE "getfattr"     (harmless ncdump attr warnings)
```
Uses `nereus` for unstructured plotting; mesh from `fortran_pp_2yr/fesom.mesh.diag.nc`
(has lon, lat, elem_area, face_nodes, nlevels(_nod2D), zbar(nz)/Z(nz1), nod_part).

Scripts in `scripts/` (edit the `C=`/`F=` run-dir paths at the top of each):
- `surf_diff_maps.py` — SST+SSS difference maps + regional area-weighted stats.
- `mfct_validation.py` — 3-way (Fortran / C-before / C-after) SST+SSS maps + stats.
- `reflocal_validation.py` — Arctic-focused SST/SSS + north-polar ΔSSS.
- `monthfix_validation.py` / `monthfix_figures.py` — marginal-sea seasonal cycles + maps.
- `baltic_amazon_probe.py` — regional maps (Baltic, Amazon) + node localization + evolution.
- `kattegat_kv_clean.py` — Kv/N² profile compare on the nz/48 grid (post IO fix).
- `ic_c_vs_fortran.py` — **C step-0 IC vs Fortran IC** at the same partition (set FIC to the
  Fortran IC dump; vars are `S`/`T` in the C snapshot, `salt` in the Fortran).
- `ic_rank_compare.py` — C IC at two rank counts (rank-dependence of extrap).

Output figures: `docs/validation_eps_2yr_dt1800/*.png`.
Annual-mean recipe: read `var.fesom.1958.nc` (Fortran) / `var.fesom.1958.monthly.nc` (C),
mask `|x|<1e29`, `np.nanmean(axis=0)`; area-weight with node areas from `elem_area`+`face_nodes`.
**Gotcha**: Kv/bvfreq/Av are now on `nz` (48 interfaces); temp/salt/u/v on `nz_1` (47 layers).

---

## 7. Key paths (quick reference)

- C source: `/home/a/a270088/port2/fesom2_port/src/` ; build: `…/build/fesom_port`
- Fortran source (port ref): `/home/a/a270088/port2/fesom2/src/`
- Fortran exe + setup: `/home/a/a270088/fesom27/fesom2/{bin/fesom.x,work_core}`
- Mesh: `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2` (dist_8,16,32,…,512,864)
- PHC IC: `/home/a/a270088/FESOM_port/fesom2/tests/data/INITIAL/phc3.0/phc3.0_winter.nc`
  (== `/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc`, md5-identical)
- Forcing/runoff/SSS-restoring: `/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/`
  (`CORE2_runoff.nc` Time=1 constant; `PHC2_salx.nc` 12-month seasonal — see `reference_paths`)
- Reference Fortran run: `/scratch/a/a270088/fortran_pp_2yr`
- C runs this session: `/work/ab0995/a270088/port/{mfct_1yr,reflocal_1yr,monthfix_1yr,
  kvfix_1yr,avfix_1yr,ic_864,ic_8,fort_ic864}`  (best climate = `monthfix_1yr`)
- Memory index: `/home/a/a270088/.claude/projects/-home-a-a270088-port2/memory/MEMORY.md`
- SLURM: account `ab0995`, partition `compute`, 864 ranks = 7 nodes.

---

## 8. PROMPT FOR THE NEXT SESSION

> Continue the FESOM2→C port (CORE2, dt=1800, linfs, PP). Read
> `docs/NEXT_SESSION_HANDOFF.md` first, plus the auto-memory in
> `~/.claude/projects/-home-a-a270088-port2/memory/`. HEAD is `6d65bd3`; the build is
> current. This session's goal: **run the definitive climate validation with all fixes and
> confirm long-term stability.**
>
> 1. Build (`source /sw/etc/profile.levante && source env.sh && make -C build fesom_port`),
>    then run a **2-year C run at 864 ranks, dt=1800, year 1958** (clone `job_monthfix_1yr`
>    → new OUT dir). This is the first full run that includes the cyclic-seam IC fix.
> 2. Compare to `/scratch/a/a270088/fortran_pp_2yr` with the `scripts/` plotters
>    (`surf_diff_maps.py`, `monthfix_validation.py`): confirm global SST/SSS RMS stays
>    ~0.04/0.016, and check the **western Mediterranean** SSS now tracks the Fortran
>    (validates `6d65bd3`).
> 3. Check **multi-year stability**: does the C survive past year ~4.9 at dt=1800 without
>    tripping the `uv>5` guard (`fesom_main.c:1100`)? See memory `project_5yr_milestone`.
> 4. Leave the Baltic/Kattegat ~1 PSU residual alone unless asked — it's a documented
>    resolution limit (dynamical, IC proven faithful), not a port bug.
>
> Work one-to-one with the Fortran (no simplifications/inventions). Dual-instrument C vs
> Fortran with output fields when diagnosing — and isolate IC interpolation from evolution
> by running the Fortran at dt=1s for apples-to-apples IC compares.
