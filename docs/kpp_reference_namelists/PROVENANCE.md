# CORE2 KPP reference namelist provenance (Phase K0)

Archived 2026-05-24 for the KPP vertical-mixing port
(`docs/plans/20260524-kpp-vertical-mixing.md`).

## Source

Copied verbatim from `/home/a/a270088/port2/fesom2/work_core/`.

The end-to-end KPP reference run is `/scratch/a/a270088/fortran_2yr_dt1800`
(KPP, dt=1800, years 1958-1959). Its work-dir lineage is `work_core` /
`work_linfs_d1800`. The run dir itself does not retain its namelists
(`ResultPath = /scratch/a/a270088/fesom27/`, later copied/renamed), so the
namelists are archived here from the work dirs.

**`work_core/namelist.oce` and `work_linfs_d1800/namelist.oce` are byte-identical;
likewise `namelist.tra`** (`diff -q` reported IDENTICAL for both), so the KPP
configuration is unambiguous regardless of which work dir produced the run.

## Confirmed values (plan "CORE2 reference KPP configuration" table)

| Symbol | Value | Source | Status |
|---|---|---|---|
| `mix_scheme` | `'KPP'` (→ `mix_scheme_nmb=1`) | namelist.oce:67 | ✓ |
| `Ricr` | `0.3` | namelist.oce:70 | ✓ |
| `A_ver` | `1.e-4` | namelist.oce:19 | ✓ |
| `visc_sh_limit` | `5.0e-3` | namelist.oce:66 | ✓ |
| `diff_sh_limit` | `5.0e-3` | namelist.tra:97 | ✓ |
| `K_ver` | `1.0e-5` | namelist.tra:101 | ✓ |
| `Kv0_const` | `.true.` | namelist.tra:100 | ✓ (→ Kv0_background_qiang NOT used) |
| `double_diffusion` | `.false.` | namelist.tra:105 | ✓ (→ ddmix gate only) |
| `use_kpp_nonlclflx` | **absent** → `.false.` | oce_modules.F90:168 | ✓ (→ nonlocal flux gate only) |
| `use_instabmix` | `.true.` | namelist.tra:88 | ✓ (mo_convect kept after KPP) |
| `use_sw_pene` | `.true.` | namelist.config:108 | ✓ (Bo/bfsfc/bldepth sw branch) |
| `step_per_day` | `48` → dt = 86400/48 = **1800 s** | namelist.config:29 | ✓ |
| `ref_sss_local` | `.true.` | namelist.tra:116 | ✓ (already in C, b8e1d36) |

Module-parameter / in-routine constants (not in namelist; from
`oce_ale_mixing_kpp.F90`): `Riinfty=0.8` (:737), `smooth_blmc=.true.` (:73),
`smooth_hbl=.false.` (:74), `minmix=3.0e-3` (:408), `cekman=0.7`/`cmonob=1.0`
(:478-479). These are verified against the Fortran source in K1.
