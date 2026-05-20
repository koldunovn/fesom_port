# FESOM2 Fortran→C Port — Experience Report

**Author:** assistant, distilled from the full memory store + docs/ + plans/  
**Source span:** 2026-04-23 (FRESH_START) → 2026-04-25 (Phase G/Sea-ice complete)  
**Port size at end:** ~15,300 LoC C (vs. ~73,000 LoC Fortran source)  
**Validated configurations:**
- 1, 8, 16, 32, 256-rank CORE2 with PHC IC + JRA55-do forcing.
- Full ocean, sea ice (EVP + FCT + thermo + coupling), GM/Redi mesoscale eddy parameterisation.
- 256-rank model month completes; T.min held at −2.06 °C by sea-ice insulation.

This report is the digest the *next* port needs in front of it on day 1. It is
written to be **the second document a porter reads after `FRESH_START.md`** —
FRESH_START tells you *what* to port; this tells you *how to avoid the
specific traps we already paid for*. Where the two overlap, FRESH_START wins
on physics; this wins on process.

---

## Table of Contents

1. Executive summary
2. The single most important rule we learned
3. The four bug classes that ate the most debug time
4. Confusing parts of FESOM2 itself (not C-port-specific)
5. Confusing parts of the *port* (C/MPI/SLURM-specific)
6. Process discipline that paid off
7. The phase order that worked
8. Validation infrastructure that earned its keep
9. Things tried and failed (don't redo)
10. Suggestions for FESOM2 Fortran-side changes that would make the next port faster
11. A starter checklist for day 1 of the next attempt

---

## 1. Executive summary

The port succeeded on the second serious attempt, after several earlier
attempts that died from accumulated approximations. The decisive change:
**we stopped trying to be clever and ported line-by-line, with a per-substep
dump-shim driving validation from line 1 of C**. Every attempt to
"simplify" or "improve" Fortran physics cost weeks of debugging downstream.
Every `myDim_*` loop bound that "looked fine" in C against a `myDim+eDim`
Fortran original turned out to be a multi-rank halo bug.

The port is structured in five phases: ocean dynamics → MPI → sea ice →
GM/Redi (and KPP/zlevel/restart still pending). Each phase ended in a
documented validated state with bit-identity off-switches where possible.

---

## 2. The single most important rule we learned

> **Trust the Fortran code exactly. Port it line-by-line. Every shortcut
> you take will become a multi-week bug.**

We have a memory file (`feedback_port_fidelity.md`) dedicated to this
because it was repeatedly the difference between a working session and a
lost session. Every "small simplification" in earlier attempts caused a
silent bug that took weeks to localise:

- linearising α/β → broke OBL depth (stuck at 5 m).
- using `dt/areasvol` instead of `dt*area/areasvol` → imprinted mesh
  geometry on tracer evolution.
- swapping a lookup-table `ghats` for a closed-form `(1−σ)²` → wrong
  magnitudes everywhere.
- a single sign flip on `edge_vflux` (vertical advective flux) → no symptom
  for ~80 steps, then blowup at the deepest CORE2 element.

The mechanical changes that ARE allowed: 1-based → 0-based indexing,
column-major → row-major, `USE`-globals → struct passing. **Nothing else.**
Even when `area/areasvol == 1` for `linfs`, port the ratio so the code
survives the day someone enables `zlevel`.

---

## 3. The four bug classes that ate the most debug time

These are the patterns to grep for *before* you run, not after. Each has
its own memory file; this section names them and tells you the symptom.

### 3.1 Halo write-loop coverage (highest blast radius)

**Pattern:** Fortran writes an array with `do n=1, myDim+eDim`, the C port
writes with `for n=0; n<myDim`. Halo entries stay calloc-zero (or stale).
Any reader that touches halo gets bad data; the bug is silent on 1 rank,
appears at some rank counts and not others.

**Real instances (all in tree as fixes):**

| Commit | Site | Symptom |
|---|---|---|
| `7db862e` | `init_tracers_AB` valuesAB halo | 8/16/32/144-rank tracer drift→blowup at step ~95–150 |
| `20715c6` | FCT step a1 `fct_ttf_max/min` halo | 16-rank S.min drift |
| `028277d` | PHC + SSS/runoff halo | small residual `pgf` signature on 16 ranks |
| `55f6b6c` | JRA forcing arrays sized myDim only | Sea-ice EVP `m_ice→7e6 m` at halo, `ice_strength→1e9` |
| *(2026-05-20)* | `fesom_bulk_compute` Ch/Ce over myDim | **dt capped at 500**; SSH non-conservation, slow blow-up at dt=1800 |

**How to catch proactively:**

```bash
grep -nE 'for.*<.*myDim_(nod|elem)2D;' src/*.c | grep -v eDim
```

For every match, open the Fortran subroutine and check the loop bound.
If Fortran says `myDim+eDim`, the C must too. We learned to do this audit
**at write time**, not after a multi-rank divergence appears.

**The most insidious instance (the dt=1800 blocker, 2026-05-20).** This one
hid for weeks behind a *partial* fix and demands its own write-up because it
breaks two comforting assumptions:

- *"The field is allocated `myDim+eDim`, so it's fine."* It was. The bug was
  the **write loop** (`fesom_bulk_compute` computed `Ch_atm_oce`/`Ce_atm_oce`
  on `myDim` only). Fortran `ncar_ocean_fluxes_mode` (gen_bulk_formulae.F90:181)
  loops `myDim+eDim`. Allocation size is necessary, not sufficient.
- *"bulk_compute already halo-exchanges its outputs."* It exchanged
  `stress_node_surf`/`heat_flux`/`water_flux` — but **not** the exchange
  coefficients `Ch/Ce`, which sea-ice `thermodynamics` reads directly at eDim.
  An exchange audit that checks "is there a `call exchange_*` for the outputs"
  passes while the real consumer's input stays stale. (See the now-corrected
  "bulk_compute halo exchange … No observable effect" line in
  `MPI_PORT_REPORT.md` — that *was* this bug, fixed halfway.)

**Consequence chain (why the symptom was a time-step limit, not a crash):**
stale eDim `Ch/Ce` → sea-ice thermo diverges `a_ice` across ranks (444 nodes)
→ divergent wind `stress_surf` on **replicated boundary elements** (2D elements
that legitimately appear in `myDim` on >1 rank; their per-rank redundant
computation must be bit-identical or the SSH edge-flux assembly stops
telescoping) → `Σ ssh_rhs ≠ 0` (grew to 5.6e-6 relative) → non-conservative
free surface → velocity blow-up whose onset scales with `dt²` (the stiffness
matrix `M = area·(1+g·dt²·α²·Δ)`). At dt=500 the leak was tolerable; at
dt=1800 it blew up. Fortran was stable at dt=1800, which (per the literal-port
rule) *proved* it was a C bug, not a physics/linfs choice.

**Diagnostic recipe that cracked it (re-add in minutes if multi-rank drifts again):**
1. `Σ ssh_rhs` over **owned** nodes after assembly must be machine-zero
   (telescoping). Print `sum`, `max|ssh_rhs|`, ratio. 1-rank ~1e-16 vs
   multi-rank growing = a cross-rank assembly leak. (env `FESOM_DIAG_SSHRHS`)
2. Per-cut-edge dump (edges with a halo endpoint): same global node-pair on
   two ranks should give identical `(c1+c2)`. Identical geometry/`uv_rhs` but
   divergent `(c1+c2)` ⇒ the differing input is an element field (here:
   `stress_surf`).
3. **Exchange-and-compare probe** — the general "which node field is stale at
   eDim" finder: copy field → `exchange_nod(copy)` → diff `copy` vs original
   over `[myDim, myDim+eDim)`. Nonzero ⇒ inconsistent at eDim. Run it over the
   whole candidate-input list at once; it named `Ch/Ce` as the only stale ones.
   This probe is the fastest tool for this entire bug class.

### 3.2 Array allocation size vs. reader's loop bound (cross-module)

**Pattern:** Module A (e.g. JRA forcing) writes a per-node array sized
`myDim` only. Module B (e.g. sea-ice thermo) reads it in a `myDim+eDim`
compute-over-halo loop. The halo read is OOB into heap garbage; values
look insane on multi-rank only. Single-rank survives because halo is empty.

**Concrete failure:** sea-ice port crashed at step 2 only on multi-rank.
Always-on stderr probe (`if (m_ice[n]>1e3) fprintf …`) revealed
`m_ice = 7×10⁶ m` at a halo node where `Tair = 0 K` (heap garbage). Root
cause: `jra->Tair` allocated `myDim`, sea-ice thermo loops `myDim+eDim`.

**Rule baked in:** any per-node forcing array defaults to
`(myDim+eDim) * stride`. Any per-node array that *might* be read at halo
by *any* downstream routine is sized `myDim+eDim`. Halo memory cost is
trivial; OOB reads are catastrophic.

### 3.3 Tracer storage stride: `nl`, not `nl-1`

**Pattern:** tracers are conceptually mid-layer values (count `nl-1`), but
storage is `calloc(N * nl, …)`. Index with `nl` stride; the trailing slot
at `nz=nl-1` is unused but is part of the layout. Using `nl-1` stride
silently shifts every read by one slot per vertex.

**Symptom (Phase G6b):** `sigma_xy` 1000× too large → bolus velocity
21 m/s → tracer FCT overshoot to S=59.69 PSU at step 1 → CG NaN at
step 2. Lost a session before the print revealed magnitudes 1000× off.

**Rule:** every per-node 3D field in this port uses stride `nl`. The
arrays explicitly sized `[N * (nl-1) * …]` (e.g. `neutral_slope`,
`slope_tapered`) are the *only* exception. Iteration range and storage
stride are independent — use the allocation shape to determine stride.

### 3.4 Unported producers surface as new consumers land

**Pattern:** an aux/forcing field is allocated and halo-exchanged, but its
Fortran *producer* was never ported. Calloc-zero sits there until a new
phase introduces the first reader. Phase output looks "bit-identical to
baseline" — masquerading as success.

**Real instances:**
- `sw_alpha`/`sw_beta` (oce_ale_pressure_bv.F90:2751-2846) — exposed by
  GM Phase G5/G2b; ported in `0d321b4`.
- JRA halo coverage — exposed by sea-ice Phase D; producer fixed `55f6b6c`.
- `MLD1_ind` — explicitly deferred until G2a needed it.

**Rule:** when you port a phase that introduces a *reader* for `aux->Y`,
grep the C port for any *writer* of `aux->Y`. If none, the producer is a
pre-existing gap; port it first or alongside. Don't trust "alloc + exchange
exist" as evidence the field is populated.

---

## 4. Confusing parts of FESOM2 itself (not specific to C porting)

These are the parts of FESOM2 that *consistently* misled us. They are
sources of bugs even for someone porting in the same language, and they
are the parts a Fortran maintainer should consider documenting better.

### 4.1 Node-vs-cell vs. depth-per-node

**Confusion:** `aux3d.out` lists *one depth per node*. The model geometry
is *cellwise*. Confusing the two leads to: `helem` computed from depth-per-node
(wrong), partial-cell thickness mis-derived, `areasvol` scaling off.

**Truth:**
- `aux3d.out` depth is *input metadata*, used at mesh-prep time to derive
  `elvls.out` (`K_c = nlevels(elem)`).
- After mesh prep, the operative quantities are
  `nlevels(elem)` (cell-based, K_c, from `elvls.out`),
  `nlevels_nod2D(node)` = MAX over surrounding cells (K_v⁺, from `nlvls.out`),
  `nlevels_nod2D_min(node)` = MIN over surrounding cells (K_v⁻).
- Cell thickness `helem(nz, elem) = (1/3) Σ_v hnode(nz, v)` is the mean of
  the three vertex thicknesses — *not* derived from depth-per-node.

(See `feedback_node_vs_cell_geometry.md` and the `geometry.html` page in
the FESOM2 docs.)

### 4.2 Edge length units: `l_e` is in radians, `d_ec` is in metres

`l_e` (edge length) is stored in *radian* angular measure; `d_ec`
(cell-to-edge distance) is stored in *physical* metres. To convert
between them you need the cos(latitude) of the relevant cell's centroid.
This bites because both are "lengths" by name and both appear in the
same edge-flux formulas.

### 4.3 Mesh element orientation must be normalised to CW

Fortran's `test_tri` check at startup swaps element node order to
*clockwise* (~244,654 of 244,659 CORE2 elements get swapped). Without
this, half the stiffness-matrix terms have the wrong sign and SSH
explodes. The orientation check is invisible to a casual reader of
`gen_modules_read_mesh.F90`; it's easy to skip.

### 4.4 PHC IC is in-situ T, model uses potential T

`phc3.0_winter.nc` provides *in-situ* temperature. FESOM2 internally uses
*potential* temperature. Conversion via Bryden 1973 polynomial
(`subroutine insitu2pot` in `oce_ale_pressure_bv.F90:~141`) is required
after horizontal interp + vertical fill, before storing. Skipping is
silent: ~0.1–1 °C deep-water bias. Long-term parity with Fortran is then
gone with no obvious symptom.

### 4.5 Runoff routing depends on which thermodynamics path is active

There are two contradictory conventions in FESOM2 for getting runoff into
`water_flux`:

- **Non-ICEPACK / standard path** (`ice_thermo_oce.F90:therm_ice`):
  `prec = rain + runo + snow*(1-A)` → `fw = prec + evap + fwice + fwsnw`.
  `oce_fluxes` writes `water_flux = -fresh_wa_flux`. **No** `-runoff` term;
  runoff is already inside `fresh_wa_flux`.
- **ICEPACK path** (`ice_oce_coupling.F90:378`):
  `water_flux = -(fresh_wa_flux*inv_rhowat) - runoff`. Runoff is
  subtracted *separately* because ICEPACK's freshwater output doesn't
  include it.

If the port mirrors one path and the auxiliary `fesom_sss_runoff_step`
mirrors the other, runoff is double-counted (or missing). We hit exactly
this when wiring sea ice: pre-ice the C port silently mirrored ICEPACK
semantics (because we had no `therm_ice`), then sea ice landed and we had
to remove the manual subtraction in lock-step (`feedback_runoff_fold_when_ice.md`).

### 4.6 PHC `extrap_nod3D` is non-deterministic across rank counts

The Gauss-Seidel "fill from neighbours" algorithm depends on iteration
order, which depends on partition. Surface T can disagree by O(1 °C) and
S by O(10 PSU) across rank counts at coastal / inland-sea / shallow nodes
(~2,700 of 126,858 on CORE2). **The Fortran does the same thing the same
way.** Same code, same non-determinism — accept it; do not burn cycles
diffing PHC at 1-rank vs N-rank.

The signature: identical step-1 stats *except* small `pgf` and `rs`
deltas; per-node fields differ at exactly the dummy-fill nodes. Rank
groupings appear ({1, 8} vs {16, 32, 144}) because the iteration order
crosses partition-size thresholds. Looks like a port bug; isn't.

### 4.7 `compute_pressure_bv` runs a sequential downward integral

`hpressure` is a column-cumulative integral. With MPI, even when T and S
are *bit-identical* between rank counts, the *order of summation* across
the column differs imperceptibly per FP rounding rules and can produce
~1e-15 differences in `hpressure`. This is *not* a port bug — Fortran
does the same — but it does mean expecting "bit-identical" beyond raw
state vectors is wrong.

### 4.8 Driver order between ice and ocean

Fortran calls ice at the *end* of `oce_timestep_ale`; its outputs feed the
*next* ocean step. That means:

- `oce_fluxes_mom` (ice-aware wind stress) runs at the **top of the next
  step**, consuming `uice` from the previous EVP run.
- `oce_fluxes` (heat/water/virt/relax) is written at the **end of this
  step**, so the next step's tracer TDMA sees the ice-aware values.

We had to move `fesom_ice_step` from inside `fesom_timestep` to
`fesom_main.c` *before* `fesom_timestep` to match this ordering — the
C port was placing ice last and it ran a step out of phase.

### 4.9 The aEVP/mEVP dispatcher is a trap for newcomers

`MOD_ICE.F90:743-754` allocates `uice_aux/vice_aux/alpha_evp_array/
beta_evp_array` *only* when `whichEVP /= 0`. Standard `whichEVP=0` does
not need them, but the standard reader/comprehender of `MOD_ICE` would
allocate them all and burn weeks porting `ice_maEVP.F90`. Our port skips
those entirely with a loud `MPI_Abort` if `whichEVP != 0`.

### 4.10 The `dyngr/dyngrsn/dyngra` CMIP6 fields are output-only

`MOD_ICE.F90:901-904` and `ice_setup_step.F90:203-205,307-309` allocate
and update three diagnostic fields whose only consumer is the CMIP6
output writer. They have *zero* physics impact. Skip them.

### 4.11 KPP comes in two flavours: FESOM1.4 KPP vs CVMix wrapper

`mix_scheme_nmb=1` is the ORIGINAL FESOM KPP (`oce_ale_mixing_kpp.F90`,
1185 lines, self-contained, no CVMix). `mix_scheme_nmb=3` is the CVMix
wrapper (`gen_modules_cvmix_kpp.F90`). They are *different code* with
different lookup tables and OBL recipes. CORE2 uses the original; we
spent time confused which to port before identifying mix_scheme_nmb at
the namelist.

---

## 5. Confusing parts of the *port* (C/MPI/SLURM-specific)

### 5.1 Build environment on Levante (mambaforge poisons everything)

If `mambaforge` or `miniconda` is in `PATH`, its `nc-config`, `mpicc`, and
related tools shadow the spack-loaded gcc-11.2 / openmpi-4.1.2 stack and
silently break MPI and NetCDF detection. We always use `bash -l configure.sh`
which sources `env.sh` and `module --force purge`s the conda variants.

Symptoms: "no MPI_C" CMake error, `nc-config --has-parallel` returning
empty, link failures referencing wrong netcdf paths. If you see this,
the cause is mambaforge in PATH. (See `feedback_levante_build.md`.)

Also: the OpenMPI runtime knobs in `env.sh` assume IB hardware, so on the
login node you must `unset OMPI_MCA_pml OMPI_MCA_btl OMPI_MCA_osc
OMPI_MCA_coll`. Compute nodes (SLURM allocation) keep them.

**Always clean-build** (`bash -l configure.sh --clean`). Incremental
builds across `git checkout` boundaries can link stale `.o` files
together with newly-compiled ones, producing a binary that doesn't match
any source-tree state. Lost a session on this once.

### 5.2 SLURM OUT_DIR clobbering

Two SLURM jobs queued at the same time, both writing to
`/work/.../core2_16r_fix/`, will clobber each other's `run.log`,
`run.err`, `slurm.*`, and `snap_*.nc` files. The headers come from one
job, the per-step lines from another. The result is hours of debugging
phantom crashes.

**Rule:** every SLURM job gets a unique OUT_DIR (encode the experiment in
the name: `core2_16r_phaseD_evp_a`, `core2_16r_g6b`, etc.). For diff
diagnostics, also encode a `STAMP=$(date +%s)_${SLURM_JOB_ID}` so even
re-runs of the same experiment don't collide.

(See `feedback_unique_outdir_per_job.md`.)

### 5.3 Halo exchange semantics: `myList_*[i]` ordering is convention

`myList_nod2D[1..myDim_nod2D]` — interior nodes *this rank owns*.
`myList_nod2D[myDim_nod2D+1 .. myDim_nod2D+eDim_nod2D]` — halo nodes
(owned by another rank, but this rank reads them).

For *element* arrays with full halo, the layout is
`[1..myDim_elem2D]` interior, `[myDim_elem2D+1 .. +eDim_elem2D]` halo,
`[+eDim_elem2D+1 .. +eXDim_elem2D]` *extended* halo. The extended halo
is used by 2-ring viscosity stencils and by inverse CSR
(`nod_in_elem2D`). It must be allocated and populated; do not skip it
even if your immediate kernel doesn't read it.

### 5.4 Halo packing requires explicit pack/unpack buffers

The Fortran halo exchange uses `MPI_TYPE_INDEXED` to scatter directly
into halo positions because `rlist[]` entries from a single sender are
*not* contiguous in the local index space. We use the simpler equivalent:
pack send buffer from `slist`, recv into a separate buffer, then scatter
from `rlist` into `field[]`. Direct `MPI_Irecv` into the halo positions
*does not work* without a custom datatype.

### 5.5 Allreduce summation order is not bit-deterministic across rank counts

`MPI_Allreduce(MPI_SUM)` with default options has implementation-defined
summation order. At 8 ranks we got lucky and matched 1-rank; at 16+ the
tree-reduction order changes and produces ~1e-15 differences per
reduction. Ocean dynamics is chaotic; small perturbations grow.

The fact that {1, 8} and {16, 32, 144} fell into two discrete groups for
`max(rs)` is the smoking-gun signature of FP-summation noise *amplified
by chaos*, not a halo bug.

This *is* the only legitimate source of non-bit-identical behaviour
between rank counts. Tolerances:
- Bit-identical for fields not touched by Allreduce.
- ulp-relative (`< 1e-12` relative) for fields downstream of CG dot
  products and `integrate_nod` sums.

If you need bit-identity for testing, OpenMPI supports
`MPI_REPRODUCIBLE` reductions — slower but deterministic.

### 5.6 Snapshot writer must use rank-0 gather, not parallel HDF5

The Fortran reference does serial gather → rank-0 NetCDF write. We use
**serial netcdf-c** for this reason — there is no need for parallel
HDF5, and trying to use it on Levante with conflicting modules is its
own debugging session. The gather protocol is documented in
`docs/plans/completed/20260424-mpi-port.md` Task 30g.

### 5.7 The Fortran `T_PARTIT` reader vs ASCII partition files

The partition files (`my_list*.out`, `com_info*.out`, `rpart.out`) are
ASCII line-oriented integers. They are read with `fscanf` against
integer tokens. The binary restart format `READ_T_COM_STRUCT` adds
extra fields (`nreq`, etc.) that are *not* in ASCII; don't expect them.
Verified format is documented in the MPI plan.

---

## 6. Process discipline that paid off

These are the working agreements that turned debugging hours into
debugging minutes. Each was hard-earned.

### 6.1 Per-substep dump shim drives validation from line 1

Before any C code was written, we built a Fortran **reference-dump shim**
that records named arrays at probe nodes after each substep of
`oce_timestep_ale`. The shim writes `evp_dump_*`, `phc_dump_*` files
with consistent format; the C port writes the same format (gated by an
env var `FESOM_EVP_DUMP_DIR`).

A diff script (`scripts/evp_dump_diff.py`) reads down the table and
flags the first row above a 1e-9 threshold. This turns "step-2
divergence" debugging from days into ~10 minutes: the offending substep
is named directly. (See `reference_evp_dump_diagnostic.md`.)

The dump infrastructure cost ~2 days to build. It paid back at the very
first multi-rank divergence and has paid back every time since.
**Build this on day 1.**

### 6.2 Always-on cheap stderr probes

`if (ice_strength[el] > 1e7) fprintf(stderr, …)` is always on. It surfaces
the JRA-halo bug class instantly without any env knob. Cost: zero outside
the bug case. Benefit: bug names itself.

Generalise: for any kernel with a known physical magnitude, add a
sanity-threshold print as a permanent guard. Don't gate it; gating
means it's off when you need it most.

### 6.3 Env-knob bisection ladder

Each subsystem has a `FESOM_NO_*=1` knob that disables it. They compose:
- `FESOM_FREEZE_TS=1` (reset T, S to IC each step) — pins tracer-path bugs.
- `FESOM_NO_WIND=1`, `FESOM_NO_HFLUX=1` — zero forcing.
- `FESOM_NO_TRADV=1`, `FESOM_NO_TRDIFF=1` — skip tracer advection/diffusion.
- `FESOM_NO_ICE_DYN`, `FESOM_NO_ICE_ADV`, `FESOM_NO_ICE_THERMO` — sea ice.
- `FESOM_NO_GMREDI=1` — entire GM/Redi block (also a bit-identity off-switch).

Three short SLURM runs typically localise a multi-rank bug to one
subsystem in ~10 minutes. Keep these knobs in the tree forever; they
cost nothing and pay back every regression.

### 6.4 Bit-identity off-switches

Wherever a phase introduces new physics, a master env-switch makes it a
no-op. The acceptance test for the phase includes "with switch on,
output is byte-identical to pre-phase HEAD." This catches gating bugs
(`FESOM_NO_GMREDI` byte-matched the pre-G commit on G9 — a strong
correctness guarantee for the gating logic).

### 6.5 Commit at every working state

Without a git commit anchored at the verified-good run, false-positive
crashes (e.g. SLURM log corruption — see 5.2) are indistinguishable from
real regressions. Lost a working Phase D state to this once. Now:
**commit the moment a multi-rank validation passes**.

### 6.6 Multi-rank validation on every phase end

The end of each phase is "build + halo audit + 16-rank 200-step
validation run + check acceptance criterion". We treat this triplet as a
unit-test suite. No phase starts until the triplet passes.

### 6.7 Plan files live in `docs/plans/`, archived to `docs/plans/completed/`

Each major phase has its own plan with checkboxes that get marked `[x]`
as they land. Archive when complete. The plan is the durable artefact;
the conversation is ephemeral.

### 6.8 Memory store for cross-session knowledge

We use a structured `memory/` directory with `feedback_*`, `project_*`,
`reference_*`, `user_*` types. Each memory file leads with its rule, then
**Why** (the incident), then **How to apply** (the trigger). This works
because every memory is actionable; it is not lore.

(For another porter, the equivalent is a `docs/lessons/` folder with one
file per non-obvious pattern.)

---

## 7. The phase order that worked

Concrete order, retroactively. Use this for the next port.

### Phase 1 — Mesh + serial baseline (1 rank)
- Mesh I/O, areas, edge dx/dy, gradient_sca.
- State allocation, struct shapes mirroring Fortran modules.
- Constant T=10/S=35 IC. Zero forcing.
- linfs only. Skip `zlevel`, `zstar`, `wsplit`, until later.
- JM-EOS + hydrostatic pressure.
- Linear free surface.
- CG SSH solver (single-rank).
- Coriolis (AB2!) + PGF + explicit vertical viscosity + bottom drag.
- Upwind tracer advection (no FCT yet).
- PP mixing (KPP later).
- Should pass: rest-state stays at rest. SSH bump propagates as gravity wave.

### Phase 2 — Pi-mesh stable (1 rank)
- FCT advection (watch `edge_vflux` sign).
- Horizontal viscosity (opt_visc=5 with backscatter).
- `wsplit` (implicit/explicit vertical-velocity split).
- Analytical forcing (constant wind stress, SST restoring).
- 1000 steps at dt=100 must not blow up.

### Phase 3 — CORE2 1-rank
- PHC IC bilinear interp + `extrap_nod3D` + vertical fill + insitu2pot.
- CORE2 mesh specifics: rotation auto-detect, CW orientation, partial cells.
- `zlevel` ALE (still not zstar).
- JRA55 reader (bilinear interp + L&Y09 bulk formulae).
- SSS restoring + runoff.
- 172 steps at dt=500 must be stable.

### Phase 4 — MPI (CORE2 16+ ranks)
- `fesom_partit` reader (rpart, my_list, com_info — ASCII).
- Generic halo exchange with explicit pack/unpack buffers.
- Halo identity test (fill interior with rank-id, exchange, halo == owner-rank-id).
- All `_alloc()` use `myDim+eDim` (or `+eXDim` for full halo elements).
- Audit every loop — `mesh->nod2D` → `myDim_nod2D` etc.
- Halo exchange after every kernel that writes a field other ranks read.
- Parallel CG (Allreduce dot products, exchange before SpMV).
- `integrate_nod` global reductions in SSS/runoff.
- FCT mid-pipeline exchange of `fct_plus`/`fct_minus`.
- Snapshot writer = gather-to-rank-0 + serial NetCDF.

**This is where halo bugs hide.** Audit every C write-loop against
Fortran *as you write*, not after divergence appears.

### Phase 5 — Sea ice (the missing physics for stability)
- Phase A: scaffold (struct, allocator, no-op driver, snapshot fields).
- Phase B: thermodynamics (`cut_off`, `obudget`, `budget`, `therm_ice`,
  `thermodynamics`). Folds runoff into `flx_fw`.
- Phase C: ice→ocean coupling (heat/water/virt/relax fluxes). **Removes**
  the runoff subtraction in `fesom_sss_runoff_step` (lock-step).
- Phase D: EVP dynamics (`stress_tensor`, `stress2rhs`, `EVPdynamics`).
  120-subcycle solver. Halo discipline critical.
- Phase D5: `oce_fluxes_mom` (ice-aware surface stress).
- Phase E: FCT advection of (a_ice, m_ice, m_snow). Highest halo-bug risk.

The whole point of sea ice is to fix `T.min → −27 °C` drift in long
runs, by insulating the high-latitude ocean. Acceptance: T.min stays
above ≈ −2 °C in 256-rank model month.

### Phase 6 — GM/Redi
- G0: mesh prereqs (`mesh_resolution`, `nlevels_nod2D_min`,
  `ulevels_nod2D_max`, `zbar_3d_n`).
- G1: state allocation.
- G2a: `MLD1_ind` in pressure_bv (unconditional).
- G2b: `compute_sigma_xy`, `compute_neutral_slope`.
- G3: `init_Redi_GM`.
- G4: `fer_solve_Gamma`, `fer_gamma2vel`.
- G5: Redi K33 in vertical-implicit tracer diffusion.
- G6a: `fer_w` accumulator inside `vert_vel_linfs` (highest risk —
  invasive mod to validated kernel; bit-identity gate required).
- G6b: bolus add/sub around tracer block.
- G7a: vertical-explicit Redi.
- G7b: horizontal Redi (5 partial-cell branches).
- G8: env knobs + driver wiring.
- G9: validation (off-switch must byte-match pre-G HEAD).

Note for next time: the divergence balance between `fer_uv` and `fer_w`
is **continuity-preserving** (built from the same `fer_solve_Gamma`).
Never clamp them per-component. Scale `fer_gamma` instead. (See
`feedback_bolus_divergence_balance.md`.)

### Phase 7 — Still pending
- Full FESOM1.4 KPP (deferred — start with PP, validate, then upgrade).
- `zstar` ALE.
- Restart I/O.
- Multi-month / multi-year qualitative comparison vs. Fortran reference.

---

## 8. Validation infrastructure that earned its keep

These are the artefacts in tree that you should *not* delete; they pay
back at every regression.

### 8.1 EVP dump diagnostic (`fesom_ice_evp.c`, `scripts/evp_dump_diff.py`)
Per-substep dumps of (a_ice, m_ice, m_snow, srfoce_*, u_ice, v_ice, sigma11/12/22,
ice_strength, rhs, …) keyed by 1-based gid, owner-only. Env-gated by
`FESOM_EVP_DUMP_DIR`. Comparison script flags first row with relative
diff above 1e-9. This is the tool that found the JRA halo bug in
~10 minutes after a 2-day suspicion.

### 8.2 Always-on `HUGE ice_strength` stderr trace
Cheap, surfaces the array-size-vs-reader-loop bug class instantly.
Permanent (no env gate). Leave it on.

### 8.3 Env-knob bisection ladder (Section 6.3 above)
`job_core2_8_freeze`, `job_core2_8_nowind`, `job_core2_8_nohflux`,
`job_core2_8_notradv`, `job_core2_8_notrdiff` — short SLURM runs that
disable a single subsystem each. Three runs typically localise.

### 8.4 Job-script templates encoding rank+phase
`job_core2_{1,8,16,32}_{phase}` with unique OUT_DIRs and
`STAMP`-prefixed dump dirs. Each is committable; rerunning gives
matching dirs.

### 8.5 PHC post-load and pre-extrap dumps
`phc_dump_postload_rank{R}.txt`, `phc_dump_preextrap_rank{R}.txt` —
surface T,S, bilin_i/j, lon, lat. Lets us decide whether a multi-rank
delta is *PHC IC noise* (accepted) or *real* (port bug).

### 8.6 Plan files with `[x]` checkboxes
Each phase has a checked plan in `docs/plans/`. The diff between
"plan submitted for review" and "plan executed" is itself instructive
for the next porter — what was missed at design time, what was added
on the fly.

---

## 9. Things tried and failed (don't redo)

From `FRESH_START.md §21` and our own logbook:

- **Linearising α/β** (constant thermal/haline expansion) → broke OBL
  depth (stuck at 5 m), wrong density gradient. Use full JM-EOS.
- **Adding 500 m cap on OBL depth** ("safety bound") → not in Fortran,
  caused wrong physics. Don't add bounds Fortran lacks.
- **Simplifying `ghats`** to a closed-form shape function → wrong
  magnitudes, didn't fix mosaic. Port the lookup table.
- **Applying horizontal diffusion directly to T then dividing by hnode** →
  double-scales the diffusion tendency. Must accumulate into `del_ttf`
  first, then ALE-reconstruct.
- **Using `dt/areasvol` instead of `dt*area/areasvol`** in TDMA w_i terms
  → imprinted mesh geometry on tracer evolution. Even when ratio = 1
  for zlevel, write the full expression.
- **Starting with KPP** → too complex to debug alongside other physics.
  Start with PP, then upgrade.
- **CVMix KPP wrapper instead of FESOM1.4 KPP** → CORE2 uses the
  ORIGINAL FESOM KPP (`oce_ale_mixing_kpp.F90`, `mix_scheme_nmb=1`),
  not the CVMix wrapper.
- **AB_order=2 with AB_order=3 history-shift** → Coriolis weight 2.5
  instead of 1.0 → SSH runaway. Match the AB-history slot count to
  AB_order *exactly*.
- **Per-component clamps on `fer_uv` / `fer_w`** during GM debugging →
  breaks divergence (continuity), tracer FCT overshoots even at small
  magnitudes. Scale `fer_gamma` uniformly instead.

---

## 10. Suggestions for FESOM2 Fortran-side changes that would make the next port faster

These are *suggestions* to the Fortran maintainers, not requirements for
the port. Each addresses a confusion that consumed real debugging time.
None are physics changes; all are organisational / API / documentation.

### 10.1 Document loop-bound conventions per write-loop

Every Fortran write-loop of the form
`do n=1, partit%myDim_nod2D + partit%eDim_nod2D` is *intentionally* filling
halo entries to avoid a downstream `call exchange_*`. This is invisible
to a porter or maintainer reading the routine — there is no comment,
the bound looks unusual but not load-bearing. It *is* load-bearing.

**Suggestion:** add a one-line comment to each such loop, e.g.

```fortran
do n=1, partit%myDim_nod2D + partit%eDim_nod2D   ! halo-cover write — replaces post-loop exchange
   ...
```

This single discipline change would have eliminated the three halo-write
bugs we found (`init_tracers_AB`, `fct_ttf_max/min`, PHC + SSS/runoff).

### 10.2 Standardise allocation comments

When `T_ICE.data(1:3) is allocated `(myDim+eDim) * stride`, the comment
at allocation should declare the reader contract:

```fortran
allocate(ice%data(t)%values(myDim_nod2D + eDim_nod2D)) ! read at halo by ice_thermo (myDim+eDim loop)
```

Knowing *who reads it at what bound* at the allocation site catches the
"sized myDim only, read at halo" bug class at code review, not at
runtime.

### 10.3 Single source of truth for partition-file format

`gen_modules_partitioning.F90` documents the format implicitly via the
reader. The ASCII format is identical across all FESOM versions but only
specified in code. A short standalone document
(`docs/PARTITION_FILE_FORMAT.md`) listing the file layout, units, and
sentinel values would have saved a day on the C side.

### 10.4 Distinguish in-situ vs. potential T at the API boundary

`subroutine read_ic_tracer` accepts a NetCDF and a tracer index but does
not encode (in the type system or in a comment) whether the input is
in-situ or potential. The PHC-as-in-situ-T trap (silent ~0.1–1 °C bias)
is structural: any porter who skips the `insitu2pot` step gets the same
bug. A wrapper `read_in_situ_temperature` vs `read_potential_temperature`
would make the conversion mandatory.

### 10.5 Standardise the runoff-routing convention

The two contradictory paths (ICEPACK vs non-ICEPACK) for getting runoff
into `water_flux` are inviting double-counting. A single helper
`apply_runoff_to_water_flux(mode)` that does the right thing for the
selected ice scheme would eliminate the `fesom_sss_runoff_step` ↔
`therm_ice` lock-step contract entirely.

### 10.6 Unify the two KPP implementations

The original FESOM KPP (`oce_ale_mixing_kpp.F90`) and the CVMix wrapper
(`gen_modules_cvmix_kpp.F90`) coexist with different tables and OBL
recipes. New porters / new users routinely confuse them; two mix schemes
selected by `mix_scheme_nmb` is a lurking footgun. Pick one (CVMix,
modern) and deprecate the other; or rename the original to make the
distinction loud (`oce_kpp_fesom14.F90`).

### 10.7 Mark output-only diagnostics

`MOD_ICE.F90:901-904` allocates `dyngr/dyngrsn/dyngra` for CMIP6 output
*only*. There is no annotation distinguishing "physics state" from
"output diagnostic". A simple comment block (or a section header
`!=== output-only diagnostics ===`) would let porters skip them
confidently.

### 10.8 Replace `extrap_nod3D` with a deterministic algorithm

The Gauss-Seidel fill is *partition-dependent by construction* — an
`O(1 °C)` per-rank-count IC noise on coastal nodes. In a port, it
masquerades as a port bug for hours. A deterministic alternative
(distance-weighted ANN fill, e.g.) would eliminate the rank-grouping
signature that misleads porters. The Fortran has the same artefact;
fixing it Fortran-side benefits both maintenance paths.

### 10.9 Document the `test_tri` orientation normalisation

The fact that ~all CORE2 elements get swapped to CW orientation at
startup is *not visible* in `gen_modules_read_mesh.F90` to a casual
reader; the orientation check is buried. Adding a section header
`!=== element-orientation normalisation (REQUIRED for stiffness sign) ===`
would make it impossible to skip when reimplementing.

### 10.10 Make `whichEVP=0` the only path, or split files

`MOD_ICE.F90:743-754` with its conditional aux allocations,
`ice_maEVP.F90` (mEVP/aEVP), and the `whichEVP` dispatcher all live
together. A porter has to discover that `whichEVP=0` is the only
production path. Splitting the standard-EVP code into its own file
(without the aEVP/mEVP entanglement) would make the supported subset
obvious.

### 10.11 Document edge-length-units once, prominently

`l_e` is in radians; `d_ec` is in metres. Stating this once at the top
of `gen_modules_read_mesh.F90` (where the edge fields are read) and once
at the top of any module that uses them (e.g. flux computation in
`oce_dyn.F90`) would eliminate the unit-confusion bug class.

### 10.12 Make the ice/ocean call ordering explicit

The fact that `oce_fluxes_mom` consumes the *previous* step's `uice` and
that ice runs *at the end* of `oce_timestep_ale` is a constraint on the
driver order. A short ASCII comment block at the top of the driver
showing the temporal flow:

```
!  step n:    ocean dynamics  →  oce_fluxes_mom (using uice from step n-1)
!             ...
!             ice step          →  uice^n stored
!  step n+1:  ocean dynamics  →  oce_fluxes_mom (using uice from step n)
```

would have saved us the move of `fesom_ice_step` from inside
`fesom_timestep` to before it.

### 10.13 Prefer `MPI_REPRODUCIBLE` reductions in solver paths

The CG dot products use `MPI_Allreduce(MPI_SUM)` with default ordering.
Switching to `MPI_REPRODUCIBLE` eliminates rank-count non-determinism in
SSH (and therefore everything downstream of SSH) at a small performance
cost. Worth offering as a runtime flag.

---

## 11. A starter checklist for day 1 of the next attempt

A literal starter for the next porter, distilled from this report.

### Day 1
- [ ] Read `FRESH_START.md` cover to cover. Take notes.
- [ ] Read this file (Sections 2, 3, 4 minimum).
- [ ] Build the Fortran reference under `/home/a/a270088/fesom27/fesom2/work_core/`
      and run a 1-day CORE2 reference. Save snapshots.
- [ ] Build the dump shim *first*. The C code does not exist yet; you need
      to validate against named arrays at probe nodes after each substep
      of `oce_timestep_ale` from the very first kernel.
- [ ] Decide your struct mapping: one C file per Fortran file, mirror
      subroutine names with `fesom_` prefix. Document this in the README.
- [ ] Set up `bash -l configure.sh --clean` build script that strips
      mambaforge.

### Phase 1 (mesh + serial baseline)
- [ ] Mesh I/O. Print nod2D/elem2D/nl. CORE2 must read 126,858 / 244,659 / 48.
- [ ] Element-orientation normalisation (CW). ~244,654 of 244,659 swap.
- [ ] State allocation per Section 7 above.
- [ ] Constant-T/S IC, zero forcing, linfs only, JM-EOS, hydrostatic
      pressure, CG SSH solver, AB2 Coriolis with 1-slot history.
- [ ] Acceptance: rest-state stays at rest. SSH bump propagates as
      gravity wave at sqrt(gH).

### Halo bug class — proactive audits
- [ ] `grep -nE 'for.*<.*myDim_(nod|elem)2D;' src/*.c | grep -v eDim`
      after every new C file lands. For each match, open the matching
      Fortran subroutine and check the bound.
- [ ] Document each audited loop with an inline comment citing the
      Fortran line and the agreed bound:
      `/* matches oce_dyn.F90:NNN — Fortran loop bound is myDim_nod2D */`

### Validation rituals
- [ ] Build + halo audit + 16-rank 200-step validation run is the
      "test suite" at every phase end. No phase starts until the triplet
      passes.
- [ ] Commit at every working state. Without commits, false-positive
      crashes from log corruption are indistinguishable from real
      regressions.
- [ ] Always-on stderr probes for any kernel with a known physical
      magnitude (e.g. `if (m_ice[n]>1e3) fprintf …`). Permanent. No
      env gate.

### Forever-on env knobs
- [ ] Each phase introduces a `FESOM_NO_*=1` knob disabling its physics.
      They compose. Do not delete them after the phase ships.
- [ ] Bit-identity off-switch (e.g. `FESOM_NO_GMREDI=1`) byte-matches the
      pre-phase HEAD on a 1-rank smoke. Make this part of the
      acceptance test for every phase that adds a major subsystem.

### Don'ts
- [ ] Don't simplify Fortran code. Don't linearise α/β. Don't cap OBL
      depth. Don't replace lookup tables with closed forms.
- [ ] Don't trim per-node arrays to `myDim` to "save memory". Default
      to `myDim+eDim`.
- [ ] Don't use mambaforge. Always `bash -l configure.sh`.
- [ ] Don't share OUT_DIR across SLURM jobs.
- [ ] Don't ignore PHC IC noise across rank counts. It's accepted
      Fortran behaviour. (Burned-in: the smell is multi-rank diffs at
      coastal/inland-sea/shallow nodes only.)
- [ ] Don't start with KPP. Start with PP. Don't start with sea ice.
      Get the ocean stable first.

---

## End notes

This port is at HEAD `30b59b3` (2026-04-25). 15,300 LoC C against
73,000 LoC Fortran. The compression ratio mostly comes from collapsing
exchange variants (50+ Fortran shape-specific exchange routines → 4
generic C functions) and from skipping inactive code paths (mEVP/aEVP,
ICEPACK, cavity, icebergs, `__oifs`/`__yac`/`__oasis`, multiple ALE modes,
CVMix wrappers, the unported-by-design diagnostic CMIP6 fields).

Three roadmap items remain: full FESOM1.4 KPP, `zstar` ALE, and restart
I/O. The strategy for each is the same as the strategy that worked for
ocean / MPI / sea-ice / GM-Redi: find the matching Fortran source, build
a per-substep dump shim, port line-by-line, audit halo writes, validate
1/8/16/32-rank with bit-identity off-switches at the phase boundary.

The rule that this port lives by — *trust the Fortran, port line-by-line,
ease-of-porting is the criterion* — should be the first thing the next
porter writes on the wall.
