# Restart I/O for the C port (and then for Kokkos)

> **Revision 3, 2026-08-18, gate green.** §2 was written from a static reading of the code and was
> incomplete: running the gate found four further pieces of persistent state (§2b) and then one
> thing that is not state at all — elements on a partition boundary are computed redundantly on
> several ranks and their copies drift apart at roundoff, so no partition-independent file can
> reproduce an uninterrupted run (§2c). §5 restates the gate accordingly and §7 records where it
> now stands: exact at four (N, M) pairs.

Paper plan: `~/port_paper/docs/plans/20260818-port-paper-v2-63yr-and-optimisation.md`, Task 2.

`io_restart.F90` was a deliberate non-port. Without it a 63-year hindcast has to fit inside one
12 h job, which is why the two earlier attempts stopped at 42 and 53 years. This document
enumerates the state a restart has to carry, states where the port deliberately departs from the
Fortran, and defines the gate. It is written before the code, because the interesting part of this
job is the enumeration and not the netCDF plumbing.

## 1. What the Fortran actually writes

Read off the reference run itself (`/work/ab0995/a270088/fesom2_core2/fesom.1958.*.restart/`), not
off the source, so the list is the one our configuration produces — `cvmix_TKE` + zstar + mEVP,
`AB_order = 2`, two tracers, no ICEPACK, no icebergs, no isotopes.

    fesom.1958.oce.restart/   18 files, 11 GB     dims (time: 12, node: 126858, nz_1: 47)
      ssh  hbar  ssh_rhs_old  hnode
      u  v  urhs_AB  vrhs_AB
      temp  temp_AB  temp_M1     salt  salt_AB  salt_M1
      tke
      w  w_expl  w_impl

    fesom.1958.ice.restart/    5 files, 59 MB
      area  hice  hsnow  uice  vice

    fesom.clock                 the model clock

`time: 12` is the monthly restart cadence kept inside the year, not a physical dimension.

The registration lives in `ini_ocean_io` / `ini_ice_io` (`io_restart.F90:109` and `:236`). Fields
registered there but absent from our run, and therefore out of scope: `urhs_AB3`/`vrhs_AB3`
(`AB_order == 3`), `iwe` (IDEMIX), `uke`/`uke_rhs` (`opt_visc == 8`), `temp_M2`/`salt_M2`, the
passive and isotope tracers, and the `__oifs`/`__oasis` ice fields.

Mapping onto the C port's variables:

| restart name | C port | shape |
|---|---|---|
| `ssh` | `dyn.eta_n` | node |
| `hbar` | `mesh.hbar` | node |
| `ssh_rhs_old` | `dyn.ssh_rhs_old` | node |
| `hnode` | `mesh.hnode` | node × nl |
| `u`, `v` | `dyn.uv` (comp 0, 1) | elem × nl |
| `urhs_AB`, `vrhs_AB` | `dyn.uv_rhsAB` (comp 0, 1) | elem × nl |
| `temp`, `salt` | `tracers.data[T\|S].values` | node × nl |
| `temp_AB`, `salt_AB` | `tracers.data[T\|S].valuesAB` | node × nl |
| `temp_M1`, `salt_M1` | `tracers.data[T\|S].valuesold` | node × nl |
| `tke` | `tke.tke` | node × nl |
| `w`, `w_expl`, `w_impl` | `dyn.w`, `dyn.w_e`, `dyn.w_i` | node × nl |
| `area`, `hice`, `hsnow` | `ice.data[AICE\|MICE\|MSNOW].values` | node |
| `uice`, `vice` | `ice.uice`, `ice.vice` | node |

## 2. What the Fortran restart leaves out — and why that matters here

**The Fortran restart is lossy: a restarted run is not the bitwise continuation of an
uninterrupted one.** Two pieces of state survive a timestep in memory and are never written:

- **The EVP stress tensor.** `sigma11`, `sigma12`, `sigma22` are updated incrementally inside the
  subcycling loop and are *not* zeroed on entry — `ice_maEVP.F90:752` reads `sigma11(el)` on the
  right-hand side, and the C port carries the same behaviour with the comment
  "sigma11/12/22 NOT zeroed on entry" (`fesom_ice_maevp.c:25`). They are elements, so the loss is
  small and local, but it is not zero.
- **The ice skin temperature.** `ice.thermo.t_skin` is read at `fesom_ice_thermo.c:470` before it is
  written at `:508` — a genuine read-before-write across steps.

This is not a defect we introduced and it is not a bug in the reference; FESOM restarts have never
claimed bit-reproducibility. But the paper's Task 2 gate demands a round trip with a maximum
absolute difference of **0**, and a literal port of the field list above cannot deliver that.

**Decision — the port's restart is a strict superset.** We write the Fortran list plus every
remaining piece of persistent state, so that a restart is an exact continuation. The deviation from
"literal port" is deliberate, is recorded here, and is itself a result worth a sentence in the
paper: enumerating the state to port the restart is what exposed that the original's restart drops
some of it.

Added to the file, beyond the Fortran list:

| field | C port | shape | why |
|---|---|---|---|
| `sigma11`, `sigma12`, `sigma22` | `ice.work.sigma*` | elem | EVP stress, carried between steps |
| `t_skin` | `ice.thermo.t_skin` | node | read before write in the thermodynamics |
| `hbar_old` | `mesh.hbar_old` | node | so `dhe` reproduces exactly instead of `restart_thickness_ale`'s first-step approximation |

## 2b. What the gate found that reading the code did not

Each of these was found by running the round trip and bisecting, not by inspection. Every one is
state that a step consumes before it produces it, and none is in the Fortran restart. They are
listed in the order they were found, which is also the order of how much they mattered.

1. **The Adams–Bashforth cold-start branch.** `fesom_step.c` passed
   `is_first_step = (step_n == 1)` to `fesom_compute_vel_rhs`, which uses `ff_step = 1.0` instead of
   `ab2` on the first step. A resumed run starts again at `step_n == 1`, so it took the cold-start
   branch a second time. Fortran guards exactly this with `if (lfirst .and. (.not. r_restart))`
   (`oce_ale_vel_rhs.F90:287`). Fixed by adding `restarted` to `fesom_step_ctx`.
   Effect: SSH differed by 1.2e-1 m after 20 resumed steps; the fix took it to 1.7e-3 m.

2. **The ALE free-surface stiffness matrix.** `fesom_update_stiff_mat_ale` adds one `dhe`-weighted
   increment per step into `S->values`, so the matrix is the integral of the run. A resumed run
   that rebuilds the base matrix and applies only the last increment has skipped all the earlier
   ones. The Fortran covers this by having `restart_thickness_ale` set `dhe` to the *absolute*
   elemental sea surface height for the first step after a restart, which telescopes to the right
   matrix in exact arithmetic but not in floating point. We store the matrix instead: per owned row,
   columns sorted by global node id, so the file stays partition-independent. The column ids travel
   with the values and are checked on read.
   Note the coupling: because the matrix is restored exactly, `dhe` must be saved as the true
   increment. Saving the matrix and using Fortran's absolute-`dhe` trick together would double-count.

3. **`dyn->d_eta`.** `fesom_ssh_solve_cg` takes `dyn->d_eta` as the CG initial guess and never zeroes
   it, so the solve warm-starts from the previous step. Starting from zero lands on a different
   point inside the solver tolerance — a difference of order the tolerance, every step.

4. **`dyn->uvnode`.** The subtlest one. It is written in the middle of the step (`fesom_step.c:136`)
   but read at the *top* of the next iteration by the bulk formula (`fesom_bulk.c:261`), so at the
   moment a restart is written it holds the velocity from the step *before* the one being saved.
   Recomputing it from the saved `uv` gives a value one step too new, and the surface stress — hence
   TKE, hence everything — is wrong from the first resumed step. The driver's pre-loop
   `fesom_compute_vel_nodes` call is therefore skipped on a restart.

**The preconditioner is not on this list, and that is worth recording.** `fesom_ssh_preconditioner`
is called once at initialisation (`fesom_main.c:516`), before any `dhe` update, so both a cold start
and a resumed run build it from the base matrix and no extra state is needed. Had it been built
lazily at the first solve, as the Fortran comment at `fesom_ssh.c:233` describes, it would have had
to be saved too.

## 2c. The thing that is not state: replicated elements

The four items above are all *state*, and adding them to the file removed most of the difference.
What was left was 3.6e-12 on `Z_3d_n` after a single resumed step, and no missing field explains it,
because the last piece is not a field.

`gen_comm.F90:264-270` puts an element into `myList_elem2D` of **every** rank that owns one of its
nodes, and each of those ranks computes it. For a quantity assembled from an element's own three
nodes that is harmless: same inputs, same arithmetic, same bits. It is not harmless for a quantity
assembled by an edge loop. `visc_filt_bcksct` sums into `U_b` over `myDim_edge2D + eDim_edge2D`
(`oce_dyn.F90:330-360`, mirrored in `fesom_momentum.c:535-568`), and two ranks reach the same
element's three edges in a different order, so they get the same sum with different rounding.
`exchange_elem` cannot repair it — the element is *owned* on both ranks, so it is in neither one's
halo list, and a halo exchange never touches it.

Measured on CORE2 at 128 ranks with `FESOM_IO_GATHER_DUPCHECK=1`, which counts, per gathered field,
how many global element ids are written more than once and how many of those writes disagree:

| after step | replicated slots | `uv` disagreeing | max &#124;Δ&#124; | `uv_rhsAB` disagreeing | max &#124;Δ&#124; |
|---|---|---|---|---|---|
| 1  | 15348 | 0    | 0      | 0    | 0 |
| 2  | 15348 | 167  | 6.9e-18 | 0    | 0 |
| 5  | 15348 | 688  | 1.1e-16 | 489  | 1.8e-12 |
| 20 | 15348 | 3311 | 4.2e-16 | 2994 | 2.9e-11 |

This is inherited from the Fortran, not introduced by the port: the Fortran loop has the same shape
and the same missing reduction.

The consequence is that **the state of a running model is not a function of the global fields**. It
also includes which rank happens to hold which rounding of those 15348 slots, and a
partition-independent file cannot record that. "A resumed run reproduces a run that never stopped"
is therefore unattainable, and no amount of extra fields would have reached it. It also explains why
the first gates passed at N = 1 and failed at N = 20: after one step the copies still agree.

What *is* attainable is the case a checkpointed hindcast actually needs. Writing a restart pushes the
gathered element values back into the live arrays (`rst_canonicalise`, `fesom_io_restart.c`), so the
writer leaves the model in exactly the state the reader will produce, and every owner of a replicated
element agrees from then on. Cost: one extra scatter per element field per checkpoint, and a
perturbation of one ulp. `FESOM_RESTART_CANONICALIZE=0` turns it off.

### What one ulp costs, measured

`tests/restart_roundtrip/job_restart_cadence` runs 200 steps three ways: checkpointing every 10
(P), checkpointing only at the end (Q), and checkpointing every 10 with the canonicalisation off
(R). **R and Q are bit-identical**, so writing a file changes nothing by itself and every P–Q
difference is the canonicalisation and nothing else. Comparing P against R at each of the 20
checkpoints gives what one ulp buys:

| step | `temp` max &#124;Δ&#124; | nodes differing (of 6.09 M) | `u` max &#124;Δ&#124; |
|---|---|---|---|
| 10  | 0        | 0         | 0        |
| 20  | 1.3e-05  | 3 641 408 | 3.5e-05  |
| 40  | 6.5e-04  | 3 698 823 | 4.3e-03  |
| 100 | 5.5e-03  | 3 705 549 | 3.8e-02  |
| 200 | 4.6e-01  | 3 705 863 | 4.9e-02  |

⚠️ **Read that carefully, because the obvious reading is wrong.** Essentially the whole wet domain
differs ten steps after the first kick, and the maximum then wanders around 1e-2–1e-1 instead of
growing exponentially. That is not chaos spreading from 15348 boundary elements. It is the SSH
solve: `FESOM_PHASE1_SOLTOL = 1e-5` (the Fortran's `solver.F90` value), so the CG is converged only
to 1e-5 relative and is global, and *any* perturbation makes it land on a different point inside its
own tolerance, everywhere, in one step. A one-ulp kick therefore produces a difference of order the
solver tolerance immediately, not eventually. The same is true of the Fortran and of any other
perturbation of the same kind.

The flip side is what makes the gate worth having: an **exact** round trip means the CG took a
bitwise identical path at every step, which happens only if the whole state — `d_eta` included — was
restored exactly. There is no "close enough" version of this test.

Two things follow that have to be said out loud rather than buried:

- **The trajectory depends on the checkpoint cadence**, and at the level of the SSH solver tolerance
  rather than of roundoff. Every production run in the campaign must use the same cadence, and the
  cadence belongs in the run manifest.
- **Turning the canonicalisation off does not avoid the cost, it only moves it.** Without it the
  *write* perturbs nothing but the *resume* does, by the same mechanism and the same amount — that
  is what the 3.6e-12 on `Z_3d_n` at the top of this section was. Production always resumes, so the
  number of kicks is the same either way. What changes is that with it on there is a reference run
  the chain is bit-identical to, and a gate that can find a missing field; with it off there is
  neither. That is why it is the default.
- **The alternative is open, not rejected.** Making the model itself keep replicated elements
  consistent — an owner-wins broadcast after every update of `uv` and `uv_rhsAB` — would remove the
  dependence and make the element state genuinely partition-consistent. It is a per-step
  communication and a change to a validated model's arithmetic everywhere rather than only at
  checkpoints, so it is a decision about the model and not about restart I/O. Not taken here; see
  §8.4.

## 3. What is *not* saved, and the argument for each

Everything below is rebuilt from the saved state, the mesh, or the calendar. Each was checked, not
assumed.

- `aux.*` — `density_m_rho0`, `hpressure`, `bvfreq`, `sw_alpha`, `sw_beta`, `Kv`, `Av`, `dbsfc`,
  `pgf_x/y`, `MLD1_ind`: all diagnosed from T, S and `uv` at the top of the step.
- `tke.tke_Av`, `tke.tke_Kv` — the header states prior-step values are never read; only `tke` itself
  is prognostic, which is exactly why the Fortran restart carries `tke` and nothing else from TKE.
- `gm.*` — slopes, tapers and `fer_*` are rebuilt each step from the density field.
- `forcing.*` — the fluxes are recomputed by `fesom_bulk_compute` from JRA55 and the ice state;
  `Ssurf`, `runoff` and `chl` are functions of the calendar month.
- `mesh.hnode_new`, `mesh.helem`, `mesh.zbar_3d_n`, `mesh.Z_3d_n`, `mesh.dhe` — re-derived on read,
  following `restart_thickness_ale` (`oce_ale.F90:1449`): bottom layer to `bottom_node_thickness`,
  `zbar_3d_n` and `Z_3d_n` cumulated upward from `hnode`, `helem` as the mean over the three
  vertices, bottom to `bottom_elem_thickness`. We set `dhe` from the saved `hbar`/`hbar_old`
  rather than from `hbar` alone, which is where we improve on the Fortran.
- `ice.data[*].values_old`, `ice.uice_old`, `ice.vice_old` — written before read within the step;
  under mEVP the `*_old` velocities are never touched at all.
- `ice.thermo.thdgr_old` — written every step, never read in the C port (a diagnostic in Fortran).
- The JRA55 reader and the SSS/runoff reader — both seek by calendar, so a restart at an arbitrary
  date reloads the same records. **This must be tested, not asserted** (see the gate).
- The output accumulators in `fesom_io` — a restart is written on a calendar boundary, where the
  monthly accumulators are empty by construction. The writer refuses any other cadence.

## 4. Format

One **CDF-5** file per restart, `fesom.<YYYYMMDD>_<seconds-of-day>.restart.nc`, written through the
existing gather path and read through a new scatter path, so the file is partition-independent: a run on 128
ranks can restart on 864, and a restart written by the CPU build can be read by the CUDA build.
That is the point of not using the per-rank raw format that `write_all_raw_restarts` provides.

Not netCDF-4/HDF5: the file is a flat set of large 1-D and 2-D double arrays needing no groups,
chunking or compression, and the HDF5 layer failed with `NetCDF: HDF error` on the second 1.2 GB
write inside one run. CDF-5 carries 64-bit sizes and `NC_INT64` attributes, which is all this format
uses.

Dimensions `node`, `elem`, `nz` (= nl), `stiff_nnz`; every field `double`. All `nl` levels of a 3-D
field are written, including the slots below the bottom, so the round trip is exact by construction
rather than by an argument about which slots are live. The clock and the
provenance ride as global attributes: `restart_format_version`, model year/month/day/second, step
index, `dt`, git SHA, and the knob string.

Env-gated, off by default, so an existing run is unchanged:

    FESOM_RESTART_OUT=<dir>     write restarts into <dir>
    FESOM_RESTART_EVERY=<n>     write every n steps (0 = only at the end)
    FESOM_RESTART_IN=<file>     start from this file instead of the PHC initial condition

## 5. The gate

**Round trip, maximum absolute difference 0.** Run N + M steps continuously, checkpointing at N and
at N + M. Separately run N steps, write a restart, start a second executable from it, run M more,
write a restart. The two final restart files must agree bit for bit, field by field. The restart
file is the right object to compare because it *is* the model's persistent state — a snapshot
comparison would only cover the subset that gets written to the output streams.

The uninterrupted leg checkpoints at N as well, and that is load-bearing, not cosmetic: for the
reason in §2c the comparison against a leg that never wrote anything at step N is not attainable by
any file. The gate as stated is the operationally meaningful one — a 63-year hindcast checkpoints
anyway, and what has to be true is that resuming from a checkpoint continues the run that wrote it.
`FESOM_RESTART_AT=<n1,n2,...>` exists for this: a fixed interval cannot express "at N and at N + M".

Then, in order:

1. **Serial byte gate for Kokkos.** Kokkos-Serial restarted must be byte-identical to the C
   reference restarted, the same way the two are already certified without restarts.
2. **Cross-backend read.** A restart written by the CPU build read by the CUDA build, and the
   reverse. If a difference appears it is stated as a restriction rather than papered over.
3. **File-format unit tests** on fixtures: round trip, a truncated file, a wrong field count, a
   wrong node count, a missing field.
4. **Calendar continuity**: the run restarted at a year boundary must load the same JRA55 records
   and the same monthly SSS as the uninterrupted run at that step.

## 6. Order of work

1. `fesom_io_restart.{c,h}` — write, then read, then the post-read re-derivation.
2. Wire into `fesom_main.c`; keep every existing default unchanged.
3. Unit tests (`tests/test_restart_format.c`, run by CTest).
4. The round-trip gate on core2 at a small rank count, N = M = 20.
5. Only then the Kokkos port.

## 7. Where the gate stands (2026-08-18)

**Green.**

**Control first.** `tests/restart_roundtrip/job_restart_control` runs the identical configuration
twice, no restart anywhere, and the two restart files agree bit for bit. The model is reproducible
run to run at 128 ranks, so any round-trip difference is real and attributable to the restart.

| case | result |
|---|---|
| unit tests, `test_restart_format` | 89 checks green (round trip, component order, optional TKE, cadence incl. `FESOM_RESTART_AT`, four guard cases) |
| control: the same run twice, no restart | exact |
| write twice in one run vs a second run's file at the same step | exact |
| round trip, N = 1, M = 20 | exact, all 38 fields |
| round trip, N = 20, M = 1 | exact |
| round trip, N = 20, M = 20 | exact |
| round trip, N = 137, M = 91 | exact |
| round trip, N = 200, M = 200 | exact |
| probe: owned state at the top of the first resumed step | identical |
| probe: owned state at entry to the timestep | identical |
| cadence: 200 steps checkpointing every 10, canonicalisation off, vs checkpointing once | exact |
| cadence: the same with it on | differs at the SSH solver tolerance — see §2c |

`tests/restart_roundtrip/job_restart_sweep` runs the four pairs. N = 1 is the trivial case that
passed before any of the fixes and is kept as a floor; N = 137 / 91 puts the interruption at a step
with no relation to the output or forcing cadence; N = 200 / 200 crosses several JRA55 3-hourly
brackets.

The remaining probe differences are all in arrays the step rebuilds before it reads them — at the
top of a resumed step the forcing has not been read yet and the scratch still holds the previous
step. They are printed under their own tags so they cannot be mistaken for state. Inside the step
the stage probes show `pgf_x`, `Kv`, `Av` differing at `1-eos` and `uv_rhs`, `ssh_rhs` at `1-eos`
and `3-mixing` — in every case before the stage that writes them — and nothing beyond that.

### The diagnosis tooling, and two traps in it

`FESOM_RESTART_PROBE=<step>` dumps a bitwise checksum (XOR of the IEEE patterns reduced with
`MPI_BXOR`, so it is independent of partitioning and reduction order) of every array the restart
carries plus the forcing and the scratch, at the entry to that step.
`FESOM_RESTART_PROBE_TOP=<step>` does the same at the top of the driver iteration, before the
forcing readers and the ice step, which separates "the restored state is wrong" from "the forcing
computed at that time is wrong". Each runs twice, once over owned entries and once over the full
local extent including halo, under separate tags. `tests/restart_roundtrip/job_restart_probe` runs
both legs and diffs all eight blocks.

`FESOM_IO_GATHER_DUPCHECK=1` counts, per gathered field, how many global element ids the gather
writes more than once and how many of those writes disagree. That is what identified §2c.

⚠️ **Trap 1, already paid for once.** `srun -l` prefixes every line with the rank, so
`grep '^\[probe\]'` matches nothing and produces two empty files — and an empty diff reads exactly
like "no differences". The job now uses `grep -o` and refuses to continue on an empty probe file.
No conclusion of the form "everything matches" is worth anything unless the probe files are
non-empty.

⚠️ **Trap 2.** The full-extent probe flags every stale halo, and most stale halos are read by
nothing. Reading it as a list of bugs sends you after arrays that are fine. Read the owned block
first; the full-extent block only says whether a difference is in the halo, and whether that matters
is a question about the reader, not about the checksum.

## 8. The Kokkos port (2026-08-18)

Done, and its byte gate is green. The code lives in `port_kokkos_int` on `m14-integrate`:
`src/fesom_io_restart.{cpp,h}`, `scatter_node`/`scatter_elem` + the duplicate check added next to the
gathers in `src/fesom_io.cpp`, and the same three driver hooks as the C.

**What is the same, on purpose.** The format, the field list, the guards and the canonicalisation are
the C reference's, so the two trees write interchangeable files — which is what lets
Kokkos-Serial-restarted be compared byte for byte against C-restarted rather than merely "closely".
`tests/restart_roundtrip/check_table_drift.py` compares the two tables directly and fails on any
divergence; it takes a second and should be run after touching either tree.

**What is new: the sync discipline.** Every array is a `fesom::Field`, so the host pointer the gather
reads is one of two copies. The table therefore carries the Field beside the host pointer:
`sync_host()` before the gather on write, `modify_host()` + `sync_device()` after the scatter on
read. Pairing them in the table rather than at the call site is deliberate — a field added without
its Field is then a visible hole rather than a CUDA-only wrong answer, and
`check_table_drift.py` also asserts that every entry's Field is the one belonging to its host
pointer. Nothing else can see that: on Serial the two views alias and every sync is a no-op.

**The case that makes it load-bearing** is the stiffness matrix. `fesom_update_stiff_mat_ale_kk` is a
device kernel ending in `modify_device()`, so under zstar the accumulated matrix lives on the device
from step 1 and the host mirror still holds the base matrix. Without the pull, a CUDA checkpoint
would carry the cold-start matrix and a resumed run would lose every accumulated increment — the
same failure §2b.2 describes, arriving by a different route. It also fixes the order of the restart
read in the driver: it must come after `fesom_ssh_preconditioner`, which pushes the base matrix to
the device.

**The C's cold-start bug was here too.** `fesom_step.cpp` passed `is_first_step = (step_n == 1)`;
`fesom_step_ctx` now carries `restarted` and the guard is `(step_n == 1 && !ctx->restarted)`.

### Gate: `jobs/job_restart_gate_kk`, 8 ranks, dt=1800, N = M = 10, zstar + cvmix_TKE + mEVP

Four legs, four comparisons, one job — the first two make a failure of the last two attributable.

| comparison | result |
|---|---|
| C round trip, A@20 vs B@20 | exact |
| Kokkos round trip, C@20 vs D@20 | exact |
| C == Kokkos-Serial, no restart, A@20 vs C@20 | exact |
| **C == Kokkos-Serial across a restart, B@20 vs D@20** | **exact** |

8 ranks and `OMPI_MCA_btl_vader_single_copy_mechanism=none` to match the established Kokkos-Serial
certification (`jobs/job_m14_gate_serial`), physics the paper's rather than that gate's default.

### CUDA: `jobs/job_restart_gate_cuda`, 4 GPUs, N = M = 10 — green

CUDA is not run-to-run reproducible on Levante A100, so "the restart is exact" is not a testable
statement there. What is testable: **the restarted legs are not distinguishable from the continuous
ones.** Four continuous legs and four restarted legs, six pairwise differences each — the same
statistic and the same count on both sides, so neither arm is favoured by taking a maximum.

| field | continuous vs continuous | restarted vs restarted |
|---|---|---|
| `temp` | 2.345e-04 | 1.815e-04 |
| `ssh`  | 5.479e-06 | 4.372e-06 |
| `u`    | 9.585e-05 | 8.248e-05 |

The restarted arm spreads *less* than the continuous arm in all three. Cross-backend: a
Serial-written restart is read by the CUDA build and a CUDA-written one by the Serial build, both
clean through the reader's own guards — format version, field list, node/element/level counts, mesh
signature. Those two legs are pass/fail on the guards and not on numbers, because a CUDA trajectory
and a Serial one are not the same trajectory and comparing them would be the sort of claim this
paper is about not making.

⚠️ **The first version of this job failed, and the failure was in the job.** Two errors, both worth
keeping because both are easy to repeat:

1. It estimated the noise floor from **one pair** of control legs. A single sample of a noisy
   quantity is not a floor — and the corrected six-pair floor for `temp` (2.35e-04) turns out to sit
   *above* the 2.01e-04 the one-pair version had flagged as a failure. The plan's "≥ 3 control legs"
   exists for exactly this.
2. Worse: the two arms did not share a checkpoint schedule. The round-trip legs checkpointed at N and
   at N+M, the controls only at N+M — so the round-trip legs took the §2c one-ulp canonicalisation
   kick at step N that the controls never took. An extra divergence source in one arm and not the
   other, compared as though the arms were alike.

**Task 2 is complete.** Both ports have restart I/O; the C round trip is exact at five (N, M) pairs;
Kokkos-Serial restarted is byte-identical to C restarted; the CUDA restart is indistinguishable from
continuous integration; and a restart file crosses backends in both directions.

## 9. Still open

1. **Calendar continuity** (§5.4): a run restarted at a year boundary must load the same JRA55
   records and the same monthly SSS as the uninterrupted run at that step. The 400-step sweep does
   not cross a month, so this is untested. Note that `io.prev_calendar` is deliberately not in the
   file — see the comment at the restart patch site in `fesom_main.c` for why that is safe and what
   would make it unsafe.
2. **Decide about §2c's alternative.** Keeping replicated elements consistent inside the model would
   remove the checkpoint-cadence dependence at the cost of a per-step communication. It is a
   decision about the model, not about restart I/O.
