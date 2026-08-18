# Restart I/O for the C port (and then for Kokkos)

> **Revision 2, 2026-08-18, after the first round-trip gates.** §2 below was written from a static
> reading of the code and was incomplete: the gate found four further pieces of persistent state,
> all of them consumed before they are written and none of them in the Fortran restart. They are
> listed in §2b. The gate is not yet green — see §7 for exactly where it stands.

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

**Round trip, maximum absolute difference 0.** Run N + M steps continuously and write a restart at
the end. Separately run N steps, write a restart, start a second executable from it, run M more,
write a restart. The two final restart files must agree bit for bit, field by field. The restart
file is the right object to compare because it *is* the model's persistent state — a snapshot
comparison would only cover the subset that gets written to the output streams.

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

**Control first.** `tests/restart_roundtrip/job_restart_control` runs the identical configuration
twice, no restart anywhere, and the two restart files agree bit for bit. The model is reproducible
run to run at 128 ranks, so any round-trip difference is real and attributable to the restart.

| case | result |
|---|---|
| unit tests, `test_restart_format` | 79 checks, green (round trip, component order, optional TKE, four guard cases) |
| write twice in one run, compare against a second run's file at the same step | exact |
| restart after 1 step, then 1 step | exact |
| restart after 1 step, then 20 steps | exact |
| restart after 20 steps, then 1 step | **differs at roundoff**: max 3.6e-12 on `Z_3d_n`, 1.2e-14 on `d_eta`, 2.3e-14 on `hbar` |
| restart after 20 steps, then 20 steps | differs at ~1e-3, i.e. the roundoff above amplified |

So something small is still missing, and it is present at step 20 but not at step 1 — which is what
a state variable that is still zero after one step looks like.

**Diagnosis tooling in place.** `FESOM_RESTART_PROBE=<step>` dumps a bitwise checksum (XOR of the
IEEE patterns, reduced with `MPI_BXOR`, so it is independent of partitioning and reduction order) of
every array the restart carries plus the forcing and the scratch, over the FULL local extent —
owned *and* halo — at the top of that step. `FESOM_RESTART_PROBE_TOP=<step>` does the same at the
top of the driver iteration, before the forcing readers and the ice step, which separates "the
restored state is wrong" from "the forcing computed at that time is wrong".
`tests/restart_roundtrip/job_restart_probe` runs both legs and diffs them.

⚠️ **A trap in that tooling, already paid for once.** `srun -l` prefixes every line with the rank,
so `grep '^\[probe\]'` matches nothing and produces two empty files — and an empty diff reads
exactly like "no differences". The job now uses `grep -o` and refuses to continue on an empty probe
file. Any conclusion of the form "everything matches" is worthless unless the probe files are
non-empty.

## 8. Next

1. Run `job_restart_probe` with `FESOM_RESTART_PROBE_TOP` and read the *first* differing array at
   the top of the first resumed iteration. That single answer decides everything below.
2. Candidates not yet excluded, in order of suspicion:
   - `io.prev_calendar`. On a restart it is set equal to `io.calendar`; in the uninterrupted run it
     is one step behind. It feeds `fesom_sss_runoff_step_cal` and the monthly-crossing tests. It
     cannot explain a difference at 10:00 on 1 January — no boundary is crossed either way — but it
     is wrong and should be saved in the file regardless.
   - The JRA55 reader's bracket. `fesom_jra55_step` refreshes coefficients only when `rdate` leaves
     `[t_indx, t_indx_p1]`, so a resumed run picks its bracket from scratch. The coefficients depend
     on the two records and not on `rdate`, so this *should* be identical — verify rather than assume.
   - Anything reached through `mesh` that the step writes and that is not in the file.
3. Then the Kokkos port, the cross-backend read, and the calendar-continuity check (§5.1, 5.2, 5.4).
