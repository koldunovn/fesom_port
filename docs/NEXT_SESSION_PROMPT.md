# Next session: debug Phase D EVP multi-rank crash

## TL;DR

Phase A-C of the FESOM2 sea-ice port are committed and working. Phase D
(standard EVP) is on branch `debug-evp` but crashes at multi-rank. 1-rank
passes 200 steps; 16-rank crashes at step 2 (BC on) or step 70 (BC off).
A code-level halo-loop audit against Fortran ice_EVP.F90 is clean. Need to
add runtime NaN detection to find the first bad node and trace from there.

## Read these memory files first (in order)

These give you full background before touching code:

```
~/.claude/projects/-home-a-a270088-port2/memory/MEMORY.md
~/.claude/projects/-home-a-a270088-port2/memory/project_sea_ice_port_state.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_write_loops_halo.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_unique_outdir_per_job.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_levante_build.md
~/.claude/projects/-home-a-a270088-port2/memory/feedback_port_fidelity.md
```

The first two are the most important — `project_sea_ice_port_state.md` has
the full debug status from the last session.

## Repo state

Current working directory: `/home/a/a270088/port2/fesom2_port`

```
main          e60a418  Revert "Phase D EVP" — current Phase C baseline (works at 16 ranks)
              163090d  Phase D EVP commit (cherry-pick exists on debug-evp)
              edc9e65  Phases A-C
              3da152b  Phase 4 wrap-up (ocean only, pre-sea-ice)

debug-evp     0490c13  Phase D EVP cherry-pick (BROKEN at multi-rank)
              + stash{0}  BC-disabled debug edit (for bisecting)
```

Always do clean builds:
```
bash -l configure.sh --clean
```

Always use a unique OUT_DIR per SLURM job (concurrent jobs writing to the
same dir destroyed several hours of debug last session — see
`feedback_unique_outdir_per_job.md`).

## The bug

**Symptom:** Phase D EVP code (commit 163090d, on branch `debug-evp`):
- 1-rank @ 200 steps: PASS, T.min=−2.06°C, identical to Fortran-style result
- 16-rank @ 200 steps: CRASH at step 2 (BC on) or step 70 (BC off), CG NaN

**Already-tried bisects (all multi-rank, branch `debug-evp`):**

| Test | Result | Conclusion |
|---|---|---|
| Full Phase D | crash step 2 | Reproducible |
| Subcycle disabled (skip stress_tensor + stress2rhs + velocity update) | PASS 200 | Bug is in subcycle body |
| evp_rheol_steps=1 (1 iter instead of 120) | crash step 2 | Bug appears in single iter, not from accumulation over 120 |
| Subcycle on, BC loop disabled | crash step 70 (slow buildup) | BC isn't the root, but accelerates blowup |
| 1-rank, subcycle on, BC off | PASS 200 | Confirms multi-rank-only bug |

**Code-level audit done — all CLEAN against Fortran:**
- Every C `for ... < myDim_*2D` loop bound matches Fortran loop bound in `ice_EVP.F90`
- Only halo exchange in Fortran EVPdynamics is `exchange_nod(U_ice, V_ice)` at end of subcycle (line 827) — present in C
- Thermo only halo-exchanges `ustar` (line 238 in F) — present in C
- ocean2ice halo-exchanges `u_w`, `v_w` (line 244) — present in C
- Verified empirically that the partition is overlap-based (rank 0 and rank 1 share 130 elements in their myDim_elem2D), so Fortran's "Pattern A" scatter from `myDim_elem2D` only is correct
- mesh `metric_factor`, `elem_cos`, `elem_area` halo-exchanged (mirrors Fortran lines 2521-2522) — present in C
- `gradient_sca` only allocated for myDim_elem2D in Fortran too — matches port

So static reading isn't finding the bug. Need runtime instrumentation.

## What to try next

**Step 1: NaN detection at the source**

Add diagnostic prints to `src/fesom_ice_evp.c` in `fesom_ice_evp_dynamics` —
after each subcycle iteration, scan `u_ice` and `v_ice` for NaN/Inf and
`fprintf(stderr, ...)` the first bad local node index + global node id +
which rank.

Suggested code (insert just after `fesom_exchange_nod2D(v_ice, partit);` at
end of subcycle):

```c
for (int n = 0; n < N; ++n) {
    if (!isfinite(u_ice[n]) || !isfinite(v_ice[n])) {
        fprintf(stderr, "[rank %d] NaN at iter %d local n=%d gid=%d "
                        "(myDim=%d, eDim=%d) u=%g v=%g a_ice=%g\n",
                partit->mype, sub, n, partit->myList_nod2D[n],
                mesh->myDim_nod2D, mesh->eDim_nod2D,
                (double)u_ice[n], (double)v_ice[n], (double)a_ice[n]);
        break;
    }
}
```

Then submit a 16-rank job (NSTEPS=2, FESOM_PRINT_EVERY=1, fresh OUT_DIR)
and look at the first NaN report.

**Step 2: Trace inputs at the bad node**

Once you know which node first goes NaN, find what input value at that node
is wrong. Likely candidates:
- `inv_mass[n]` — if NaN, drag → NaN, velocity → NaN
- `u_w[n]` or `v_w[n]` — if NaN at halo (ocean2ice halo coverage)
- `inv_areamass[n]` — if NaN
- `stress_atmice_x[n]` (should be 0; if not, something corrupted)

**Step 3: Backwards trace**

If e.g. `inv_mass` is NaN at a halo node where Step 2 of EVPdynamics didn't
write it (myDim only), then maybe something earlier wrote NaN to it. Or
maybe `a_ice` at that node is NaN.

The hint from the user during the last session was Russian:

> "скорее всего где то забыты коммуникации, или же цикл не включает гало
> (в некоторых местах мы делаем вычисления и по гало, чтобы избежать
> коммуникации"

Translation: "most likely there are forgotten communications somewhere, or
a loop that doesn't include halo (in some places we compute over halo too
to avoid communication)". Keep this lens — even though static audit is
clean, NaN detection will find the actual missing piece.

## Files to focus on

```
/home/a/a270088/port2/fesom2_port/src/fesom_ice_evp.c    # the suspect
/home/a/a270088/port2/fesom2_port/src/fesom_ice_evp.h
/home/a/a270088/port2/fesom2_port/src/fesom_ice.c         # dispatcher + step driver
/home/a/a270088/port2/fesom2/src/ice_EVP.F90              # Fortran ground truth (read in full)
```

Reference (do NOT modify Fortran source — it's the ground truth):
```
/home/a/a270088/port2/fesom2/src/ice_EVP.F90
/home/a/a270088/port2/fesom2/src/MOD_ICE.F90
/home/a/a270088/port2/fesom2/src/oce_mesh.F90        # mesh setup, halo exchanges
```

## Job script template

There's a `job_core2_16_fix` template at the repo root. Always:
1. `cp job_core2_16_fix /tmp/job_my_test`
2. Edit OUT_DIR to a UNIQUE name like `/work/ab0995/a270088/port/core2_16r_evp_nantrace_$(date +%s)`
3. `mkdir -p $OUT_DIR`
4. `sbatch /tmp/job_my_test`

For 1-rank tests:
```
sed -i 's|--ntasks=16|--ntasks=1|; s|--ntasks-per-node=16|--ntasks-per-node=1|' /tmp/job_my_test
```

Note 1-rank takes ~14 min wallclock for CORE2 200 steps; 16-rank is ~3 min.

## Workflow rules from last session

1. **ALWAYS clean build** (`bash -l configure.sh --clean`). Stale objects
   from previous source state will silently link with new source — symptom
   looks like flaky behavior. Burned hours on this.
2. **ALWAYS unique OUT_DIR per SLURM job.** Concurrent jobs sharing the
   same `/work/.../` directory clobber each other's `run.log` and
   `snap_*.nc`, producing logs that look like crashes but aren't. Burned
   hours on this too.
3. **2× repeat tests** before declaring a phase "validated". One PASS is
   not enough — the EVP commit was based on a single PASS that turned out
   not to reproduce.
4. **Commit at every working state** with a clear message about what's
   verified. The Phase C commit (`edc9e65`) saved us today; without it the
   session would have been a total loss.
5. **Follow Fortran literally.** When tempted to invent a structure
   Fortran doesn't have (like a `edge_is_boundary` per-edge array), stop
   and find the Fortran convention instead (`myList_edge2D > edge2D_in`
   was the right answer).

## What's already saved

- This prompt: `/home/a/a270088/port2/fesom2_port/docs/NEXT_SESSION_PROMPT.md`
- Sea-ice port plan: `/home/a/a270088/port2/fesom2_port/docs/plans/20260425-sea-ice-port.md`
- Phase D EVP code: branch `debug-evp` in this repo
- D5 attempt patch (`oce_fluxes_mom`): saved at `/tmp/sea_ice_modified.patch` and `/tmp/sea_ice_work/` (may be cleared by /tmp policy — back up if needed)

## When EVP is fixed

1. Verify with **3 fresh runs** at 16-rank (different OUT_DIRs) all passing
2. Verify at 1-rank still passing
3. Verify at 32-rank passing (different partition exposes different halo
   patterns)
4. Commit to `debug-evp`
5. `git checkout main && git merge debug-evp`
6. Move to D5 (oce_fluxes_mom) — saved attempt is mostly correct, just
   needs to ride on a working EVP

## Then Phase E

FCT advection of ice tracers. Plan calls it the highest-risk phase
(`init_tracers_AB` and `fct_ttf_max/min` halo bugs in ocean FCT are
expected to recur). Do the proactive halo audit per
`feedback_write_loops_halo.md` as you write it.
