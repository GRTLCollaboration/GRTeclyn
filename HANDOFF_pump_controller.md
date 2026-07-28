# GRTeclyn — controller / pump work, handoff

Repo: `/home/jovyan/nachevsky/test/simulation/GRTeclyn`
Branch: `feature/interstellar`   Remote: `myfork` (github Nikchik-coder/GRTeclyn)
Binary: `Examples/RadialRecipe/main3d.gnu.CUDA.ex` (built from `d6a0c350`, 2026-07-28 00:18 —
verified to postdate the governor-fix headers; later commits are docs-only, no rebuild needed)

---------------------------------------------------------------------------
## 1. TL;DR

The plan was "keep the PD pump on to t=30, pay for it with an explicit
energy-momentum reservoir so the Einstein constraints stay defensible".

Outcome: **the reservoir approach is dead, and it was never needed.**

* The reservoir (as specified: zero spatial stress) is mathematically dust.
  Dust with momentum and no pressure grows without bound. It diverged.
* Worse, in mode >= 1 the reservoir fed the pump's SAFETY GOVERNOR, which
  throttled the pump to exactly zero at t~10 in every run of the first
  campaign. That produced a fake "confinement saturates at t_pump=8" result.
* With the reservoir off (mode 0), the pump runs the full 30 time units at
  constant strength, and the constraints are FINE. The always-on run has
  *lower* constraint violation than the pump-free run.

Constraint defence is now: (a) pump-free run is worse than pumped, (b) a tight
Duhamel bound. No reservoir required.

---------------------------------------------------------------------------
## 2. Commits (all pushed to myfork/feature/interstellar)

* `fe5ef9f8` fix: correct controller reservoir conservation law and force diagnostics
* `d6a0c350` fix: never feed the controller ledger into the pump safety governor
* `2ddaa87a` docs: this handoff
* `a9376f77` docs: findings log — `research/neuralspacetime/Debug.md`
* `b97868de` docs: node-local plotfile scratch is now MANDATORY (see §10)

(parent was `aaba6b71`)

**RULE: never `scp` code to the cluster.** Edit locally -> commit -> push ->
`git pull` on the cluster -> run from the in-repo path. An scp'd file is
untracked, so the cluster silently diverges and the change is lost on the next
pull. This was violated once (the queue script) and has been corrected; the
duplicate at `runs/ladder_queue.py` was deleted.

---------------------------------------------------------------------------
## 3. What was broken and what was fixed

### 3.1 Pump-force diagnostic read ~1e-17 instead of ~1e-4  [FIXED, verified]
Root cause: the diagnostic used ABSOLUTE grid coordinates while the evolution
RHS uses CENTRE-RELATIVE coordinates (`Coordinates(iv, dx, center)` subtracts
the centre). So the pump envelope was evaluated in the wrong frame and
underflowed. This is why the paper's `E_pump ~ 1.6e-17` and the "soft
trajectory guide" narrative are artifacts — WITHDRAW THEM.
Fixed in `Examples/RadialRecipe/RadialRecipeLevel.cpp`.
Verified: pre-fix vs post-fix binaries give BIT-IDENTICAL physics columns and
differ only in `pump_force_L2` (1.56e-17 -> 1.16e-4).

### 3.2 PD law duplicated in 4 places  [FIXED]
New single source of truth: `Source/GRTeclynCore/RL/RLPumpForce.hpp`.
Call sites rewired: ComplexScalarField, BiComplexScalarField,
ComplexExoticScalarField, and the RadialRecipeLevel diagnostic.

### 3.3 Reservoir had 3 physics bugs  [FIXED, but see 3.4 — approach abandoned]
`Source/Matter/ControllerReservoirMatter.hpp`
1. Momentum source had the WRONG SIGN. This codebase defines
   `j_i = sum_A s_A (-Pi_A d_i phi_A)`, so the pump's contribution to matter
   momentum is `-f_i` and the reservoir must absorb `+f_i` — opposite to the
   energy case. Verified by direct flat-space computation of
   `d_t j_i + d_j S^j_i` for the scalar EMT.
2. Both sources carried a SPURIOUS LAPSE. `rhs[c_Pi] += S_A` is a
   coordinate-time rate, so `d_t rho|pump = f_perp` exactly with no alpha.
3. ALL TRANSPORT TERMS WERE MISSING (shift advection, `D_i j_c^i` flux,
   `j_ck d_i beta^k`, `alpha K j_ci`). It was a per-cell ODE.

NOTE: the plan's own densitised energy equation has a spurious `alpha K rho_c`
on the RHS — it fails the FLRW-dust check. The correct form for zero stress is
`d_t(sqrt(g) rho) + d_i[sqrt(g)(alpha j^i - beta^i rho)] = -sqrt(g) j^i d_i alpha`.

### 3.4 Reservoir is STRUCTURALLY UNSOUND  [NOT FIXABLE as specified]
Even fully corrected, `S_c^ij = 0` makes the reservoir pressureless dust.
The momentum equation has no `grad rho` term, so nothing pushes back; leftover
momentum is frozen in and keeps feeding the energy equation forever.
Evidence: pump force is EXACTLY 0 from t=6 (pump off) yet the ledger keeps
growing 1.9e-2 -> 9.8e-1 by t=30. Growth with zero input = self-divergence.
This would happen in exact arithmetic. It is the ANSATZ that is wrong.
In mode 2 (reservoir gravitates) it killed a run with a NaN in K at t=9.7.

### 3.5 Ledger fed the safety governor  [FIXED — this is the big one]
`compute_radial_recipe_constraint_norms` publishes `cached_L2_Ham`, which
drives the pump governor (tanh, centre 0.035, width 0.003). It called
`fill_active_constraints`, which in mode >= 1 INCLUDES the reservoir. So the
controller's own bookkeeping was wired into its safety interlock.
Consequence: reservoir diverged -> reported L2_Ham crossed 0.035 at t~7.5 ->
governor cut pump to zero by t~10 in EVERY mode-1 run.
Governor minima: mode-0 baseline 1.000, pump-free 1.000, every mode-1 run
with an active pump 0.000.
Fix: `fill_active_constraints` gained `reservoir_mode_override`; the governor
path passes 0 so it always sees the true physical constraint violation.
Files: `RadialRecipeMatterDispatch.hpp`, `RadialRecipeConstraintNorms.hpp`.

### 3.6 Diagnostics added  [DONE]
`constraint_norms.dat` is APPEND-ONLY. Column layout (1-indexed):
```
1 time  2 L2_Ham  3 L2_Mom  4 min_rho_req  5 max_rho_req  6 integral_neg_rho
7 L2_Ham_rel  8 L2_Mom_rel  9 pump_force_L2  10 governor  11 pump_fi_L2
```
`pump_fi_L2` (col 11) is new — real measurement of the momentum force density
via a 4th-order stencil (state_new is FillPatch'd with 2 ghosts, so the
earlier "needs a ghost-filled D1 pass" comment was wrong).
`pump_work` in `collapse_diagnostics.dat` (col 15) is now the CORRECTED
`integral f_perp sqrt(gamma) dV`. Same column name, NEW meaning — do not
compare against pre-`fe5ef9f8` runs.

### 3.7 Geodesic emit gate  [DONE]
`GEODESIC_EMIT_MIN_TIME` now takes precedence over `RL_PUMP_STOP_TIME` in
`evolving_geodesic.py` / `ftl.py`; default preserves old behaviour. Tests updated.

---------------------------------------------------------------------------
## 4. Verification performed

* mode 0 is BIT-IDENTICAL before/after all changes.
* Null test (pump off + mode 2 vs pump off + mode 0): BIT-IDENTICAL. Passes.
* Mode 1 does not perturb the evolution (min/max rho_req agree to 10 sig figs).
* CAUTION LEARNED: I first validated at t=2 and it looked good. The ledger
  only turns bad after t~2.5. DO NOT validate this system on short runs.

---------------------------------------------------------------------------
## 5. Data status

### CONTAMINATED — do not use constraint columns (2,3,7,8)
`runs/always_on_pump/hq146_m1_tp4_t30`, `..._m1_tp8...`, `..._m1_tp16...`,
`..._m1_tp30...`  (mode 1; ledger diverged AND governor closed at t~7-8;
effective pump duration was ~7 in all of them regardless of setting, which is
why tp8/tp16/tp30 are bit-identical)
`runs/always_on_pump/hq146_m2_tp30_t30` — CRASHED, NaN in K at t=9.69.

Their *physics* outputs (confinement.dat, collapse_diagnostics.dat,
ftl_timeseries.dat) are still valid up to the point the governor closed.

### VALID
`runs/always_on_pump/hq146_m0_tp4_t30`  — mode 0, governor 1.0 throughout.
`runs/always_on_pump/hq146_m1_tp0_t30`  — pump never on, reservoir identically 0.

Key valid numbers:
* pump-free final L2_Ham 3.79e-3  vs  pumped-to-t=4 final L2_Ham 2.92e-3
  -> THE PUMP IS NOT DRIVING CONSTRAINT GROWTH (direct evidence, not a bound)
* Duhamel bound on m0_tp4: measured peak L2_Ham 3.14e-3 vs bound 3.59e-3
  (1.15x, tight). Momentum 7.9e-4 vs 1.75e-3 (2.2x). Governor never engaged.
* f_i is ~3.5x LARGER than f_perp — the momentum forcing dominates.

---------------------------------------------------------------------------
## 6. What is running RIGHT NOW

The **fast ladder is FINISHED** (5 runs, mode 0, all reached t=30.02). Its
directory `runs/pump_ladder_fast/` has had its plotfiles pruned; `data/` and
`small_data/` are intact. Headline fast-tier results are in
`research/neuralspacetime/Debug.md`.

The **HQ ladder is RUNNING**: `runs/pump_ladder_m0/`, 6 runs, ALL SIX
CONCURRENT on GPUs 0-5. This is only possible because plotfiles now go to
node-local NVMe — see §10.

```
launcher : grteclyn-wrapper/scripts/campaigns/rl/pump_ladder_queue.py
log      : runs/ladder_queue.out   (also runs/pump_ladder_m0/queue.log)
launched : 01:40   ETA ~02:50   (measured 35 rows/min = 0.35 t/min to t=30)
grid     : L=128 N=256 max_level=3, stop_time=30, plot_interval=144
```

| run | `rl_pump_stop_time` | purpose |
|---|---|---|
| `lad_m0_tp0`  | 0.0  | pump-free control rung |
| `lad_m0_tp4`  | 4.0  | bit-identity regression vs `runs/always_on_pump/hq146_m0_tp4_t30` |
| `lad_m0_tp8`  | 8.0  | |
| `lad_m0_tp16` | 16.0 | last rung healthy at fast tier |
| `lad_m0_tp24` | 24.0 | brackets the 16->30 turnover |
| `lad_m0_tp30` | -1.0 | always on; headline + lapse-collapse test |

Health at t~5-6 (20 min in), all confirmed good:
* **governor = 1.000 in all six** — the `d6a0c350` fix confirmed in production
  (previously it clamped to 0 at t~10).
* **determinism check PASSES exactly.** tp8/tp16/tp24/tp30 differ only in stop
  time, all > 5, so at t=5.01 they must be identical — and they are to every
  printed digit: `3.4418205697e-03 / 2.6949606572e-04 / 3.0898006783e-06`.
* consumers clean, no NaN, scratch flat at 3 plotfiles/run.

`lad_m0_tp30` may blow up near t~29.8 as it did at the fast tier. That is a
GAUGE pathology (min_lapse -> 1e-10 while min_chi stays 0.337), NOT collapse,
and is itself a result. Do not treat it as a failed run.

**ONE THING TO WATCH — may change the headline.** At t=5.01 the pumped runs
sit at L2_Ham 3.44e-3 against 1.03e-3 pump-free (3.3x), and tp4 has relaxed
only to 3.01e-3 a full time unit after switch-off. The fast-tier claim was
"statistically indistinguishable", but that was a WHOLE-RUN MEAN (2.84e-3
unpumped) — at t=5 the unpumped run is still far below its own average, which
inflates the ratio. Unresolved until the curves reach t=15-30. **If the gap
holds at late times, the "pump does not affect constraint growth" claim must be
restated at 256^3.** Check this first when the runs finish.

---------------------------------------------------------------------------
## 7. Next steps

1. **When the HQ ladder finishes**, in order:
   a. Resolve the open question above — plot L2_Ham(t) for all six rungs and
      compare whole-run mean AND peak, not instantaneous values.
   b. Confirm `lad_m0_tp4` is bit-identical to
      `runs/always_on_pump/hq146_m0_tp4_t30` (the governor fix must be a no-op
      in mode 0).
   c. Re-run `research/neuralspacetime/analysis/pump_constraint_budget.py` for
      the Duhamel bound at HQ resolution. It reads `pump_fi_L2` (col 11) for
      the momentum bound.
   d. Confinement dose-response from `*/small_data/confinement.dat`
      (col 5 = confined_frac). Fast tier gave 1.9% -> 20.3% (10.7x), monotonic,
      no saturation — but confinement DECAYS everywhere; the pump only slows
      the leak. Check whether that survives at 256^3.
   e. `f_geo`: fast tier showed NO dose-response (26.9 / 27.5 / 31.0 / 23.4 /
      23.4) — the pump-free run beat tp16 and tp30. This is the weakest part of
      the story. Verify at HQ before claiming anything about FTL geometry.
2. Resolution ladder on whichever t_pump wins.
3. research.tex: NOT to be touched incrementally — user's decision is that a
   validated fix means a full re-run and a total rewrite.

---------------------------------------------------------------------------
## 8. Traps / gotchas (learned the hard way)

* `amr.plot_file` and `amr.check_file` in the RM baseline params are ABSOLUTE
  paths into the baseline dir. If you clone params without redirecting them,
  your run writes plotfiles into the baseline AND the `--delete` consumer
  prunes them there. Both queue scripts redirect them.
* Plotfiles are deleted after extraction, so ANY extraction not enabled during
  the run is unrecoverable. `--areal-radius` and `--shell-fields` were NOT
  enabled in these runs (matching the previous campaign's config), so no
  `areal_radius.dat` / `shell_profiles.dat`. Horizon info is still available:
  `collapse_diagnostics.dat` is written by the simulation itself.
* The consumer's drain window must be generous. At 180 s the final plotfile
  was cut off in 3 runs and had to be recovered manually. **Now 600 s.**
* Scratch is TRANSIENT. A plotfile deleted from `/tmp` before extraction is
  gone for good (on NFS it would have survived a stalled consumer). Purge a
  run's scratch only when every resident plotfile appears in
  `small_data/consume_state.json`; otherwise keep it and log `[KEEP-SCRATCH]`.
* The pump does NOT track the matter. `trajectory_mode = 1` overwrites the
  lump centres every coarse step from a prescribed parametric path
  (`TrajectoryEvaluator`), and sets amplitude from `well_depth` (constant).
  See `Main_RadialRecipe.cpp:92-117`.
* Build (single GPU, per wrapper README):
  `cd Examples/RadialRecipe && PATH="/usr/local/cuda/bin:$PATH" \
   NO_MPI_CHECKING=TRUE make USE_MPI=FALSE USE_CUDA=TRUE -j 32`
  Do NOT put the grtresna conda env on PATH (its g++ 15 breaks nvcc).
* Gridinit does NOT need rebuilding for lower resolution — the reader
  trilinear-interpolates. Only a change of BOX SIZE (L) would need a re-solve.

---------------------------------------------------------------------------
## 9. Disk

Cleaned: all plotfiles/checkpoints removed from `runs/`, plus 5 contaminated /
crashed HQ runs from the first campaign (was 145 G, now **62 G**).
All `data/` and `small_data/` results preserved.
Kept deliberately: `hq146_m0_tp4_t30` (bit-identity reference for
`lad_m0_tp4`) and `hq146_m1_tp0_t30`.
Remaining bulk: 26 G `runs/grtresna_promote/_cache` (gridinits, NEEDED) plus
~20 G of DUPLICATE gridinit copies inside published runs
(`bcma_rf_freefall_corrected_*/initial_data.gridinit` is 14 G and duplicates
`_cache/bcma_rf_L128_N384_*`). Deduping those would free ~18 G — not done,
needs a byte-identity check first.

---------------------------------------------------------------------------
## 10. Plotfile scratch is NODE-LOCAL — mandatory for every run

This is the change that unblocked 6-way concurrency. Full documentation:
`grteclyn-wrapper/README.md` and `grteclyn-wrapper/scripts/campaigns/README.md`.
Reference implementation:
`grteclyn-wrapper/scripts/campaigns/rl/pump_ladder_queue.py`.

```
GPU sim ──plotfiles──▶ /tmp/grteclyn_scratch/<run>/RadialRecipePlt*  (local NVMe)
   │                        │ consumer extracts (~16 s), then deletes (keep-last 3)
   │                        ▼
   └──.dat, KB/s──▶ <output_path>/data/  +  small_data/          (NFS, ~500 KB)
```

`amr.plot_file` / `amr.check_file` are INDEPENDENT of `output_path`, so:

```
output_path    = "<NFS>/runs/pump_ladder_m0/<run>"
amr.plot_file  = "/tmp/grteclyn_scratch/<run>/RadialRecipePlt"
amr.check_file = "/tmp/grteclyn_scratch/<run>/RadialRecipeChk"
```

Consumer: `--data /tmp/grteclyn_scratch/<run>` `--out <NFS>/.../small_data`.

**Why.** A 256^3 ml=3 plotfile is ~3.2 GB, written AND read back every ~288 s
per run. On NFS that capped concurrency at 2 — consumers blocked in NFS I/O
(`D` state), backlogs grew, plotfiles accumulated. Measured after the move
(2026-07-28, 6 concurrent HQ runs): local scratch 8.8 GB/run, NFS run dirs
~500 KB each, extraction 15.7 s against a 288 s cadence (18x headroom), zero
plotfiles leaked to NFS, 53 GB of 1.1 TB overlay used.

**Keep every write inside your own directories.** On this cluster
`~/.local/bin` and `~/.local/lib` are owned by `nobody:nogroup` and are NOT
writable — admin policy is "write only to your own". An agent previously
tampered with them and admins went looking for the culprit. So:
* call `grteclyn-wrapper/.venv/bin/python -m ...` DIRECTLY, never `uv run`
  (which writes `~/.cache/uv`);
* pin `XDG_CACHE_HOME`, `UV_CACHE_DIR`, `MPLCONFIGDIR`, `TMPDIR`,
  `PYTHONPYCACHEPREFIX` into `/tmp/grteclyn_scratch/_cache`.

The complete set of write targets is `/tmp/grteclyn_scratch/` (transient) and
the NFS run directory (`data/`, `small_data/`, `run.log`, `params.txt`).
Nothing else. Do not write anywhere under `~/.local`.

**Not yet enforced in code.** `grteclyn-wrapper/src/grteclyn_wrapper/core/
params.py` still defaults `amr.plot_file` / `amr.check_file` to the episode
directory, so on an NFS-backed `RUNS_DIR` the DEFAULT still lands on NFS.
Until that default is changed this is a launcher-level convention that every
campaign script must apply explicitly. Fixing the default would be a good
piece of work.
