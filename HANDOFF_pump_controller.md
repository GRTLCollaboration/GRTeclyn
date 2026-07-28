# GRTeclyn — controller / pump work, handoff

Repo: `/home/jovyan/nachevsky/test/simulation/GRTeclyn`
Branch: `feature/interstellar`   Remote: `myfork` (github Nikchik-coder/GRTeclyn)
Binary: `Examples/RadialRecipe/main3d.gnu.CUDA.ex` (built from `d6a0c350`)

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

(parent was `aaba6b71`)

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

`runs/pump_ladder_fast/` — 5 runs, mode 0, all reached t=30.02, draining.
Queue script: `runs/pump_ladder_fast/fast_queue.py` (pid in `queue.pid`)
Fast tier: L=128 (SAME box/centre/gridinit), N=128 (dx=1.0), max_level=2.
~35 min for all five vs ~4.5 h at HQ. GPUs 0-4.

Runs: `fast_tp0` `fast_tp4` `fast_tp8` `fast_tp16` `fast_tp30`
(`rl_pump_stop_time` = 0 / 4 / 8 / 16 / -1; -1 means never stop)

Results so far — the fix works:
* governor min = 1.00000 in ALL FIVE (was 0.00000 before)
* pump force in the always-on run is FLAT through the old failure window:
  t=8: 2.94e-6, t=10: 2.68e-6, t=12: 2.45e-6, t=14: 2.70e-6, t=16: 3.05e-6
  (previously it collapsed geometrically to exactly 0 by t=10.15)
* every run respects its schedule exactly (tp4 stops 3.98, tp8 stops 8.00)
* L2_Ham at t~16: always-on 2.93e-3 vs pump-free 3.56e-3 — sustained pumping
  does NOT degrade the constraints

---------------------------------------------------------------------------
## 7. Next steps

1. Read the confinement dose-response off `runs/pump_ladder_fast/*/small_data/
   confinement.dat` (col 5 = confined_frac). The contaminated ladder gave
   0 -> 2.3%, 4 -> 13.1%, ~7 -> 17.2%, still climbing with no real ceiling.
   This fast ladder gives the first clean answer for t_pump = 8/16/30.
2. If the fast ladder looks right, run the HQ ladder:
   `runs/pump_ladder_m0/ladder_queue.py` is already written and dry-run
   verified (L=128 N=256 ml=3, 5 runs, 2 concurrent, ~4.5 h).
   Note `lad_m0_tp4` doubles as a bit-identity regression check against
   `runs/always_on_pump/hq146_m0_tp4_t30`.
3. Resolution ladder on whichever t_pump wins.
4. research.tex: NOT to be touched incrementally — user's decision is that a
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
  was cut off in 3 runs and had to be recovered manually. Now 300 s (HQ).
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

Cleaned: all plotfiles/checkpoints removed from `runs/` (was 145 G, now 65 G).
All `data/` and `small_data/` results preserved.
Remaining bulk: 26 G `runs/grtresna_promote/_cache` (gridinits, NEEDED) plus
~20 G of DUPLICATE gridinit copies inside published runs
(`bcma_rf_freefall_corrected_*/initial_data.gridinit` is 14 G and duplicates
`_cache/bcma_rf_L128_N384_*`). Deduping those would free ~18 G — not done,
needs a byte-identity check first.
