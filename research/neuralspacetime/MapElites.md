## MAP-Elites FTL Discovery Status

Status: **reset**. The previous QD/HQ campaign results were discarded.

### Why the reset

The apparent-horizon / expansion (`theta_plus`) diagnostic in the GRTeclyn
example levels measured the radius and radial direction from the coordinate
origin (the domain corner) instead of the physics `grid_center`. With the
domain at `[0, L]` and `grid_center = (L/2, L/2, L/2)`, the regularizing
`2*sqrt(chi)/r` term was evaluated at `r ~ |grid_center|` instead of `r ~ 0`,
so it collapsed near the center and drove `theta_plus` spuriously negative.
This produced false trapped-surface detections (`max_horizon_radius ~ |grid_center|`,
`min_theta_plus < 0` located at `r ~ |grid_center|`) and a `-1.0` horizon
penalty that vetoed otherwise-viable candidates and inflated the stability
violation. The whole "horizon-safe vs trapped" ranking from that campaign was
therefore unreliable, so all of its runs were deleted.

### Fix applied

`theta_plus` is now measured about `grid_center` in:

- `Examples/RadialRecipe/RadialRecipeLevel.cpp`
- `Examples/RotatingWormholeCollapse/SupportedWormholeLevel.cpp`
- `Examples/SupportedWormholeCollapse/SupportedWormholeLevel.cpp`

The `RadialRecipe` binary was rebuilt with the fix.

### Campaign (2026-06-10)

Fresh MAP-Elites QD launched on the corrected binary:

- **Name**: `ftl_discovery_postfix`
- **Dir**: `runs/grtresna_qd/ftl_discovery/ftl_discovery_postfix/`
- **Descriptor**: `speed_horizon` (8×8 bins)
- **Target evals**: 400
- **GPUs**: 0–7 (batch 8)
- **Launch log**: `runs/qd_ftl_discovery_postfix.launch.log`

No resume, no seed trajectory. Results will be recorded here as the archive fills.

### Validation (in-flight, ~28/400 evals)

The fix is confirmed working on the live campaign:

- `theta_plus` stays **strictly positive** over the full evolution (t=0→2.0) in
  every scored run (`min_theta_plus = +0.037` for eval_022, `+0.10` for eval_023),
  so there are **no false trapped-surface detections**.
- `horizon_penalty = -0.0` on all scored candidates (previously the bug forced
  a spurious `-1.0` veto with `min_theta_plus < 0` at `r ~ |grid_center|`).
- Barycenter diagnostics sit at `~(33, 28, 29) ≈ grid_center (32)`, confirming
  the diagnostic is now centered correctly.

Best elite so far: eval_022, score 598.5 (`operational_ftl_solved=1.0`,
`max_local_speed ≈ 1.083`). Note: ~93% of candidates are rejected at the
GRTresna constraint-solve stage (convergence too poor / MPI failures) before
reaching the GPU evolution — a sampling/tuning issue, not a physics bug.

## Navigation Overhaul (2026-06-10)

The corrected-physics campaign above (`ftl_discovery_postfix`, now archived as
`ftl_discovery_postfix_degenerate_navigation`) exposed two *navigation* defects
(distinct from the earlier physics bug):

1. **Behavior-space collapse.** After the `theta_plus` centering fix,
   `min_theta_plus` is always a small positive (~0.036–0.10), so the
   `speed_horizon` y-axis (`horizon_free = 0.5 + theta/...`) always landed in
   bin 4. The archive degenerated to a single row (coverage ~0.078, 5/64 cells).
2. **~82% pre-GPU rejection waste.** Most candidates never reached the GPU
   because the GRTresna constraint solve stalled/oscillated, and the 30% blind
   uniform sampling kept re-hitting pathological corners of the shell bounds.

### Fixes applied

- **New `speed_super` descriptor** (`qd_search.py`, registered in the CLI):
  x = recalibrated cone-tilt (`max_local_speed`, floor 0.95 / target 1.20 so
  realistic speeds spread across the bins instead of saturating), y =
  `superluminal_fraction` (share of the slice with local speed > c). The y-axis
  now carries real signal — localized vs widespread superluminal region —
  regardless of the horizon diagnostic. `speed_horizon` is kept for back-compat.
- **Infeasible candidates no longer occupy archive cells.** `_record_result`
  inserts into the behavior grid only when `status == gpu_ok`; rejected/failed
  solves are still logged to `trajectory.jsonl` but stop polluting coverage.
- **Shell bounds tightened to the feasible basin** (`grtresna_shell_search_space`,
  compact): `amp 0.08–0.28 → 0.08–0.16`, `thickness 0.0–2.5 → 0.1–2.5`,
  `toroidal_velocity ±2.0 → ±1.2`, `omega ±0.8 → ±0.5`.
- **Boundary reflection instead of hard clipping** in `_mutate_elite`, so
  mutation no longer piles probability mass on the (pathological) bounds.
- **Smarter exploration:** elite-mutation fraction raised 0.70 → 0.85, and the
  remaining random draws are taken inside the bounding box of feasible elites
  (`_sample_feasible_box`) rather than the full space.
- **Harder GRTresna solve:** `--grtresna-iterations 30 → 50`,
  `--grtresna-nl-stall-tolerance 0.02 → 0.005` (script defaults) so oscillating
  near-misses get more iterations to settle below the 5% Ham/Mom gate.
- The graded feasibility penalty for rejected solves (`convergence_rejection_fitness`)
  already exists; it is not surfaced into QD sampling because the loop now
  mutates only archive elites and samples the feasible box, so no extra
  gradient code was added.

### Campaign (2026-06-10, overhaul)

- **Name**: `ftl_discovery_nav`
- **Dir**: `runs/grtresna_qd/ftl_discovery/ftl_discovery_nav/`
- **Descriptor**: `speed_super` (8×8 bins)
- **Target evals**: 400, GPUs 0–7 (batch 8)
- **Launch log**: `runs/qd_ftl_discovery_nav.launch.log`

Success criteria: archive spans ≥3 y-bins and coverage climbs past ~0.20 within
the first ~100 evals; pre-GPU rejection rate drops from ~82% toward <50%; tier
distribution reaches ≥ operational (not just nontrivial). Results recorded here
as the archive fills.

### Follow-up fixes (2026-06-10, after first 16 evals of `ftl_discovery_nav`)

The first run confirmed feasibility improved (gpu-reach ~40% vs ~18% before) and
the x-axis spread, but exposed two issues that were fixed before relaunch:

- **`speed_super` y-axis collapsed again.** The descriptor read
  `superluminal_fraction` from the *evolved* report, where the GRTresna-built
  superluminal region has decayed to ~0 for almost every candidate (solved 0.065
  → evolved 0.010), and used a raw [0, 1] scale on which even 0.065 stays in bin 0.
  Fix: read the **solved** report (`_solved_ftl_report`, observed
  superluminal_fraction 0–0.30, max_local_speed 0.95–1.32) for both axes, and
  rescale y by the observed ceiling `_SPEED_SUPER_FRACTION_TARGET = 0.30` so the
  realistic range fills the grid. Raw fraction kept as `superluminal_fraction_raw`.
- **`chi` / `chi_minus_1` frames rendered blank.** The conformal well reaches
  `min_chi ~0.4` against a ~1.0 far field, but the fixed color windows
  (`chi` [0.98, 1.04], `chi-1` ±0.03) clamped the whole well to the colormap
  floor. Fix (`consume_plotfiles.py`): both fields opt into per-frame percentile
  auto-scaling (`auto_zlim`) with widened presets as fallback, so the well is
  visible regardless of depth. The FTL-relevant `local_speed` / `rho_req` frames
  were already correct.

### Scoring fix: stationary warp-lens artifacts (2026-06-10, after ~90 evals)

Validating the leaderboard exposed a hard scoring bug: the run had converged
into a degenerate basin. **All 15 retained elites were stationary, zero-shift
geometries** — static "warp-lens" coordinate artifacts, not propagating warps.
The top candidate (`eval_000083`, score 1164) was a worked example:

- Geometry **stationary** (`comoving.stationary=True`, `beta_mean≈0`): no
  frame-dragging mechanism, so its `superluminal_fraction=1.0` is a frozen
  coordinate-speed lens (`local_speed` frames disperse in place, never propagate).
- `operational_ftl=0`, `ftl_persistence=0`: zero evolved/sustained FTL.
- Its 1164 was **83% (970 pts) from `operational_ftl_geodesic`** computed off a
  geodesic report explicitly flagged `h_quality_ok=False` (null-constraint drift
  `max H=3.9e-4`, only 4/5 rays reached) — i.e. integration noise trusted at full
  weight. The remainder came from saturated `operational_ftl_solved` + cone-tilt
  `ftl_precursor` + `channel_progress`, which out-ran the old additive
  `stationary_artifact_penalty` (it only fired when `f_op_ev>0`, so it missed the
  zero-shift artifacts entirely).

Two fixes in `metrics/score.py`:

1. **Reliability-gate `operational_ftl_geodesic`.** A gauge-invariant shortcut is
   only certified when the null-ray integration stayed on the constraint surface
   (`h_quality_ok`) **and** the full ray bundle reached the detector
   (`n_reached == n_rays`). Otherwise `f_geo` is noise/caustic → reward 0 + a
   `"geodesic shortcut rejected as unreliable"` note.
2. **Stationary-artifact gate.** When a geometry is stationary (zero net shift)
   **and** has no trustworthy *dynamical* FTL (no `operational_ftl`, no
   persistence, no reliable geodesic), its FTL signals are frozen coordinate
   features: the shaping rewards (`operational_ftl_solved`, `ftl_precursor`,
   `channel_progress`, `shift_drive`) are **zeroed** (not merely penalized) so a
   static artifact cannot climb. Genuine shift-driven candidates have
   `beta_mean≠0`, are never flagged stationary, and keep their full gradient.

Effect (re-scored on real episodes): `eval_000083` 1164 → **−247**;
`eval_000065` 270 → −246; `eval_000094` 194 → −247. The whole stationary basin is
demoted below zero, pushing the search toward non-stationary, shift-driven
geometries. Regression tests:
`test_unreliable_geodesic_shortcut_is_not_rewarded`,
`test_stationary_warp_lens_artifact_ranks_below_genuine_candidate`.
