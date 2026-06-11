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

## Matter model — reference & future directions (2026-06-10)

This is a map of *what the lumps actually are* and how to extend them, so we
don't have to re-trace the code each time. The `shell`/MAP-Elites campaign
evolves **N independent massive real scalar fields** ("lumps"), not a single
field and not massless matter.

### What the campaign actually runs

Each eval's `params.txt` sets:

```
recipe_matter_model     = grtresna_independent_scalars
recipe_num_scalar_fields = 5
recipe_scalar_mass       = 0.1
recipe_exotic_matter     = 0     # per-lump exotic flags are used instead
```

So GRTeclyn dispatches to `GRTresnaIndependentScalars` (5 fields, mass m=0.1),
**not** the `ScalarField<DefaultPotential>` / `ExoticScalarField` paths.

### Where the matter lives (file map)

**GRTeclyn (evolution side):**
- `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp` — runtime model
  selection. `uses_independent_scalars()` picks `GRTresnaIndependentScalars`
  when `recipe_matter_model == "grtresna_independent_scalars"`; otherwise
  `ExoticScalarField<DefaultPotential>` (if `recipe_exotic_matter`) or
  `ScalarField<DefaultPotential>`.
- `Source/Matter/GRTresnaIndependentScalars.{hpp,impl.hpp}` — the per-lump
  fields. `compute_emtensor` sums `sign[k]·T_ab(φ_k)` over lumps;
  `add_matter_rhs` is the curved-space Klein–Gordon RHS per lump.
- `Source/Matter/GRTresnaScalarPotential.hpp` — the **only** non-trivial
  potential currently wired in: `V = ½ m² (Σφ_k)²` (massive, shared sum).
- `Source/Matter/DefaultPotential.hpp` — `V=0` (massless); used by the legacy
  `ScalarField`/`ExoticScalarField` paths. `recipe_scalar_mass` is **ignored**
  on those paths — only `GRTresnaIndependentScalars` honors the mass.
- `Source/Matter/{ScalarField,ExoticScalarField}.{hpp,impl.hpp}` — canonical
  (ρ≥0) and phantom (ρ≤0, T scaled by `-recipe_support_strength`) single-field
  models. Templated on the potential class.
- `Examples/RadialRecipe/StateVariables.hpp` — state layout: `c_phi, c_Pi`
  plus `phi_lump0..4 / Pi_lump0..4` (`2 * GRTRESNA_MAX_INDEPENDENT_SCALARS`).
- `Examples/ScalarField/Potential.hpp` — reference `½ m² φ²` potential class
  (template signature to copy when adding new potentials).

**GRTresna (initial-data side, sibling repo `../GRTresna`):**
- `Examples/ScalarFieldBH/MatterParams.hpp` — the **lump definition**. `lump_t`
  = {amp, width, center, velocity (boost→linear momentum), omega
  (rotation→L_z), mode (0=axisym, 1=dipole, 2=quadrupole), exotic flag}.
  `lump_phi` = `amp·angular_factor(mode)·exp(-r²/2w²)`; `lump_pi` =
  `-(v·∇φ) - omega·∂_azimuth φ` (so each cloud carries net P_i ~ v and L_z ~ omega).
  Exotic lumps are damped by `EXOTIC_AMP_SCALE = 0.25` to stay in the
  Lichnerowicz/York convergent regime.
- `Examples/ScalarFieldBH/MyMatterFunctions.cpp` — paints initial `(phi, Pi)` =
  background + `Σ_k lump_phi/lump_pi`. `my_potential_function` uses
  `½ (scalar_mass·φ)²` (must match GRTeclyn's potential).

### How lumps interact / change shape (answers we keep re-deriving)

- **Shape:** live wave fields under Klein–Gordon — they oscillate, disperse,
  deform. Not rigid. Initial shape is the analytic Gaussian cloud above.
- **Interaction:** only two indirect channels — **(a) gravity** (all lumps
  source one shared metric) and **(b) a shared mass potential** `½m²(Σφ_k)²`
  whose gradient `m²·Σφ_j` is applied to *every* lump's `Pi` RHS
  (`GRTresnaIndependentScalars.impl.hpp`), so overlapping lumps cross-couple.
  There is **no** direct gradient-gradient or self-interaction term.
- **Why they "fly away":** initial boosts are O(1) (`lumpK_velocity ~ |1|`),
  amplitudes are tiny (`amp ~ 0.09` → negligible self-gravity), and the mass is
  light (m=0.1 → Compton wavelength 1/m=10 ≈ box). Nothing binds them, so they
  free-stream out. This is physically expected, not a bug — there is no soliton
  / bound-state mechanism in the current matter sector.

### Future directions (ranked: leverage vs. effort)

1. **Search the mass + cap boosts (config-only, do first).** *Done
   (2026-06-11).* `m=0.1` and O(1) velocities were why lumps dispersed/flew
   away. `grtresna_scalar_mass` is now a searched QD dimension (range
   `0.3–1.5`, default `0.6`) wired into `grtresna_shell_search_space` and
   applied to `cfg.scalar_mass` in `_apply_grtresna_overrides`
   (`search/optimize.py`); it already propagated to both the GRTresna solve
   (`solver.py`) and the GRTeclyn evolution (`matter_wiring.py`), so the
   consistency rule holds with no new physics code. The toroidal current (warp
   motor) keeps its full range, while the net-outflow velocities are capped
   (`poloidal ±1.5→±0.8`, `radial ±0.8→±0.3`) so matter stays bound. Heavier
   mass binds within the shell width (Compton wavelength `1/m`), letting the
   search settle on persistent, near-static matter.
2. **Add `λφ⁴` (or φ⁶) self-interaction ⇒ oscillons / Q-balls** (genuinely bound,
   long-lived). Extend `GRTresnaScalarPotential` (GRTeclyn) **and**
   `my_potential_function` (GRTresna `MyMatterFunctions.cpp`) identically, plus
   one `lambda` param threaded through both param files.
3. **Complex scalar field ⇒ boson stars / Q-balls** with conserved U(1) charge —
   the textbook *stationary, non-dispersing* lump. Larger lift: two real
   components per field, and constraint-solver support on the GRTresna side.
4. **Per-lump independent mass** (currently one shared `scalar_mass`) to mix a
   heavy "anchor" lump with light mobile ones.

### Hard consistency rule (do not violate)

`T_ab` used in the GRTresna **constraint solve** must equal `T_ab` used in the
GRTeclyn **evolution**, or the run starts off-constraint and any "FTL" is a
constraint-relaxation transient (root cause documented in
`Examples/RadialRecipe/Debug.md`). The `grtresna_independent_scalars` path
exists precisely to keep them identical. Any new matter (mass search, `λφ⁴`,
complex field) must be added to **both** sides with matching analytic forms.

## Campaign `ftl_discovery_v2` — first healthy run + a scoring concern (2026-06-10)

Fresh 8-GPU campaign (`runs/grtresna_qd/ftl_discovery_v2/`, `speed_super` 8×8,
`ftl_first`) launched on the vectorized-conversion + corrected-stationarity +
moderated-penalty fixes. This is the first run that scores sanely: gpu_ok
candidates get gradient-rich positives (14→405), `stationary_artifact_penalty`
is 0 on survivors (graded, e.g. −0.955 on the near-stationary eval_011), the
horizon veto fires correctly (eval_026/038 → −1.0), and the archive spreads
across both axes (y-bins 0→7, cells up to [7,7]) — the earlier y-collapse is
fixed.

### Validation of the top elite (eval_000036, score 405, cell [7,7])

Checked in detail because it carries the run's first nonzero
`operational_ftl_solved` (0.717). **Verdict: best *physical precursor* so far,
but NOT confirmed gauge-invariant FTL.** Evidence from `score.json`:

- **Genuinely physical (unlike the old eval_083 lens):** `comoving.stationary =
  false`, `beta_mean = 0.371` (real net shift), final Ham/Mom L2 = 8.2e-4 /
  1.8e-4 (stays on-constraint — not a relaxation transient), `min_theta_plus =
  0.118` (no horizon), survives to t=2.
- **But the 405 is driven by a *coordinate*-speed metric.** `operational_ftl_solved`
  derives from `general_ftl_solved.max_local_speed = 1.287`, explicitly noted as
  a "superluminal **coordinate** channel" — gauge-dependent, a precursor not a proof.
- **The one gauge-invariant test failed its reliability gate:** `geodesic_ftl`
  has `h_quality_ok = False` (null-constraint drift max H=5.4e-4, only 4/5 rays
  reached), so `operational_ftl_geodesic` was correctly zeroed.
- **No sustained signal:** evolved `operational_ftl = 0`, `ftl_persistence = 0`
  (max local speed decays 1.287 → 1.092 over the evolution).
- **Frames** (`shift1_z`, `local_speed_z`, t=0 vs t=2): a real ±0.05 shift dipole
  and a ~1.05–1.10 (on z-slice) superluminal region, both **localized near
  center and weakening in place** — they do not propagate as a warp channel.

So the scoring behaved honestly here: it labeled the result as precursor +
channel progress, *not* certified FTL.

### Concern + plan: gauge-invariant signal is under-weighted

`operational_ftl_solved` (a **coordinate** speed, gauge-dependent) currently
dominates the `ftl_first` objective, while the only gauge-invariant measure
(`geodesic_ftl`, null rays) is frequently rejected for unreliability and
contributes 0. This risks ranking gauge/shift-channel artifacts above genuine
shortcuts. Plan (future work, in priority order):

1. **Fix null-geodesic reliability** so `h_quality_ok` can pass on good
   candidates: tighten ray integration (smaller step / higher-order),
   investigate the constraint-drift source, and ensure all rays reach the
   detector. Until this passes, we have *no* trustworthy gauge-invariant verdict.
   Code: `metrics/null_geodesic.py` (`GeodesicFtlReport`, `h_quality_ok`,
   `n_reached/n_rays`, `max_h_drift`).
2. **Re-weight once geodesic is reliable:** gauge-invariant `operational_ftl_geodesic`
   should outweigh the coordinate-based `operational_ftl_solved`, so the
   leaderboard tracks real shortcuts rather than coordinate channels. Code:
   `metrics/score.py` (`ftl_first` component weights).
3. **Keep `operational_ftl_solved` as a precursor/shaping term only** (lower
   weight), gated as now by non-stationarity, so it guides the search toward
   shift channels without certifying them as FTL.
