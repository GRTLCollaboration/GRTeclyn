# Campaign results: `qball_traj_fgeo_v1` — evolving-geodesic FTL shortcut search

**Status:** complete — 400/400 evals (run 1: 0–200, warm-started extension: 201–400).
**Finished:** 2026-08-04 22:48 UTC.
**Best score:** 2541.6 (eval 322) — a genuine evolving-geodesic shortcut (24 % path saving — saturating the depth reward) with **25.0 % of the lump matter still confined at the end of the run**.
**Interactive 3D animation of the matter dynamics** (eval 161 recipe, the breakthrough choreography): https://claude.ai/code/artifact/1d2bd619-1cd2-4a70-b32f-3d6a37f60265

---

## 1. What this campaign searched for

A MAP-Elites quality-diversity search over the 39-dimensional space of Q-ball
lump trajectories (5 lumps × per-lump orbit parameters + 4 shared
breathing/bobbing parameters + per-lump exotic-matter dial), scoring each
candidate by how long a *genuine* faster-than-light geodesic shortcut stays
open in the fully evolved spacetime.

- **Objective (`f_geo_max` mode):** score ≈ 10000 × `ftl_geo_evolving`, where
  `ftl_geo_evolving = min((f_geo − 0.001) / 0.199, 1.0) × structural_persistence`.
  `f_geo` is the shortcut depth (fractional path saving of a light ray traced
  through the *evolving* metric, best over the emission sweep) and
  **saturates at 20 %**; `structural_persistence` is the fraction of lump
  matter still confined and coherent at the end of the run
  (density retention × coherence × confined mass fraction).
  Plus small shaping terms; exotic-matter penalty recorded but weight 0;
  pump-energy tax weight 40; graded horizon penalty.
  **Consequence: any depth beyond 20 % earns nothing** — see §3.2, §3.3, §5.
- **Shortcut test:** rays traced through the time-dependent metric between
  antipodal points on the r = 8 measurement shell; `f_geo` = fractional
  saving vs. the flat-space crossing time.
- **Search:** MAP-Elites, 8×8 archive, descriptor mode `ftl_lifetime`,
  batch = 4 (one eval per H100 GPU), pipeline mode; GRTresna solves each
  candidate's constraint-satisfying initial data on CPU, GRTeclyn evolves
  it (CCZ4) on GPU to t = 26.
- **Trajectory model** (`Source/GRTeclynCore/RL/TrajectoryEvaluator.hpp`):
  per lump `r_k(t) = max(R_MIN, R0 + v_rad·t + A_breath·sin(ω_breath·t))`,
  `φ(t) = phase0 + ω_rot·t`, orbital plane rotated by
  R_z(tilt_φ)·R_y(tilt_θ), shared axial bob `z_amp·sin(ω_z·t)`.
  `exotic ≥ 0.5` → phantom (negative-energy) scalar field for that lump.

Run 1 was interrupted once by infrastructure (evals 29–40 lost, recorded as
12 `pipeline_interrupted` markers in `trajectory.jsonl`) and warm-restarted;
the 200→400 extension warm-started from the 9-elite archive
(`QD_RESUME=1 QD_TARGET_EVALS=400`).

---

## 2. Headline results

| eval | score | matter retained at end | shortcut depth f_geo (best emit) | note |
|---|---|---|---|---|
| **322** | **2541.6** | **25.0 %** | 24.0 % (emit t = 10) | champion — best matter retention |
| 226 | 2529.0 | 24.8 % | 28.1 % (emit t = 12) | |
| 272 | 2488.0 | 24.4 % | 25.7 % (emit t = 12) | |
| 161 | 2455.1 | 24.1 % | — (dir pruned) | run-1 champion, breakthrough recipe |
| 179 / 330 | 2397.5 | 23.5 % | 33.3 % (emit t = 12) | best depth × retention compromise |
| 370 | 607.5 | 5.8 % | **38.3 %** (emit t = 12) | **deepest shortcut ever seen**, but matter disperses fast |

(All these depths exceed the 20 % reward saturation, so each row's score is
effectively 10000 × its matter-retention column — see §3.2.)

All top scorers sit in archive cell [7,7], tier 2 ("operational").
Final archive: 9 elites, coverage 14.1 % of the 8×8 grid, 32 near-misses.
Tier census of elites: 7 operational, 1 nontrivial, 1 constructed, 0 rejected.

**Record progression** (running best over the campaign):

| eval | score | matter retained |
|---|---|---|
| 1 | 77.6 | 0.5 % |
| 2 | 141.2 | 2.2 % |
| 3 | 1873.6 | 18.3 % |
| 6 | 1965.0 | 19.3 % |
| 161 | 2455.1 | 24.1 % |
| 226 | 2529.0 | 24.8 % |
| 322 | 2541.6 | 25.0 % |

The search plateaued at ~19 % for 155 evals (evals 6–160), broke through at
161, and the warm-started extension refined that recipe three more times.

---

## 3. Key findings

### 3.1 Exotic-matter quantity gates admission; choreography decides the outcome

Bucketing all 388 scored evals by how many of the 5 lumps are phantom
(negative energy):

| phantom lumps | n | max `ftl_geo_evolving` (saturated depth × retention) | mean |
|---|---|---|---|
| 0/5 | 3 | 0.0 % | 0.0 % |
| 1/5 | 15 | 0.0 % | 0.0 % |
| 2/5 | 77 | 7.9 % | 0.7 % |
| 3/5 | 140 | 18.2 % | 2.9 % |
| 4/5 | 82 | 19.3 % | 3.5 % |
| 5/5 | 71 | **25.0 %** | 6.5 % |

With fewer than two phantom lumps a shortcut *never* appears. But quantity
alone is far from sufficient: among all-phantom configurations, scores span
~80× on identical matter content — the difference is purely how the lumps
move. Zero-shortcut runs also die young (mean structure survival ≈ 3 % vs
≈ 8.5 % for shortcut-bearing runs).

### 3.2 Above 20 % depth, the score *is* matter retention (by construction)

The scoring component saturates the depth signal:
`ftl_geo_evolving = min((f_geo − 0.001)/0.199, 1) × structural_persistence`
(`GEO_FTL_TARGET = 0.2` in `metrics/score/ftl.py`), where
`structural_persistence` = density retention × coherence × final confined
mass fraction. Every top-10 scorer has depth above 20 %, so its saturated
depth factor is exactly 1.0 and its score component is *bit-exact equal* to
`survival` (for eval 322 both are 0.25022202502458873; 67 of 330 nonzero
rows show this equality, and the two correlate at r = 0.752 overall —
partly mechanical, via the saturation).

Consequence: once a configuration crossed 20 % depth, **all further score
progress in this campaign was matter-retention progress** — the search had
zero incentive to deepen the corridor, and any depth beyond 20 % was
invisible to the objective. (An earlier reading of these numbers as
"fraction of the run the shortcut stays open" was wrong: the component
contains no lifetime measurement at all.)

### 3.3 Depth and matter retention trade off

Eval 370 achieves the deepest shortcut on record — a light signal's
effective path is **38.3 % shorter** than flat space — but its matter
disperses fast (only 5.8 % still confined at the end, score 607). The
champion (322) is shallower (24 %) but retains 4× more matter. Eval 179 is
the best compromise (33.3 % deep, 23.5 % retained). Because the depth
reward saturates at 20 % (§3.2), eval 370's extra 18 points of depth earned
*nothing* — the objective paid only for retention. Depth was therefore
never optimized, and the 38 % figure is likely not the ceiling (see §5).

Shortcut depth vs. emission time (later emission → deeper corridor, i.e.
the inward squeeze keeps concentrating negative energy):

| emit t | 322 | 226 | 179 | 370 |
|---|---|---|---|---|
| 0 | 0.183 | 0.202 | 0.193 | 0.202 |
| 4 | 0.210 | 0.226 | 0.222 | 0.260 |
| 8 | 0.232 | 0.254 | 0.279 | 0.323 |
| 12 | 0.224 | **0.281** | **0.333** | **0.383** |

(3 rays per emit time, all reached, none captured. For 226/179/370 the depth
is still *rising* at the last emission time sampled.)

### 3.4 The winning choreography: a slow coordinated inward squeeze

All record-holders share the same qualitative recipe, discovered at eval 161
and refined through 226 → 322:

- **All five lumps phantom** (negative energy) — `recipe_scalar_field_signs = -1 -1 -1 -1 -1`.
- **Nested shells:** inner lumps at r ≈ 2.3–3.6, middle at r ≈ 4–5.3, one
  outer at r ≈ 7–7.8.
- **Barely rotating:** |ω_rot| ≈ 0.005–0.07 (a few tens of degrees of arc
  over the whole run) — rotation is nearly irrelevant.
- **Everything drifts inward.** Eval 161 had 4/5 lumps falling inward at
  mixed rates; eval 226 and 322 push all five inward at a brisker, more
  uniform pace (v_rad ≈ −0.08 to −0.17). This is the main refinement axis
  from 161 → 322 (+0.95 % lifetime).
- **Axial bob** z_amp ≈ 2.6–3.0, ω_z ≈ 1.5–1.6 in all champions.
- **Breathing is optional:** eval 322 nearly switched it off
  (A_breath = 0.013 vs 0.25–0.29 in 161/226) — the shared radial breathing
  mode is apparently not load-bearing.

Physical picture: the coordinated inward drift keeps the negative-energy
corridor between the antipodal shell points fed and concentrated as the
structure slowly implodes; a quarter of the lump matter is still confined
and coherent when the run ends, and the corridor keeps deepening for as
long as the implosion stays coherent (the frozen-probe depth is still
rising at the final frame).

Watch it: the interactive 3D animation of this choreography (eval 161
parameters, play/pause, orbit/top/side views, coherent-window timeline) is
at https://claude.ai/code/artifact/1d2bd619-1cd2-4a70-b32f-3d6a37f60265

### 3.5 Champion time series (eval 322)

From `eval_000322/small_data/ftl_timeseries.dat` (10 coarse samples):
peak f_geo 0.323 at t = 26 (final snapshot, frozen-metric measure), peak
operational fraction f_op 0.217 at t = 16, max local coordinate speed
1.914 c at t = 26, superluminal fraction reaches 1.0 by t = 6.4. No horizon
formed (horizon_penalty = 0); ANEC flag clean in the scoring sense
(anec_condition = 1.0) while the pointwise energy-condition component is
≈ 0 (as expected — the matter is phantom by construction); exotic penalty
recorded at −1.6 (weight 0 in this mode); pump-energy tax −0.36.

---

## 4. Champion configurations (exact values)

Shared parameters (breathing / axial bob):

| eval | A_breath | ω_breath | z_amp | ω_z |
|---|---|---|---|---|
| 161 | 0.253 | 2.774 | 2.985 | 1.585 |
| 226 | 0.288 | 2.728 | 2.600 | 1.605 |
| 322 | 0.013 | 2.877 | 2.793 | 1.513 |

Per-lump (R0, ω_rot, v_rad, tilt_θ, tilt_φ, exotic):

**Eval 322 (score 2541.6):**

| lump | R0 | ω_rot | v_rad | tilt_θ | tilt_φ | exotic |
|---|---|---|---|---|---|---|
| 0 | 5.32 | −0.0301 | −0.1671 | 1.57 | 0.72 | 0.80 |
| 1 | 3.32 | −0.0574 | −0.1216 | 1.34 | 0.43 | 0.56 |
| 2 | 2.28 | −0.0323 | −0.0834 | 1.85 | 5.86 | 0.85 |
| 3 | 4.13 | −0.0248 | −0.1275 | 2.67 | 1.06 | 0.76 |
| 4 | 7.76 | −0.0350 | −0.1276 | 2.14 | 5.51 | 0.62 |

**Eval 226 (score 2529.0):**

| lump | R0 | ω_rot | v_rad | tilt_θ | tilt_φ | exotic |
|---|---|---|---|---|---|---|
| 0 | 6.73 | −0.0049 | −0.1369 | 2.66 | 1.71 | 0.93 |
| 1 | 3.61 | −0.0338 | −0.1238 | 2.14 | 0.30 | 0.86 |
| 2 | 3.64 | −0.0385 | −0.1252 | 2.52 | 6.21 | 0.74 |
| 3 | 4.55 | −0.0172 | −0.1599 | 3.00 | 1.71 | 0.76 |
| 4 | 7.26 | −0.0388 | −0.0089 | 2.45 | 5.63 | 0.61 |

**Eval 161 (score 2455.1, run-1 champion — eval dir pruned, parameters preserved here and in `trajectory.jsonl`):**

| lump | R0 | ω_rot | v_rad | tilt_θ | tilt_φ | exotic |
|---|---|---|---|---|---|---|
| 0 | 4.81 | −0.0227 | −0.1565 | 1.94 | 0.24 | 0.87 |
| 1 | 2.52 | −0.0735 | −0.0832 | 1.59 | 1.54 | 0.61 |
| 2 | 2.45 | −0.0427 | −0.0806 | 2.18 | 6.08 | 0.92 |
| 3 | 4.73 | −0.0202 | −0.1684 | 2.92 | 1.96 | 0.59 |
| 4 | 7.15 | −0.0393 | −0.0201 | 2.33 | 5.45 | 0.62 |

Full override lists (all 39 genome values, full precision) are in the
matching rows of `trajectory.jsonl` (`overrides` field).

---

## 5. Caveats & open questions (for the paper)

1. **Depth was never optimized — worse, it saturates.** The reward caps the
   depth signal at 20 % (`GEO_FTL_TARGET = 0.2`) and multiplies it by
   matter retention, so eval 370's 38.3 % earned no more than a 20 % 
   corridor would have. Depth is also still rising at the last emission
   time sampled (t = 12) and at the final frozen-probe frame — a
   depth-targeted campaign (uncapped reward, no retention multiplier,
   later emission window) would likely find substantially deeper corridors.
   **Update 2026-08-05:** the cap was removed from the scorer for all
   modes (`_geo_magnitude` is now linear and uncapped, as are the
   operational arrival-time components), and a dedicated uncapped
   depth campaign (`qball_traj_fgeo_depth_v1`, objective `f_geo_depth`,
   stop_time 32, emission sweep to t = 18) was launched seeded from this
   campaign's elites. Numbers in this document reflect the capped v1
   scoring that actually governed the search.
2. **Coarse ray statistics.** The evolving-geodesic test uses 3 rays per
   emission time and 7 emission times; the per-timestep series uses 5 rays
   at 10 coarse samples. Adequate for ranking, thin for a publication-grade
   depth measurement.
3. **Single measurement shell** (r = 8) and single resolution; no
   convergence study yet on the champions.
4. **Score = matter retention above 20 % depth** (§3.2) means further
   progress under `f_geo_max` is a structure-retention problem, not a
   geometry problem. Candidate next steps: a pure-depth objective
   (caveat 1), longer evolutions, or per-lump breathing/feeding to prolong
   the coherent implosion.
5. **Phantom matter** violates pointwise energy conditions by construction;
   the ANEC component reported 1.0 under the scorer's definition — worth an
   explicit averaged-condition analysis before making claims.
6. **Late-time shell extraction warnings** (`RadialRecipePlt02560/02600`:
   "no valid samples at any radius") on the final eval's post-processing —
   harmless for scoring but check before reusing late-time shell profiles.

---

## 6. Data availability

| artifact | where |
|---|---|
| Scoreboard, all 400 rows (12 interrupted markers; filter `'score' in row`) | `trajectory.jsonl` |
| Kept eval dirs (full data, params, time series) | `eval_000002, 000141, 000179, 000226, 000272, 000322, 000370` |
| Per-eval shortcut summary | `eval_*/small_data/evolving_geodesic.json`, `ftl_timeseries.dat` |
| Run logs | `../_logs/qball_traj_fgeo_v1.log`, `../qball_traj_fgeo_v1_resume2.log` |
| Interactive 3D animation of eval 161's lump motion | https://claude.ai/code/artifact/1d2bd619-1cd2-4a70-b32f-3d6a37f60265 |
| Static trajectory plotter (regenerate figures for any kept eval) | `grteclyn-wrapper/tests/visualisation/plot_lump_trajectories.py` |

Note: MAP-Elites prunes non-elite eval dirs automatically; eval 161's
directory (and its `lump_trajectories.png`) was pruned at campaign end.
Its full parameter set is preserved above and in `trajectory.jsonl`; the
figure can be regenerated by writing those overrides into a `params.txt`
and running the plotter.
