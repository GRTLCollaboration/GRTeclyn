**NEXT STEP — `qball_traj_spiral_v2` is COMPLETE (200/200). Improve the gated search.**

The dispersion-gated campaign finished cleanly (~1 h 45 m, 200/200, no stalls, GPUs
drained). Best **603.39** (eval 118, cell 7,7, operational), 9 elites, **14% coverage**.
The gate works: all operational elites now hold `anec`=1.0, `tidal`=1.0, high
`constraint_growth` (0.91–0.98) — no more disperse-then-warp cheats. But the search
**converged early** (best hit by iteration 15/24, coverage frozen at 14%) and 65 of 200
evals were rejected (50 postload + 15 GRTresna). Priorities:

1. **Break the coverage plateau (14% → higher).** All non-trivial elites sit in the top
   lifetime row (y=7); the low-x cells and the interior of the archive are empty. The
   spiral generator keeps proposing the same persistent-lump basin. Try (a) larger
   MAP-Elites mutation σ / restart injection once coverage stalls 3 iters, (b) widen the
   `ftl_peak_strength` descriptor binning, or (c) seed a few interior cells from v1's
   high-op_ftl (dispersed) configs to force exploration off the y=7 ridge.

2. **Cut the 25% postload rejection rate.** 50/200 evals died at the GRTresna postload
   gate (Ham/Mom ≤ 5%). Either pre-filter candidates by predicted constraint quality in
   `candidates.py`, or relax the gate to ~7–8% and let the scorer's `constraint_growth`
   term demote the marginal ones — 25% wasted GPU-lease time is the biggest throughput leak.

3. **Chase op_ftl, not just ftl_peak.** The score leader (eval 118) has op_ftl only 0.347
   while the `max_local_speed` champion hit **2.006 c** (eval 44) at a much lower total
   score. Consider raising the `operational_ftl` weight in `_general_ftl_total` now that
   it's gated, so the archive optimum tracks a genuinely traversable channel.

4. **HQ-promote eval 118.** Re-run the top elite at 256³ / max_level=3 / t_stop=30 to
   confirm the FTL channel + confinement survive at production resolution (as done for
   `trajectory_5lump_v1` eval 122). Run dir kept: `runs/grtresna_qd/qball_traj_spiral_v2/eval_000118`.

Optional (carried from v1, still open): `confinement_min_frac` floor on the operational
tier; raise `survival` weight in `_general_ftl_total`.

---

**~~Relaunch QD as `qball_traj_spiral_v2` with the dispersion gate~~ — DONE (2026-07-01, 200/200)**

`qball_traj_spiral_v1` stopped at **73/200** evals. Top hits (evals **56/46/40**, scores
~1100, tier operational) saturated `operational_ftl` + `ftl_persistence` but kept only
**28–40% confinement** — matter disperses while a late coordinate FTL channel opens.
The dispersion gate below now fixes the scoring; the campaign is ready to relaunch.

Run: `bash scripts/campaigns/qball_trajectory/run.sh` (defaults to `qball_traj_spiral_v2`,
39D spiral space, `SCORE_EXOTIC_PENALTY_WEIGHT=0.2`, `SCORE_FTL_DISPERSION_GATE=1.0`).
Compare the archive: it should shift toward **persistent lumps + FTL**, not dissolving
clouds. Soften with `SCORE_FTL_DISPERSION_GATE=0.5` if the full gate starves the search.

Note: 5-lump confinement uses one barycenter shell — interpret **time-drop** in
confined fraction (e.g. 56%→28%), not absolute level alone.

Optional follow-ups (not yet done): `confinement_min_frac` floor on the operational
tier; raise `survival` weight in `_general_ftl_total`.


**~~Strengthen matter-dispersion penalty in `general_ftl` scoring~~ — DONE (2026-07-01)**

`qball_traj_spiral_v1` leaders banked most of their score on `operational_ftl` +
`ftl_persistence` at only 28–40% confinement.  Root cause: `structural_persistence`
(`confinement × density_retention × coherence`) gated **`survival`** and
**`ftl_geo_evolving`**, but **not** those two dominant coordinate-FTL terms.

Fix: `compute_ftl_components` now scales `operational_ftl` and `ftl_persistence` by
`structural_persistence` (same as the geodesic terms), controlled by
`SCORE_FTL_DISPERSION_GATE` (via `scorer.py`): `1.0` = full gate (default), `0.0` =
legacy ungated, partial strengths interpolate
`(1-s) + s·structural_persistence`.  Defaults to 1.0 when the matter series is
unavailable, so legacy episodes are untouched.  Guarded by
`tests/metrics/score/test_ftl_dispersion_gate.py`.  Exposed in
`qball_trajectory/run.sh`.

Payoff: MAP-Elites no longer rewards FTL channels that only exist after the Q-balls
fly apart; the archive shifts toward persistent lumps + FTL.


**~~Reduce the exotic-matter penalty in the Q-ball MAP-Elites objective~~ — DONE (2026-07-01)**

`SCORE_EXOTIC_PENALTY_WEIGHT` env (default **0.2** in `qball_trajectory/run.sh`) scales
`exotic_penalty` in `general_ftl` / `ftl_first` / `robust_ftl` via `objectives.py`.
Campaign `qball_traj_spiral_v1` explores exotic-rich configs; eval 5 uses 4/5 exotic lumps
at tier operational.


**~~Add radial-spiral trajectory mode (`v_rad`)~~ — DONE (2026-07-01)**

Per-lump `trajectory_lump{k}_v_rad` in C++ (`TrajectoryParams.hpp`,
`TrajectoryEvaluator.hpp`: `r_k = max(r_min, R0 + v_rad*t + breath)`), Python
search space (`spaces.py`, 39D pinned), joint speed cap + inward spiral clamp
(`candidates.py`, `trajectory_r_min=0.1`).


---

Here is the plain English breakdown of why this idea is so powerful and how it
changes the game.

To understand it, we need to think of our FTL (Faster-Than-Light) effect not as
a solid pipe or a metal bridge, but as a wave in the ocean.

The Problem: We’re only sending one surfer at dawn

Right now, the simulation works like this: the moment we hit "start" on the
experiment (which we call t=0), the computer immediately fires a single test
light-ray (our "probe") across the grid to see if it reaches the finish line
faster than normal light.

But remember what we just discovered: the FTL effect is caused by the matter
dispersing and spreading out. At t=0, the matter hasn't spread out yet! The
"warp" hasn't formed. Our test light-ray is being forced to travel through the
space before the FTL highway is actually built, and it only gets the FTL boost
for the tail-end of its journey.

We are basically sending a surfer out into the ocean at exactly 6:00 AM, and
measuring how good the waves are. If the perfect FTL wave actually swells up
at 6:10 AM, our surfer misses the best part of it.

The Solution: Send a surfer every few minutes

Instead of firing just one test ray at the very beginning, the new idea is to
fire a continuous stream of them—say, one every 2 seconds.

Here is exactly how this helps us turn a random discovery into a working FTL
engine:

1. Finding the "Sweet Spot" By firing a ray at 0 seconds, 2 seconds, 4
seconds, 6 seconds, etc., we can watch the FTL highway open and close in
real-time. We might find that the ray fired at t=6 arrives massively faster than
the one fired at t=0. This tells us exactly when the gravitational tailwind is
at its absolute strongest.

2. Writing the "Train Schedule" of the Warp Highway Right now, we know the FTL
effect lasts for about 16 "code units" (think of these as seconds) before dying
out. But we don't know the exact shape of that lifespan. Does it slowly ramp up
for 8 seconds and fade for 8 seconds? Does it violently snap open at 2 seconds
and slowly bleed out? Firing continuous rays acts like a radar ping, giving us a
perfect, second-by-second map of exactly how long the window stays open.

3. The Ultimate Goal: Pacing the Engine (The "Pulse" Drive) This is the most
exciting part. If we map this out and discover that a single burst of exotic
matter creates a "surfable" FTL wave that lasts exactly 16 seconds... we now
know how to build a permanent highway.

Instead of trying to keep one lump of matter perfectly stable forever (which
makes the simulation crash), we can build an engine that simply "pulses" or
injects a new burst of exotic matter every 15 seconds.

  - Pulse 1 creates a wave.
  - Just before Pulse 1 dies, Pulse 2 fires, creating the next wave.
  - Just before Pulse 2 dies, Pulse 3 fires.

By knowing exactly when to fire the next pulse, a spaceship (or a light beam)
could just surf from one dispersing wave to the next, forever.

In short: Firing continuous test rays allows us to stop treating the FTL effect
like a static, frozen object, and start treating it like a rhythm. Once we know
the rhythm, we can build an engine that beats to it, creating a permanent,
stable FTL corridor out of transient, dying waves.


**4. ~~Rethink the FTL Probes for Transients~~ — DONE (2026-07-01)**
Continuous null-ray emission sweep implemented and validated.  Fires a ray fan
every Δt=2 code units (7 launches) and maps f_geo(t_emit).  First result on eval
122 compact: monotonic decay 9.38%→5.89%, peak at t_emit=0, 100% FTL lifetime.
See MapElitesDynamics.md top section.


**5. ~~Q-ball trajectory QD campaign~~ — STOPPED (2026-07-01, 73/200 evals)**
Launched as `qball_traj_spiral_v1`; stopped manually. Top-3 kept: evals **56, 46, 40**
(scores 1140 / 1102 / 1097). See [MapElitesDynamics.md](./MapElitesDynamics.md) top section.
Relaunch blocked on dispersion-penalty fix (see top of this file).


**6. Pivot from "Warp Ship" to "Warp Highway" (or Stargate)**
Currently, the MAP-Elites scoring function heavily penalizes dispersion: the overall score is multiplied by `structural_persistence`. If the matter drops to 10% density, the FTL score is crushed. 
*   **Action:** Fork the objective function. Keep one search looking for compact "ships," but create a new objective (`standing_channel_ftl`) that **removes the structural persistence penalty**. Reward configurations that turn an initial dense seed into a stable, grid-spanning "highway" of diffuse FTL medium.


**7. Fix AMR (Adaptive Mesh Refinement) Tracking**
Because the smooth boson-star fields do not trigger the default `regrid_threshold=0.02`, the Adaptive Mesh Refinement is failing to engage. Tune the regrid tagger to tag on matter directly, lower the threshold, or increase the base resolution ($N$).


**8. ~~Self-gravitating soliton trajectories~~ — DEPRIORITIZED**
Originally proposed to replace the pump-spotlight with freely orbiting boson
stars.  Deprioritized because: (a) the pump/trajectory potential is the control
mechanism — without it, matter follows geodesics with zero steerability; (b) the
two main dispersal causes (sech seed radiation, superluminal orbits) are now
fixed by the ODE profile + speed cap; (c) the compact-seed run already shows
31% confinement and 100% FTL lifetime with the pump model.  The productive path
is improving confinement within the pump framework (stronger binding, tuned orbit
radii, AMR tagger) rather than discarding the control mechanism.