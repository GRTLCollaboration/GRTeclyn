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


**5. Q-ball trajectory QD campaign — NEXT**
New MAP-Elites campaign with compact Q-ball solitons on retrograde orbits.
*   **Campaign script:** `scripts/campaigns/qball_trajectory/run.sh`
*   **Matter:** boson_star + trajectory ansatz, compact preset (m=1, λ=640,
    μ=85333, ω=0.8), ODE profile seeding, equilibrium amplitude cap.
*   **Constraints:** all-retrograde orbits (`trajectory_retrograde_only=1`,
    implemented in `candidates.py`); sub-luminal speed cap (v_max=0.3c, already
    default).  Together these remove 50%+ of the search space — HQ-validated
    that counter-rotation and superluminal orbits are false-positive generators.
*   **FTL probe:** multi-ray emission sweep (Δt=2, 7 launches) + xyz geodesic
    directions; `general_ftl` objective.
*   **Goal:** find less-dispersive Q-ball orbital configurations that sustain
    FTL longer.  The compact ODE seed + speed cap fix the two known dispersal
    causes (relaxation radiation + superluminal pump mismatch); the search can
    now explore orbit geometry without those confounders.


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