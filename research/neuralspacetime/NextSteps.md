Matter-First Metric Discovery: Revised Technical Plan

Based on trajectory_5lump_v1 campaign results (60 evals, best stable score 1039.5,
peak geodesic shortcut 10.82%).

This revision supersedes the prior blueprint. Several originally proposed changes
are abandoned based on empirical evidence showing they would degrade or destroy
the discovered FTL mechanism. The pipeline is restructured around what the data
actually reveals: independent per-lump orbital freedom, static initial data, and
direct HQ promotion of top elites.

Discovered FTL Mechanism (Summary)

The MAP-Elites search has converged on a robust pattern:

  - 5 independent scalar lumps in counter-rotating (all-retrograde) orbits
  - Nested shell structure: tight inner lump (R0~2.5) + wide outer lumps (R0~6-8)
  - Mixed orbital tilts: some near-equatorial, some near-inverted (~165 deg)
  - 3/5 lumps carrying exotic matter (ANEC-violating)
  - Strong z-oscillation (breathing mode) coupling all lumps vertically

This creates a 3D frame-dragging vortex (shift_drive = 0.906) that sustains a
10.82% superluminal geodesic shortcut for the full 16M evolution. The mechanism
is NOT gravitational lensing -- the null-ray probe measures genuine end-to-end
transit faster than flat-space light across identical coordinate separation.

Key performance metrics (trajectory vs previous SH ansatz):
  - FTL hit rate per GPU eval: 54% vs 1.3% (40x improvement)
  - GRTresna solver pass rate: 100% vs 54%
  - Best stable f_geo_peak: 10.82% vs 2.12% (5x improvement)
  - Best stable score: 1039.5 vs 470.6

Revised Pipeline

         [Optimizer Proposes Trajectory Genome (spaces.py, 40D)]
                                │
                                ▼
         [Soft Pre-GPU Filter: density cap only]
           Reject if sum(well_depth) > 0.50
                                │
                                ▼  (Passes Filter)
         [GRTresna Static Initial Data (t = 0, no momentum)]
           Trivial constraint solve, 100% pass rate
                                │
                                ▼
         [Postload Gate: Ham_L2 < 0.01, Mom_L2 < 0.01]
                                │
                                ▼  (Passes Gate)
         [GRTeclyn GPU Evolution (128^3, max_level=1)]
           Smooth-start ramp (tau=2.5), dt_mult=0.02, t_stop=16
                                │
                                ▼
         [4D Evolving Null-Geodesic Probe + Scoring]
                                │
                                ▼
         [HQ Promotion of Top Elites to 256^3]


Phase 1: Stage 0 Search -- Keep Running (NEXT LOGICAL STEP)
============================================================

>>> IMMEDIATE ACTION: Complete the trajectory_5lump_v1 campaign to 200 evals.
>>> Then HQ-promote the top 3 elites.

The current Stage 0 configuration is already highly effective:

  - Search space: 40 dimensions (5 lumps x 8 params each)
  - Batch size: 8 (one eval per GPU)
  - Pipeline: 3 concurrent GRTresna solves feeding 8 GPU slots
  - Throughput: ~2 evals/minute
  - Archive: 8x8 grid on (ftl_geo_evolving, ftl_lifetime) descriptors

No changes to the search space, ansatz, or solver are warranted at this stage.
The 40D independent-lump parameterization is the mechanism -- reducing it would
destroy what makes it work.

What NOT to change (and why):

  1. Do NOT add collision filtering (hard abort).
     - Eval 008 (best score 1166.8) has overlapping lumps at t=0 (margin = -1.73).
     - Close-proximity lumps generate the strongest frame-dragging shear.
     - The postload gate already handles configurations that are genuinely unstable.

  2. Do NOT cap R0 at 6.5.
     - All top-5 configs have at least one lump with R0 > 6.5 (up to 8.0).
     - The nested inner/outer shell structure IS the mechanism.
     - Current [1.5, 8.0] range is validated by the data.

  3. Do NOT switch to co-orbital rings.
     - The FTL mechanism requires independent tilts, radii, and omega per lump.
     - A co-orbital model forces shared parameters that would collapse the
       diversity the optimizer exploits (varied tilts 5-175 deg, varied omega
       -0.25 to -0.99, varied R0 2.5 to 8.0 within a single config).
     - Dimensionality reduction here is mechanism destruction.

  4. Do NOT introduce momentum-solved initial data.
     - Static t=0 gives 100% GRTresna pass rate (vs 54% with momentum solve).
     - The smooth-start ramp builds frame-dragging dynamically during evolution.
     - Reintroducing S^i sources would re-create the 46% failure mode.

  5. Do NOT halve dt_multiplier for Stage 0.
     - Would double wall-clock time per eval (currently ~45s GPU phase).
     - gpu_failed rate is 10% -- not worth halving throughput to save 10% of evals.
     - Reserve CFL reduction for HQ promotion where accuracy matters more.

What CAN be tuned (low-risk):

  - Soft density cap: reject if sum(well_depth) > 0.50. This catches the
    over-massive configurations (postload_rejected mean = 0.473) without
    blocking the top performers (mean = 0.34). Saves ~15% of wasted GRTresna
    CPU time. Implement as early abort, not penalty score.

  - Postload threshold relaxation: Currently Ham_L2 < 0.01. The median
    rejected value is 0.0155, with 68% within 2x threshold. Relaxing to 0.015
    would admit ~40% more candidates to GPU phase. Trade-off: slightly higher
    constraint violation at evolution start, but the CCZ4 damping handles it.

  - Near-miss seeding: Feed the MAP-Elites near-miss pool with mutations of
    the top-3 elites (especially eval 050's all-retrograde pattern). The
    archive has only 6.2% coverage (4/64 cells); biased seeding from known
    winners would accelerate filling.


Phase 2: HQ Promotion (NEXT AFTER 200 EVALS COMPLETE)
======================================================

Once the 200-eval campaign finishes, promote the top 3-5 elites directly to
high-resolution verification. CMA-ES refinement is bypassed -- the trajectory
ansatz already produces configurations that are qualitatively correct; they need
resolution, not parameter tweaking.

Promotion Candidates (current):

  1. Eval 050 (score 1039.5, gpu_ok, survival=1.0)
     - f_geo_peak = 10.82%, ftl_lifetime = 1.0
     - All-retrograde, 3 exotic, nested R0 structure
     - PRIMARY promotion target

  2. Eval 008 (score 1166.8, gpu_failed, survival=0.592)
     - Highest raw score, crashed with NaN in h11 at 59% evolution
     - Counter-rotating (mixed prograde/retrograde), overlapping lumps
     - Hypothesis: higher resolution prevents the shear-zone NaN

  3. Eval 038 (score 734.9, gpu_ok, survival=1.0)
     - Different topology (no geodesic FTL but strong operational FTL)
     - Tests whether the mechanism generalizes

  4-5. Best new elites from evals 60-200 (TBD)

HQ Configuration:

  - Grid: N=256, L_full=128, max_level=3 (effective dx = 0.0625)
  - CFL: dt_multiplier = 0.01 (half of Stage 0, for stability in shear zones)
  - KO dissipation: enhanced (sigma_KO = 0.5) near orbital boundaries
  - Evolution window: t_stop = 32 (2x Stage 0, tests long-term persistence)
  - Geodesic probe: hq mode, directional sweep x/y/z, n_rays=17
  - Frames output: enabled (for visualization of the vortex structure)
  - Keep all plotfiles (no consumer-delete)

Success Criteria:

  - f_geo_peak confirmed > 5% at 256^3 (rules out resolution artifact)
  - ftl_geo_evolving > 0.3 sustained through t=32
  - Constraint growth < 2x over evolution window
  - No NaN/crash through full t_stop

If eval 008 survives at HQ: it becomes the strongest FTL result to date and
validates the close-proximity exotic lump mechanism.


Phase 3: Post-Verification Analysis
====================================

After HQ promotion confirms the geodesic shortcut at high resolution:

  1. Gauge-invariance check: Repeat with different gauge conditions (1+log
     slicing variants, Gamma-driver shift damping) to confirm the shortcut is
     not a coordinate artifact.

  2. Extended geodesic survey: Fire null rays from a dense grid of emission
     points and directions to map the full "FTL volume" -- the spatial region
     where superluminal transit is possible.

  3. Tidal tensor extraction: Compute the Weyl tensor components along the
     fastest geodesic to characterize passenger survivability quantitatively
     (not just the scalar tidal_comfort metric).

  4. Energy budget: Integrate the stress-energy tensor over the domain to
     quantify total exotic matter required and compare to known energy-condition
     violation bounds (Ford-Roman quantum inequality).

  5. Parameter sensitivity: Small perturbations around eval 050's genome to
     map the local fitness landscape -- is this a sharp peak or a broad basin?
     This informs whether the mechanism is fine-tuned or robust.


Phase 4: Future Search Directions (After Verification)
======================================================

Only after HQ verification confirms the physical reality of the shortcut:

  1. Increase lump count (7 or 9 lumps): More lumps may strengthen the vortex.
     Expect diminishing returns but worth a 50-eval test campaign.

  2. Longer evolution (t_stop = 64): Test whether the frame-dragging vortex is
     truly stationary or slowly decaying. Self-bound solitons (boson stars)
     could replace the pump-maintained lumps for genuinely eternal configurations.

  3. Larger domain (L=128, N=256 base): Test whether the shortcut scales with
     propagation distance or is confined to the near-field of the matter
     distribution.

  4. Asymmetric exotic assignment: The search found 3/5 exotic consistently.
     A focused sub-campaign fixing n_exotic=3 and searching only over placement
     could accelerate convergence.

  5. Omega sign constraint: All top stable results are all-retrograde. A
     constrained search fixing all omega < 0 would eliminate half the search
     space and accelerate coverage.
