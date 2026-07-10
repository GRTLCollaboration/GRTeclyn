Based on the provided internal reviews, lab journals, and roadmap documents,
here is a synthesis of the Directions for Research Improvements and the
Publishable Results.

The project has achieved remarkable engineering milestones (MAP-Elites searching
over constraint-solved initial data in full numerical relativity), but the
physics instrumentation needs to catch up to the search capabilities before
making definitive physical claims.

Part 1: Directions for Research Improvements

The immediate research directions are structured into four main pillars:
Instrumentation/Calibration, Matter Confinement, Active Control (RL), and
Geometry Expressivity.

1. Instrumentation & Calibration (The "P1-P4" Physics Baseline)

Before any claims of "gauge-invariant FTL" can be published, the 4D
null-geodesic probe must be calibrated and made rigorously gauge-independent.

  - Probe Calibration (Positive & Negative Controls): The probe must be tested
    against a Schwarzschild baseline to confirm it correctly reports a Shapiro
    delay (negative f_{geo}), and against an analytic Alcubierre/Krasnikov
    metric to confirm it recovers a known transit-time advantage.
  - Gauge-Honest Baseline (Gao-Wald standard): The current metric compares
    transit time against a flat-space coordinate distance. This must be upgraded
    to compare against a reference null ray routed through the far-field region
    of the same spacetime to ensure pure gauge-invariance.
  - Wave-Zone GW Extraction: Current GW extraction at radii r=8, 12, 24 for
    sources of size 10-18 is in the "near-zone." The extraction domain must be
    pushed further out (or mesh-refined) to ensure \Psi_4 amplitude falls as 1/r
    before quoting beaming ratios or doing LIGO matched-filtering.
  - Rigorous Convergence: The claim that FTL configs hit a "resolution ceiling"
    must be replaced with a true 3-resolution ladder (e.g., 128^3, 256^3, 384^3)
    to demonstrate Richardson extrapolation and formal convergence order for the
    f_{geo} metric.

2. Matter Confinement & Boundary Physics

The search proved that FTL shortcuts in GR are transient because the exotic
matter always disperses or collapses. Stabilizing the matter is the next
physical bottleneck.

  - Pump-Free Self-Gravitating Boson Stars: The current boson star requires an
    unphysical "pump" to hold together. The focus must shift to truly
    stable-branch stars. This involves fixing the Kaup limit discrepancy (0.62
    vs textbook 0.633), defaulting to a safely stable amplitude
    (\phi_c \approx 0.06-0.07), and stabilizing high-resolution AMR runs with
    Kreiss-Oliger dissipation and CCZ4 tuning.
  - Wormhole Stabilization: Massless phantom scalar wormholes evaporate because
    they lack a confining potential. The fix is moving to a bound-state profile
    (Q-ball or self-gravitating boson star) with a non-zero mass (phantom_mass
    > 0) to prevent dispersion.
  - Boundary Reflection Mitigation (The "Sponge Zone"): To solve the
    accumulation of boundary reflections over long runs, implement a Sponge
    Layer (Kreiss-Oliger dissipation ramped up in the outer 20% of the grid).
    This slows and decays outgoing waves before they can reflect off the
    Sommerfeld boundaries.

3. Active Control (Reinforcement Learning)

Moving from MAP-Elites (fire-and-forget ballistics) to RL requires an agent that
obeys Einstein's constraints (i.e., you cannot artificially "respawn" matter
without breaking the Hamiltonian/Momentum constraints).

  - The "Tractor Beam" (Active Trajectory Control): Use a Python RL library
    (like PPO via Stable-Baselines3) to continuously "pilot" the C++
    TrajectoryEvaluator. The agent observes center-of-mass and dispersion rates,
    and actively tweaks the velocity vector (v_{rad}, \omega_{rot}) and pump
    amplitude to shepherd the matter through centripetal balance.
  - Bicomplex Scaffolding: Train the agent to control only Normal (canonical)
    matter to create a stable gravitational "cage." The agent learns to squeeze
    and rotate the normal matter to trap the Exotic (phantom) matter, milking it
    for an FTL wake without letting it disperse.
  - Physics Accounting: Any active RL "pump" injects a non-conservative source
    term. These runs must be strictly tagged as actuated (engineering) rather
    than physical (pure GR solutions), unless the pump's stress-energy is
    properly accounted for in the CCZ4 matter sources.

4. Geometry Expressivity

  - Puncture Free Data: Expand the GRTresna search space from matter-only to
    matter + Bowen-York punctures. This unlocks genuine Einstein-Rosen bridges,
    multi-throat configurations, and Morris-Thorne style
    exotic-shells-around-throats.

Part 2: Publishable Results

The project has generated highly defensible and novel results. Rather than
publishing a sensationalist "FTL exists" paper (which is already known
analytically to be possible with exotic matter), the documents outline four
distinct, high-impact papers:

Paper 1: Methods — AI-Driven Discovery in Numerical Relativity

  - The Claim: Demonstrating a Quality-Diversity (MAP-Elites) search pipeline
    operating directly over constraint-solved initial data in full 3D numerical
    relativity.
  - The Highlight: The "validation ladder." The paper will detail how the
    pipeline successfully decoupled gauge artifacts from physical geometry. It
    will feature the "false-positive table" (e.g., catching a 24.6\% FTL signal
    that was a low-res artifact, or a 2.26c coordinate speed that was just a
    lapse collapse). This proves the robustness of the methodology.

Paper 2: A Negative Result — The Canonical Matter Barrier

  - The Claim: A systematic, high-dimensional search found zero superluminal
    shortcuts using canonical (positive-energy) matter (0/94 in boson shell
    runs, 0 in Lentz arm).
  - The Highlight: This provides rigorous, empirical numerical relativity
    support for the Gao-Wald energy-condition theorems. It serves as a direct
    empirical challenge to recent literature proposing positive-energy warp
    drives, proving that without negative energy, the optimizer cannot cheat the
    topology of GR.

Paper 3: Quantitative Physics — ANEC Violation vs. Shortcut Depth

  - The Claim: Establishing the first empirical "exchange rate" between Averaged
    Null Energy Condition (ANEC) violation and shortcut depth (f_{geo}) in full,
    non-linear GR dynamics.
  - The Highlights:
      - Universal Transience: Demonstrating that all discovered FTL shortcuts in
        dynamic GR are transient—they eventually close via matter dispersion or
        horizon formation.
      - Ford-Roman Quantum Inequality (QI) Margins: Quantifying exactly how many
        orders of magnitude the necessary negative-energy pulses exceed
        quantum-inequality bounds, converting "warp drive" hype into an honest
        physical limitation constraint.

Paper 4: Rotating Phantom Wormhole Dynamics

  - The Claim: The first robust dynamical simulation of a rotating
    phantom-scalar wormhole collapsing/dispersing, mapping the outcome phase
    space against rotation (\omega), azimuthal modes (m), and amplitude
    (\kappa).
  - The Highlights:
      - Extracting the natural-quadrupole (\ell=2, m=\pm2) gravitational wave
        burst that occurs when a rotating massless-phantom throat undergoes
        dynamic instability and radiates its support away.
      - A methodological breakthrough for the field: proving that using a
        constraint-clean trigger (\kappa re-solve) keeps the Hamiltonian
        residual bounded until horizon formation, directly fixing the "steady
        constraint growth" defect that plagued previous Ellis-Bronnikov collapse
        literature.
