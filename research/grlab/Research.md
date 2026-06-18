That bath or pool analogy is exceptionally accurate. In physics, this is known
as constructive wave focusing (or spherical/cylindrical implosion). It is
governed by some of the most profound non-linear mathematics in wave mechanics
and General Relativity.

When you disturb water in a circular pool, the waves travel inward. Because the
circumference of the wave circle shrinks as it approaches the center, the energy
is compressed into a smaller and smaller area. The wave amplitude mathematically
grows as 1/\sqrt{r} (in 2D) or 1/r (in 3D). When they meet in the exact middle,
they interfere constructively to produce a massive, localized, non-linear peak—a
giant "splash" [1].

In General Relativity, using ordinary (positive-energy) matter, we can construct
this exact "gravitational splash" [2].

Designing the "Gravitational Bath" in Your Code

To do this, we do not need to invent unphysical boundaries or reflecting walls
on your numerical grid (which would crash the simulation). Instead, we use your
Toroidal Ring Layout (matter_layout = 3) [4]:

                     [ THE GRAVITATIONAL SPLASH ]

                 Lump 1 (Re)             Lump 2 (Re)
                  (ρ > 0)                 (ρ > 0)
                  ┌────┐                  ┌────┐
                  │ Φ₁ │                  │ Φ₂ │
                  └────┘                  └────┘
                     \                      /
                      \   Ingoing Waves    /
                       ▼                  ▼
                              * r = 0  ◄─── Perfect Constructive
                       ▲                  ▲     Interference Peak
                      /   Ingoing Waves    \
                     /                      \
                  ┌────┐                  ┌────┐
                  │ Φ₃ │                  │ Φ₄ │
                  (ρ > 0)                 (ρ > 0)
                 Lump 3 (Re)             Lump 4 (Re)

1.  The Ring as the Transmitter: The ring layout arranges your positive-energy
    boson stars in a symmetric circle around the coordinate origin (r=0) [4].
2.  Coherent Breathing: You configure the scalar momentum (\Pi) so that all the
    stars breathe in phase [2]. They act as a circular array of synchronized
    transmitters, sending gravitational and scalar waves traveling inward toward
    the center of the ring [1, 4].
3.  The Focus at r = 0: As the waves travel inward, their geometric wavefront
    area shrinks. Their amplitude grows. When they meet at the exact coordinate
    origin (r=0), they interfere in perfect phase.
4.  Non-Linear Gravity Takes Over: Because the energy density of the scalar
    field is proportional to the square of the field derivatives
    (\rho \sim \dot{\phi}^2 + |\nabla \phi|^2), the energy density at r=0 spikes
    violently [2].

This extreme, localized energy concentration dramatically warps the surrounding
spacetime. If you look at the visualization frames you just shared, you would
see the conformal factor well (\chi - 1, Group 12) and the local coordinate
speed (c, Group 10) plunge dynamically at r=0 as the waves peak, creating a
massive, transient gravitational potential well solely from wave focusing [4].

Programming This into MAP-Elites

You can configure MAP-Elites to find and optimize the configurations that
produce this "gravitational splash" [4].

1. The Dynamic Scorer (general_ftl Modification)

Instead of searching for FTL, write a custom objective profile that rewards Wave
Concentration:

  - The Fitness: Maximize the maximum central energy density reached during the
    run:
    \text{Score} = \max_t \left( \rho_{\rm req}(t, \vec{0}) \right) \times \text{survival} - \text{penalties}
    To prevent the optimizer from just using brute-force positive mass to create
    a static black hole, you add a penalty if the initial mass is already
    collapsed. The waves must start dilute and physically travel inward to form
    the peak.

2. The MAP-Elites Descriptor Axes

Map the discovered "splashes" over two physical axes [4]:

  - X-Axis — Wave Chromaticity (Frequency): The oscillation frequency of the
    focusing waves.
  - Y-Axis — Initial Ring Radius (R): The distance the waves must travel before
    they meet at the center [4].

What the Search Will Discover

By running this campaign, MAP-Elites will find the exact combinations of scalar
mass (m), self-coupling (\lambda), and ring radius (R) that act as a perfect
gravitational lens, focusing the ordinary scalar waves into a highly
concentrated, non-linear spike [2, 4].

  - The Sub-Critical Splash: If the initial wave energy is just below the
    threshold, the waves will focus at r=0, produce a massive, non-linear peak
    in curvature, and then cleanly bounce and disperse back out to the
    boundaries. This is the exact gravitational analog of your circular pool
    splash [1, 4].
  - The Super-Critical Collapse: If the initial wave energy exceeds the critical
    threshold, the constructive interference at the center will become so dense
    that the self-gravitation cannot be resisted [1, 4]. The minimum lapse will
    collapse to zero (\alpha_{\rm min} \to 0), and your Apparent Horizon finder
    will detect a black hole forming at r=0 purely out of focused,
    positive-energy waves [2].

This is a beautiful, highly dynamic, and 100% physical experiment that can be
run on your current 128³ gridinit pipeline [2, 4].


The short answer is: you can run the initial wave-focusing "splash" and critical
collapse experiments with your current real scalar setup right now. You do not
need to wait for the complex boson star implementation to get these results
[2, 4].

However, if you want to explore sustained, long-term resonance and stable
standing wave cavities, the complex scalar (boson star) implementation will be
necessary [2].

Here is a breakdown of why both models have their place in this research, and
how you can transition between them systematically:

I. What You Can Do Right Now with Real Scalars

The physical "splash" (constructive wave focusing) and critical collapse are
fundamentally wave propagation phenomena. Real scalar fields (\phi) are governed
by the Klein-Gordon wave equation and are fully capable of this [2]:

  - Historical Precedent: Matthew Choptuik’s original discovery of critical
    collapse and discrete self-similarity was achieved using a real, massless
    scalar field [2].
  - The "One-Off" Splash: You can use your current 5-lump real scalar setup
    (grtresna_independent_scalars) in the Ring Layout (matter_layout = 3)
    [2, 4]. If you configure the initial scalar momentum (\Pi) to point inward,
    the real scalar waves will travel toward the center, meet at r=0, and
    produce a massive, non-linear peak in energy density (\rho) and
    gravitational curvature [2, 4].
  - The Result: This will trigger either a dynamic bounce (the wave packet
    disperses) or a lapse collapse (\alpha \to 0), forming an Apparent Horizon
    [2]. This entire lifecycle can be run and visualized on your current
    codebase today [2].

II. Why You Need Complex Scalars (Boson Stars) for "Resonance"

While real scalars are excellent for a "one-off" splash, they suffer from a
major physical limitation when it comes to long-term resonance: radiative
dispersion [2].

  [ REAL SCALAR WAVE RETREAT ]                  [ COMPLEX SCALAR CAVITY (BOSON STAR) ]
  
        r = 0 (The Splash)                            r = 0 (Trapped Energy)
             *                                             *
            / \  (Waves bounce and                        / \  Conserved charge Q
           /   \  leak to boundaries)                    /   \ prevents evaporation;
          ▼     ▼                                       ▲     ▲ phase rotation locks
       Waves disperse;                                  └─────┘ energy in a stable
       energy is lost                                           oscillating cavity

1.  The Energy Leak: When a real scalar wave packet focuses at the center and
    bounces, it naturally spreads out. As the waves travel outward, they leak
    past your matter configuration and escape to the boundaries of your
    simulation grid. The system rapidly loses energy and dampens out.
2.  No Conservation Law: Real scalar fields lack a conserved global charge.
    There is no physical law preventing the field from simply evaporating into
    flat space [2].
3.  The Boson Star Cure: Complex scalar fields introduce a global U(1) symmetry,
    yielding a conserved Noether charge Q (particle number) [2]. This charge
    acts as a physical barrier that prevents the matter from dispersing to
    infinity.
4.  Sustained Resonance: With complex fields, the time-dependence is locked into
    the phase (e^{-i\omega t}), allowing you to construct stable, non-linear
    standing waves (gravitational geons/cavities) [2]. The waves can bounce and
    focus at the center periodically for hundreds of code units without losing
    their energy to the boundaries, allowing you to study long-term orbital and
    wave resonance [2].

III. The Strategic Recommendation

To maximize your research throughput, I recommend an iterative, two-step
approach:

Step 1: The Real-Scalar "Splash" Campaign (Immediate)

Use your existing v21/v22 pipeline with real scalars to run the wave-focusing
sweeps [2].

  - The Run: Lock your layout to the Ring (matter_layout = 3) with
    grtresna_exotic_fraction=0 (no exotic matter) [2, 4].
  - The Objective: Set the optimizer to maximize the central energy density peak
    (\max \rho_{\rm req}(t, \vec{0})) [2, 4].
  - The Payload: This is a low-friction way to generate your first set of highly
    visual, dynamic results—watching the "splash" form, the lapse collapse, and
    the apparent horizon outline the newly formed singularity [2].

Step 2: The Complex-Scalar "Resonance" Campaign (Following Phase 1)

While your real-scalar runs are generating data, complete the compilation of
your Complex Scalar (ComplexScalarField) model in GRTresna and GRTeclyn [2].

  - The Run: Once compiled, transition the ring search to complex fields.
  - The Objective: Evolve the configurations to t = 50M or longer. You will be
    able to observe true, non-linear standing wave resonance, where the complex
    field traps its own radiation in a stable, breathing gravitational cavity
    [2].
