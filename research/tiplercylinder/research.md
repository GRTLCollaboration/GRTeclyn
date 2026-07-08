Extending the existing grteclyn-wrapper framework to model Tipler cylinders—and
the associated multi-cylinder configurations like "spacetime railguns" or
"spacetime blenders"—presents a structured path for investigating rotational
numerical relativity.

To transition these concepts from pen-and-paper theoretical physics into your
existing CCZ4/AMR pipeline, we must address the primary lesson of your prior
campaigns: ansatz and matter-sector design are the dominant levers for physical
discovery, while numerical validation is critical to weeding out gauge
artifacts.

Below is a proposed five-phase research plan designed to integrate cylindrical,
spinning configurations into the current wrapper ecosystem.

Phase 1: Boundary Conditions & Initial Data (GRTresna)

The "infinite" cylinder is a mathematical abstraction. In a 3D numerical grid,
we can model this using a clever coordinate trick, or construct a compact
(finite) rotating cylinder.

1.  The Infinite Cylinder Approximation (Periodic Boundary Conditions):
      - Mechanism: Configure the spatial grid with periodic boundary conditions
        along the z-axis in both GRTresna and GRTeclyn. This mathematically
        mimics an infinitely long cylinder of radius R aligned with the z-axis,
        while keeping the computational domain finite.
      - Action: Modify the params.txt generation in
        src/grteclyn_wrapper/grtresna/solver.py to enforce periodic boundaries
        in the z-direction while leaving x and y open (or using radiative
        boundary conditions).
2.  The Matter Sector (Rotating Boson Stars / Q-Rings):
      - Mechanism: To generate rotation, we cannot use static, real scalars. We
        must use a complex scalar field \Phi with an azimuthal phase winding
        (e^{i m \phi}, where m is the integer angular momentum).
      - Action: Build a cylindrical/toroidal initial data profile in
        grtresna/profiles/ (e.g., modifying the existing boson_star_ode.py or
        shell profiles). Arrange multiple spinning scalar lumps along the
        z-axis, or implement a continuous rotating complex scalar cylinder.
3.  Multi-Cylinder (Railgun/Blender) Layout:
      - Mechanism: Define a twin-cylinder ansatz where two parallel cylindrical
        matter arrays are placed at coordinates x = -d and x = +d.
      - Action: Extend the independent scalar lump bridge
        (GRTresnaIndependentScalars in C++) to allocate matter channels for two
        separate cylinders, assigning opposite velocities or phase winds to
        model counter-rotating (railgun) or co-rotating (blender/shear) setups.

Phase 2: Evolution & Gauge Stability (GRTeclyn)

Extreme rotation leads to extreme frame-dragging (non-zero shift vectors
\beta^i). This places massive stress on the CCZ4 system.

1.  Managing Shift Runaway:
      - Your previous findings showed that extreme runs often suffer from gauge
        collapse or coordinate artifacts when the shift vector runaways (e.g.,
        SH eval 151/101).
      - Action: Monitor the \beta^i evolution closely. You may need to tune the
        Gamma-driver shift condition parameters in GRTeclyn's CCZ4 evolution to
        prevent coordinate singularities from forming outside the matter
        cylinder.
2.  Handling the CCZ4 Damping Terms:
      - High-resolution runs in your self-gravitating boson star tests still
        developed NaNs at t \approx 6\text{--}9 due to stability limits.
      - Action: Calibrate the CCZ4 damping parameters \kappa_1 and \kappa_2
        specifically for rotational shear. Introduce these as configurable
        hyperparameters in the wrapper's search defaults
        (scripts/campaigns/lib/search_common.sh).
3.  The Post-Load Gate adaptation:
      - Action: Ensure the projection/postload_gate.py (which runs a short
        t=0.01 evolution to check L_2 Hamiltonian and Momentum constraints) is
        adjusted to handle the rotational momentum of the cylinders without
        falsely rejecting valid, constraint-satisfying rotating initial data.

Phase 3: Diagnostics & 4D Geodesic Probe Customization

To prove the existence of Closed Timelike Curves (CTCs), frame-dragging, or
railgun acceleration, you must measure them in a gauge-invariant manner.

1.  Measuring Frame-Dragging (\omega):
      - Action: Implement a diagnostic parser in
        src/grteclyn_wrapper/visualisation/ to extract the off-diagonal metric
        components (specifically g_{t\phi} or g_{tx}, g_{ty}) and plot the
        frame-dragging "whirlpool" profiles.
2.  Adapting the 4D Null-Geodesic Probe for CTC Precursors:
      - Current State: The EvolvingMetricField probe traces rays along straight
        lines to find FTL shortcuts.
      - Modification: Create a rotational geodesic tracer. Instead of firing
        rays end-to-end along the z-axis, fire null rays in a circular path
        around the spinning cylinder.
      - Goal: Measure the "arrival time asymmetry" between co-rotating and
        counter-rotating rays. A precursor to a CTC occurs when the coordinate
        transit time for a co-rotating ray becomes negative or drastically
        smaller than the flat-space baseline.
3.  The "Railgun" Particle Acceleration Metric:
      - Action: For the counter-rotating dual cylinder setup, write a metric
        that integrates the acceleration of a test mass dropped into the "river"
        of space between x = -d and x = +d. Compare its coordinate velocity and
        proper acceleration to verify the zero-G catapult effect.

Phase 4: Optimization Campaigns (QD → CMA-ES)

With the parameters and diagnostics established, launch targeted search
campaigns to find the most stable and physically extreme configurations.

1.  Define the Objectives (spacetime_shear and spacetime_railgun):
      - Build new objective functions in
        src/grteclyn_wrapper/metrics/score/objectives.py:
          - spacetime_railgun: Maximizes the velocity boost of test particles in
            the central channel while penalizing matter dispersion and black
            hole collapse (using the existing horizon_penalty to prevent the
            cylinders from merging into a single black hole).
          - ctc_precursor: Maximizes the tipping of light cones (i.e.,
            maximizing the region where g_{tt} < 0 outside the event horizon)
            while enforcing stability.
2.  Stage 0: MAP-Elites (QD) Exploration:
      - Launcher: Create scripts/campaigns/tipler_cylinder/run_qd.sh.
      - Descriptors: Use spacetime_shear and frame_dragging_strength (the peak
        value of the shift vector) as the behavior axes for the MAP-Elites grid.
      - Run with intermediate grid settings
        (N=128, L=64, \text{stop\_time}=16.0) to map out the parameter space of
        cylinder radii, distances, and spin rates.
3.  Stage 1: CMA-ES Basin Tightening:
      - Take the most stable, highest-performing "railgun" or "time-cone tipper"
        elites and run local hill-climbing to find the exact thresholds where
        the physics is maximized before numerical breakdown.

Phase 5: High-Resolution Validation & Falsification

As proven by your previous trajectory campaigns (where a Stage 0 leader like
eval 008 was revealed to be a low-resolution artifact), a rigorous validation
path is non-negotiable.

1.  Resolution-Ladder Replay (Stage 2 - HQ):
      - Promote promising cylinder candidates to N=256, L=128, and
        \text{stop\_time}=30 (or greater) to confirm the effects survive grid
        refinement.
2.  Stability and Geodesic Trust Flags:
      - Ensure that any claim of light-cone tipping or zero-G acceleration
        passes the geodesic trust flags (low coordinate drift, correct step
        sizing, and lack of code NaNs/crashes).
3.  Falsification Ladder Integration:
      - Update scripts/search/validate_tiers.py to include specific checks for
        the cylindrical runs, assigning tiers from T0 (constructed rotating
        initial data) up to T5 (converged rotating geometry).

Summary Checklist for Implementation

| Task                          | File to Create/Modify                  | Purpose                                           |
| :---------------------------- | :------------------------------------- | :------------------------------------------------ |
| **Z-Periodic Boundaries**     | `grtresna/solver.py` & GRTeclyn setup  | Mimics infinite cylinder on a finite grid         |
| **Rotating Matter Profile**   | `grtresna/profiles/cylinder_lumps.py`  | Generates spinning complex scalar cylinders       |
| **Frame-Dragging Diagnostic** | `src/grteclyn_wrapper/visualisation/`  | Visualizes the spacetime whirlpool                |
| **CTC Geodesic Probe**        | `src/grteclyn_wrapper/metrics/probes/` | Traces curved/circular null-rays around cylinders |
| **Railgun Objective**         | `metrics/score/objectives.py`          | Rewards acceleration, penalizes merger/collapse   |
| **Campaign Launchers**        | `scripts/campaigns/tipler_cylinder/`   | Coordinates the Stage 0 & 1 search runs           |

