# Spinning-Configuration FTL Hunt (Tipler-Class Cylinders) — Grounded Research Plan

> Rewritten 2026-07-08 against the verified state of the GRTresna→GRTeclyn
> pipeline. Supersedes the earlier speculative five-phase sketch. Companion
> journals: [`../neuralspacetime/MapElitesDynamics.md`](../neuralspacetime/MapElitesDynamics.md),
> [`../rotatingwormhole/RotatingWormholePlan.md`](../rotatingwormhole/RotatingWormholePlan.md),
> [`../grlab/LabJournal.md`](../grlab/LabJournal.md),
> [`../../grteclyn-wrapper/README.md`](../../grteclyn-wrapper/README.md).

---

## 1. Goal and physical framing

**Headline goal:** find **FTL/superluminal configurations** built from spin —
rapidly rotating matter whose frame dragging lets a co-rotating signal beat the
flat-space light-travel time — using the existing CCZ4/AMR evolution +
GRTresna constraint-solve + QD/CMA-ES search pipeline.

This is the same objective as the trajectory/wormhole campaigns (a
gauge-invariant geodesic shortcut, the `f_geo` measurement), attacked with a
different lever: **frame dragging** — the "river of space" around a spinning
column flows azimuthally, and anything riding that river gets carried. Two
exploitation modes, both searched MAP-Elites-style (the MapElitesDynamics
pattern, but with *rotating* lumps/columns as the ansatz):

- **Mode A — Sagnac shortcut:** a co-rotating light signal rides the dragged
  space and beats the flat-space transit time (signal FTL, azimuthal analog of
  `f_geo`).
- **Mode B — frame-dragging accelerator ("spacetime railgun"):** two
  counter-rotating columns at x = ±d drag space in the *same* direction through
  the central channel; a test mass dropped in the channel is carried and
  accelerated by geometry alone — it free-falls (zero proper acceleration)
  while gaining coordinate velocity. The payoff metric is the exit velocity of
  a timelike geodesic through the channel vs its entry velocity.

We are not hunting time machines; CTCs appear below only as the mathematical
extreme end of the same measurement, useful for calibrating how far up the
ladder a configuration got.

### The correct FTL measure for spin (fixing the old plan)

The earlier draft proposed "maximize the region where g_tt < 0 outside the
horizon." That is wrong: with signature (−,+,+,+), g_tt < 0 is the *normal*
state of spacetime. The quantity that actually encodes the spin-assisted
shortcut is the **co-rotating null transit time vs the flat-space baseline** —
the azimuthal analog of the existing end-to-end `f_geo` shortcut fraction.

The graded, measurable FTL ladder (weakest → strongest):

1. **Frame dragging:** nonzero g_tφ (equivalently, azimuthal shift β^φ). Any
   spinning source produces this; magnitude is the cheap descriptor.
2. **Sagnac shortcut (the FTL signal):** co-rotating vs counter-rotating null
   rays around the axis arrive at different coordinate times; the co-rotating
   ray beating the flat-space circumference time by a fraction f is the direct
   spin analog of the `f_geo` shortcut. Gauge-honest, cheap to measure, scales
   with the frame-dragging integral.
3. **Strong tipping:** the ratio g_tφ² / (g_tt · g_φφ) approaching 1 — light
   cones tilted far enough that the co-rotating transit time collapses toward
   zero (arbitrarily large effective signal speed).
4. **Extreme limit (g_φφ < 0):** the closed azimuthal circle becomes timelike —
   the Tipler/CTC regime. Not the goal; flagged only as the asymptote of the
   ladder and an automatic red flag for numerical breakdown checks.

### Reality checks (do not skip)

- **Tipler's infinite cylinder is a mathematical idealization.** Hawking's
  chronology-protection result: CTCs generated in a *compact* region require
  weak-energy-condition violation. A finite cylinder of ordinary matter cannot
  do it. The pipeline *does* have phantom/exotic scalar sectors (used
  throughout the wormhole campaigns), so the compact hunt is not a priori
  impossible in-simulation — but any extreme-ladder claim inherits the full
  exotic-matter penalty and falsification burden. Modes A and B do **not**
  need the extreme regime: frame dragging and channel acceleration exist at
  any spin and simply scale up.
- **Expected regime at survivable parameters** (ω ≈ 0.05–0.2, m = 1–2, field
  amplitudes below collapse): measurable frame-dragging, Sagnac scaling, and
  railgun channel boost — levels 1–2 of the ladder. Level 3 is pursued by
  climbing the amplitude/spin ladder (κ-family style) as far as numerics
  allow. Any level-3+ signal must survive the resolution ladder before being
  claimed.
- **The dominant failure mode is matter dispersal, not gauge collapse.** Every
  campaign to date says so: Q-ball trajectory HQ confinement fell 53% → 23%
  (t=30); the rotating-wormhole phantom cloud dispersed at t ≈ 4.5–5 across all
  κ. Both modes need *sustained* rotation — a shortcut or accelerator that
  exists only while the matter flies apart scores as the eval-118 failure did.
  Matter design is the whole game (§2).

---

## 2. Matter sector — the decisive lever

Angular momentum can be carried two ways. This choice dominates persistence.

### (a) PRIMARY: internal phase winding — matter at rest, field rotating

Complex (phantom or canonical) scalar Φ = f(ρ,z) e^{i(mφ − ωt)}:
|Φ|² is axisymmetric and *stationary*; J_z is carried entirely by the phase
winding; no bulk matter motion at all. This sidesteps the boosted-lump
dispersal problem — there is no kinetic energy trying to unbind the
configuration.

This machinery **already exists and is validated**:

- RotatingWormholeCollapse evolves (φ₁,Π₁,φ₂,Π₂) with `wormhole_azimuthal_m`,
  `wormhole_rotation_omega`, phantom sign flip — built, smoked, and run.
- GRTresna solves the winding ID to Ham < 1%, Mom < 1%
  (`params_rotating_wormhole_test.txt`, κ-family driver
  `scripts/wormhole/solve_kappa_family.{sh,py}`).
- **Headline B5 result:** analytic rotating ID NaN'd at t ≈ 2.8; the *same*
  configuration with GRTresna-solved ID ran stably with J_z conserved
  (77.5677 at t=0 and t=3.5) and constraints *decreasing*. The rotating
  instability was the O(ω) K_ij=0 ID defect, not physics.
- Lesson L3 stands: a single *real* scalar cannot rotate smoothly — cos(mφ)
  modulation is the four-lobe dispersal pattern. Winding complex scalar is
  mandatory for any rotating profile.

For the cylinder: replace the wormhole throat profile f(r,θ) with a
**cylindrical envelope** f(ρ, z) — a column (or z-stack of tori) of winding
field centered on the z-axis. New profile module
`grtresna/profiles/cylinder_winding.py` + envelope painting; the GRTresna
matter model and evolution code are reused as-is.

### (b) SECONDARY: self-gravitating stable-branch boson star / Q-ball stack

The [`SELFGRAV_HANDOFF`](../../grteclyn-wrapper/docs/SELFGRAV_HANDOFF.md)
four-bug fix (profile loader, curved-space Π painting with real α(r), pump
target amplitude, stable-branch φ_c < Kaup) is the best persistence result the
MAP-Elites program has produced: confinement ~0.90–0.97 to t=16 vs 0.58 for
the broken seed. Honest caveats carry over verbatim:

- late-time RMS spreading (2.1 → 4.2 over the last third — "held but leaking");
- **max_level=3 still NaNs at t ≈ 6–9** — the clean result is coarse-grid only;
- validated for a star **at rest**; nothing is known about it under boost.

Use case here: a z-stack of stable-branch stars, each given spin by *internal
winding* (m ≥ 1 rotating boson star / Q-torus), not by orbital motion.
Known physics risk: m ≥ 1 rotating boson stars are unstable to
non-axisymmetric perturbations unless self-interactions stabilize them —
Q-ball-type potentials (already in `grtresna/profiles/qball_ode.py` /
`qball_couplings.py`) are exactly the stabilizing family. Treat "does the
winding Q-torus hold axisymmetry" as an explicit early experiment, watched in
the frames, not assumed.

### (c) FALLBACK ONLY: boosted-lump trajectories

Lumps orbiting at v_t ≈ 0.3–0.5c (the trajectory-ansatz campaigns) are the
*documented dispersal failure*: eval 118 HQ showed the channel peak at t ≈ 19
then decay as the lumps sheared apart; the two gpu_failed evals were
dispersal-driven blow-ups. Given the tangential speeds a cylinder-scale spin
would demand, bulk-motion matter is used only as a comparison arm — never as
the primary carrier of angular momentum.

**Bottom line:** the "spin" of the Tipler column comes from ω·m phase winding
in matter that is spatially at rest. That is the only configuration class the
pipeline has ever shown to both rotate and persist.

---

## 3. What the infrastructure verifiably provides (reuse, don't rebuild)

| Capability | Where | Status |
|---|---|---|
| Winding complex phantom scalar evolution | `Examples/RotatingWormholeCollapse` (`SupportedWormholeLevel.cpp`) | built, smoked, production runs done |
| GRTresna rotating exotic constraint solve | κ-family driver `scripts/wormhole/solve_kappa_family.{sh,py}` | converges Ham/Mom < 1%; dx=0.5 ID produced |
| Amplitude-ladder (κ) ID families | same | ready-made spin/amplitude sweep tool |
| Collapse/spin diagnostics (`min_chi`, `min_θ₊`, barycenter, **J_z**) | `collapse_diagnostics.dat` | wired |
| `spacetime_shear` objective + descriptor | `src/.../metrics/score/objectives.py`, `search/qd_search/descriptors.py` | exists |
| 4D evolving null-geodesic probe | `src/.../metrics/probes/ftl/evolving_geodesic.py` | straight rays only — Sagnac variant is the extension (§5) |
| QD → CMA-ES → HQ staging, postload gate, tier ladder | `scripts/campaigns/`, `projection/postload_gate.py`, `scripts/search/validate_tiers.py` | production-proven |
| Disk-safe launcher + streaming consumer | `scripts/wormhole/wormhole_case.sh` pattern (sidecar, `--keep-last`, frame drain) | L6-compliant |
| Per-lump independent scalars (counter-rotating pairs) | `GRTresnaIndependentScalars` C++ | exists; opposite winding signs = railgun arm |

Paid-for lessons that constrain every design choice (from Debug.md / the plans):

- **L1** frozen prescribed sources are fatal — matter must co-evolve.
- **L2** never evolve finer than the ID file's native dx (solve ID at target dx).
- **L3** real scalars cannot rotate; **L4** winding complex scalar is the cure.
- **L6** always run the consumer sidecar; never the bare binary.
- **Analytic rotating ID is fatal** (t≈2.8 NaN) — GRTresna solve is mandatory,
  not optional.
- **Reward hacking is real** (gw_beam eval 51: score 336 from a constraint
  blow-up, → −116 only after *multiplicative* health gates). Any curvature /
  off-diagonal-metric objective ships with hard gates from day one.

---

## 4. Phase plan

### Phase 1 — Initial data: compact winding column (finite cylinder)

**Primary arm: no boundary-condition surgery.** Build a finite z-extended
winding column inside the existing open-boundary box and let GRTresna solve it:

1. `grtresna/profiles/cylinder_winding.py`: envelope f(ρ,z) =
   A · sech(ρ/σ_ρ)^p · (z-window), painted onto the complex (optionally
   phantom) scalar with phase winding m; two-channel Π from ω, using the
   *solved* lapse where applicable (SELFGRAV lesson — never paint Π with α=1).
2. Solve with the κ-family pattern at the **evolution's level-0 dx** (L2);
   verify the gridinit structurally (axisymmetric |Φ|², both winding channels,
   J_z ≠ 0) with regression tests mirroring
   `tests/grtresna/test_rotating_wormhole_kappa.py`.
3. **Multi-column (railgun/blender) arm:** two parallel columns at x = ±d via
   `GRTresnaIndependentScalars`, opposite (counter-rotating) or equal
   (co-rotating) winding signs. Only after the single column is stable.

**Deferred stretch goal: periodic-z "infinite" cylinder.** Grounded caveats
that demoted it from Phase 1:

- `is_periodic = 0 0 0` is **hardcoded** in
  `src/grteclyn_wrapper/grtresna/solver/params.py`; GRTeclyn/AMReX supports
  periodic geometry, but the wrapper plumbing (solver params, postload gate,
  consumer slicing) all assumes open boxes.
- Physics: an infinite rotating cylinder's exterior (Levi-Civita/Lewis class)
  is **not asymptotically flat**. GRTresna's momentum solve (compact V_i
  ansatz) and the Sommerfeld x/y faces both assume decay — a mixed
  periodic/open elliptic solve with a non-decaying exterior is genuine methods
  work with unclear convergence, not a params tweak.

Attempt only if the finite column produces a tipping signal that is visibly
end-effect-limited.

### Phase 2 — Evolution: dispersal is enemy #1, gauge is enemy #2

1. **Dispersal watch (primary).** The rotating-wormhole high-res result is the
   warning: all κ arms dispersed at t ≈ 4.5–5, κ-independently — signature of a
   numerical/boundary artifact, still unresolved. Before any campaign, run the
   single winding column through the same discrimination protocol: vary box
   size and resolution; if the dispersal time moves, it is numerical. Add an
   AMR-robust volume-integrated ρ diagnostic (the `rho_sum` refined-region sum
   is confounded by de-refinement).
2. **Pump as stabilizer (optional arm).** The corrected RL/PD pump
   (`rl_pump_target_amp`, profile-matched width) is available if the column
   leaks; account for pump energy injection honestly in scoring (NextSteps
   pump-accounting item).
3. **Gauge/damping.** Strong frame dragging stresses the Gamma-driver; expose
   CCZ4 κ₁/κ₂ and shift-driver η as campaign hyperparameters in
   `scripts/campaigns/lib/search_common.sh`. Watch β^φ for runaway — but note
   the B5/C evidence that with constraint-clean ID, rotation at ω ≈ 0.05 is
   *not* gauge-hostile.
4. **Postload gate:** verify `projection/postload_gate.py` tolerances against a
   known-good solved winding column (nonzero momentum is legitimate here);
   adjust Mom thresholds only with evidence, not preemptively.

### Phase 3 — Diagnostics: Sagnac probe (Mode A) + railgun probe (Mode B)

1. **Frame-dragging map (cheap, every run).** Parser in
   `src/grteclyn_wrapper/visualisation/` extracting the azimuthal shift
   β^φ = (x β^y − y β^x)/ρ² on the mid-plane; peak and radial profile logged to
   `small_data/`. Peak |β^φ|·ρ is the `frame_dragging_strength` descriptor.
2. **Sagnac azimuthal geodesic probe (gauge-honest, the Mode A FTL meter).** Extend
   `metrics/probes/ftl/evolving_geodesic.py` with an azimuthal launch mode:
   integrate co- and counter-rotating null geodesics on a ring ρ = R around the
   axis through the evolving 4D metric stack; report
   Δt = t_counter − t_co and the normalized asymmetry vs the flat-space
   circumference time. Ladder levels: Δt > 0 (dragging), the co-rotating ray
   beating flat space (the FTL signal), t_co → 0 (strong tipping). Reuses the
   existing metric-stack cache and trust flags (step-size, drift, NaN checks)
   unchanged.
3. **Tipping-ratio field.** From plotfile metric components, compute
   g_tφ²/(g_tt g_φφ) on the mid-plane (g_tφ = γ_φj β^j etc. from ADM
   variables); log max over the ring family. This is a *coordinate* quantity —
   it ranks candidates but only the Sagnac probe substantiates a claim.
4. **Railgun probe (gauge-honest, the Mode B accelerator meter).** Integrate a
   timelike geodesic dropped at rest in the central channel between the
   counter-rotating columns through the evolving 4D metric stack; report exit
   vs entry velocity (the boost), the proper acceleration along the worldline
   (should be ~0 — the "zero-G catapult" signature), and the boost per unit
   channel length. Explicitly label raw coordinate-velocity numbers as
   gauge-dependent; the geodesic integration with trust flags is the honest
   measurement. A family of drop points and initial velocities gives the
   channel's "acceleration map".

### Phase 4 — Campaigns (QD → CMA-ES), hard-gated from day one

1. **Two objectives** in `metrics/score/objectives.py`, one per mode:
   - **`spin_ftl` (Mode A):**
     `score = (w₁·sagnac_shortcut + w₂·tipping_ratio + w₃·J_z_persistence)
     × health_multiplier + penalties`;
   - **`spacetime_railgun` (Mode B):**
     `score = (w₁·channel_boost + w₂·boost_per_length + w₃·J_z_persistence)
     × health_multiplier + penalties`, with a dispersion gate on the columns
     (the `SCORE_FTL_DISPERSION_GATE` pattern — a boost measured while the
     columns fly apart does not count).

   In both, `health_multiplier = 0` on constraint spike (gw_beam Fix A/C
   pattern: truncate time series at t_spike, zero the physics terms — additive
   penalties provably lose to unbounded exploits). Graded `horizon_penalty`
   retained (columns must not merge/collapse to a BH); exotic penalty applies
   whenever the phantom sector is on.
2. **Descriptors:** `spacetime_shear` (exists) × `frame_dragging_strength`
   (new, from the Phase 3 parser) on the 8×8 archive. Search dimensions:
   column radius σ_ρ, z-extent, amplitude A (κ-style), winding m ∈ {1,2},
   ω ∈ [0, 0.2], potential couplings (Q-ball family), phantom on/off, and for
   the pair arm: separation d and relative winding sign.
3. **Stage 0** at the standard search grid (N=128, L=64, t_stop=16),
   `scripts/campaigns/tipler_cylinder/run_qd.sh` cloned from the qball launcher
   with shared `search_common.sh`. **Stage 1** CMA-ES warm-start from archive
   elites, same config (scores must stay comparable).

### Phase 5 — Validation and falsification (non-negotiable)

1. **Resolution ladder:** promote leaders to HQ (N=256, L=128, t_stop ≥ 30,
   max_level per L2 discipline). The trajectory campaign precedent is binding:
   eval 008 was a low-res artifact; eval 118's channel *decayed* under the
   longer HQ window. Assume the same until shown otherwise.
2. **Boundary/box study** for any persistent signal (the t≈5 dispersal protocol
   from RotatingWormholePlan applies verbatim).
3. **Geodesic trust flags** required on every Sagnac number (drift, step
   sizing, no NaNs); ω=0 control run must show zero asymmetry (built-in null
   test), and a slow-ω run should match the linear frame-dragging expectation
   (built-in calibration).
4. **Tier ladder:** extend `scripts/search/validate_tiers.py` with the
   cylinder checks — T0 constructed rotating ID up through converged geometry;
   no FTL or accelerator claim below the convergence tier.

---

## 5. First milestones (concrete, ordered)

1. `cylinder_winding.py` profile + GRTresna solve of a single winding column
   (m=1, ω=0.05, canonical-first, then phantom) to Ham/Mom < 1%; structural
   regression tests on the gridinit.
2. Coarse evolution smoke (64³, t=8) via a `wormhole_case.sh`-style launcher
   with the consumer sidecar; confirm J_z conservation and no dispersal cliff —
   including the box-size discrimination run.
3. β^φ frame-dragging parser + `frame_dragging_strength` descriptor.
4. Sagnac probe mode in the evolving-geodesic tracer, with ω=0 null test and
   small-ω linear calibration.
5. Counter-rotating column pair ID (`GRTresnaIndependentScalars`, opposite
   winding) + railgun probe (drop-in timelike geodesic through the channel).
6. `spin_ftl` and `spacetime_railgun` objectives with multiplicative gates;
   clone the QD launcher; 50-eval shakeout before any full campaign.

## 6. Honest expectations

At the parameter ranges that survive evolution, the deliverables of the first
campaign generation are:

1. a **map of frame-dragging strength, Sagnac shortcut fraction, and railgun
   channel boost versus (m, ω, amplitude, geometry, separation d)** with
   validated scaling — levels 1–2 of the ladder;
2. the **persistence answer** for winding columns: does phase-winding matter
   hold its rotation where boosted lumps did not? This alone decides whether
   spin is a viable FTL lever in this pipeline;
3. the best accelerator elite: how much velocity can geometry alone impart to
   a free-falling payload before the columns destabilize.

Strong tipping (level 3 — co-rotating transit time collapsing toward zero) is
chased by riding the amplitude ladder toward the collapse threshold with the
graded horizon penalty as the guardrail. The g_φφ < 0 extreme (level 4) is not
a target; if it ever appears it is treated first as a numerical artifact and
must survive the full falsification ladder — resolution, boundary, gauge, and
probe trust — before being reported at all.
