# Research Roadmap

## Purpose

This document is the top-level scientific roadmap for the repository. It
separates demonstrated results from hypotheses, identifies the gates that must
be passed before stronger claims are made, and links each direction to its
detailed working plan. Implementation details and dated campaign logs belong in
the linked documents, not here:

- [Matter-first spacetime search](./neuralspacetime/NextSteps.md)
- [Trajectory campaign results](./neuralspacetime/MapElitesDynamics.md)
- [Static wormhole dynamics article](./wormholedynamics/article.md)
- [Rotating wormhole program](./rotatingwormhole/RotatingWormholePlan.md)
- [Orbital-pump / collapse-on-command lab record](./rotatingwormhole/OrbitalPumpPlan.md)
- [Rotating wormhole manuscript](./rotatingwormhole/article/research.tex) and its
  [referee-style review](./rotatingwormhole/NextSteps.md)
- [Natural \(m=2\) campaign](./rotatingwormhole/PRDM2CampaignPlan.md)
- [Spin, Sagnac, and finite-cylinder program](./tiplercylinder/research.md)
- [Active-control research](./RL/research.md)
- [GW-beam lab journal](./grlab/LabJournal.md)
- [Pipeline capabilities and operations](../grteclyn-wrapper/README.md)

The repository now contains a capable discovery and falsification pipeline:
GRTresna solves the constraints, GRTeclyn evolves the resulting spacetime,
MAP-Elites and CMA-ES search matter configurations, and 4D geodesic and
\(\Psi_4\) probes test the outcome. The next phase should use that pipeline to
answer sharply defined physics questions rather than treating a high search
score as a discovery.

## Verified state

| Area | Best supported result | Current limit |
|---|---|---|
| Matter-first null shortcut | `eval 144` gives \(f_{\rm geo}=7.96\%\) in the HQ 4D trace, with 5/5 rays and acceptable null-Hamiltonian error | Transient, approximately 91% exotic matter, structural persistence about 0.63, and late-time trapping |
| Trajectory shortcut | `eval 118` gives a 13.0% end-to-end HQ 4D shortcut; frozen slices peak near 22.8% | Confinement falls from about 53% to 23% by \(t=30\); the channel peaks and decays |
| Static Ellis--Bronnikov dynamics | 3D evolution reproduces expansion/collapse branches and a coherent \(\ell=2\) signal with inner-radius propagation near \(c\) | Trigger, horizon, convergence, and finite-radius waveform systematics still limit the strongest manuscript claims |
| Rotating wormhole collapse-on-command | Early support cut of a rotating (m=1) Q-torus throat produces deep collapse (\(\min\alpha\sim10^{-3}\)), a trapped-surface proxy, phantom rebound, and an \(m=0\) \(\ell=2\) burst with **no imposed kick** — reproduced at three resolutions and in an \(L=128\), \(t=80\) box; \(M_{\rm ADM}=0.209\); manuscript drafted | Engineered (support coupling is a nonconservative source); only the *early* cut works (no-ramp / late-ramp arms NaN at \(t\approx19\)); near/intermediate-zone extraction; coordinate \(\theta_+\) proxy, no AH finder; non-rotating control was pathological (killed \(t\approx13.5\)); global \(Q/J\) budgets do not close |
| Rotating passive confinement | Genuine 2D spinning Q-torus eigenstate (profile 4) roughly **doubles** charge-retention half-life (\(\approx27\text{--}28\) vs 13--16) and removes the \(t\approx13.5\) blow-up | Still a flat-space eigenstate — it slowly spreads once gravity is on; phantom self-gravity actively unbinds it |
| Active pump support | Strong PD pumping converts a passive blow-up into a bounded, controlled evolution; measured pump work is negligible (\(\sim10^{-6}\)) — the pump is a transient *igniter*, not an energy prop | The earlier "stable config B" claim was **retracted**: at \(L=128\), \(t=50\) it disperses (confinement \(0.75\to0.004\)) with or without support; supercritical constellations collapse to a BH at every pump setting |
| Directional GW search | Canonical Q-ball configurations show directional preference that persists at HQ | HQ power is weak and the \(1/r\) wave-zone gate fails at the available extraction radii |
| Search methodology | Multiplicative health gates reject collapse-driven \(\Psi_4\) reward hacking; QD \(\rightarrow\) CMA-ES \(\rightarrow\) HQ is operational | Search-resolution leaders frequently fail or weaken under HQ replay |

These results support a **transient warp-precursor** claim, an **engineered
collapse-on-command GW experiment** on a rotating throat, and a productive
numerical-relativity methodology. They do not yet support a stable warp drive,
a new general class of FTL metric, a passively stable rotating wormhole, a GW
laser, or a macroscopic spacetime accelerator.

## Central scientific finding: confinement is the bottleneck

The independent campaigns now point to the same result:

1. Exotic trajectory lumps produce a real null-timing advantage while their
   matter distribution spreads.
2. A mass term delays phantom-field dispersal but does not turn a painted
   profile into an eigenstate.
3. A spherical Q-ball twisted into an \(m=1\) torus is not a stationary
   solution of the rotating Einstein--Klein--Gordon system; increasing the
   Q-ball frequency did not improve throat confinement.
4. Solving the **genuine 2D spinning eigenstate** (rather than twisting a 1D
   profile) roughly doubled charge retention and removed a violent blow-up —
   direct evidence that a self-consistent bound state, not more forcing, is
   what confines.
5. The PD pump does negligible net work: it can *shape* a configuration into a
   basin (igniter) but cannot substitute for gravitational binding — tighter
   gains do not delay collapse, and more pumped matter collapses *earlier*.
   The "stable pumped throat" (config B) was retracted once run longer in a
   bigger box: it disperses with or without support.
6. Canonical matter can be held longer by the trajectory pump, but the resulting
   gravitational field is weak and the pump makes the system engineered rather
   than self-supported.

The next decisive question is therefore not “can the optimizer increase the
score?” It is:

> Can a constraint-clean, dynamically evolved matter configuration retain the
> source of a verified spacetime effect for many crossing times without an
> unaccounted external source?

There are two scientifically honest responses if the answer remains no:

- characterize the transient as the result, including its lifetime and
  radiation; or
- introduce a physical or explicitly accounted active support mechanism and
  label the result accordingly.

Dispersion should not be hidden by a coordinate observable or a score. It is a
primary measured outcome.

**Cross-cutting objective: remove the artificial pump.** The forcing pump is a
useful control-system probe, but as a stabilizer it injects an unaccounted
source and undermines every downstream claim — and the July 2026 measurements
show it *cannot* stabilize anyway: its net work is negligible and no gain
setting prevents a supercritical configuration from collapsing or a
sub-critical one from dispersing. The most valuable single result
this program can produce is a *self-consistent bound state* that confines its
own source with no pump. Two routes lead there, both under Priority B: a
stationary rotating eigenstate (B1) and a composite normal-plus-exotic bound
object (B2). Either would convert transient results across the FTL and wormhole
lines into stable ones simultaneously, which is why passive confinement is the
program's highest-leverage physics goal. The binding caveat below is that a
bound object must be *solved as one state*, not assembled from separately-solved
pieces.

## Priority A: complete the strongest publishable results

### A1. Static wormhole collapse and gravitational-wave signal

**Question:** What is the converged 3D gravitational-wave signature of the two
Ellis--Bronnikov instability branches?

**Next action:** finish the controls and numerical validation identified in
[the article](./wormholedynamics/article.md): use a constraint-clean
out-of-equilibrium trigger where possible, distinguish the throat minimal
surface from a trapped surface, add a genuine apparent-horizon measurement,
strengthen the resolution study, and quantify finite-radius extraction error.

**Decision criterion:** submit only when the waveform and branch assignment
survive the trigger control, resolution comparison, and wave-zone checks.

**Why first:** the core simulation and manuscript already exist. Closing its
validation gaps has higher near-term value than opening another search space.

### A2. Matter-first transient warp precursor

**Question:** Which parts of the measured timing advantage survive changes in
resolution, box size, endpoint placement, emission time, and matter
confinement?

**Next action:** center the paper on `eval 144` as the reproducible HQ result and
use `eval 118` as a complementary transient-dynamics case. Report the full
world-tube history \(f_{\rm geo}(t_{\rm emit})\), null-constraint error,
confinement, energy-condition cost, and trapped-surface status. Do not promote
the 17.7% search-tier emission-sweep value as the headline result.

**Decision criterion:** claim a transient precursor only if the end-to-end
timing advantage remains positive across a resolution/box control and is
measured before numerical or causal boundary contamination.

**Deliverable:** a paper about matter-first discovery and 4D falsification, not
a claim of a functional transport drive.

### A3. Rotating wormhole collapse-on-command manuscript

**Status: this is now the most advanced rotating result — a drafted
manuscript, not a pending pilot.** The
[lab record](./rotatingwormhole/OrbitalPumpPlan.md) closed most of the
referee-style blockers in July 2026: causality controls (early support cut
selects the resolvable deep-collapse branch; no-ramp/late-ramp arms NaN),
attribution of the early \(\theta_+\) crossing to the throat minimal surface,
\(m\)-spectrum (\(\gtrsim99.998\%\) of \(\ell=2\) power in \(m=0\)), ADM mass
and angular momentum, a three-resolution Richardson ladder (convergence is
quantity-dependent), an \(L=128\), \(t=80\) intermediate-zone waveform check
(overlap 0.993, 12.3% residual), and an energy budget showing the drift is not
pump injection. The [draft](./rotatingwormhole/article/research.tex) already
carries the required hedges.

**Next action:** finish, do not extend. The named remaining gaps are:
(i) a clean non-rotating control — the first attempt was pathological and was
killed at \(t\approx13.5\); retry only with a *different* setup (lower AMR
depth, tighter tagging, or a better-matched ID), not more GPUs on the same
trajectory; (ii) a production apparent-horizon finder before any horizon mass
or spin is quoted; (iii) wave-zone or characteristic extraction before any
precision energy, quasinormal-frequency, or detector claim.

**Decision criterion:** submit with the current claim boundary — an engineered
dynamical experiment, early-cut trigger only, proxy horizon, intermediate-zone
waveform. Do not let the remaining gaps inflate into a new campaign; they are
either cheap (documentation of the AH-finder gap) or explicitly out of scope
for this paper (CCE).

**Follow-up, gated:** the [natural \(m=2\) hunt](./rotatingwormhole/PRDM2CampaignPlan.md)
(a genuinely non-axisymmetric burst) and the Phase-7 scale-matched redesign
(bare-mass throat \(b_0\approx2\text{--}3\) hugged by a ring of compact,
*bound* Q-tori, tuned to the critical edge so a support cut collapses instead
of dispersing) are the next experiments on this line — start them only after
the manuscript is out.

## Priority B: attack the confinement problem

### B1. Stationary rotating Einstein--Klein--Gordon solve

This is the highest-value fundamental direction, and its first half is already
done: the **flat-space** 2D spinning eigenstate
\(\Phi=f(\rho,z)e^{i(m\varphi-\omega t)}\) is solved (bordered Newton +
amplitude continuation, residual \(\sim10^{-9}\), grid-convergent), painted via
the profile-4 pipeline, and verified to roughly double charge retention with
no blow-up ([lab record, Phase 8](./rotatingwormhole/OrbitalPumpPlan.md)). The
residual slow spread has two identified causes: the eigenstate ignores
gravitational back-reaction, and phantom self-gravity is repulsive.

**Next action:** the two levers named in the lab record — solve the
**self-gravitating** rotating eigenstate (\(f(\rho,z)\), metric, lapse,
frame-dragging shift, and eigenfrequency together, by continuation from the
non-rotating solution in small rotation steps), and evolve a **normal**
(non-exotic) torus as a control to isolate the phantom-gravity contribution
(needs a non-exotic complex evolution branch).

**Decision criterion:** require a small scalar-equation residual at \(t=0\),
bounded constraints, and flat \(Q_{\rm sphere}(t)\) through \(t=30\). A slower
leak is not equilibrium.

**Outcome if negative:** a documented non-existence or instability boundary
over the searched family is still publishable and is more valuable than another
painted-profile search.

### B2. Composite normal-plus-exotic bound object ("sandwich" cage)

This is the flagship attempt at pump-free confinement: a stable canonical
Q-ball/boson-star shell caging an exotic core, so that a physical equilibrium —
not a forcing term — holds the source together. If it works it converts the
transient results across the FTL and wormhole lines into stable ones at once, so
its importance is high. It is a multi-week methods bet, not a quick win, and it
must survive one decisive caveat.

**The caveat that decides success.** A normal shell painted around an exotic
core, from two separately-solved profiles, is generically *not* a solution of
the coupled two-sector-plus-gravity system. That is the same failure that made
the painted rotating torus disperse: a state that is bound in isolation still
radiates its non-eigenstate modes when it is not an eigenstate of the *combined*
system. Two related warnings: opposite gravitational signs do not imply
attraction (a negative-energy core is gravitationally repelled, not pulled in),
and even a genuine equilibrium may be a saddle (like the Ellis--Bronnikov
throat), i.e. exists but is unstable.

**Next action:** solve the smallest concentric canonical-plus-phantom
configuration as **one self-consistent state** (both matter sectors and the
constraints together, in the bicomplex model), not as superposed profiles. Then
run a one-dimensional radial perturbation scan before any CMA-ES or MAP-Elites.
Track the two charge budgets separately, the lowest radial normal mode, and the
total ADM mass.

**Decision criterion:** proceed to 3D optimization only if the object (i) is
constraint-clean at \(t=0\) as a single solve, (ii) retains both sectors with
flat \(Q_{\rm sphere}(t)\) over several crossing times, and (iii) shows no
growing radial mode — all with the pump off. Only then test whether it also
opens a 4D shortcut.

**Avoid:** assuming a canonical shell "solves" stability, superposing two
solved profiles and calling the result a bound state, or labeling a tuned
transient a stable wormhole.

### B3. Active support as a controlled experiment

RL is useful only after a deterministic controller and an honest source budget
are established. Overwriting matter or spawning a fresh lump during evolution
breaks the Einstein constraints unless the geometry and source terms are
updated consistently.

The PD baseline has now been characterized and its result reframes this
direction: the pump does **negligible net work** (\(\sim10^{-6}\); the
`pump_work` diagnostic exists, column 23) and acts as a *shaping* force, not
gravitational support. It cannot delay collapse of a supercritical
configuration at any gain, and adding pumped matter makes collapse *earlier*.
Its one demonstrated success is as a **transient igniter** — shepherding a
sub-critical constellation into a longer-lived basin, after which it can be
switched off. Confinement itself is set by the mass balance of the
configuration, not by the controller.

**Next action:** treat the controller question as "what basins can an igniter
reach?", not "how much can a pump hold?". Train RL only if it beats the PD
igniter under the same action limits and the same `pump_work`/charge
accounting. See [the RL plan](./RL/research.md).

**Decision criterion:** report the Pareto frontier among shortcut duration,
confinement, injected work/charge, and constraint growth. The correct endpoint
is an **actively supported spacetime**, not a self-sustaining one.

### B4. Transient relay without illegal matter injection

The “pulse drive” remains worth testing, but a second pulse cannot simply be
pasted onto an evolved grid. A clean first experiment places all relay stages
in constraint-solved initial data and uses their initial positions and phases
to stagger the transient response.

**Next action:** infer each candidate's response kernel from the continuous
emission sweep, then optimize a two-stage relay before attempting a longer
chain. Compare the combined result with the linear superposition predicted by
the single-stage kernels.

**Decision criterion:** require a continuous interval of trusted
\(f_{\rm geo}(t_{\rm emit})>0\), bounded constraints, and no loss of source
confinement beyond the known single-stage behavior. If the second stage merely
adds more exotic matter without extending normalized lifetime, stop.

## Priority C: exploit the platform in lower-risk directions

### C1. Frame dragging and Sagnac transport

Use finite, asymptotically flat winding configurations rather than an infinite
Tipler cylinder. The natural building block is now the validated profile-4
spinning Q-torus eigenstate from the rotating-wormhole program: a z-stack of
3--5 tori, with inter-torus spacing, per-torus \((\omega, m)\), and the
winding-sign pattern as search dimensions. Implement the azimuthal
null-geodesic controls and frame-dragging diagnostics in
[the spin plan](./tiplercylinder/research.md).

The original “railgun” claim needs a physics correction: in a stationary,
localized spacetime, a freely falling particle cannot acquire arbitrary net
asymptotic energy because the timelike Killing energy is conserved. Coordinate
exit speed alone is gauge dependent. Net extraction requires a time-dependent
source, a Penrose-like split, or an exchange of energy/angular momentum with
the rotating matter.

**Decision criterion:** first recover the zero-spin null result and slow-spin
linear Sagnac scaling. Then report invariant timing, conserved energy, and
angular-momentum transfer. Treat a simple exit-speed increase as a diagnostic,
not as proof of acceleration.

### C2. Canonical-matter compact-object astrophysics

The canonical sector avoids the exotic-matter interpretation problem and can
reuse the search and waveform infrastructure.

Priorities:

1. self-gravitating boson-star/Q-ball encounter waveforms;
2. compact rotating configurations and ergoregion-instability boundaries;
3. directional GW emission with denser sources and extraction radii in the
   actual wave zone;
4. proper-time differences between specified worldlines, replacing
   gauge-dependent “minimum lapse” stasis claims.

Every waveform campaign should optimize a wave-zone-validated observable, not
raw local \(\Psi_4\) power.

### C3. Exotic compact-object waveform library

Build a library only from runs that pass propagation, \(1/r\), resolution, and
boundary checks. Store source parameters, ADM quantities, constraints, waveform
uncertainties, and detector-scaled strain alongside each mode.

The first useful entries are the static wormhole collapse controls and the
rotating collapse-on-command \(m=0\) burst (with its documented finite-radius
and convergence caveats). FTL candidates should enter only after their
radiation is separated from matter crossing the extraction sphere.

### C4. Lifetime and no-go atlas

A new synthesis direction is to treat failed confinement as data. Across
canonical and phantom families, map:

- \(Q_{\rm sphere}\) half-life and matter RMS growth;
- shortcut strength and duration;
- ADM mass, angular momentum, exotic fraction, and injected work;
- dominant failure mode: dispersal, horizon formation, constraint failure, or
  boundary contamination.

Use a Pareto archive rather than one weighted total. This can establish whether
the data indicate a smooth trade-off or a sharp stability boundary. A robust
no-go region for the tested scalar models would be a substantive result.

## Reconciliation of the original five ideas

| Original idea | Grounded status | Home |
|---|---|---|
| Sequential “pulse drive” | Testable only with all stages in constraint-solved initial data or with explicitly accounted external injection; continuous emission sweep already exists | Priority B4 and [trajectory results](./neuralspacetime/MapElitesDynamics.md) |
| Bicomplex pressure vessel ("sandwich") | Flagship pump-free confinement bet; high value but must be solved as one self-consistent bound state and pass a radial stability scan, not superposed profiles | Priority B2 |
| Spacetime railgun | Reframe as Sagnac/frame-dragging and energy-transfer study; net free-geodesic energy gain is constrained in a stationary spacetime | Priority C1 and [spin plan](./tiplercylinder/research.md) |
| GW signatures of warp collapse | Partially delivered: the rotating collapse-on-command \(m=0\) burst is drafted with three-resolution and big-box checks; precision claims still gated on wave-zone extraction and an AH finder | Priorities A1, A3, and C3 |
| RL tractor beam | Reframed by measurement: the PD pump does negligible net work and acts as a transient igniter, not support; RL must beat that igniter under the same accounting | Priority B3 and [RL plan](./RL/research.md) |

## Repository-wide falsification rules

These lessons have already been paid for and should be treated as project
requirements:

1. **Co-evolve matter.** Frozen prescribed sources are not valid dynamical
   solutions.
2. **Match initial-data and evolution resolution.** Do not evolve finer
   structure than the native initial-data solve without a demonstrated
   prolongation/convergence argument.
3. **Use complex winding for rotation.** A single real scalar with
   \(\cos(m\phi)\) modulation is a density pattern, not a smooth rotating state.
4. **Solve rotating momentum constraints.** Analytic \(K_{ij}=0\) rotating data
   have already produced a false instability.
5. **Use multiplicative validity gates.** Unbounded curvature or \(\Psi_4\)
   observables cannot be made safe with additive penalties.
6. **Promote before claiming.** Search-tier score is a hypothesis. HQ replay,
   longer duration, and a resolution/box control are evidence.
7. **Measure confinement with fixed physical regions.** Use
   \(Q_{\rm sphere}\), \(\rho_{\rm sphere}\), and matter RMS; `rho_sum` can remain
   nearly constant while matter leaves the source.
8. **Prove the wave zone.** Retarded-time alignment alone is insufficient;
   require \(1/r\) scaling, multiple radii, causal boundary exclusion, and
   uncertainty estimates.
9. **Separate coordinate diagnostics from observables.** Lapse, shift, local
   coordinate speed, and coordinate exit velocity are useful search features,
   not standalone physical claims.
10. **Account for external control.** Record pump work, charge/momentum
    injection, and constraint response. Never describe a driven solution as
    passive equilibrium.

No stable-transport, horizon, natural-mode, or detector claim should be made
below the relevant convergence and falsification tier.

## Recommended sequencing

Two tracks run in parallel and should not block each other:

- **Track 1 — harvest** (mostly writing and control runs; low GPU cost):
  finish the papers the existing data already supports.
- **Track 2 — confinement** (methods work and pilot solves; the physics bet):
  pursue the pump-free bound state via B1 and B2.

Each step below names a deliverable and a pass/fail check. A step without its
check completed does not count as done.

### Weeks 1--2

1. **Rotating manuscript finalization (A3, Track 1).** The runbook is done and
   the draft carries the hedges; this is now mostly writing plus one decision:
   whether a redesigned non-rotating control (different AMR/tagging/ID, per the
   lab record) is cheap enough to include, or whether the \(m=0\) power
   dominance stands alone as the rotation attribution.
   *Deliverable:* a submittable manuscript with the AH-finder and wave-zone
   gaps documented as limitations.
   *Check:* no claim in the abstract exceeds the early-cut / proxy-horizon /
   intermediate-zone evidence tier.
2. **Static wormhole controls (A1, Track 1).** Launch the cheap control pair
   already identified: pump-on/no-trigger run and shifted-trigger run, plus the
   throat-vs-trapped-surface disambiguation.
   *Deliverable:* control-run diagnostics folded into the manuscript.
   *Check:* the collapse time tracks the trigger, and the early \(\theta_+\le0\)
   surface is attributed (throat minimal surface vs genuine trapping).
3. **Sandwich feasibility pilot (B2, Track 2).** Before committing to the full
   two-sector methods work, attempt the smallest concentric
   canonical-plus-phantom constraint solve in the existing bicomplex model —
   static, spherical, one shell, one core.
   *Deliverable:* a solved (or demonstrably unsolvable) `.gridinit` with
   Ham/Mom below the postload gate.
   *Check:* if the coupled solve converges, proceed to the radial scan; if it
   diverges at every amplitude pairing, record the existence boundary and
   re-weight toward B1 — that is signal, not wasted time.
4. **Common campaign record.** Freeze one schema for every future run:
   \(Q_{\rm sphere}\) per sector, ADM mass and angular momentum, injected pump
   work, world-tube \(f_{\rm geo}(t_{\rm emit})\), failure mode, validation tier.
   *Check:* the next campaign writes it without hand-editing.

### Weeks 3--6

5. **Precursor paper consolidation (A2, Track 1).** Rebuild the matter-first
   manuscript around the verified 7.96% HQ result, with `eval 118` as the
   transient-dynamics case and the confinement decay reported as a primary
   measurement.
   *Check:* every number in the abstract traces to an HQ, gated artifact.
6. **Sandwich evolution + radial scan (B2).** If the pilot solve converged,
   evolve it unperturbed to \(t\ge30\) at matched resolution, then repeat with
   small inward/outward radial kicks.
   *Check (the go/no-go for the whole direction):* both sectors'
   \(Q_{\rm sphere}(t)\) flat within ~10% over several crossing times, pump off,
   no growing radial mode. Pass -> open the (shell mass, core amplitude) family.
   Fail -> write up the mechanism of failure and shift the confinement budget
   to B1.

### Months 2--3

7. **Self-gravitating rotating eigenstate (B1).** The flat-space eigenstate and
   the profile-4 painting pipeline exist; couple the solve to the metric by
   continuation from the ω=0 equilibrium. Run regardless of the sandwich
   outcome — the two attack different sectors (rotating single-field vs static
   two-field) and either can win. Add the normal-torus control when the
   non-exotic complex evolution branch lands.
   *Check:* scalar-EOM residual small at \(t=0\) and \(Q_{\rm sphere}\) flat to
   \(t=30\); a slower leak is a fail (the flat-space eigenstate's half-life
   \(\approx27\text{--}28\) is the baseline to beat).
8. **Phase-7 scale-matched collapse design (A3 follow-up).** Only after the
   rotating manuscript is out: bare-mass throat \(b_0\approx2\text{--}3\)
   hugged by a ring of compact bound Q-tori, tuned to the critical edge. The
   fast pre-test named in the lab record — re-solve the existing ID with
   \(b_0=2\) and cut support — decides whether matter binding or bare mass is
   the blocker.
   *Check:* a support cut produces `max_ah_r > 0` (real trapped surface), not
   dispersal to flat.
9. **Igniter characterization (B3).** The PD-pump work accounting already
   exists (`pump_work`, budget scripts); what remains is mapping which basins
   an igniter can reach across the (κ, R0, ω_orb) family.
   *Check:* a Pareto table (post-ignition lifetime vs injected work) exists
   before any RL training is proposed.
10. **Sagnac probe calibration (C1).** Zero-spin null test and slow-spin linear
    scaling first; physics runs only after both pass.

### Decision points (end of month 3)

- **Sandwich holds:** it becomes the program flagship — map its stability
  boundary, then put an FTL probe and a GW probe on it. This is the outcome
  that retires the pump everywhere.
- **Only the rotating eigenstate holds:** proceed to the (ω, m, κ) collapse
  grid, the [natural \(m=2\) campaign](./rotatingwormhole/PRDM2CampaignPlan.md),
  and the Phase-7 collapse redesign on that branch.
- **Only active control holds:** publish controlled spacetime engineering with
  the explicit igniter/energy budget from step 9 — the honest framing is
  "ignition into a longer-lived basin", since the pump has been shown not to be
  a continuous prop.
- **Nothing holds:** the lifetime/no-go atlas (C4) and the validated ECO
  waveform catalog (C3) become the primary output — a defensible negative
  result over the tested scalar families.

The central strategic shift is simple: finish the strongest existing evidence,
make confinement a measured research question with early cheap pilots, and let
failed validation change the claim rather than merely changing the score.
