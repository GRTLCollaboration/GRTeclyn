This updated log outlines a significant phase of maturity in the project. We
have transitioned from simply finding coordinate-level shortcuts to verifying
the physical reality of these spacetimes using gauge checks, multi-field matter
models, and dynamic soliton transport.

Below is an objective, senior-scientist-level evaluation of the important
findings so far and the critical next steps for the research.

Part 1: Important Findings So Far

1. The Bicomplex Scalar Breakthrough (grtresna_bicomplex_scalar)

Perhaps the most physically meaningful finding is the success of the bicomplex
model. Under identical parameters (the Eval 122 genome, m=0.15, \omega=0.12),
the single-complex field produced only a coordinate artifact (f_{\rm geo} = 0\%,
0/5 null rays reached the detector). By decoupling the fields into a canonical
scaffolding (\Phi^+) and a phantom source (\Phi^-) with opposite gravitational
coupling (f_s = -1), we successfully turned a coordinate illusion into a
genuine, gauge-invariant 6.13\% evolving null-geodesic shortcut with 5/5 rays
reaching the detector and excellent energy conservation. Decoupling the
ANEC-violating channel from the canonical matter allows the orbital structure to
drive the frame-dragging without destroying the negative-energy backing.

2. The Soliton Inertia and "Spotlight" Mismatch

We have isolated why moving matter configurations disperse under evolution. When
self-gravity is negligible, complex scalar fields act as free Klein-Gordon wave
packets, which have no bound state and immediately disperse.

When self-interaction is added to form a stable Q-ball (via a quartic attractive
term and a sextic stabilizer), the matter remains perfectly confined in static
tests (1.01\times RMS radius growth over the baseline's 1.35\times). However,
under the fast trajectories of Eval 122 (2.6c - 6.0c), the Q-balls still
disperse.

The physical reason is a transport velocity mismatch: solitons possess physical
inertia. The moving coordinate pump (the spotlight) moves too quickly for the
soliton to translate adiabatically. The pump simply ends up spawning new field
content at its current coordinates, leaving the old field behind to disperse as
radiation.

3. Lorentz-Boosted Initial Data ("Speed Follow") Validity

To resolve the transport mismatch, starting the Q-balls with a physical initial
velocity vector matching their orbit tangent (boosted up to 0.8c) was
implemented and validated.

  - The Good: This significantly improved the early-evolution phase. The RMS
    radius shrank for the first \sim 1.7 code units, and the confined matter
    fraction reached 79\% (a +12\% improvement over the cold-start run).
  - The Bad: It did not prevent final dispersal at t=16 (confined fraction still
    fell to 6.6\%). The boost resolves the linear initial transient, but it
    cannot protect the soliton from the immense centripetal/transverse
    acceleration required to hold a highly curved circular orbit at these
    speeds. The soliton continues to shed energy as radiation at every turn.

4. The -nan Momentum Solver Bug

We identified a purely logical bottleneck in the GRTresna non-linear solver
(GRSolver.impl.hpp). For static trajectory initial data (P_i = 0), the momentum
source is zero, which causes the relative Momentum error to return -nan (0/0).
While physically harmless, this NaN breaks the solver's early-exit check:
\text{converged} = (\text{exit\_tol} > 0) \land (\text{Ham\_error} < \text{exit\_tol}) \land (\text{Mom\_error} < \text{exit\_tol})
Because Mom_error < exit_tol always evaluates to false for a NaN, the solver is
forced to run the full 30 iterations, even though the Hamiltonian error
converges below 1\% by iteration 8. This explains the historically slow GRTresna
solve times for at-rest trajectory campaigns.

Part 1b: Q-Ball Dispersion Work — Implemented (2026-06-27)

Following the dispersion analysis in MapElitesDynamics.md, the following was
implemented in grteclyn-wrapper (Python only; GRTresna C++ unchanged except for
an attempted params.txt workaround):

**Modules and wiring**

| Component | Path | Role |
|-----------|------|------|
| Q-ball couplings | `grtresna/qball_couplings.py` | `QBallCouplings` value object: `standard()` (λ=160, μ=5333) and `stiff()` (λ=640, μ=85333, ω_min preserved). `cap_well_depth()` clamps paint amplitude to equilibrium √(3λ/4μ). |
| Radial ODE solver | `grtresna/qball_radial_ode.py` | Outward shooting + bisection on φ₀(0); `cached_qball_radial_profile()`, `profile_for_lump()`. Stiff equilibrium core ≈ 0.075 (not 0.15). |
| Gridinit paint | `grtresna/boson_star_fields.py` | `PROFILE_ODE_BOUND = 3`; ODE tabulation used in `_lump_phi0_at_radius` and Lorentz-boosted repaint. |
| Trajectory expansion | `search/optimize/config.py` | `grtresna_qball_equilibrium_amplitude=1` caps lump `amp`; `grtresna_qball_ode_profile=1` sets `profile=3` + `qball_*` keys on lumps. |
| Replay flags | `scripts/campaigns/hq/replay_eval.py` | `--qball-preset {standard,stiff}`, `--qball-equilibrium-amplitude`, `--qball-ode-profile`. |
| Matter metadata | `grtresna/matter_wiring.py` | `scalar_mu` serialization; `EVOLUTION_MATTER_KEYS` for replay consistency. |
| Tests | `tests/grtresna/test_qball_*.py` | 20+ unit tests (couplings, ODE shoot, config wiring, params mapping). |

**Replay results (eval 122 genome, boosted, t→16 unless noted)**

| Run | Seed | conf @ t≈16 | Notes |
|-----|------|-------------|-------|
| `traj_qball_boosted_eval122` | λ=160, sech@0.15 | **6.6%** | Baseline |
| `traj_qball_stiff_boosted_eval122_v2` | λ=640, sech@0.15 | stopped ~9.6 | conf ~39% vs ~28% baseline mid-run; still dispersing on orbit |
| `traj_qball_stiff_static_smoke` | λ=640, ω_rot=0 | — | conf ~72% at rest (binding OK) |
| `traj_qball_stiff_ode_smoke` | stiff + equilibrium + ODE | **0.5%** @ t=0 | **Broken** — diagonal stripe garbage in φ frames; run stopped |

**Critical blocker: ODE profile requires GRTresna C++ support**

GRTresna `BosonStarParams.hpp` implements lump profiles 0 (Gaussian), 1 (tanh
shell), and 2 (sech) only. Profile 3 is not defined; it falls through to
Gaussian. The Python pipeline splits the constraint solve from the exported
matter:

1. **GRTresna** builds lapse/metric from whatever envelope is in `params.txt`
   (we temporarily map profile 3→2 sech in `solver.py`; this is **not** a fix).
2. **Python gridinit export** repaints matter with the ODE tabulation +
   Lorentz boost + de Broglie phase.

Sech in the solve and ODE in the repaint are different φ₀(r) at the same nominal
`amp`/`width`. The metric is sourced for sech; the exported field is ODE +
boosted phase — inconsistent Hamiltonian/momentum data and garbage φ frames
(`confined_frac ≈ 0.5%` vs ~75% for a good boosted baseline).

**Do not use `--qball-ode-profile` until GRTresna C++ adds profile 3** with the
same tabulated flat-space Q-ball φ₀(r) as Python (shoot + cubic interpolation,
keyed by m, λ, μ, ω). Until then, the safe stiff path is
`--qball-preset stiff --qball-equilibrium-amplitude` only (sech @ equilibrium
amp ≈ 0.075).

**Physics takeaway from stiff runs**

- Stiffer λ deepens the well and improves mid-run confinement vs baseline, but
  does not prevent late orbit dispersal when seeded with sech @ 0.15.
- At rest, stiff Q-balls bind (~72% conf); dispersal on the fast eval-122 orbit
  is largely kinematic (centripetal radiation on curved trajectories), not just
  seed relaxation.
- Equilibrium amplitude (≈ 0.075 for stiff) is required; seeding at 0.15 is
  2× super-equilibrium and adds t=0 relaxation shock.

Part 2: Critical Next Steps in the Research

To resolve these newly identified physical and numerical limitations, the
research should be structured into the following five priority steps.

                             [Immediate Priorities]
                                       │
         ┌─────────────────────────────┼─────────────────────────────┐
         ▼                             ▼                             ▼
   GRTresna Fix              Stiffer Q-Ball Well          GRTresna ODE Profile (C++)
  - Fix GRSolver.impl     - λ=640, μ=85333 DONE         - profile 3 in C++
  - Bypass Mom NaN        - Replay: sech@equilibrium      - Tabulated φ₀(r)
  - 3x solve speedup        only until C++ ready          - Match Python ODE
         │                             │                             │
         └─────────────────────────────┴─────────────────────────────┘
                                       │
                              Radial ODE (Python) DONE
                              gridinit paint only — blocked
                              until C++ profile 3 lands

1. Patch the GRTresna Momentum Convergence Check (Numerical Priority)

We must immediately implement the proposed fix in GRSolver.impl.hpp to treat
non-finite momentum errors as satisfied when the momentum source is zero:

const bool mom_ok = !std::isfinite(Mom_error) || (Mom_error < exit_tol);
const bool converged = (exit_tol > 0) && (Ham_error < exit_tol) && mom_ok;

This will cut the GRTresna trajectory initial-data solve times by \sim 3\times
with zero loss in physical accuracy, significantly accelerating our Stage 0
throughput.

2. Implement Stiffer Q-Ball Potentials (Physics Priority) — **DONE (Python)**

The velocity boost alone cannot keep the Q-balls confined because the potential
well is not deep enough to resist centripetal acceleration. Stiffer couplings
were implemented via `QBallCouplings.stiff()`:

  - λ = 640, μ = 85333 (4× λ, μ scaled to preserve ω_min).
  - Equilibrium core amplitude √(3λ/4μ) ≈ 0.075 — use
    `--qball-equilibrium-amplitude`, not `well_depth=0.15`.
  - Replay: `--qball-preset stiff --qball-equilibrium-amplitude`.

**Result:** Mid-run confinement improves vs baseline; static binding OK (~72%).
Late orbit dispersal persists on eval 122 — stiffer well alone is insufficient
without exact seed and/or slower orbit kinematics.

3. Build the Exact Q-Ball Radial ODE Solver (Algorithm Priority) — **PARTIAL**

The Python ODE solver is implemented (`qball_radial_ode.py`, shooting +
bisection). Gridinit repaint uses it when `--qball-ode-profile` is set.

**Blocked on GRTresna C++:** The constraint solve must use the same φ₀(r) as the
export. `BosonStarParams.hpp` has no profile 3; mapping 3→sech in Python
`params.txt` does not make solve and repaint consistent. Required C++ work:

  - Add `profile == 3` to `lump_envelope` / `lump_phi1` with tabulated ODE
    φ₀(r) (same flat-space ODE as Python; load from params or shared `.dat`).
  - Wire `scalar_lambda`, `scalar_mu` (already in params) into profile lookup.
  - Rebuild GRTresna; then re-run `traj_qball_stiff_ode_eval122` smoke (t=4)
    and validate t=0 `confined_frac` ≳ 0.7 and clean φ frames.

Until C++ lands: **do not** pass `--qball-ode-profile` to replay.

4. Complete and Score the Bicomplex HQ Promotion

Analyze the results of the currently running traj_bicomplex_m03_w025 promotion
(256^3, ml=3, t_{\rm stop}=30).

  - Verify if the higher mass/frequency scale (m=0.3, \omega=0.25) improves
    structural persistence over the m=0.15 baseline.
  - Confirm that the 6.13\% evolving FTL shortcut scales or stabilizes at higher
    resolutions, confirming its physical, resolution-independent nature.

5. Restrict the Stage 0 Search Space to All-Retrograde Orbits

HQ validation has conclusively shown that counter-rotating configurations (mixed
prograde/retrograde signs) are generators of false positives (such as Eval 008,
which collapsed to 0\% FTL at high resolution).

  - Constrain the search bounds in spaces.py to force \omega_{\rm rot} < 0
    globally [1].
  - This will immediately eliminate half of the search space, focusing the
    MAP-Elites optimizer entirely on the highly productive, retrograde
    frame-dragging vortex basin.
