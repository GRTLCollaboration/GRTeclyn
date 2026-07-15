# Orbital-Pump Wormhole Support — Implementation Plan

Companion to [`RotatingWormholePlan.md`](./RotatingWormholePlan.md) (this is its
**Rung 2 — active support**, fused with the multi-lump constellation idea) and to
the trajectory-pump machinery journaled in
[`../neuralspacetime/MapElitesDynamics.md`](../neuralspacetime/MapElitesDynamics.md).

---

## ✅ FINAL RESULT (2026-07-14) — collapse-on-command produces a clean GW burst, no l=2 kick

**Headline.** A rotating (m=1, ω=0.25067, κ=1.0) Q-torus wormhole, held open by the
PD support pump, **collapses to an apparent horizon and radiates a gravitational-wave
burst when the support is ramped off — with *zero* imposed perturbation.** This is the
rotating analogue of the static Ellis–Bronnikov `S_support→0` collapse branch, and the
GW signal is a genuine consequence of the axisymmetric throat pinch, not of any seeded
mode.

**Run.** `runs/rotating_wormhole/evo_omega_p0p25_m1_kappa_1p00_dx0p5_ml3_mass0p5_qball_lam170_mu614450_torus_wh_collapse_ml3_gw`
(N=128, dx=0.5, `max_level=3`, `t→40`, dense in-code Weyl4 + Python Ψ₄ extraction).

**No kick.** `wormhole_phi_perturbation_amplitude = 0.0`, `wormhole_phi_monopole_amplitude = 0.0`.
The *only* trigger is the support ramp `support_ramp_t_start=8.0 → support_ramp_t_end=10.0`,
`support_ramp_floor=0.0` (pump turned off between t=8 and t=10). No l=2 kick, no Gaussian K
perturbation.

**Collapse — confirmed.**
- `min_lapse`: 1.000 → **0.033** (deep collapse well).
- `min_chi`: → **0.0185** (throat pinch).
- `max|K|`: → **0.62**.
- Apparent horizon / trapped surface **forms**: `max_ah_r = 9.18 (>0)`, `min_theta_plus = −1.90 (≤0)`.

**Constraints — clean, ran to completion, no NaN.**
- Ham L2 ≤ **3.25e-2**, Mom L2 ≤ **9.67e-3**, flat to t=40 (STEP 4000).

**GW radiation — confirmed (Python l=2,m=0 extraction).**
- Smooth coherent waveform: near-zone response ~t=6–8, **main burst ~t=13–24**, then decay —
  causally consistent with support-off at t=8–10 → throat collapse → horizon → radiation
  reaching r=12.
- Extraction is physically valid: `m=0` imaginary part RMS ~**1e-5** (≈0, as required).
- Physical scaling: peak strain frequency ≈ **169 Hz** at M=30 M☉ (Advanced-LIGO band),
  ≈ **5 Hz** at M=1000 M☉ (decihertz / DECIGO band).

**⚠️ Caveats / open items.**
- **J_z drift ~141%, Q_total drift ~116%** (both change sign) over the run.
  Decomposed 2026-07-15: ∫pump_work=0 (PD pump off); drift is **not** pump
  injection. Q_sphere retained (~5%). Do not quote global conservation.
- **The dense in-code C++ Weyl4 extraction (`data/Weyl4_mode_2{0,1,2}.dat`) is buggy**
  (spurious `O(1)` `m=0` imaginary part + ~1.8 M⁻¹ junk oscillation, uncorrelated with the
  trusted Python signal). Physics plots use `small_data/psi4_mode_l2m0.dat` (Python post-hoc
  spherical-harmonic extraction). C++-side fix tracked in the progress log below.
- **2026-07-15 runbook:** causality controls, ADM mass, m-spectrum, and coarse
  convergence revise several headline phrasings — see
  **§ Lab journal — Article Results Runbook** below (lab record; manuscript
  edits are separate).

**Reproduce the figures** (see `grteclyn-wrapper/README.md` §Visualization):
```bash
RUN=runs/rotating_wormhole/evo_omega_p0p25_..._torus_wh_collapse_ml3_gw/output
bash grteclyn-wrapper/scripts/plot/plot_diagnostic.sh "$RUN" 12 18 24   # constraints, collapse, GW
bash grteclyn-wrapper/scripts/plot/make_movies.sh   "$RUN" --framerate 10
```

**Goal.** Keep the rotating-wormhole throat supported by a constellation of
exotic Q-ball lumps driven along prescribed orbits by the PD "spotlight" matter
pump. Then trigger collapse **on command** by ramping the pump off — the
rotating analogue of the static Ellis–Bronnikov article's `S_support = 0.5`
collapse branch. MAP-Elites searches the orbit/pump space for maximally stable,
minimal-intervention configurations.

**Why (one paragraph of context for the executor).** All passive rotating IDs
disperse: a spherical Q-ball twisted into an `m=1` torus is not a stationary
eigenstate, so `Q_sphere` (the throat-sphere Noether charge — the ONLY valid
confinement metric; `rho_sum` is misleading) drains with half-life t≈13–16 and
the spacetime relaxes to flat instead of collapsing (RotatingWormholePlan
"Rung 1a′ DONE"). The pump replaces the missing stationarity with active
control. This is an *engineered* (actively supported) wormhole, not a solution
of the pure Einstein–Klein–Gordon system — constraint drift from the pump's
Bianchi violation must be monitored and penalized throughout.

---

## Verified file map (all paths relative to repo root `GRTeclyn/`)

| What | Where | Status |
|---|---|---|
| Wormhole evolution example | `Examples/RotatingWormholeCollapse/` (`Main_SupportedWormhole.cpp`, `SupportedWormholeLevel.cpp`, `SimulationParameters.hpp`) | builds; `complex_scalar` branch at `SupportedWormholeLevel.cpp` ~line 241 uses `ComplexExoticScalarField<ComplexScalarPotential>` |
| Matter class used by wormhole | `Source/Matter/ComplexExoticScalarField.{hpp,impl.hpp}` | **NO pump hook yet** |
| Matter class WITH pump (reference) | `Source/Matter/ComplexScalarField.{hpp,impl.hpp}` | has `RLMatterPumpParams m_pump`; PD target-soliton drive at `impl.hpp` lines ~146–215 (`k_p>0` closed-loop; `k_p<=0` legacy open-loop) |
| Trajectory machinery | `Source/GRTeclynCore/RL/{TrajectoryParams.hpp, TrajectoryEvaluator.hpp, RLMatterPumpParams.hpp, RLLumpState.hpp}` | complete, example-agnostic |
| Wiring reference (copy from here) | `Examples/RadialRecipe/`: params parsing `SimulationParameters.hpp` ~lines 94–123 (`trajectory_mode`, `trajectory_lump<k>_*`, shared `trajectory_*`); per-step update `Main_RadialRecipe.cpp` ~lines 38–110 (seed at t=0, update lump centres each coarse step); pump construction `RadialRecipeMatterDispatch.hpp` `build_rl_pump()` ~lines 70–130 | working, validated in QD campaigns |
| GRTresna multi-lump ID | `../GRTresna/Source/Matter/BosonStarParams.hpp`: `num_lumps` + per-lump `lump<k>_{amp,width,center,velocity,omega,mode,exotic,winding,profile,profile_path}` (parse at ~lines 499–540) | multi-lump, per-lump velocity + exotic flag already supported |
| ID driver | `grteclyn-wrapper/scripts/wormhole/id/solve_kappa_family.{sh,py}` | single-lump today; extend for constellation |
| Run CLI | `grteclyn-wrapper/scripts/wormhole/run/wormhole_case.{sh,py}` | renders params from template + flags, launches with plotfile-deleting sidecar |
| Diagnostics | `output/data/collapse_diagnostics.dat`, 22 cols; **`Q_sphere` = col 21** (confinement), `rho_sum` = col 18 (DO NOT use for confinement), `min_chi`, `min_theta_plus`, `max_ah_r`, `J_z`, `Q_total` | AMR-correct (level-0 synchronised integrals) |
| QD campaign pipeline | `grteclyn-wrapper/scripts/campaigns/{qd/run.sh, lib/, hq/}` | validated 8-GPU pipeline, 200-eval runs, HQ promotion |

Baseline numbers to beat (passive, κ=1, m=1, N=128/dx=0.5, ml=1, μ=0.5, t→30):
`Q_sphere` retention 0.99/0.92/0.56/0.27/**0.05** at t=6/10/15/20/30
(half-life ≈ 13–16). Reproduce commands are in RotatingWormholePlan §"Rung 1a′ DONE".

---

## Lab journal — Article Results Runbook (2026-07-15)

Campaign to harden the 2026-07-14 headline claim before treating the manuscript as
submission-ready. Driven by referee-style blockers in [`NextSteps.md`](./NextSteps.md)
(Tier 1 + Tier 2). This section is the **lab record** (runs, numbers, failures) —
not a draft of the paper.

Headline ID (all arms unless noted):
`runs/rotating_torus_id/torus_m1_om0p250_kappa1p00_dx0p5_L64_lam170_mu614450_exotic_throat1/initial_data.gridinit`
(m=1, ω≈0.25067, κ=1, exotic, throat bare mass 1, L=64 unless noted).
Launch always via `wormhole_case.py` with `--no-frames` + plotfile delete sidecar
(`--delete --keep-last 3`); disk still ballooned when consumers lagged behind AMR
plotfiles (~1–2 GB each) — manually trimmed live backlogs; crashed runs pruned to
~13–18 MB keeping diagnostics only.

### Tooling added (postrun)

| Script | Purpose |
|--------|---------|
| `grteclyn-wrapper/scripts/wormhole/postrun/adm_quantities.py` | ADM surface integrals on gridinit |
| `…/charge_energy_budget.py` | ΔJ_z / ΔQ / ∫pump_work / E_GW proxy |
| `…/psi4_extrapolate.py` | 1/r fit of peak \|rΨ₄\| + relative constraint norms |
| `consume_plotfiles` | now also writes `psi4_mode_l2_all.dat` (m∈{−2..2} per radius) |

### Done — Tier 1

**A1 causality controls — DONE (NaN, not t→40).**
| Run suffix | Ramp | Outcome |
|------------|------|---------|
| `noramp_ctrl` | none (`t_start=-1`) | NaN in `K` @ **t≈18.82**; min_a never &lt;0.15; max_ah→15.6; θ₊→−O(100) |
| `ramp_t16` | [16,18]→0 | NaN in `h11` @ **t≈18.85**; same shallow-lapse / runaway-proxy pattern |
| headline (prior) | [8,10]→0 | clean t→40, min_a→9.8e−4 |

Bit-identical to headline until the early ramp would have acted. Early support cut
**selects the resolvable deep-collapse branch**; no-ramp / late-ramp do **not** hold
a quiet throat — they abort without deep lapse. “Collapse anytime on command” is
**not** demonstrated; “early support cut → deep collapse + GW” is.

**A2 t≈2.3 θ₊ — DONE.** Identical crossing in noramp and headline at t=2.26 while
still in the gauge transient → throat minimal-surface proxy, not collapse. Robust
collapse markers after support cut (t≳12).

**A3 modes + non-rotating — MOSTLY DONE.**
- Python directional / all-m: burst power **≳99.999% in m=0** (`conv_dx067`
  `psi4_mode_l2_all.dat`; headline directional beam_ratio ~1e−5).
- Non-rotating ID: same torus \|Φ\| table, m=0, ω=0, exotic+throat1
  (`torus_m0_om0p000_*_exotic_throat1`, Ham~0.51%, J_ADM=0, M_ADM≈−0.24).
- `norot_m0_ramp`: **Arena OOM** @ t≈9.6 (GPU memory), not a physics verdict.
  Attribution → toroidal geometry, not rotation; full norot waveform still open.

**A4 ADM — DONE.** On exotic_throat1 gridinit, R∈{16,20,24,28}:
**M_ADM = 0.209 ± 0.001**, **J_ADM_z ≈ −1.93** (volume J_z(t=0)≈−1.47, same sign).
Prior Hz scalings assumed code mass=1; correct physical f scales with M_ADM
(≈×0.21): ~35 Hz @ 30 M☉ / ~1 Hz @ 10³ M☉ if one identifies M_ADM with the
astrophysical mass.

**A5 / A6 text+budget — DONE (analysis).**
- ∫pump_work dt = **0** (PD pump off in headline; only support-strength ramp).
- Q_sphere retained (~5% drift); Q_total +116%, J_z +141% (sign change) — **do not
  claim conservation**; drift is not pump injection.
- E_GW (l=2,m=0 NP proxy) ≈ 0.025 ≈ 0.12 M_ADM (order-of-mag, near-zone).
- Propagation: v = 1.00 ± 0.08 c (sampling-limited), not “exactly 1.00”.

### Done — Tier 2 (partial)

**B1 convergence — ONE RUNG DONE.**
`conv_dx067` (dx≈0.667, N=96, ml=3) finished **cleanly t→40**:
- min_a → 0.0020, recovers to 0.052; min_chi → 0.033; max_ah 9.11 @ t≈15.5 → 1.97
- Ψ₄ peak \|Re\|(R=12) ≈ 6.1e−2 (headline 7.7e−2); Ham max 3.1e−2
- Phantom-bounce-like (ah shrink + lapse recover) **present at coarse dx**
- Two-point vs dx=0.5: waveform overlap ~0.98, ~20% amplitude diff; constraint
  max **does not** scale as h^p (obs. p≈0) → **Richardson order not yet shown**
  (need finished dx=0.4 / 0.333)

**B2 extrapolation — analysis DONE; large-box RUN.**
Linear lim \|rΨ₄\|_∞ ≈ 1.29 from {12,18,24}; quadratic unstable (near-zone).
`bigbox_L128` (L=128, R_ext=50/60, sponge 54–62, stop 70) still evolving —
heavy (~80 GB GPU).

**B3 relative constraints — DONE enough.** Ham L2 max ≈ 10× t=0 residual;
Mom similarly elevated. No Θ/Z_i or L∞ columns in `constraint_norms.dat`.

**B4 AH finder — GAP.** No production AH finder in GRTeclyn; keep coordinate
θ₊ proxy; documented.

**B5 energy budget — DONE** at proxy level (see A6).

### Still going (live GPUs, 2026-07-15 ~morning)

| Suffix | GPU load (approx) | Progress | Notes |
|--------|-------------------|----------|-------|
| `rerun_allm` | ~53 GB | t~23/40 | headline physics + all-m file |
| `conv_dx040` | ~77 GB | t~17/40 | fine ladder rung |
| `conv_dx033` | ~49 GB | t~15/40 | finest rung (may grow mem) |
| `bigbox_L128` | ~81 GB | t~19/70 | wave-zone extract; near OOM |

Idle GPUs = finished/crashed arms (noramp, ramp_t16, norot, conv_dx067).

### Results summary (lab, not paper prose)

1. **Deep collapse + GW + late bounce** reproduce at dx=0.5 and dx≈0.67.
2. **Burst is axisymmetric (m=0)** — geometry, not spin modes.
3. **Early support ramp** is what yields the clean deep-collapse branch; no-ramp
   / late-ramp **NaN ~t=19** without min_a→10⁻³.
4. **M_ADM≈0.209** must anchor any physical-unit claim; old 169/5 Hz assumed M=1.
5. **No global Q/J conservation**; pump_work=0 so prior “pump-injected” guess for
   the drift was wrong for this headline config.
6. **Richardson order: not yet** — waiting on finer finished rungs.
7. **Non-rotating control: OOM** — re-run at lower ml or memory still TODO.

### Open follow-ups

- Finish `conv_dx040`, `conv_dx033`, `bigbox_L128`, `rerun_allm`; then 3-point
  Richardson on Ψ₄ / min_a / constraints + wave-zone cross-check.
- Re-launch norot at `max_level≤2` or smaller box to complete A3 waveform compare.
- Optional: tighter plot_interval / keep-last=1 on heavy arms to protect disk.

---

## Progress log

### Phases 1–3 DONE (2026-07-13) — pump wired, constellation ID solves, passive constellation is unstable

**Phase 1 (pump hook + CLI + diagnostic) — COMPLETE, GPU-validated.**
- `ComplexExoticScalarField` gained the `RLMatterPumpParams m_pump` hook + the
  7-arg `add_matter_rhs(...,coords,time)` overload (PD trap `k_p>0` + legacy
  open-loop), copied verbatim from `ComplexScalarField` (no extra phantom-sign
  flip — the exotic class flips only the stress tensor, not the field EOM).
- `CCZ4RHSWithMatter` now routes to the time/coords overload via a `void_t`
  detection trait (`detail_matter::has_time_rhs_v`) instead of the hard-coded
  two-class `if constexpr`, so any matter class with the overload gets the pump.
- Wormhole `SimulationParameters` parses `trajectory_*` + `rl_pump_*` +
  `pump_ramp_*`; `Main_SupportedWormhole` seeds/updates the `TrajectoryEvaluator`
  each coarse step and applies the ramp; `SupportedWormholeLevel` builds the pump
  (`build_wormhole_pump`) in the `complex_scalar` RHS branch and publishes
  `L2_Ham` for the governor.
- **`pump_work` diagnostic** added (col 22): level-0 integral of
  `alpha (Pi1 S1 + Pi2 S2) sqrt(gamma) dV` (the pump's Bianchi-violating power);
  0 when the pump is off.
- **CLI (no repo params files):** `wormhole_case.py` gained
  `--num-lumps --orbit-radius --orbit-omega --well-depth --well-width --pump-kp
  --pump-kd --pump-frequency --pump-target-amp --pump-ramp-*` and
  `--id-num-lumps/--id-orbit-radius/--id-orbit-omega` (ID selection). All params
  are generated into the run dir; `--num-lumps 0` emits zero trajectory keys
  (clean no-op).
- **L6 fix:** `--no-frames` now still runs a delete-only sidecar, and plotfile
  pruning runs even when the binary NaN/aborts (was leaking tens of GB on both
  `--no-frames` and crashed runs).
- GPU build (`USE_CUDA=TRUE USE_MPI=TRUE CUDA_ARCH=90`, local OpenMPI 5.0.8 for
  gcc-11/nvcc compat) links clean. Pump-on vs pump-off smoke differ measurably
  (min_chi 0.8988 vs 0.8964, Q_sphere retention slightly better with pump),
  constraints bounded, `pump_work` populates.

**Phase 2 (constellation ID) — COMPLETE.** `solve_kappa_family.py` gained
`NUM_LUMPS/ORBIT_RADIUS/ORBIT_OMEGA`: N Q-ball lumps on a circle, tangential
boost velocity `v=omega_orb*R0` (coherent orbital J_z), `winding=0` (orbital AM
replaces phase winding), `exotic=1`, per-lump internal U(1) frequency
`LUMP_OMEGA` → Noether charge via `pi2=-(omega/alpha)phi1`. A 2-lump solve
(R0=6, omega_orb=0.1, mass=0.5, lam=170, mu6=14450, N=128) **converged Ham
0.62% / Mom 0.89%** (both < 1%). Tags match between solver and run CLI
(`..._nlump2_R6_worb0p1`).

**Phase 3 (passive constellation) — DONE, negative as expected.** 2-lump
constellation, pump OFF, t→30, at **both** max_level=1 and unigrid:

| t | Q_sphere | J_z | min_chi | outcome |
|---|----------|-----|---------|---------|
| 0 | −2.19 | 26.1 | 1.00 | — |
| 8 | −2.20 | 9.0 | 0.991 | charge fully held |
| 12 | −2.01 (~92%) | 9.6 | 0.987 | smooth decay |
| 13.5 | **NaN blowup (h11)** | — | — | dynamically unstable |

Q_sphere is retained ~92% to t=12 (comparable-to-better than the single-torus
baseline's t=10 value), then a **violent dynamical blowup at t≈13.5** — NOT the
single-torus's smooth dispersal. Reproduces at ml=0 and ml=1 ⇒ **physical, not
an AMR artifact**: passive phantom lumps do not orbit stably (self-repulsion /
non-stationary orbit). This is the quantitative green light for the pump
(Phase 4): the question is whether active support beats the t≈13 instability.

> **Column-map correction (found in Phase 4).** `collapse_diagnostics.dat`
> columns are (1-based, from `SupportedWormholeLevel.cpp` header):
> 1 time, 2 min_lapse, **3 min_chi**, 4 max_abs_K, 8 max_ah_r,
> **9 min_theta_plus**, 18 rho_sum, **19 J_z**, 20 Q_total, **21 Q_sphere**,
> 22 rho_sphere, **23 pump_work**. The Phase 3 table above was read with an
> off-by-one/older map (its "min_chi 1.00" was min_lapse; "Q_sphere −2.19" was
> Q_total; "J_z 26.1" was rho_sum). Correct t≈0 values: min_chi 0.513 (the ID
> already has a deep throat), Q_sphere −1.919, Q_total −2.20, J_z −0.73.

### Phase 4 DONE (2026-07-13) — pump beats the blow-up, but the native constellation collapses to a BH

k_p sweep at the ID-matched orbit (2-lump, R0=6, ω_orb=0.1, ω_U(1)=0.05,
target_amp=0.1 Gaussian σ=5, k_d=0.01, governor 0.03/0.005, ml=0, t→30):

| k_p | outcome |
|-----|---------|
| 0.05, 0.1, 0.5 | numerical **blow-up at t≈13** (Q,K→1e30), identical to passive |
| **1.0** | **bounded to t=30** — survives the t≈13.5 barrier |

- **pump_work is negligible** (~1e-8 at k_p=0.1, ~1e-6 at k_p=1.0) at every k_p:
  the PD trap does almost no *net* energy work because the field tracks the
  target closely (error→0); its value is the *restoring force*, not energy
  injection. L2_Ham stays bounded (≤0.028), so blow-ups are **not** constraint
  runaway — they are the code failing to resolve an unsupported collapse.
- Weak pumps (k_p≤0.5) are too feeble to change the passive outcome.
- **k_p=1.0 changes the physics qualitatively:** min_chi rises to 0.92 (throat
  opening) through t≈5–8, then a **genuine apparent horizon forms at t≈8.8**
  (min_theta_plus crosses +0.026 → −0.046, deepening to −5.5; max_ah_r≈13→9),
  i.e. the coherent rotating exotic constellation **collapses into a spinning
  black hole**, and the moving-puncture-like evolution follows it stably
  (min_chi pinned ≈0.03–0.05) out to t=30. Q_sphere oscillates about 0
  (not retained).

**Decision:** the active-pump mechanism *works* — strong pumping converts the
passive numerical blow-up into a resolved, bounded evolution (control
achieved). But at its native amplitude/orbit this constellation is
gravitationally **supercritical**: it collapses to a BH rather than holding an
open throat, so the strict Phase-4 PASS (open throat, no horizon, Q retention
≥0.5) is **not** met.

### Phase 4b — "can more pump support prevent collapse?" (NO) — same ID, pump-only sweep

Horizon-formation time (min_theta_plus first < 0) as a function of pump knobs,
same 2-lump gridinit, ml=0:

| config | horizon at t |
|--------|--------------|
| k_p=1.0, target_amp=0.1 (baseline) | 8.8 |
| k_p=2.0, amp=0.1 | 8.57 |
| k_p=4.0, amp=0.1 | 8.51 |
| k_p=1.0, amp=0.15 | 7.61 |
| k_p=1.0, amp=0.2 | 6.94 |

- **Tighter PD tracking (k_p 1→4) does NOT delay collapse** (8.8→8.5, slightly
  *earlier*): the field already sits on the target (pump_work~1e-6), so raising
  the gain just tracks a supercritical target more faithfully. The pump is a
  *shaping* force, not gravitational repulsion.
- **More matter (target_amp 0.1→0.2) collapses ~2 t EARLIER** (8.8→6.9): the
  rotating cloud is net-**attractive** (orbital-boost KE + bare-mass μ=0.25
  dominate the phantom's negative potential energy), so adding "support" mass
  accelerates the BH. The config is gravitationally supercritical and no pump
  setting avoids it.

**Conclusion:** collapse resistance is not a pump knob. To find a sub-critical
open throat try: mass reduction (kappa<1), larger R0 / more lumps (spread the
mass), faster ω_orb (centrifugal support). Cheap pre-check done in Phase 4c.

### Phase 4c DONE (2026-07-13) — sub-critical open throat FOUND (mass reduction wins)

Four hand-tuned IDs re-solved (Ham/Mom all <1%), each evolved with k_p=1.0:

| ID | knob | v=ωR0 | outcome |
|----|------|-------|---------|
| A: R0=6, ω_orb=0.13, κ=1.0 | more spin | 0.78 | **horizon @ t=8.29** (earlier than baseline — extra orbital KE dominates) |
| C: R0=10, ω_orb=0.07, κ=1.0 | spread out | 0.70 | **horizon @ t=7.29** (spreading removed throat support → bare-mass throat collapses) |
| B: R0=6, ω_orb=0.1, **κ=0.6** (amp 0.06) | less mass | 0.60 | **STABLE open throat** — no horizon, min_chi plateau ≈0.86 held t≈8→30 |
| D: R0=8, ω_orb=0.1, **κ=0.7** (amp 0.07) | spread + less mass | 0.80 | survives (no horizon) but min_chi drifts 0.43→0.97 = slow **dispersal** (throat opening flat) |

- **Spin and spread both make it WORSE** (A, C collapse earlier): more orbital
  KE promotes collapse; moving lumps off-centre starves the throat of exotic
  support so the bare mass (μ=0.25) collapses on its own.
- **Mass reduction is the lever.** B (κ=0.6) is gravitationally sub-critical:
  the pump holds a **genuine persistent rotating throat** — min_chi rises from
  0.44 to a flat **≈0.86 plateau** and stays there (0.868→0.864→0.856) from
  t≈8 through the passive blow-up time (13.5) and beyond, with
  min_theta_plus=+0.025 (no trapped surface) the entire time. Q_sphere/J_z
  oscillate about 0 (orbital modulation through the fixed r<10 sphere), not
  decay. **This is the Phase-4 PASS: an actively-supported stable rotating
  wormhole throat.**
- D shows the *other* edge: slightly too little central support → the throat
  slowly disperses toward flat (min_chi→1). So there is a Goldilocks band in
  (amplitude, R0): B sits in it, D just above, A/C below.

**Result: a stable actively-supported rotating throat exists.** The critical
control axis is per-lump amplitude / total mass (κ); orbit radius and spin
trade off support-vs-collapse. Winning config: **B = 2-lump, R0=6, ω_orb=0.1,
amp=0.06 (κ=0.6), k_p=1.0** (min_chi≈0.86, no horizon, bounded to t=30). See
"Do we still need MAP-Elites?" below for why this makes the campaign optional.

**Surviving artifacts (evolution run dirs were pruned 2026-07-13; the ID
gridinits were kept).** B's initial data is at
`runs/rotating_wormhole_id/rotwh_omega_p0p05_m1_kappa_0p60_dx0p5_mass0p5_qball_lam170_mu614450_nlump2_R6_worb0p1/initial_data.gridinit`.
Reproduce B's evolution (auto-loads that gridinit, streams+prunes plotfiles):

```
grteclyn-wrapper/.venv/bin/python \
  grteclyn-wrapper/scripts/wormhole/run/wormhole_case.py \
  --kappa 0.6 --omega 0.05 --m 1 --dx 0.5 --mass 0.5 --lambda 170 --mu6 14450 \
  --max-level 0 --stop-time 30 --num-lumps 2 --orbit-radius 6 --orbit-omega 0.1 \
  --well-width 5.0 --well-depth 0.1 --pump-kp 1.0 --pump-kd 0.01 \
  --pump-frequency 0.05 --pump-target-profile 0 --pump-target-amp 0.06 \
  --governor-center 0.03 --governor-width 0.005 --no-frames --plot-interval 200 \
  --run-suffix pk1_B_lowmassk06 --gpu 0
```

Other surviving gridinits: κ=0.7/R8 (D, disperses), κ=1.0 nlump2 R6/R10/R6w0.13
(baseline, C, A — all collapse), κ=1.0 single-throat, and the ω=0.3/0.4 winding
IDs. Re-solve any missing one with `NUM_LUMPS=2 ORBIT_RADIUS=.. ORBIT_OMEGA=..
MASS=0.5 LAMBDA=170 MU6=14450 RES_N=128 solve_kappa_family.sh <kappa> 2`.

**NOTE (CLI change this phase):** `wormhole_case.py` gained `--run-suffix` so
parallel pump/κ sweeps sharing one ID land in distinct run dirs (case_tag does
not encode pump gains). Used by all Phase 4/4b/4c runs.

### Phase 6 attempt (2026-07-13) — pump-ramp is NOT a collapse trigger; the pump is a transient IGNITER

Ran the pump-ramp on config B (κ=0.6). The ramp scales `rl_pump_amplitude`;
in PD mode (`k_p>0`) amplitude only *gates* the pump, so a floor of 0 fully
cuts it while any floor>0 leaves the PD drive unchanged (graded floors would
need the ramp to scale `k_p` — deferred). Floor=0 cut at varying times, ml=0,
t→30. Outcome axis: does the throat pinch to a BH (min_theta_plus<0) or hold?

| run | pump | min_chi outcome | verdict |
|-----|------|-----------------|---------|
| control | on, never cut | 0.86 plateau (slow drift up) | stable wormhole |
| cut@12–14 (floor0) | cut late | 0.86, unchanged vs control | **no collapse** — self-sustaining |
| cut@8–10 | cut | 0.86 held to t≈16+ | **no collapse** |
| cut@4–6 | cut | opened to 0.86 then re-pinches (0.86→0.63→…) | **delayed collapse** |
| passive (cut@0) | never on | peaks 0.75, re-pinches | **collapse, horizon @ t≈10** |

**Key findings.**
- **The PD pump does negligible net work (~1e-6)** and is *not* a continuous
  gravitational prop. Cutting it from the *established* stable state (t≳8–10)
  does essentially nothing — the wormhole has become **self-supporting**.
- **The pump is a transient "igniter":** it is essential only for roughly the
  first t≈7–9 to shepherd the sub-critical constellation into the stable
  open-throat basin. Remove it before that window → the throat re-pinches and
  collapses to a BH (passive → horizon @ t≈10; cut@4–6 → delayed collapse).
  Remove it after → the throat persists autonomously.
- The pump's other visible effect is churning the U(1) charge: with the pump
  on, Q_sphere oscillates at the drive frequency; cut it and Q_sphere settles
  to a conserved value (passive Q_sphere is flat, as Noether requires).
- **So "collapse on command by ramping the pump" only works in the narrow
  ignition window and produces a *delayed* collapse — it is not a clean,
  arbitrary-time trigger for an established wormhole.**

**Consequence / recommended real trigger.** To collapse an *established* stable
rotating wormhole on command at an arbitrary time, ramp the quantity that
actually gravitates: **`wormhole_support_strength`** (it scales the exotic
stress-energy `ρ, S_ij` directly — the rotating analogue of the static
Ellis–Bronnikov `S_support` cut), or equivalently ramp the field mass/amplitude.
This needs a small code change (make `support_strength` time-dependent with the
existing ramp schedule, applied in the matter class) + a GPU rebuild. That is
the clean path to the collapse + Ψ₄ GW-burst result. The pump-ignition +
self-sustaining-wormhole finding is itself a strong result (an engineered
rotating wormhole that becomes autonomous after transient support).

### Phase 6 result (2026-07-13) — support-strength cut → DISPERSAL, not collapse (config B is sub-critical)

Ran the `wormhole_support_strength` ramp on config B at the larger box: L=128,
N=256 (dx=0.5), sponge ON (inner 28 / outer 62), extraction r=12,24, k_p=1.0
igniter, to t=50. Two arms, **bit-identical until t=12.02** (confirms the ramp is
the only difference):
- **collapse:** `support_ramp t=12→14, floor=0` (full exotic-support cut).
- **control:** support held at 1.0.

Diagnostics (`output/data/collapse_diagnostics.dat`, cols: min_chi=3, max_ah_r=8,
min_theta_plus=9, Q_sphere=21):

| quantity | collapse (cut) | control (held) |
|----------|----------------|----------------|
| horizon `max_ah_r` | **0 for the whole run** | 0 |
| `min_theta_plus` | stays +0.0127 (never <0) | stays +0.0127 |
| `min_chi` t=0 (ID throat) | 0.45 | 0.45 |
| `min_chi` post-cut dip | 0.854 @ t≈19.5 | 0.877 @ t≈19.8 |
| `min_chi` late | → 0.99 by t≈39 | → 0.99 by t≈28 |

**Findings.**
- **No black hole in either arm** (`max_ah_r=0`, `min_theta_plus>0` throughout).
  The `max_level=0` NaN worry never triggered — because nothing collapsed.
- **The throat is DEEPEST at t=0 (min_chi=0.45) and only ever *opens*.** Even the
  *control* (support held) drifts to flat (min_chi→0.99 by t≈28); the "0.86–0.88
  plateau" is a temporary pause in an overall opening trend, not a stable state.
- **Cutting the exotic support produced only a modest transient** — a slightly
  deeper re-pinch (0.854 vs 0.877) around t≈15–20 — then it, too, relaxed to flat.
  So the sub-critical throat **disperses on command; it does not collapse.**
- **GW (Ψ₄ ℓ=2 m=0) is weak and nearly identical** between arms (peak |Ψ₄| at
  R=24 ≈ 3e-3, both at t≈25, from the orbiting lumps + junk), **no collapse burst**.

**Interpretation.** Config B is light and exotic-dominated: removing the
NEC-violating support gives the deformation nothing to hold and it *heals toward
flat* (the "D"/dispersal outcome), rather than leaving a positive-mass remnant to
collapse. A genuine collapse-on-command + Ψ₄ burst needs a configuration poised
**near the critical (collapse) boundary** — i.e. more compact/heavier support
(the supercritical A/C arms collapse on their own; B is well on the sub-critical
side). This is a direct argument for the **Phase 7 scale-matched redesign**
(bigger throat + compact hugging lump ring), tuned near the stability edge so a
support cut tips it into collapse instead of dispersal.

**Artifacts:** `runs/rotating_wormhole/evo_..._collapse_supcut_t12_L128/` and
`..._control_nosupcut_L128/` (both ran to t=50, no NaN). Embedding-diagram script
`grteclyn-wrapper/src/grteclyn_wrapper/visualisation/wormhole_embedding.py`.

### Matter analysis of the L=128 t=50 runs (2026-07-13) — it DISPERSES; no stable-to-t=50 case

Pulled the matter diagnostics (`confinement.dat` cols confined_frac/rms_radius;
`collapse_diagnostics.dat` Q_sphere=21, rho_sphere=22). **Both arms behave almost
identically** (the support cut barely matters):

| t | confined_frac (r<10) | matter rms radius | Q_sphere |
|---|----------------------|-------------------|----------|
| 0  | 0.75  | 8.6  | −0.84 |
| 10 | 0.53  | 11.1 | ~0 (osc) |
| 20 | 0.23  | 16.3 | ~0 |
| 30 | 0.09  | 22.0 | ~0 |
| 40 | 0.02  | 28.6 | ~0 |
| 50 | **0.004** | **35.7** | ~0 |

- **Unambiguous dispersal:** ~99.5% of the matter has left the throat region by
  t=50; the cloud's rms radius grows ~linearly (8.6→36, ≈ the lump speed) — the
  lumps + debris fly outward. `Q_sphere` drains from −0.84 → ~0.
- **NO stable-to-t=50 case.** The "min_chi ≈ 0.86 plateau" holds only t≈8→18;
  by t≈28 even the **control** (support held) has `min_chi → 0.99` (throat opened
  to ~flat) and stays there. The earlier "config B is stable" claim was an
  artifact of stopping at t≈30 on the smaller L=64 box. **Run to t=50 at L=128,
  config B disperses with or without support.**

**Why it disperses instead of collapsing (the crux).** Two things block collapse:
1. **The support matter is not bound.** A Q-ball twisted onto an orbit is not a
   stationary eigenstate in the curved+rotating background, so the exotic cloud
   spreads and leaves regardless of support — nothing remains at the throat to
   collapse. (This is the "all passive rotating IDs disperse" problem from the
   intro; k_p=1.0 does not actually confine it — confined_frac still → 0.4%.)
2. **There is almost no positive gravitating mass.** A throat *stays open* via
   NEC-violating (negative-energy) matter and *collapses* only if enough
   *positive* mass-energy is concentrated. Config B is light and exotic-dominated:
   remove the exotic and the deformation just **heals to flat** — no positive-mass
   remnant to implode. (Contrast: the static `SupportedWormholeCollapse` DID
   collapse on an `S_support` cut because it had a substantial **bare-mass throat**
   held open by exotic matter. That mass balance is what the rotating case lacks.)

### Collapse-trigger recipe (what Phase 7 must deliver)

To get collapse-on-command + Ψ₄ burst, the ID must satisfy "enough positive bare
mass + just enough exotic support to hold it open; then cut the support":
1. **Give the throat real bare mass — enlarge `b0`** (0.5 → ~2–3;
   `wormhole_throat_radius` in the GRTresna base params). Same change as the
   Phase-7 "bigger throat", but the *reason* is collapse: bigger throat = more
   positive bare mass = something that implodes when unsupported. Reuses the
   proven static-case mechanism.
2. **Make the support matter bound so it does not disperse first** (compact
   Q-balls in a deeper well, or a central bound exotic core + orbiting ring). If
   the matter leaves before t_hold there is nothing to trigger.
3. **Tune to the critical edge:** ID stable *with* full support, supercritical
   *without* it. A 1-D sweep over `b0` (or the support floor) brackets the
   boundary. Diagnose success by `max_ah_r > 0` / `min_theta_plus < 0` (real
   trapped surface), not just `min_chi`.
4. **Collapse run needs AMR** (`max_level ≥ 4`, ChiTagger) to resolve the forming
   puncture; keep sponge + Ψ₄ extraction.

**Fast next test (before the full Phase-7 redesign):** re-solve config B's ID with
`b0 = 2` (keep the 2 lumps as support), establish with the pump igniter, then ramp
`support_strength → 0`. If it forms a horizon (`max_ah_r > 0`), the trigger works
and we proceed to the scale-matched Phase-7 version for the clean GW burst; if it
still disperses, the matter-binding fix (step 2) is the blocker.

### Phase 8 DONE (2026-07-13) — genuine rotating-eigenstate support: Q-torus retains charge ~2× longer, no blow-up

Addresses the root cause named in the intro ("a spherical Q-ball twisted into an
m=1 torus is not a stationary eigenstate, so Q_sphere drains"). Instead of
twisting a 1D radial soliton by `(sin θ)^m`, we now solve the **genuine 2D
spinning Q-ball eigenstate** `Φ = f(ρ,z) e^{i(m φ − ω t)}` and paint `f(ρ,z)`
directly.

**Solver (new).** `grteclyn_wrapper.grtresna.profiles.qball_torus`: the flat-space
PDE `∂²_ρ f + (1/ρ)∂_ρ f + ∂²_z f − (m²/ρ² + κ²)f + λf³ − μ₆f⁵ = 0`
(κ²=mass²−ω²) solved by a **bordered Newton** that pins the off-axis peak
amplitude and lets ω float (so it cannot collapse to the f≡0 vacuum — the failure
mode of fixed-ω Newton and Newton–Krylov), wrapped in **amplitude continuation**
to hit the target ω. Verified: residual ~1e-9, grid-convergent (f_max→0.09173,
Q→3.098 at N=80/120/160), vanishes on axis, localized (edge/peak ~1e-6).
File I/O (`write/read_torus_profile`) round-trips to ~1e-12.

**GRTresna (new profile 4).** `BosonStarParams.hpp` gained a 2D torus table
loader (`qball_torus_load_table`) + bilinear interp (`qball_torus_interp`);
`lump_winding_modulus` special-cases `profile == 4` to return `amp·f(ρ,z)/f_max`
with NO `(sin θ)^m` factor (the 2D profile already carries the toroidal
structure). `read_params` inherits `qball_profile_path` for profile 3 **or** 4.
Rebuilt clean.

**ID driver (new).** `scripts/wormhole/id/solve_torus.{py,sh}`: throat-free,
centered, full-box (all-Sommerfeld ⇒ GRTresna centres at L/2, `bh1_bare_mass=0`
⇒ flat Bowen-York psi=0), single profile-4 winding lump. Solves the eigenstate,
tabulates it, GRTresna re-solves the constraints, exports a centered gridinit.
Compact couplings (mass 0.5, λ170, μ₆14450), **ω=0.25** (thick-wall; ω=0.05 is
extreme thin-wall and does not fit the box — the solver detects box-escape and
errors). ID converged **Ham 0.095 % (normal) / 0.62 % (exotic), Mom 0.48 %**.

**Evolution.** `wormhole_case.py` gained `--full-box` (center z=L/2,
`lo_boundary 1 1 1`) for a z-symmetric centered object; reuses `--gridinit`
override. Passive isolated evolution of the EXOTIC torus (matches the
`ComplexExoticScalarField` evolution sign), ml=0, L=64/dx=0.5, t→30.

| t | **Q-torus** `|Q_sphere/Q0|` | twisted-sphere baseline | notes |
|---|------|------|-------|
| 6 | 0.98 | 0.99 | tie |
| 10 | 0.93 | 0.92 | tie |
| 15 | **0.89** | 0.56 | torus pulls ahead |
| 20 | **0.80** | 0.27 | ~3× |
| 30 | **0.40** | 0.05 | ~8× |

- **Charge-retention half-life ≈ 27–28 vs the baseline's 13–16 — roughly
  doubled**, and the violent **t≈13.5 dynamical blow-up is GONE** (the twisted
  2-lump constellation NaN'd at 13.5; the torus runs clean to t=30). `min_chi`≈1,
  `max_ah_r=0` (no horizon), Ham L2 bounded ≤3.5e-3 throughout — physical, not
  numerical.
- **Honest caveat:** the torus still slowly spreads (confined_frac 0.81→0.43 by
  t=20; rms 8.3→13.3). Two reasons: (1) it is a **flat-space** eigenstate, only
  approximately stationary once gravity is on; (2) it is **exotic/phantom**, whose
  repulsive self-gravity actively unbinds it. So this is a large improvement, not
  a permanent soliton. Next levers: solve the **self-gravitating** rotating
  eigenstate (metric + field together), and/or evaluate a **normal** (non-exotic)
  torus (needs a non-exotic complex evolution branch) to isolate the phantom-gravity
  contribution.
- **Verdict:** confirms the intro's thesis — a genuine stationary eigenstate is
  the right object; twisting spheres was the root cause of the drain. The
  eigenstate + profile-4 painting pipeline is now the foundation for the
  scale-matched Phase-7 throat (swap the twisted-lump support for a hugging
  **ring of Q-tori** or a central rotating Q-torus core).

**Artifacts:** ID `runs/rotating_torus_id/torus_m1_om0p250_kappa1p00_dx0p5_L64_lam170_mu614450{,_exotic}/`
(gridinit + `qball_torus.dat` kept); evolution
`runs/rotating_wormhole/evo_..._torus_iso{,_t30}/output/data/collapse_diagnostics.dat`.
Reproduce: `EXOTIC=1 TORUS_OMEGA=0.25 solve_torus.sh 2` then
`wormhole_case.sh --gridinit <…_exotic/initial_data.gridinit> --full-box
--omega 0.25067 --m 1 --dx 0.5 --box-size 64 --max-level 0 --stop-time 30
--mass 0.5 --lambda 170 --mu6 14450 --no-frames --run-suffix torus_iso_t30 --gpu 1`.

---

## Phase 1 — Wire the trajectory pump into RotatingWormholeCollapse

### 1.1 Add the pump hook to `ComplexExoticScalarField`

- Mirror `ComplexScalarField`: add `RLMatterPumpParams m_pump{};` member, a
  constructor argument (defaulted so existing call sites compile unchanged),
  and copy the pump source block from
  `Source/Matter/ComplexScalarField.impl.hpp` lines ~146–215 into
  `ComplexExoticScalarField.impl.hpp`'s RHS (`add_matter_rhs` /
  `matter_rhs`-equivalent), including both the PD target-soliton mode
  (`k_p > 0`) and the legacy open-loop mode, and the `governor` scaling.
- **Phantom-sign caveat:** the pump drives `Pi`/`Pi2` toward the target; the
  exotic class flips the *stress tensor* sign, not the field EOM sign — copy
  the source term verbatim, do NOT add an extra sign flip. Verify by the smoke
  test in 1.5.
- Keep it opt-in: `m_pump.num_sites < 1` ⇒ exact no-op (this check already
  exists in the copied block).

### 1.2 Parse trajectory + pump params in the wormhole example

In `Examples/RotatingWormholeCollapse/SimulationParameters.hpp`, copy from
`Examples/RadialRecipe/SimulationParameters.hpp` (~lines 94–123 plus the
`rl_pump_*` block it references):

- `trajectory_mode` (0=off default, 1=on), `trajectory_num_lumps`,
  per-lump `trajectory_lump<k>_{R0,omega_rot,phase0,tilt_theta,tilt_phi,well_depth,v_rad}`,
  shared `trajectory_{A_breath,omega_breath,z_amp,omega_z,well_width}`,
  `trajectory_pump_frequency`.
- Pump gains/limits: `rl_pump_kp`, `rl_pump_kd`, `rl_pump_width`,
  `rl_pump_target_profile`, `rl_pump_target_width`, `rl_pump_target_amp`,
  `rl_pump_max_amplitude`, governor `rl_l2_ham_governor_{center,width}`.
- **NEW (needed for Phase 6):** pump ramp schedule
  `pump_ramp_t_start` (default −1 = never), `pump_ramp_t_end`,
  `pump_ramp_floor` (default 0). Multiply every site amplitude by
  `r(t) = 1` for `t < t_start`, linear → `floor` over `[t_start, t_end]`,
  `= floor` after. Implement in the per-step update (1.3), NOT in the kernel.

### 1.3 Per-step trajectory update + pump construction

Mirror `Main_RadialRecipe.cpp` ~lines 38–110 in
`Main_SupportedWormhole.cpp` / `SupportedWormholeLevel.cpp`:

1. At init, if `trajectory_mode == 1`: seed `RLRuntime::g_lump_state` with the
   t=0 `TrajectoryEvaluator` positions; set per-lump pump amplitudes from
   `well_depth`; set `rl_pump_width = trajectory_params.well_width`.
2. Each coarse timestep (level 0, before RHS): evaluate
   `TrajectoryEvaluator` at the new time, write centres + amplitudes
   (× ramp factor from 1.2) into the global lump state.
3. In the `complex_scalar` branch of `specificEvalRHS`
   (`SupportedWormholeLevel.cpp` ~line 241): build the pump via a local copy of
   `build_rl_pump()` (port from `RadialRecipeMatterDispatch.hpp` ~lines 70–130;
   it reads the global lump state + governor) and pass it into the
   `ComplexExoticScalarField` constructor.

### 1.4 Pump-work diagnostic (required for scoring)

Add one column to `collapse_diagnostics.dat`: `pump_work` = level-0 volume
integral of the pump source term contracted with the field momentum
(instantaneous power), plus its running time-integral `pump_energy`. Follow the
existing `Q_sphere` integral pattern (level-0 synchronised grid). Update the
column-count note in the README and the parser
`grteclyn-wrapper` collapse plotting script if it validates column counts.

### 1.5 Build + no-op regression

- Build CPU + GPU: `main3d.gnu.ex`, `main3d.gnu.MPI.CUDA.ex`
  (see `grteclyn-wrapper/scripts/wormhole/build/`).
- **Regression:** re-run the passive κ=1 smoke with `trajectory_mode = 0` and
  diff `collapse_diagnostics.dat` against a pre-change run — must be
  bit-identical (or numerically identical) since the pump is a no-op.

**Acceptance (Phase 1):** builds pass; no-op regression identical; a 64³ t=2
smoke with `trajectory_mode=1`, 1 lump, tiny `well_depth=0.01` runs NaN-free
and shows nonzero `pump_work`.

---

## Phase 2 — Multi-lump orbital initial data (GRTresna)

### 2.1 Extend `solve_kappa_family.py` for constellations

GRTresna already parses `num_lumps` and per-lump `lump<k>_*` including
`velocity` and `exotic` — the extension is params *generation* only:

- New env vars: `NUM_LUMPS` (default 1 = current behaviour),
  `ORBIT_RADIUS` (R₀), `ORBIT_OMEGA` (ω_orb), `ORBIT_PLANE` (default z=0).
- For k = 0..N−1: place lump k at angle `2πk/N` on the circle of radius R₀ in
  the orbital plane (centred on the throat), with **tangential velocity**
  `v_k = ω_orb · R₀ · t̂_k` (t̂ = unit tangent). Set `lump<k>_exotic = 1`,
  `lump<k>_winding = 0` (orbital angular momentum replaces winding — that is
  the point), `lump<k>_mode = 0`, Q-ball profile (`profile = 3` +
  `profile_path`) with the stiff naturally-small couplings from Rung 1a′
  (φ_c ≈ 0.1; e.g. μ=0.5, λ=170, μ₆=14450 for the ω=0.05 eigenstate — reuse
  the tabulated profile the ODE wrapper already emits).
- Keep the singular-ψ throat term exactly as in the single-lump solve.
- Tag directories `..._nlump<N>_R<R0>_worb<ω>` so families coexist.
- Total exotic amplitude must respect the Lichnerowicz bound: start with
  per-lump amp 0.05–0.1; if the solve diverges (Ham ≫ 1%), halve per-lump amp.

### 2.2 Solve and verify

- Solve 2-lump ID at N=128 (dx = 0.5), R₀ ∈ {6, 10}, ω_orb = v/R₀ with
  v ∈ {0.05, 0.1}.
- **Acceptance:** Ham < 1%, Mom < 1%; `.gridinit` produced; structural check
  (adapt `tests/grtresna/test_rotating_wormhole_kappa.py`): two |Φ|² maxima at
  ±R₀, two-channel Π populated, throat χ dip present.
- Nonzero net `J_z` check: the two-lump tangential momenta should source
  frame dragging — verify `J_z ≠ 0` in the t=0 diagnostics of Phase 3.

---

## Phase 3 — Passive constellation control run (pump OFF)

One run per ID from 2.2, `trajectory_mode = 0`, t → 30, using
`wormhole_case.sh` (extend it with `--nlump/--orbit-radius/--orbit-omega` to
locate the new gridinit tags; ALWAYS via the script — never the bare binary,
lesson L6):

- Measure: `Q_sphere(t)` half-life, `min_chi`, `max_ah_r`, `J_z(t)`,
  lump barycenters (do the lumps orbit at all, fly apart from phantom
  self-repulsion, or fall through the throat?).
- This is the **science control**: it answers "do phantom lumps orbit a
  throat passively?" — expected NO (repulsive self-gravity), but the failure
  mode (unbind vs infall vs disperse) sets the pump's job description.

**Decision gate:** record the passive half-life t½_passive. Everything later
is judged against it AND against the single-torus baseline (t½ ≈ 13–16).

---

## Phase 4 — Hand-tuned pumped smoke (the falsification test)

Cheapest decisive experiment — do this BEFORE any campaign:

- Same 2-lump ID; `trajectory_mode = 1` with orbits **matched to the ID**:
  `R0` = ID's R₀, `omega_rot` = ID's ω_orb, `phase0` = lump placement angles,
  `tilt = 0`, `v_rad = 0`, `well_width ≈` lump width.
- PD mode: `rl_pump_kp > 0` (start 0.1), `rl_pump_kd` (start 0.01),
  `rl_pump_target_profile` = the Q-ball profile, `rl_pump_target_amp` = ID amp.
  Governor centred at ~3× the run's initial L2_Ham.
- t → 30, N=128/dx=0.5, ml ≤ 1. Sweep k_p ∈ {0.05, 0.1, 0.5} if needed.

**Acceptance / decision:**
- **PASS:** `Q_sphere(30)/Q_sphere(0) ≥ 0.5` (vs ≤ 0.05 passive), L2_Ham
  bounded (< 5× initial, no runaway), `min_chi` holds a plateau < 1 (throat
  exists), `pump_energy` finite and recorded. ⇒ proceed to Phase 5.
- **FAIL after gain sweep:** the pump cannot beat the instability timescale ⇒
  falsifies the approach cheaply; fall back to Rung 1½ (stationary solve) /
  Rung 3 (GW burst characterisation) per RotatingWormholePlan.

---

## Do we still need the MAP-Elites campaign? (decision, 2026-07-13)

**Short answer: NO — it is not on the critical path to the headline result.**
Phase 4c already produced a *working* stable configuration (**B**: 2-lump,
R0=6, ω_orb=0.1, amp=0.06 [κ=0.6], k_p=1.0), which holds an open throat with
no horizon out to t=30. The publishable payoff is **Phase 6** (collapse on
command + GW signature), and that needs only *one* good "hold" config to then
ramp the pump down. So the plan is now:

1. **Small hand-tuned refinement around B** (cheap, ~4–6 runs, no campaign) to
   flatten B's slow secular drift — B holds a clean min_chi≈0.86 plateau over
   t≈8→18 but then gradually opens to ≈0.98 by t=30. For a clean Phase-6
   "hold, then trigger" experiment the hold phase should be genuinely
   stationary so the collapse is unambiguously caused by the ramp, not by
   pre-existing drift. Sweep κ∈{0.55,0.6,0.65}, target_amp match, and a
   governor/k_p tweak; pick the flattest, longest plateau.
2. **Go straight to Phase 6** on that refined elite.

**What MAP-Elites would add (now OPTIONAL, paper-enrichment only):** a
*stability map* — the boundary surface in (per-lump amplitude, R0, ω_orb, k_p)
between the three fates we have already seen by hand: collapse-to-BH
(supercritical: baseline, A, C), stable open throat (B), and dispersal-to-flat
(D). That map is a nice figure ("which engineered rotating throats can be held
open, and at what pump cost") and would locate the *most* stable / least-pump
corner, but it is **not required** to demonstrate support or collapse. Run it
only if we want the map for the paper, or if the hand-refinement in step 1
fails to find a stationary-enough hold. The campaign spec below is kept for
that contingency.

---

## Phase 5 — MAP-Elites stability campaign (OPTIONAL — see decision above)

Reuse the QD pipeline (`scripts/campaigns/qd/run.sh` + `lib/`); the deltas:

- **Evaluation target:** the RotatingWormholeCollapse binary via
  `wormhole_case.py` (not RadialRecipe). Per-eval: GRTresna multi-lump solve
  (Phase 2 generator) with **postload constraint gate Ham/Mom ≤ 5%**, then
  evolution to t = 30 at N=128/dx=0.5 (search res), scoring from
  `collapse_diagnostics.dat` — no Ψ₄/frames needed at search res.
- **Genome (~4 + 7·N_lumps dims):** N_lumps ∈ {2,3,4}; per-lump
  `{R0, omega_rot, phase0, tilt_theta, tilt_phi, well_depth, v_rad}`; shared
  `{k_p, k_d, well_width, A_breath}`. ID placement derives from the same genes
  (lump k at phase0_k, radius R0_k, tangential velocity ω·R0) so the pump
  corrects rather than fights the ID.
- **Fitness:**
  `score = w1·Q_retention(t=30) + w2·throat_persistence − w3·constraint_growth − w4·pump_energy_norm − w5·horizon_flag`
  where `Q_retention = Q_sphere(30)/Q_sphere(0)` clipped [0,1];
  `throat_persistence` = fraction of the run with `min_chi ∈ [0.4, 0.95]`
  (i.e. a throat exists — neither collapsed nor relaxed flat);
  `constraint_growth = max(L2_Ham)/L2_Ham(0)` soft-capped;
  `pump_energy_norm` = `pump_energy` / initial matter energy;
  `horizon_flag` = 1 if `max_ah_r > 0` (a horizon here means failure of
  *support*, not the desired triggered collapse). Start
  (w1..w5) = (1.0, 0.5, 0.3, 0.3, 1.0); tune after 1 batch.
- **Descriptors (8×8 archive):** axis 1 = total orbital angular momentum
  proxy `Σ ω_rot·R0²` (binned); axis 2 = `pump_energy_norm` (binned, log).
  The interesting corner is high-retention/low-pump-work — configs approaching
  a *passive* equilibrium (input to Rung 1½).
- **Budget:** 200 evals, 8 GPUs (~2 h at RadialRecipe-like throughput; wormhole
  evals are similar cost). Then HQ promotion of the top elite
  (256³ or dx=0.25, t→60, Ψ₄ ON via `GRTECLYN_PSI4=1` in `promote_common.sh`).

**Acceptance:** archive coverage > 10%; ≥ 1 elite with
`Q_retention ≥ 0.7` at HQ replay. Log everything in a new journal
`research/rotatingwormhole/OrbitalPumpJournal.md`.

---

## Phase 6 — Pump-ramp collapse experiment (the paper arm)

On the HQ-validated elite(s):

1. **Hold:** evolve with full pump to t_hold (quasi-stationary: `Q_sphere`
   flat over ≥ 10 t, `min_chi` plateau).
2. **Trigger:** `pump_ramp_t_start = t_hold`, ramp over Δt = 2, to
   `pump_ramp_floor ∈ {0.5, 0.25, 0}` — the rotating analogue of the static
   article's `S_support` cut. One run per floor value.
3. **Measure:** branch outcome per floor — throat pinch (trapped surface:
   `min_theta_plus < 0`, `max_ah_r > 0`) vs rarefactive opening
   (`min_chi → 1`); `J_z` through the transition; Ψ₄ ℓ=2, m=0/±2 from both
   the orbiting phase (binary-like periodic signal at 2ω_orb) and the
   collapse/dispersal transient. Extraction radii 12/16/20/24 as in the
   existing wormhole `plot_run.sh`.
4. Convergence pair (dx vs dx/2) + box-size pair (L=64 vs 96) on the headline
   case, per the earlier boundary-study protocol.

**Deliverable:** "An actively supported rotating wormhole: minimal-control
stabilization map (MAP-Elites), and collapse-on-command GW signatures vs pump
floor" — with the passive Rung 1½ solve as the exact companion arm if the
low-pump-work corner of the map points to one.

---

## Phase 7 — Scale-matched throat + orbiting lump ring (NEXT DESIGN, 2026-07-13)

**Motivation (from the embedding-diagram review).** Config B is geometrically a
*tiny neck buried inside two huge clouds*, not a clean "throat + orbiting
supports": throat areal radius `b0 = 0.5`, but each Q-ball lump has half-max
radius ≈ 5.7 (10% radius ≈ 8.4) centred at `R0 = 6`. So (i) the lumps do not
actually *orbit* anything — most of their mass sits on top of the centre; and
(ii) the object is better described as one exotic blob with a coordinate pinch.
Physically we want the throat and the support lumps to be of **comparable size**,
ideally a **larger throat encircled by many small lumps** as the support/engine
ring — an interpretable, genuinely orbiting configuration.

Visualisation of the current (mismatched) state:
`grteclyn-wrapper/src/grteclyn_wrapper/visualisation/wormhole_embedding.py`
(3-panel: 3-D funnel zoomed on the neck, side-on cross-section with the exotic
matter's radial reach, equatorial |Φ| map — all to scale, code units).

**Hard constraint (do NOT design around it).** A traversable throat stays open
only if NEC-violating (exotic) stress exists **at the throat itself**
(Morris–Thorne flaring-out). This is a theorem, not a tunable. Our own data
confirms it: Phase 4c **config C** (lumps spread out to `R0 = 10`) *starved the
throat and collapsed to a BH at t ≈ 7.3*. Therefore **the support lumps' fields
must still blanket the neck** — small lumps orbiting *far* from a *big empty*
throat WILL collapse.

**"Swallowed" caveat.** A wormhole has no horizon, so matter is not captured like
a BH — but matter not in a bound orbit **plunges through the throat to the other
universe**. Supports must be (a) in genuinely bound orbits AND (b) close enough
that their fields overlap the neck.

**Target design (satisfies both).** A **hugging ring**: throat `b0 ≈ 2–3`
(comparable to lump size); `N = 6–12` **compact** Q-balls on a ring at
`R0 ≈ b0 + r_lump` so neighbours nearly touch → an exotic ring that collectively
blankets the neck while each lump is individually compact and orbiting.
Optional hybrid: keep a **central exotic core** (a static/winding throat field)
for the flaring-out support and use the ring purely for spin/dynamics.

**Build steps.**
1. **Shrink the Q-ball → compact profile.** Re-solve the flat-space Q-ball ODE
   (`grteclyn_wrapper.grtresna.profiles.qball_ode`) with a **larger field mass
   `m`** and **weaker quartic `λ`** so the soliton loses its fat flat top and
   gets a half-max radius comparable to the target `b0` (order 1–2, not ~6).
   `MASS`, `LAMBDA`, `MU6`, `LUMP_OMEGA` must match on the solve and evolution
   sides (gotcha 5). Sanity-check the new half-max/10% radii by re-running the
   embedding script (it reads `qball_profile.dat`).
2. **Bigger throat in the base ID params.** Set `wormhole_throat_radius` (b0) in
   the GRTresna base params (`params_rotating_wormhole_test.txt`) to the target
   (2–3). Larger throat = more curvature/"mass" to support → expect to need more
   total exotic amplitude, but respect the amplitude ceiling (gotcha 7): raise
   `N` rather than per-lump amp.
3. **New constellation ID solve.** `NUM_LUMPS=N ORBIT_RADIUS=R0 ORBIT_OMEGA=ω
   MASS=.. LAMBDA=.. MU6=.. RES_N=256 EVO_L=128 solve_kappa_family.sh <kappa> 8`.
   Postload constraint gate **Ham/Mom ≤ 5%** (gotcha 7: if Ham stuck ≫1%, reduce
   per-lump amp / N, don't fight the solver).
4. **Stability hunt (mirror Phase 4c).** Evolve with `k_p = 1.0` igniter;
   sweep `{b0, N, R0, ω_orb, kappa}` for a flat `min_chi` plateau, no horizon to
   t = 30, `Q_retention ≥ 0.5`. Mass reduction (kappa) remains the main
   sub-critical lever (Phase 4c lesson).
5. **Then Phase 6 collapse-on-command** on the scale-matched winner (support-
   strength ramp; AMR ON — see the AMR note below).

**Acceptance (Phase 7):** lump half-max radius ≈ `b0` (comparable scales);
`R0 ≈ 1.5–2 b0` (a ring that hugs, not engulfs, the neck); ≥ 6 lumps visibly
orbiting outside the throat; stable open throat (min_chi plateau, no horizon) to
t = 30; embedding diagram shows a clean "big neck + small orbiting ring."

**AMR note (carried into Phase 6/7 collapse runs).** The collapse arm forms a BH;
at `max_level = 0` (uniform dx = 0.5) the puncture is unresolved → NaN. The
static `SupportedWormholeCollapse` that DID extract GW through collapse used
`max_level = 5` (ChiTagger refinement) with `puncture_tracking.enabled = 0`
(no tracker needed — collapse is central; 1+log slicing handles the singularity).
Our `RotatingWormholeCollapse` has the same gauge and the ChiTagger; the collapse
run must set `--max-level ≥ 4`. Caveat: the rotating case has **no x/y/z symmetry**
(spin breaks reflection symmetry — the static octant run computed 1/8 the box), so
the full-domain N=256 base grid + refined puncture patches is memory-heavy
(~60 GB already at ml0); use the fewest AMR levels that resolve the puncture and
tag tightly on the centre.

**Open knobs to fix before the solve:** support architecture (hugging ring vs
central core + ring vs hybrid); target `b0`; lump count `N`; the compact-Q-ball
couplings `(m, λ, μ₆)`.

---

## Standing gotchas (bind ALL phases — from hard lessons)

1. **`Q_sphere` (col 21) is the confinement metric. NEVER `rho_sum` (col 18).**
2. **L2 (gridinit dx):** evolve at the gridinit's native level-0 dx; AMR below
   is OK only with `regrid_interval` given **one entry per level** (AMReX
   aborts otherwise).
3. **L6 (disk):** never launch the bare binary; always `wormhole_case.sh` /
   `run_rotating_wormhole.sh` so the plotfile-deleting consumer sidecar runs
   *during* the evolution. `plot_interval ≥ 10` at search res.
4. **Gauge kick:** maximal-slicing (K=0, α=1) ID into 1+log evolution spikes
   `max_K` ≈ 2.2 at t ≈ 0.5. Keep it (all baselines share it) but note it; if
   the pump fights it, try `--initial-lapse-type 1`.
5. **Param consistency:** `phantom_mass`(evolution) == `scalar_mass`(solve);
   same λ, μ₆ on both sides; shared convention
   `V = ½μ²|Φ|² − ¼λ|Φ|⁴ + ⅙μ₆|Φ|⁶`.
6. **Coupled potential only:** the complex branch must use
   `ComplexScalarPotential` on |Φ|² — never per-component potentials (breaks
   U(1)/Noether charge at t=0).
7. **Exotic amplitude ceiling:** total painted exotic amp near 0.1-equivalent;
   Lichnerowicz divergence shows up as Ham stuck ≫ 1% — reduce per-lump amp,
   don't fight the solver.
8. **Pump = Bianchi violation by design:** L2_Ham/Mom drift proportional to
   pump work is expected; it is a *scored penalty*, not a bug — but runaway
   growth (> ~5× initial) means gains too high / governor mis-set.
9. Commit style: `feat:`/`docs:` prefixes as in recent history; journal every
   run outcome in `OrbitalPumpJournal.md` (tables: params, Q_retention,
   half-life, constraint growth, pump energy, verdict).

---

## C++ in-code Weyl4 extraction — root-cause analysis + DEFERRED (2026-07-14)

**Status: DEFERRED — upstream feature is incomplete.** The GRTL collaboration lists
**"Weyl scalar / CCE extraction" as 🔧 In progress** in GRTeclyn's development status
(alongside particle-based diagnostics, apparent-horizon finder, ADM quantities). So
the in-code C++ Ψ₄ extraction here is expected to be unreliable — it is a partial
port, not a bug in our wiring. **We skip the C++ fix and use the validated Python
post-hoc extraction (`small_data/psi4_mode_l2m0.dat`) for all physics.** Revisit the
port below once upstream marks Weyl/CCE extraction ✅ Ported.


**Symptom.** `data/Weyl4_mode_2{0,1,2}.dat` (dense in-code C++ Ψ₄) disagrees with
the validated Python post-hoc extraction (`small_data/psi4_mode_l2m0.dat`):
- ~8× larger RMS at r=12, and **uncorrelated** (corr ≈ 0.07) with the Python signal.
- Large `O(0.6)` value at t≈0 where a near-stationary ID should radiate ≈0.
- `O(1)` `m=0` imaginary part (Python gives ~1e-5).
- ~80 per-shell `1e51`-scale blowups at r=18/24, **clustered at violent/regrid
  phases** (t≈6.1, 7.8, 9.0–9.2 = support-off collapse onset, and t≈37.7).

**Verified NOT the cause (all standard/correct, identical to the working BinaryBH
example):**
- Extraction center: `extraction_center = 32 32 32` is read into
  `extraction_params.center` (`SimulationParametersBase.hpp:123`) — correct throat
  center. Python sidecar uses the same 32,32,32 (`wormhole_case.py:618,624`, full-box).
- Modes `2 0 / 2 1 / 2 2`, `num_points_phi=24`, `num_points_theta=37` — parsed OK.
- Component map: `WeylExtraction` reads derived `Weyl4` comps {0,1} = {Re,Im};
  `Weyl4WithMatter::compute_mf` writes Re,Im from `out_comp` — correct.
- Projection math (`SphericalExtraction::add_mode_integrand`): canonical
  `∫ (r·Ψ₄)·conj(sYlm)/r² · r²sinθ dθdφ` with `spin_Y_lm` — the same code BinaryBH
  uses successfully. Both Python and C++ ultimately consume the *same*
  `Weyl4WithMatter` field.

**Root cause — GRChombo comparison (the decisive finding).** GRChombo's
`Examples/BinaryBH/BinaryBHLevel.cpp::specificPostTimeStep` (lines 142–169) does a
sequence that GRTeclyn's port **dropped**:
```cpp
fillAllGhosts();                                             // fill CCZ4 state ghosts, all levels
BoxLoops::loop(Weyl4(center,m_dx,formulation),
               m_state_new, m_state_diagnostics,
               EXCLUDE_GHOST_CELLS);                         // Weyl4 -> STORED diagnostic, VALID cells only
if (m_level == min_level) {
    m_gr_amr.m_interpolator->refresh(false);
    m_gr_amr.fill_multilevel_ghosts(                         // fill diagnostic Weyl4 ghosts across levels
        VariableType::diagnostic, Interval(c_Weyl4_Re,c_Weyl4_Im), min_level);
    my_extraction.execute_query(m_gr_amr.m_interpolator);    // interpolate STORED, consistent field
}
```
GRTeclyn instead computes Weyl4 **on the fly** inside the interpolator via
AMReX `derive("Weyl4", time, s_num_ghosts=2)` (see
`ParticleInterpolator.impl.hpp:367`), with no stored diagnostic and no explicit
multilevel diagnostic-ghost fill.

**Ruled out:** AMReX's `derive` is *not* under-filling input ghosts —
`AMReX_AmrLevel.cpp:1638` sets `ngrow_src = ngrow + 2` (from the Weyl4 `boxMap`
`grow(box,2)`), so the CCZ4 source is FillPatched to 4 ghosts, enough for the
4th-order Weyl4 stencil. Center/modes/components/projection math all verified
identical to the working BinaryBH. So this is **not** a one-line ghost-count patch.

**Remaining suspects (need runtime localization).** The `1e51` spikes cluster at the
violent, heavily-regridding phases (support-off collapse onset t≈6–9; t≈37). Likely
the quartic *particle* interpolation onto the derived Weyl4 near freshly-created
coarse–fine boundaries, and/or `getData(time)` time-interpolation against a
just-regridded level. GRChombo sidesteps both by extracting from a stored,
multilevel-ghost-filled diagnostic rather than an on-the-fly per-query derive.

**Fix (chosen: diagnostic-MultiFab port of the GRChombo path) — substantial:**
GRTeclyn currently has **no** diagnostic state or `fill_multilevel_ghosts` (grep of
`Source/GRTeclynCore` finds none), so the fix is a real port, not a patch:
1. Add a diagnostic state (or reuse spare state comps) holding `Weyl4_Re/Im`.
2. In `SupportedWormholeLevel::specificPostTimeStep`, before extraction: fill CCZ4
   ghosts, compute Weyl4 into the diagnostic on valid cells for every level,
   average-down + fill multilevel ghosts of the diagnostic.
3. Point `WeylExtraction` at the stored diagnostic (state-var interp path) instead
   of the derived-var `derive` path.
Recommended pre-step: a short (t≈2) run with `write_extraction = 1` to dump the raw
Weyl4 spheres and localize *where* on the sphere the `1e51` originates — cheap
runtime evidence to confirm the mechanism before the full port + rerun.

**Interim policy (in force).** Physics uses the Python extraction
(`plot_diagnostic.sh` → `psi4_mode_l2m0.dat`). The dense C++ files are debugging-only
until the port above is implemented and re-validated (Im(m=0)≈0, matches Python,
no `1e51`).
