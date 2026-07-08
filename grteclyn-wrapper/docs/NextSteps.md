# Next Steps — Critical Review & Implementation Plan

Senior-scientist review of the FTL / GW / self-grav campaigns (2026-07-07), plus a
prioritized, actionable implementation plan. Companion to
[`README.md`](README.md) (§Main results) and [`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md).

---

## 1. Verdict

The **engineering** is strong: constraint-solved initial data in the loop, a
validation ladder that has caught five distinct classes of false positives,
honest negative results (boson 0/94), and documented reward-hacking closure
(gw_beam v3 → v4). That is more methodological discipline than most
computational-physics projects have.

The **physics framing has a soft core**: the headline claim ("genuine
gauge-invariant FTL shortcuts exist") rests on a baseline definition that is
not actually gauge-invariant, the probe has never been calibrated against a
known analytic answer, and the most theorem-relevant quantity (ANEC along the
shortcut ray) is never computed. Fixing those three things makes the results
defensible; skipping them invites a one-paragraph referee rejection.

The README's own finding #3 ("search design is the dominant lever") is
correct, and should be sharpened: **the search is done discovering; the
instrumentation now has to catch up to the claims.**

---

## 2. Validity assessment

### What is solid

| Strength | Why it matters |
|----------|----------------|
| Constraint-solved initial data + post-load gate | Every `gpu_ok` candidate is a genuine Einstein solution at t=0 — many "warp drive" papers lack this |
| Evolving 4D null-ray tracing | Correct *class* of probe (vs coordinate Dijkstra, correctly demoted after the wormhole paradox) |
| Falsification tiers, HQ ladder, geodesic trust flags | Caught eval 008 (low-res artifact), SH 151/101 (gauge), gw_beam v3 (numerical bomb) |
| Controlled boson-vs-scalar comparison (0/94 vs 32/92) | A real negative control inside the search itself |
| Consistency with theory | Exotic matter → transient shortcuts, canonical matter → none. Matches Olum 1998 / Gao–Wald 2000: superluminal shortcuts require ANEC violation |

### Core problem 1 — the `f_geo` baseline is coordinate-dependent

`evolving_geodesic.py` defines `t_flat = |x_end - x_start|` — a **coordinate**
distance (see module docstring, `f_geo = max(0, (t_flat - t_travel)/t_flat)`).
The null ray itself is gauge-honest; the *comparison* is not. "Beats flat
space" only means something if:

- the ray endpoints sit where the metric is ≈ Minkowski in these coordinates
  (unchecked — box L=64–128 containing structure of RMS radius up to ~18), and
- coordinate time at the endpoints ≈ proper time of an asymptotic observer
  (unchecked).

If ψ, α, β deviate from (1, 1, 0) at the box edges, percent-level `f_geo` can
be manufactured or masked purely by the spatial gauge. The ~10–20% signals are
probably above that floor — but the floor has never been *measured*. The
theorem-clean formulation (Gao–Wald) compares arrival against a reference ray
between the **same endpoints routed through the far field of the same
spacetime** (same ADM mass), not literal Minkowski. Additionally,
`max(0, ...)` clips the Shapiro-*delay* regime, discarding the built-in sanity
check that positive-mass configs should show negative would-be `f_geo`.

### Core problem 2 — no positive or negative control has ever run through the pipeline

The probe has only ever measured its own discoveries. Until it demonstrably
(a) reports a Shapiro **delay** on Schwarzschild-like data and (b) recovers the
**known** transit-time advantage of an analytic Alcubierre/Krasnikov metric
(cf. Clough et al. 2024 warp-collapse evolutions), every champion number is
produced by uncalibrated instrumentation.

### Core problem 3 — the pump breaks Einstein-consistency wherever it is active

The RL/feedback pump injects a non-conservative source term into the scalar
field equation. Unless its stress-energy is included in the constraints and in
the CCZ4 matter sources, pumped runs are **not solutions of the Einstein
equations** (Bianchi identities violated; constraint growth guaranteed). Fine
for "actuated engineering" exploration; fatal for any physics claim. Note the
self-grav "fixed" star in `SELFGRAV_HANDOFF.md` still runs with the pump on — a
truly correct stable-branch boson star needs **no** pump, and the residual leak
(RMS 2.1 → 4.2 by t=16) with the pump active suggests remaining inconsistency
in the initial-data painting, not a missing trap.

### Secondary concerns

- **Exotic matter = ghost scalar** — classically unstable, quantum-forbidden in
  sustained form. The honest framing is "how much NEC violation buys how much
  shortcut in full nonlinear GR," not "FTL exists."
- **"All FTL configs converge to ~20% peak f_geo (resolution ceiling)"** — if
  the observable has a resolution ceiling, that is a resolution-dependent
  measurement, not convergence. Need ≥3 resolutions + demonstrated convergence
  order for `f_geo` itself.
- **Global L2 constraint gates (1e-2) hide local violation** in the channel
  where the shortcut actually lives.
- **GW extraction at r = 8/12/24 for sources of size ~10–18 is near-zone**, not
  radiation. gw_beam ratios and LIGO matched-filter templates built on this Ψ₄
  are not asymptotic quantities.
- **Linear temporal interpolation of the metric stack** has unquantified error
  vs plotfile cadence.
- **Kaup limit 0.62 vs textbook 0.633** — 2% error in the quantity that decides
  stable-vs-unstable branch selection.
- **Uncommitted work across two repos** (self-grav fix) is a reproducibility
  time bomb.

---

## 3. Implementation plan

Ordered by (value ÷ cost). Each workstream lists tasks, files, and acceptance
criteria. P1–P4 are instrumentation and require **no new campaigns** — they
re-analyze existing champions (eval 122, wormhole 046, spiral 118) or run
single cheap evolutions.

### P1 — Control suite for the FTL probe *(highest value, cheapest)*

Calibrate the instrument before trusting any more readings.

| Task | Detail |
|------|--------|
| Negative control | Schwarzschild (or boosted single star) initial data → probe must report Shapiro **delay** (negative unclipped `f_geo`) of analytically predictable magnitude, at every resolution |
| Positive control | Analytic Alcubierre / Krasnikov metric stack (via `evolving_field_from_analytic_stack`, already exists) with a *known* transit-time advantage → probe must recover it to stated tolerance |
| Cadence convergence | Re-run a champion probe at 2× and 4× plotfile cadence → quantify temporal-interpolation error on `f_geo` |
| Permanent gate | Wire all three into a pytest module so every future probe change is recalibrated automatically |

- **Files:** `src/grteclyn_wrapper/metrics/probes/ftl/evolving_geodesic.py`
  (reuse `evolving_field_from_analytic_stack`), new
  `tests/metrics/ftl/test_probe_controls.py`, analytic metric generators under
  `metrics/probes/ftl/`.
- **Acceptance:** Schwarzschild delay within tolerance of analytic Shapiro
  value; analytic warp `f_geo` recovered within ±X% (fix X from cadence study);
  tests run in CI/pytest preflight.

### P2 — Gauge-honest baseline; re-measure the champions

| Task | Detail |
|------|--------|
| Unclip `f_geo` | Report signed value; keep `max(0,·)` only at scoring time |
| Far-field reference ray | Replace/augment `t_flat = |x_end - x_start|` with a reference null ray between the same endpoints routed through the far-field region of the **same** spacetime (Gao–Wald-style comparison); report both `f_geo_flat` and `f_geo_reference` |
| Endpoint weak-field check | At launch/arrival points, record ψ, α, β deviations from (1,1,0); flag the report `endpoint_weakfield_ok` and gate scoring on it |
| Re-measure champions | eval 122 (`trajectory_5lump_v1`), wormhole 046 (HQ), spiral 118 — offline re-probe from retained metric stacks, no rerun needed where stacks exist |

- **Files:** `metrics/probes/ftl/evolving_geodesic.py`,
  `metrics/probes/ftl/geodesic.py`, report dataclasses, scoring hookup in
  `metrics/score/ftl.py`; docs update in `metrics/README.md`.
- **Acceptance:** champions re-reported with reference-ray baseline + endpoint
  diagnostics; any champion whose signal collapses under the new baseline is
  reclassified (update README results tables honestly).

### P3 — ANEC integrals + quantum-inequality margin *(the novel result)*

The theorem-relevant quantity, and the most publishable number in the project:
an empirical **exchange rate between ANEC violation and shortcut depth** in
full nonlinear GR.

| Task | Detail |
|------|--------|
| ANEC along champion rays | Integrate ∫ T_μν k^μ k^ν dλ along the traced shortcut geodesics (ray path already produced by the probe; T_μν available from matter fields on the grid) |
| Exchange-rate table | For all confirmed champions + several sub-threshold candidates: (ANEC integral) vs (peak/evolving `f_geo`) |
| Ford–Roman QI margin | Compute how many orders of magnitude the champions' negative-energy pulses exceed quantum-inequality bounds — converts "FTL evidence" into an honest physical statement |

- **Files:** new `metrics/probes/ftl/anec.py`; hook into
  `EvolvingGeodesicFtlReport`; analysis notebook/script under
  `scripts/search/`.
- **Acceptance:** ANEC < 0 along every confirmed shortcut ray (theorem
  consistency check); exchange-rate plot; QI margin quoted per champion.

### P4 — Real convergence + channel-local constraints

| Task | Detail |
|------|--------|
| 3-resolution ladder for eval 122 | 128³ / 256³ / 384³ (or 512³ if feasible), same physics; Richardson-extrapolate `f_geo` and quote a convergence order — retire the "resolution ceiling" phrase |
| Channel-local constraint monitor | Max/L∞ Ham & Mom in a tube around the champion ray path per plotfile (not just global L2); add to `small_data/` outputs |
| Constraint-triggered checkpointing | For the self-grav high-res NaN: checkpoint on constraint-spike so blow-ups can be bisected instead of rerun blind |

- **Files:** `scripts/campaigns/hq/` ladder launcher; new diagnostic in
  `metrics/diagnostics/`; consumer hookup in `visualisation/`.
- **Acceptance:** quoted convergence order for `f_geo(eval 122)`; channel-local
  constraint violation demonstrably small relative to the curvature scale
  producing the shortcut.

### P5 — Pump physics accounting / run-class split

| Task | Detail |
|------|--------|
| Audit pump stress-energy | Determine whether the pump term enters CCZ4 matter sources + constraint diagnostics. If not, either add its T_μν or — simpler — formally split run classes |
| Run-class split | Tag every run `physical` (no pump, claims allowed) vs `actuated` (pump on, engineering only); propagate the tag through `EpisodeMetrics`, trajectories, and results tables |
| Documentation | README + lab journals updated so no pumped number is cited as a GR result |

- **Files:** `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp`,
  `ComplexScalarField.impl.hpp` (audit); `grtresna/matter/wiring.py`,
  `metrics/episode_metrics.py` (tagging).
- **Acceptance:** every result table row carries a run class; no `actuated` row
  is presented as physics.

### P6 — Self-grav boson star: pump-free, consistent, committed

| Task | Detail |
|------|--------|
| Pump-free stable star | A correct stable-branch star must sit without a trap. Chase the residual leak in the initial-data painting (profile shape vs generic `sech`, lapse/momentum consistency) rather than compensating with the pump |
| Kaup discrepancy | Resolve 0.62 vs 0.633 (ODE tolerance / isotropic conversion) before trusting branch selection; default `phi_c` to a safely stable 0.06–0.07 (current default 0.08 is marginal) |
| High-res stability | KO dissipation, CCZ4 κ₁/κ₂ tuning, 1+log gauge start, fixed (non-regridding) refinement around the star — per the handoff's own list, now with P4's constraint-triggered checkpoints |
| Commit both repos | GRTeclyn `feature/interstellar` + sibling GRTresna; tag campaign configs and seeds |

- **Acceptance:** pump-free star holds RMS radius within ~10% over t=16 at
  `max_level≥2`; no NaNs to t=16 at `max_level=3`; both repos committed.

### P7 — Wave-zone GW extraction *(prerequisite for any further gw_beam / LIGO work)*

| Task | Detail |
|------|--------|
| Larger extraction domain | Push extraction radii to r ≳ several × source size (bigger box or mesh-refined far zone); extrapolate Ψ₄ r→∞ across radii |
| Full multipole set | Beyond l=2,m=0 — beaming characterization needs the (l,m) spectrum |
| Pause LIGO matched-filter search | Until templates come from wave-zone waveforms with defensible physical scaling |

- **Files:** `visualisation/process_wave/`, `scripts/plot/`,
  `gw_search/README.md` (add prerequisite note).
- **Acceptance:** Ψ₄ amplitude falls as 1/r across extraction radii
  (radiation-zone check) before any beam ratio or template is quoted.

---

## 4. Scoring & search improvements (cross-cutting)

1. **Drop `operational_ftl` (coordinate Dijkstra, weight 400) from scoring** —
   it has produced both false positives and false negatives; keep as a
   diagnostic only. Headline scoring = evolving geodesic + trust gates
   exclusively.
2. **Move toward constrained optimization** — hard feasibility gates + a single
   physical objective, instead of ~15 weighted components per mode. Every
   reward-hacking episode so far came from weighted-sum leakage.
3. **Hold-out validation metrics** — keep at least one physical diagnostic out
   of the objective entirely and check champions against it post hoc
   (adversarial validation against metric-hacking).
4. **ADM mass/momentum tracking** — free integration-quality diagnostic; ghost
   matter permits negative ADM mass, worth knowing per champion.
5. **Proper apparent-horizon finder** (outgoing null expansion) if the current
   trapped-surface veto is a proxy.

### P8 — Geometry expressivity via puncture free data

Extend the forward search from matter-only to matter + Bowen–York puncture
free data (masses, spins, momenta — already exposed in `GRTresnaConfig` but
never searched): genuine Einstein–Rosen bridges, multi-throat configurations,
and exotic-shell-around-throat (Morris–Thorne) campaigns. Full phased plan:
[`PuncturePlan.md`](PuncturePlan.md). Depends on P1/P2 (probe calibration)
before its discovery phase.

---

## 5. What the publishable claims actually are

Not "FTL exists" — exotic-matter shortcuts have been known possible since
Alcubierre/Olum/Visser. The defensible papers in this work:

1. **Methods** — quality-diversity search over constraint-solved initial data
   with a gauge-artifact-rejecting validation ladder ("AI-driven discovery in
   numerical relativity"), with the false-positive table as the centerpiece.
2. **Negative result** — systematic search finds *zero* shortcuts from
   canonical matter (boson 0/94, Lentz arm): numerical support for the
   energy-condition theorems and a direct empirical challenge to
   positive-energy warp proposals.
3. **Quantitative** — the ANEC-violation-vs-shortcut exchange rate (P3) and the
   universal transience result (every shortcut closes via dispersal or horizon
   formation) — dynamics no constructed-metric paper can show.

All three require P1–P4 to be complete first.

---

## 6. Suggested execution order

```
P1 (controls)          ~ days      ── unlocks trust in everything downstream
P2 (baseline)          ~ days      ── re-scores existing champions offline
P3 (ANEC/QI)           ~ 1 week    ── the novel result; needs P1+P2 semantics
P4 (convergence)       ~ 1 week GPU── runs alongside P3
P5 (pump split)        ~ days      ── mostly bookkeeping + audit
P6 (self-grav)         ~ weeks     ── the physics bottleneck; parallel track
P7 (GW wave zone)      ~ 1 week    ── only if gw_beam/LIGO remains a goal
P8 (punctures)         ~ 3–4 weeks ── geometry expressivity; phases 0–1 can
                                      start now, discovery phase after P1+P2
                                      (see PuncturePlan.md)
```

P1 + P2 + P3 together are the minimum bar before writing up any FTL claim.
