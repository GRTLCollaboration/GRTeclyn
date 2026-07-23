# Puncture Campaign — Implementation Plan

Geometry-expressive search via **Bowen–York puncture free data** instead of (and
later combined with) scalar-lump matter. Companion to [`NextSteps.md`](NextSteps.md)
(referenced there as **P8**) and [`README.md`](README.md).

---

## 1. Motivation and physical scope

The matter-first campaigns are expressivity-limited by the Gaussian-lump basis
and by conformally trivial geometric free data: every candidate is "lumps on a
conformally flat slice." But in the conformal (CTT/CTTKHybrid) method, matter
is **not** the only well-posed input — the analytic ψ/A_ij seed terms are
*geometric* free data. GRTresna already ships the full two-puncture Bowen–York
sector (bare mass, spin, momentum, offset per hole), and the wrapper's
`GRTresnaConfig` already exposes all of it. The search has simply never used
these dimensions.

**What punctures buy:**

- **Nontrivial slice topology** — each puncture is an Einstein–Rosen bridge
  (extra asymptotic end). Two punctures = two bridges ("multi-throat").
- **Vacuum geometry knobs** — spin (frame dragging) and momentum (boosted
  geometry) with zero matter, which the lump sector cannot produce.
- **The Morris–Thorne question in full nonlinear GR** — can an exotic scalar
  shell around a puncture throat hold it open / open an exterior shortcut, and
  for how long before horizon recollapse? (The wormhole v21/v22 campaigns found
  throat-*like* conformal channels without topology; this is the real thing.)

**What punctures do NOT buy (be honest in all write-ups):**

- No same-universe traversable handle. Geroch's theorem forbids topology change
  in evolution without singularities/CTCs; an R³-chart grid stays R³-chart. The
  ER bridge itself is non-traversable in vacuum. The searchable science is
  *transient exterior/bridge shortcuts and horizon-delay dynamics with exotic
  shells* — same epistemic class as the existing transient-FTL results.
- Vacuum Bowen–York alone is just a BBH: positive energy, Shapiro **delay**, no
  FTL. Phase 2 (puncture-only) is therefore *infrastructure + calibration*;
  the discovery campaign is Phase 3 (puncture + exotic shell).

**Paradigm note:** this is not "geometry-first" in the ill-posed inversion
sense. It is *forward* search over an enlarged cause space (matter **+**
geometric free data). Every candidate remains a constraint-solved Einstein
solution, so the entire existing gate/score/validation ladder applies.

---

## 2. What already exists (verified in-tree)

| Piece | Location | Status |
|-------|----------|--------|
| Bowen–York ψ and A_ij for two punctures (mass, spin, momentum, offset) | `../GRTresna/Source/Core/PsiAndAijFunctions.cpp` (`compute_bowenyork_psi`, `compute_bowenyork_Aij`) | Built, used by ScalarFieldBH example |
| Constraint methods handle puncture seeds | `../GRTresna/Source/Methods/CTTK.impl.hpp`, `CTTKHybrid.impl.hpp` | Built |
| Wrapper config fields `bh1_/bh2_{bare_mass,spin,momentum,offset}` | `src/grteclyn_wrapper/grtresna/solver/config.py` (lines ~41–48) | Exposed, defaults benign |
| params.txt emission of all 8 puncture params | `src/grteclyn_wrapper/grtresna/solver/params.py` (lines ~198–205) | Wired |
| Search-space plug point (`grtresna_ansatz` dispatch) | `src/grteclyn_wrapper/search/optimize/spaces.py` (`build_search_space`) | Add a `"puncture"` branch |
| Pre-evolution gates (GRTresna reject, solved-FTL, postload) | `search/grtresna_evaluation_gates.py`, `search/solved_ftl_gate.py`, `projection/postload_gate.py` | Reused as-is (thresholds may need puncture variants) |
| Graded trapped-surface penalty | `src/grteclyn_wrapper/metrics/score/horizon.py` | **Must be regionalized** (§Phase 1) — punctures contain horizons by construction |
| `psi_floor`, `maximal_slicing`, exotic-safe solver | `grtresna/solver/` | Interplay with punctures to be validated in Phase 0 |

**Not yet verified (Phase 0 targets):** `.gridinit` ψ→χ conversion near the
puncture (χ→0), GRTeclyn RadialRecipe moving-puncture gauge settings (1+log
lapse + Γ-driver shift with advection, η parameter), AMR tagging around
puncture locations, probe-ray capture accounting.

---

## 3. Phase 0 — Feasibility smokes (no search, ~days)

Goal: prove one puncture survives the full pipeline. All runs single-GPU,
coarse, frames on.

| # | Task | Acceptance |
|---|------|-----------|
| 0.1 | **GRTresna vacuum puncture solve:** `bh1_bare_mass=1.0`, all lumps zero, solve on 128³/L=64 with `regrid_radius` centered on the puncture | Converged solve; Ham/Mom within existing thresholds away from the puncture; ψ profile matches Brill–Lindquist analytic to a few % at r > 2m |
| 0.2 | **gridinit conversion:** inspect χ = ψ⁻⁴ near puncture; confirm no NaN/negative χ; check `psi_floor` behavior | Loadable `.gridinit`; postload gate passes with a puncture-specific threshold if needed |
| 0.3 | **GPU evolution smoke (t=16):** verify/enable moving-puncture gauge in RadialRecipe (1+log lapse, Γ-driver shift + advection, η≈1/m); AMR tags follow the puncture | No NaN to t=16; lapse collapses at puncture (expected); constraint norms bounded away from puncture |
| 0.4 | **Probe calibration (ties to NextSteps P1/P2):** run the 4D evolving geodesic probe on the puncture run with rays *offset* from the hole | Probe reports Shapiro **delay** (negative unclipped f_geo) of plausible magnitude; rays that fall into the hole are flagged `captured`, excluded from f_geo, counted in the report |
| 0.5 | **Two punctures:** `bh1/bh2` offset ±4–6 on z, small masses | Solve + short evolution clean; head-on infall visible in frames |

Deliverable: `scripts/campaigns/puncture_ftl/smoke.sh` running 0.1–0.5;
findings appended to this file.

**Probe change required (0.4):** `metrics/probes/ftl/evolving_geodesic.py` /
`geodesic.py` — add capture detection (ray enters χ < χ_min region or proper
step collapses) and a `n_captured` field on the reports. Without this, a
puncture on the ray path silently degrades `n_reached` and poisons trust
statistics.

---

## 4. Phase 1 — Scoring & gate adaptation (~days, parallel with Phase 0)

Punctures break three scoring assumptions. Fix them **before** any search or
the optimizer will exploit the mismatch (see gw_beam v3 lesson).

### 1.1 Regionalize the trapped-surface veto

`horizon.py` grades a single centered horizon proxy; a puncture candidate has a
horizon at t=0 *by construction* and would be floored instantly.

- Add an **expected-horizon mask**: spheres of radius `k·m_i` (k≈1.5) around
  each puncture's (advected) location, sourced from the genome via
  `matter.json`-style metadata.
- Penalty applies only to trapped surfaces **outside** the mask or
  **intersecting the probe channel** (the ray-fan corridor already known to the
  probe). A horizon swallowing the channel = full graded penalty; a horizon
  sitting quietly at the puncture = none.
- Files: `metrics/score/horizon.py`, `metrics/diagnostics/collapse.py` (needs
  per-region θ₊ minima, not just the global minimum), objective hookups in
  `metrics/score/objectives.py`.

### 1.2 Persistence for zero-matter genomes

`structural_persistence = density_retention × morphological_coherence` is
matter-defined and degenerates to junk on vacuum candidates.

- Add a **geometry-persistence** alternative for genomes with zero matter
  amplitude: retention of the χ-well depth / throat-channel profile (reuse the
  solved-slice channel metrics) instead of ρ retention.
- Files: `metrics/score/survival.py` (or equivalent), gated by a
  `has_matter` flag from the genome.

### 1.3 Gate thresholds

- Postload `L2_Ham/Mom ≤ 1e-2` will be dominated by the puncture neighborhood.
  Either mask the expected-horizon spheres out of the postload norm or add a
  puncture-specific threshold env (`POSTLOAD_HAM_MAX_PUNCTURE`).
- Solved-FTL gate: confirm it doesn't reject puncture slices as "degenerate"
  (deep χ wells are exactly what it should keep).

Acceptance: a hand-built puncture candidate passes gates and scores finitely;
a hand-built "horizon eats the channel" candidate is penalized.

---

## 5. Phase 2 — Genome + campaign wiring (~1 week)

### 2.1 Search space

New `grtresna_puncture_search_space()` in `search/optimize/spaces.py` and a
`grtresna_ansatz == "puncture"` branch in `build_search_space`. Production
space pinned to **12-D** (axis conventions matching the probe: channel along z):

| Dim | Bounds | Notes |
|-----|--------|-------|
| `grtresna_bh1_bare_mass` | 0.3 – 2.5 | horizon r≈m/2; keep ≥ ~4 finest cells (dx_fine=0.125–0.25 at ml=2) |
| `grtresna_bh1_offset_z` | −8 – 8 | throat placement on the probe axis |
| `grtresna_bh1_spin_z` | −0.5 – 0.5 | in units of m² (enforce |J| ≤ 0.5 m² in wiring) |
| `grtresna_bh1_momentum_x/z` | −0.4 – 0.4 | units of m; transverse + axial boosts |
| `grtresna_bh2_bare_mass` | 0.0 – 2.5 | 0 disables the second puncture (single-throat subspace) |
| `grtresna_bh2_offset_z` | −8 – 8 | |
| `grtresna_bh2_spin_z` | −0.5 – 0.5 | |
| `grtresna_bh2_momentum_x/z` | −0.4 – 0.4 | |
| `grtresna_bh_separation_min` | — | not a dim: wiring-level reject if |offset₁−offset₂| < m₁+m₂ (merged punctures = one hole, wasted eval) |

Wiring: map `grtresna_bh*` param keys → `GRTresnaConfig` fields alongside the
existing `grtresna_lump{k}_*` mapping (same module that builds the config from
the params dict); emit puncture metadata into the eval's matter/geometry JSON
for the horizon mask (§1.1).

### 2.2 Campaign launcher

`scripts/campaigns/puncture_ftl/run.sh` modeled on
`scripts/campaigns/general_ftl/`:

- `GRTRESNA_ANSATZ=puncture`, matter lumps pinned to zero (Phase 2 is
  geometry-only), `OBJECTIVE_MODE=general_ftl`, descriptor `ftl_lifetime`,
  `GRTECLYN_GEO_DIRECTIONS="x y z"`, `GRTECLYN_EVOLVING_GEODESIC=1`.
- `maximal_slicing` forced on (punctures + K=0 is the standard, and the exotic
  wedge in Phase 3 requires it anyway).
- Search defaults from `search_common.sh` unchanged (128³, L=64, t=16) but
  `max_level=2` minimum and regrid centered on punctures.

### 2.3 Expected outcome (calibration, not discovery)

Vacuum punctures should produce **no FTL** (Shapiro delay everywhere) and a
healthy BBH-like dynamics distribution. This is the negative-control campaign
at scale: if *any* vacuum candidate shows f_geo > 0, that is a probe bug, not a
discovery — file it against NextSteps P1/P2.

Acceptance: 100-eval QD run completes ≥40% `gpu_ok`; zero vacuum FTL
positives; archive populated across descriptor bins.

---

## 6. Phase 3 — Puncture + exotic shell (the discovery campaign, ~1–2 weeks)

The Morris–Thorne attempt: exotic scalar shells co-located with a puncture
throat, searched for exterior/bridge shortcuts and horizon-delay.

- **Genome:** Phase 2's 12-D puncture block **+** 2–3 lumps restricted to a
  shell parameterization centered on `bh1_offset` (amp, width, radius, exotic
  flag ≈ +10–12 D). Reuse the existing shell ansatz machinery rather than free
  lumps — the physics target is a shell *around the throat*.
- **Solver:** exotic lumps ⇒ `apply_exotic_safe_solver` path (K=0 maximal) —
  already the puncture default, so no mode clash.
- **Objective:** `general_ftl` with the §1.1 channel-aware horizon penalty. The
  interesting gradient is *horizon-delay*: the graded penalty already rewards
  later collapse; verify it composes with the mask (a throat that stays open to
  t=16 with an open channel should dominate the archive).
- **Probe:** ray fans offset laterally from the throat axis at several impact
  parameters (the shortcut, if any, lives in the exterior channel near the
  shell, not down the hole). Add impact-parameter sweep to the probe options if
  not already covered by the z-fan.
- **HQ promotion:** unchanged ladder (256³, t=30, hq probe). Key comparison
  row: wormhole v22 eval 046 (throat-like channel, no puncture, horizon killed
  it at t≈21) vs the best puncture+shell candidate — same table, same metrics.

Acceptance: campaign completes with gates holding (no reward-hack signatures);
either (a) a confirmed transient exterior shortcut near a real bridge — new
result class, or (b) a clean negative — also publishable against the v22
baseline ("exotic shells around genuine throats do not outperform matter-only
channels").

---

## 7. Phase 4 — Extensions (opportunistic)

| Item | What | Cost |
|------|------|------|
| Multi-throat | bh1+bh2 both massive with independent shells; descriptor axis = throat count / separation | genome + descriptor tweak |
| Spin ladders | frame-dragging assist on the channel (spin_z sweep at fixed shell) | pinned-subspace CMA-ES |
| Boosted throat | `bh1_momentum_z` ≠ 0 — a *moving* throat vs the static v22 wormholes | included in Phase 2 dims |
| ADM bookkeeping | ADM mass/momentum per candidate (NextSteps §4.4) — punctures make this cheap to validate against m₁+m₂ | small; do during Phase 0 |

---

## 8. Risks & open questions

| Risk | Mitigation |
|------|-----------|
| GRTresna nonlinear solve diverges with exotic shell *at* the puncture (ρ<0 in a strong-ψ region) | Phase 0.1/0.2 sweep amp×radius before the campaign; reject-and-log in gates (counts toward funnel stats) |
| RadialRecipe gauge not puncture-ready (no Γ-driver advection / wrong η) | Phase 0.3 explicitly verifies; if missing, it's a small C++ params change, budget 1–2 days rebuild |
| Under-resolved punctures at search res (m=0.3 ⇒ r_h≈0.15 ≈ 1 fine cell) | mass lower bound tied to dx_fine in the space definition; wiring-level reject below 3 cells |
| Probe rays captured → biased f_geo | capture accounting (Phase 0.4) **before** any scored run |
| Optimizer exploits the horizon mask (hides a channel-eating horizon inside a mask sphere by inflating m) | mask radius from *genome* mass, penalty still applies to any trapped surface intersecting the channel regardless of mask; add a mask-abuse check to the validation tiers |
| Vacuum genomes break matter-based metrics silently | §1.2 `has_matter` gating; smoke test includes a vacuum candidate through full scoring |
| Frames/consumer divide-by-χ blowups near puncture | clamp χ in visualization only (never in physics) |

---

## 9. Execution order & effort

```
Phase 0  smokes + probe capture accounting        ~3–5 days   (1 GPU)
Phase 1  horizon mask + persistence + gates       ~3–4 days   (parallel, CPU)
Phase 2  genome, launcher, vacuum QD (100 evals)  ~1 week     (8 GPUs, calibration)
Phase 3  exotic-shell discovery campaign          ~1–2 weeks  (8 GPUs)
Phase 4  extensions                               opportunistic
```

Dependencies: Phase 2 requires 0+1 complete. Phase 3 additionally benefits from
NextSteps **P1/P2** (probe controls + gauge-honest baseline) being done first —
a new geometry family measured with an uncalibrated probe repeats the old
mistake at higher stakes.

Commit discipline: land Phases 0–1 as separate reviewed commits (probe capture
accounting and horizon regionalization change *existing* scoring behavior —
they need regression tests on the current champions to prove old campaigns'
scores are unchanged where no puncture metadata is present).
