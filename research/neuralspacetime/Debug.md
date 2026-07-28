# Pump controller / constraint investigation — findings log

Status as of 2026-07-28. Branch `feature/interstellar`.
Commits: `fe5ef9f8`, `d6a0c350`, `2ddaa87a`.

This file records what was broken, what was fixed, what was abandoned, and
which numbers are safe to quote in the manuscript. Written so the results are
not lost between sessions. **Nothing here has been written into
`research.tex`** — the standing decision is that a validated fix means a full
re-run and a total rewrite, not incremental edits.

---

## 0. One-paragraph summary

The paper's claim that the PD pump is a "soft trajectory guide" with
`E_pump ~ 1.6e-17` is an **artifact of a coordinate-frame bug** in the
diagnostic and must be withdrawn. The real control effort is O(1). The plan to
account for that effort with an explicit energy–momentum reservoir was
implemented, debugged, and then **abandoned** — the reservoir ansatz is
mathematically unstable. It turned out not to be needed: a direct ablation
shows the constraints barely notice the pump. Sustained driving does buy a
large, monotonic gain in **matter retention**, but buys **nothing** in
`f_geo`, and past `t_pump ~ 16` it destabilises the run.

---

## 1. WITHDRAW from the manuscript

### 1.1 `E_pump ~ 1.6e-17` and the "soft trajectory guide" narrative
Sections `sec:pd_pump` (~lines 435–465), plus 1548, 1559, 1570–71.

**Root cause.** The `pump_work` diagnostic evaluated the pump envelope in
ABSOLUTE grid coordinates, while the evolution RHS uses CENTRE-RELATIVE
coordinates (`Coordinates(iv, dx, center)` subtracts the grid centre). The
Gaussian envelope was therefore sampled ~64 code units away from where the
pump actually acts, and underflowed.

**Evidence.** Two binaries, identical config, 200-step replay:
every physics column **bit-identical**; `pump_force_L2` changed
`1.56e-17 -> 1.16e-4` (13 orders of magnitude).
Old diagnostic was nonzero in only 64/3000 rows; the fixed one is nonzero in
**399/399** rows of the driven window.

**Also withdraw:** `pump_work` in `collapse_diagnostics.dat` (col 15) now
means the CORRECTED `∫ f_perp √γ dV`. Same column name, new meaning — never
compare it against runs built before `fe5ef9f8`.

### 1.2 "Confinement saturates at t_pump ≈ 8"
This was produced by the contaminated campaign (see §3) and is false. See §5.

---

## 2. Bugs found and fixed

| # | Bug | Where | Status |
|---|-----|-------|--------|
| 1 | Pump diagnostic used absolute, not centre-relative, coordinates | `RadialRecipeLevel.cpp` | FIXED `fe5ef9f8` |
| 2 | PD law duplicated in 4 places, diagnostic copy had drifted | new `RLPumpForce.hpp` | FIXED `fe5ef9f8` |
| 3 | Reservoir momentum source had the wrong SIGN | `ControllerReservoirMatter.hpp` | FIXED `fe5ef9f8` |
| 4 | Reservoir sources carried a spurious lapse factor | same | FIXED `fe5ef9f8` |
| 5 | Reservoir had NO transport terms (per-cell ODE) | same | FIXED `fe5ef9f8` |
| 6 | `f_i` was never measured at all | `RadialRecipeLevel.cpp` | FIXED `fe5ef9f8` |
| 7 | **Controller ledger fed the pump SAFETY GOVERNOR** | `RadialRecipeConstraintNorms.hpp` | FIXED `d6a0c350` |

### Detail on #3 (sign)
This codebase defines `j_i = Σ_A s_A (−Π_A ∂_i φ_A)`. So the pump's
contribution to matter momentum is `−f_i`, and the reservoir must absorb
`+f_i` — **opposite in sign to the energy case**. Verified by direct
flat-space computation of `∂_t j_i + ∂_j S^j_i` for the scalar EMT.
Symptom of the bug: `L2_Mom` got 2.16× WORSE with the ledger on.

### Detail on #4 (lapse)
`rhs[c_Pi] += S_A` is a coordinate-time rate, so `∂_t ρ|_pump = f_perp`
exactly, with no α. The α in the plan's "α f_⊥" is already absorbed by
`S_A = α J_A`. The extra factor caused under-subtraction wherever α<1.

### Detail on #7 — the one that wrecked a whole campaign
`compute_radial_recipe_constraint_norms` publishes `cached_L2_Ham`, which
drives the pump governor (tanh, centre 0.035, width 0.003). It called
`fill_active_constraints`, which in `controller_reservoir_mode >= 1` INCLUDES
the reservoir. **The controller's own bookkeeping was wired into its own
safety interlock.** When the ledger diverged, the governor throttled the pump
to zero on the strength of a number describing the ledger, not the spacetime.

Governor minima, first campaign:

| config | governor min |
|--------|--------------|
| mode 0 baseline | 1.000 |
| pump-free | 1.000 |
| **every mode-1 run with an active pump** | **0.000** |

Closure times: t=6.90–7.99 (tp8/tp16/tp30), t=7.33 (tp4).
Consequence: tp8/tp16/tp30 came out **bit-identical** because their effective
pump duration was ~7 regardless of the configured stop time.

### NOTE — an error in the original plan's own equations
The plan's densitised energy equation carries a spurious `α K ρ_c` on the RHS;
it fails the FLRW-dust check. Correct form for zero spatial stress:
`∂_t(√γ ρ) + ∂_i[√γ(α j^i − β^i ρ)] = −√γ j^i ∂_i α`.

---

## 3. ABANDONED: the controller energy–momentum reservoir

Implemented in `Source/Matter/ControllerReservoirMatter.hpp` (modes 0/1/2),
fully debugged, then abandoned. **The ansatz is unsound, not the code.**

`S_c^{ij} = 0` makes the reservoir *pressureless dust*. The momentum equation
has no `∇ρ` term, so nothing pushes back; leftover momentum is frozen in and
keeps feeding the energy equation forever. Density grows without bound by
construction — this would happen in exact arithmetic, in flat space, with a
perfect discretisation.

**Decisive evidence.** Pump force is EXACTLY 0 from t=6 (pump off at t=4), yet
the ledger keeps growing `1.9e-2 -> 9.8e-1` by t=30. **Growth with zero input
is self-divergence.**

Head-to-head, mode 0 vs mode 1, identical physics (HQ, t_pump=4):

| t | no ledger | with ledger | ratio |
|---|-----------|-------------|-------|
| 1.0 | 1.28e-3 | 4.54e-4 | 0.36 ✓ |
| 2.0 | 2.13e-3 | 1.18e-3 | 0.55 ✓ |
| 3.0 | 2.72e-3 | 3.18e-3 | 1.17 ✗ |
| 10.0 | 2.32e-3 | 4.42e-2 | 19 ✗ |
| 30.0 | 2.92e-3 | 9.79e-1 | **335** ✗ |

Helps for ~2.5 time units, then diverges. In mode 2 (reservoir gravitates) it
killed a run with a **NaN in K at t=9.69**.

**METHODOLOGICAL WARNING.** The fix was first validated at t=2 and looked
excellent (mode 2: Ham ×0.22, Mom ×0.31). The ledger only turns bad after
t≈2.5. **Never validate this system on short runs.**

What still holds about the reservoir machinery: mode 0 is bit-identical to
pre-reservoir behaviour; the pump-off/mode-2 null test passes bit-identically;
mode 1 does not perturb the evolution (min/max `rho_req` agree to 10 sig figs).

---

## 4. HEADLINE RESULT — the constraints do not care about the pump

This is the defensible claim, and it needs **no reservoir**.

Fast tier (L=128, N=128, ml=2, mode 0), `controller_reservoir_mode = 0`:

| t_pump | peak L2_Ham | mean L2_Ham | governor min |
|--------|-------------|-------------|--------------|
| 0 (none) | **4.4465e-3** | 2.8379e-3 | 1.0000 |
| 4 | 3.8868e-3 | 2.8449e-3 | 1.0000 |
| 8 | 4.1716e-3 | 3.0209e-3 | 1.0000 |
| 16 | 4.6754e-3 | 3.0894e-3 | 1.0000 |
| 30 (t<29) | 1.9579e-2 | 3.4793e-3 | 0.0000 |

**Claim, for `t_pump ≤ 16`:** constraint violation is statistically
indistinguishable from the pump-free run.
* mean varies 2.84e-3 → 3.09e-3 over 16 units of continuous driving: **9% spread**
* peak is **non-monotonic** — the UNPUMPED run (4.45e-3) sits *between*
  tp4 (3.89e-3) and tp16 (4.68e-3)
* a referee cannot claim the pump drives constraint growth when driving for
  4–8 units gives a *lower* peak than not driving at all
* the safety governor never engaged for any `t_pump ≤ 16`

**Independent supporting evidence (HQ, valid runs):**
* pump-free final L2_Ham **3.787e-3** vs pumped-to-t=4 **2.917e-3** — the
  pump-free run is WORSE
* Duhamel bound is tight: measured peak 3.136e-3 vs bound 3.593e-3 (**1.15×**);
  final 2.917e-3 vs 3.594e-3 (1.23×); momentum 7.946e-4 vs 1.748e-3 (2.2×)
* late-time constraint growth is dominated by dispersing matter + regridding
  and is present identically in the pump-free run

**HARD BOUNDARY.** The claim holds to `t_pump = 16`, NOT to 30. tp30 reaches
1.96e-2 (4.4× baseline) and then destabilises. Do not over-claim.

Analysis script: `research/neuralspacetime/analysis/pump_constraint_budget.py`

---

## 5. Pump retention — real, large, monotonic

`confined_frac` (col 5 of `small_data/confinement.dat`), fast tier:

| t | tp0 | tp4 | tp8 | tp16 | tp30 |
|---|-----|-----|-----|------|------|
| 8.64 | 0.655 | 0.616 | 0.601 | 0.599 | 0.599 |
| 17.28 | 0.410 | 0.403 | 0.404 | 0.451 | 0.440 |
| 25.92 | 0.097 | 0.206 | 0.239 | 0.246 | 0.296 |
| 28.80 | 0.028 | 0.144 | 0.184 | 0.228 | 0.261 |
| 30.02 | 0.019 | 0.124 | 0.157 | **0.203** | *0.668 — blowup, discard* |

**Frame this as RETENTION, not confinement growth.** Confinement *decays* in
every run (0.60 → 0.02–0.20). Nothing accumulates matter; the controller only
resists dispersion. What sustained driving buys is a slower leak:

* at t=30.02: **1.9% (no pump) → 20.3% (t_pump=16)** — a **10.7× gain**
* monotonic in pump duration at every late epoch; **no saturation** in 0→16
* early on (t=8.64) the pump marginally HURTS (unpumped 0.655 is highest)

### Lapse health is also monotonic — until it isn't
`min_lapse` at t=30: **0.074 → 0.112 → 0.228 → 0.421 → floor(1e-10)**
for t_pump 0/4/8/16/30. Longer driving keeps the spacetime *healthier*, right
up to tp16, then tp30 collapses.

---

## 6. The tp30 failure — a driven instability, not a black hole

At the fast tier, tp30 dies at t≈29.8. The governor engaged at t=29.10 for a
**legitimate** reason (true L2_Ham crossed 0.035) and did its job; the run
died 0.7 units later anyway.

**It is NOT gravitational collapse:**

| t | min_lapse | min_chi | max_Pi |
|---|-----------|---------|--------|
| 27.8 | 2.17e-1 | 3.94e-1 | 9.4e-2 |
| 29.0 | 4.66e-2 | 3.70e-1 | 1.47 |
| 29.8 | 6.81e-5 | 3.43e-1 | 2.71 |
| 30.0 | 1.00e-10 | **3.37e-1** | 1.44e+4 |

`min_chi` barely moves (0.394 → 0.337) while `min_lapse` falls 5 orders. In
real collapse the conformal factor collapses *with* the lapse. This is lapse
collapse with a healthy metric — a **gauge / driving pathology**. `max_Pi`
stays flat at ~2.4 until the final step, so 1.4e4 is a consequence, not a
cause.

**Mechanism, and it is measurable.** Control effort ESCALATES:

| t | pump force | governor |
|---|-----------|----------|
| 20 | 3.85e-6 | 1.0000 |
| 24 | 6.56e-6 | 1.0000 |
| 28 | **1.11e-4** | 1.0000 |

A **30× rise** from t=20 to t=28. The field drifts further from target, the PD
controller pushes harder, and that drives the lapse collapse. This is exactly
the failure the `RLMatterPumpParams.hpp` comment (lines 55–62) warns about:
the trap fighting gravity and exciting a breathing mode.

**Open question:** whether this survives at HQ resolution, or is a `dx=1.0`
artifact. The HQ ladder answers it. `tp24` brackets the 16→30 turnover.

---

## 7. f_geo — measured, but does NOT respond to the pump

Column 3 of `small_data/ftl_timeseries.dat`. Peak per rung (fast tier):

| t_pump | peak f_geo |
|--------|-----------|
| 0 (none) | 26.9% |
| 4 | 27.5% |
| 8 | **31.0%** |
| 16 | 23.4% |
| 30 | 23.4% |

**No dose–response.** tp8 is best, and the PUMP-FREE run (26.9%) beats both
tp16 and tp30 (23.4%). f_geo also decays to 0 at late times in every run.

**This is awkward for the sustained-control narrative:** what the pump
demonstrably improves is matter retention; the paper's headline metric is the
thing it does not improve. Confront this directly rather than burying it.

Caveats: fast tier only (dx=1.0), and geodesic integration is
resolution-sensitive — quote HQ. `geo_trustworthy` (col 4) is 0 at several
epochs; filter those rows.

**GAP:** `evolving_geodesic.json` is MISSING in all fast-tier runs. The queue
runs the plotfile consumer but not the scoring layer, so `f_geo_evol` was
never computed. Per `research.tex` line 172, `f_geo^evol = 0` for every
accepted evaluation anyway. If it is wanted for these runs, it needs a
separate scoring pass BEFORE the plotfiles are pruned.

---

## 8. Data provenance

**VALID**
* `runs/always_on_pump/hq146_m0_tp4_t30` — HQ, mode 0, governor 1.0 throughout
* `runs/always_on_pump/hq146_m1_tp0_t30` — HQ, pump never on, reservoir ≡ 0
* `runs/pump_ladder_fast/fast_tp{0,4,8,16}` — complete, governor 1.0
* `runs/pump_ladder_fast/fast_tp30` — valid to t≈29, garbage after

**CONTAMINATED — do not use constraint columns (2,3,7,8)**
* `runs/always_on_pump/hq146_m1_tp{4,8,16,30}_t30` — ledger diverged AND
  governor closed at t≈7–8; effective pump duration ~7 in all of them
* `runs/always_on_pump/hq146_m2_tp30_t30` — CRASHED, NaN in K at t=9.69

Their *physics* outputs (confinement, collapse diagnostics, ftl timeseries)
remain valid up to the point the governor closed.

### `constraint_norms.dat` column layout (APPEND-ONLY)
```
1 time  2 L2_Ham  3 L2_Mom  4 min_rho_req  5 max_rho_req  6 integral_neg_rho
7 L2_Ham_rel  8 L2_Mom_rel  9 pump_force_L2  10 governor  11 pump_fi_L2
```
`pump_fi_L2` (col 11) is new — the first real measurement of the momentum
force density. **`f_i` is ~3.5× LARGER than `f_perp`**: the controller's
momentum forcing dominates its energy forcing. Notable because that is exactly
the component the old code never measured and got the sign wrong on.

---

## 9. In flight

`runs/pump_ladder_m0/` — HQ ladder, 6 rungs (`t_pump` = 0/4/8/16/24/30),
mode 0, 2 concurrent, ~5 h total. Queue: `ladder_queue.py`.
`lad_m0_tp4` doubles as a bit-identity regression check against
`runs/always_on_pump/hq146_m0_tp4_t30`.

Purpose: (a) does the tp30 lapse collapse survive at HQ, or was it a `dx=1.0`
artifact; (b) where between 16 and 30 does driving turn destructive; (c) HQ
values of f_geo and retention for the manuscript.

---

## 10. What the paper can and cannot say

**CAN say**
* the pump does not drive constraint growth for `t_pump ≤ 16` (§4), with a
  tight Duhamel bound as backup
* sustained control gives a 10.7× improvement in matter retention (§5)
* the safety governor never engaged in any physically valid run
* control effort is measurable and O(1), not 1e-17

**CANNOT say (yet)**
* anything about `t_pump > 16` — tp30 destabilises
* that the pump improves f_geo — it does not (§7)
* that confinement *grows* — it decays everywhere; the pump slows the decay
* any HQ number until `runs/pump_ladder_m0` completes

**MUST NOT say**
* `E_pump ~ 1.6e-17`, or "soft trajectory guide" (§1.1)
* "confinement saturates at t_pump = 8" (§1.2)
* anything relying on the controller reservoir (§3)
