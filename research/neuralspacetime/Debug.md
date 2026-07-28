# Pump controller / constraint investigation — findings log

Status as of 2026-07-28. Branch `feature/interstellar`.
Commits: `fe5ef9f8`, `d6a0c350`, `2ddaa87a`, `6daea1a3`, `fbdbd740`.

**READ §15 AND §16 FIRST.** §15: the metric cache feeding the geodesic probes
was sampling at the coarsest AMR level and silently erasing refined geometry —
every ladder `f_geo_evol` number is void. §16: the full audit of every other
metric in the paper and the campaign that §15 forced; nine more defects found
(none voiding the paper's headline — the RM stack passes the §15 fidelity
check and the 20.15% trace reproduces bit-exactly), several paper-level
corrections required. The simulations themselves are unaffected — constraints
(§4), retention (§5) and the tp30 instability (§6) all stand.

The 2026-07-28 HQ ladder is complete and validated (§11); its data now lives at
`runs/pump_ladder_m0_v1`. §12 is a separate retracted claim.

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

### HQ CONFIRMATION (256³, `runs/pump_ladder_m0`) — the claim extends to `t_pump ≤ 24`

| t_pump | peak L2_Ham | mean L2_Ham | final | governor min |
|--------|-------------|-------------|-------|--------------|
| 0 (none) | 4.2723e-3 | 2.7116e-3 | 3.7868e-3 | 1.0000 |
| 4 | **3.9018e-3** | 2.8828e-3 | 2.9172e-3 | 1.0000 |
| 8 | 4.1969e-3 | 3.0880e-3 | 3.9745e-3 | 1.0000 |
| 16 | 4.6773e-3 | 3.1984e-3 | 4.6773e-3 | 1.0000 |
| 24 | 4.1474e-3 | 3.4027e-3 | 3.9829e-3 | 1.0000 |
| 30 | *4.7606e-2* | *5.0860e-3* | *4.7606e-2* | **0.0002** |

* peak is again **non-monotonic**: the UNPUMPED run (4.27e-3) sits *between*
  tp4 (3.90e-3) and tp16 (4.68e-3) — the same shape as the fast tier
* mean 2.71e-3 → 3.40e-3 over **24** units of continuous driving: 26% spread
* governor never engaged for any `t_pump ≤ 24`
* **the hard boundary moves from 16 to 24.** tp24 is completely healthy; only
  tp30 destabilises (§6)

Duhamel bound at HQ (all tight, all satisfied):

| run | final L2_Ham | ham_bound | ratio |
|-----|--------------|-----------|-------|
| tp4 | 2.917e-3 | 3.594e-3 | 1.23× |
| tp8 | 3.974e-3 | 4.189e-3 | 1.05× |
| tp16 | 4.677e-3 | 5.235e-3 | 1.12× |
| tp24 | 3.983e-3 | 6.892e-3 | 1.73× |

**Claim, for `t_pump ≤ 16` (fast tier) / `≤ 24` (HQ):** constraint violation is
statistically indistinguishable from the pump-free run.
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

**HARD BOUNDARY.** The claim holds to `t_pump = 24` at HQ (16 at the fast
tier), NOT to 30. tp30 reaches 4.76e-2 (11× baseline) and destabilises. Do not
over-claim.

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

### HQ CONFIRMATION — larger gain, still monotonic, still no saturation

`confined_frac`, 256³:

| t | tp0 | tp4 | tp8 | tp16 | tp24 | tp30 |
|---|-----|-----|-----|------|------|------|
| 25.92 | 0.0852 | 0.1945 | 0.2286 | 0.2382 | 0.3075 | 0.2866 |
| 28.80 | 0.0225 | 0.1312 | 0.1699 | 0.2161 | 0.2191 | 0.2462 |
| 30.02 | **0.0157** | 0.1111 | 0.1431 | 0.1903 | **0.2042** | *0.2390* |

* at t=30.02: **1.6% (no pump) → 20.4% (t_pump=24)** — a **13.0× gain**
  (fast tier gave 10.7× at tp16)
* monotonic in pump duration at every late epoch, `t_pump` 0→24, no saturation
* tp30's 0.2390 is inside the destabilised window (§6) — quote tp24, not tp30
* the early-time harm reproduces: below t≈14 every pumped rung sits *below* the
  pump-free run (ratio 0.89–0.99)

### Lapse health is also monotonic — until it isn't
Fast tier, `min_lapse` at t=30: **0.074 → 0.112 → 0.228 → 0.421 → floor(1e-10)**
for t_pump 0/4/8/16/30.

HQ at t=30, `min_lapse` / `min_chi` for t_pump 0/4/8/16/24/30:

| t_pump | min_lapse | min_chi | max_abs_K | max_Pi |
|--------|-----------|---------|-----------|--------|
| 0 | 2.607e-2 | 5.682e-4 | 4.916e-1 | 6.61e-2 |
| 4 | 1.111e-1 | 5.956e-2 | 2.746e-1 | 5.13e-2 |
| 8 | 2.254e-1 | 1.720e-1 | 2.491e-1 | 8.08e-2 |
| 16 | 4.236e-1 | 4.026e-1 | 2.090e-1 | 1.79e-2 |
| 24 | **6.282e-1** | **5.826e-1** | 7.648e-2 | 3.36e-2 |
| 30 | *4.953e-5* | *2.914e-1* | 3.853e-1 | *5.24* |

Monotonic all the way to **tp24** at HQ (not 16), then tp30 collapses. Note the
pump-free run has `min_chi` = 5.7e-4 — its central region compacts hard. That
is a real, large effect; see §12 for what it does NOT license you to say.

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

**ANSWERED AT HQ — it is NOT a `dx=1.0` artifact.** `runs/pump_ladder_m0/lad_m0_tp30`
reproduces the failure at 256³:

| quantity | fast tier (t=30) | HQ (t=30) |
|----------|------------------|-----------|
| min_lapse | 1.00e-10 | **4.95e-5** |
| min_chi | 0.337 | **0.291** |
| max_Pi | 1.44e+4 | 5.24 |
| peak L2_Ham | 1.96e-2 | 4.76e-2 |
| governor min | 0.0000 | 0.0002 |

Same signature: `min_lapse` collapses five orders while `min_chi` barely moves.
Lapse collapse with a healthy conformal factor — a gauge / driving pathology,
at both resolutions. The governor engaged legitimately in both.

**The turnover is between 24 and 30, not 16 and 30.** `tp24` is completely
healthy at HQ (governor 1.0000, min_lapse 0.628, min_chi 0.583, peak L2_Ham
4.15e-3 — *below* the pump-free run). This ladder cannot localise the turnover
further: `tp24` and `tp30` are bit-identical up to t=24 by construction, so a
finer bracket needs new rungs at t_pump = 26/28.

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

### HQ CONFIRMATION — no dose–response, confirmed at 256³

Peak `f_geo` per rung (all 22/22 rows `geo_trustworthy=1`, better than the
fast tier where several epochs were 0):

| t_pump | peak f_geo | at t |
|--------|-----------|------|
| 0 (none) | 27.0% | 30.00 |
| 4 | 29.5% | 10.08 |
| 8 | **30.9%** | 14.40 |
| 16 | 23.5% | 12.96 |
| 24 | 23.5% | 12.96 |
| 30 | 23.5% | 12.96 |

Identical verdict to the fast tier: tp8 is best, and the PUMP-FREE run (27.0%)
beats tp16/24/30 (23.5%). **The pump does not improve `f_geo`, at either
resolution.** tp16/24/30 tie exactly because their peak falls at t=12.96, before
any of them differ.

**GAP — NOW CLOSED. See §13.** `f_geo_evol` was never computed for the fast
tier *or* the HQ ladder, because no campaign queue runs the scoring layer at
all. A post-hoc pass now computes it from the `metric_stack` cache; the values
are **~12–13% and trustworthy**.

CORRECTION to the sentence previously here: it cited `research.tex` line 172 as
saying `f_geo^evol = 0` "for every accepted evaluation anyway", dropping the
qualifier. Line 172 is about the **canonical-only (positive-energy)** search,
not these exotic runs. The paper's claim stands; it was never in tension with
this measurement.

---

## 8. Data provenance

**VALID**
* `runs/pump_ladder_m0_v1/lad_m0_tp{0,4,8,16,24}` — HQ, 256³, t=30, complete,
  validated (§11). Manuscript numbers for constraints / retention / lapse.
  **Its `metric_stack` was deleted and `f_geo_evol` reset — see §15.**
* `runs/pump_ladder_m0_v1/lad_m0_tp30` — valid to t≈29; governor engaged
  legitimately at the end, late-time values are inside the destabilised window
* `runs/pump_ladder_m0/` — RERUN IN PROGRESS (launched 04:46), 257³ cache
* `runs/always_on_pump/hq146_m0_tp4_t30` — HQ, mode 0, governor 1.0 throughout
  (bit-identity reference for `lad_m0_tp4`)
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

## 9. HQ ladder — COMPLETE

`runs/pump_ladder_m0/` — 6 rungs (`t_pump` = 0/4/8/16/24/30), mode 0, L=128
N=256 ml=3, t=30, `plot_interval=144`. Launched 01:40, finished 03:29
(2026-07-28). **All six ran 6-way concurrent on one node** — possible only
because plotfiles now go to node-local NVMe (wrapper README).

Queue: `grteclyn-wrapper/scripts/campaigns/rl/pump_ladder_queue.py`.

Its three purposes, all answered:
* (a) does the tp30 lapse collapse survive at HQ → **yes**, §6
* (b) where does driving turn destructive → **between 24 and 30**, not 16 and 30
* (c) HQ f_geo + retention → §5, §7

Throughput, measured against the surviving NFS-era campaign
(`runs/always_on_pump/queue.log`, max 2 concurrent, 7 runs, 5 h 40 min):

| | per-run rate | extraction mean / max | 6 runs to t=30 |
|---|---|---|---|
| NFS plotfiles, 2 concurrent | 0.333 t/min | 74.0 s / **310.3 s** | ~4 h 33 min |
| node-local, 6 concurrent | 0.335–0.378 t/min | 15.7 s / 21.9 s | **1 h 30 min** |

**3.2× on campaign wall-clock, and per-run speed is unchanged** — 3× the jobs
at the same cost each, which is the proof the I/O wall is gone. The old worst
case (310 s) exceeded the 288 s plotfile cadence; that is how backlogs formed.

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
* that `f_geo^evol` is unmeasured or zero **for the exotic pump runs** — it is
  ~12–13% and trustworthy (§13). This leaves `research.tex` line 172 intact:
  that line is about the *canonical-only positive-energy* search, a different
  sector. Do not conflate them.
* that a black hole / apparent horizon forms in any run (§12)

### Updated for HQ

**CAN now say**
* the constraint claim holds to `t_pump ≤ 24` at 256³, not just 16
* retention gain is **13.0×** at t_pump=24 (1.6% → 20.4% at t=30.02)
* the tp30 lapse-collapse instability is resolution-independent (§6)
* the governor fix is a verified no-op in mode 0 — `lad_m0_tp4` is
  bit-identical to `hq146_m0_tp4_t30` on every physics column

**STILL CANNOT say**
* anything about `t_pump > 24`
* **anything at all about `f_geo_evol`** — every value is void (§15)
* that the pump does not improve `f_geo` — the §7 verdict is *plausible but
  unproven*; that probe is a milder case of the same resample problem and must
  be audited against the rerun (§15)
* where between 24 and 30 the turnover sits — needs rungs at 26/28

---

## 11. Validation record — HQ ladder

| check | result |
|-------|--------|
| completion | 6/6 reached t=30.00, rc=0, 3000 rows each |
| NaN / Inf | none (only expected `structure_coherence=nan`) |
| **bit-identity** `lad_m0_tp4` vs `hq146_m0_tp4_t30` | **identical, cols 1–8** |
| determinism tp16/tp24/tp30 for t ≤ 16 | exact at t=8.00 / 12.00 / 16.00 |
| governor | 1.0000 for tp0–tp24; 0.0002 for tp30 (legitimate) |
| Duhamel bound | satisfied, 1.05–1.73× |
| plotfile GC | 0 deleted-before-extracted across all six; no empty skeletons |

The bit-identity result is the important one: it proves `d6a0c350` (the
governor fix) changes nothing in `controller_reservoir_mode = 0`, so every
mode-0 number from the first campaign remains valid.

---

## 12. RETRACTED — "the pump delays black hole formation"

**A claim was made from `collapse_diagnostics.dat` cols 8/9 and then withdrawn
after checking. Do not re-derive it.**

The claim was that horizons form at t = 15.82 (pump-free) → 24.50 (tp16), i.e.
the pump delays collapse by ~9 time units. It is not supported.

`min_theta_plus` (col 9) and `max_ah_r` (col 8) are **pointwise proxies**, not
a surface-integrated expansion ([`RadialRecipeLevel.cpp:838`](../../Examples/RadialRecipe/RadialRecipeLevel.cpp)):

```
theta_plus = 2*sqrt(chi)/r - dchi_dr/sqrt(chi) + A_rr - (2/3)*K
```

`max_ah_r` is a `ReduceRealMax` of `r` over any cell with `theta_plus <= 0` —
the OUTERMOST such cell, not a horizon radius. Four checks kill the claim:

* `r_at_min_theta_plus` stays pinned at **8.2–8.6** and never migrates inward;
  a forming horizon moves inward
* `min_theta_plus` tracks `max_abs_K` almost exactly (tp0: 0.026 → −1.70 as
  K goes 0.12 → 0.91)
* at r ≈ 8.3 with chi ≈ 1, `2*sqrt(chi)/r ≈ 0.24`, so `K >~ 0.36` alone flips
  the sign — no trapped surface required
* `tp24` later returns to `min_theta_plus = +0.041`. **A horizon does not
  un-form.**

This is the `-(2/3)K` term in a gauge where K grows, nothing more.

**What survives:** `min_chi` at t=30 (§5) — 5.7e-4 pump-free vs 0.583 at
tp24 — is a real and large difference, monotonic to tp24. Central compaction
is genuine. Calling it collapse requires a real apparent-horizon finder, which
this codebase does not have.

---

## 13. `f_geo_evol` — computed at last, and it is NOT zero

### Why it was never computed
`consume_plotfiles --evolving-geodesic` does **not** run the 4D trace. Its only
effects are setting `GRTECLYN_EVOLVING_GEODESIC=1` inside the consumer's own
process (where nothing reads it) and enabling the `small_data/metric_stack`
cache at the right resolution. The trace lives in
`metrics.aggregation.collector`, reached only via `core/evaluation.py` →
`runner.py` — the QD/CMA-ES evaluation path. **A hand-rolled campaign queue
bypasses it entirely**; `pump_ladder_queue.py` has no post-run step at all.

Worse, `consume_plotfiles/extraction/ftl.py` writes literal `0.0  0` into cols
13/14 on every row, so the columns look *computed-and-zero* rather than
*never-computed*. That is how this went unnoticed.

### The fix
`grteclyn-wrapper/scripts/campaigns/rl/score_evolving_geodesic.py` — a post-hoc
pass reading the `metric_stack` cache (on NFS, survives plotfile pruning), so
it can be run any time after a campaign finishes:

```
grteclyn-wrapper/.venv/bin/python \
  grteclyn-wrapper/scripts/campaigns/rl/score_evolving_geodesic.py \
  runs/pump_ladder_m0/lad_m0_tp*
```

### Results — **ALL RETRACTED, see §15**

The numbers this pass produced (12–13%, and a later emission sweep giving up to
26.4%) were computed from a cache that could not represent the geometry. They
are deleted, not corrected. **§15 is the section that matters.** The machinery
below is sound; the input it was fed was not.

### Three consequences — (1) stands, (2) and (3) are void

Only the first survives §15. The other two were readings of retracted numbers.

1. **This does NOT contradict `research.tex` line 172.** That line reports
   `f_geo^evol = 0` for a *canonical-only (positive-energy) search*
   (`sec:canonical_bound`) — a different matter sector from these exotic pump
   runs. §7 of this file previously quoted it without the "canonical-only"
   qualifier; that was a misreading. The two numbers do not compete.

   What IS worth knowing: `_post_pump_emit_ok` (`metrics/score/ftl.py`) rejects
   any `t_emit` below `GEODESIC_EMIT_MIN_TIME` / `RL_PUMP_STOP_TIME`. Every
   trace here emits at `t_emit = 0.00`, so **in a scoring context these six
   values would be rejected and recorded as 0** — the raw measurement is real,
   but the scorer would not credit it. That is a property of this pass's
   configuration, not a defect in the paper.
2. **No dose–response here either.** tp4 is best (13.22%); pump-free (12.26%)
   beats tp16/24/30. Same verdict as frozen `f_geo`.
3. **This configuration cannot resolve pump duration.** tp16/24/30 are
   bit-identical (0.12036) because a ray emitted at t=0 arrives before t=16 and
   never samples the spacetime where those rungs differ. A meaningful test needs
   an emission sweep with `GEODESIC_EMIT_MIN_TIME` set past each rung's stop time.

### Caveat
The cache is **33³** (search-mode `GRTECLYN_METRIC_STACK_N_SPACE`), fixed at
write time. It cannot be raised after the fact — the plotfiles are gone. Future
campaigns wanting HQ-resolution 4D traces must set
`GRTECLYN_EVOLVING_GEODESIC_MODE=hq` **before launch**.

---

## 14. Bug fixed while validating — five disagreeing trust bars (`6daea1a3`)

Both geodesic probes treat a captured ray (fell into a puncture throat /
horizon) as physics, not an integration failure:
`n_reached == n_rays - n_captured`. Every *consumer* of their reports still
used `n_reached == n_rays`. Any run where a ray was captured was certified by
the probe and simultaneously marked untrusted downstream.

| site | gated | was |
|------|-------|-----|
| `metrics/score/ftl.py:204` | `geo_trustworthy` (frozen) | stale |
| `metrics/score/ftl.py:256` | `evo_trustworthy` (4D) | stale |
| `metrics/aggregation/collector.py:279` | `f_geo_evol_ok`, col 14 | stale |
| `search/ftl_peak_metrics.py:54` | QD/MAP-Elites retention gate | stale |
| `consume_plotfiles/extraction/ftl.py:82` | `geo_trustworthy`, col 4 | stale |

The scoring ones had teeth: they zeroed `ftl_geo_evolving` /
`operational_ftl_geodesic`, so a certified shortcut earned no FTL credit.

Fixed: `geodesic.rays_complete(n_rays, n_reached, n_captured)` is now the single
definition and all five call it. `EvolvingGeodesicMetrics` gained `n_captured`,
which it never carried — so the scorer and the QD gate *could not* have applied
the correct bar even in principle. Regression test:
`tests/metrics/ftl/test_ray_bundle_trust_bar.py`.

**Behaviour is unchanged whenever `n_captured == 0`, which covers every run
recorded to date** — including this ladder (0 captured in all six). No existing
number moves.

---

## 15. THE BIG ONE — the metric cache silently erased the geometry (`fbdbd740`)

**Every geodesic number derived from `metric_stack` is void. The simulations
themselves are fine. A rerun is in progress.**

### How it was caught
Not by a test. By a physics objection from the user:

> *"cause for the pump free matter collapsing and light rays should have more
> trouble reaching detector compare to pump runs with more stable matter"*

That is exactly right, and it is the check that should have been applied
before reporting anything. A deep gravitational well produces **Shapiro
delay** — light through it arrives *later*, not earlier. The pump-free run has
the deepest well by three orders of magnitude (`min_chi` = 5.7e-4 vs 0.58 for
tp24) and was reporting the **largest** apparent shortcut. That is backwards,
and being backwards it falsified the measurement, not the physics.

### The bug
The evolving-geodesic probe does not read plotfiles. It reads
`small_data/metric_stack/`, a cache of `g_{mu nu}` resampled onto a **uniform**
grid — `n_space³` points over `±half_width`. Defaults: `n_space` = 33 in search
mode, 65 in hq mode; `half_width` = `ftl_L` = 8.

So the cache spacing is `dx = 2*8/32 = 0.5`.

The runs use `amr.n_cell = 256` over `L = 128` (level-0 `dx` = 0.5) with
`amr.max_level = 3`, i.e. finest `dx` = **0.0625**.

**The cache sampled at the coarsest level and discarded all three levels of
refinement.** There is no error, no warning, and no visible artifact: the
resampled metric is smooth, positive-definite, and null geodesics integrate
through it perfectly happily. It is simply a different spacetime — one with the
sharp features removed.

### Why it corrupted a comparison rather than shifting it

The cache is faithful while the geometry is smooth, and fails exactly when a
sharp feature appears. `min_chi` as the simulation reports it vs what the
33³ cache can represent, `lad_m0_tp0`:

| t | true `min_chi` | cache can represent | error |
|---|---------------|---------------------|-------|
| 0.00 | 9.419e-1 | 9.428e-1 | 1.0× |
| 17.28 | 5.300e-1 | 5.453e-1 | 1.0× |
| 23.04 | 1.175e-1 | 1.498e-1 | 1.3× |
| 25.92 | 7.271e-3 | 6.319e-2 | 8.7× |
| 27.36 | 1.790e-3 | 4.200e-2 | 23.5× |
| 28.80 | 1.212e-3 | 8.796e-2 | 72.6× |
| 30.00 | **5.682e-4** | **5.541e-2** | **97.5×** |

And across the ladder at t=30:

| run | true `min_chi` | cached | error |
|-----|---------------|--------|-------|
| **tp0** | 5.682e-4 | 5.622e-2 | **99.0×** |
| tp4 | 5.955e-2 | 6.806e-2 | 1.1× |
| tp8 | 1.720e-1 | 1.717e-1 | 1.0× |
| tp16 | 4.026e-1 | 4.106e-1 | 1.0× |
| tp24 | 5.826e-1 | 5.930e-1 | 1.0× |
| tp30 | 2.914e-1 | 3.690e-1 | 1.3× |

**Only the pump-free run develops a deep sharp well, so only its cache was
wrong.** The pumped runs stay smooth and their caches are accurate to 1.0×.
The error therefore did not shift all six equally — it inflated exactly one,
and that one was the control. Its rays never paid the Shapiro delay, so it
posted the biggest shortcut *because* it was the worst-resolved.

Secondary red flag in the same direction: tp0's cache has `max(g_tt)` = **+1.79**
(positive — shift exceeding lapse), while tp24 stays at −0.40. The tp0 cache is
not a physical metric.

### What it invalidates

* every `f_geo_evol` number in §13 — deleted, not corrected
* the emission sweep that followed (peak 26.4% for tp0 at `t_emit`=26, with the
  pumped runs' ray bundles apparently failing at late launch times). The reading
  drawn from it — *"the pump-free spacetime stays transparent while pumped ones
  go opaque"* — was an artifact of the missing well and is withdrawn.
* **NOT** the frozen `f_geo` (§7). That probe reads plotfiles directly via yt at
  `n=65` over `±8` → `dx` = 0.25, twice as fine, and yt traverses the AMR
  hierarchy. Less exposed — but it is still a uniform resample coarser than
  `dx`=0.0625, so §7's numbers need the same audit before being quoted.
* **NOT** anything that comes from the simulation itself: constraints (§4),
  retention (§5), `min_chi`/`min_lapse` (§5), the tp30 instability (§6). Those
  are written by the C++ at full AMR resolution and never pass through the cache.

### The fix (`fbdbd740`)

1. **`metric_stack_cache.cache_fidelity()`** — compares each cached slice's
   representable `min_chi`, computed exactly as `(det gamma)^(-1/3)`, against
   the value the simulation itself wrote to `collapse_diagnostics.dat` col 3.
   The simulation is the ground truth and it was always sitting right there.
   Verified: flags 5/22 tp0 slices, passes tp4/tp8/tp16/tp24 clean.
2. **`metric_stack_cache.required_n_space()`** — derives the needed resolution
   from `n_cell`, `box_length`, `max_level`, `half_width` instead of a guessed
   default. For this ladder it returns **257**.
3. **`score_evolving_geodesic.py` REFUSES** to report `f_geo_evol` from an
   unfaithful cache (`--force` to override).
4. **`pump_ladder_queue.py`** exports `GRTECLYN_METRIC_STACK_N_SPACE=257` and
   `GRTECLYN_EVOLVING_GEODESIC_MODE=hq` to the consumer.

Cost of 257³: ~650 MB/slice compressed, ~85 GB for a 6-run ladder, against
6.2 TB free. Cheap. The 33³ default was never a considered trade-off — it was
a default nobody had tied to the AMR depth.

### How this affects the paper

**Unaffected — these do not touch the cache:**
* the `E_pump ~ 1.6e-17` withdrawal (§1.1) — a C++ diagnostic bug, independent
* the constraint claim, `t_pump ≤ 24` (§4), and the Duhamel bound
* the retention result, 13.0× at tp24 (§5)
* the tp30 lapse-collapse instability being resolution-independent (§6)
* `research.tex` line 172 (canonical-only positive-energy sector) — untouched

**Must not be written until the rerun lands:**
* any `f_geo_evol` value
* any claim about whether the pump helps or harms 4D geodesic FTL

**Needs an audit before it is quoted:**
* §7's frozen `f_geo` — same class of resample, milder. The verdict there
  ("no dose–response", pump-free beats tp16/24/30) is *plausible but unproven*
  until checked against the well-resolved rerun. Note the frozen probe's
  worst-case exposure is the same run for the same reason.

**The methodological point, which belongs in the paper if any geodesic number
does:** a null-geodesic FTL measurement is only as good as the metric it is
integrated through, and a uniform resample of an AMR grid is lossy in a way
that is *invisible downstream*. Any reported `f_geo` should be accompanied by
evidence that the sampled metric reproduces the run's own `min_chi`. That check
is now `cache_fidelity()` and costs nothing.

### Data status after this

* `runs/pump_ladder_m0_v1/` — the completed 2026-07-28 ladder, **physics valid**
  (§11 validation stands), metric_stack **deleted**, `f_geo_evol` columns reset
  to placeholder. Retained as the bit-identity reference for the rerun.
* `runs/pump_ladder_m0/` — rerun launched 04:46, identical simulation config,
  only the consumer's cache resolution differs. **It must reproduce
  `constraint_norms.dat` and `collapse_diagnostics.dat` bit-identically** —
  that is the regression check, since nothing touching the sim changed.
* `runs/reservoir_fix_check/` — deleted (abandoned ansatz, t=2 smoke tests).

### The lesson

Two bugs in one session (§14, §15) had the same shape: **a value was derived
from a stale or lossy intermediate while the authoritative source sat unused
next to it.** §14 re-derived a trust bar the probe already computed. §15
resampled a metric whose true `min_chi` the simulation was already writing to
disk every step. In both cases the fix was to compare against the source rather
than trust the copy — and in both cases the check was cheap enough that there
was never a reason not to.

---

## 16. FULL METRIC AUDIT (2026-07-28) — every metric in the paper and this campaign

After §15, every other measurement pipeline was audited before restarting the
campaign: the C++ in-situ diagnostics, every Python probe/extractor the paper
or the ladder quotes, the analysis scripts, and the paper's own surviving data
artifacts (`runs/grtresna_promote/bcma_*`). Verdict first, details after.

### 16.0 Verdict table

| metric | measured by | verdict |
|---|---|---|
| evolving `f_geo` (headline 20.2%) | metric_stack + 4D tracer | **SOUND** — RM stack passes `cache_fidelity`; trace reproduces bit-exactly (16.2). Two probe defects found (frozen tail, overshoot; 16.4) that do NOT touch the t_emit=4 headline |
| frozen `f_geo` (29.5%) | `geodesic.py`, n=65/dx=0.25 | SOUND for candidate-146 (fidelity passes); resample caveat stands for collapsing runs (§15) |
| `f_ff` (7.6%) | `observer_timing.py` | math audited clean (mass-shell projection, tetrad, proper-length D0); reads the same fidelity-passing stacks |
| `f_op` | 2D Dijkstra, y-midplane slice | machinery exact per edge; **2D-slice-only** caveat (16.6) |
| `L2_Ham` / `L2_Mom` / governor | C++, `constraint_norms.dat` | **LEVEL-0 ONLY** (16.5). Within-protocol comparisons stand; absolute values carry the caveat |
| `collapse_diagnostics` reductions | C++ | **finest-level footprint only** + `min/max_Pi` reads 1 of 8 matter components (16.5) |
| `energy_conditions.dat` | C++ | **BUG — wrong potential (μ=0)**, fixed (16.3) |
| `curvature_invariants.dat` | C++ | formulas correct; finest-footprint caveat |
| confinement / dispersion | `extraction/confinement.py` | **3 BUGS, fixed** (16.3): weight dropped `phi`/`Pi`, unit-inconsistent `total`, coordinate-volume fractions |
| `structure_coherence` gate | `general.py` | **BUG — canonical-Re-only weight**, fixed (16.3) |
| `boundary_flux.dat` | `probes/boundary.py` | **BUG — "radial" derivative was ∂φ/∂x**, fixed (16.3) |
| Duhamel bound | `pump_constraint_budget.py` | integrator fine; §4's ratio table mixes categories (16.6) |
| trapped-surface claim | θ+ proxy | RM's signal is more collapse-like than §12's ladder case, but still a pointwise proxy; paper's t≈27 does not match data (16.7) |
| Alcubierre control (32%) | analytic stack | 5-slice ±0.4 stack for a ~7-unit flight → validates shortcut detection, mostly frozen bubble (16.6) |

### 16.1 The paper's stacks PASS the §15 check

`cache_fidelity` (min_chi representability, tol 1.5×) run on all five surviving
candidate-146 stacks — RM, RC, and the three freefall twins (65³, dx=0.25,
`max_level=3`): **every slice passes, at every time, in every run.**
Candidate-146 never develops a feature the probe grid cannot hold; the deep
sharp well that voided the ladder's tp0 (99× error) simply never forms here
(RM `min_chi` at t=30 is 5.96e-2, and the headline flight ends at t=15.5 where
`min_chi` ≈ 0.8). **The §15 bug does not invalidate the paper's headline.**

Correction to §15's *secondary* flag: `max(g_tt) > 0` (≈ +1.5) appears at t=30
in **all** stacks, including every fidelity-passing one. β²>α² late in these
gauges is a real feature, not per-se evidence of cache corruption. The §15
primary evidence (99× min_chi mismatch) is untouched.

### 16.2 The headline number reproduces bit-exactly

Current code (post §14 trust-bar fix), stored RM stack, t_emit=4, x-axis:
`f_geo = 0.201480`, `t_arrival = 15.498689474416453` — **identical to the
stored `evolving_geodesic.json`** and within 0.2pp of the paper's table. The
pipeline that produced the paper is deterministic and still reproducible.

### 16.3 Bugs found and FIXED in this audit

1. **Confinement weight dropped `phi`/`Pi` for bicomplex runs.** The extractor
   fell through to the independent-lump branch, so the canonical field's Re
   component and momentum — ∫|φ|dV=6.9, ∫|Π|dV=15.3 of a ~62 total — never
   entered any dispersion number, and Re/Im of one complex field were summed as
   `|Re|+|Im|` (a stationary Q-ball reads as breathing at ω). Aggregate impact
   measured live at t=10: rms 6.768→6.851, conf_frac 0.6197→0.6122 (~1.2%) —
   BUT canonical-only vs phantom-only differ by 5.4pp (0.588 vs 0.642): **the
   sectors disperse at different rates and the old weight could not see it.**
   Fixed: model-aware `_matter_sectors` (tag > sniff), per-sector columns
   12–17 appended to `confinement.dat`.
2. **`total_activity` changed units row-to-row.** The legacy covering-grid
   branch summed cell values with no dV; the AMR branch integrates w·dV. A run
   momentarily de-refining emitted rows 8× its neighbours (v1 tp0 t=1.44/2.88:
   644.95 vs ~78). Ratios (rms, frac) were unaffected — dV cancels — and
   nothing scored reads `total`. Fixed: both branches now integrate.
3. **`confined_frac` was a coordinate-volume fraction.** The physical mass
   element is w·√γ d³x = w·χ^(-3/2) d³x. At t=30 the ladder compares tp0
   (χ_min=5.7e-4 → χ^(-3/2)=7.4e4) against tp24 (0.58 → 2.2): the control's
   core was weighted ~3e4× less per unit proper volume — same error class as
   §15, biased against the same run. At t=10 (χ_min 0.82) the shift is already
   0.620→0.594. Direction: proper weighting RAISES tp0, so **§5's "13.0×
   retention gain" is an upper bound**; the monotonic ordering is expected to
   survive (matter that has left R_conf is in χ≈1 territory) but must be
   re-measured. Fixed: proper-volume columns (`total/rms/frac_proper`) written
   alongside; scoring still gates on the coordinate value for continuity.
4. **The coherence gate saw 1 of 8 matter components.**
   `matter_coherence_from_plotfile` (→ `structure_coherence` →
   `structural_persistence` → scales EVERY FTL reward in Eq. (score)) used
   `sqrt(phi²+Pi²)` — half the canonical field, none of the phantom sector.
   For candidate-146 (2 canonical / 3 phantom lumps) fragmentation of the
   phantom sector was invisible to the gate. Fixed: shares `_matter_sectors`.
   (Same fallthrough fixed in the `scalar_activity` frame field.)
5. **`energy_conditions.dat` used the wrong potential.** The EC diagnostic
   constructed `BiComplexScalarField(mass, lambda)` — **μ silently defaulted
   to 0** — while the evolution passes `recipe_scalar_mu = 85333`. At
   candidate amplitudes (|φ|~0.08) the sextic term is comparable to the mass
   and quartic terms, so every NEC/WEC/SEC/DEC margin and
   `matter_integral_NEC_violation` ever written is quantitatively wrong (the
   *sign structure* — phantom-sector violation — survives, since the kinetic
   terms dominate the margins' signs). Same omission in the complex-scalar
   branch. Fixed in `RadialRecipeLevel.cpp` (+ rebuild).
6. **`boundary_flux.dat`'s "radial" derivative was ∂φ/∂x.** Sample points were
   displaced by a constant (dr,0,0) instead of along each point's radial unit
   vector — correct at the ±x poles, a transverse derivative elsewhere. Fixed.
   (Remaining caveat: the flux still uses `phi`/`Pi` only.)
7. **Rays completed through FROZEN geometry past the stack's end.**
   `EvolvingMetricField` clamps its time bracket, so a ray still in flight
   after the last slice silently finishes in a static spacetime. Demonstrated
   on the paper's own data: the RM t_emit=12 launch "arrived" at t=32.75
   against a stack ending at 30.0. Fixed: strict by default (such rays are
   integration failures); `allow_frozen_tail=True` retained for analytic
   controls and stationary-stack tests only.
8. **Detector-plane overshoot.** Arrival time was read at the first RK sample
   PAST the detector plane (up to one step late), biasing `f_geo` LOW by up to
   ~0.4pp — conservative, but sloppy. Fixed: linear interpolation back to the
   crossing. NOTE: future `f_geo` values shift up by ≤0.4pp relative to
   published traces; within-campaign comparisons are unaffected.
9. (§15's fixes — 257³ cache, `cache_fidelity` refusal, auto-scoring — carried
   forward; the restarted campaign has all of them from t=0.)

### 16.4 The paper's t_emit=12 sweep points are artifacts

Stored sweeps: at t_emit=12 the bundles are INCOMPLETE — RM 2/5, **RC 0/5**,
RF 4/5 rays — and RM's two "arrivals" ran 2.75 units past the stack end
(defect 7). The caption of the emission-sweep figure claims "Every launch has
5/5 rays reaching and passes the null-Hamiltonian quality gate; at
t_emit=12 the rays arrive no earlier than the flat baseline". **That is false
at t_emit=12** — the plotted zeros there are bundle failures, not measured
zero advantage. t_emit ≤ 10 launches are clean 5/5 everywhere, and the peak
(t_emit=4) is untouched. The manuscript must drop or annotate the t=12 point;
"the window brackets the channel's life" survives via the 20%→12% decay over
t_emit=4→10.

### 16.5 C++ diagnostics: what they actually measure

* **`constraint_norms.dat` is computed on LEVEL 0 ONLY** (`Level()==0`,
  level-0 state, single cell_vol; `RadialRecipeLevel.cpp:437`). So are the
  pump-force norms and the governor's input. Constraint violation living on
  the refined levels (dx 0.25→0.0625) is invisible to all of them. Within-
  protocol comparisons (the entire §4 ladder argument) are unaffected — every
  run is measured the same way — but the absolute numbers are the level-0
  restriction, and the paper's "volume-weighted L2 norms" should say so. The
  Richardson order p≈3.3 is the convergence of that restriction.
* **`collapse_diagnostics` / `energy_conditions` / `curvature_invariants`
  reduce over the FINEST level's boxes only** — which cover 0.01–0.3% of the
  domain mid-run (measured live: level 2 = 0.01% at t=14.4). Spot checks show
  `min_chi`/`min_lapse` DO match the global values (refinement tracks the
  well), but `max_Pi` is already wrong: on tp24 the footprint's canonical-Pi
  range was [-1.4e-2, -3.1e-3] while the global max was +5.4e-4 and the true
  activity peak sat in `Pi_lump0` (+8.4e-2) — **a component the diagnostic
  never reads** (`c_phi`/`c_Pi` only). §6's max_Pi numbers are the canonical
  Re component on a partial footprint: fine as a blow-up flag, not quotable
  as "the field momentum".
* `theta_plus` / `max_ah_r`: §12 stands. See 16.7 for the paper's version.

### 16.6 Caveats that stay caveats (documented, not code-fixed)

* **`f_op` is a 2D y-midplane probe** (yt FRB slice). A 3D channel off the
  plane is invisible. Design choice; state it wherever f_op is quoted.
* **§4's Duhamel table compares the wrong pair.** `ham_bound` bounds the
  PUMP-ATTRIBUTED violation; comparing it against the TOTAL `L2_Ham` (which
  includes the free-evolution residual ~3e-3 present with no pump) mixes
  categories. The meaningful check — pumped-minus-pump-free difference
  (≲1e-3) vs the bound (3.6e-3) — passes with more room, so the conclusion
  survives; the table's framing should change.
* **Resolution-matrix caveat:** the probe grid is FIXED at 65³/dx=0.25 while
  the sim resolution varies (RC/RM/RF), so part of the f_geo ladder agreement
  reflects the common probe grid, not pure sim convergence. Defensible
  because fidelity passes at all three; the paper should state it.
* **Temporal cadence:** the paper's stacks store Δt=0.72 (RM) / 0.24 (RF);
  the ladder rerun stores Δt=2.88 — 4× coarser time interpolation for
  `f_geo_evol`. The 257³ spatial fix did not change the cadence.
* **Alcubierre control:** the analytic stack spans t∈[−0.4,+0.4] for a
  ~7-unit flight — the 32% recovery validates shortcut detection through a
  (mostly) frozen bubble, not the tracer's time interpolation. Tests now
  build honest full-span stacks on the production cache path.
* **Canonical-only bound (research.tex Eq. canonical_bound):** computed
  through 33³ search-mode caches. §15's failure mode (missing wells) biases
  TOWARD spurious shortcuts, so `max f_geo_evol = 0` over 105 evals survives
  that bias; but smoothing can also erase a genuinely narrow corridor, so the
  bound reads "none found within this ansatz/budget/resolution", which is how
  the paper already phrases it. No action needed beyond the resolution note.
* Pre-existing failing test (not from this audit's changes):
  `test_horizon_finder_guard.py::test_charge_retention_columns_are_parsed`
  writes a 22-column wormhole fixture; `read_collapse_metrics` requires ≥23
  for that layout. RotatingWormholeCollapse column drift — reconcile against
  that example's current writer before trusting its j_z/charge fields.

### 16.7 Paper corrections required (for the rewrite; research.tex untouched)

1. §1.1's E_pump withdrawal (already recorded).
2. Emission-sweep figure caption: t_emit=12 points are bundle failures (16.4).
3. **Trapped-surface claim**: the abstract and discussion say "collapses to an
   apparent horizon by t≈27" with "θ+≈−0.9, r_AH≈10.8". In RM's own data
   θ+ first goes negative at **t=18.71**, and unlike the ladder case the
   proxy behaves collapse-like (r_at_min migrates inward 9.6→5.7, no
   recovery, θ+=−0.9 not attributable to the −⅔K term at K=0.27). So the
   signal is *plausible* — but it is still a pointwise radial proxy on the
   finest-level footprint, there is no AH finder, "corroborated" overstates
   it, and t≈27 matches no obvious feature of the data. The headline ray
   (arrival 15.5) precedes even the earlier onset, so the transport
   certificate is safe either way.
4. Confinement numbers ("72%→11%" RM, "→2%" pump-free; spread 3.27): produced
   by the old weight (defect 1) and coordinate volume (defect 3). RM's
   min_chi(t=30)=0.0596 → χ^(-3/2)≈69 at the core, and the pump-free twin
   collapses deeper — the same one-sided bias as the ladder. Dispersal itself
   is not in doubt; the specific percentages are weight- and volume-choice
   dependent and cannot be recomputed (plotfiles gone). Quote from the
   restarted campaign's proper-volume columns instead.
5. State the level-0 norm caveat (16.5) and the fixed-probe-grid caveat
   (16.6) where those numbers appear.

### 16.8 What was checked and found CLEAN

For the record, because "audited" should mean audited: confinement AMR
integral hygiene (Σ cell_volume = V_domain exactly; yt child-masking
correct); code-unit handling and R_conf plumbing (col 9 = 6.668 everywhere);
no radiation-halo contamination of the moments (99.57% of activity inside
r<15; r<32 restriction changes nothing to 6 s.f.); Dijkstra edge speeds solve
the exact per-edge null condition; both probes' flat baselines agree (14.4);
trilinear interpolation, null re-projection (direction-preserving root
choice), and the scale-free h_rel drift measure; `future_null_cov`
future-directed root selection; the freefall probe's mass-shell projection,
Gram–Schmidt tetrad, and proper-length D0 integral; `rays_complete` (§14)
applied consistently at all five consumer sites; Duhamel trapezoid and 16π/8π
factors; `remove_duplicate_time_data` on restart paths.

### 16.9 Campaign restart

The 04:46 rerun was killed at t≈20 (its confinement stream carried defects
1–3 and its consumer predated the fixes; §15's caches were already fine). The
campaign restarts from scratch with: fixed confinement extractor (sector +
proper-volume columns), fixed coherence gate, fixed boundary flux, strict
frozen-tail geodesics, interpolated detector crossing, 257³ metric_stack,
auto-triggered end-of-run scoring, and the rebuilt binary whose only change
is the EC μ fix. The bit-identity check against `runs/pump_ladder_m0_v1`
(physics cols of `constraint_norms.dat` / `collapse_diagnostics.dat`) now
also validates that the rebuild changed diagnostics only.
