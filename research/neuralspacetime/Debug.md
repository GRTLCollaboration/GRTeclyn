# Pump controller / confinement — findings log

Status 2026-07-29. Branch `feature/interstellar`.
Commits: `fe5ef9f8` `d6a0c350` `2ddaa87a` `6daea1a3` `fbdbd740` `64c89be4`
`4f31f33a` `e01ec730` `856606c7`.

**Standing rules.**
1. Nothing here has been written into `research.tex`. A validated fix means a
   full re-run and a total rewrite, not incremental edits.
2. A config is "held" only on a run reaching **t ≥ 40**. t=30 proves nothing:
   `pcb_match` looked healthy at 30 and collapsed at ≈32.
3. Collapse is called on central `min_lapse` + `min_chi` + `max|K|`. Never on
   the `theta_plus` proxy. There is no apparent-horizon finder in this repo.

Bracketed tags like `[§15]` are the pre-2026-07-29 section numbers, kept so
older commit messages and notes still resolve. Full-detail prose version:
`git show 856606c7:research/neuralspacetime/Debug.md`.

**Where it stands** [§0]. The pump's original headline — a near-free "soft
trajectory guide" — was a diagnostic bug; real control effort is O(1). The
constraints genuinely do not care about the pump (`t_pump ≤ 24` at HQ), and that
claim needs no reservoir; the reservoir ansatz was abandoned as unstable. The
confinement story then went through two measurement corrections: `confined_frac`
is a ratio with a growing denominator [§19.1], and a params-parsing bug meant
**the phantom sector was never actually pumped** in any run before `19.8`. With
both fixed, the pump does hold matter — the lever is the **target amplitude**,
not gain: aim 0.15 is lethal, 0.10 conserves both sectors with the healthiest
geometry, 0.08 holds then dies at t≈32. **Every configuration so far collapses
eventually; only the timing differs.** Whether a stable window exists between
0.10 and 0.15, or the 2+3 lump geometry is simply unstable at this mass, is what
campaign F (§6) is running to decide. Nothing is in `research.tex` yet.

---

# PART A — FOR THE PAPER

## 1. MUST NOT say [§1, §10]

| claim | why it's dead |
|---|---|
| `E_pump ~ 1.6e-17`, "soft trajectory guide" | diagnostic evaluated the pump envelope in absolute, not centre-relative, coords — sampled ~64 units from where the pump acts, and underflowed. Corrected: `1.56e-17 → 1.16e-4` (13 orders). Real control effort is O(1) [§1.1] |
| "confinement saturates at `t_pump ≈ 8`" | product of the governor-contaminated campaign [§1.2, §3] |
| anything resting on the controller reservoir | ansatz mathematically unstable, abandoned [§3] |
| "the pump delays black hole formation" (horizons t=15.8 → 24.5) | `theta_plus`/`max_ah_r` are pointwise proxies, not surface-integrated expansion [§12] |
| a black hole / apparent horizon forms in **any** run | no AH finder exists here [§12, §20]. The shell-averaged expansion (§16b) is a better instrument, still unvalidated and still not a finder |
| "the pump moves mass from the phantom to the canonical sector" | zero transfer; 100% numerator growth from injection [§18.1] |
| "the pump acts on / reaches the phantom sector" (pre-`19.8` runs) | a params-parsing bug routed **all five** spotlights to canonical; phantom received zero pump force in every bicomplex campaign to date [§19.8] |
| `f_geo^evol` is unmeasured or zero **for the exotic pump runs** | it is 12.35–13.13%, trustworthy at 257³ [§17]. `research.tex` line 172 is the *canonical-only positive-energy* sector — a different sector. Do not conflate |
| retention scaling with pump duration below `t_pump = 16` | the coordinate-volume ladder's smooth rise across tp4–tp16 is mostly a χ-weighting artifact [§16.3-3] |

Also: `pump_work` (`collapse_diagnostics` col 15) changed meaning at `fe5ef9f8`.
Never compare it against runs built before that commit.

## 2. Corrections required before the rewrite

1. **E_pump withdrawal** (§1 above).
2. **Emission-sweep figure caption.** At `t_emit=12` the bundles are INCOMPLETE
   — RM 2/5, RC 0/5, RF 4/5 rays — and RM's two "arrivals" ran 2.75 units past
   the stack end. The caption's "every launch has 5/5 rays reaching" is **false
   at t=12**; the plotted zeros are bundle failures, not measured zero
   advantage. `t_emit ≤ 10` is clean 5/5; the peak (t_emit=4) is untouched.
   "The window brackets the channel's life" survives via the 20%→12% decay over
   `t_emit` 4→10. [§16.4]
3. **Trapped-surface claim.** Abstract/discussion say "collapses to an apparent
   horizon by t≈27" with "θ+≈−0.9, r_AH≈10.8". In RM's data θ+ first goes
   negative at **t=18.71**. The signal is *plausible* (r_at_min migrates inward
   9.6→5.7, no recovery, not attributable to −⅔K at K=0.27) but it is still a
   pointwise proxy on a partial footprint; "corroborated" overstates it and
   t≈27 matches no feature of the data. The headline ray (arrival 15.5)
   precedes even the earlier onset, so the transport certificate is safe.
   [§16.7]
4. **Confinement percentages** ("72%→11%" RM, "→2%" pump-free; spread 3.27).
   Produced by the old sector weight and coordinate volume. Cannot be
   recomputed — plotfiles gone. Quote the restarted campaign's proper-volume
   columns instead. Dispersal itself is not in doubt. [§16.7-4]
5. **Constraint norms must be re-measured, not merely caveated.** `L2_Ham` /
   `L2_Mom` are level-0-only, unmasked and domain-diluted — a control-loop
   number, not an accuracy figure (§4 item 0). **No absolute constraint value
   may be quoted until the composite AMR norm exists**; comparisons between
   runs are unaffected. Two caveats that do stay caveats: the geodesic probe
   grid is **fixed at 65³/dx=0.25** while sim resolution varies, and `f_op` is a
   **2D y-midplane** probe. [§16.5, §16.6]
6. **Methodological point worth including**: a null-geodesic FTL measurement is
   only as good as the metric it is integrated through, and a uniform resample
   of an AMR grid is lossy in a way that is *invisible downstream*. Any reported
   `f_geo` should carry evidence that the sampled metric reproduces the run's
   own `min_chi`. That check is `cache_fidelity()` and costs nothing. [§15]

## 3. Safe to quote [§10]

| result | value | caveat it must carry |
|---|---|---|
| constraints ignore the pump | `t_pump ≤ 24` at 256³ (≤16 fast tier); mean `L2_Ham` 2.71e-3 → 3.40e-3 over 24 units = 26% spread; peak **non-monotonic** (unpumped 4.27e-3 sits between tp4 3.90e-3 and tp16 4.68e-3) | quote as a **comparison only** — the norm is mis-measured (§4 item 0), so the ratios stand but the absolute values must not be printed [§16.5] |
| Duhamel bound | satisfied, 1.05–1.73× | §4's table compares the wrong pair — bound is on the *pump-attributed* violation; the meaningful check (pumped − pump-free ≲1e-3 vs bound 3.6e-3) passes with more room [§16.6] |
| governor never engaged | 1.0000 for every physically valid run | — |
| control effort | measurable, O(1), not 1e-17 | — |
| tp30 lapse-collapse instability | resolution-independent (fast tier + 256³) | — |
| headline evolving `f_geo` | **20.148%**, `t_arrival` 15.4987 — reproduces bit-exactly from the stored RM stack | RM stack passes `cache_fidelity` at every slice [§16.1-2] |
| `f_geo_evol`, exotic pump runs | 12.35–13.13% at 257³, all 5/5 rays | no dose–response; tp16/24/30 tie by construction [§17] |
| the pump does **not** improve `f_geo` | pump-free 27.0% beats tp16/24/30 (23.5%); tp8 best at 30.9% | *plausible but unproven* — same resample class as §15, milder [§7, §15] |
| retention gain (fractions) | 13.0× at tp24 coordinate-volume; **16.7× proper-volume** — 13.0× was a LOWER bound | `confined_frac` is a ratio with a growing denominator — prefer absolute [§19.1] |

**Coordinate vs proper volume, HQ ladder at t=30** [§16.3-3]. The physical mass
element is `w·√γ d³x = w·χ^(−3/2) d³x`, not `w d³x`. Proper weighting *lowers*
tp0 and *raises* every pumped run, so the gain goes UP — **13.0× was a lower
bound**:

| run | frac (coord) | frac (proper) | gain (coord) | gain (proper) | min_chi |
|---|---|---|---|---|---|
| tp0 | 0.0169 | 0.0159 | 1.00 | 1.00 | 5.68e-4 |
| tp4 | 0.0965 | 0.2205 | 5.72 | 13.89 | 5.96e-2 |
| tp8 | 0.1415 | 0.2242 | 8.38 | 14.13 | 1.72e-1 |
| tp16 | 0.2014 | 0.2245 | 11.93 | 14.15 | 4.03e-1 |
| tp24 | 0.2683 | 0.2656 | 15.89 | 16.74 | 5.83e-1 |
| tp30 | 0.2875 | 0.2956 | 17.03 | 18.63 | 2.91e-1 |

**And the ladder's shape changes.** Under proper weighting the tp4→tp16
dose–response **disappears** (0.2205 / 0.2242 / 0.2245 — 1.8% spread, against a
2.1× spread in coordinate volume). A real dose–response survives only at
tp24/tp30. The coordinate ladder's smooth monotonic rise across tp4–tp16 is
**mostly a χ-weighting artifact** driven by how deeply each run's core collapsed
— the same error class as the metric cache, biased against the same run. tp0's
χ_min=5.7e-4 core is nearly empty by t=30, so up-weighting it adds almost
nothing to the numerator while the same weight applies to the populated pumped
cores.

**Preferred confinement metric — absolute per-sector activity** [§19.1]:
```
abs_total = total_activity × confined_frac                                (cols 2 × 5)
abs_canon = total_activity × canon_mass_frac × confined_frac_canon        (2 × 17 × 14)
abs_phan  = total_activity × (1 − canon_mass_frac) × confined_frac_phantom(2 × 17 × 16)
```
`confined_frac` is meaningful only when `total_activity` is flat — true of no
run in this campaign (68.6 → 171.1 unpumped, 68.6 → 279.0 always-on).
Caveats that stay attached: activity is an **amplitude** measure (`A√(1+ω²)`),
not an energy; `canon_mass_frac` is an amplitude *share* and a misnomer — quote
it as "canonical share of total field activity", never as a mass fraction, and
never across configs with different ω [§18.1].

---

# PART B — OPEN WORK

## 4. Code fixes, priority order

0. **BUG — constraint norms were measured on the wrong grid** (found 2026-07-29;
   **code fixed `7b57f2a1`+, analysis still open**). `constraint_norms.dat`
   cols 2–3 are a *control-loop input* that had been quoted as a physics result.
   Three compounding defects, all at
   [RadialRecipeLevel.cpp:438](../../Examples/RadialRecipe/RadialRecipeLevel.cpp#L438):
   * gated to `Level() == 0`, so the three refined levels — finest `dx=0.0625`,
     8× level 0, and where all the matter is — never contribute;
   * **no covered-cell mask**: coarse cells underneath the refinement are still
     summed, so the active region enters at 8×-too-coarse resolution and its
     real residual never appears at all;
   * normalised by the **whole** domain volume, ~99.9% of which is near-vacuum
     (level 1 covers 0.116% of the domain), so the mean is diluted toward zero.

   For its actual job — a cheap, consistent, once-per-coarse-step scalar
   feeding the pump governor — it is fine. As a published accuracy figure it is
   not, and it was never meant to be one.

   **DONE — six new append-only columns 12–17** (`L2_Ham_amr`, `L2_Mom_amr`,
   `L2_Ham_amr_rel`, `L2_Mom_amr_rel`, `Linf_Ham_amr`, `L2_Ham_amr_ref`). The
   composite evaluates the constraints on **every** level, masks coarse cells
   covered by finer ones (`amrex::makeFineMask`), weights by each level's own
   `cell_vol`, drops 4 level-0 cells of outer boundary (BC violation is not
   physics), and additionally reports the **undiluted peak** `|Ham|` and the
   norm **restricted to the refined region** (levels ≥ 1) so the vacuum dilution
   is visible rather than hidden.

   **Cols 2–3 and the governor's input are deliberately untouched.** Changing
   them would be a *physics* change, not a readout change — `L2_Ham` feeds
   `publish_cached_L2_Ham` → the pump governor (tanh, centre 0.035), and the
   composite is a different magnitude, so recentring the governor would
   invalidate every prior run's comparability. Verified: cols 1–11 **bit-identical**
   to stored `pcc_t010` over a 200-step clone; `collapse_diagnostics`,
   `energy_conditions`, `curvature_invariants` identical to the pre-fix binary
   run under the same conditions. Cost is **free** — paired A/B on idle GPUs,
   1.1653 s vs 1.1649 s per coarse step (+0.03%). The naive version was +13%;
   all of it was recomputing level 0 twice, removed by passing the caller's
   level-0 constraint MultiFab in. Old binary kept as
   `main3d.gnu.CUDA.ex.pre_amrnorm`.

   **Still open — the analysis, which is the part that gates the paper:**
   * **Redo the Richardson test.** p≈3.3 certifies the convergence of the
     level-0 restriction, not of the solution. Needs 2–3 resolutions scored on
     `L2_Ham_amr`.
   * **Explain `L2_Ham_rel` ≈ 0.64–0.81.** All four campaign F arms
     (0.65/0.72/0.64/0.64 at t=10 → 0.74/0.81/0.74/0.69 at t=19): the residual
     is comparable to the terms it must cancel. First composite reading on a
     `pcc_t010` clone at t≤2 gives `L2_Ham_amr_rel` ≈ 0.22–0.40, i.e. **lower**
     but still large. Cannot go into the paper unexplained.
   * **Quote `Linf_Ham_amr` and `L2_Ham_amr_ref` beside the mean.** The t≤2
     clone already shows peak `|Ham|` ≈ **218×** the composite L2, and the
     refined-region norm ~30× the domain mean — the dilution is not a rounding
     detail. These are early-time, barely-refined figures; the honest numbers
     need a full run on the new binary.

   **Scope — rankings survive, absolutes do not.** Every run to date is measured
   the same wrong way, so all within-protocol comparisons (§9, and the whole
   "constraints ignore the pump" argument) stand. What is not usable is any
   absolute claim of the form "the constraint violation of this run is 4e-3".
   **Campaigns A–F predate the new columns and have 11, not 17.**

0b. **BUG — the `theta_plus` fix treated the wrong cause** (found 2026-07-29 by
   its own field A/B, §16a). `e01ec730` masked refinement-edge ghosts; the false
   trapped surface at r≈10.3 is unchanged, because the minimum sits at an
   *interior* cell at the footprint's outer radius and re-pins one cell in.
   Theta tracks global `max|K|`, not local geometry. **Nothing published moves**
   — collapse calls were always made on lapse/chi/max|K| — but the three theta
   columns are not evidence of anything and must not be quoted. Cheapest real
   fix: restrict the reduction to r < `recipe_basis_radius_max` so the answer is
   not simply "the edge of the box I looked in". Proper fix: surface-integrated
   expansion on a located surface, i.e. an actual horizon finder.

0c. **Metric stack under-resolves collapse-epoch slices.** The geodesic scorer's
   fidelity gate refused `pcf_t010_t60` (1/28 slices, true `min_chi` 4.889e-02
   read as 7.744e-02 at 257³). The gate is working; the cache is too coarse.
   **Every t≥40 arm will hit this**, so raise `GRTECLYN_METRIC_STACK_N_SPACE`
   for endurance runs — never `--force`, which is exactly how a spurious FTL
   shortcut would enter the record.

1. **Profile-matched target.** Both the overlap strip [§19.10] and the
   bridge-feeding [§19.11] are *shape* errors: the controller aims at an
   analytic sech whose skirts do not match the real lump. Measure the actual
   radial profile from a `pcc_t010` snapshot (retained for this) and drive
   toward it. Subsumes both defects; makes the pump a pure stabiliser rather
   than a source. **Blocker:** retained snapshots are metric-only — needs a
   radial-profile extraction added to the consumer first.
2. **Adaptive aim.** Scale target amplitude down as `min_lapse`/`min_chi` fall
   — closes the loop on the measured stability-vs-warp tension [§19.13].
3. **Deferred hygiene** (none blocking): `wiring.py`
   `EVOLUTION_MATTER_KEYS`/`_pump_controller_overrides` lack the phantom and
   superpose knobs, so generated configs *silently cannot set them*;
   `RadialRecipeLevel.cpp` `reservoir_mode_override` race;
   `structure_coherence` NaN in every row [§18.7]; failing
   `test_horizon_finder_guard.py` (22-col wormhole fixture vs ≥23 required).

## 5. Open scientific questions

* **Does any configuration hold to t ≥ 40?** **PARTLY ANSWERED (campaign F, §6).**
  `pcf_t012_t60` (aim 0.12) cleared t=40 with margin — lapse 0.378, chi 0.231 —
  the first arm ever to do so, then collapsed at t=43.86. So the bar is
  reachable, but *nothing has yet held indefinitely*: every configuration still
  collapses and only the timing moves. **Open:** where the improvement turns
  over. Kill time is monotonic in aim across four runs (0.06 → 26.7, 0.08 → 32,
  0.10 → 39.6, 0.12 → 43.9), so the next rungs are **0.14 and 0.16**. If it does
  not turn over by 0.16 the aim is not the limit and the lever is the
  *configuration* (lump count, separation, exotic fraction, rotation).
* **Stability vs warp strength — the assumed trade did not appear over 0.10→0.12.**
  The arm with the strongest warp signature (`pcb_match`, f_geo 0.279, local
  speed 1.78) ends deepest (lapse 0.28); the healthiest (`pcc_t010`, lapse 0.63)
  carries 65% of that signature; the flattest (`tw2`, 115 units diffuse) the
  least (0.089) — that is the concentration-buys-warp-costs-stability reading.
  But campaign F's 0.12 arm held **more** matter than 0.10 (68.7 vs 56.1 at
  t=20) **and** survived longer with ~2× the lapse/chi margin. Over that
  interval both improve together, so the tension is not universal; either it
  reappears above 0.12 or the earlier trade was confounded by the pre-fix
  over-feeding law. The 0.14/0.16 rungs decide it.
* **Sectors are not independently tunable.** `pcf_split_t60` (canonical 0.08,
  exotic 0.12) died at t=29.7 — earlier than uniform aiming at *either*
  constituent value, while holding the richest exotic sector of the four. Why a
  mixed-aim geometry is worse than both its endpoints is unexplained.
* **Where between `t_pump` 24 and 30 does driving turn destructive?** Needs
  rungs at 26/28; tp24 and tp30 are bit-identical to t=24 by construction.
* **`f_geo_evol` vs pump duration** is unresolved by design: a `t_emit=0` ray
  arrives at t≈12.6, before tp16/24/30 differ. Needs an emission sweep with
  `GEODESIC_EMIT_MIN_TIME` past each rung's stop time.
* **The §16.3-4 coherence-gate fix is unexercised** — committed, correct by
  inspection, never run on bicomplex data end-to-end [§18.7].

## 6. Campaign F — COMPLETE (`runs/pump_confine_f`, 2026-07-29 00:32–02:41)

All arms clone `pcc_t010` (previous best), all to **t=60**, binary `e01ec730`,
routing echo `parsed: 1 -1 -1 1 -1` verified in all four logs. Every arm was
ended by the watchdog; none reached t=60. **No arm crossed governor 0.5 at any
time** (`governor_min = 1.0000` in all four) — no result here was bought with
constraint violation.

| arm | aim | killed at | trigger | t_end | end canon | end exotic |
|---|---|---|---|---|---|---|
| `pcf_t012_t60` | 0.12 flat | **43.86** | `min_chi` 0.0495 | 43.2 | 32.4 | 35.7 |
| `pcf_t010_t60` | 0.10 flat | 39.62 | `min_lapse` 0.1432 | 38.9 | 27.0 | 29.1 |
| `pcf_split_t60` | 0.08 / 0.12 | 29.67 | `min_chi` 0.0231 | 28.8 | 22.3 | 37.0 |
| `pcf_sup_a06_t60` | 0.06 superposed | 26.66 | `min_chi` 0.0475 | 25.9 | 20.1 | 31.4 |

**`pcf_t012_t60` is the first arm in any campaign to clear the t ≥ 40 bar.**
At t=40: total 65.0, canonical 33.1, exotic 36.0, `min_lapse` 0.378,
`min_chi` 0.231 — not marginal, comfortably above both kill floors. Score with
`score_confine_abs.py --endurance` (`79e965a8`; the old late window [20.16,
27.36] was built for t=30 runs and is blind to this regime).

**Aiming higher won on both axes at once.** 0.12 held *more* matter than 0.10
(68.7 vs 56.1 total at t=20) *and* survived 4.2 time units longer, with roughly
twice the lapse and chi margin throughout. This contradicts the
stability-vs-warp tension assumed in §19.13: over 0.10 → 0.12 there is no
trade, both improve. It does not establish monotonicity — 0.15 is untested
since aim 0.08 died at t≈32 and the cliff was inferred, not measured.

**Split aiming failed, and that is the informative negative.** Aiming canonical
low (0.08) and exotic high (0.12) died at t=29.7 — *earlier than flat aiming at
either constituent value*. The sectors are not independently tunable; the
mixed-aim geometry is worse than either uniform one. Its exotic sector was the
richest of the four at kill time (37.0), so the failure is geometric, not a
matter shortage.

**The superposed-target fix (`64c89be4`) hit its predicted number exactly**:
48.3 total at t=20 against the ~67 the over-feeding law produced. The law is
right; 0.06 is simply too low an aim to survive — it died earliest of the four.

**Answer to the open question: no stable window was found, but the trend now
points somewhere.** Every configuration still collapses; only the timing moves.
The ordering (0.06 → 26.7, 0.08 → ~32, 0.10 → 39.6, 0.12 → 43.9) is monotonic
in aim across four independent runs, so the next test is **0.14 and 0.16**, to
find where the improvement turns over. If it does not turn over by 0.16 the
limit is not the aim and the lever is the configuration (lump count/separation
/exotic fraction), per the §19.13 decision tree.

**Geodesic scores — the two short arms scored, both endurance arms were
refused.** `sup_a06` 8.16e-02 and `split` 8.41e-02 (5/5 rays, fidelity PASS).
The fidelity gate **refused** both survivors rather than emit a spurious
shortcut:

| arm | bad slices | worst | true `min_chi` | cached | error |
|---|---|---|---|---|---|
| `pcf_t010_t60` | 1 / 28 | t=38.88 | 4.889e-02 | 7.744e-02 | 2× too shallow |
| `pcf_t012_t60` | 2 / 31 | t=43.20 | 7.227e-02 | 1.929e-01 | 2.7× too shallow |

This is the gate working, not a failure — a well read 2–3× too shallow means
rays do not feel it and arrive early, which is exactly how a fake FTL result
would enter the record. **The pattern is structural: both refusals are in the
collapse epoch, and every arm that survives to t≥40 will hit it**, because the
deepest slices are precisely the ones 257³ cannot resolve. The consequence is
that *the arms that survive longest are currently the ones we cannot score for
FTL* — so raise `GRTECLYN_METRIC_STACK_N_SPACE` for endurance runs (§4 item 0c).
Never `--force`. Campaign F therefore has **no FTL number for its two best
arms**; the 8.2–8.4e-02 figures above are from the two that died early and are
not evidence about the endurance regime.

**Watchdog** (new, in `pump_confine_queue.py`): aborts on `min_lapse < 0.15`,
`min_chi < 0.05`, NaN, or a pre-t=27 breach of `L2_Ham > 0.035` / governor <0.5.
Offline replay kills `pcd_match_t60` at t=30.62 and `pcd_match_p0` at t=26.36,
and leaves `pcc_t010`/`pcc_match_tw2` untouched (they sit 4–10× clear).
Column indices are width-guarded so a layout change fails the watchdog *open*.
**Field record: fired correctly on all four arms**, and `min_lapse` led on
`pcf_t010_t60` (0.1432 while `min_chi` was 0.0287, still above its own floor),
confirming it as the earlier trigger.

---

# PART C — THE RECORD

## 7. Bugs found and fixed [§2]

| # | bug | where | fixed |
|---|---|---|---|
| 1 | pump diagnostic used absolute, not centre-relative, coords | `RadialRecipeLevel.cpp` | `fe5ef9f8` |
| 2 | PD law duplicated in 4 places; diagnostic copy had drifted | new `RLPumpForce.hpp` | `fe5ef9f8` |
| 3 | reservoir momentum source had the wrong SIGN | `ControllerReservoirMatter.hpp` | `fe5ef9f8` |
| 4 | reservoir sources carried a spurious lapse factor | same | `fe5ef9f8` |
| 5 | reservoir had NO transport terms (per-cell ODE) | same | `fe5ef9f8` |
| 6 | `f_i` was never measured at all | `RadialRecipeLevel.cpp` | `fe5ef9f8` |
| 7 | **controller ledger fed the pump SAFETY GOVERNOR** | `RadialRecipeConstraintNorms.hpp` | `d6a0c350` |
| 8 | five disagreeing ray trust bars | 5 consumer sites | `6daea1a3` |
| 9 | metric cache sampled below AMR depth | `metric_stack_cache` | `fbdbd740` |
| 10 | confinement weight dropped `phi`/`Pi` for bicomplex runs | `extraction/confinement.py` | audit |
| 11 | `total_activity` changed units row-to-row | same | audit |
| 12 | `confined_frac` was a coordinate-volume fraction | same | audit |
| 13 | coherence gate saw 1 of 8 matter components | `general.py` | audit |
| 14 | `energy_conditions.dat` used the wrong potential (μ=0) | `RadialRecipeLevel.cpp` | audit |
| 15 | `boundary_flux.dat` "radial" derivative was ∂φ/∂x | `probes/boundary.py` | audit |
| 16 | rays completed through FROZEN geometry past stack end | `EvolvingMetricField` | audit |
| 17 | detector-plane overshoot (bias ≤0.4pp LOW) | geodesic probe | audit |
| 18 | scoring pass recomputed a number it already had | `score_evolving_geodesic.py` | `4f31f33a` |
| 19 | **`recipe_scalar_field_signs` parsed only its first token** | `SimulationParameters.hpp` | §19.8 |
| 20 | superposed-target PD law (overlap strip) | `RLPumpForce.hpp` | `64c89be4` |
| 21 | `theta_plus` proxy read the refinement edge | `RadialRecipeLevel.cpp` | `e01ec730` — **partial: stencil contamination gone, false horizon REMAINS**, §16a. Shell-averaged replacement now built but **not field-validated**, §16b |
| 22 | constraint norms level-0-only, unmasked, domain-diluted | `RadialRecipeLevel.cpp` | cols 12-17 added; **analysis still open**, §4 item 0 |

### The ones with teeth

**#3 sign.** This codebase defines `j_i = Σ_A s_A(−Π_A ∂_i φ_A)`, so the pump's
contribution to matter momentum is `−f_i` and the reservoir must absorb `+f_i`
— **opposite in sign to the energy case**. Symptom: `L2_Mom` got 2.16× WORSE
with the ledger on.

**#4 lapse.** `rhs[c_Pi] += S_A` is a coordinate-time rate, so
`∂_t ρ|_pump = f_perp` exactly, with no α — already absorbed by `S_A = α J_A`.

**#7 governor.** `compute_radial_recipe_constraint_norms` publishes
`cached_L2_Ham`, which drives the pump governor (tanh, centre 0.035, width
0.003). It called `fill_active_constraints`, which in
`controller_reservoir_mode ≥ 1` **includes the reservoir**. The controller's own
bookkeeping was wired into its own safety interlock: when the ledger diverged
the governor throttled the pump to zero on the strength of a number describing
the ledger, not the spacetime. Every mode-1 run with an active pump hit governor
min 0.000, closing at t=6.90–7.99 — so tp8/tp16/tp30 came out **bit-identical**,
effective pump duration ~7 regardless of configured stop time.

**#9 metric cache — the big one.** The evolving-geodesic probe reads
`small_data/metric_stack/`, `g_{μν}` resampled onto a **uniform** grid
(`n_space³` over `±half_width`; defaults 33 search / 65 hq, `half_width=8`), so
`dx = 0.5`. The runs use `n_cell=256` over `L=128` with `max_level=3`, finest
`dx = 0.0625`. **The cache sampled at the coarsest level and discarded all three
levels of refinement** — no error, no warning, and the resampled metric is
smooth, positive-definite, and integrates fine. It is simply a different
spacetime.

Why it corrupted a *comparison* rather than shifting it — cached vs true
`min_chi` at t=30:

| run | true `min_chi` | cached | error |
|---|---|---|---|
| **tp0** | 5.682e-4 | 5.622e-2 | **99.0×** |
| tp4 | 5.955e-2 | 6.806e-2 | 1.1× |
| tp8 | 1.720e-1 | 1.717e-1 | 1.0× |
| tp16 | 4.026e-1 | 4.106e-1 | 1.0× |
| tp24 | 5.826e-1 | 5.930e-1 | 1.0× |
| tp30 | 2.914e-1 | 3.690e-1 | 1.3× |

Only the pump-free run develops a deep sharp well, so only its cache was wrong
— it inflated exactly one run, and that one was the control. Its rays never paid
the Shapiro delay, so it posted the biggest apparent shortcut *because* it was
the worst-resolved. **Caught by a physics objection, not a test**: a deep well
means light arrives *later*, so the deepest-well run reporting the largest
shortcut is backwards. Fix: `cache_fidelity()` compares each cached slice's
representable `min_chi` against what the simulation itself wrote to
`collapse_diagnostics.dat`; `required_n_space()` derives resolution from
`n_cell`/`box_length`/`max_level`/`half_width` (**257** for this ladder); the
scorer **refuses** an unfaithful cache (`--force` to override).
Correction to a secondary flag: `max(g_tt) > 0` appears at t=30 in **all**
stacks including fidelity-passing ones — β²>α² late in these gauges is real, not
per-se cache corruption. The 99× `min_chi` mismatch is untouched.

**#14 EC potential.** The diagnostic built `BiComplexScalarField(mass, lambda)`
with **μ silently defaulting to 0** while the evolution passes
`recipe_scalar_mu = 85333`. At candidate amplitudes (|φ|~0.08) the sextic term
is comparable to the mass and quartic terms. Measured old-vs-new on
bit-identical evolutions:

| column | changed? | magnitude |
|---|---|---|
| `matter_min_NEC` | **no** | kinetic-dominated, potential cancels |
| `matter_min_WEC` | yes, all rows | ≤137% rel.; **151 rows flip sign** |
| `matter_min_SEC` | yes, all rows | ≤141% rel.; **50 rows flip sign** |
| `matter_min_DEC` | yes, all rows | ≤102% rel.; **18 rows flip sign** |
| `matter_integral_NEC_violation` | **no** | NEC-derived |

Every flip runs the same way — old margin negative, corrected positive. **The
μ=0 bug fabricated energy-condition violations that do not exist.** Any paper
statement resting on a WEC/SEC/DEC *violation* in the canonical sector must be
re-derived. NEC and the phantom sector's NEC violation are unaffected.

**#19 field-signs parsing — invalidates every phantom-sector claim before it.**
`recipe_scalar_field_signs = 1 -1 -1 1 -1` is tokenized by ParmParse into five
values, but the loader read the key into a single `std::string`, and a scalar
string query returns **only token 0**. The parsed array was `[1,0,0,0,0]`, and
the router treats 0 as canonical. **All five spotlights drove Phi+; Phi− received
zero pump force in every bicomplex campaign to date.** (`load_rl_lump_seed_axis`
had the same bug; both fixed via `countval`+`getarr`, and the parsed signs are
now echoed to the log at startup.) Caught because `pca_pg2` (phantom gains 24/14)
and `pca_pg4` (48/28) produced `confinement.dat` rows **bit-identical to each
other and to `lad_m0_tp30`** — two different gains cannot produce identical bits
unless the gain multiplies a force that is exactly zero.
Retro-explains at one stroke: the "injection asymmetry" [§19.7], "pure canonical
injection" [§18.1], the `pca_pwd_up` NaN (it tripled the *canonical* drive), and
the phantom retention gains (real measurements, wrong attribution — the benefit
was **indirect**, geometry sourced by the pumped canonical sector).

**#21 theta_plus refinement edge.** See §16/§20 below — and §16a, where the
field A/B shows the fix removed the secondary cause only. The proxy still
reports a trapped surface at r≈10.3 with the interior at chi≈0.6.

## 8. Abandoned / retracted — do not re-derive

**The controller energy–momentum reservoir** [§3]. Implemented (modes 0/1/2),
fully debugged, then abandoned — **the ansatz is unsound, not the code.**
`S_c^{ij} = 0` makes the reservoir pressureless dust: the momentum equation has
no `∇ρ` term, so nothing pushes back, leftover momentum is frozen in and keeps
feeding the energy equation forever. Density grows without bound *by
construction* — in exact arithmetic, in flat space, with a perfect
discretisation. Decisive evidence: pump force is EXACTLY 0 from t=6 (pump off at
t=4) yet the ledger grows `1.9e-2 → 9.8e-1` by t=30. Mode 0 vs mode 1 ratio goes
0.36 → 335 over t=1→30; mode 2 died with NaN in K at t=9.69.
**Methodological warning:** validated at t=2 it looked excellent (mode 2: Ham
×0.22, Mom ×0.31); it only turns bad after t≈2.5. **Never validate this system
on short runs.** It turned out not to be needed — a direct ablation shows the
constraints barely notice the pump [§4].
Also: the original plan's densitised energy equation carries a spurious `α K ρ_c`
and fails the FLRW-dust check. Correct form for zero spatial stress:
`∂_t(√γ ρ) + ∂_i[√γ(α j^i − β^i ρ)] = −√γ j^i ∂_i α`.

**Per-sector PD gain as the phantom lever** [§18.6]. Dead. `pcb_match_pg4`
(gains ×4) is identical to `pcb_match` to ~3 digits at every slice. With the
target matched the PD error is small, so quadrupling a multiplier on a near-zero
error changes nothing. What the phantom sector needed was **routing** [§19.8] and
a **correct target**, not more force.

**"Over-driving the target" as a confinement lever** [§19.5]. Fatal, not weak.
Raising `well_depth` 0.15 → 0.30 refined **12,255,232 level-3 cells vs 110,592**
(111×), ran 3× slower, and aborted with `NaN in K` at t=4.31. Re-attributed to
the canonical sector after §19.8.

**"The PD controller is mis-signed for phantom matter"** [§18.4]. False.
`GRTresnaIndependentScalars.impl.hpp:83` gives `∂_t φ = α Π + β·∂φ` for every
field regardless of sign; the phantom nature enters only the stress-energy
coupling. The PD target `tPi1 = ω·tphi2/α` is sign-correct in both sectors.

**"`env < 1e-8` makes it a local trap"** [§18.4]. Not with this config.
`rl_pump_target_profile = 2` is `sech(r/1.667)`, so the cutoff sits at r ≈ 32 —
far outside the matter region. The gate essentially never fires. What the
controller has is a *soft, exponentially decaying authority* (~5e-2 at r=6,
~1.5e-3 at r=12), not a hard boundary.

**"3-spotlight phase conflict"** [§19.4]. Falsified. `rl_pump_phase` is never
written in trajectory mode; all five sites share frequency 0.8 (= the initial
data's `bs_omega`, with `Pi_im = −ω·φ/α` painted, so pump and matter rotate in
phase by construction). `trajectory_lump<k>_phase0` is *orbital azimuth*, not
field phase. The real overlap effect is **summed per-site PD errors** [§19.10].

**`f_geo_evol` from the 33³ cache** [§13]. All values deleted, not corrected —
they came from a cache that could not represent the geometry. The machinery was
sound; the input was not. Replaced by §17's 257³ values. (Corollary: the void
values were numerically close to correct — the corruption lived at t≳26, which
`t_emit=0` rays never sample. They were still rightly retracted: that they came
out close could only be known by redoing them on a faithful cache.)

**`trajectory_lump<k>_exotic` is dead config** [§18.5]. Written into every
`params.txt` and never read — there is no `exotic` member on
`TrajectoryLumpParams`. Affects no result, but makes `params.txt` look as though
per-lump exotic fractions were in play when they were not.

---

# PART D — RESULTS

## 9. Constraints — the pump does not drive constraint growth [§4]

HQ 256³ (`runs/pump_ladder_m0_v1`), mode 0:

| t_pump | peak L2_Ham | mean | final | governor min |
|---|---|---|---|---|
| 0 (none) | 4.2723e-3 | 2.7116e-3 | 3.7868e-3 | 1.0000 |
| 4 | **3.9018e-3** | 2.8828e-3 | 2.9172e-3 | 1.0000 |
| 8 | 4.1969e-3 | 3.0880e-3 | 3.9745e-3 | 1.0000 |
| 16 | 4.6773e-3 | 3.1984e-3 | 4.6773e-3 | 1.0000 |
| 24 | 4.1474e-3 | 3.4027e-3 | 3.9829e-3 | 1.0000 |
| 30 | *4.7606e-2* | *5.0860e-3* | *4.7606e-2* | **0.0002** |

Duhamel: tp4 1.23×, tp8 1.05×, tp16 1.12×, tp24 1.73× — all satisfied, all
tight. A referee cannot claim the pump drives constraint growth when driving for
4–8 units gives a *lower* peak than not driving at all.
**Hard boundary:** holds to tp24 at HQ (16 fast tier), NOT to 30.
Analysis: `research/neuralspacetime/analysis/pump_constraint_budget.py`.

## 10. Geometry and collapse

**Lapse/chi at t=30, HQ ladder** — monotonic to tp24, then tp30 collapses [§5]:

| t_pump | min_lapse | min_chi | max_abs_K | max_Pi |
|---|---|---|---|---|
| 0 | 2.607e-2 | 5.682e-4 | 4.916e-1 | 6.61e-2 |
| 4 | 1.111e-1 | 5.956e-2 | 2.746e-1 | 5.13e-2 |
| 8 | 2.254e-1 | 1.720e-1 | 2.491e-1 | 8.08e-2 |
| 16 | 4.236e-1 | 4.026e-1 | 2.090e-1 | 1.79e-2 |
| 24 | **6.282e-1** | **5.826e-1** | 7.648e-2 | 3.36e-2 |
| 30 | *4.953e-5* | *2.914e-1* | 3.853e-1 | *5.24* |

**tp30 is a driven instability, not gravitational collapse** [§6]. `min_lapse`
falls 5 orders while `min_chi` barely moves (0.394 → 0.337) — in real collapse
the conformal factor collapses *with* the lapse. Lapse collapse with a healthy
metric = a gauge/driving pathology. Mechanism is measurable: pump force rises
**30×** from t=20 to t=28 (3.85e-6 → 1.11e-4) as the field drifts from target
and the PD pushes harder. Reproduces at 256³, so not a `dx=1.0` artifact. The
governor engaged legitimately in both.

**Central compaction is genuine** — pump-free `min_chi` = 5.7e-4 vs 0.583 at
tp24 is a real, large, monotonic difference. Calling it collapse requires a real
AH finder, which this codebase does not have [§12].

## 11. Confinement [§18, §19]

**Absolute per-sector activity, tp0 (off) vs tp30 (always on)** [§19.2]:

| t | tp0 total | canon | phan | tp30 total | canon | phan |
|---|---|---|---|---|---|---|
| 0.00 | 48.1 | 18.6 | 31.6 | 48.1 | 18.6 | 31.6 |
| 1.44 | 49.4 | 18.5 | 33.1 | **79.2** | **46.8** | 33.0 |
| 12.96 | 47.2 | 13.4 | 34.6 | 78.8 | 45.3 | 33.9 |
| 18.72 | 40.4 | 8.0 | 32.7 | 78.7 | 46.3 | 33.2 |
| 24.48 | 19.4 | 4.9 | 15.5 | 73.0 | 48.4 | 26.4 |
| 30.00 | 2.9 | 2.7 | 2.0 | 80.2 | 62.2* | 22.9 |

\* t>27.4 canonical rise sits inside the L2_Ham runaway — contaminated, not
recapture.

Three phases, not one decay: **injection** (t<1.44, canonical ×2.5, phantom
×1.0 exactly), **hold** (t≈1.4–19, both sectors flat in both runs — the fraction
was falling only because totals grew), **late divergence** (t>19; unpumped both
collapse, λ_phan ≈ 0.174/unit; pumped canonical holds and phantom declines at
0.025/unit — **7× slower**). End-state gains, absolute: total 27.8×, canonical
~23×, phantom ~12×.

**The exotic strip is geometric** [§19.10]. The legacy PD law takes one error
*per spotlight*, each against that site's own target lump alone, and sums the
forces. Where two same-sector lumps overlap, each site reads the neighbour's
contribution as *excess above its own target* and drives it down — and the
down-drives add. Geometry (centres from the elite's `initial_data.matter.json`,
signs `1 -1 -1 1 -1`, envelope `sech(r/1.667)`):

| pair | sector | separation | overlap at midpoint | sector outcome by t=1.44 |
|---|---|---|---|---|
| C0–C3 | canonical | 8.93 | 9.4e-3 (isolated) | **+11%** |
| P1–P2 | exotic | 4.64 | 0.12 | **−19%** |
| P2–P4 | exotic | 4.97 | 0.10 | (same chain) |

C0–P1 sit 1.98 apart but are cross-sector — different field components, no
crosstalk. Confirmed by `pcd_match_p0` (exotic pump literally OFF): exotic holds
33.05 at t=1.44, **no strip**. The legacy pump is what removes the matter.

**Two fixes for the strip, and both have a cost.**
* *Superposed-target law* [§19.11] (`rl_pump_superpose_targets`, `64c89be4`):
  sites of a sector accumulate ONE summed target and envelope, one PD error
  against the superposition, weight capped at `governor·min(Σenv, 1)`. Reduces
  to the legacy law for an isolated site; off by default (bit-identical).
  Strip eliminated (exotic 42.87 at t=1.44 vs 25.67 legacy) — **and
  over-corrected**: the summed sech skirts are fatter than the real three-lump
  field, so the controller feeds the bridges. 67 units of mass instead of 44 →
  collapse.
* *Aiming slightly ABOVE the matter's amplitude* (`pcc_t010`, aim 0.10) removes
  the strip under the **legacy** law with a one-line config change: with the
  target above the local field, a neighbour's contribution no longer reads as
  excess and the down-drive never fires. Exotic at t=1.44 is 31.75 — *above* its
  initial value.

**Direct exotic pumping is NET ESSENTIAL** [§19.13]. `pcd_match_p0` ends at
exotic 11.9 (−62%), `min_lapse` 0.015, `min_chi` 0.014 — the worst end-state of
any completed arm. Its `matter_min_NEC` goes −4.4e-3 (t=1) → **+2.8e-4 (t=18)**
and the integrated violation to **exactly 0.0**: the negative-energy content is
entirely gone, and the geometry collapses immediately after. "Legacy exotic
pumping is net harmful" is correctly scoped to **t ≲ 20 only**.

## 12. FTL / geodesic [§7, §17]

Frozen `f_geo` peak, HQ ladder (22/22 rows trustworthy): tp0 27.0%, tp4 29.5%,
tp8 **30.9%**, tp16/24/30 23.5%. `f_geo_evol` at 257³: tp0 12.35% (truncated
stack, `--max-time 27.5`), tp4 **13.13%**, tp8 12.61%, tp16/24/30 12.39%.
**No dose–response in either.** tp16/24/30 tie to the last digit *by
construction* — bit-identical until t=16, and a `t_emit=0` ray arrives at t≈12.6.

Confinement-campaign arms [§19.13]:

| arm | peak f_geo | end-to-end | max local speed | superluminal frac |
|---|---|---|---|---|
| pcb_match | **0.279** | 0.0798 | **1.78** | 0.988 |
| pcb_match_pg4 | 0.278 | 0.0789 | 1.78 | 0.988 |
| pcd_match_p0 | 0.266 | — | 1.51 | 0.963 |
| pcd_match_t60 | 0.257 | — | 1.74 | 0.949 |
| pcc_t010 | 0.182 | 0.0748 | 1.62 | 0.957 |
| pce_sup | 0.114 | — | 1.41 | 0.855 |
| pcc_match_tw2 | 0.089 | — | 1.34 | 0.939 |

## 13. Campaign finals [§19.13]

Start is always total 48.05 / canonical 18.62 / exotic 31.56.

| arm | aim | law | t_end | total | canon | exotic | min_chi | min_lapse(29) |
|---|---|---|---|---|---|---|---|---|
| pcb_match | 0.08 | legacy | 30 | 42.8 | 22.1 | 23.6 | 0.075 | 0.280 |
| pcb_match_pg4 | 0.08 | legacy | 30 | 43.5 | 22.0 | 24.4 | 0.085 | — |
| **pcc_t010** | **0.10** | legacy | 30 | 54.2 | 27.9 | **29.9** | **0.499** | **0.630** |
| pcc_match_tw2 | 0.08, w×2 | legacy | 30 | 115.0 | 61.1 | 58.0 | 0.427 | 0.488 |
| pcd_match_p0 | 0.08, exotic OFF | legacy | 30 | 25.2 | 13.3 | 11.9 | 0.014 | 0.015 |
| pce_sup | 0.08 | superposed | 20† | 66.7 | 26.9 | 43.4 | 4.1e-3 | 0.109 |
| pcb_base / pcb_pg4 | 0.15 | legacy | 12† | 87 | 43 | 48 | 1.6e-3 | — |
| pcd_match_t60 | 0.08 | legacy | 32† | — | — | — | 0.0153 | 0.108 |

† killed on the pre-registered `min_chi < 0.05` criterion / NaN.

**`pcc_t010` is the best configuration found so far** on every axis except warp
strength: healthiest geometry by 6× in `min_chi`, exotic sector essentially
conserved (−5%), no early strip.
**`pcd_match_t60` settles endurance: aim 0.08 dies at t ≈ 32.** Constraints
stayed clean throughout (L2_Ham ~5e-3, governor 1.000), so this is the geometry
running away, not the solver failing — and it means `pcb_match`'s "clean finish"
at t=30 **was collapse already in progress**.
**Target amplitude is the lever, not gain.** 0.15 is lethal (injects 48→87 by
t=11.5, `NaN in K` at t≈12.3); 0.10 conserves; 0.08 holds then dies at 32;
0.06+superposed is untested until campaign F.

## 14. Data provenance [§8, §9]

**VALID**
* `runs/pump_ladder_m0_v1/lad_m0_tp{0,4,8,16,24}` — HQ 256³, t=30, validated
  [§11]. Manuscript numbers for constraints/retention/lapse. `metric_stack`
  deleted, `f_geo_evol` reset [§15].
* `runs/pump_ladder_m0/` — 257³ rerun; the campaign the fixed pipeline produced.
* `runs/always_on_pump/hq146_m0_tp4_t30` — bit-identity reference for
  `lad_m0_tp4`. `hq146_m1_tp0_t30` — pump never on, reservoir ≡ 0.
* `runs/pump_ladder_fast/fast_tp{0,4,8,16}`; `fast_tp30` valid to t≈29.
* `runs/pump_confine_{b,c,d,e}` — see §13. **Ran the pre-fix binary**, so their
  `theta_plus` columns carry the §20 artifact — as do the post-fix ones, §16a.

**CONTAMINATED — do not use constraint columns (2,3,7,8)**
* `runs/always_on_pump/hq146_m1_tp{4,8,16,30}_t30` — ledger diverged AND
  governor closed at t≈7–8; effective pump duration ~7 in all.
* `runs/always_on_pump/hq146_m2_tp30_t30` — CRASHED, NaN in K at t=9.69.
Their *physics* outputs remain valid up to the point the governor closed.

**VOID / DELETED**
* `runs/pump_confine_a` — void as a sector experiment (§19.8 routing bug); dirs
  deleted, findings preserved here.
* `runs/reservoir_fix_check/` — deleted (abandoned ansatz, t=2 smoke tests).
* Pruned 2026-07-28 (141 G → 113 G): `pcb_base`, `pcb_pg4`, `pce_sup_t60`
  deleted entirely. Kept in full with snapshots: `pcc_t010`, `pcb_match`,
  `pcd_match_t60`. Measurement files kept but `metric_stack/*.npz` dropped:
  `pcb_match_pg4`, `pcc_match_tw2`, `pcd_match_p0`, `pce_sup` — the 4D geodesic
  re-scorer cannot be re-run on these without repeating the simulation.
* Per-sector *spatial* analysis is no longer possible from stored data
  (`/tmp/grteclyn_scratch` plotfiles purged); `confinement.dat`'s 22-row moments
  are all that remain [§18.8].

### Column layouts

`constraint_norms.dat` (17, APPEND-ONLY; runs before 2026-07-29 have 11):
```
1 time  2 L2_Ham  3 L2_Mom  4 min_rho_req  5 max_rho_req  6 integral_neg_rho
7 L2_Ham_rel  8 L2_Mom_rel  9 pump_force_L2  10 governor  11 pump_fi_L2
12 L2_Ham_amr  13 L2_Mom_amr  14 L2_Ham_amr_rel  15 L2_Mom_amr_rel
16 Linf_Ham_amr  17 L2_Ham_amr_ref
```
**Cols 2-3 are the GOVERNOR'S INPUT, not an accuracy figure. Quote 12-17.**
Cols 12-17 are the whole-hierarchy composite (all levels, covered cells masked,
4 level-0 cells of boundary dropped); 16 is the undiluted peak `|Ham|` and 17 is
restricted to the refined region. See §4 item 0.
`pump_fi_L2` is the first real measurement of the momentum force density, and
**`f_i` is ~3.5× LARGER than `f_perp`** — the controller's momentum forcing
dominates its energy forcing, and that is exactly the component the old code
never measured and got the sign wrong on.

`collapse_diagnostics.dat` (15):
```
1 time  2 min_lapse  3 min_chi  4 max_abs_K  5-7 min_lapse_{x,y,z}
8 max_ah_r  9 min_theta_plus  10 r_at_min_theta_plus
11 min_phi  12 max_phi  13 min_Pi  14 max_Pi  15 pump_work
```

### What the C++ diagnostics actually measure [§16.5]

* **`constraint_norms.dat` is LEVEL 0 ONLY** (`Level()==0`, level-0 state,
  single cell_vol), **unmasked** (coarse cells under the refinement are still
  summed, at 8×-too-coarse resolution) and **domain-diluted** (normalised by the
  full volume, ~99.9% near-vacuum). So are the pump-force norms and the
  governor's input. Within-protocol comparisons (the entire §4 argument) are
  unaffected — every run is measured the same way — but absolute numbers are the
  level-0 restriction, and Richardson p≈3.3 is the convergence of *that
  restriction*. **This is a bug to fix, not a caveat to state — see §4 item 0.**
* **`collapse_diagnostics` / `energy_conditions` / `curvature_invariants` reduce
  over the FINEST level's boxes only** — 0.01–0.3% of the domain mid-run.
  `min_chi`/`min_lapse` do match global values (refinement tracks the well), but
  `max_Pi` is already wrong: on tp24 the footprint's canonical-Pi range was
  [−1.4e-2, −3.1e-3] while the global max was +5.4e-4 and the true activity peak
  sat in `Pi_lump0` (+8.4e-2) — **a component the diagnostic never reads**
  (`c_phi`/`c_Pi` only). Fine as a blow-up flag; not quotable as "the field
  momentum".

### Operational notes

* **Plotfiles must go to node-local NVMe, never NFS.** Measured: NFS 2-concurrent
  gave 0.333 t/min with extraction max 310.3 s (over the 288 s plotfile cadence
  — that is how backlogs formed); node-local 6-concurrent gives 0.335–0.378
  t/min with extraction max 21.9 s. **3.2× on campaign wall-clock at unchanged
  per-run speed**, which is the proof the I/O wall is gone.
* **Launch hygiene:** `export VARS && nohup python … & for k in …; do …; done`
  backgrounds the ENTIRE `&&` list, so exports live in that subshell and every
  loop launch runs WITHOUT them — silently falling back to search mode while
  still writing plausible result files. Export on its own line, then launch.
  Detection: first log line `mode=hq`, result line `rays=5/5`.
* **Scoring is ~3 min/run at ~25 GB** post-`4f31f33a` (was >60 min unfinished at
  90–110 GB). The dominant cost had been recomputing the frozen peak from the
  float32 cache — a number the consumer already measured at full AMR fidelity
  into `ftl_timeseries.dat` col 3.
* **`--max-time <t>`** truncates a stack before late under-resolved slices.
  Conservative by construction: the strict no-frozen-tail guard fails any ray
  still in flight past the last kept slice, so truncation can only lose rays,
  never fabricate arrivals.
* Future campaigns wanting HQ 4D traces must set
  `GRTECLYN_EVOLVING_GEODESIC_MODE=hq` and `GRTECLYN_METRIC_STACK_N_SPACE` (257
  here) **before launch** — the cache resolution is fixed at write time and
  cannot be raised afterwards.

## 15. Metric audit verdicts [§16]

Every measurement pipeline was audited after the cache bug. Verdicts:

| metric | verdict |
|---|---|
| evolving `f_geo` (headline 20.2%) | **SOUND** — RM stack passes `cache_fidelity`; trace reproduces bit-exactly. Two probe defects found (frozen tail, overshoot) that do NOT touch the t_emit=4 headline |
| frozen `f_geo` (29.5%) | SOUND for candidate-146; resample caveat stands for collapsing runs |
| `f_ff` (7.6%) | math audited clean; reads fidelity-passing stacks |
| `f_op` | machinery exact per edge; **2D y-midplane slice only** |
| `L2_Ham`/`L2_Mom`/governor | **BUG — level-0 only, unmasked, domain-diluted** (§4 item 0); within-protocol comparisons stand, absolutes not quotable |
| `collapse_diagnostics` | **finest-level footprint only**; `max_Pi` reads 1 of 8 matter components |
| `energy_conditions.dat` | **BUG — μ=0**, fixed |
| `curvature_invariants.dat` | formulas correct; finest-footprint caveat |
| confinement / dispersion | **3 BUGS, fixed** |
| `structure_coherence` gate | **BUG — canonical-Re-only weight**, fixed (unexercised) |
| `boundary_flux.dat` | **BUG — "radial" derivative was ∂φ/∂x**, fixed |
| Duhamel bound | integrator fine; §4's ratio table mixes categories |
| trapped-surface claim | pointwise proxy; paper's t≈27 matches no feature of the data |
| Alcubierre control (32%) | 5-slice ±0.4 stack for a ~7-unit flight — validates shortcut *detection*, mostly frozen bubble |

**The paper's stacks PASS the fidelity check.** All five surviving candidate-146
stacks (RM, RC, three freefall twins; 65³, dx=0.25, `max_level=3`): every slice
passes, at every time, in every run. Candidate-146 never develops a feature the
probe grid cannot hold — RM `min_chi` at t=30 is 5.96e-2 and the headline flight
ends at t=15.5 where `min_chi` ≈ 0.8. **The cache bug does not invalidate the
paper's headline.**

**Checked and CLEAN**, for the record: confinement AMR integral hygiene
(Σ cell_volume = V_domain exactly; yt child-masking correct); code-unit handling
and R_conf plumbing; no radiation-halo contamination (99.57% of activity inside
r<15); Dijkstra edge speeds solve the exact per-edge null condition; both
probes' flat baselines agree; trilinear interpolation and null re-projection;
`future_null_cov` future-directed root selection; the freefall probe's mass-shell
projection, Gram–Schmidt tetrad and proper-length D0 integral; `rays_complete`
applied at all five consumer sites; Duhamel trapezoid and 16π/8π factors;
`remove_duplicate_time_data` on restart paths.

**Caveats that stay caveats** (documented, not code-fixed): `f_op`'s 2D slice;
the Duhamel table's category mix; the **fixed 65³/dx=0.25 probe grid** while sim
resolution varies (so part of the f_geo ladder agreement reflects the common
probe grid, not pure sim convergence); temporal cadence Δt=0.72 (RM) / 0.24 (RF)
vs the ladder rerun's 2.88 — the 257³ spatial fix did not change cadence; the
canonical-only bound reads "none found within this ansatz/budget/resolution",
which is how the paper already phrases it.

## 16. The theta_plus proxy read the refinement edge (`e01ec730`) [§20]

**Diagnostics-only: no simulation number moves. Removes a false-horizon signal;
does NOT make the proxy a horizon finder.**

**Bug.** `theta_plus` is reduced over the finest level's valid cells, but its
radial derivative is a centred stencil. At a cell on that level's outer
boundary, `chi[i±1]` is a ghost interpolated down from the coarse level — so the
stencil differences across a resolution jump, and since `theta_plus` carries
`−dchi_dr/√chi` with no smoothing, the jump lands in the reported minimum. The
reduction takes the min over all cells, so it selects the most contaminated cell
in the domain.

**Evidence it was an artifact** (`pcc_t010`, pre-fix binary):

| t | min_theta_plus | r_at_min | max_ah_r | min_lapse | min_chi |
|---|---|---|---|---|---|
| 15.01 | +0.1190 | 10.26 | 0.00 | 0.852 | 0.840 |
| 21.01 | +0.0657 | 10.24 | 0.00 | 0.770 | 0.741 |
| 27.01 | **−0.0282** | 10.35 | 11.16 | 0.676 | 0.601 |
| 30.00 | **−0.1161** | 10.40 | 11.92 | 0.601 | 0.499 |

The radius does not move for fifteen time units (a forming horizon migrates);
the interior stays healthy (a trapped surface at r≈11.9 enclosing everything
forces prompt collapse); and the radius is where refinement ends, just outside
`recipe_basis_radius_max = 8.0`, not where the matter is.

**Fix.** An `iMultiFab` mask on the finest level — 1 on valid cells, 0 on ghosts,
then `FillBoundary()` so ghosts covered by a *sibling box* recover 1 (an interior
box face is not a resolution jump). The proxy is evaluated only where all six
stencil neighbours are valid, in **both** reductions (the min, and the second
pass locating `r_at_min_theta_plus`). Exclusion cost is negligible: one cell at
finest dx = 0.0625, and refinement tracks the well, so a real horizon forms
interior to the finest level.

**Verified physics-neutral.** 200-step clone of `pcc_t010`, new binary vs stored
run: `constraint_norms.dat` and `energy_conditions.dat` **bit-identical**;
`collapse_diagnostics.dat` bit-identical on every column except
`min_theta_plus` (0.1636 → 0.1683 at t=2 — the contaminated cell read *lower*
than the true minimum, i.e. biased toward a false horizon) and `r_at_min`
(11.79 → 11.38). One `curvature_invariants` row differs by 1 ULP; the **pre-fix
binary reproduces the same difference**, so it is pre-existing GPU
reduction-order nondeterminism, not a regression. Old binary preserved as
`Examples/RadialRecipe/main3d.gnu.CUDA.ex.pre_thetafix`.

**CAVEAT.** `theta_plus`'s leading term is `2√chi/r`, so in a near-flat region it
is just `2/r`: the minimum over any footprint sits at the **largest radius in
that footprint by construction**, no pathology required. Campaign F reports
`r_at_min ≈ 109.55` early — the domain corner (√3·64 ≈ 110.85), one cell in, as
the mask now enforces. The fix removes stencil contamination; it does not remove
the fact that this is a pointwise radial proxy on a reduction footprint rather
than a surface-integrated expansion on a located surface. Note the score-side
guard `HORIZON_OFFCENTER_FRACTION = 0.5` (`metrics/score/horizon.py`) is
calibrated for corner-origin miscentering and would not have caught this; it is
unchanged and still does not cover it.

**Scope.** Campaigns B–E ran the pre-fix binary, so their theta columns carry the
artifact — their collapse conclusions are unaffected because they were already
derived from lapse/chi/max|K| by hand. The paper is untouched: §16.7's
correction concerns RM data where the proxy migrates *inward*, the opposite
signature. `pcf_t010_t60` is the field A/B.

### 16a. FIELD VALIDATION: THE FIX DID NOT ACHIEVE ITS PURPOSE

`pcf_t010_t60` (fixed binary) vs `pcc_t010` (pre-fix), same config, same seed,
3000 rows joined on time. **The false-horizon signal is unchanged.**

| window | `min_theta_plus` | `r_at_min` |
|---|---|---|
| t < 14 | differs in ~50% of rows | differs (e.g. 11.792 → 11.376) |
| **t = 16 → 30** | **bit-identical** | **bit-identical, pinned 10.22 → 10.40** |

The fix does what it was written to do — it removes stencil contamination, and
early rows move. But at the epoch it was built for it changes *nothing*:

| t | theta | r_at_min | min_lapse | min_chi |
|---|---|---|---|---|
| 25.64 | −0.0001 | 10.29 | 0.702 | 0.639 |
| 27.91 | −0.0500 | 10.35 | — | 0.573 |
| 36.00 | −0.5861 | 10.50 | — | — |

Still claiming a trapped surface at r≈10.3 with the interior at chi ≈ 0.6.

**Why the diagnosis was wrong.** The pinned minimum is not at a ghost cell — it
is at an *interior* cell, at the largest radius the finest level's footprint
reaches. Masking one layer moves it one cell and it re-pins. The §16 CAVEAT
already said this ("the minimum over any footprint sits at the largest radius in
that footprint by construction") but it was filed as a caveat *to* the fix
rather than recognised as the dominant mechanism, so the fix was written for the
secondary cause. **The radius is essentially static (10.26 → 10.50 over twenty
time units) while theta plunges +0.11 → −0.59** — theta tracks global
`max|K|` growth (0.049 → 0.374), not a migrating surface. A pointwise proxy with
a `−(2/3)K` term evaluated on a bounded footprint goes negative at the footprint
edge whenever global K grows, regardless of local geometry.

**Consequence.** The proxy is **structurally unfit for collapse calls** and
masking cannot fix it. Do not quote `min_theta_plus`, `r_at_min_theta_plus` or
`max_ah_r` as horizon evidence in any campaign, pre- or post-fix. Collapse calls
stay on `min_lapse` / `min_chi` / `max|K|`, as they already were — so no
published conclusion moves. Options, in order of cost: drop the three columns;
restrict the reduction to r < `recipe_basis_radius_max` so the footprint edge is
not the answer; or replace it with a surface-integrated expansion on a located
surface (the only real fix). The `e01ec730` mask is harmless and stays.

### 16b. The third option is now implemented — and NOT yet field-validated

Status 2026-07-29. `RadialRecipeLevel::specificPostTimeStep` no longer takes a
pointwise minimum. It averages `theta_plus` over concentric coordinate spheres
centred on the **lapse minimum** (the one place a horizon can form first), and
counts a sphere only when the finest level covers ≥ 90% of its volume with
≥ 32 cells. That is §16a's option 3, and it attacks the actual mechanism: the
answer can no longer be "the largest radius in the footprint", because a
partially sampled sphere is discarded rather than averaged.

The three columns keep their positions (8, 9, 10 — see the column map above), so
nothing reading `collapse_diagnostics.dat` by index moves. `min_theta_plus` is
**NaN** when no sphere qualified: unmeasured, and it must never read as zero.

*This code had been written but never connected.* The output line still named the
three deleted pointwise variables, so the file could not compile and the binary
on disk predated it. Reconnecting it was a rename, not a physics change.

**What has been checked:** it builds, and on an Ellis–Bronnikov smoke run
(t=1, ml=1) it reports `max_ah_r = 0` with `min_theta_plus = +0.149` — no
trapped surface on data that has none, which is the trivially necessary
condition and nothing more.

**What has NOT been checked, and it is the one that matters:** the §16a field
A/B. The old proxy's signature was a radius pinned at 10.2–10.5 for twenty time
units while theta plunged with global `max|K|`. Whether the shell average breaks
that pinning can only be seen on a run that actually collapses — i.e. re-run the
`pcf_t010_t60` arm against this binary and join on time. **Until that A/B is
done, §16a stands unchanged: the three columns remain unquotable as horizon
evidence.** A new implementation is not a validated one — that is §17's seventh
habit applied to this very section, which earned its own entry for exactly this
mistake.

One diagnostic is computed and then discarded: `n_shells_ok`, the number of
spheres that passed the coverage test. It is the natural companion to a NaN
(*why* was nothing measured), and appending it would not move any existing
column. Left out deliberately — adding a column is a decision for whoever runs
the validation.

## 17. The recurring failure mode

Five bugs had the same shape: **a value was derived across a fidelity boundary
and then read as physics, while the authoritative source sat unused next to it.**

| bug | derived from | authoritative source that was ignored |
|---|---|---|
| trust bars [§14] | a re-derived `n_reached == n_rays` | the probe's own `rays_complete` |
| metric cache [§15] | a 33³ resample of an AMR grid | the run's own `min_chi` in `collapse_diagnostics.dat` |
| scoring cost [§17] | a recomputed frozen peak from the float32 cache | `ftl_timeseries.dat` col 3, measured at full AMR fidelity |
| confined_frac [§19.1] | a ratio with a growing denominator | `total_activity`, sitting in column 2 |
| theta_plus [§20] | a pointwise proxy at the edge of its own reduction footprint | the interior lapse/chi in the same file |

In every case the fix was to compare against the source rather than trust the
copy, and in every case the check was cheap enough that there was never a reason
not to.

A sixth, found 2026-07-29, is a *different* shape and worth naming separately:
**a number built as a control-loop input was promoted to a published result**
without anyone re-deriving what it measures (§4 item 0). It is correct for its
job and wrong for the paper, which is why no test would ever have flagged it.
Habit to adopt: before quoting any diagnostic, check what it was written *for*.

A seventh, found 2026-07-29 (§16a), is different again: **a fix was verified for
the wrong property.** `e01ec730` was checked exhaustively for physics-neutrality
— bit-identity across four diagnostic files, ULP-level accounting, a preserved
pre-fix binary — and it passed all of it. Nobody checked whether it *worked*.
The field A/B shows the false horizon is bit-identical over t=16–30, the exact
window the fix existed to repair. Worse, the correct explanation was already
written down in the same section, filed as a CAVEAT *to* the fix rather than
recognised as the dominant cause. Habits: **"physics unchanged" is not
"defect removed"** — every fix needs a test that fails before it and passes
after; and **when you write a caveat that explains the symptom, stop** — you
have probably just found the real cause and mislabelled it.

Two further habits earned the hard way: **never validate this system on
short runs** (the reservoir looked excellent at t=2 and diverges after t≈2.5;
`pcb_match` looked healthy at t=30 and collapsed at 32), and **when a result is
backwards, suspect the measurement** — the cache bug was caught by a physics
objection (a deep well means light arrives *later*), not by a test.
