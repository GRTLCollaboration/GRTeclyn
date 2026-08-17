# Bondi dipole — review fixes & rerun plan

Updated 2026-08-17, after two external reviews of `bondi_dipole.tex`.

---

## 0. State of the results — snapshot 2026-08-17 12:10

Live numbers, not final. Everything below is read from
`runs/bondi_rerun/<cell>/bondi_sg_*/small_data/sector_barycenters.dat`; the
three resolutions record on different time grids (frames every 40 *steps*, and
dt ∝ dx), so **every cross-resolution number here is linearly interpolated onto
a common time** — comparing raw rows silently compares different instants and
was the first thing that made the fit look broken.

### Item A convergence matrix — running, 4/6 past the science window

All six cells: `L_full = 64`, unigrid (`max_level 0`), `dt_multiplier = 0.02`,
`stop_time 60`, separation 8, lumps at rest. Only `N_full` differs.

| cell | N | dx | t reached | frames |
|---|---|---|---|---|
| `convA_pm_n128` | 128 | 0.500 | 60.01 done | 152 |
| `convA_pp_n128` | 128 | 0.500 | 60.01 done | 152 |
| `convA_pm_n192` | 192 | 0.333 | 39.5 | 149 |
| `convA_pp_n192` | 192 | 0.333 | 49.1 | 185 |
| `convA_pm_n256` | 256 | 0.250 | 12.0 | 61 |
| `convA_pp_n256` | 256 | 0.250 | 19.8 | 100 |

Measured GPU memory: 8.8 GB (N=128), 25 GB (N=192), **58 GB of 81 GB**
(N=256) — unigrid N=256 fits one H100 with ~23 GB spare, confirming the
`BONDI_MAXLEVEL=0` warning in LAUNCH.md. Evolution speed 180 / ~30 / 13.5 code
units per hour.

**The runaway is resolution-independent.** Mixed-pair centroid drift, all three
grids, common times:

| t | N=128 | N=192 | N=256 | spread |
|---|---|---|---|---|
| 5 | 0.00252 | 0.00260 | 0.00258 | 2% |
| 10 | 0.03422 | 0.03415 | 0.03377 | 1.3% |
| 12 | 0.05277 | 0.05225 | 0.05166 | 2.1% |

**The control's drift is a numerical artefact and it converges away.** `PP`
centroid displacement at t=10: **−0.01615 (128), −0.01042 (192), −0.00550
(256)** — falling with h at an apparent order between 1.1 and 2.2, i.e. → 0.
So the effect-to-artefact ratio *improves* with refinement: 2.1× at N=128,
3.3× at N=192, **6.1× at N=256**. This is the single most useful thing the
series has produced so far: the mixed-pair drift sits still under refinement
while the same-sector control's drift dies, which is exactly the signature a
referee asks for and cannot be faked by a converging discretisation.

**Two-grid error over the whole science window** (128 vs 192, mixed pair) —
smooth, single-signed, sub-percent:

| t | 10 | 15 | 20 | 25 | 30 | 35 |
|---|---|---|---|---|---|---|
| drift 128 | 0.0342 | 0.0897 | 0.2152 | 0.4337 | 0.7739 | 1.2721 |
| drift 192 | 0.0342 | 0.0887 | 0.2133 | 0.4309 | 0.7705 | 1.2639 |
| rel. diff | +0.2% | +1.1% | +0.9% | +0.7% | +0.5% | +0.7% |

So the published drift numbers are not going to move materially. Coordinate
separation over the same window: 7.96 → 7.33 (N=128).

**Do not fit the convergence order before t ≈ 15.** Attempted at t=10–12 it
returns p ≈ 0.1 and a 64% "extrapolation error" — a failed fit, not a result.
Two reasons, both structural: (i) the drift there is 10⁻²–10⁻³, so the
grid-to-grid differences (10⁻⁴) are the same size as other error sources, and
(ii) **each resolution solves its own initial data independently** to a 0.1%
non-linear tolerance, so the three cells do not start from the same state
(final Ham residual 8.3227e-2 / 8.3188e-2 / 8.3185e-2 — converging; Mom
residual 7.420e-2 / 7.517e-2 / 7.499e-2 — *not* monotone). Below t≈12 the
inter-grid differences also flip sign, which alone kills a three-point fit.
From t≈15 the ordering is stable (see the table above) and the fit is
well-posed. Expect the fitted order to come out **below** the formal scheme
order regardless: the drift is read from a threshold-weighted centroid, and
blunt diagnostics converge more slowly than the field they are built from.
Three grids give a rate *and* an extrapolation with no spare point to validate
the fit — if the rate lands far from expectation, add a fourth resolution
rather than forcing a curve through three points.

### Box and wave extraction, as actually configured (all six cells)

| quantity | value |
|---|---|
| evolution box | width 64, centre 32 32 32, outer wall 32 from centre |
| outer BC | Sommerfeld outgoing, non-periodic (`hi/lo_boundary = 1 1 1`) |
| **constraint-solve box** | **width 128** — already 2× the evolution box |
| Ψ₄ shells | **R = 8 and R = 16** (confirmed in the `.dat` headers) |
| built-in extraction | `activate_extraction = 0` — the in-code module is OFF |

The Ψ₄ numbers come from the plotfile consumer (`--radii`, wired from
`BONDI_RADII`), not from the evolution executable. **Neither shell is in the
wave zone**: the inner shell sits *at* the star separation (8) and the outer at
twice it, with the reflecting wall at only 32. Amplitudes at these radii carry
near-field terms, the two shells disagree by construction (r·ψ₄ falls 5–10×
between them), and no radiated-energy statement can be built on them. The
drift result is unaffected — it is a matter-distribution measurement and never
touches the shells. This is item C, now staged (below).

### Launch-path fixes made while getting the matrix running

* **All six cells initially aborted at parse time** with
  `ParmParse::addDefn(): no values for definition regrid_interval`. Cause:
  `regrid_intervals_for_max_level()` returned `[]` for `max_level ≤ 0`,
  emitting a key with no tokens; ParmParse rejects that when it *parses* the
  file, before any code reads the key. GRTeclyn `getarr`-reads exactly
  `max_level` entries and appends its own terminal 0
  (`Source/GRTeclynCore/AMReXParameters.hpp`), so unigrid must still emit one
  benign `0`. Fixed in the **shared** helper
  `grteclyn-wrapper/src/grteclyn_wrapper/core/params.py`, which every launcher
  (`search/optimize/driver.py`, `search/qd_search/driver.py`,
  `campaigns/hq/replay_eval.py`) goes through — not in any one campaign script.
* **Parallel plotfile extraction validated and in use.** `CONSUMER_JOBS=8`:
  A/B on 7 real plotfiles gave byte-identical output to serial (9/9 streams,
  all time-ordered) at 3.2 → 12.4 files/min. Serial consumption loses to
  production and idles the GPU 25–35 min per cell during the end-of-run drain.
  A 256³ plotfile is 2.9 GB and takes ~8× the serial processing of a 128³ one,
  so `-j 8` is mandatory, not an optimisation, at N=256.
* **Chained cells do not need to wait for the predecessor's post-processing.**
  `chain_next.sh` waits for the whole launcher to exit, but the GPU is free as
  soon as evolution ends — the remaining backlog is CPU-only. Two N=192 cells
  were started ~20 min early by killing the watcher and launching directly.
  If you do this, **kill the watcher first**, or it fires a duplicate later.

### Not started / still open

Items **D–H** untouched. Item **C** staged (see below). `bondi_dipole.tex`
edits from items B and the Sec. VI additions are **uncommitted**, and cannot be
compile-checked on this node (no `pdflatex`/`latexmk` present).

---

**Environment note (2026-08-17): the runs moved to a GPU build node, and the
old split between "laptop" and "cluster" no longer applies.** The working tree
now lives on the same machine as the GPUs (4× H100, all four free), with
AMReX, Chombo and the GRTresna solver as siblings of the repo — so edit, build
and run happen in place. There is nothing to push and pull between machines,
and nothing to copy. All campaign output lives under `runs/bondi_rerun/` (see
its `README.md`); the original paper cells sit in `runs/bondi_rerun/published/`.

Two things that bite on this node:

* **Build with the system toolchain, not the wrapper's.** Do *not* source
  `grteclyn-wrapper/scripts/lib/env.sh` before `make`: it prepends the conda
  `grtresna` environment, which hides the system `g++ 11.4` behind a conda
  `gcc 15.2` and leaves no `nvcc` on `PATH`. `make` then reports the target
  *up to date* while compiling nothing. `env.sh` is for running, not building.
  Exact working recipe: `UPSTREAM_MERGE_PLAN.md` §10.
* **Single rank still.** Unchanged from before — launch everything with one
  rank.
* **`env` on this machine is not `env` — spell out `/usr/bin/env`.** Two other
  users' bin directories on this shared box sit ahead of `/usr/bin` on `PATH`,
  and each holds an executable named `env` that is really a PATH-setup snippet
  meant to be sourced (run `type -a env` to see the current set).
  Invoked as `env VAR=x cmd` it prepends its own directory to `PATH`
  and **exits 0 without ever running `cmd`**. LAUNCH.md's documented detach
  pattern (`setsid nohup env ... bash run_*.sh`) therefore launches *nothing*,
  silently, leaving an empty log and a success exit code. Verify any detached
  launch with `pgrep -af bondi_sg_` before believing it started.

The upstream merge from the collaboration is **deferred until after
submission**; see the decision note at the top of `UPSTREAM_MERGE_PLAN.md`.
Nothing on the list below needs it.

---

## 1. Answered: why the max level only triggers at t ≈ 48

Not a schedule — a consequence of the tagging criterion. The run executable
(`Examples/RadialRecipe`) tags with `ChiTagger`
(`Source/Tagging/ChiTagger.hpp`):

```
criterion = dx * sqrt( sum_ij (d_i d_j chi)^2 )  >= regrid_threshold (0.02)
```

i.e. it tags on the **curvature of χ**, not on the matter. The stars are
weak-field objects: |χ − 1| ≲ 0.02 spread over several m⁻¹, so at t = 0 the
criterion sits an order of magnitude below 0.02. `max_level = 1` was allowed
from the start, but nothing crossed threshold until the deepening canonical
well (χ_min 0.98 → 0.44 as the pair nears contact and radiation falls back)
steepened past it — at t ≈ 47.5 in the mixed runs.

Correction (2026-08-17): "never in the singles" was imprecise. The singles do
refine, but only a **fixed, negligible 32 768-cell patch** appearing at
t ≈ 3.4 and never growing — verified identical in July's published `single_p`
and in the fresh reproduction on the GPU node. It is small enough not to
disturb the unigrid reading of the science window, but the tagger is not
inactive there. The genuinely unigrid statement is about the *science window
resolution*, not about level 1 being absent.

Consequence: the entire science window (t ≤ 30) ran **unigrid at Δx = 0.5**,
~10 grid points per stellar rms radius. Both reviews flag this as the top
blocking issue. Fix is item 3-D below (field-based tagging from t = 0).

The paper now states this explicitly (Sec. Setup + Appendix).

## 2. Text fixes applied to bondi_dipole.tex (2026-08-17, no rerun needed)

- **Eq. (3)**: added the missing factor of 2 (`□Φ± = 2 ∂V/∂|Φ|² Φ±`), with a
  clause tying it to the ½-normalized kinetic term; `dV/dφ₀ ≡ 2φ₀ ∂V/∂|Φ|²`
  now defined next to Eq. (5).
- **Eq. (6) sign**: the flipped element was the Π *definition*, not the seed.
  The code convention is `∂ₜΦ = αΠ` (`complex_scalar.py: phase_velocity_pi2`,
  seeds `pi2 = −ω φ₀/α`), so the text now reads `Π = α⁻¹∂ₜΦ` and Eq. (6)
  keeps its minus sign — all three (ansatz, definition, seed) now consistent.
- **New theory paragraph** (Sec. II): the action
  S = ∫√−g[R/16π + ℒ(Φ₊) − ℒ(Φ₋)], Bianchi consistency, inherited strong
  hyperbolicity, and the named quantum-ghost/vacuum-decay objection
  (Carroll–Hoffman–Trodden 2003; Cline–Jeon–Moore 2004).
- **Forward link** added at "held together by something other than gravity …
  its own gravity pushes it apart" (the self-repulsion Forward identified).
- **"First" claims scoped**: "asymptotically flat, constraint-solved,
  dynamically evolved"; static phantom stars (Dzhunushaliev et al. 2014) and
  de Sitter negative-mass bubbles (Mbarek–Paranjape 2014) acknowledged;
  positive-mass-theorem citations (Schoen–Yau, Witten) added to the intro.
- **GPU-first claim** now shared with our own wormhole paper
  (arXiv:2604.00071), cited as companion work `\cite{shirokov2026}`.
- **Conclusion (3) demoted** to a theoretical expectation with the explicit
  P_signed(t) integral stated; abstract softened to match ("should carry …;
  direct momentum-budget measurement deferred").
- **Headline numbers**: 50–150× tied to the single-star noise floor with the
  ≥20× same-sign figure alongside (abstract, results, conclusions); 0.6% →
  "0.6–3.4% across the clean window"; 0.2% → "0.2–3%".
- **Table I** d₀ = 16 row: the wrong-signed canonical −0.34 now printed with a
  footnote explaining the boundary-loss bias.
- **Barycentre defined as an equation** (new Eq., Sec. Setup): field-weighted
  centroid with w_s = (|Φ_s|² + |Π_s|²)^{1/2} over the full domain
  (= accounting region), noted as coordinate/gauge-dependent.
- **Ψ₄ caveat**: recorded but not used quantitatively — both shells (R = 8, 16)
  sit inside the matter.
- **Superposition/phase paragraph** (Sec. III): lump phases 0 and π documented
  (inert for the mixed pair, physical for same-sector pairs); Helfer et al.
  CQG 39, 074001 (2022) cited; phase-offset PP rerun flagged as future work.
- **PM-eq defense**: why matching *ADM* mass (not bare mass/charge) is the
  right choice for Bondi's point-mass limit.
- **"Pull exceeds M/r²"** reworded to the point-mass value M/d² at the
  barycentric separation.
- **Table V** χ row relabeled "extremal χ (min for +, max for −)".
- **Δt = 0.02Δx** acknowledged as deliberately conservative and untuned.
- **Tagging wording** corrected ("χ-gradient" → curvature-of-χ with the exact
  criterion) + the t ≈ 48 explanation in Setup and Appendix.
- **AMR/gauge citations**: CTTK method paper (Aurrekoetxea–Clough–Lim 2023),
  1+log (Bona et al. 1995), Gamma-driver (Campanelli et al. 2006; Baker et al.
  2006), González–Guzmán–Sarbach 2009 (prior nonlinear ghost-scalar NR).
- **Title** trimmed: dropped "with GRTeclyn" (code stays in abstract/intro).
- **Draft box** expanded to 6 items and now points at this file.

Compile verified: 3 passes, zero errors, zero undefined citations/references.

## 3. Blocking items — need GPU reruns (priority order)

Launch commands and env knobs are in `results/bondi-dipole-runaway/LAUNCH.md`
(single rank only — MPI segfaults on the node). Runs now happen in place on
the GPU build node; no commit → push → pull round-trip is needed to get code
onto the hardware.

- **A. Convergence series** (blocks *every* headline number): rerun `PM` and
  one control at N = 192 and 256. The configuration is symmetric under y→−y
  and z→−z, so quadrant symmetry buys ~4× and funds N = 256 at roughly
  current cost (~2–3 h/cell at N = 128 now). Richardson-extrapolate drifts,
  quote convergence order + error bars for 50–150×, 0.6–3.4%, ~10%.
- **B. Momentum-balance diagnostic** — **MEASURED 2026-08-17** (`PM` and
  `PM-eq` reruns to t=60 with `GRTECLYN_SECTOR_DYNAMICS=1`; streams in
  `results/bondi-dipole-runaway/data/pair_pm{,_eqm}_v2/sector_dynamics.dat`,
  analysis in `analysis/momentum_balance.py` → `momentum_balance.csv`).
  Findings, science window t ≤ 30:
  - Sector momenta grow to ±(5–7)×10⁻³ by t=30 and cancel in the signed sum
    to a residual 3–5× smaller — the bookkeeping largely cancels, but
    P_signed(t) is **not** zero and does **not** close to truncation error.
  - The residual is physics, not error: it matches the point-mass sum
    M₊v₊ + M₋v₋ built from the *independently measured* core velocities and
    the dressed ADM masses to 3–10% (t=10: −1.96e-4 vs −1.89e-4; t=30:
    −1.98e-3 vs −2.16e-3). A pair with net mass −0.013 moving at v *must*
    carry P = M_net·v; Forward's "near-zero" is the equal-|M| idealization.
  - Per sector, integral/(M·v_core) climbs 0.72 → 0.95 (canonical) and
    0.81 → 0.94 (phantom) over t=10→30: early on part of each star's cloud
    lags its core; by t=30 the lumps carry their momentum near-rigidly.
  - `PM-eq` does **not** null the residual (−1.72e-3 vs −1.98e-3 at t=30):
    with M_net = 0 the residual is M₊(v₊−v₋), and the gap closes (phantom
    catches up; sep 8 → 7.79 by t=30, → 1.1 by t=60) — a finite-size /
    field-overlap effect absent from the point-particle analysis. Reciprocity
    M₊a₊ = |M₋|a₋ is violated at sep=8 where the lumps overlap (R₉₀ ≈ 5–7.6).
  - Late times are contaminated: past t≈45 the cores approach merger and
    px_canon collapses to ~0 (canonical momentum shed/absorbed); quote
    nothing beyond t=30 (consistent with the published window).
  - Conclusion (3) can be upgraded to a *matter-sector* measurement (done in
    the tex). Full closure — matter + gravitational field = 0 — still needs
    the ADM surface integral, i.e. item C; the volume integral alone
    exchanges momentum with the field by construction.
  - **Bonus, 2026-08-17:** the residual has a closed form in `PM-eq`. With
    M₊+M₋ = 0 the signed mass dipole is D = M₊(x₊−x₋) = M₊d, so Ḋ = M₊ḋ —
    the residual *is* the mass-dipole rate, driven by the closing gap.
    Verified numerically: measured px_total vs M₊ḋ = −1.78e-4/−1.78e-4
    (t=10), −4.46e-4/−4.75e-4 (t=20), −1.72e-3/−2.00e-3 (t=30),
    −1.43e-2/−1.26e-2 (t=45): agreement 1–14%. Two independent
    diagnostics (field volume integral, core tracker) on one number.
    Written into Sec. VI of the tex.

  **Energy-budget / dipole-radiation proposal — assessed 2026-08-17.**
  An external suggestion proposed adding: (i) the radiated energy drives the
  binary's internal energy negative without bound ("Ostrogradsky
  instability"), and (ii) because the gap closes, D̈ ≠ 0 and the pair
  therefore emits **ℓ=1 dipole bremsstrahlung**.
  - **(i) accepted, reframed.** The no-floor argument is right and worth
    stating: Bondi mass loss is geometric (news-squared flux positive-definite
    regardless of matter), the pair starts at M ≈ 0, and the positive-mass
    theorem's DEC hypothesis fails, so nothing bounds the descent — no ground
    state, no terminus. Added to Sec. VI. **Do not call it Ostrogradsky**:
    that theorem concerns higher-derivative Lagrangians, and both sectors here
    have canonical second-order dynamics (the sign enters only the metric
    variation, Sec. II). The correct anchor is the positive-mass theorem
    (already cited, Schoen–Yau/Witten) plus the vacuum-decay cites
    (carroll2003, cline2004) the paper already carries.
  - **(ii) REJECTED — it is wrong, and it would have contradicted the paper.**
    Ψ₄ has spin weight −2, so its multipole expansion *starts* at ℓ=2: there
    is no ℓ=1 harmonic to populate, whatever D̈ does. The tex already says
    this correctly (Sec. VI, "no ℓ=1 mode to populate"); inserting the
    proposal verbatim would have put a flat self-contradiction in the paper
    and handed a reviewer an instant kill. Separately, the proposal misreads
    Ḋ: it is a *momentum*, not a source. For the whole system Ḋ = P_ADM = 0
    for data released at rest, so total D̈ ≡ 0; only the *matter* dipole
    accelerates, and its rate of change is momentum handed to the
    gravitational field — item C, not a radiation channel.
  - What the closing gap *does* contribute is extra time dependence in the
    **quadrupole**, on top of the linear acceleration. Written as such.
- **C. Boundary problem + wave-zone Weyl** (blocks separation-scaling finals
  AND the radiated-energy question raised in the new Sec. VI): implement a
  **sponge layer** (Kreiss–Oliger dissipation profile ramping up near the
  boundary) *and/or* enlarge the domain to L ≥ 128 (ideally 256 so boundaries
  stay causally disconnected past t = 60). Code change in RadialRecipe +
  **rebuild** (GPU build on cluster after pull). Massive-field modes reflect
  perfectly off Sommerfeld boundaries — this is why the window is capped at
  t = 60. Weyl program on the enlarged domain (the near-zone analysis of
  2026-08-17 cannot settle any of these):

  **Enlarged-box cell STAGED 2026-08-17** — `boxC_pm_L128_n256`, queued behind
  `convA_pp_n192` on its GPU. **No code change and no rebuild is needed for the
  domain part**, which is why this could be staged today rather than after a
  sponge-layer implementation:

  | knob | value | why |
  |---|---|---|
  | `BONDI_LFULL` | 128 | wall moves 32 → **64** from centre |
  | `BONDI_NFULL` | 256 | keeps **dx = 0.5**, the published resolution |
  | `BONDI_MAXLEVEL` | 0 | unigrid; 256³ = 58 GB, one GPU, measured |
  | `BONDI_RADII` | `16 24 32 40` | four shells, two of them at R ≥ 32 |
  | `BONDI_STOP_TIME` | 90 | see the causal window below |
  | `BONDI_PSI4_HIGHER_L` | 1 | l = 3,4 stream — only meaningful now |
  | `BONDI_SCRUTINY` | 1 | momentum stream, for the item B closure |
  | `CONSUMER_JOBS` | 8 | mandatory at 2.9 GB/plotfile |

  Same cell count as the `convA_*_n256` cells, so the same 58 GB — but spent on
  *volume* instead of resolution. dt = 0.01 (twice the `convA_*_n256` step), so
  ~27 code units/h and ~3.3 h to t = 90.

  **Why these radii and this stop time.** With the wall at 64, junk from the
  initial data reaches it at t ≈ 64 and returns inward to shell R at
  t ≈ 128 − R. So the clean window per shell is R < t < 128 − R:
  R = 16 → 16–112, R = 24 → 24–104, R = 32 → 32–96, R = 40 → 40–88. All four
  shells are simultaneously clean over **t = 40–88**, which is what a
  multi-shell r⁰ peel test of r·ψ₄ needs. `stop_time 90` buys that window;
  going further only adds time past R = 40's contamination. As a bonus the
  *drift* window also grows — the massive-field bath now needs t ≈ 128 to
  return to the centre instead of t ≈ 64, so this cell can also test whether the
  published t ≤ 30 cap was conservative.

  **The constraint solve needs no change: it already runs on a box of width
  128** (`--grtresna-domain-l 128`, hardcoded in the launcher, 64³ base + 3 AMR
  levels), while the evolution box has been the *smaller* one at 64. Setting
  the evolution box to 128 makes the two coincide, so the initial data covers
  the new box exactly and its configuration stays **identical to the published
  cells** — the enlarged box is then the only variable that changed. What is
  lost is the current 2× margin between solve box and evolution box; solving on
  the evolution domain with a 1/r outer condition is standard practice, but if
  junk turns out to contaminate the shells anyway, the knobs to raise are
  `--grtresna-domain-l 256` with `--grtresna-domain-nx 128` (both exist in the
  CLI; the launcher hardcodes them, so they would need promoting to
  `BONDI_*` env knobs like the rest of §3's rerun set).

  **Verified before staging, not assumed:** the extraction centre is derived as
  `L_full/2` in `campaigns/hq/replay_eval.py` and flows to the consumer's
  `--center`/`--frames-center`, so it becomes 64 64 64 automatically — a
  hardcoded 32 would have silently extracted spheres around the wrong point.
  `--consumer-radii` is `nargs="+"`, so four shells are fine and the `.dat`
  header labels each column group `R=`. Star placement is set by
  `trajectory_lump{k}_R0 = 4` in code units, independent of the box, so the
  separation stays 8.

  **Control not staged yet, deliberately.** The radiated-energy *sign* wants a
  same-sector comparison, but this box configuration has never been run. Verify
  the PM cell past t ≈ 2 (centre at 64 64 64, t = 0 barycentre fingerprints per
  LAUNCH.md §4, four `R=` column groups present) before chaining
  `boxC_pp_L128_n256`.

  - Move the Ψ₄ shells into the wave zone (R ≥ 32, several shells to verify
    the r⁰ peel of r·ψ₄; today's shells at 8/16 show r·ψ₄ falling 5–10×
    between them — pure near zone).
  - **Extend the extraction to l = 3 (and l = 4)** — **DONE 2026-08-17.** New
    module `extraction/psi4_higher_l.py`, its own stream
    `psi4_mode_higher_l.dat`, off by default behind `BONDI_PSI4_HIGHER_L=1`
    (`GRTECLYN_PSI4_HIGHER_L`); the l = 2 column contracts are untouched.
    Choose the multipoles with `BONDI_PSI4_ELLS` (default `3 4`).
    l = 1 does NOT exist for ψ₄ (spin weight −2), which the paper states
    explicitly against the "dipole GW" expectation.

    **Bug found while doing it — read before trusting any old l ≥ 3 work.**
    The l = 3 and l = 4 spin-weighted harmonics already present in
    `consume_plotfiles/sphere.py` were **not normalised**: ⟨|Y|²⟩ ranged over
    0.055–0.56 instead of 1, and they were not mutually orthogonal. They had
    never been exercised because only l = 2 was wired into the extraction path,
    so the error was latent — and item C is precisely the task that would have
    used them, silently producing confident, wrong multipole amplitudes. They
    are now computed from the general Wigner-d formula, verified orthonormal to
    quadrature accuracy for l = 2, 3, 4. The published l = 2 helpers were left
    byte-identical (regression-checked), so no existing number moves.

    One convention note that matters for this item specifically: the legacy
    l = 2 closed forms carry the opposite overall sign on m = ±1 relative to the
    general formula. That cancels in |amplitude| and in power fractions, so
    nothing published changes — but odd-l/even-l *interference* is a relative
    phase measurement, so `psi4_higher_l.py` deliberately uses the general
    harmonic for every l, including l = 2, rather than mixing conventions.
  - **Net radiated energy sign** for the mixed pair (positive, negative, or
    inter-sector cancellation) from wave-zone dE/dt — the headline follow-up
    number; Sec. VI promises it.
  - **ADM mass + momentum surface integrals** at the outer boundary vs time
    (should stay ≈ const; PM-eq should give ADM ≈ 0 globally) — decouples the
    conservation statement from scalar-field centroids.
  - **Bonnor–Swaminarayan benchmark**: compute the B&S acceleration parameter
    for (M₊, M₋, d₀) and compare with the measured early-time coordinate
    acceleration, per run — the strong-field analogue of Table II.
  - Optional figure for the follow-up: shift vector β^x profile along the
    axis (frame dragging of the runaway — matter-free velocity proxy).
- **D. AMR from t = 0** (fixes §1): `ChiTagger` cannot see weak-field stars.
  Add a field-based tagger (tag on scalar amplitude/gradient of both sectors,
  or a fixed-region tagger around the cores) or drop the threshold to ~0.001
  (test for noise-tagging first). Target core resolution Δx ≤ 0.125 — two
  extra levels over the N = 128 base, or one over N = 256. Code change in
  `RadialRecipeLevel::tag_cells` + **rebuild**. Expect the ±8% breathing in
  Table V to shrink; re-quote it.
- **E. Core-weighted centroid rerun**: the diagnostic is "built and waiting"
  (paper, Next steps). Re-run the separation series (d₀ = 12, 16) with it;
  repair the d₀ = 16 canonical −0.34. Add a proper-distance separation
  between χ extrema along the axis as a gauge-robust check of Fig. 4c.
- **F. Noether-charge accounting**: per-sector conserved U(1) charge inside
  spheres — clean core-vs-halo bookkeeping the field-weighted rms stream
  can't give (it tracks the radiation bath).
- **G. Extra controls**: a phase-offset `PP` rerun (current pairs: phases 0/π)
  to bound the superposition systematic, and an **unequal-mass PP** control
  (mixed pair is unequal; an asymmetry-correlated diagnostic bias is
  currently undetectable).
- **H. Robustness point** at a second frequency or coupling (draft-box item).
- **I. CFL experiment — DONE 2026-08-17: REJECTED, keep `dt_multiplier = 0.02`.**
  Ran `single_p` at 0.2 against a fresh 0.02 baseline on the same machine and
  binary (`runs/bondi_rerun/experiments/cfl_dt0p2`, `…/experiments/cfl_base`). The coarse
  step fails the star-health gate outright — rms 5.045 → 13.32 by t = 32 and
  19.22 by t = 40, against a gate of ±10 % — i.e. the star disperses. It also
  drives runaway noise-tagging: the refined region grows from the baseline's
  fixed 32 768 cells to 12 513 280 (380×), taking the cell from 8.8 GB to 41 GB
  of GPU memory. So 0.2 is not a 10× saving, it is an instability being
  resolved into existence.
  Useful side effect: this is the evidence for the "Δt = 0.02Δx is deliberately
  conservative and untuned" caveat — it is now measured, not asserted, and the
  paper can say a 10× coarser step was tried and fails.

Constraint targets after C+D: L∞(H) must stay ≪ 0.05 in the window (currently
0.05 at t ≤ 30, 0.9 by t = 60 — reviewers will reject on that curve; the
constraints figure (`fig_constraints.pdf`, appendix — Fig. 9 since the Weyl
figure took Fig. 8) needs to be re-made from the improved runs).

## 4. Local tasks (CPU, this machine — no GPU needed)

- **R(ω) check**: DONE 2026-08-17 — **the reviewer was right, and the sentence
  was wrong in a second way too.** New script
  `analysis/star_radius_scan.py`, own output `stars/star_radius.csv` (27
  frequencies × both sectors, 0.530–0.900; `--reuse` re-derives the verdict from
  the CSV without redoing the ~15 min of shooting solves). The published
  `star_family.csv` column contract is untouched.
  - **R(ω) is non-monotonic**, as suspected. Canonical: `R_90` falls
    7.567 → 5.108 with a minimum at **ω ≈ 0.80**, then inflates (thick-wall);
    `R_99` turns over much earlier, at **ω ≈ 0.63** (7.57 → 8.30 → 8.30 → 13.40
    at ω = 0.90), because the exponential tail starts spreading while the core
    is still contracting. Phantom is the same shape, `R_99` minimum at ω ≈ 0.67.
  - **But ω = 0.55 and ω = 0.56598 both sit on the shrinking side**, so
    "higher ω → smaller star" is directionally correct. Bounded, though: only
    −32 % in `R_90` (−5 % in `R_99`) is available before the turnover.
  - **The word "compact" was flatly wrong.** `2|M|/R_99` decreases
    *monotonically* across the entire family — 1.456e−2 → 2.282e−3 between
    ω = 0.550 and 0.80 (factor 6.4), factor 14.8 over the full sweep — because
    the ADM mass falls much faster than the radius. Higher ω gives *smaller but
    puffier* stars, i.e. a better separation-to-size ratio bought with a weaker
    drift signal.
  - The other half of the sentence checks out: same-ω sector mass asymmetry
    `|M_-|/|M_+| − 1` shrinks 25.1 % (ω = 0.53) → 20.3 % (0.55) → 4.3 % (0.80)
    → 2.5 % (0.90).
  - Free bonus verification of the campaign's equal-|ADM| pairing: phantom at
    ω = 0.56598 gives |M| = 0.063950 against canonical at ω = 0.550 with
    0.063951 — agreement to 5 significant figures, from an independent code
    path.
  - **Discussion sentence rewritten** in `bondi_dipole.tex` accordingly (no
    longer claims more compact stars; states the turnover, the bounded 32 %,
    the monotonic compactness *decrease*, and that ω cannot substitute for the
    larger domain). **Not compile-checked — there is no LaTeX on the GPU node**;
    only pre-existing macros and the pre-existing `fig:family` label are used.
  - Detector caveat worth knowing if this is re-run: `R_99` is read off an ODE
    *grid point*, so it is quantised and a naive turning-point test fires on
    ~0.01 % wiggles (it first reported a spurious min/max pair at ω = 0.575/0.580,
    displacing the real minimum by 0.04 in ω). `TURN_REL_TOL = 2e-3` suppresses
    that.
- **Figures**: DONE 2026-08-17. `make_article_figures.py` restyled so line
  style, not gray shade, distinguishes runs: Fig. 4b — PM-eq now dash-dot
  (Φ₊) / dash-dot-dot (Φ₋) at 0.35 gray, controls dotted at 0.45 (was 0.78
  light gray); Fig. 4c — four distinct styles (solid / dash-dot / dashed /
  dotted); Fig. 5 — point-mass model dotted (canonical) / dash-dotted
  (phantom) at 0.4 (was 0.78); Fig. 8 — controls dotted / tightly dashed at
  0.45. Captions updated to match; regenerated and visually verified.
  The reviewer's "stray '1' artifacts in Figs. 1–8": **not reproducible** in
  the current PDFs (all eight inspected panel by panel) — likely fixed by the
  earlier "update figure annotations" commit or a viewer artifact on their
  side; re-check once more on the submission build.
- **Error bars**: derive preliminary scatter-based error bars from the four
  control runs for every headline number (superseded by A when it lands).
- **Weyl analysis**: DONE 2026-08-17. New Sec. VI ("The gravitational-wave
  signature") + `fig_weyl.pdf` (new `fig_weyl()` in
  `make_article_figures.py`), from the packed `psi4_mode_l2_all.dat` streams
  (l = 2 spin-weighted modes, R = 8/16). Headline numbers (all reproducible
  by rerunning `fig_weyl()`, which prints them):
  - A_l2 at R = 16 grows 297× (PM), 832× (PM-eq), 481× (MP mirror) above the
    t ≤ 10 floor; PP bounded near 1e-3 with no secular ramp; MM irregular
    (its own spreading).
  - corr(A_l2, |Ẋ|) = 0.995 over 5 ≤ t ≤ 55 (pair-centre speed, ±2 central
    difference); fitted κ = 0.016.
  - Mirror reproduces the PM amplitude history to 4.7% (max-norm).
  - m = ±1 carry 0.6% of the l = 2 power (x-axisymmetric pattern about the
    z-decomposition axis).
  - Spectrum of A_l2 is secular: 99.7% below ω = 0.4; 4e-4 in the ω = 0.55
    band, 8e-6 at 2ω — NOT scalar-bath leakage.
  - Near-zone verdict: A(R=8)/A(R=16) = 5.0 at t = 30, 9.9 at t = 60 →
    r·ψ₄ ~ R⁻²…⁻³, nothing wave-zone; hence no energy/flux claims in the
    paper — that's rerun item C.
  - Physics note recorded in the paper (dipole question): signed source is
    conserved ⇒ signed-dipole current = conserved signed momentum ⇒ no
    mass-dipole radiation at leading order (point-mass reciprocity
    M₊a₊ = |M₋|a₋); ψ₄ (spin −2) has no l = 1 mode at all. Counterpoint to
    the naive "negative masses ⇒ dipole GW" reading of Trivedi–Loeb
    (arXiv:2605.10976), which itself notes universal coupling evades dipole
    bounds.
  - Raw waveform added as `fig_weyl` panel (a) (2026-08-17): monotone
    chirpless ramp in PM vs oscillatory wander in PP/MM at R = 16.
  - **WARNING — wormhole pipeline on this data**: running
    `process_wave/plot_extracted_psi4.py --combined` on
    `data/pair_pm/psi4_mode_l2m0.dat` works mechanically but most panels
    are physically invalid on near-zone data: the strain/LIGO panel reports
    "445 Hz, SNR 133 at 10 Mpc, horizon 166 Mpc" and E_rad = 1.9e-4 M —
    these assume wave-zone radiation and MUST NOT be quoted anywhere; the
    propagation-speed panel returns v = ∞ (no burst front to track — the
    secular ramp peaks at end-of-run at both shells simultaneously); the
    spectrogram is empty (all power below its f = 0.1 floor). Only the
    waveform and PSD panels are meaningful. After rerun item C (wave-zone
    shells), the full pipeline becomes valid — rerun it then, including the
    strain/LIGO comparison the wormhole paper made.

## 5. Bookkeeping for submission

- Remove the red draft box.
- Verify at proof stage: Hoffmann 1966 chapter title/page (1966 Hlavatý
  volume vs the 1965 King's College proceedings version); page numbers of
  Bonnor–Swaminarayan and Bičák et al. New bib entries added today were
  web-verified where possible (dzhunushaliev2014 = "Boson stars with
  nontrivial topology", arXiv:1409.6978); double-check gonzalez2009 II page
  015011 and mbarek2014 101502(R) formatting.
- "traversible" in the Shinkai–Hayward title is original — warn the
  copyeditor not to "fix" it.
- Sec. III D failure narrative **removed from the article** (2026-08-17, user
  decision: debug text, not article text). Full content preserved in §6 below
  and in git history. Cross-references repaired: contribution (v) now points at
  Caveats only; the couplings, dressed-star, and seeding paragraphs keep the
  physics conclusions ("first couplings at which a dressed ω=0.55 star
  exists", "any rescaling de-tunes the eigenstate") without the campaign
  story. The two load-bearing residual facts survive elsewhere: the 30×
  tighter-solve probe lives in the Caveats resolution paragraph; the t = 0
  seed verification lives in contribution (iii).
- 2026-08-17: Sec. VI (GW signature) added with `fig_weyl.pdf`; abstract,
  contributions (iv), roadmap, conclusion (1), and Next steps updated to
  match; Sec. II strengthened (Gleiser–Dotti exponential modes + no-puncture
  argument; Bonnor's inverted fluid interior, verified against the PDF;
  companion wormhole e-folding 0.11 M). New bib entry `trivedi_loeb2026`
  (arXiv:2605.10976, web-verified 2026-08-17); Table V labels shortened to
  fix a 20.6 pt overfull. Figure numbering shifted: Weyl figure is Fig. 8,
  constraints figure is now Fig. 9.
- After reruns A–E land: update Tables I/III/IV, the constraints figure
  (`fig_constraints.pdf`), the clean-window ranges, and re-promote
  conclusion (3); only then consider the full-rewrite
  rule that applies to validated-fix campaigns.

## 6. Removed from the article: the Sec. III D failure narrative (for the record)

Cut 2026-08-17. Kept here because two of the four root causes are the reason
the t = 0 verification gate exists, and future reruns must not regress them.

Before any runaway could be claimed, a single lump had to hold its own size;
four campaigns failed that gate. Heavy lumps collapsed (min χ → 0.009 by
t ≈ 25), light lumps evaporated (rms 5 → 30 by t ≈ 50), and intermediate
lumps did both — indicting the *shape* of the seed, not its weight. Four root
causes, in order:

1. **The amplitude clamp.** The seeding script clamped the amplitude to the
   flat-space thin-wall value √(3λ/4μ), 4.5% below the solved eigenstate: an
   exact Q-ball shape at the wrong height is not a Q-ball — it breathed, shed,
   and evaporated. *Fix:* seed at the solved φ_c.
2. **No dressed equilibrium existed.** The corrected flat-space seed still
   dispersed: a flat-space eigenstate is an equilibrium only in flat space,
   and the dressed family showed the target star did not exist at the original
   couplings. *Fix:* weaken the couplings to the working values and seed the
   dressed profile with its own lapse (paper Eq. 6).
3. **Two seeding paths disagreed silently.** The constraint solve used the
   tabulated star while the evolution's matter reconstruction fell through to
   a Gaussian: the solved metric backed one object, the evolution started
   another. *Fix:* both paths resolve the same star, enforced to 1e−10 by a
   regression test.
4. **The phantom path had never run.** An old veto downgraded exotic lumps to
   canonical matter whenever self-gravity was on — correct for a
   gravity-*bound* star, irrelevant for one bound by self-interaction.
   *Fix:* lift the veto. The first phantom star evolved matched its predicted
   profile, with the predicted mirror-image metric signature (χ just *above*
   unity).

Residual effects attributed rather than fixed: the ±8% radial breathing is
resolution-intrinsic (re-solving with 30× tighter momentum residuals
reproduced it crest for crest); the rms-radius stream tracks the shed
radiation bath, not star health.

The exportable lesson: **verify the seed before believing the evolution.**
Every run was checked at t = 0 against the independently integrated star
solution; two of the four root causes were caught by that check, not by
watching movies.
