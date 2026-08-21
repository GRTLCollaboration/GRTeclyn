# `campaign/` — every run behind the paper, one directory per cell

This directory holds all 32 evolutions behind the Bondi-dipole work.  Three
families live here side by side, told apart by their names:

* **`convA_*` and `boxC_*`** — the **uniform-grid convergence campaign**
  (2026-08-17/18).  This is the quantitative anchor: every number quoted in
  the article comes from these cells.
* **`pair_*`, `single_*`** (no prefix) — four earlier **production runs at
  Δx = 0.5**, plus two `*_v2` cells carrying the sector-momentum streams.
  These are used *only* where the convergence campaign has no twin: the
  mirror test, the single-star controls, and the momentum bookkeeping.  Every
  figure that touches them says Δx = 0.5 on its face.
* **`next_*`** — the **follow-up campaign** (2026-08-20), six cells that answer
  questions the published campaign could not.  **Nothing in the current article
  depends on these**; they exist to close the open items in
  `research/bondi_dipole/docs/GPU.md`, and their results are summarised below.

Every cell is **strictly uniform-grid** (`max_level = 0`): there is no mesh
refinement anywhere in this campaign, so a cell's resolution is fully
described by its grid spacing.

> These two families used to live in separate `campaign/` and `data/` trees.
> They were merged here on 2026-08-20 and `data/` was deleted, along with five
> cells the `convA_*` ladder had superseded (`pair_pp`, `pair_mm`,
> `pair_pm_eqm`, `pair_pm_sep12`, `pair_pm_sep16`).  The merge was verified by
> regenerating every figure and every printed number, unchanged.  (The paper's
> 2026-08-20 shortening pass later dropped two figures and merged two tables;
> the mapping below reflects the current nine-figure, four-table layout.)

---

## The reference cell

**`convA_pm_sep12_n256`** is the reference configuration for the whole paper:
the mixed canonical–phantom pair released at rest at x = ±6 (d₀ = 12), on the
finest grid, Δx = 0.25.

It is *not* the loudest cell — the d₀ = 8 runs drift four times further — but
it is the cleanest.  The dressed stars have rms radius ≈ 5 and R₉₀ ≈ 7.6, so
at d₀ = 8 the two envelopes overlap from t = 0 and the pair reaches near
contact (gap 1.86) by t = 60; that is a merger, not a two-body force law.  At
d₀ = 12 the gap only closes from 12.00 to 10.14 over the whole run, and the
drift still sits far above every control.  See Sec. IV of the article.

---

## The convergence campaign — `convA_*` (L = 64)

Six configurations × three resolutions.  Name pattern
`convA_<config>_n<N>`, with L = 64 throughout:

| N | Δx | points per stellar rms radius |
|---|---|---|
| 128 | 0.50 | ≈ 10 |
| 192 | 0.33 | ≈ 15 |
| 256 | 0.25 | ≈ 20 |

| config | sectors | d₀ | what it is |
|---|---|---|---|
| `convA_pm_sep12_*` | canonical + phantom | 12 | **the reference pair** |
| `convA_pm_*` | canonical + phantom | 8 | strong-signal illustration; contact-contaminated |
| `convA_pm_eqm_*` | canonical + phantom | 8 | equal-\|M\| variant, phantom retuned to ω = 0.56598 |
| `convA_pm_sep16_*` | canonical + phantom | 16 | widest pair; **not grid-stable**, quoted as a bound only |
| `convA_pp_*` | canonical × 2 | 8 | null control: same-sign pair must not drift |
| `convA_mm_*` | phantom × 2 | 8 | null control: same-sign pair must not drift |

The point of running all three grids is Sec. V C: a physical drift should be
insensitive to the grid, and a numerical one should converge away.  Both
happen — the mixed-pair drifts sit on plateaus (spread 6.3 % at d₀ = 12,
0.5 % at d₀ = 8) while the `pp`/`mm` residuals fall by factors of 3–4.  The
`sep16` family is the exception that shows the test has teeth: its 53 % spread
is why the article quotes d₀ = 16 as a bound rather than a measurement.

## The double box — `boxC_*` (L = 128)

| cell | what it is |
|---|---|
| `boxC_pm_L128_n256` | mixed pair, d₀ = 8, domain doubled to L = 128 at Δx = 0.5, extraction shells at r = 16/24/32/40 |
| `boxC_pp_L128_n256` | the same for the `PP` control (diagnostic only; no figure uses it) |

Doubling the *box* rather than refining the grid is what buys wave-zone
extraction: it is the light-crossing time to the boundary that matters, not
the resolution.  Inside the retarded-time window u = 10–25 the four shells
agree on one outgoing wave to 1.1–1.4×.  See Sec. VI B.

## The Δx = 0.5 production cells

| cell | what it is | why it survives |
|---|---|---|
| `pair_pm` | mixed pair, d₀ = 8 | the same-numerics partner the mirror is compared against |
| `pair_mp_mirror` | the same pair with the sectors **swapped in space** | the mirror control; the campaign has no mirrored twin |
| `single_p` | one canonical star alone | noise floor + stability record |
| `single_m` | one phantom star alone | noise floor + stability record |
| `pair_pm_v2` | mixed pair, d₀ = 8 | carries `sector_dynamics.dat` — the only sector-momentum stream |
| `pair_pm_eqm_v2` | equal-\|M\| pair | the same, for the equal-\|M\| case |

A mirror test has to hold *everything* fixed except the reflection, so it is
correctly run at the resolution its partner was run at — which is why this
comparison lives at Δx = 0.5 rather than on the finest grid.

---

## The follow-up campaign — `next_*` (2026-08-20)

Six cells, ~33 GPU-hours, all reaching their full stop time.  They answer four
open questions; two answers came out the opposite of what was expected — and
the longest cell turned up a fifth thing nobody asked for, which is the most
important result here and is written up as §6 below.  **Read §6 before quoting
§1–§5.**

| cell | d₀ | N | L | stop | what it tests |
|---|---|---|---|---|---|
| `next_tolB_n128` / `_n192` / `_n256` | 12 | 128/192/256 | 64 | 60 | does a *tighter initial-data solve* buy a convergence order? |
| `next_spongeA_boxC` | 8 | 256 | 128 | 90 | is the boundary sponge safe? (twin of `boxC_pm_L128_n256`) |
| `next_spongeB_sep16_long` | 16 | 256 | 128 | **120** | the constant-gap runaway, past the old t = 60 wall |
| `next_lever_omega0615` | 8 | 256 | 64 | 60 | force scaling at a wider mass contrast (ω = 0.615) |

### 1. The tolerance ladder — the hypothesis was wrong

Every rung of the published ladder stopped on the *same* solver tolerance and
so floored at the same `0.0832 %` Hamiltonian violation, and Sec. V C blamed
that floor for the missing convergence order.  These three cells scale the
tolerance as Δx⁴ instead, and the solve does refine as intended:

| rung | exit tolerance | Newton steps | achieved |
|---|---|---|---|
| `next_tolB_n128` | 0.1 | 7 | 0.0832 % |
| `next_tolB_n192` | 0.019 | 9 | 0.0134 % |
| `next_tolB_n256` | 0.00625 | 10 | **0.0054 %** |

**It changed nothing.**  Fifteen times cleaner initial data moved the drift by
≈ 0.1 % and left the grid spread at 6.3 %, exactly where the published ladder
sits.  The initial-data tolerance was never the limiting factor.

### 2. The peak tracker found the real culprit — the diagnostic

The same three cells carry `sector_dynamics.dat`.  Comparing the two ways of
locating the pair on the identical ladder:

| diagnostic | spread at t = 20 | at t = 60 |
|---|---|---|
| barycentre (what the article quotes) | 12.5 % | 6.3 % |
| halo-free core tracker | 3.3 % | 5.0 % |
| **core-to-core separation** | **0.21 %** | **1.68 %** |

The barycentre integrates over the whole domain and drags in the radiation
bath; the tracker does not.  The separation is a genuinely well-converged
quantity.  A single convergence *order* still does not fall out — some times
give ≈ 4, others do not, and one is non-monotonic — because the drift
differences sit near the tracker's own noise floor at 301 samples.

### 3. The sponge is safe, and buys ~20 time units

`next_spongeA_boxC` is `boxC_pm_L128_n256` with the sponge on and nothing else
changed, so the two are directly subtractable:

| t | negative energy in the box, sponge ON | OFF | change |
|---|---|---|---|
| 20 | 0.6145 | 0.6147 | −0.03 % |
| 60 | 2.1663 | 2.1684 | −0.10 % |
| 90 | 2.2791 | 2.2821 | **−0.13 %** |

The Hamiltonian constraint is unchanged to five decimals at *every* time.  The
sponge removes phantom radiation and does not perturb the solution.

### 4. The constant-gap runaway — the gap holds, but the star does not

`next_spongeB_sep16_long` runs d₀ = 16 to t = 120.  The coordinate gap does hold
— it closes by only 2.8 % out to t = 80 — but the object whose acceleration this
cell was built to measure turns into a black hole partway through, and two
further checks say the surviving motion is not the textbook runaway either:

| t | pair drift | coordinate separation | proper separation | canonical peak |
|---|---|---|---|---|
| 0 | 0.00 | 16.00 | 14.99 | 0.0226 |
| 20 | 0.01 | 15.94 | 16.96 | 0.0227 |
| 40 | 0.05 | 15.81 | 16.82 | 0.0233 |
| 45 | 0.11 | 15.78 | — | 0.0237 |
| 60 | 0.44 | 15.68 | 20.74 | 0.0290 |
| 80 | 1.79 | 15.56 | 23.08 | **0.0000** |

Three separate problems, in order of severity:

* **The canonical star collapses.**  Minimum lapse 0.98 → 0.047 and χ → 0.024;
  its peak amplitude first *rises* 28 % as the core squeezes down, then goes to
  zero by t ≈ 75.  Everything past t ≈ 45 is a collapsing star, not a
  self-accelerating one.  §6 shows this is intrinsic to the star, not caused by
  the pairing.
* **"Constant gap" is a coordinate statement only.**  Proper separation grows
  15.0 → 23.1, a 54 % increase, over the same window in which the coordinate
  separation looks flat.
* **The acceleration has the wrong shape.**  The phantom's speed roughly doubles
  every 10 time units (0.0018, 0.0055, 0.0157, 0.0292, 0.0504, 0.0906 over
  successive windows), taking off exactly as the canonical star dies.  This is
  *not* coordinate drag — the shift at the phantom is 0.003–0.005 against a
  coordinate speed reaching 0.12 — but a genuine runaway at fixed separation
  needs speed rising *linearly*, not compounding.

**Defensible window at this separation: t ≲ 45** — but d₀ = 16 is the worst case
for signal.  At d₀ = 12 the same tracker *does* see the Bondi signature: both
stars move +x with displacements in the ratio 1.13 / 1.16 / 1.25 at t = 40 / 45
/ 50, against the 1.205 the mass ratio predicts, while the gap holds to 1 %.
The barycentre the article quotes gives 7.11 / 4.46 for the same two stars and
converges on nothing — the signature is real, and invisible to the published
diagnostic.  What still does not match is the time dependence: the drift grows
as t³·⁹, not the t² of a constant-acceleration runaway.

`movies/` in that cell shows all of it: matter (`scalar_activity_*`, `phi_z`),
geometry (`chi_z` — the collapse), negative energy (`rho_req_z`) and
`local_speed_z`, 301 frames each at 12 fps over the whole t = 0…120 run.  On
`movie_scalar_activity_z.mp4` the right-hand (canonical) blob visibly brightens
and then vanishes while the left-hand (phantom) one is still barely moving.

### 5. Force scaling is not monotonic — because all three pairs have merged

| \|M₋\|/M₊ | cell | pair drift at t = 60 | separation 8.00 → |
|---|---|---|---|
| 0.62 | `next_lever_omega0615` | 5.66 | 0.84 |
| 0.83 | `convA_pm_n256` (published) | **6.41** | 1.86 |
| 1.00 | `convA_pm_eqm_n256` | 6.23 | 1.42 |

The lever arm is now wide, and the drift peaks near the published ratio rather
than rising with contrast — but the last column is why.  All three pairs have
essentially merged by t = 60.  What is being compared is merger dynamics, not a
force law, so the non-monotonicity should not be quoted as a scaling result.

### 6. The canonical star was never stable — it collapses on its own

This was not on the queue.  It came out of running one cell past t = 60 for the
first time, and it applies to **every cell in this directory**, published and
follow-up alike.

*Each star run alone, no partner anywhere in the box* (`single_p`, `single_m`):

| alone | lapse t = 0 → 40 | rms radius t = 0 → 40 |
|---|---|---|
| canonical | 0.976 → 0.867, falling monotonically | 5.05 → 13.40 |
| phantom | 1.001 → 0.998, flat | 5.43 → 6.32 |

The canonical star destabilises with nothing near it, and its lapse at t = 40
alone (0.867) matches the d₀ = 16 pair (0.824) and the d₀ = 8 pair (0.853) — the
partner changes almost nothing.

*Stacking the deck one way or the other:*

| pair | worst lapse | worst χ | outcome |
|---|---|---|---|
| canonical + canonical (`*_pp_*`) | 0.014 | 0.003 | horizon almost at once (lapse < 0.5 by t = 4) |
| canonical + phantom (`*_pm_*`) | 0.047 | 0.024 | horizon at t ≈ 57–75 |
| phantom + phantom (`*_mm_*`) | 0.974 | 0.989 | nothing happens, ever |

More canonical matter collapses faster; more phantom matter is more stable; two
phantoms never collapse at all.  Within the mixed runs the same ordering holds —
d₀ = 16 crosses lapse < 0.5 at t = 57 but d₀ = 8 not until t = 62, so bringing
the negative mass *closer* delays the collapse.

*And `../stars/star_radius.csv` said so before any of this was run.*  Along the
canonical sequence, ADM mass falls as the star gets denser at the frequency
every cell uses:

| ω | φ_c | M_ADM | dM/dφ_c | branch |
|---|---|---|---|---|
| 0.900 | 0.0159491 | 0.007963 | — | stable |
| 0.800 | 0.0191264 | 0.011209 | +2.4 | stable |
| 0.775 | 0.0194775 | 0.012612 | +4.0 | stable |
| **0.550** | **0.0196947** | **0.063951** | **−182** | **unstable ← used in every cell** |
| 0.615 | 0.0199119 | 0.035175 | −110 | unstable |

A negative slope is the standard turning-point signature of the unstable branch.
The stable branch is at ω ≈ 0.775–0.90, where the stars are 3–5× lighter and
~50 % larger.  The phantom sits on the equivalent negative-slope branch and is
nonetheless perfectly stable — which is the physical point: a negative-energy
object's own gravity pushes *outward*, so it has no collapse channel.  The
instability is one-sided by construction.

**One more geometric fact worth checking before quoting any separation.**  Each
star's 99 % radius is 8.72 (canonical) and 9.16 (phantom) at ω = 0.55:

| d₀ | centres apart | sum of r₉₉ | sum of r₉₀ | verdict |
|---|---|---|---|---|
| 8 | 8 | 17.9 | 15.8 | centres closer than *one* radius — a single blob |
| 12 | 12 | 17.9 | 15.8 | heavily overlapping |
| 16 | 16 | 17.9 | 15.8 | envelopes touch; cores just clear |

d₀ = 8 — the configuration behind the headline drift — is not two objects, and
does not hold a gap: it merges, 8.00 → 1.81 by t = 60 and through zero by t = 70.
d₀ = 16 is the first separation where the bulk of the two stars are apart, and
there the drift at t = 60 is 0.44 against 7.11 at d₀ = 8 — a factor 16 smaller,
where an inverse-square force between separated masses predicts a factor 4.

Full write-up, including what it costs the article and the three ways forward:
`../../../research/bondi_dipole/docs/GPU.md`.

---

## Where each cell is used in the paper

Figure and table numbers refer to the compiled `research/bondi_dipole/bondi_dipole.tex`.

| cell | used in |
|---|---|
| `convA_pm_sep12_n256` | Fig. 3 (frames), Fig. 4a/b/c, Fig. 7a, Table I, Table III |
| `convA_pm_sep12_n128`, `_n192` | Fig. 6a (ladder), Fig. 9 (constraints), Sec. V C |
| `convA_pm_n256` | Fig. 4b/c, Fig. 5, Fig. 7a/b, Table I, Table II, Table III |
| `convA_pm_n128` | Fig. 8 (the L = 64 box-independence check), Fig. 6b, Sec. VII (box comparison against `boxC_pm_L128_n256`) |
| `convA_pm_n192` | Fig. 6b (refinement sequence) |
| `convA_pm_eqm_n256` | Fig. 7a, Table I, Table II |
| `convA_pm_eqm_n128`, `_n192` | Sec. V C (ladder spread 0.3 %) |
| `convA_pm_sep16_n*` | Fig. 6b, Table I, Table III — always as a bound |
| `convA_pp_n*` | Fig. 4b, Fig. 5, Fig. 6a/b, Fig. 7a, Fig. 9, Table I |
| `convA_mm_n*` | Fig. 4b, Fig. 5, Fig. 6b, Fig. 7a, Fig. 9, Table I |
| `boxC_pm_L128_n256` | Fig. 8, Sec. VI B, Sec. VII (box comparison and the peak-tracker cross-check) |
| `boxC_pp_L128_n256` | `analysis/wave_check.md` only |
| `pair_pm` | Fig. 5, Fig. 7a (same-numerics mirror reference) |
| `pair_mp_mirror` | Fig. 5, Fig. 7a, Table I |
| `single_p`, `single_m` | Table I and Sec. V F (single-star stability records) |
| `pair_pm_v2`, `pair_pm_eqm_v2` | Sec. VI C and Discussion item (3) — the signed-momentum balance; `pair_pm_v2` also carries the peak tracker (`core_x_*`) that bounds the halo bias in Sec. VII |

Table IV (computational cost) is aggregate: it comes from the run logs, not
from any one cell's streams.

---

## What is inside a cell

Every `.dat` stream carries its own `#` header naming each column; the data
dictionary is in [`../README.md`](../README.md).  The streams the article
actually reads are:

| file | read by |
|---|---|
| `sector_barycenters.dat` | every drift, gap and velocity number in the paper |
| `constraint_norms.dat` | Fig. 9 and Appendix A (the composite norms are the `L2_Ham_amr`/`L2_Mom_amr`/`Linf_Ham_amr` columns named in the `#` header — cols 1/2 after `time` are a different, level-0-only diagnostic) |
| `psi4_mode_l2_all.dat` | Fig. 7, the ℓ = 2 Weyl amplitudes at R = 8 and 16 |
| `psi4_mode_l2m0.dat` | Fig. 8, the wave-zone shells (`boxC` header names r = 16/24/32/40) |
| `confinement.dat` | extremal χ and rms radii (last column is `min_chi`) |
| `sector_dynamics.dat` | the `*_v2` and `boxC_*` cells — signed sector momenta, and the `core_x_*` peak tracker Sec. VII compares against the barycentre |
| `metadata.json`, `evolution_params.txt`, `grtresna_params.txt` | provenance: what was solved, how well, what was evolved |

`convA_*`/`boxC_*` cells additionally carry `frames/` — matter
(`scalar_activity_z`) and geometry (`chi_minus_1_z`) on the z slice, one every
Δt = 10 plus the final state, named by simulation time.

## Regenerating

```bash
# derived tables (analysis/)
python3 analysis/{make_tables,momentum_balance,separation_scaling}.py
python3 analysis/{newtonian_reference,constraint_check,convergence_check,wave_check}.py

# the article's figures
cd research/bondi_dipole && SCRATCH=/tmp python3 make_article_figures.py
```

The figure script also prints, to stdout, every number quoted in the article's
captions and text.  Re-packing from the run tree is
`research/bondi_dipole/pack_campaign.sh` (owns `convA_*`/`boxC_*` and the
`next_*` follow-up cells) and `pack_results.sh` (owns the unprefixed cells);
both write here, under disjoint names.

The follow-up cells are launched by
`grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_next_campaign.sh`, and the
movies above by `grteclyn-wrapper/scripts/plot/make_movies.sh`.
