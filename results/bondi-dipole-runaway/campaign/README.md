# `campaign/` — every run behind the paper, one directory per cell

This directory holds all 26 evolutions the Bondi-dipole article is built on.
Two families live here side by side, told apart by their names:

* **`convA_*` and `boxC_*`** — the **uniform-grid convergence campaign**
  (2026-08-17/18).  This is the quantitative anchor: every number quoted in
  the article comes from these cells.
* **`pair_*`, `single_*`** (no prefix) — four earlier **production runs at
  Δx = 0.5**, plus two `*_v2` cells carrying the sector-momentum streams.
  These are used *only* where the convergence campaign has no twin: the
  mirror test, the single-star controls, and the momentum bookkeeping.  Every
  figure that touches them says Δx = 0.5 on its face.

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

## Where each cell is used in the paper

Figure and table numbers refer to the compiled `research/bondi_dipole/bondi_dipole.tex`.

| cell | used in |
|---|---|
| `convA_pm_sep12_n256` | Fig. 3 (frames), Fig. 4a/b/c, Fig. 7a, Table I, Table III |
| `convA_pm_sep12_n128`, `_n192` | Fig. 6a (ladder), Fig. 9 (constraints), Sec. V C |
| `convA_pm_n256` | Fig. 4b/c, Fig. 5, Fig. 7a/b, Table I, Table II, Table III |
| `convA_pm_n128` | Fig. 8 (the L = 64 box-independence check), Fig. 6b |
| `convA_pm_n192` | Fig. 6b (refinement sequence) |
| `convA_pm_eqm_n256` | Fig. 7a, Table I, Table II |
| `convA_pm_eqm_n128`, `_n192` | Sec. V C (ladder spread 0.3 %) |
| `convA_pm_sep16_n*` | Fig. 6b, Table I, Table III — always as a bound |
| `convA_pp_n*` | Fig. 4b, Fig. 5, Fig. 6a/b, Fig. 7a, Fig. 9, Table I |
| `convA_mm_n*` | Fig. 4b, Fig. 5, Fig. 6b, Fig. 7a, Fig. 9, Table I |
| `boxC_pm_L128_n256` | Fig. 8, Sec. VI B |
| `boxC_pp_L128_n256` | `analysis/wave_check.md` only |
| `pair_pm` | Fig. 5, Fig. 7a (same-numerics mirror reference) |
| `pair_mp_mirror` | Fig. 5, Fig. 7a, Table I |
| `single_p`, `single_m` | Table I and Sec. V F (single-star stability records) |
| `pair_pm_v2`, `pair_pm_eqm_v2` | Sec. VI C and Discussion item (4) — the signed-momentum balance |

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
| `sector_dynamics.dat` | the `*_v2` cells only — signed sector momenta |
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
`research/bondi_dipole/pack_campaign.sh` (owns `convA_*`/`boxC_*`) and
`pack_results.sh` (owns the unprefixed cells); both write here, under disjoint
names.
