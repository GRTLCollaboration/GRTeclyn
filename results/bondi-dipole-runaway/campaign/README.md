# Campaign extracts — the Bondi dipole runaway

Publication-level data only. One folder per cell. Movies, plotfiles and the
550 MB initial-data file stay in the gitignored run tree at
`runs/bondi/staging/`, because they carry no number an analysis reads. Cells
that were run with frames carry a `frames/` folder holding one still per field
every `dt = 10` — enough to show what the movie shows, without the ~250 frames
per field the movie is made from. The `t = 400` cell is the exception: it exists
for one moving picture, is twice as long as every other cell, and no number in
the analysis is read off its stills, so only its movie was kept.

Packed by `research/bondi_dipole/pack_runaway.sh`, which is safe to re-run at
any time and scrubs absolute paths.

**The result these cells support.** A positive-mass canonical boson star paired
with a negative-mass phantom star accelerates as a unit, at constant rate, with
zero total momentum and constant separation — the Bondi runaway — and the rate
falls off as one over the separation squared. Nothing here is corrected,
subtracted or fitted around an artefact.

## How to read a cell name

`runaway_pair_d16_L64_N128_lev0` = the runaway pair, stars `16` apart, evolved
in a box of side `L = 64` divided into `N = 128` cells per side, `lev0` = no
mesh refinement. Cell size `dx = L/N = 0.5` everywhere, for the whole run.

Two boxes are involved and they are easy to confuse:

| | box side `L` | cells `N` | refinement | cell size |
|---|---|---|---|---|
| **evolution** — what the name records | `64` | `128` | none | `0.5` |
| **elliptic solve** — 32 MPI ranks | `128` | `256` | none | `0.5` |

**The two cell sizes are equal, and that is the whole point.** When they differ,
the solved metric lands on the evolution grid displaced by a fraction of a cell
while the matter is repainted exactly, so every star is born off the centre of
its own gravitational well — canonical falls back down it, phantom is pushed up
it. That artefact was the same size as the signal being measured. Matching the
spacings makes the transfer a straight copy, and the metric-minus-matter offset
is `0.0000` on every cell here.

Note that `evolution_params.txt` records `N_full` and `max_level` but **not** the
evolution `L`; the only `L` in it is `128`, which belongs to the solve box.
`launch_config.sh` is the authoritative record of every grid setting, which is
why the cell name carries all three.

## Cells

Shared by every cell: matter rung `m = 1`, `λ = 10240`, `μ = 21845333`;
canonical star at `ω = 0.75` (`M = +0.014350`); phantom at `ω = 0.7603`
(`M = −0.014350`, mass-matched); maximal slicing.

### Controls — does a lone star sit still?

One star, alone in an empty box, nothing to pull on it.

| cell | star | drift over `t = 200` | peak amplitude |
|---|---|---|---|
| `control_lone_canonical_L64_N128_lev0` | canonical, box centre | `1.8e-03` | `−0.48%` |
| `control_lone_phantom_L64_N128_lev0` | phantom, off-centre at `x = 37` | `1.6e-03` | `−0.48%` |

Drift is the largest displacement on any axis. The identical canonical run on
the old initial data moved `−0.328` on all three axes — a `269×` reduction, and
what remains is no longer diagonal, so it is noise rather than drift. The
off-centre phantom is the sharper test: no symmetry forces it to stay put.

### The runaway pairs — the result

`a` is fitted to the pair's midpoint from `t ≥ 5`, before which the gauge is
still settling. `a·d²` should return the star's mass at every separation.

| cell | `d` | `a` measured | `a·d²` | separation at `t = 200` | total momentum |
|---|---|---|---|---|---|
| `runaway_pair_d08_L64_N128_lev0` | 8 | `2.3210e-04` | `0.01485` | `8.187` | `5.8e-05` |
| `runaway_pair_d10_L64_N128_lev0` | 10 | `1.4541e-04` | `0.01454` | `10.107` | `3.7e-05` |
| `runaway_pair_d12_L64_N128_lev0` | 12 | `1.0003e-04` | `0.01440` | `12.063` | `2.3e-05` |
| `runaway_pair_d16_L64_N128_lev0` | 16 | `5.5968e-05` | `0.01433` | `16.025` | `9.8e-06` |

True `GM = 0.014350`. A power law through the four accelerations gives
`d^−2.051` against `−2` exact. `d = 10` is the headline cell: its acceleration
is constant across four disjoint fits of the run
(`1.473 / 1.442 / 1.456 / 1.469 e−04`).

The trend in `a·d²` is physical and should not be hidden — the closest pair
overshoots by `3.5%` and the agreement tightens monotonically with separation,
which is where finite-size and post-Newtonian corrections both push.

### The resolution ladder — does the runaway survive refinement?

Same physics, same box, three cell sizes. `a` is fitted to the pair's midpoint
over the last two thirds of the run.

| cell | cells | drift at `t = 200` | separation | `a` |
|---|---|---|---|---|
| `runaway_pair_d10_L64_N128_lev0` | 128³ | `+2.8815` | `10.000 → 10.003` | `1.463e-04` |
| `runaway_pair_d10_L64_N192_lev0` | 192³ | `+3.0139` | `10.000 → 9.930` | `1.611e-04` |
| `runaway_pair_d10_L64_N256_lev0` | 256³ | `+3.0016` | `10.000 → 9.915` | `1.596e-04` |

The two finest rungs agree to `0.4%` on drift and `0.9%` on acceleration; the
coarsest sits about `4%` low. That is a converging quantity, and it is the
opposite of what a discretisation artefact does — refining does not drive the
effect toward zero, it settles it onto a value. **The N=256 rung is the
headline number**, with N=192 as its convergence partner.

**Error bars, honestly.** The triple is non-monotone (N=192 lands slightly
above N=256), so no formal convergence order is quotable. The cause is visible
in `constraint_norms.dat`: the `t = 0` violation *rises* with resolution
(`1.1e-04 → 1.5e-04 → 1.8e-04` — initial-data interpolation noise scales as
`1/dx²`) while the evolution's own violation *falls* (`1.1e-05 → 6.6e-06 →
5.6e-06`). Two error sources with opposite grid behaviour never leave a single
asymptotic power. Richardson on the fine pair alone gives corrections smaller
than the pair's own spread at any assumed order, so the quote is
**drift `3.00 ± 0.01`, `a = (1.60 ± 0.02)e-04`**, the fine-pair spread as the
error bar.

Total momentum sharpens the claim rather than diluting it: `px_total` is
`3.7e-05` at N=128 and `3.8e-06` at N=256 — the pair displaces by 3 while its
momentum converges *toward* zero. Displacement without momentum is the Bondi
signature.

One nuance the coarse grid hides: N=128 says the pair is rigid (separation
back to `10.003` at `t = 200`), but both fine rungs agree it actually closes
slightly (`9.930` / `9.915`). Rigidity is a `1%`-level statement; the
fine-grid value is the honest one.

Two independent checks agree: doubling the box moves the drift by `4%`
(`wavezone`, below), and running four times longer leaves the acceleration flat
(`longrun`, below).

### The mirror — swap the sectors, the runaway reverses exactly

`control_mirror_mp_d10_L64_N128_lev0` is the headline cell with the two stars
swapped in place (phantom left, canonical right). Everything must come out
sign-flipped and nothing else may change:

| | headline | mirror | ratio |
|---|---|---|---|
| drift | `+2.88147` | `−2.88153` | `−1.000022` |
| `a` | `+1.46336e-04` | `−1.46340e-04` | `−1.000025` |

A reversal exact to two parts in 10⁵. No grid artefact, boundary effect or
gauge drift knows which side the phantom is on; only the physics does.

### Sustained acceleration — is it steady or is it running away with itself?

| cell | `t_end` | drift | separation | `a` (last ⅔) | `a` (final ¼) |
|---|---|---|---|---|---|
| `longrun_pair_d10_t400_L64_N128_lev0` | `400` | `+11.5177` | `10.000 → 11.741` | `1.418e-04` | `1.444e-04` |

The acceleration is constant to `2%` between the middle of the run and its final
quarter. A late-time uptick visible at `t = 200` does not survive to `t = 400`,
so it was a transient.

One thing does change over the longer window: the separation *opens*,
`10.000 → 11.741`, having stayed flat all the way to `t = 200`. The acceleration
is unaffected, but any statement that the pair moves rigidly belongs to
`t ≤ 200` only.

### The same-sign nulls — the pairs merge, and the centroid still does not move

`control_pair_pp_*` and `control_pair_mm_*` put two stars of the *same* mass
sign side by side. Their report cards say *matter dispersed, only 15–17% still
confined*, and an earlier version of this section blamed a wave from the outer
boundary arriving at `t ≈ 32` and declared the cells unreadable after that.
**That was wrong.** The boundary flux measurement retires it: net inward flux in
the PP cell is at round-off for the whole run — nothing arrives from the wall,
the sponge absorbs as designed — and the only significant flux is *outward*,
near `t = 100`, which is the story leaving the box, not entering it.

What the frames of the PP N=256 cell actually show, tracking the two
gravitational wells image by image (the full series is packed as that cell's
`well_tracking.dat`; first merged frame `t = 33.6`):

| `t` | 0 | 13 | 20 | 26 | 32 | ~35 |
|---|---|---|---|---|---|---|
| gap between wells | `8.75` | `8.25` | `8.00` | `7.25` | `5.25` | **merged** |

Two positive masses attract, fall together, and **merge at `t ≈ 35`**. The ×7
activity growth (onset `t ≈ 40`, peak `t ≈ 97`) is merger ejecta, which then
drains out through the sponge; the surviving field peak is the remnant core
ringing down — the collide-and-bounce visible in the chi movie. Its flatness
across resolution (`7.1× / 7.0× / 6.9×` at N = 128/192/256) is what converged
physics looks like.

| | runaway (mixed) | null (same-sign) |
|---|---|---|
| total scalar field in box | `7.9 → 9.2` | `7.8 → 54.7 → 27.6` (ejecta, then drain) |
| **peak** field strength | `0.0247 → 0.0246` | `0.0247 → 0.0233` (remnant) |
| net inward boundary flux | round-off | round-off |
| constraint error | `3.7e-06 → 6.6e-06` | `2.1e-05 → 2.3e-05`, bounded |

**The null verdict is strengthened, not weakened.** Over the full `t = 200` —
through infall, merger and ringdown — the pair's centroid moves by:

| grid | PP centroid drift | MM centroid drift | min χ (PP / MM) |
|---|---|---|---|
| 128 | `+0.00073` | `−0.00026` | `0.97947` / `1.00000` |
| 192 | `+0.00048` | `−0.00035` | `0.97920` / `1.00000` |
| 256 | `+0.00052` | — | `0.97914` / — |

Four orders of magnitude below the runaway's `+2.88`. A configuration that
cannot manufacture net momentum even while merging is a stronger null than one
that merely sits still.

**One measurement the same-sign cells do not contain.** The sector splitter
assigns matter by field sign, so both stars of a same-sign pair land in one
sector: the tracker reports a single core at the pair midpoint and
`coord_sep = nan`. Per-star trajectories exist only where frames do — PP at
N=256. For MM, which shows the same ×7 growth with the same timing, whether the
two phantom stars merged like the PP pair or pushed apart as Bondi predicts for
two negative masses is **not yet measured**; the flux diagnostic is blind there
too (it integrates canonical-sector energy and reads zero in a phantom-only
box). The `control_pair_mm_d10_L64_N128_lev0_frames` re-run (physics identical
to the archived MM cell, frames on) exists to answer it.

What does survive from the old caveat is a statement about `t = 0`, not the
evolution: the elliptic solve's outer boundary condition is genuinely wrong for
a same-sign pair (the box carries net `2M`), which is why the solve residual
floors near `5.4e-04 %` and the gate is held flat at `0.002` on every same-sign
rung.

### The wave zone — a null with a number

`wavezone_pair_d10_L128_N256_lev0` doubles the box to `L = 128` at the same cell
size, with extraction shells at `R = 16/24/32/40`. Its drift, `+2.7606`, is
within `4%` of the `L = 64` cell, so the boundary is not driving the runaway.

`r·ψ₄` is meant to be flat in the wave zone. At `t = 199`:

| `R` | 16 | 24 | 32 | 40 |
|---|---|---|---|---|
| `ψ₄` (l=2) | `3.23e-05` | `6.52e-06` | `2.17e-06` | `7.99e-07` |
| `r·ψ₄` | `5.16e-04` | `1.57e-04` | `6.93e-05` | `3.20e-05` |

It falls by a factor of `16`. A power-law fit gives `ψ₄ ∝ r^−4.0`; radiation
would give `r^−1`. Every shell in this box is in the near zone, and no radiative
tail is measurable out to `R = 40` — consistent with the setup's own prediction
that a pair with zero total momentum has a static mass dipole and therefore no
`l = 1` growth. Report it as a null with a number, not as a missing measurement.

### Mass scaling — does the pull follow the partner's mass?

Each star is accelerated by the *other* star's mass, so changing one star's mass
must change its partner's acceleration and leave its own alone. Four cells, one
knob, everything else held fixed:

| cell | phantom | ratio to matched | gap at `t = 200` | drift |
|---|---|---|---|---|
| `massratio_heavyphantom_d10_L64_N128_lev0` | `ω = 0.7603` vs a **lightened canonical** at `ω = 0.81` | `1.333` | **`+0.6028` — opens** | `+2.4306` |
| `runaway_pair_d10_L64_N128_lev0` | `ω = 0.7603` | `1.000` | `+0.0030` — rigid | `+2.8815` |
| `massscale_pair_d10_w0804_L64_N128_lev0` | `ω = 0.8040` | `0.7995` | `−0.5919` — closes | `+2.6732` |
| `massratio_w088_r060_d10_L64_N128_lev0` | `ω = 0.88` | `0.5974` | `−1.3816` — closes more | `+2.5362` |

**The gap changes sign exactly where the masses cross.** Monotone across a
factor of 2.2 in mass ratio, and a straight line through the four points puts
the zero crossing at ratio `1.06` against `1.00` predicted. That is the test no
artefact survives: an artefact can make a gap close, but it cannot make the
closing reverse depending on which of two stars is heavier.

Getting there required lightening the *canonical* rather than fattening the
phantom, because the phantom branch has a floor — see below.

On the quantitative side, the three *changed-partner* pull ratios — the
scaling law itself — land at `0.809 / 0.597 / 0.739` against
`0.7995 / 0.5974 / 0.7471` predicted (`1.2%`, `0.05%`, `1.1%`). Each cell's
*unchanged*-partner ratio doubles as an internal control, and it degrades as
the pair deforms: `0.973` (gap `−0.59`), `1.052` (gap `+0.60`), `0.725` (gap
`−1.38`). The separation-corrected fit is trustworthy while the pair stays
near-rigid and visibly is not once the gap moves by more than `~1`; quote the
changed-partner ratios, and **never the constants** — on the matched cell the
same fit returns `0.0161` and `0.0119` for two numbers both known to be
`0.014350`.

**Read the gap late, not early.** Every cell dips around `t = 50` and then
recovers, and the dip is larger than the signal until about `t = 100`:

| `t` | 25 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|
| phantom heavier | `−0.002` | `−0.212` | `−0.062` | `+0.201` | `+0.603` |
| equal | `−0.008` | `−0.255` | `−0.173` | `−0.110` | `+0.003` |

Judged at `t = 50` the reversed cell looks like every other cell. It crosses
zero near `t = 118`. Any statement about the gap belongs to `t ≳ 150`.

See `stars/star_family_massratio.csv` for the branch scan the frequencies were
chosen from — and for the bound it turned up: **both branches have a floor on
how light their stars get**, the canonical at `M+ ≈ 0.00776` near `ω = 0.92` and
the phantom at `|M−| ≈ 0.00791` near `ω = 0.94`, about `0.54` and `0.55` of the
matched mass. Past those the stars grow heavier again, and at `ω = 1.0` neither
branch has a bound star at all. A rung at `0.40` of matched cannot be built, on
either side.

### Mesh refinement — does the uniform-grid choice shape the result?

The campaign runs uniform grids so the convergence ladder means exactly one
thing, the cell size. `amrcheck_pair_d10_L64_N128_lev1` is the headline cell
with refinement switched on, run once to answer the referee question.

| | drift | gap | acceleration |
|---|---|---|---|
| `runaway_pair_d10_L64_N128_lev0` (uniform) | `+2.881465` | `+0.0030` | `1.463362e-04` |
| `amrcheck_pair_d10_L64_N128_lev1` (refinement on) | `+2.881489` | `+0.0031` | `1.463392e-04` |

Agreement to six digits — `0.00%` on both drift and acceleration. **Level 1 was
never created**: the log records 213 regrid *checks* at level 0 and not one step
on a finer level. The tagger triggers at `|χ − 1| = 0.02` and these runs peak
near `0.005`, so there is nothing for it to refine. Mesh refinement changes
nothing here, and the uniform-grid choice costs the result nothing.

### Stability survey — is the star family stable at all?

Supporting material, `t = 120`, lone canonical star at four frequencies.

| cell | lapse at birth | lapse at `t = 120` |
|---|---|---|
| `stability_canonical_w075_L64_N128_lev0` | `0.99051` | `0.99038` |
| `stability_canonical_w080_L64_N128_lev0` | `0.99209` | `0.99171` |
| `stability_canonical_w085_L64_N128_lev0` | `0.99343` | `0.99275` |
| `stability_canonical_w090_L64_N128_lev0` | `0.99464` | `0.99343` |

These four cells carry no `Ham_and_Mom_errors.txt` or `launch_config.sh`: they
reused the cached single-star initial data (no per-cell elliptic solve) and
were launched from the since-retired job queue, so `metadata.json` and
`evolution_params.txt` are their provenance record.

All four static. On the binary that subtracted the potential twice in the
stress-energy trace, these same initial-data bytes drove the lapse to
`0.13`–`0.44` and on to a horizon.

## Files in a cell

| file | what it is |
|---|---|
| `sector_dynamics.dat` | the main stream: per-sector core position, peak, momentum, gauge |
| `sector_barycenters.dat` | whole-sector centroid and rms radius — an independent tracker |
| `confinement.dat` | how much matter stays in the lumps |
| `collapse_diagnostics.dat` | lapse, `chi`, `K` — downsampled to `dt = 0.5` |
| `constraint_norms.dat`, `energy_conditions.dat`, `curvature_invariants.dat` | likewise downsampled |
| `psi4_*.dat` | extracted wave content |
| `well_tracking.dat` (PP N=256 only) | two-well positions vs time, derived from the frame slice cache — the merger record |
| `Ham_and_Mom_errors.txt` | elliptic solve convergence history |
| `evolution_params.txt`, `grtresna_params.txt` | what was evolved, what was solved |
| `launch_config.sh` | the exact environment the cell was launched with |

Every `.dat` carries its own `#` header naming the columns.

**One caveat on `sector_dynamics.dat`.** Column 10, `proper_sep`, was quantised
to whole cells until 2026-08-21 — it integrated between the two stars' integer
peak cells, so it could not resolve anything finer than `dx = 0.5`. It is fixed
now, but **files written before that date must not be compared with files
written after it**. Nothing in the results above uses it; `coord_sep` and the
continuous core positions are what every number here is measured from.

## What is not here

Two earlier generations of extracts were removed from this folder once they were
superseded. Both remain in git history and can be restored with
`git checkout d0f1a301 -- <path>`:

| removed | why it was wrong |
|---|---|
| `_corrupted_runs/` | ran on the binary that double-subtracted the potential in the stress-energy trace; canonical stars collapse |
| `_biased_initial_data/` | solve grid finer than the evolution grid, so the metric arrived displaced; includes `single_p_x32`, `single_m_x37` and `single_p_x37`, the "before" half of the control comparison |

The full history of both bugs, their diagnosis and their fixes is in
`research/bondi_dipole/docs/MatterDebugg.md`.
