# Campaign extracts — the Bondi dipole runaway

Publication-level data only. One folder per cell, text-only: PNG frames,
movies, plotfiles and the 550 MB initial-data file stay in the gitignored run
tree at `runs/bondi/staging/archive/`, because they carry no number an analysis
reads.

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

### Stability survey — is the star family stable at all?

Supporting material, `t = 120`, lone canonical star at four frequencies.

| cell | lapse at birth | lapse at `t = 120` |
|---|---|---|
| `stability_canonical_w075_L64_N128_lev0` | `0.99051` | `0.99038` |
| `stability_canonical_w080_L64_N128_lev0` | `0.99209` | `0.99171` |
| `stability_canonical_w085_L64_N128_lev0` | `0.99343` | `0.99275` |
| `stability_canonical_w090_L64_N128_lev0` | `0.99464` | `0.99343` |

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
