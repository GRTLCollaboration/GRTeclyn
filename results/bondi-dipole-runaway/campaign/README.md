# Campaign extracts — the ω = 0.75 runaway series

Packed from `runs/bondi/runaway/` by `research/bondi_dipole/pack_runaway.sh`.
One folder per cell.  Everything here is small and text-only: PNG frames,
movies, plotfiles and the gridinit stay in the gitignored run tree, because
they are large and carry no number any analysis reads.

This is the first series produced *after* two initial-data/matter bugs were
found and fixed — the potential double-counted in the trace of the matter
stress tensor, and the slicing asymmetry that gave canonical stars a collapsing
birth slice.  Every extract produced before those fixes has been moved to
[`_corrupted_runs/`](_corrupted_runs/) and is kept as the evidence for them.
See `research/bondi_dipole/docs/MatterDebugg.md` for the full history.

## Cells

All cells share one grid — `N = 128`, `L = 64`, `dx = 0.5`, no refinement —
one matter rung (`m = 1`, `λ = 10240`, `μ = 21845333`) and one canonical star
frequency, `ω = 0.75`.  `d` is the initial coordinate separation.

| cell | pair | `d` | `t` | what it is |
|---|---|---|---|---|
| `pm_eqm`       | canonical + phantom | 10 | 200 | headline cell: the two stars carry equal `\|ADM\|` mass (phantom at `ω = 0.7603`) |
| `pm_eqm_sep8`  | canonical + phantom |  8 | 200 | separation scan |
| `pm_eqm_sep12` | canonical + phantom | 12 | 200 | separation scan |
| `pm_eqm_sep16` | canonical + phantom | 16 | 200 | separation scan, widest — the cleanest point-mass comparison |
| `pm`           | canonical + phantom | 10 | 200 | both stars at `ω = 0.75`, so the phantom is 5% heavier |
| `pm_sep8`      | canonical + phantom |  8 | 200 | mismatched separation scan |
| `pm_sep12`     | canonical + phantom | 12 | 200 | mismatched separation scan |
| `pp`           | two canonical       | 10 | 200 | control: no runaway is possible |
| `mm`           | two phantom         | 10 | 200 | control: no runaway is possible |
| `pm_eqm_t400`  | canonical + phantom | 10 | 400 | long baseline; **still evolving when this was packed** — re-run the packer to refresh |

## What each folder holds

| file | contents |
|---|---|
| `sector_dynamics.dat` | per-sector core position, peak, separation, momentum, shift |
| `sector_barycenters.dat` | independent whole-sector barycentre and rms radius |
| `confinement.dat` | total and peak activity, rms radius, confined fraction — read this before believing any geometry number |
| `collapse_diagnostics.dat` | lapse, `chi`, `max abs K`, horizon search |
| `constraint_norms.dat` | Hamiltonian and momentum violation, whole grid and interior |
| `energy_conditions.dat`, `curvature_invariants.dat` | NEC/WEC/SEC/DEC minima; Ricci and `K_ij` invariants |
| `psi4_*.dat`, `boundary_flux.dat`, `areal_radius.dat`, `shell_profiles.dat`, `ftl_timeseries.dat` | wave extraction, outgoing flux, radial profiles |
| `evolution_params.txt`, `grtresna_params.txt` | the full parameter sets actually used |
| `Ham_and_Mom_errors.txt` | elliptic solve residual history |
| `metadata.json`, `initial_data.matter.json` | provenance of the initial data |
| `launch_config.sh` | the exact environment that produced the cell |

The four evolution streams (`collapse_diagnostics`, `constraint_norms`,
`energy_conditions`, `curvature_invariants`) are downsampled to `dt = 0.5` from
the every-step originals and carry an added header line naming their columns;
the columns themselves are untouched.  Everything in `small_data` is copied
verbatim and already self-documenting.

## What the series measures

`ratio` is the measured midpoint displacement over the Newtonian prediction
`0.5 * (G M / d^2) * t^2`; `null` is the transverse displacement, forbidden by
symmetry, as a fraction of the forward one — it measures the tracker's own
systematic error.  `gap` is how much the coordinate separation closed.

| cell | `d` | ratio | null | gap | peak drift | interior `L2 Ham` |
|---|---|---|---|---|---|---|
| `pm_eqm_sep8`  |  8 | 1.089 | -2.0% | 0.722 | +0.46% | 9.7e-06 |
| `pm_eqm`       | 10 | 1.042 | -1.4% | 0.724 | +0.01% | 1.1e-05 |
| `pm_eqm_sep12` | 12 | 1.012 | -0.8% | 0.698 | -0.02% | 1.2e-05 |
| `pm_eqm_sep16` | 16 | 0.990 | +0.2% | 0.701 | -0.02% | 1.3e-05 |
| `pm_sep8`      |  8 | 1.078 | -1.9% | 0.455 | +0.19% | 9.9e-06 |
| `pm`           | 10 | 1.039 | -1.2% | 0.563 | +0.04% | 1.1e-05 |
| `pm_sep12`     | 12 | 1.012 | -0.5% | 0.590 | -0.00% | 1.2e-05 |
| `pp`           | 10 | -0.009 | — | 0 | +0.25% | 1.7e-05 |
| `mm`           | 10 | 0.036 | — | 0 | +0.66% | 1.3e-05 |

The mixed pairs track the prediction across a factor of four in the predicted
force and agree to 1% at the widest separation, where treating each star as a
point mass is most defensible.  The two same-sign controls, built by the same
code on the same grid, move about a hundred times less relative to the same
prediction.  No cell lost its stars: peak amplitude holds to better than 0.7%
everywhere.

Two caveats travel with the numbers.  The coordinate gap closes at a rate that
does not depend on separation, so it is a fixed drift set at birth rather than
gravity — see the debug document for the argument.  And **every cell here is a
single resolution**; the convergence test is not yet run.

## Refreshing

```
bash research/bondi_dipole/pack_runaway.sh
```

Safe at any time, including while cells are still evolving: each cell folder is
rebuilt from scratch, and cells with no time series yet are skipped.

## The older extracts

[`_corrupted_runs/`](_corrupted_runs/) holds everything packed before the two
fixes, with its own README.  Two things to know before reaching for the old
tooling: `research/bondi_dipole/pack_results.sh` and `pack_campaign.sh` read a
`runs/bondi/rerun` tree that no longer exists, so they now serve only as the
provenance record for those extracts; and the scripts in `../analysis/` address
cells by name, so they resolve into `_corrupted_runs/`, not into this campaign.
