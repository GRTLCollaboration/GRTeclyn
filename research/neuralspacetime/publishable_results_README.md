# Matter-first automated discovery of transient spacetime shortcuts — PRD article audit pack

Paper: *Matter-first automated discovery of transient spacetime shortcuts*
([`research/neuralspacetime/article/research.tex`](../../research/neuralspacetime/article/research.tex)).

This directory is a **lightweight publishable mirror** of the run artifacts
and article tables needed to audit every plotted number and HQ claim in the
manuscript. It is **not** a full `/runs` dump: frames, AMReX plotfiles,
`initial_data.gridinit`, checkpoints, and `small_data/metric_stack/*.npz`
are intentionally omitted (hundreds of MB–GB each).

Authoritative provenance map (also copied here):
[`PLOT_DATA_SOURCES.txt`](PLOT_DATA_SOURCES.txt)
(same content as
[`research/neuralspacetime/article/PLOT_DATA_SOURCES.txt`](../../research/neuralspacetime/article/PLOT_DATA_SOURCES.txt)).

Regenerate this pack after new HQ / search runs:

```bash
bash research/neuralspacetime/pack_publishable_results.sh
```

Re-packing is **append-only** under `runs/`: existing run extracts are kept, so
a re-pack can never silently swap out the data the published figures were made
from (RC/RM were re-run on 2026-07-21, after the article figures were
generated). To deliberately re-pull them, run with `PACK_REFRESH_RUNS=1`.
Extracts whose source run has since been pruned from `/runs` (RF, DS, DL,
pump-free) survive here as the last remaining copy.

Paths below are relative to
`results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/`.

---

## Layout

```
results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/
├── README.md                          ← this file
├── PLOT_DATA_SOURCES.txt              ← figure/table → artifact map
├── article/data/                      ← PGFPlots tables used by research.tex
├── figures/                           ← raster figures used by research.tex
├── manifests/                         ← promote manifests (e.g. pump-free)
├── validation/                        ← matrix / launch certificates (JSON)
└── runs/                              ← light extracts mirroring /runs paths
    ├── grtresna_qd/
    │   ├── qball_traj_bicomplex_v1/   ← MAP-Elites (ME-87)
    │   └── qball_traj_canonical_v1/   ← canonical-only bound
    ├── grtresna_cmaes/
    │   └── qball_traj_bicomplex_cmaes_v1/  ← CMA-ES (CMA-146)
    ├── grtresna_promote/
    │   ├── sources/.../eval_000146/   ← frozen champion genome
    │   ├── bcma_rc_* … bcma_rf_*      ← resolution ladder
    │   ├── bcma_ds_* / bcma_dl_*      ← domain ladder
    │   ├── bcma_pfrm_*                ← pump-free twin
    │   └── bcma_*_freefall_corrected_*  ← f_ff companion runs (Table V)
    └── geometry_atlas/                ← stationary atlas (Sec. VI.D, Table VII)
        ├── geometry_atlas_topologies_*   ← topology-complete, 1000 evals
        ├── geometry_atlas_maxfgeo_*      ← pure-f_geo, 500 evals
        ├── geometry_atlas_lowexotic_*    ← E_- <= 0.1 ban, 500 evals
        └── geometry_atlas_breadth_*      ← breadth sweep
```

Each HQ run directory contains only:

| Path | Role |
|---|---|
| `score.json` | Objective summary |
| `params.txt` | Evolution / pump knobs |
| `small_data/evolving_geodesic.json` | \(f_\mathrm{geo}\) emit sweeps / peaks |
| `small_data/confinement.dat` | Confinement timeseries |
| `small_data/freefall_observer_timing.json` | \(f_\mathrm{ff}\) (freefall-corrected runs only) |
| `data/constraint_norms.dat` | Hamiltonian / momentum norms |
| `data/collapse_diagnostics.dat` | Collapse / pump-work inputs |

Search eval directories contain only `score.json`, `metadata.json`,
`params.txt`, and `initial_data.matter.json` (no gridinit).

Each atlas campaign contains `archive.json`, `metadata.json`, `state.json`,
`summary.json`, `trajectory.jsonl`, and `elites/*.json`. The per-eval
`trajectory.jsonl` rows carry `f_geo`, the `(f_shift, E_-)` descriptors and the
reject flag, which is everything needed to recompute the quoted atlas
statistics; the 74 MB-per-eval `*.gridinit` grids, the `evals/*.json` dumps and
the elite galleries are intentionally omitted.

---

## Article figure / table → pack path

| Paper item | Pack path |
|---|---|
| `fig:qd_progress` | `article/data/campaign_progress.txt` ← `runs/grtresna_qd/qball_traj_bicomplex_v1/trajectory.jsonl` |
| `fig:cmaes_progress` | `article/data/cmaes_progress.txt` ← `runs/grtresna_cmaes/.../trajectory.jsonl` |
| `fig:fgeo_validation` | `article/data/fgeo_sweeps.txt` ← RC/RM/RF `small_data/evolving_geodesic.json` |
| `fig:constraints` | `article/data/constraints_rm.txt` ← RM `data/constraint_norms.dat` |
| `fig:constraint_convergence` | `article/data/constraints_resolution.txt`, `constraint_order.txt` |
| `fig:lump_dynamics` | `figures/fig_eval146_lump_trajectories.png` + frozen `initial_data.matter.json` |
| `fig:rf_field_maps` | `figures/fig_rf_{rho,shift}_t{0,4,16,30}.png` (frames not shipped) |
| `fig:null_constraint` | `article/data/null_constraint_ray*.txt` |
| `tab:endpoint_gauge` | `article/data/endpoint_gauge_rm.txt` |
| `tab:pump_compare` | RM vs `runs/.../bcma_pfrm_*` + `article/data/pump_work_budget.txt` |
| `tab:matrix` / domain ladder | RC/RM/RF + DS/DL HQ dirs under `runs/grtresna_promote/` |
| Canonical bound | `runs/grtresna_qd/qball_traj_canonical_v1/` |
| Casimir bounds | `article/data/casimir_bounds.txt` |

Promote-run abbreviations (candidate 146):

| Tag | Directory |
|---|---|
| RC | `runs/grtresna_promote/bcma_rc_L128_N192_t30_hq_eval000146/` |
| RM | `runs/grtresna_promote/bcma_rm_L128_N256_t30_hq_eval000146/` (headline) |
| RF | `runs/grtresna_promote/bcma_rf_L128_N384_t30_hq_eval000146/` |
| DS | `runs/grtresna_promote/bcma_ds_L96_N192_t30_hq_eval000146/` |
| DL | `runs/grtresna_promote/bcma_dl_L160_N320_t30_hq_eval000146/` |
| PF | `runs/grtresna_promote/bcma_pfrm_L128_N256_t30_hq_eval000146/` |

---

## What is deliberately missing

- `/runs/**/frames/**` (movie PNGs; selected RF stills live under `figures/`)
- `/runs/**/plt*`, plotfiles, checkpoints
- `initial_data.gridinit` and other binary ID dumps
- `small_data/metric_stack/*.npz` (needed only to *regenerate*
  endpoint-gauge / null-constraint tables via
  `research/neuralspacetime/analysis/endpoint_gauge_and_null_constraint.py`)

Those stay on the cluster under gitignored `/runs`. The reduced tables under
`article/data/` are what `research.tex` actually `\input`s.

---

## Related source tree (not duplicated here)

| Path | Contents |
|---|---|
| `research/neuralspacetime/article/research.tex` | Manuscript |
| `research/neuralspacetime/analysis/` | Table regenerators |
| `research/neuralspacetime/validation/` | Same certificates as `validation/` here |
| `/runs` (gitignored) | Full machine-local campaign outputs |
