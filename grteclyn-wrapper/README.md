# grteclyn-wrapper

Run isolated GRTeclyn episodes from the repo root (`GRTeclyn/`). For RadialRecipe GPU smoke tests, the wrapper can stream plotfiles into small data during the run and delete heavy HDF5 dirs afterward.

## Prerequisites

From the GRTeclyn repo root:

```bash
cd /path/to/GRTeclyn
uv sync   # Python deps including yt for plotfile extraction
```

First build (single GPU, no MPI):

```bash
BUILD=1 bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```

## Single GPU run (one guessed shape)

Pick **one** initial-data source:

| Env var | Example |
|---------|---------|
| `SEED_NAME` | `ellis_bronnikov` |
| `CANDIDATE_ID` | `bubble_wall_016`, `random_000` |
| `NONSPHERICAL_ID` | `quadrupole_bubble_001`, `dipole_lopsided_000` |

```bash
# Known seed
BUILD=0 SEED_NAME=ellis_bronnikov CUDA_VISIBLE_DEVICES_OVERRIDE=0 \
  bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh

# Spherical guesser candidate
BUILD=0 CANDIDATE_ID=bubble_wall_016 CUDA_VISIBLE_DEVICES_OVERRIDE=1 \
  bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh

# Non-spherical guessed shape
BUILD=0 NONSPHERICAL_ID=quadrupole_bubble_001 CUDA_VISIBLE_DEVICES_OVERRIDE=2 \
  bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```

Outputs go to `runs/radialrecipe_gpu_smoke/<name>_gpu_t<stop_time>_<stamp>/`.

## Scripts index

All helper scripts live in [`scripts/`](scripts/README.md), which documents each
one's purpose. The most important groups:

- **Closed-loop FTL search (current):** `run_ftl_search_cmaes.sh` (9-D radial),
  `run_ftl_search_nonspherical.sh` (13-D, gauge angular modes),
  `run_ftl_search_directional.sh` (21-D, full-z). These launch the CMA-ES loop
  across all GPUs and stream plotfiles (frames on the fly, heavy data deleted).
- **Tier-2 validation:** `run_tier2_hq_188.sh` (high-res streaming validation of
  the non-spherical winner).
- **Smoke / batch infrastructure:** `run_radialrecipe_gpu_smoke.sh` (single
  episode), `run_nonspherical_gpu_batch.sh`, `run_radialrecipe_gpu_promote.sh`,
  `run_optimize_loop.sh`, `run_subset.sh`, `validate_campaign.sh`.
- **Post-processing:** `make_movies.sh`, `summarize_scores.py`,
  `plot_run_radial.sh`.

```bash
# Example: launch the non-spherical FTL campaign
bash grteclyn-wrapper/scripts/run_ftl_search_nonspherical.sh
# Validate the winner at high quality (streaming frames, plotfiles deleted)
bash grteclyn-wrapper/scripts/run_tier2_hq_188.sh 0 val16hq_nonsph_eval188
```

## In-situ diagnostics & matter sector

Each RadialRecipe run now emits three diagnostic tables under `data/` (read back
automatically by `grteclyn_wrapper.metrics.read_episode_metrics`):

| File | Columns | Probes |
|------|---------|--------|
| `constraint_norms.dat` | `L2_Ham L2_Mom min_rho_req max_rho_req integral_neg_rho` | constraint satisfaction; `min_rho_req < 0` flags geometries needing exotic matter |
| `energy_conditions.dat` | `matter_min_{NEC,WEC,SEC,DEC} matter_integral_NEC_violation` | observer-sampled energy conditions of the **evolved matter** |
| `curvature_invariants.dat` | `max_abs_ricci_scalar max_ricci_tensor_sq max_Kij_sq L2_ricci_scalar` | coordinate-invariant geometry |

A general, mechanism-agnostic operational-FTL measure (`ftl_general.py`,
Dijkstra shortest-coordinate-time vs. flat baseline) is also computed per episode
and exposed as `EpisodeMetrics.general_ftl`. It is not warp-specific: any
geometry whose coordinate light cones open a faster channel scores `f_op > 0`.

### Spacetime → matter: exotic matter is now evolved when needed

A wormhole/warp geometry generally requires exotic (phantom, `rho <= 0`) matter
to satisfy the Hamiltonian constraint. The constrained recipe (`--phantom`)
already solves for the scalar profile under phantom coupling. The C++ level now
evolves the **matching** matter: when `--phantom` is set the wrapper injects
`recipe_exotic_matter = 1`, and `RadialRecipeLevel` evolves an `ExoticScalarField`
(`T_munu = -recipe_support_strength * canonical`) in the RHS, the constraints,
and the energy-condition diagnostic. With a canonical seed it falls back to the
ordinary `ScalarField`. Verified on Ellis–Bronnikov: the evolved matter
NEC/WEC/SEC/DEC go negative (NEC `-0.07`, integrated violation `~2.1`) exactly
where the geometry demands it, instead of the `~0` null result a canonical field
gives.

### Two findings worth keeping in mind

1. **Matter-sector EC is a null result *only* with a canonical field.** A canonical
   `ScalarField` has `rho >= 0`, so its NEC/WEC are `~0` by construction; `--phantom`
   used to shape only the initial data. With `recipe_exotic_matter = 1` the evolved
   matter is genuinely exotic and the `matter_*` columns reveal the violation. The
   curvature invariants and the general FTL measure read the geometry directly and
   are meaningful regardless. The geometry-sourced effective stress energy
   (`T^eff = G / 8pi`) is what the `matter_*` columns *cannot* see — that is
   evaluated post-hoc from plotfiles by `warpfactory.py`.
2. **Evolved-data FTL / effective-EC needs more plot vars.** Plotfiles currently
   store `chi`, `K`, `lapse` — not the full `h_ij`/shift — so the general FTL runs
   on the `t=0` reconstructed slice. To run it (and the effective EC) on evolved
   spacetimes, add the metric components to `amr.plot_vars` and feed the plotfile
   grid into `operational_ftl_on_grid` / `warpfactory.py`.

## Batch: 7 non-spherical shapes on GPUs 0–6

```bash
BUILD=0 bash grteclyn-wrapper/scripts/run_nonspherical_gpu_batch.sh
```

Outputs: `runs/radialrecipe_nonspherical/`. One log per candidate: `<id>_<stamp>.log`.

## Plotfile consumer (default on)

With `CONSUME_PLOTFILES=1` (default), each run:

1. Starts a sidecar `consume_plotfiles` process while the GPU simulation runs.
2. Appends `small_data/shell_profiles.dat`, `small_data/areal_radius.dat`, optional PNG frames.
3. Deletes processed plotfile dirs (`--keep-last 1`).
4. Runs a **post-sim drain** for any backlog.

Useful env vars:

| Variable | Default | Meaning |
|----------|---------|---------|
| `CONSUME_PLOTFILES` | `1` | Enable streaming extraction |
| `CONSUMER_DELETE` | `1` | Delete HDF5 plot dirs after extract |
| `CONSUMER_RADII` | `4 8` | Extraction radii |
| `PLOT_INTERVAL` | `10` if consumer on, else `1` | Plotfile write cadence |
| `STOP_TIME` | `2.0` | Simulation stop time |
| `N_FULL` | `64` | Grid resolution |
| `CUDA_VISIBLE_DEVICES_OVERRIDE` | `0` | GPU index for single run |

Disable consumer (keep all plotfiles):

```bash
CONSUME_PLOTFILES=0 bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```

## Post-run plots (no GW / Psi4)

```bash
bash src/scripts/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode_dir>
```

Writes EPS/PNG to `src/visualisation/plots/radial/` (constraints, collapse, shell profiles).

## Manual plotfile drain (if needed)

If a run finished before extraction completed:

```bash
bash src/scripts/plot_run_radial.sh runs/radialrecipe_nonspherical/<episode_dir> --no-delete
# or one-shot batch consume (no watch):
uv run python -m src.visualisation.process_wave.consume_plotfiles \
  --data runs/radialrecipe_nonspherical/<episode_dir> \
  --out runs/radialrecipe_nonspherical/<episode_dir>/small_data \
  --radii 4 8 --no-psi4 --shell-fields chi lapse K --areal-radius \
  --delete --keep-last 1 -j 4
```
