# `grteclyn-wrapper/scripts`

Helper scripts for the closed-loop GRTeclyn metric-discovery pipeline. Run them
from the repo root (`GRTeclyn/`). All paths below are relative to that root.

Most scripts pick up GRTeclyn via `env.sh` (sourced automatically) or hard-code
the repo path; the FTL-search campaign scripts call the Python CLI directly
(`uv run python -m grteclyn_wrapper ... optimize`).

---

## Environment

| Script | Purpose |
|--------|---------|
| `env.sh` | Sets `GRTECLYN_ROOT`, `PYTHONPATH`, and optional local OpenMPI on `PATH`, then `cd`s to the repo root. Sourced by the smoke / promote / batch / optimize scripts; not run directly. |

## Closed-loop FTL search campaigns (CMA-ES)

These are the **current** search drivers. Each launches the custom CMA-ES loop
across all 8 H100s, streams plotfiles through the consumer (frames extracted on
the fly, heavy plotfiles deleted, `--consumer-keep-last 5` retained so the
*evolved* operational FTL `F_op^ev` and effective energy conditions guide the
search), and writes to `runs/<name>/{trajectory.jsonl, result.json, eval_*/}`.

| Script | Search space | Purpose |
|--------|--------------|---------|
| `run_ftl_search_cmaes.sh` | 9-D radial (`chi`/`alpha`/`beta` basis) | Baseline spherically-symmetric FTL search, normal matter (`--no-phantom`); least-exotic shift-driven warp channels. |
| `run_ftl_search_nonspherical.sh` | 13-D (`--nonspherical`) | Adds ℓ=1/ℓ=2 Legendre angular modes on the **gauge** fields (lapse + shift) → directional channels at zero extra exotic-matter cost. |
| `run_ftl_search_directional.sh` | 21-D (`--nonspherical`, full-z) | As above plus searchable radial center/width of every angular mode, and a **full z-domain** (Sommerfeld `lo_boundary`, centered) so true odd-parity (fore/aft dipole) channels are allowed. Roughly 2× the cost. |

Usage: `bash scripts/run_ftl_search_nonspherical.sh [NAME] [MAXGEN] [SEED]`.

## Tier-2 validation (single best candidate, high quality)

| Script | Purpose |
|--------|---------|
| `run_tier2_hq_188.sh` | **Current** high-quality validation of the non-spherical winner `eval_000188`. High base resolution (`N=192`), modest refinement, `t=16`; uses the streaming consumer so frames are extracted on the fly and plotfiles are deleted (`--consumer-keep-last 2`). Output stays ~1 GB. |
| `run_tier2_validation_long16.sh` | *Legacy.* Re-validates four specific warp candidates from the earlier spherical-phantom campaign `ftl_search_cmaes_01` (hard-coded coefficients) at `N=96`, `t=16`, retaining all plotfiles. Kept only as a reproducibility record of that campaign; safe to delete if those runs are no longer of interest. |

## Smoke / batch / promotion infrastructure

The building blocks behind the validation campaigns and the paper's smoke tables.
These are env-var driven (see the top-level `README.md`).

| Script | Purpose |
|--------|---------|
| `run_radialrecipe_gpu_smoke.sh` | Single-episode driver: build (optional), generate `params.txt`, `check_params`, short GPU evolution, parse diagnostics + score. Selects initial data via `SEED_NAME` / `CANDIDATE_ID` / `NONSPHERICAL_ID`. The core unit the batch scripts call. |
| `run_nonspherical_gpu_batch.sh` | Runs the seven `accepted_ray_sane` non-spherical candidates in parallel across GPUs via the smoke driver. |
| `run_radialrecipe_gpu_promote.sh` | Resolution + long-time promotion ladder (`N=64/96/128`, `t=5`) for smoke survivors. |
| `run_optimize_loop.sh` | Generic env-var-driven CMA-ES loop (one candidate per GPU per generation). The `run_ftl_search_*.sh` scripts are the campaign-specific successors that call the Python `optimize` CLI directly; this remains as a generic/standalone entry point. |
| `run_subset.sh` | Utility: run an arbitrary `ENTRIES` subset (`label:seed|nonsph:id`) in parallel via the smoke driver. Used for ad-hoc multi-candidate runs. |
| `validate_campaign.sh` | One-off campaign: the four paper seeds (flat / Ellis–Bronnikov / Alcubierre / Schwarzschild) + four non-spherical candidates, in parallel across GPUs 0–7. |

## Post-processing & utilities

| Script | Purpose |
|--------|---------|
| `make_movies.sh` | Stitch pre-rendered PNG frames into mp4 movies (`ffmpeg` glob demuxer, robust to gapped/step-numbered frame names). `make_movies.sh EPISODE_DIR [...] [--framerate N] [--only chi_z K_z]`. |
| `summarize_scores.py` | Re-score every episode under a runs dir and print a per-candidate diagnostics table (score, gate, survival, FTL, energy conditions, curvature). `python scripts/summarize_scores.py RUNS_DIR [target_stop_time]`. |
| `plot_run_radial.sh` | Plotfile drain / diagnostic plots on a wrapper episode dir (streams plotfiles, optional frames). |
| `plot_diagnostic_radial.sh` | Post-run RadialRecipe diagnostics (constraints, collapse, shell profiles). |
| `plot_run.sh` | Wormhole-collapse sidecar: stream plotfiles, extract $\Psi_4$, render frames. |
| `plot_diagnostic.sh` | Post-run wormhole figures (constraints, collapse, $\Psi_4$ panels, evolution strips). |
| `move_files.sh` | Archive run outputs + figures into `grteclyn-wrapper/output/SimResults/`. |
| `rollback` | Roll back a run directory to a checkpoint/plot index. |

---

### Notes on obsolete scripts

`run_tier2_nonspherical_188.sh` was **removed**: it validated the same winner
`eval_000188` but retained all (multi-GB) plotfiles at `N=96`. It is fully
superseded by `run_tier2_hq_188.sh`, which runs at higher resolution and uses the
streaming frame-extraction + plotfile-deletion workflow.
