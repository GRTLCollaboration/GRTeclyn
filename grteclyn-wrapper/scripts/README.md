# `grteclyn-wrapper/scripts`

Helper scripts for the GRTeclyn wrapper. Run from the repo root (`GRTeclyn/`).
Paths below are relative to that root unless noted.

## Layout

| Folder | What lives here |
|--------|-----------------|
| [`lib/`](lib/env.sh) | Shared `env.sh` (sets `GRTECLYN_ROOT`, `PYTHONPATH`, optional OpenMPI/GRTresna env). Sourced by other scripts; do not run directly. |
| [`plot/`](plot/) | **Plotfile streaming and post-run figures** — shared across examples (wormhole collapse, RadialRecipe episodes, Teo runs, etc.). |
| [`radial/`](radial/) | RadialRecipe smoke tests, batch runs, promotion, validation campaigns. |
| [`search/`](search/) | FTL / GRTresna CMA-ES search campaigns, tier-2 replay, offline scoring. |
| [`wormhole/`](wormhole/) | Wormhole-collapse helpers: `.gridinit` generators, rollback, archive to `SimResults/`. |
| [`build/`](build/) | Rebuild GRTresna MPI binary. |

There are **no** duplicate scripts at the `scripts/` root — only `README.md`, `lib/env.sh`, and the folders below.

---

## `plot/` — visualization (any example)

| Script | Use for |
|--------|---------|
| `plot_run.sh` | **Live** plotfile consumer during a run (Ψ₄, frames, optional delete). Default data dirs: `data_2gpu`, `data_supported`, `data`. |
| `plot_diagnostic.sh` | **Post-run** wormhole-style figures (constraints, collapse, Ψ₄ panels). |
| `plot_run_radial.sh` | Plotfile drain on a wrapper **episode** dir under `runs/`. |
| `plot_diagnostic_radial.sh` | Post-run RadialRecipe diagnostics (constraints, shell profiles). |
| `make_movies.sh` | Stitch PNG frame folders to MP4 (`ffmpeg` glob demuxer). |

```bash
# Supported / WormholeCollapse / RotatingWormholeCollapse output
./grteclyn-wrapper/scripts/plot/plot_run.sh /path/to/simulation/output
./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh /path/to/simulation/output

# RadialRecipe search episode
./grteclyn-wrapper/scripts/plot/plot_run_radial.sh runs/radialrecipe_gpu_smoke/<episode>
./grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode>
```

---

## `search/` — FTL campaigns

| Script | Purpose |
|--------|---------|
| `run_grtresna_search.sh` | **Production** matter-first search (GRTresna → `.gridinit` → GPU). |
| `run_ftl_search_cmaes.sh` | 9-D radial CMA-ES (geometry-first). |
| `run_ftl_search_nonspherical.sh` | 13-D gauge-angular search. |
| `run_ftl_search_directional.sh` | 21-D full-z directional search. |
| `run_optimize_loop.sh` | Generic env-driven CMA-ES loop. |
| `run_tier2_hq_188.sh` | HQ replay of non-spherical winner `eval_000188`. |
| `run_tier2_validation_long16.sh` | Legacy long validation (four spherical candidates). |
| `validate_tiers.py` | Offline falsification-tier assessment. |
| `rescore_grtresna_solved_ftl.py` | Re-score solved-geometry FTL on a campaign. |
| `summarize_scores.py` | Per-episode diagnostics table for a `runs/` tree. |

---

## `radial/` — RadialRecipe smoke & batch

| Script | Purpose |
|--------|---------|
| `run_radialrecipe_gpu_smoke.sh` | Single GPU episode (`SEED_NAME` / `CANDIDATE_ID` / `NONSPHERICAL_ID`). |
| `run_nonspherical_gpu_batch.sh` | Seven non-spherical candidates on GPUs 0–6. |
| `run_radialrecipe_gpu_promote.sh` | Resolution + long-time promotion ladder. |
| `run_subset.sh` | Ad-hoc parallel subset via `ENTRIES` + `RUNS_DIR`. |
| `validate_campaign.sh` | Paper seeds + four non-spherical candidates (GPUs 0–7). |

---

## `wormhole/` — collapse & Teo initial data

| Script | Purpose |
|--------|---------|
| `make_teo_wormhole_gridinit.py` | Teo rotating wormhole `.gridinit` CLI. |
| `make_rotating_wormhole_id.py` | Rotating wormhole ID helper. |
| `scan_rotating_wormhole_support.py` | Support-parameter scan. |
| `rollback` | Roll back run dir to a checkpoint/plot index. |
| `move_files.sh` | Archive outputs to `output/SimResults/`. |

---

## `build/`

| Script | Purpose |
|--------|---------|
| `rebuild_grtresna_mpi.sh` | Rebuild GRTresna `Main_ScalarFieldBH3d` MPI executable. |
