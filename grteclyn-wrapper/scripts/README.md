# `grteclyn-wrapper/scripts`

Helper scripts for the GRTeclyn wrapper. Run from the repo root (`GRTeclyn/`).
Paths below are relative to that root unless noted.

## Layout

| Folder | What lives here |
|--------|-----------------|
| [`campaigns/`](campaigns/README.md) | **Three-stage FTL pipeline** — QD, CMA-ES, HQ promotion launchers + shared config. Start here for production searches. |
| [`lib/`](lib/env.sh) | Shared `env.sh` (sets `GRTECLYN_ROOT`, `PYTHONPATH`, optional OpenMPI/GRTresna env). Thin redirects to `campaigns/lib/` for search/promote defaults. |
| [`plot/`](plot/) | **Plotfile streaming and post-run figures** — shared across examples (wormhole collapse, RadialRecipe episodes, Teo runs, etc.). |
| [`radial/`](radial/) | RadialRecipe smoke tests, batch runs, promotion, validation campaigns. |
| [`search/`](search/) | FTL offline analysis (tier validation, campaign reports, artifact cleanup). |
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
| `make_movies.sh` | Stitch PNG frame folders to MP4; writes `<run>/movies/movie_<field>_<axis>.mp4`. |

```bash
# Supported / WormholeCollapse / RotatingWormholeCollapse output
./grteclyn-wrapper/scripts/plot/plot_run.sh /path/to/simulation/output
./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh /path/to/simulation/output

# RadialRecipe search episode
./grteclyn-wrapper/scripts/plot/plot_run_radial.sh runs/radialrecipe_gpu_smoke/<episode>
./grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode>
```

---

## `campaigns/` — three-stage FTL pipeline

See **[`campaigns/README.md`](campaigns/README.md)** for the full guide.

| Stage | Script | Purpose |
|-------|--------|---------|
| 0 — QD | `campaigns/qd/run.sh` | MAP-Elites survey (N=128, t=16, 4D search profile) |
| 1 — CMA-ES | `campaigns/cmaes/run.sh` | Local refine from QD trajectory (same grid as stage 0) |
| 2 — HQ | `campaigns/hq/run_batch.sh` | Full-res promotion (N=256, t=30, frames, 4D HQ verify) |
| HQ single | `campaigns/hq/replay_eval.py` | One eval replay |

Plotfile **consume + on-the-fly delete** (`--consumer-delete`, keep last 3) is default for all stages. PNG **frame rendering** is off for QD/CMA-ES (`GRTECLYN_FRAMES=0`) and on for HQ (`GRTECLYN_FRAMES=1` in `promote_common.sh`).

---

## `search/` — FTL analysis utilities

Launchers live under `campaigns/` only. This folder keeps offline helpers:

| Script | Purpose |
|--------|---------|
| `project_geometry_motif.py` | Scout → GRTresna projection → post-load gate → optional evolve. |
| `run_matter_geometry_smokes.sh` | End-to-end matter–geometry consistency smokes. |
| `validate_tiers.py` | Offline falsification-tier assessment. |
| `report_campaign_ftl.py` | Rank evals by per-frame FTL peaks (`f_geo`, speed, lifetime) from `ftl_timeseries.dat`. |
| `rescore_grtresna_solved_ftl.py` | Re-score solved-geometry FTL on a campaign. |
| `cleanup_heavy_sim_artifacts.sh` | Delete heavy plotfile/checkpoint trees under a QD campaign. |

HQ promotion defaults: [`campaigns/lib/promote_common.sh`](campaigns/lib/promote_common.sh).

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
