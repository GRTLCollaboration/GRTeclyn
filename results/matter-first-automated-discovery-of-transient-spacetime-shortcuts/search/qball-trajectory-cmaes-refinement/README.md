# Q-ball trajectory CMA-ES depth refinement — git-friendly extract

(Campaign internal name: `qball_traj_fgeo_depth_cmaes_v1`.)

Light copy of the completed CMA-ES campaign
`runs/grtresna_cmaes/qball_traj_fgeo_depth_cmaes_v1` (200 evaluations, 13
generations, finished 2026-08-06 21:25 UTC). It is the local-refinement
follow-up to the MAP-Elites campaign in
[`../qball-trajectory-geodesic-depth-search/`](../qball-trajectory-geodesic-depth-search/),
using the identical objective (`f_geo_depth`), physics, gates and search space —
only the search strategy differs. Full analysis:
**[CAMPAIGN_RESULTS.md](CAMPAIGN_RESULTS.md)**.

**Headline: 48.38 % evolving-geodesic path saving** (eval 195), up from 45.90 %
for the MAP-Elites champion and 46.28 % for the warm-start seed it began from.
The best result of the entire campaign came from the last, half-sized
generation — the search was still improving when the eval budget ran out.

Two caveats travel with that number:

- About **+1.7 of the +2.1 points is clean**; the remaining +0.4 coincides with
  a jump to a different initial-data solution branch that requires 40 % more
  exotic matter and sits at 95 % of the constraint acceptance gate. The
  conservative refined figure is **47.97 %** (eval 185). Both genomes are kept
  in `genomes/`.
- The depth is still **clipped by the measurement window** — the emission sweep
  rises monotonically to its final launch at t = 18 with no turnover, so 48.38 %
  is a lower bound.

## Layout

- `CAMPAIGN_RESULTS.md` — write-up: headline numbers, per-generation table,
  findings, winner configuration, caveats, run commands, next steps.
- `analysis/`
  - `per_eval.csv` — all 200 rows: generation, status, depth, emission and
    arrival times, constraint residuals, exotic-matter integral, branch label.
  - `per_generation.csv` — 13 rows: best and median depth, per-branch bests and
    counts, crash and gate-rejection counts.
- `genomes/` — the three genomes that matter, each with its measurements:
  - `winner_eval195_familyB.json` — 48.38 %.
  - `best_clean_branch_eval185_familyA.json` — 47.97 %, conservative branch.
  - `warm_start_seed.json` — 46.28 %, the starting point.
- `logs/cmaes_progress_evals_1-200.log` — per-eval and per-generation progress
  lines extracted from the launch log.
- `run/` — the campaign directory minus heavy binaries (no
  `initial_data.gridinit`, no `small_data/metric_stack`, no plotfiles):
  - `trajectory.jsonl` — all 200 scoreboard rows with full genomes.
  - `metadata.json` — campaign configuration and search-space bounds.
  - `warm_start_gen1_seed.jsonl` — the trajectory this run was warm-started from.
  - `ftl_champions.json`, `ftl_retention.jsonl` — per-metric champions.
  - `eval_000193` … `eval_000198` — the six deepest evals: `params.txt`,
    `score.json`, `small_data/` (emission sweep in `evolving_geodesic.json`,
    `ftl_timeseries.dat`, …), `grtresna/Ham_and_Mom_errors.txt` (solver residual
    history), `data/` diagnostics, per-eval `run.log`.

Pack size: 27 MB.

## Reproducing

Launcher, exact environment, monitoring, resume and stop commands are in
[CAMPAIGN_RESULTS.md §7](CAMPAIGN_RESULTS.md#7-run-commands--properties).
Replaying the winner from this pack:

```bash
python grteclyn-wrapper/scripts/campaigns/hq/replay_eval.py \
  --trajectory results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/trajectory.jsonl \
  --eval 195 --objective-mode f_geo_depth
```
