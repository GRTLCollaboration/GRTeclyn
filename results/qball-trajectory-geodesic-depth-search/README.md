# Q-ball trajectory pure-depth hunt — git-friendly extract

(Campaign internal name: `qball_traj_fgeo_depth_v1`.)

Light copy of the completed QD campaign
`runs/grtresna_qd/qball_traj_fgeo_depth_v1` (200 evals, finished 2026-08-05),
which re-ran the `qball_traj_fgeo_v1` search with the objective switched from
saturated-depth-times-retention to **raw evolving-geodesic depth**
(`f_geo_depth`), warm-started from the seven v1 elites. Full analysis and
champion configurations: **[CAMPAIGN_RESULTS.md](CAMPAIGN_RESULTS.md)**.

Headline: deepest honest shortcut so far — **45.9 % path saving**, still
rising at the last emission time. The record genome is the v1 seed
(eval 370) unchanged — 200 evals never improved on it.

## Layout

- `CAMPAIGN_RESULTS.md` — results write-up (headline numbers, findings,
  champion parameters, caveats, next steps).
- `run/` — the campaign directory minus heavy binaries (no
  `initial_data.gridinit`, no `small_data/metric_stack`):
  - `trajectory.jsonl` — all 200 scoreboard rows (full genomes/overrides).
  - `archive.json`, `validation.json`, `pre_gpu_archive.json` — final
    MAP-Elites archive state and convergence report.
  - `ftl_champions.json`, `ftl_retention.jsonl` — per-metric champions.
  - `metadata.json` — campaign configuration.
  - `eval_000005 / 7 / 14 / 26 / 63 / 80` — kept champion eval dirs:
    `params.txt`, `score.json`, diagnostics in `data/`, shortcut
    measurements in `small_data/` (`evolving_geodesic.json`,
    `ftl_timeseries.dat`, …), per-eval `run.log`.
    `eval_000007` is the depth champion (= v1 eval 370 recipe).
- `logs/` — `[qd]` extract of the launch log (batch reports, prunes,
  record announcements).
