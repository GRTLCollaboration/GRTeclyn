# Q-ball trajectory search for evolving-geodesic shortcuts — git-friendly extract

(Campaign internal name: `qball_traj_fgeo_v1`.)

Light copy of the completed QD campaign
`runs/grtresna_qd/qball_traj_fgeo_v1` (400 evals, finished 2026-08-04),
which searched Q-ball lump trajectories for long-lived evolving-geodesic
FTL shortcuts. Full analysis and champion configurations:
**[CAMPAIGN_RESULTS.md](CAMPAIGN_RESULTS.md)**.

## Layout

- `CAMPAIGN_RESULTS.md` — results write-up (headline numbers, findings,
  exact champion parameters, caveats for the paper).
- `run/` — the campaign directory minus heavy binaries:
  - `trajectory.jsonl` — all 400 scoreboard rows (full genomes/overrides).
  - `archive.json`, `validation.json` — final MAP-Elites archive state.
  - `eval_000002 / 141 / 179 / 226 / 272 / 322 / 370` — the kept elite
    eval dirs: `params.txt`, `score.json`, diagnostics in `data/`,
    shortcut measurements in `small_data/` (`evolving_geodesic.json`,
    `ftl_timeseries.dat`, …), per-eval `run.log`.
    `eval_000322/lump_trajectories.png` is the champion trajectory figure.
- `logs/` — `grep '[qd]'` extracts of the two campaign logs (search
  progress, per-batch reports, record announcements).

## What was stripped (still on disk in `runs/`)

| item | size | why |
|---|---|---|
| `eval_*/initial_data.gridinit` | 531 MB × 7 | raw solver initial-data grids |
| `eval_*/small_data/metric_stack/` | 14 MB × 7 | metric snapshots (npz) used for ray tracing |
| full campaign logs (`runs/_logs/qball_traj_fgeo_v1.log`, `runs/grtresna_qd/qball_traj_fgeo_v1_resume2.log`) | 248 MB × 2 | replaced by the `[qd]` extracts in `logs/` |

Everything needed to re-run any configuration is retained: each eval's
`params.txt` (and the `overrides` in `trajectory.jsonl`) fully specifies
the candidate; initial data is regenerated deterministically by GRTresna.
Trajectory figures for any eval:
`grteclyn-wrapper/tests/visualisation/plot_lump_trajectories.py <eval_dir>`.
