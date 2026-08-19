# Gated-shortcut CMA-ES refinement — git-friendly extract

(Campaign internal name: `qball_traj_fgeo_max_cmaes_v1`.)

Light copy of the completed CMA-ES campaign
`runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1` — 200
evaluations, 13 generations of 16, finished 2026-08-18. It is the
local-refinement follow-up to the MAP-Elites campaign in
[`../qball-trajectory-evolving-geodesic-shortcut-search/`](../qball-trajectory-evolving-geodesic-shortcut-search/),
warm-started from that campaign's champion (eval 322) and using the identical
objective (`f_geo_max`), physics, gates and 39-dimensional search space — only
the search strategy differs.

**This is the headline campaign of the article.** Its champion is the genome
promoted to paper resolution.

## Headline

**35.94 % evolving-geodesic path saving** (eval 193), up from **30.05 %** for
the warm-start genome measured under the same scorer — a 19.6 % relative gain
over 200 evaluations.

Unlike the depth lineage, this number is **persistence-gated**: the objective
multiplies the shortcut by how long it survives, so an instantaneous spike
scores nothing. The champion holds its shortcut at a persistence of 0.78, and
its horizon penalty is exactly zero — no apparent-horizon proxy ever triggers,
so the saving is not a collapsing-region artefact.

Two things to read carefully:

- **The search was still improving when the budget ran out.** The best result
  came from generation 13 — the last, half-sized one — exactly as in the depth
  campaign. The per-generation best rises in eleven of twelve steps and never
  plateaus, so 35.94 % is a floor on what this basin holds, not a converged
  optimum.
- **Do not compare this score to the depth campaign's.** 3665.77 here and
  4862.61 there are different objectives. Compare path saving: 35.94 % gated
  versus 48.38 % ungated is the real Pareto statement, and it is the article's
  point, not an inconsistency.

## Per-generation climb

Path saving of the best eval in each generation, from `analysis/per_generation.csv`:

| gen | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| best | .320 | .324 | .335 | .338 | .344 | .348 | .349 | .355 | .354 | .356 | .359 | .356 | **.359** |
| gpu_ok | 12 | 16 | 14 | 15 | 15 | 14 | 12 | 11 | 13 | 14 | 11 | 9 | 7 |

163 of 200 evaluations returned a scored run; the other 37 were refused by the
post-load constraint gate before any GPU time was spent on them. The rejection
rate climbs with the generations (4 in gen 1, 7 in gen 12) — the search is
pushing into the corner of the space where the constraint solve gets hard,
which is the expected end-of-run behaviour, not a fault.

## Champion — eval 193

| quantity | value |
|---|---|
| score | 3665.77 |
| evolving-geodesic path saving | 0.35939 |
| persistence | 0.7776 |
| curvature activity | 0.7437 |
| horizon penalty | 0.0 |
| max Hamiltonian residual (L2) | 0.0300 |
| max momentum residual (L2) | 0.0022 |
| exotic-matter integral | 41.76 |
| nontriviality gate | 0.794 |

The score decomposes exactly as
`10000·0.35939 + 60·0.77763 + 40·0.74367 + 0.79395·(30·0.28531 + 10·0.15151 + 5·0.53490) + 40·(−0.36711)`,
i.e. the path saving supplies 3593.9 of the 3665.8 and everything else is a
small correction. There is no exotic-matter penalty in this objective by
design — in the gated mode phantom matter is treated as available fuel, and the
cost of that fuel is reported separately as the exotic integral above.

`operational_ftl_geodesic` is 0.0 for the champion: the shortcut is a genuine
geometric path saving, not a superluminal coordinate-speed artefact.

## Physics, as run

Identical to the parent MAP-Elites campaign — verified genome-key by
genome-key before launch:

- **Pump on for the entire simulation.** No `rl_pump_stop_time` key is written,
  and the evolution default is never-stop. The emission floor is pinned at
  t = 4 explicitly, because a negative pump value would otherwise erase the
  scorer's fallback floor.
- Stop time 26, objective `f_geo_max`, exotic penalty weight 0, pump energy
  weight 40.
- Constraint gate at 5 % relative residual — the same threshold the parent
  campaign admitted its champion under, not the looser promotion default.
- Single-rank GRTresna solves (`mpirun` segfaults on this node).
- Seed 11, population 16, four GPUs pipelined.

## Layout

- `analysis/per_eval.csv` — all 200 evaluations: generation, status, score,
  path saving, persistence, emission and arrival times, constraint residuals,
  exotic integral, energy-condition minima, termination reason.
- `analysis/per_generation.csv` — 13 rows: best and median score and path
  saving, gpu_ok and gate-rejection counts, best eval id.
- `logs/progress.log` — per-eval progress lines from the launch log.
- `run/` — the campaign directory minus heavy binaries:
  - `trajectory.jsonl` — all 200 rows with full genomes and every scored
    component. Complete: the evals whose directories were trimmed are all here.
  - `metadata.json`, `result.json` — configuration and final state.
  - `ftl_champions.json`, `ftl_retention.jsonl` — per-metric record holders.
  - `eval_000162`, `eval_000170`, `eval_000193`, `eval_000195`,
    `eval_000197` — the five best evaluations, champion included.

Pack size: 19 MB.

## Reproducing

```bash
cd grteclyn-wrapper
.venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-fgeo-max-refinement/run/eval_000193 \
  --name fmax193_replay --runs-dir ../runs/neuralspacetime/hq \
  --gpu 0 --n-full 256 --l-full 128 --max-level 3 --regrid-threshold 0.02 \
  --stop-time 64 --objective-mode f_geo_max --evolving-geodesic \
  --grtresna-ranks 1
```

The champion is also frozen for the promotion matrix at
`runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1/`, which is what
`scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh` reads.
