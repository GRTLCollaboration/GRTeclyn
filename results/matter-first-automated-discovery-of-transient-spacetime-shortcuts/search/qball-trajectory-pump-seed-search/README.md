# Pump-era seed search — git-friendly extract

(Campaign internal name: `qball_traj_pump_v2`.)

Light copy of `runs/neuralspacetime/search/map_elites/qball_traj_pump_v2` —
MAP-Elites, 200 evaluations (started 2026-07-30, resumed from eval 80 on
2026-08-03), descriptor axes on FTL lifetime.

## Why this pack exists

This is the **ancestor of every later campaign**, and it is kept for provenance
rather than for its numbers. Its retained elites are what
`scripts/campaigns/qball_trajectory/run_fgeo.sh` warm-starts from by default:

```bash
SEED_EVAL_DIRS="$(ls -d .../qball_traj_pump_v2/eval_*)"
```

so the gated-shortcut lineage did not begin from random throws in a
39-dimensional space — it began in the productive corner this campaign found.
Deleting it would make that starting point unreproducible.

## What it found

Best evolving-geodesic path saving **22.15 %** (eval 13, score 559.59); 128 of
200 evaluations returned a scored run, 72 were refused by the constraint gate.

These numbers are **not comparable to the later campaigns** and are not quoted
in the article. This campaign optimised the pump under an earlier scorer, before
the evolving-geodesic measurement was what it is today; its value is the
genomes, not the scores. The successor campaign re-evaluated every seed under
the current physics before using it.

## Layout

- `analysis/per_eval.csv` — all 200 evaluations. No per-generation table:
  MAP-Elites is steady-state and has no generations.
- `run/trajectory.jsonl` — every evaluation with its full genome.
- `run/archive.json` — the descriptor grid and its occupants.
- `run/eval_000013`, `eval_000059`, `eval_000068`, `eval_000104`,
  `eval_000165` — the five best retained evaluations, i.e. the warm-start
  seeds.

Two further seed directories (`eval_000002`, `eval_000004`, scoring 137 and 89)
were held live until 2026-08-19 and are not in this pack. A re-run of
`run_fgeo.sh` with the default seed set will therefore start from five seeds
rather than seven; pass `SEED_EVAL_DIRS` explicitly to pin an exact set. Both
dropped genomes remain in `run/trajectory.jsonl`.

Pack size: 22 MB.
