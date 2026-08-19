#!/usr/bin/env python3
"""Derive per-eval and per-generation tables for one trajectory-search campaign.

Own module, own output files: writes ``per_eval.csv`` and ``per_generation.csv``
into ``$OUT_DIR`` and touches nothing else. Called by
``pack_search_campaigns.sh`` only when a pack has no analysis table yet, so a
hand-curated table with campaign-specific columns is never clobbered.

The column set is deliberately objective-agnostic — the same table describes an
``f_geo_max`` and an ``f_geo_depth`` campaign, so packs stay comparable — and
every value is read straight from ``trajectory.jsonl``; nothing is recomputed.

Rows carrying an eval id seen earlier (a resumed campaign re-scores the evals
that were in flight when it stopped) collapse to their best-scoring row, which
is the one the optimiser itself kept.

Usage: search_campaign_tables.py TRAJECTORY.jsonl [METADATA.json]
Env: OUT_DIR (required) — directory to write the two CSVs into.
"""

from __future__ import annotations

import csv
import json
import math
import os
import pathlib
import statistics
import sys

# (csv column, trajectory path). A path is a tuple of nested dict keys.
EVAL_COLUMNS: list[tuple[str, tuple[str, ...]]] = [
    ("status", ("status",)),
    ("score", ("score",)),
    ("f_geo_evolving", ("components", "ftl_geo_evolving")),
    ("f_geo_depth", ("components", "ftl_geo_depth")),
    ("f_geo_peak", ("components", "ftl_geo_peak")),
    ("ftl_persistence", ("components", "ftl_persistence")),
    ("operational_ftl_geodesic", ("components", "operational_ftl_geodesic")),
    ("curvature_activity", ("components", "curvature_activity")),
    ("survival", ("components", "survival")),
    ("stability", ("components", "stability")),
    ("constraint_health", ("components", "constraint_health")),
    ("horizon_penalty", ("components", "horizon_penalty")),
    ("exotic_penalty", ("components", "exotic_penalty")),
    ("pump_energy_penalty", ("components", "pump_energy_penalty")),
    ("nontriviality_gate", ("components", "nontriviality_gate")),
    ("t_emit", ("evolving_geodesic", "t_emit")),
    ("t_arrival", ("evolving_geodesic", "t_arrival")),
    ("rays_reached", ("evolving_geodesic", "n_reached")),
    ("rays_total", ("evolving_geodesic", "n_rays")),
    ("h_drift_ok", ("evolving_geodesic", "h_quality_ok")),
    ("max_ham_l2", ("constraints", "max_hamiltonian_l2")),
    ("max_mom_l2", ("constraints", "max_momentum_l2")),
    ("integral_neg_rho", ("constraints", "integral_negative_rho")),
    ("min_rho_required", ("constraints", "min_rho_required")),
    ("min_nec", ("energy_conditions", "min_nec")),
    ("stability_violation", ("stability", "violation")),
    ("final_time", ("constraints", "final_time")),
    ("termination_reason", ("termination_reason",)),
]


def dig(record: dict, path: tuple[str, ...]):
    node = record
    for key in path:
        if not isinstance(node, dict):
            return None
        node = node.get(key)
    return node


def load_rows(traj: pathlib.Path) -> list[dict]:
    """Best-scoring row per eval id, in eval order."""
    best: dict[int, dict] = {}
    for line in traj.open(encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        try:
            eval_id = int(rec.get("eval"))
        except (TypeError, ValueError):
            continue
        rec["eval"] = eval_id
        previous = best.get(eval_id)
        if previous is None:
            best[eval_id] = rec
            continue
        # A gpu_ok row always beats a rejected one; otherwise take the higher score.
        def rank(r: dict) -> tuple[int, float]:
            ok = 1 if r.get("status") == "gpu_ok" else 0
            try:
                return (ok, float(r.get("score")))
            except (TypeError, ValueError):
                return (ok, -math.inf)

        if rank(rec) > rank(previous):
            best[eval_id] = rec
    return [best[k] for k in sorted(best)]


def generation_of(eval_id: int, population: int | None) -> str:
    if not population or population <= 0:
        return ""  # steady-state search (MAP-Elites): no generations
    return str((eval_id - 1) // population + 1)


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    traj = pathlib.Path(argv[0])
    meta = json.loads(pathlib.Path(argv[1]).read_text(encoding="utf-8")) if len(argv) > 1 else {}
    out_dir = pathlib.Path(os.environ["OUT_DIR"])
    out_dir.mkdir(parents=True, exist_ok=True)

    population = meta.get("population_size")
    rows = load_rows(traj)

    per_eval = out_dir / "per_eval.csv"
    with per_eval.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["gen", "eval"] + [name for name, _ in EVAL_COLUMNS])
        for rec in rows:
            writer.writerow(
                [generation_of(rec["eval"], population), rec["eval"]]
                + [
                    "" if (value := dig(rec, path)) is None else value
                    for _, path in EVAL_COLUMNS
                ]
            )

    per_gen = out_dir / "per_generation.csv"
    if population:
        buckets: dict[int, list[dict]] = {}
        for rec in rows:
            buckets.setdefault((rec["eval"] - 1) // population + 1, []).append(rec)
        with per_gen.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                ["gen", "evals", "gpu_ok", "rejected", "best_score", "median_score",
                 "best_eval", "best_f_geo_evolving", "median_f_geo_evolving"]
            )
            for gen in sorted(buckets):
                members = buckets[gen]
                ok = [r for r in members if r.get("status") == "gpu_ok"]
                scores = [float(r["score"]) for r in ok
                          if isinstance(r.get("score"), (int, float))]
                fgeo = [v for r in ok
                        if isinstance(v := dig(r, ("components", "ftl_geo_evolving")), (int, float))]
                best = max(ok, key=lambda r: float(r.get("score") or -math.inf), default=None)
                writer.writerow([
                    gen, len(members), len(ok), len(members) - len(ok),
                    f"{max(scores):.4f}" if scores else "",
                    f"{statistics.median(scores):.4f}" if scores else "",
                    best["eval"] if best else "",
                    f"{max(fgeo):.6f}" if fgeo else "",
                    f"{statistics.median(fgeo):.6f}" if fgeo else "",
                ])
        wrote = f"{per_eval.name}, {per_gen.name}"
    else:
        per_gen.unlink(missing_ok=True)  # steady-state search has no generations
        wrote = per_eval.name

    print(f"[tables] {out_dir}: {wrote} ({len(rows)} evals)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
