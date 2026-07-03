#!/usr/bin/env python3
"""Re-score a gw_beam QD campaign and rebuild archive.json from trajectory."""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from grteclyn_wrapper.core.episode import write_json
from grteclyn_wrapper.metrics import dataclass_to_dict, read_episode_metrics, score_episode
from grteclyn_wrapper.metrics.score.helpers import domain_half_width_for_episode
from grteclyn_wrapper.search.qd_search.archive import Elite, QDArchive
from grteclyn_wrapper.search.qd_search.descriptors import _bin_index, _descriptor_details
from grteclyn_wrapper.search.qd_search.io import _load_trajectory_records
from grteclyn_wrapper.search.trajectory_log import format_trajectory_line, infer_trajectory_status
from grteclyn_wrapper.search.validation_tiers import DEFAULT_TIER_CONFIG, Tier, evaluate_tiers


def _episode_dir(campaign: Path, record: dict[str, Any]) -> Path | None:
    episode = record.get("episode")
    if episode:
        path = Path(str(episode))
        if path.is_dir():
            return path
    eval_id = record.get("eval")
    if eval_id is not None:
        path = campaign / f"eval_{int(eval_id):06d}"
        if path.is_dir():
            return path
    return None


def _rescore_record(
    campaign: Path,
    record: dict[str, Any],
    *,
    bins: int,
    target_stop_time: float,
    exotic_penalty_weight: float,
) -> dict[str, Any]:
    updated = dict(record)
    if updated.get("status") != "gpu_ok" and infer_trajectory_status(updated) != "gpu_ok":
        return updated

    ep_dir = _episode_dir(campaign, updated)
    if ep_dir is None:
        updated["rescore_status"] = "missing_eval_dir"
        return updated

    overrides = updated.get("overrides") or {}
    metrics = read_episode_metrics(ep_dir, objective_mode="gw_beam")
    score = score_episode(
        metrics,
        target_stop_time=target_stop_time,
        objective_mode="gw_beam",
        domain_half_width=domain_half_width_for_episode(ep_dir, overrides),
        exotic_penalty_weight=exotic_penalty_weight,
    )
    write_json(
        ep_dir / "score.json",
        {"score": dataclass_to_dict(score), "metrics": dataclass_to_dict(metrics)},
    )

    descriptor_details = _descriptor_details(
        score.components,
        dataclass_to_dict(metrics),
        mode="gw_beam",
        overrides=overrides,
    )
    assessment = evaluate_tiers(score.components, metrics=metrics, config=DEFAULT_TIER_CONFIG)

    d1, d2 = descriptor_details["x"], descriptor_details["y"]
    updated["score"] = float(score.total)
    updated["components"] = dict(score.components)
    updated["descriptor_details"] = descriptor_details
    updated["descriptors"] = [d1, d2]
    updated["cell"] = [_bin_index(d1, bins), _bin_index(d2, bins)]
    updated["tier"] = assessment.tier
    updated["tier_name"] = assessment.tier_name
    updated["rescore_status"] = "ok"
    return updated


def _rebuild_archive(records: list[dict[str, Any]], *, bins: int) -> QDArchive:
    archive = QDArchive(bins=bins)
    dims_params: dict[str, float] = {}
    for record in records:
        if infer_trajectory_status(record) != "gpu_ok":
            continue
        if int(record.get("tier", -1)) < int(Tier.CONSTRUCTED):
            continue
        cell_raw = record.get("cell")
        if not cell_raw or len(cell_raw) != 2:
            continue
        cell = (int(cell_raw[0]), int(cell_raw[1]))
        score = record.get("score")
        if not isinstance(score, (int, float)):
            continue
        overrides = record.get("overrides") or {}
        params = {k: float(v) for k, v in overrides.items()}
        dims_params = params
        details = record.get("descriptor_details") or {}
        d1 = float(details.get("x", record.get("descriptors", [0, 0])[0]))
        d2 = float(details.get("y", record.get("descriptors", [0, 0])[1]))
        elite = Elite(
            cell=cell,
            score=float(score),
            descriptors=(d1, d2),
            params=params,
            episode=record.get("episode"),
            tier=int(record.get("tier", -1)),
            tier_name=str(record.get("tier_name", "rejected")),
            objectives=dict(record.get("objectives") or {}),
            descriptor_details=dict(details),
        )
        archive.insert(elite)
    _ = dims_params
    return archive


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign_dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    campaign = args.campaign_dir.expanduser().resolve()
    traj_path = campaign / "trajectory.jsonl"
    meta_path = campaign / "metadata.json"
    if not traj_path.is_file():
        raise SystemExit(f"missing {traj_path}")

    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.is_file() else {}
    bins = int(meta.get("bins", 8))
    target_stop_time = float((meta.get("base_overrides") or {}).get("stop_time", 24.0))
    exotic_weight = float(meta.get("exotic_penalty_weight", 0.0))

    records = _load_trajectory_records(traj_path)
    rescored: list[dict[str, Any]] = []
    ok = 0
    for record in records:
        new_record = _rescore_record(
            campaign,
            record,
            bins=bins,
            target_stop_time=target_stop_time,
            exotic_penalty_weight=exotic_weight,
        )
        if new_record.get("rescore_status") == "ok":
            ok += 1
        rescored.append(new_record)

    archive = _rebuild_archive(rescored, bins=bins)
    best = archive.best.score if archive.best else float("nan")

    print(f"Rescored {ok}/{len(records)} trajectory rows; archive elites={len(archive.cells)} best={best:.4f}")

    scored = [
        r for r in rescored
        if isinstance(r.get("score"), (int, float)) and r.get("rescore_status") == "ok"
    ]
    scored.sort(key=lambda r: float(r["score"]), reverse=True)
    print("\nTop 5 after rescore:")
    for row in scored[:5]:
        print(f"  eval {row['eval']:4d}  score={float(row['score']):8.2f}  tier={row.get('tier_name')}")

    if args.dry_run:
        print("dry-run: no files written")
        return

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    shutil.copy2(traj_path, campaign / f"trajectory.jsonl.bak.{stamp}")
    archive_path = campaign / "archive.json"
    if archive_path.is_file():
        shutil.copy2(archive_path, campaign / f"archive.json.bak.{stamp}")

    with traj_path.open("w", encoding="utf-8") as fh:
        for record in rescored:
            fh.write(format_trajectory_line(record))
    write_json(archive_path, archive.to_dict())
    write_json(campaign / "gw_beam_rescore.json", {
        "rescored_at": stamp,
        "rows": len(records),
        "rescored_ok": ok,
        "archive_elites": len(archive.cells),
        "best_score": best,
    })
    print(f"Wrote {traj_path} and {archive_path}")


if __name__ == "__main__":
    main()
