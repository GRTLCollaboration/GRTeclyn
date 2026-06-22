#!/usr/bin/env python3
"""Offline re-score splash / critical_collapse campaign evals with updated objective."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode


def _load_old_score(eval_dir: Path) -> float | None:
    score_path = eval_dir / "score.json"
    if not score_path.is_file():
        return None
    payload = json.loads(score_path.read_text(encoding="utf-8"))
    score = payload.get("score") or {}
    total = score.get("total")
    return float(total) if total is not None else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "campaign_dir",
        type=Path,
        help="e.g. runs/grtresna_qd/spacetime_splash_v11",
    )
    parser.add_argument("--evals", nargs="*", help="specific eval dirs (default: all)")
    parser.add_argument(
        "--splash-mode",
        default="discovery",
        choices=["discovery", "threshold"],
    )
    args = parser.parse_args()

    campaign = args.campaign_dir.expanduser().resolve()
    if args.evals:
        eval_dirs = [campaign / e for e in args.evals]
    else:
        eval_dirs = sorted(campaign.glob("eval_*"))

    rows: list[dict] = []
    for ep_dir in eval_dirs:
        row: dict = {"episode": ep_dir.name}
        old_total = _load_old_score(ep_dir)
        row["old_total"] = old_total

        if not ep_dir.is_dir():
            row["status"] = "missing"
            rows.append(row)
            continue

        try:
            metrics = read_episode_metrics(
                ep_dir, objective_mode="critical_collapse"
            )
            score = score_episode(
                metrics,
                objective_mode="critical_collapse",
                splash_mode=args.splash_mode,
            )
        except Exception as exc:  # noqa: BLE001
            row["status"] = "rescore_failed"
            row["error"] = str(exc)
            rows.append(row)
            continue

        row["status"] = "ok"
        row["new_total"] = float(score.total)
        row["delta"] = (
            float(score.total) - old_total if old_total is not None else None
        )
        row["components"] = dict(score.components)
        row["notes"] = list(score.notes)
        rows.append(row)

    scored = [r for r in rows if r.get("new_total") is not None]
    scored.sort(key=lambda r: r["new_total"], reverse=True)

    out_path = campaign / "splash_rescore.json"
    out_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")

    print(f"Wrote {out_path} ({len(rows)} evals, {len(scored)} rescored)")
    print("\nTop 10 after rescore:")
    for row in scored[:10]:
        old = row.get("old_total")
        old_s = f"{old:.1f}" if old is not None else "n/a"
        print(
            f"  {row['episode']}  new={row['new_total']:.1f}  "
            f"old={old_s}  delta={row.get('delta')}"
        )


if __name__ == "__main__":
    main()
