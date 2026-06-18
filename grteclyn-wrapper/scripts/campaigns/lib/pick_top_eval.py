#!/usr/bin/env python3
"""CLI: print top eval id(s) from a campaign trajectory.jsonl."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from grteclyn_wrapper.search.campaign_pick import pick_top_eval_ids


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path, help="Path to trajectory.jsonl")
    parser.add_argument("--top-k", type=int, default=1)
    parser.add_argument(
        "--rank-by",
        default="score",
        help="Trajectory field to rank by (default: score)",
    )
    parser.add_argument(
        "--format",
        choices=("eval_id", "pairs"),
        default="eval_id",
        help="Output eval ids or eval/gpu candidate pairs",
    )
    parser.add_argument(
        "--gpu",
        default="0",
        help="GPU id used with --format pairs (default: 0)",
    )
    parser.add_argument(
        "--include-dry-run",
        action="store_true",
        help="Treat dry_run trajectory rows as pickable (campaign wiring tests)",
    )
    args = parser.parse_args(argv)
    try:
        eval_ids = pick_top_eval_ids(
            args.trajectory,
            top_k=args.top_k,
            rank_by=args.rank_by,
            include_dry_run=args.include_dry_run,
        )
    except ValueError as exc:
        print(f"[pick_top_eval] {exc}", file=sys.stderr)
        return 2
    if args.format == "pairs":
        print(" ".join(f"{eval_id} {args.gpu}" for eval_id in eval_ids))
    else:
        print(" ".join(str(eval_id) for eval_id in eval_ids))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
