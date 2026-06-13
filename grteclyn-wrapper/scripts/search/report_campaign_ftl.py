#!/usr/bin/env python3
"""Print peak per-frame FTL metrics for a saved QD / search campaign.

Peaks are computed by scanning every frame in ``ftl_timeseries.dat`` (not the
last frame alone). Includes pruned evals via ``trajectory.jsonl``.

Examples
--------
    uv run python scripts/search/report_campaign_ftl.py runs/grtresna_qd/ftl_discovery_v15
    uv run python scripts/search/report_campaign_ftl.py runs/grtresna_qd/ftl_discovery_v15 \\
        --sort max_speed --top 30 --status gpu_ok
    uv run python scripts/search/report_campaign_ftl.py runs/grtresna_qd/ftl_discovery_v15 \\
        --curve 127 --curve-top 3
"""

from __future__ import annotations

import argparse
from pathlib import Path

from grteclyn_wrapper.search.ftl_campaign_report import (
    FtlSortKey,
    format_campaign_report,
    format_ftl_curve,
    load_campaign_ftl_summaries,
    rank_summaries,
)


def _parse_sort_key(value: str) -> FtlSortKey:
    try:
        return FtlSortKey(value)
    except ValueError as exc:
        choices = ", ".join(k.value for k in FtlSortKey)
        raise argparse.ArgumentTypeError(f"unknown sort key {value!r}; choose: {choices}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("campaign_dir", type=Path, help="campaign dir with trajectory.jsonl and eval_*")
    parser.add_argument("--top", type=int, default=20, help="number of rows to print (default: 20)")
    parser.add_argument(
        "--sort",
        type=_parse_sort_key,
        default=FtlSortKey.F_GEO_PEAK,
        help="rank key: f_geo_peak, f_op_peak, max_speed, superluminal, lifetime, timeavg, score",
    )
    parser.add_argument("--status", default=None, help="filter to trajectory status, e.g. gpu_ok")
    parser.add_argument(
        "--min-f-geo",
        type=float,
        default=0.0,
        help="minimum trustworthy f_geo peak to include (default: 0)",
    )
    parser.add_argument(
        "--no-trajectory",
        action="store_true",
        help="only report eval dirs still on disk (skip trajectory.jsonl fallback)",
    )
    parser.add_argument(
        "--curve",
        type=int,
        nargs="*",
        metavar="EVAL",
        help="print per-frame f_geo table for eval id(s)",
    )
    parser.add_argument(
        "--curve-top",
        type=int,
        default=0,
        help="also print per-frame tables for top N ranked evals",
    )
    args = parser.parse_args()

    campaign_dir = args.campaign_dir.expanduser().resolve()
    all_summaries = load_campaign_ftl_summaries(
        campaign_dir,
        include_trajectory=not args.no_trajectory,
    )
    ranked = rank_summaries(
        all_summaries,
        sort_by=args.sort,
        status=args.status,
        min_f_geo_peak=args.min_f_geo,
    )
    print(
        format_campaign_report(
            ranked,
            campaign_dir=campaign_dir,
            sort_by=args.sort,
            top_n=max(args.top, 0),
            total_in_campaign=len(all_summaries),
        )
    )

    curve_ids: list[int] = list(args.curve or [])
    if args.curve_top > 0:
        curve_ids.extend(s.eval_id for s in ranked[: args.curve_top])
    if not curve_ids:
        return

    by_id = {s.eval_id: s for s in all_summaries}
    seen: set[int] = set()
    for eval_id in curve_ids:
        if eval_id in seen:
            continue
        seen.add(eval_id)
        summary = by_id.get(eval_id)
        if summary is None:
            print(f"\neval {eval_id}: not found in campaign")
            continue
        print()
        print(format_ftl_curve(summary))


if __name__ == "__main__":
    main()
