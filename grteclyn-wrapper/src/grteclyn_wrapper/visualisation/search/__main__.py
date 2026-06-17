"""CLI: QD campaign batch-progress figure from trajectory.jsonl."""

from __future__ import annotations

import argparse
from pathlib import Path

from grteclyn_wrapper.visualisation.search.qd_batch_progress import (
    _detect_optimizer,
    default_output_path,
    render_qd_batch_progress,
)
from grteclyn_wrapper.visualisation.search.trajectory_batches import (
    batch_stats_from_campaign,
    format_batch_summary,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot MAP-Elites QD score improvement per GPU batch (marginal gains + saturation). "
            "Reads trajectory.jsonl from a QD campaign directory."
        ),
    )
    parser.add_argument("campaign_dir", type=Path, help="QD campaign dir with trajectory.jsonl")
    parser.add_argument("--batch-size", type=int, default=8, help="evals per batch (default: 8)")
    parser.add_argument("--rolling", type=int, default=3, help="rolling-mean window (default: 3)")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="output PNG path (default: research/neuralspacetime/figures/qd_batch_progress_<name>.png)",
    )
    parser.add_argument("--title", default=None, help="figure title (default: campaign name)")
    args = parser.parse_args()

    campaign_dir = args.campaign_dir.expanduser().resolve()
    rows = batch_stats_from_campaign(campaign_dir, args.batch_size)
    out_path = args.out or default_output_path(campaign_dir)

    opt = _detect_optimizer(campaign_dir)
    if opt == "cmaes":
        default_title = f"CMA-ES refinement progress \u2014 {campaign_dir.name}"
        xlabel = r"CMA-ES generation batch $k$"
    else:
        default_title = f"MAP-Elites batch progress \u2014 {campaign_dir.name}"
        xlabel = r"MAP-Elites batch index $k$"

    title = args.title or default_title

    paths = render_qd_batch_progress(
        rows, title=title, out_path=out_path, rolling=args.rolling, xlabel=xlabel
    )
    print(format_batch_summary(rows))
    print()
    for p in paths:
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
