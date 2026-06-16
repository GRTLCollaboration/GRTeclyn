"""Render MAP-Elites QD batch-improvement / saturation figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from grteclyn_wrapper.visualisation.search.trajectory_batches import (
    BatchStats,
    format_batch_summary,
    rolling_mean,
)

_PLOTS_DIR = Path(__file__).resolve().parents[1] / "plots"


def default_output_path(campaign_dir: Path) -> Path:
    return _PLOTS_DIR / f"qd_batch_progress_{campaign_dir.name}.png"


def render_qd_batch_progress(
    rows: list[BatchStats],
    *,
    title: str,
    out_path: Path,
    rolling: int = 3,
) -> Path:
    if not rows:
        raise ValueError("no batch rows to plot")

    batches = [r.batch for r in rows]
    deltas = [r.delta for r in rows]
    cum = [r.cum_best for r in rows]
    trend = rolling_mean(deltas, rolling)

    fig, (ax_delta, ax_cum) = plt.subplots(
        2,
        1,
        figsize=(9.0, 6.8),
        sharex=True,
        gridspec_kw={"height_ratios": [1.1, 1.0], "hspace": 0.14},
    )

    colors = ["#2ca02c" if d > 0 else "#c7c7c7" for d in deltas]
    ax_delta.bar(batches, deltas, color=colors, edgecolor="white", linewidth=0.6, zorder=2)
    ax_delta.plot(
        batches,
        trend,
        color="#d62728",
        linewidth=2.0,
        marker="o",
        markersize=4,
        label=f"{rolling}-batch rolling mean Δ",
        zorder=3,
    )
    ax_delta.set_ylabel("Archive gain per batch")
    ax_delta.set_title(title, fontsize=11, loc="left", pad=8)
    ax_delta.grid(axis="y", alpha=0.25)
    ax_delta.legend(loc="upper right", framealpha=0.9)
    ax_delta.axhline(0.0, color="black", linewidth=0.8, zorder=1)

    last_gain = max((r.batch for r in rows if r.delta > 0), default=None)
    if last_gain is not None:
        ax_delta.axvline(last_gain, color="#9467bd", linestyle="--", linewidth=1.2, alpha=0.8)
        ymax = max(deltas) if max(deltas) > 0 else 1.0
        ax_delta.annotate(
            f"last gain (batch {last_gain})",
            xy=(last_gain, ymax * 0.55),
            xytext=(last_gain + 2.5, ymax * 0.82),
            fontsize=8,
            color="#9467bd",
            arrowprops={"arrowstyle": "->", "color": "#9467bd", "lw": 1.0},
        )

    ax_cum.plot(batches, cum, color="#1f77b4", linewidth=2.2, marker=".", markersize=5)
    ax_cum.set_xlabel("MAP-Elites batch (8 parallel evals)")
    ax_cum.set_ylabel("Cumulative best score")
    ax_cum.grid(alpha=0.25)
    ax_cum.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.0f}"))

    summary = format_batch_summary(rows)
    fig.subplots_adjust(left=0.10, right=0.96, top=0.94, bottom=0.26, hspace=0.14)
    fig.text(
        0.10,
        0.04,
        summary,
        fontsize=7.5,
        family="monospace",
        va="bottom",
        ha="left",
        transform=fig.transFigure,
        clip_on=False,
    )

    out_path = out_path.expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path
