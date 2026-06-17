"""Render MAP-Elites QD batch-improvement / saturation figure."""

from __future__ import annotations

from pathlib import Path
from typing import List

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from grteclyn_wrapper.visualisation.search.trajectory_batches import (
    BatchStats,
    rolling_mean,
)

matplotlib.use("Agg")

_FIGURES_DIR = (
    Path(__file__).resolve().parents[5] / "research" / "neuralspacetime" / "figures"
)

# ---------------------------------------------------------------------------
# Publication style (consistent with project: serif + stix mathtext)
# ---------------------------------------------------------------------------

_SCIENTIFIC_STYLE = {
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
    "mathtext.fontset": "stix",
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.linewidth": 1.0,
    "grid.alpha": 0.30,
    "grid.linewidth": 0.6,
    "legend.fontsize": 11,
    "legend.framealpha": 0.85,
    "legend.edgecolor": "0.8",
}

# Colour palette
_CLR_GAIN = "#2e7d32"  # green for positive gains
_CLR_GAIN_FILL = "#a5d6a7"  # lighter green fill
_CLR_TREND = "#c62828"  # deep red for rolling mean
_CLR_CUM = "#1565c0"  # deep blue for cumulative curve
_CLR_CUM_FILL = "#bbdefb"  # light blue fill
_CLR_LAST = "#6a1b9a"  # purple for "last gain" annotation


def _detect_optimizer(campaign_dir: Path) -> str:
    """Detect optimizer type from campaign path name."""
    name = campaign_dir.name.lower()
    parent = campaign_dir.parent.name.lower() if campaign_dir.parent else ""
    if "cmaes" in name or "cma_es" in name or "cmaes" in parent:
        return "cmaes"
    return "qd"


def default_output_path(campaign_dir: Path) -> Path:
    opt = _detect_optimizer(campaign_dir)
    prefix = "cmaes_progress" if opt == "cmaes" else "qd_batch_progress"
    return _FIGURES_DIR / f"{prefix}_{campaign_dir.name}.png"


def _despine(ax: plt.Axes) -> None:
    """Remove top and right spines for a cleaner look."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def render_qd_batch_progress(
    rows: list[BatchStats],
    *,
    title: str,
    out_path: Path,
    rolling: int = 3,
    xlabel: str | None = None,
) -> List[Path]:
    """Render a publication-quality two-panel search-progress figure.

    Returns a list of paths written (PNG and PDF).
    """
    if not rows:
        raise ValueError("no batch rows to plot")

    plt.rcParams.update(_SCIENTIFIC_STYLE)

    batches = np.array([r.batch for r in rows])
    deltas = np.array([r.delta for r in rows])
    cum = np.array([r.cum_best for r in rows])
    trend = np.array(rolling_mean(deltas.tolist(), rolling))

    # ----- Figure layout -----
    fig, (ax_delta, ax_cum) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.25},
    )

    # ===================================================================
    # Upper panel: per-batch archive gain (step + fill)
    # ===================================================================
    ax_delta.step(
        batches,
        deltas,
        where="mid",
        color=_CLR_GAIN,
        linewidth=1.2,
        zorder=2,
        label=r"$\Delta f^*_k$",
    )
    ax_delta.fill_between(
        batches,
        0,
        deltas,
        step="mid",
        color=_CLR_GAIN_FILL,
        alpha=0.55,
        zorder=1,
    )

    # Rolling-mean trend
    ax_delta.plot(
        batches,
        trend,
        color=_CLR_TREND,
        linewidth=1.8,
        marker="o",
        markersize=3.5,
        markeredgewidth=0,
        label=rf"{rolling}-batch rolling mean $\langle\Delta f^*\rangle$",
        zorder=4,
    )

    ax_delta.axhline(0.0, color="0.4", linewidth=0.7, zorder=0)
    ax_delta.set_ylabel(r"$\Delta f^*$ (archive gain)")
    ax_delta.set_title(title, fontsize=13, loc="left", pad=10)
    ax_delta.grid(axis="y", linestyle="--")
    ax_delta.legend(loc="upper left", fontsize=10)
    _despine(ax_delta)

    # "Last gain" annotation
    last_gain = max((r.batch for r in rows if r.delta > 0), default=None)
    if last_gain is not None:
        ax_delta.axvline(
            last_gain,
            color=_CLR_LAST,
            linestyle="--",
            linewidth=1.0,
            alpha=0.75,
            zorder=3,
        )
        # Place annotation below the x-axis via the axes transform to avoid overlap
        ax_delta.text(
            last_gain,
            -0.08,
            rf"last gain $k\!={last_gain}$",
            transform=ax_delta.get_xaxis_transform(),
            fontsize=8.5,
            color=_CLR_LAST,
            ha="center",
            va="top",
        )

    # ===================================================================
    # Lower panel: cumulative best score (step staircase)
    # ===================================================================
    ax_cum.step(
        batches,
        cum,
        where="mid",
        color=_CLR_CUM,
        linewidth=2.0,
        zorder=3,
        label=r"$f^*_{\mathrm{cum}}$",
    )
    ax_cum.fill_between(
        batches,
        cum.min() * 0.98 if cum.min() > 0 else 0,
        cum,
        step="mid",
        color=_CLR_CUM_FILL,
        alpha=0.40,
        zorder=1,
    )

    ax_cum.set_xlabel(xlabel or r"Batch index $k$")
    ax_cum.set_ylabel(r"$f^*_{\mathrm{cum}}$ (cumulative best)")
    ax_cum.grid(axis="both", linestyle="--")
    ax_cum.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.0f}"))
    ax_cum.legend(loc="upper left", fontsize=10)
    _despine(ax_cum)

    # Minor ticks for both panels
    for ax in (ax_delta, ax_cum):
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.tick_params(which="minor", length=3, width=0.6)
        ax.tick_params(which="major", length=5, width=1.0)

    # ----- Save -----
    out_path = out_path.expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    paths_written: List[Path] = []

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    paths_written.append(out_path)

    pdf_path = out_path.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight")
    paths_written.append(pdf_path)

    plt.close(fig)
    return paths_written
