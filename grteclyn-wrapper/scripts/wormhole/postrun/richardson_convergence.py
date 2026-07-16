#!/usr/bin/env python3
"""Three-resolution convergence and large-box waveform diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
    }
)


DX = {"coarse": 2.0 / 3.0, "medium": 0.5, "fine": 0.4}
RUN_SUFFIXES = {
    "coarse": (
        "evo_omega_p0p25_m1_kappa_1p00_dx0p667_ml3_mass0p5_qball_"
        "lam170_mu614450_conv_dx067"
    ),
    "medium": (
        "evo_omega_p0p25_m1_kappa_1p00_dx0p5_ml3_mass0p5_qball_"
        "lam170_mu614450_torus_wh_collapse_ml3_gw"
    ),
    "fine": (
        "evo_omega_p0p25_m1_kappa_1p00_dx0p4_ml3_mass0p5_qball_"
        "lam170_mu614450_conv_dx040"
    ),
    "bigbox": (
        "evo_omega_p0p25_m1_kappa_1p00_dx0p5_ml2_L128_mass0p5_qball_"
        "lam170_mu614450_bigbox_L128_ml2"
    ),
}
LABELS = {
    "coarse": r"$h=0.667$",
    "medium": r"$h=0.5$",
    "fine": r"$h=0.4$",
}
STYLES = {
    "coarse": {"color": "0.55", "ls": "--"},
    "medium": {"color": "tab:blue", "ls": "-."},
    "fine": {"color": "black", "ls": "-"},
}


def load_run(run_dir: Path) -> dict[str, np.ndarray]:
    output = run_dir / "output"
    files = {
        "collapse": output / "data" / "collapse_diagnostics.dat",
        "constraints": output / "data" / "constraint_norms.dat",
        "psi4": output / "small_data" / "psi4_mode_l2m0.dat",
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing run data:\n" + "\n".join(missing))
    return {name: np.loadtxt(path) for name, path in files.items()}


def observed_order(values: dict[str, float]) -> float | None:
    """Solve the unequal-spacing Richardson relation for positive p."""
    qc, qm, qf = (values[key] for key in ("coarse", "medium", "fine"))
    dcm, dmf = qc - qm, qm - qf
    if dcm * dmf <= 0.0 or abs(dmf) < 1.0e-15:
        return None
    target = dcm / dmf
    hc, hm, hf = (DX[key] for key in ("coarse", "medium", "fine"))

    def residual(order: float) -> float:
        ratio = (hc**order - hm**order) / (hm**order - hf**order)
        return ratio - target

    lo, hi = 1.0e-6, 20.0
    flo, fhi = residual(lo), residual(hi)
    if flo * fhi > 0.0:
        return None
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fmid = residual(mid)
        if flo * fmid <= 0.0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def richardson_result(values: dict[str, float]) -> dict[str, float | None]:
    order = observed_order(values)
    result: dict[str, float | None] = {
        "coarse": values["coarse"],
        "medium": values["medium"],
        "fine": values["fine"],
        "order": order,
        "continuum": None,
        "gci_fine_percent": None,
    }
    if order is None:
        return result
    qm, qf = values["medium"], values["fine"]
    ratio = (DX["medium"] / DX["fine"]) ** order
    continuum = qf + (qf - qm) / (ratio - 1.0)
    result["continuum"] = continuum
    if abs(qf) > 0.0:
        result["gci_fine_percent"] = (
            125.0 * abs(qf - qm) / (abs(qf) * (ratio - 1.0))
        )
    return result


def scalar_results(runs: dict[str, dict[str, np.ndarray]]) -> dict[str, dict]:
    collapse = {key: runs[key]["collapse"] for key in DX}
    constraints = {key: runs[key]["constraints"] for key in DX}
    values = {
        "max_abs_K": {key: np.max(data[:, 3]) for key, data in collapse.items()},
        "max_ah_proxy": {key: np.max(data[:, 7]) for key, data in collapse.items()},
        "time_max_ah": {
            key: data[np.argmax(data[:, 7]), 0] for key, data in collapse.items()
        },
        "min_lapse": {key: np.min(data[:, 1]) for key, data in collapse.items()},
        "min_chi": {key: np.min(data[:, 2]) for key, data in collapse.items()},
        "hamiltonian_max": {
            key: np.max(data[:, 1]) for key, data in constraints.items()
        },
        "momentum_max": {
            key: np.max(data[:, 2]) for key, data in constraints.items()
        },
    }
    results = {name: richardson_result(series) for name, series in values.items()}
    for rejected in ("min_lapse", "min_chi"):
        continuum = results[rejected]["continuum"]
        if continuum is not None and continuum <= 0.0:
            results[rejected]["order"] = None
            results[rejected]["continuum"] = None
            results[rejected]["gci_fine_percent"] = None
    return results


def plot_resolution_curves(
    ax: plt.Axes,
    runs: dict[str, dict[str, np.ndarray]],
    group: str,
    column: int,
    ylabel: str,
    log: bool = False,
) -> None:
    for key in DX:
        data = runs[key][group]
        ax.plot(
            data[:, 0],
            data[:, column],
            label=LABELS[key],
            lw=1.25,
            **STYLES[key],
        )
    if log:
        ax.set_yscale("log")
    ax.set_xlim(0.0, 40.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(ylabel)
    ax.axvspan(8.0, 10.0, color="0.9", zorder=-2)
    ax.grid(alpha=0.25)


def make_convergence_figure(
    runs: dict[str, dict[str, np.ndarray]],
    results: dict[str, dict],
    out_dir: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), sharex=True)
    plot_resolution_curves(
        axes[0, 0], runs, "collapse", 1, r"$\min\alpha$", log=True
    )
    plot_resolution_curves(
        axes[0, 1], runs, "collapse", 7, r"$r_{\rm AH}$ proxy"
    )
    plot_resolution_curves(
        axes[1, 0], runs, "constraints", 1, r"$\|\mathcal{H}\|_2$", log=True
    )
    plot_resolution_curves(
        axes[1, 1], runs, "constraints", 2, r"$\|\mathcal{M}\|_2$", log=True
    )
    axes[0, 0].set_title("(a) Lapse collapse")
    axes[0, 1].set_title("(b) Trapped-surface proxy")
    axes[1, 0].set_title("(c) Hamiltonian constraint")
    axes[1, 1].set_title("(d) Momentum constraint")
    axes[0, 0].legend(frameon=False, ncol=3, fontsize=8)
    momentum = results["momentum_max"]
    ah = results["max_ah_proxy"]
    note = (
        rf"$p_{{\max\|\mathcal{{M}}\|}}={momentum['order']:.2f}$; "
        rf"$p_{{\max r_{{\rm AH}}}}={ah['order']:.2f}$"
    )
    fig.text(0.5, 0.01, note, ha="center", fontsize=9)
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 1.0))
    save_figure(fig, out_dir / "fig_convergence")


def make_headline_figure(run: dict[str, np.ndarray], out_dir: Path) -> None:
    collapse = run["collapse"]
    panels = (
        (1, r"$\min\alpha$", r"(a) Minimum lapse", True),
        (2, r"$\min\chi$", r"(b) Minimum conformal factor", True),
        (3, r"$\max |K|$", r"(c) Maximum curvature", True),
        (7, r"$r_{\rm AH}$ proxy", r"(d) Trapped-surface proxy", False),
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.2), sharex=True)
    for ax, (column, ylabel, title, log_scale) in zip(axes.flat, panels):
        ax.plot(collapse[:, 0], collapse[:, column], color="black", lw=1.3)
        ax.axvspan(8.0, 10.0, color="0.9", zorder=-2)
        ax.set(xlim=(0.0, 40.0), ylabel=ylabel, title=title)
        if log_scale:
            ax.set_yscale("log")
        ax.grid(alpha=0.25)
    for ax in axes[1]:
        ax.set_xlabel(r"$t$")
    fig.tight_layout()
    save_figure(fig, out_dir / "fig_collapse_diagnostics")


def contiguous_segments(time: np.ndarray) -> list[slice]:
    if len(time) < 2:
        return [slice(0, len(time))]
    cadence = np.median(np.diff(time))
    breaks = np.where(np.diff(time) > 1.5 * cadence)[0] + 1
    bounds = np.concatenate(([0], breaks, [len(time)]))
    return [slice(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]


def plot_waveform_with_gaps(
    ax: plt.Axes, time: np.ndarray, wave: np.ndarray, key: str
) -> None:
    for index, segment in enumerate(contiguous_segments(time)):
        ax.plot(
            time[segment] - 12.0,
            wave[segment],
            label=LABELS[key] if index == 0 else None,
            lw=1.35,
            **STYLES[key],
        )


def normalized_overlap(x: np.ndarray, y: np.ndarray) -> float:
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / denominator) if denominator else float("nan")


def bigbox_alignment(psi4: np.ndarray) -> dict[str, float]:
    time = psi4[:, 0]
    wave36, wave48 = psi4[:, 1], psi4[:, 3]
    u36, u48 = time - 36.0, time - 48.0
    shifts = np.linspace(-2.0, 2.0, 401)
    common = np.linspace(-8.0, 30.0, 381)
    reference = np.interp(common, u36, wave36)
    overlaps = []
    for shift in shifts:
        candidate = np.interp(common, u48 + shift, wave48)
        overlaps.append(normalized_overlap(reference, candidate))
    best = int(np.nanargmax(overlaps))
    shifted = np.interp(common, u48 + shifts[best], wave48)
    residual = np.linalg.norm(reference - shifted) / np.linalg.norm(reference)
    return {
        "time_shift": float(shifts[best]),
        "overlap": float(overlaps[best]),
        "relative_l2_residual": float(residual),
    }


def make_waveform_figure(
    runs: dict[str, dict[str, np.ndarray]], out_dir: Path
) -> dict[str, float]:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    for key in DX:
        psi4 = runs[key]["psi4"]
        plot_waveform_with_gaps(axes[0], psi4[:, 0], psi4[:, 1], key)
    axes[0].axvspan(0.0, 8.0, color="tab:red", alpha=0.08)
    axes[0].set(
        xlabel=r"$u=t-R_{\rm ext}$",
        ylabel=r"$r\,{\rm Re}\,\Psi_4^{2,0}$",
        title=r"(a) Resolution ladder, $R_{\rm ext}=12$",
    )
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.25)

    bigbox = runs["bigbox"]["psi4"]
    time = bigbox[:, 0]
    axes[1].plot(time - 36.0, bigbox[:, 1], label=r"$R=36$", color="black")
    axes[1].plot(
        time - 48.0, bigbox[:, 3], label=r"$R=48$", color="tab:blue", ls="--"
    )
    alignment = bigbox_alignment(bigbox)
    axes[1].set(
        xlabel=r"$u=t-R_{\rm ext}$",
        ylabel=r"$r\,{\rm Re}\,\Psi_4^{2,0}$",
        title=rf"(b) $L=128$: overlap ${alignment['overlap']:.3f}$",
        xlim=(-15.0, 35.0),
    )
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    save_figure(fig, out_dir / "fig_waveform_systematics")
    return alignment


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    script = Path(__file__).resolve()
    simulation_root = script.parents[4]
    parser.add_argument(
        "--runs-root",
        type=Path,
        default=simulation_root / "runs" / "rotating_wormhole",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=(
            simulation_root
            / "research"
            / "rotatingwormhole"
            / "article"
            / "figures"
        ),
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    runs = {
        key: load_run(args.runs_root / suffix)
        for key, suffix in RUN_SUFFIXES.items()
    }
    results = scalar_results(runs)
    make_headline_figure(runs["medium"], args.out)
    make_convergence_figure(runs, results, args.out)
    alignment = make_waveform_figure(runs, args.out)
    summary = {
        "grid_spacings": DX,
        "richardson": results,
        "bigbox_waveform": alignment,
        "notes": {
            "psi4_quantity": "Files already store r*Psi4.",
            "waveform_order": (
                "Rejected: coarse extraction misses the principal burst interval."
            ),
            "bigbox_role": (
                "Intermediate-zone extraction check, not a resolution rung."
            ),
        },
    }
    print(json.dumps(summary, indent=2))
    print(f"Saved article figures in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
