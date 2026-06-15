#!/usr/bin/env python3
"""Plot metric-guesser radial profiles -- one candidate per figure.

Run from the grteclyn-wrapper directory:

    uv run --extra plots python tests/plot_metric_profiles.py --output-dir validation_out

Each candidate produces a clear 2x2 figure explaining the full pipeline:
  top-left:     the guessed geometry  a(r) = 1/sqrt(chi)
  top-right:    the raw conformal factor  chi(r)
  bottom-left:  the derived scalar field  phi(r)
  bottom-right: the required energy density  rho_req(r)  (the physics verdict)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from grteclyn_wrapper.initial_data.constrained_recipe import RecipeBasis, constrained_phi
from grteclyn_wrapper.initial_data.validate_guesser import (
    CandidateSpec,
    generate_candidates,
    validate_candidate,
)

SHOWCASE = {
    "accepted_phantom": [
        ("random_000", "Random bumpy geometry (phantom matter allowed)"),
        ("bubble_wall_016", "Alcubierre-style bubble wall (phantom matter allowed)"),
    ],
    "rejected_exotic_energy": [
        ("random_001", "Random bumpy geometry (normal matter only)"),
        ("bubble_wall_015", "Alcubierre-style bubble wall (normal matter only)"),
    ],
    "rejected_negative_chi": [
        ("bad_control_024", "Deliberately broken: chi goes negative"),
    ],
}


def _setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "mathtext.fontset": "cm",
            "font.family": "serif",
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "axes.unicode_minus": False,
        }
    )
    return plt


def evaluate_profiles(spec: CandidateSpec):
    """Evaluate chi, a(r), and optionally the constrained phi/rho profiles."""
    basis = RecipeBasis(
        num_bases=spec.num_bases,
        basis_width=spec.basis_width,
        basis_radius_max=spec.basis_radius_max,
    )
    r_max = 2.0 * spec.basis_radius_max
    r = np.linspace(r_max / 2048.0, r_max, 2048)
    chi = basis.evaluate(r, spec.chi_asymptotic, spec.chi_coeffs)
    chi_safe = np.clip(chi, 1.0e-12, None)
    stretch = 1.0 / np.sqrt(chi_safe)

    phi_profile = rho_required = None
    if spec.apply_constrained and float(np.min(chi)) > 0.0:
        result = constrained_phi(
            basis,
            spec.chi_asymptotic,
            spec.chi_coeffs,
            K_const=spec.K_asymptotic,
            phantom=spec.field_type == "phantom",
        )
        phi_profile = (result.r, result.phi_profile)
        rho_required = (result.r, result.rho_required)

    return r, chi, stretch, phi_profile, rho_required


def plot_single_candidate(
    spec: CandidateSpec,
    *,
    human_label: str,
    output_dir: Path,
) -> Path:
    plt = _setup_matplotlib()

    record = validate_candidate(spec)
    r, chi, stretch, phi_data, rho_data = evaluate_profiles(spec)

    is_accepted = record.classification.startswith("accepted")
    verdict_color = "#1a7a2e" if is_accepted else "#b71c1c"
    verdict_word = "ACCEPTED" if is_accepted else "REJECTED"
    verdict_text = f"{verdict_word}: {record.classification_reason}"

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"{human_label}\n"
        f"[{spec.candidate_id}  |  family: {spec.family}  |  "
        f"field: {spec.field_type}]",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.005,
        verdict_text,
        ha="center",
        fontsize=13,
        fontweight="bold",
        color=verdict_color,
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor=verdict_color,
            alpha=0.12,
            edgecolor=verdict_color,
        ),
    )

    # -- top-left: geometry stretch a(r) --
    ax = axes[0, 0]
    ax.set_title(r"Step 1:  Guessed geometry — radial stretch $a(r)$")
    ax.plot(r, stretch, color="#7b1fa2", linewidth=2.2)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.0, label="flat space")
    ax.fill_between(
        r,
        1.0,
        stretch,
        where=stretch > 1.0,
        alpha=0.10,
        color="#7b1fa2",
        label="stretched",
    )
    ax.fill_between(
        r,
        stretch,
        1.0,
        where=stretch < 1.0,
        alpha=0.10,
        color="#e65100",
        label="compressed",
    )
    ax.set_xlabel(r"radius $r$")
    ax.set_ylabel(r"$a(r) = \chi^{-1/2}$")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.2)
    ax.annotate(
        f"max stretch = {np.max(stretch):.3f}",
        xy=(r[np.argmax(stretch)], np.max(stretch)),
        xytext=(0.55, 0.85),
        textcoords="axes fraction",
        fontsize=10,
        arrowprops=dict(arrowstyle="->", color="gray"),
    )

    # -- top-right: raw conformal factor chi(r) --
    ax = axes[0, 1]
    ax.set_title(r"Internal variable: conformal factor $\chi(r)$")
    ax.plot(r, chi, color="#1565c0", linewidth=2.2)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.0, label="flat space")
    ax.axhline(0.0, color="black", linewidth=0.8)
    if float(np.min(chi)) < 0:
        ax.fill_between(
            r,
            chi,
            0.0,
            where=chi < 0.0,
            alpha=0.25,
            color="tab:red",
            label=r"$\chi < 0$ (invalid!)",
        )
    ax.set_xlabel(r"radius $r$")
    ax.set_ylabel(r"$\chi(r)$")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.2)
    ax.text(
        0.02,
        0.05,
        f"range: [{record.chi_min:.3f},  {record.chi_max:.3f}]",
        transform=ax.transAxes,
        fontsize=10,
    )

    # -- bottom-left: derived phi(r) --
    ax = axes[1, 0]
    if phi_data is not None:
        r_phi, phi_vals = phi_data
        ax.set_title(
            r"Step 2:  Scalar field $\phi(r)$ derived from constraint"
        )
        ax.plot(r_phi, phi_vals, color="#e65100", linewidth=2.2)
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_xlabel(r"radius $r$")
        ax.set_ylabel(r"$\phi(r)$")
        ax.grid(True, alpha=0.2)
        ax.text(
            0.02,
            0.05,
            rf"Gaussian fit error: {record.phi_fit_residual:.4f}",
            transform=ax.transAxes,
            fontsize=10,
        )
    else:
        ax.set_title(r"Step 2:  (skipped — $\chi$ invalid)")
        ax.text(
            0.5,
            0.5,
            r"Cannot derive $\phi$" + "\n" + r"because $\chi < 0$",
            ha="center",
            va="center",
            fontsize=14,
            color="tab:red",
        )
        ax.set_xlabel(r"radius $r$")

    # -- bottom-right: required energy density --
    ax = axes[1, 1]
    if rho_data is not None:
        r_rho, rho_vals = rho_data
        ax.set_title(
            r"Step 3:  Required energy density $\rho_{\rm req}(r)$  "
            + r"$\leftarrow$ THE PHYSICS VERDICT"
        )
        ax.plot(r_rho, rho_vals, color="#2e7d32", linewidth=2.2)
        ax.axhline(0.0, color="black", linewidth=1.0)
        ax.fill_between(
            r_rho,
            rho_vals,
            0.0,
            where=rho_vals >= 0.0,
            alpha=0.15,
            color="#2e7d32",
            label="normal matter OK",
        )
        ax.fill_between(
            r_rho,
            rho_vals,
            0.0,
            where=rho_vals < 0.0,
            alpha=0.30,
            color="#c62828",
            label="EXOTIC MATTER needed",
        )
        ax.set_xlabel(r"radius $r$")
        ax.set_ylabel(r"$\rho_{\rm req}(r)$")
        ax.legend(fontsize=10, loc="upper right")
        ax.grid(True, alpha=0.2)

        rho_min = float(np.min(rho_vals))
        if rho_min < 0:
            idx_min = int(np.argmin(rho_vals))
            ax.annotate(
                f"min = {rho_min:.2e}",
                xy=(r_rho[idx_min], rho_min),
                xytext=(0.55, 0.15),
                textcoords="axes fraction",
                fontsize=11,
                fontweight="bold",
                color="#c62828",
                arrowprops=dict(arrowstyle="->", color="#c62828", lw=1.5),
            )
    else:
        ax.set_title(r"Step 3:  (skipped — $\chi$ invalid)")
        ax.text(
            0.5,
            0.5,
            "Cannot evaluate energy\n"
            r"because $\chi < 0$",
            ha="center",
            va="center",
            fontsize=14,
            color="tab:red",
        )
        ax.set_xlabel(r"radius $r$")

    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.94))

    safe_id = spec.candidate_id.replace("/", "_")
    output_path = output_dir / f"{safe_id}.png"
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot metric-guesser radial profiles (one figure per candidate)."
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_out"),
    )
    parser.add_argument(
        "--candidate-id",
        action="append",
        dest="candidate_ids",
        default=None,
        help="Plot only these candidate ids.",
    )
    args = parser.parse_args()
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = {s.candidate_id: s for s in generate_candidates(seed=args.seed)}

    if args.candidate_ids:
        pairs = [(cid, cid) for cid in args.candidate_ids]
    else:
        pairs = []
        for group_name, items in SHOWCASE.items():
            for cid, label in items:
                pairs.append((cid, label))

    for cid, label in pairs:
        spec = candidates.get(cid)
        if spec is None:
            print(f"WARNING: {cid} not found, skipping")
            continue
        path = plot_single_candidate(
            spec,
            human_label=label,
            output_dir=output_dir,
        )
        record = validate_candidate(spec)
        print(f"  {path.name:30s}  {record.classification}")

    print(f"\nAll figures in: {output_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
