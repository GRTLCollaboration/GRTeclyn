#!/usr/bin/env python3
"""Plot first-stage non-spherical metric proposals in physical space.

Run from grteclyn-wrapper:

    uv run --extra plots python tests/plot_nonspherical_profiles.py --output-dir validation_out

Each candidate produces a 2x2 figure:
  top-left:     physical cross-section of chi(x,z) — shows actual spatial shape
  top-right:    physical cross-section of geometry stretch a(x,z) = 1/sqrt(chi)
  bottom-left:  polar "shape" of chi at radius of peak deformation
  bottom-right: angular contrast as a function of radius
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from grteclyn_wrapper.nonspherical_guesser import (
    evaluate_chi_ray,
    generate_nonspherical_candidates,
    validate_candidate_rays,
    NonSphericalCandidate,
)


SHOWCASE = {
    "dipole_lopsided_000": "Dipole (l=1): space compressed more on one side",
    "quadrupole_bubble_001": "Quadrupole (l=2): prolate/oblate bubble shape",
    "strong_quadrupole_bad_003": "Bad control: deformation so strong chi goes negative",
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


def _evaluate_xz_slice(
    candidate: NonSphericalCandidate,
    *,
    n_grid: int = 300,
    extent: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate chi on a physical (x, z) cross-section (y=0 plane)."""
    if extent is None:
        extent = 1.8 * candidate.basis.basis_radius_max
    coords = np.linspace(-extent, extent, n_grid)
    X, Z = np.meshgrid(coords, coords)
    R = np.sqrt(X**2 + Z**2)
    R = np.clip(R, 1.0e-6, None)
    THETA = np.arccos(np.clip(Z / R, -1.0, 1.0))

    chi = np.zeros_like(R)
    for i in range(n_grid):
        for j in range(n_grid):
            r_val = R[i, j]
            th_val = THETA[i, j]
            r_arr = np.array([r_val])
            chi[i, j] = evaluate_chi_ray(candidate, r_arr, th_val)[0]

    return coords, coords, chi


def _evaluate_xz_fast(
    candidate: NonSphericalCandidate,
    *,
    n_grid: int = 300,
    extent: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fast vectorized version using angular sampling."""
    if extent is None:
        extent = 1.8 * candidate.basis.basis_radius_max
    coords = np.linspace(-extent, extent, n_grid)
    X, Z = np.meshgrid(coords, coords)
    R = np.sqrt(X**2 + Z**2)
    R_safe = np.clip(R, 1.0e-6, None)
    THETA = np.arccos(np.clip(Z / R_safe, -1.0, 1.0))

    chi_base = candidate.basis.evaluate(
        R_safe.ravel(), candidate.chi_asymptotic, candidate.chi_coeffs
    ).reshape(R.shape)

    from grteclyn_wrapper.nonspherical_guesser import legendre_p

    chi = chi_base.copy()
    for mode in candidate.modes:
        radial = np.exp(
            -((R_safe - mode.radial_center) ** 2) / (2.0 * mode.radial_width**2)
        )
        mu = np.cos(THETA)
        angular = legendre_p(mode.ell, mu)
        chi += mode.amplitude * radial * angular

    return coords, coords, chi


def _polar_chi_at_radius(
    candidate: NonSphericalCandidate,
    radius: float,
    n_angles: int = 360,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate chi as a function of angle at a fixed radius."""
    thetas = np.linspace(0.0, 2.0 * np.pi, n_angles)
    r_arr = np.array([radius])
    values = np.array(
        [evaluate_chi_ray(candidate, r_arr, th % np.pi)[0] for th in thetas]
    )
    return thetas, values


def _find_peak_deformation_radius(candidate: NonSphericalCandidate) -> float:
    """Find the radius where angular contrast is largest."""
    r_max = 2.0 * candidate.basis.basis_radius_max
    r = np.linspace(r_max / 512, r_max, 512)
    thetas = [0.0, np.pi / 4, np.pi / 2, 3.0 * np.pi / 4, np.pi]
    bundle = np.vstack([evaluate_chi_ray(candidate, r, th) for th in thetas])
    contrast = np.max(bundle, axis=0) - np.min(bundle, axis=0)
    return float(r[np.argmax(contrast)])


def plot_candidate(
    candidate_id: str,
    label: str,
    output_dir: Path,
    seed: int,
) -> Path:
    plt = _setup_matplotlib()
    candidates = {c.candidate_id: c for c in generate_nonspherical_candidates(seed=seed)}
    candidate = candidates[candidate_id]
    record = validate_candidate_rays(candidate)

    is_accepted = record.classification.startswith("accepted")
    verdict_color = "#1a7a2e" if is_accepted else "#b71c1c"
    verdict_word = "ACCEPTED" if is_accepted else "REJECTED"

    x, z, chi_2d = _evaluate_xz_fast(candidate)
    chi_safe = np.clip(chi_2d, 1.0e-12, None)
    stretch_2d = 1.0 / np.sqrt(chi_safe)

    peak_r = _find_peak_deformation_radius(candidate)
    polar_theta, polar_chi = _polar_chi_at_radius(candidate, peak_r)

    r_max = 2.0 * candidate.basis.basis_radius_max
    r_line = np.linspace(r_max / 512, r_max, 512)
    thetas_line = [0.0, np.pi / 4, np.pi / 2, 3.0 * np.pi / 4, np.pi]
    angular_contrast = (
        np.max(
            np.vstack(
                [evaluate_chi_ray(candidate, r_line, th) for th in thetas_line]
            ),
            axis=0,
        )
        - np.min(
            np.vstack(
                [evaluate_chi_ray(candidate, r_line, th) for th in thetas_line]
            ),
            axis=0,
        )
    )

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        f"{label}\n[{candidate.candidate_id}]",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.005,
        f"{verdict_word}: {record.reason}",
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

    # -- top-left: physical cross-section of chi --
    ax = fig.add_subplot(2, 2, 1)
    vmin_chi = max(float(np.min(chi_2d)), -0.5)
    vmax_chi = float(np.max(chi_2d))
    im = ax.pcolormesh(
        x,
        z,
        chi_2d,
        cmap="RdYlGn",
        shading="auto",
        vmin=vmin_chi,
        vmax=vmax_chi,
    )
    if float(np.min(chi_2d)) <= 0.0:
        ax.contour(x, z, chi_2d, levels=[0.0], colors="red", linewidths=2.5)
    ax.set_aspect("equal")
    ax.set_title(r"Physical cross-section: $\chi(x,z)$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$z$")
    fig.colorbar(im, ax=ax, label=r"$\chi$", shrink=0.85)

    # -- top-right: physical cross-section of stretch --
    ax = fig.add_subplot(2, 2, 2)
    vmax_stretch = min(float(np.max(stretch_2d)), 2.0)
    im = ax.pcolormesh(
        x,
        z,
        np.clip(stretch_2d, 0.0, vmax_stretch),
        cmap="magma",
        shading="auto",
        vmin=0.9,
        vmax=vmax_stretch,
    )
    ax.set_aspect("equal")
    ax.set_title(r"Geometry stretch: $a(x,z)=\chi^{-1/2}$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$z$")
    fig.colorbar(im, ax=ax, label=r"$a$", shrink=0.85)

    # -- bottom-left: polar plot of chi at peak deformation radius --
    ax = fig.add_subplot(2, 2, 3, projection="polar")
    ax.plot(polar_theta, polar_chi, color="#1565c0", linewidth=2.5)
    ax.fill(polar_theta, polar_chi, alpha=0.15, color="#1565c0")
    circle_flat = np.full_like(polar_theta, 1.0)
    ax.plot(polar_theta, circle_flat, color="gray", linestyle="--", linewidth=1.0)
    ax.set_title(
        rf"Angular shape of $\chi$ at $r={peak_r:.1f}$",
        pad=15,
    )
    ax.set_rlabel_position(45)

    # -- bottom-right: angular contrast vs radius --
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(r_line, angular_contrast, color="#6a1b9a", linewidth=2.2)
    ax.axvline(peak_r, color="gray", linestyle=":", linewidth=1.0, label=rf"$r={peak_r:.1f}$")
    ax.set_title(r"Angular contrast: $\max_\theta\chi - \min_\theta\chi$")
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\Delta_\theta\chi$")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=10)
    ax.text(
        0.97,
        0.95,
        "\n".join(
            [
                rf"$\min\chi={record.chi_min:+.3f}$",
                rf"$\max\chi={record.chi_max:+.3f}$",
                rf"$\max\Delta_\theta\chi={record.max_angular_contrast:.3f}$",
            ]
        ),
        transform=ax.transAxes,
        fontsize=10,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.93))
    output_path = output_dir / f"{candidate_id}.png"
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot non-spherical ray validation cases.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("validation_out"))
    parser.add_argument("--candidate-id", action="append", default=None)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.candidate_id:
        requested = {cid: cid for cid in args.candidate_id}
    else:
        requested = SHOWCASE

    for cid, label in requested.items():
        path = plot_candidate(cid, label, args.output_dir, args.seed)
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
