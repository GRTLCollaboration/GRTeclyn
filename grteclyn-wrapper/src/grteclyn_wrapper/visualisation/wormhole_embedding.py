#!/usr/bin/env python3
r"""Plot the solved rotating Q-torus used as wormhole initial matter.

The profile table stores the stationary amplitude ``f(rho,z)`` in
``Phi=f exp[i(m phi-omega t)]``.  The figure separates the axisymmetric modulus
from the winding phase; it does not depict the obsolete two-orbiting-lump setup.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
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
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }
)

SIM_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_ID_DIR = (
    SIM_ROOT
    / "runs"
    / "rotating_torus_id"
    / "torus_m1_om0p250_kappa1p00_dx0p5_L64_lam170_mu614450_exotic_throat1"
)
DEFAULT_OUT = (
    SIM_ROOT
    / "research"
    / "rotatingwormhole"
    / "article"
    / "figures"
    / "fig_initial_matter"
)


@dataclass(frozen=True)
class TorusProfile:
    rho: np.ndarray
    z: np.ndarray
    amplitude: np.ndarray
    omega: float
    winding: int
    f_max: float


def header_value(lines: list[str], name: str, cast):
    pattern = re.compile(rf"(?:^|\s){re.escape(name)}=([^\s]+)")
    for line in lines:
        match = pattern.search(line)
        if match:
            return cast(match.group(1))
    raise ValueError(f"Missing {name}=... in qball_torus.dat header")


def load_torus(path: Path) -> TorusProfile:
    if not path.is_file():
        raise FileNotFoundError(f"Q-torus profile not found: {path}")
    lines = path.read_text().splitlines()
    header = [line for line in lines if line.startswith("#")]
    rows = np.loadtxt(path)
    if rows.ndim != 2 or rows.shape[1] < 3:
        raise ValueError(f"Expected rho z f columns in {path}")
    rho_values = np.unique(rows[:, 0])
    z_values = np.unique(rows[:, 1])
    expected = len(rho_values) * len(z_values)
    if len(rows) != expected:
        raise ValueError(f"Profile is not a rectangular grid: {len(rows)} != {expected}")
    amplitude = np.empty((len(rho_values), len(z_values)))
    rho_index = {value: index for index, value in enumerate(rho_values)}
    z_index = {value: index for index, value in enumerate(z_values)}
    for rho, z, field in rows[:, :3]:
        amplitude[rho_index[rho], z_index[z]] = field
    return TorusProfile(
        rho=rho_values,
        z=z_values,
        amplitude=amplitude,
        omega=header_value(header, "omega", float),
        winding=header_value(header, "m_az", int),
        f_max=header_value(header, "f_max", float),
    )


def read_parameter(path: Path, name: str, default: float) -> float:
    if not path.is_file():
        return default
    pattern = re.compile(rf"^\s*{re.escape(name)}\s*=\s*([0-9.eE+-]+)")
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1))
    return default


def symmetric_meridional(profile: TorusProfile) -> tuple[np.ndarray, np.ndarray]:
    negative = -profile.rho[:0:-1]
    rho = np.concatenate((negative, profile.rho))
    amplitude = np.concatenate(
        (profile.amplitude[:0:-1, :], profile.amplitude), axis=0
    )
    return rho, amplitude


def make_figure(profile: TorusProfile, throat_radius: float) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.15))
    rho_sym, field_sym = symmetric_meridional(profile)
    display_limit = min(12.0, float(profile.rho[-1]), float(profile.z[-1]))

    meridional = axes[0].pcolormesh(
        rho_sym,
        profile.z,
        field_sym.T,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=profile.f_max,
    )
    axes[0].plot(
        [-throat_radius, throat_radius],
        [0.0, 0.0],
        color="cyan",
        lw=2.0,
    )
    axes[0].annotate(
        "throat",
        xy=(throat_radius, 0.0),
        xytext=(1.5, 1.2),
        color="cyan",
        fontsize=7,
        arrowprops={"arrowstyle": "->", "color": "cyan", "lw": 0.8},
    )
    axes[0].set(
        xlim=(-display_limit, display_limit),
        ylim=(0.0, display_limit),
        aspect="equal",
        xlabel=r"$\rho$",
        ylabel=r"$z$",
        title=r"(a) Meridional amplitude $f(\rho,z)$",
    )
    colorbar = fig.colorbar(meridional, ax=axes[0], fraction=0.047, pad=0.03)
    colorbar.set_label(r"$f$")

    span = display_limit
    points = 400
    x = np.linspace(-span, span, points)
    y = np.linspace(-span, span, points)
    xx, yy = np.meshgrid(x, y)
    radius = np.hypot(xx, yy)
    equatorial_profile = profile.amplitude[:, 0]
    modulus = np.interp(radius, profile.rho, equatorial_profile, right=0.0)
    phase = np.arctan2(yy, xx)
    real_field = modulus * np.cos(profile.winding * phase)

    for ax, data, cmap, title in (
        (axes[1], modulus, "magma", r"(b) Equatorial modulus $|\Phi|$"),
        (
            axes[2],
            real_field,
            "RdBu_r",
            rf"(c) Winding phase: $\mathrm{{Re}}\,\Phi$, $m={profile.winding}$",
        ),
    ):
        maximum = profile.f_max if ax is axes[1] else profile.f_max
        image = ax.pcolormesh(
            x,
            y,
            data,
            shading="auto",
            cmap=cmap,
            vmin=0.0 if ax is axes[1] else -maximum,
            vmax=maximum,
        )
        throat = plt.Circle(
            (0.0, 0.0), throat_radius, fill=False, color="cyan", lw=1.5
        )
        ax.add_patch(throat)
        ax.set(
            xlim=(-span, span),
            ylim=(-span, span),
            aspect="equal",
            xlabel=r"$x$",
            ylabel=r"$y$",
            title=title,
        )
        fig.colorbar(image, ax=ax, fraction=0.047, pad=0.03)

    fig.suptitle(
        rf"Initial exotic Q-torus: "
        rf"$\Phi=f(\rho,z)e^{{i(m\varphi-\omega t)}}$, "
        rf"$m={profile.winding}$, $\omega={profile.omega:.4f}$",
        fontsize=11,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    return fig


def save_figure(fig: plt.Figure, output: Path) -> None:
    stem = output.with_suffix("") if output.suffix else output
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {stem}.png")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id-dir", type=Path, default=DEFAULT_ID_DIR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--throat-radius",
        type=float,
        default=None,
        help="Override wormhole_throat_radius (default: read params, then 0.5)",
    )
    args = parser.parse_args()
    profile = load_torus(args.id_dir / "qball_torus.dat")
    throat_radius = (
        args.throat_radius
        if args.throat_radius is not None
        else read_parameter(args.id_dir / "params.txt", "wormhole_throat_radius", 0.5)
    )
    save_figure(make_figure(profile, throat_radius), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
