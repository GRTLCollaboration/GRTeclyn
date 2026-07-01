#!/usr/bin/env python3
"""Plot the 3D trajectories (and current positions) of Q-ball lumps for an eval.

This is a faithful Python port of the C++ ``TrajectoryEvaluator`` in
``Source/GRTeclynCore/RL/TrajectoryEvaluator.hpp``.  It reads the
``trajectory_*`` keys out of an eval's ``params.txt`` and renders the full
orbit of every lump over the evolution window, plus a marker at the final
lump location.

Trajectory model (per lump k):

    r_k(t)  = max(R_MIN, R0 + v_rad*t + A_breath*sin(omega_breath*t))
    phi(t)  = phase0 + omega_rot*t
    (x_orb, y_orb) = (r_k*cos phi, r_k*sin phi)      # orbit lives in its plane
    z_center(t)    = z_amp * sin(omega_z * t)        # shared axial bob
    then the orbital plane is rotated into the lab frame by (tilt_theta, tilt_phi).

Usage:
    python plot_lump_trajectories.py <eval_dir> [--out FILE] [--tmax T] [--nt N]

If ``--out`` is omitted the figure is written into the eval dir as
``lump_trajectories.png``.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401,E402  (registers 3d proj)


def _parse_params(params_path: Path) -> dict[str, float]:
    """Parse ``key = value`` lines from a GRTeclyn params.txt into floats."""
    out: dict[str, float] = {}
    text = params_path.read_text()
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip().strip('"')
        try:
            out[key] = float(val)
        except ValueError:
            continue
    return out


def _get(params: dict[str, float], key: str, default: float) -> float:
    return params.get(key, default)


def evaluate_lump(
    t: np.ndarray,
    lump: dict[str, float],
    a_breath: float,
    omega_breath: float,
    z_amp: float,
    omega_z: float,
    r_min: float,
) -> np.ndarray:
    """Return an (N, 3) array of lab-frame positions for a single lump.

    Mirrors ``TrajectoryEvaluator::evaluate`` exactly.
    """
    breath_dr = a_breath * np.sin(omega_breath * t)
    z_center = z_amp * np.sin(omega_z * t)

    r_k = np.maximum(r_min, lump["R0"] + lump["v_rad"] * t + breath_dr)
    phi = lump["phase0"] + lump["omega_rot"] * t

    x_orb = r_k * np.cos(phi)
    y_orb = r_k * np.sin(phi)

    ct = math.cos(lump["tilt_theta"])
    st = math.sin(lump["tilt_theta"])
    cp = math.cos(lump["tilt_phi"])
    sp = math.sin(lump["tilt_phi"])

    x = cp * ct * x_orb - sp * y_orb
    y = sp * ct * x_orb + cp * y_orb
    z = -st * x_orb + z_center
    return np.column_stack([x, y, z])


def load_lumps(params: dict[str, float]) -> list[dict[str, float]]:
    n = int(_get(params, "trajectory_num_lumps", 0))
    lumps: list[dict[str, float]] = []
    for k in range(n):
        pfx = f"trajectory_lump{k}_"
        lumps.append(
            {
                "R0": _get(params, pfx + "R0", 5.0),
                "omega_rot": _get(params, pfx + "omega_rot", 0.0),
                "phase0": _get(params, pfx + "phase0", 0.0),
                "tilt_theta": _get(params, pfx + "tilt_theta", 0.0),
                "tilt_phi": _get(params, pfx + "tilt_phi", 0.0),
                "v_rad": _get(params, pfx + "v_rad", 0.0),
                "well_depth": _get(params, pfx + "well_depth", 0.05),
                # exotic >= 0.5 => phantom (negative energy) matter
                "exotic": _get(params, pfx + "exotic", 0.0),
            }
        )
    return lumps


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("eval_dir", type=Path, help="Path to the eval_XXXXXX directory")
    ap.add_argument("--out", type=Path, default=None, help="Output image path")
    ap.add_argument("--tmax", type=float, default=None, help="Max time (default: stop_time)")
    ap.add_argument("--nt", type=int, default=2000, help="Samples along each orbit")
    args = ap.parse_args()

    eval_dir: Path = args.eval_dir
    params_path = eval_dir / "params.txt"
    if not params_path.is_file():
        raise SystemExit(f"No params.txt found in {eval_dir}")

    params = _parse_params(params_path)
    lumps = load_lumps(params)
    if not lumps:
        raise SystemExit("No lumps found (trajectory_num_lumps = 0).")

    a_breath = _get(params, "trajectory_A_breath", 0.0)
    omega_breath = _get(params, "trajectory_omega_breath", 0.0)
    z_amp = _get(params, "trajectory_z_amp", 0.0)
    omega_z = _get(params, "trajectory_omega_z", 0.0)
    r_min = _get(params, "trajectory_r_min", 0.1)

    tmax = args.tmax if args.tmax is not None else _get(params, "stop_time", 16.0)
    t = np.linspace(0.0, tmax, args.nt)

    fig = plt.figure(figsize=(13, 11))
    ax = fig.add_subplot(111, projection="3d")

    cmap = plt.get_cmap("tab10")
    all_pts = []
    for k, lump in enumerate(lumps):
        pos = evaluate_lump(t, lump, a_breath, omega_breath, z_amp, omega_z, r_min)
        all_pts.append(pos)
        color = cmap(k % 10)
        is_phantom = lump["exotic"] >= 0.5
        matter = "phantom (-rho)" if is_phantom else "canonical (+rho)"
        # linestyle distinguishes matter sign at a glance
        ls = "--" if is_phantom else "-"

        ax.plot(
            pos[:, 0], pos[:, 1], pos[:, 2],
            color=color, lw=1.6, ls=ls, alpha=0.85,
            label=(
                f"lump {k}: R0={lump['R0']:.2f}, w_rot={lump['omega_rot']:+.3f}, "
                f"v_rad={lump['v_rad']:+.3f}, tilt_th={lump['tilt_theta']:.2f}  [{matter}]"
            ),
        )
        # start marker (t=0): hollow circle
        ax.scatter(*pos[0], color=color, s=60, marker="o",
                   facecolors="none", edgecolors=color, lw=1.8)
        # final position (t=tmax): filled star
        ax.scatter(*pos[-1], color=color, s=180, marker="*",
                   edgecolors="black", lw=0.6, zorder=5)
        ax.text(pos[-1, 0], pos[-1, 1], pos[-1, 2], f" {k}",
                color=color, fontsize=10, fontweight="bold")

    # origin
    ax.scatter([0], [0], [0], color="black", s=40, marker="+")

    pts = np.vstack(all_pts)
    rng = np.max(np.abs(pts)) * 1.05
    ax.set_xlim(-rng, rng)
    ax.set_ylim(-rng, rng)
    ax.set_zlim(-rng, rng)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    title = (
        f"Lump trajectories  ({eval_dir.name}),  t in [0, {tmax:g}]\n"
        f"A_breath={a_breath:.2f}, w_breath={omega_breath:.2f}, "
        f"z_amp={z_amp:.2f}, w_z={omega_z:.2f}\n"
        "hollow o = start (t=0),  * = final position;  "
        "solid = canonical, dashed = phantom"
    )
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper left", fontsize=7, bbox_to_anchor=(-0.15, 1.0))
    ax.view_init(elev=22, azim=-60)

    out = args.out if args.out is not None else eval_dir / "lump_trajectories.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
