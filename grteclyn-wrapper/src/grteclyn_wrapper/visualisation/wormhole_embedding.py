#!/usr/bin/env python3
r"""Embedding diagram of the spinning supported wormhole + its exotic lump constellation.

Draws, to scale and in simulation (code, G=c=1) units, the two ingredients of the
``RotatingWormholeCollapse`` config-B initial data:

1. **The wormhole geometry.**  The initial data is conformally flat with
   ``chi = (4 r^2 / (4 r^2 + b0^2))^2`` (isotropic radius ``r``; see
   ``SupportedWormholeInitialData.hpp``).  GRTeclyn's convention
   ``gamma_ij = chi^{-1} delta_ij`` makes the areal radius

       R(r) = chi^{-1/2} r = r + b0^2 / (4 r),

   i.e. an Einstein-Rosen / Ellis bridge whose throat (minimal areal radius) is at
   ``r = b0/2`` with ``R_min = b0``.  The equatorial (z=0) 2-slice embeds in flat
   3-space as the surface of revolution

       Z(R) = +/- b0 * arccosh(R / b0),        R >= b0,

   the classic two-sheet wormhole funnel.  ``b0 = wormhole_throat_radius = 0.5``.

2. **The exotic matter.**  Two phantom Q-ball lumps (``num_lumps = 2``) sit at
   ``(+/-R0, 0, 0)`` in the z=0 plane and orbit the throat at angular velocity
   ``omega_orb`` (tangential speed ``v = omega_orb * R0``), carrying the orbital
   angular momentum that makes the throat *spin* (real frame-dragging J_z).  Each
   lump's modulus is the tabulated flat-space Q-ball eigenstate
   ``phi(r) = amp * phi0(r)/phi0(0)`` used verbatim by the solver
   (``profile == 3`` => no width rescale; ``mode == 0`` => isotropic), so the
   clouds are genuinely *fat* (half-max radius ~ 5.7) and overlap right through the
   tiny r=0.5 neck.

Panel A: the 3-D embedding funnel with the two lump centres (to-scale half-max
spheres), the orbit ring, the rotation sense and the J_z spin axis.
Panel B: the true equatorial |Phi| map of the two overlapping Q-ball clouds, with
the throat circle, the orbit circle and the tangential velocity arrows -- honest
proportions, axes in code units.

Usage:
    python wormhole_embedding.py [--id-dir DIR] [--out PNG]
                                 [--b0 0.5] [--orbit-radius 6] [--orbit-omega 0.1]
                                 [--amp 0.06]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402
from mpl_toolkits.mplot3d import proj3d  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

# Config-B defaults (rotwh ... kappa_0p60 ... nlump2_R6_worb0p1).
DEFAULT_ID_DIR = (
    "/home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/rotating_wormhole_id/"
    "rotwh_omega_p0p05_m1_kappa_0p60_dx0p5_L128_mass0p5_qball_lam170_mu614450_"
    "nlump2_R6_worb0p1"
)


class Arrow3D(FancyArrowPatch):
    """A 3-D arrow (matplotlib has none built in)."""

    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, _ = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return float(np.min(zs3d))


def load_qball_profile(id_dir: Path):
    """Return (r, phi0/phi0(0)) from the ID's qball_profile.dat, or None."""
    dat = id_dir / "qball_profile.dat"
    if not dat.is_file():
        return None
    r, phi = [], []
    for line in dat.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            r.append(float(parts[0]))
            phi.append(float(parts[1]))
        except ValueError:
            continue
    if len(r) < 2:
        return None
    r = np.asarray(r)
    phi = np.asarray(phi)
    return r, phi / phi[0]


def profile_radius(r: np.ndarray, env: np.ndarray, frac: float) -> float:
    """Radius where the normalised envelope first drops below `frac`."""
    below = np.where(env < frac)[0]
    return float(r[below[0]]) if below.size else float(r[-1])


def embedding_Z(R: np.ndarray, b0: float) -> np.ndarray:
    """Z(R) = b0 * arccosh(R/b0) for the equatorial ER/Ellis bridge slice."""
    return b0 * np.arccosh(np.maximum(R / b0, 1.0))


def envelope_fn(r_query, prof):
    """Normalised radial envelope phi0(r)/phi0(0); Gaussian fallback if no table."""
    if prof is None:
        # Fallback shape (only used if the table is missing).
        return np.exp(-0.5 * (r_query / 5.7) ** 2)
    r_tab, env_tab = prof
    return np.interp(r_query, r_tab, env_tab, left=env_tab[0], right=0.0)


def build_panel_A(ax, b0, R0, omega_orb, amp, prof, r_half):
    """3-D wormhole funnel, zoomed on the throat so the neck is visible."""
    R_zoom = 5.0  # zoom: the funnel is only curved within a few b0 of the neck
    Rg = np.linspace(b0, R_zoom, 120)
    th = np.linspace(0.0, 2.0 * np.pi, 120)
    Rgrid, Thgrid = np.meshgrid(Rg, th)
    Xg = Rgrid * np.cos(Thgrid)
    Yg = Rgrid * np.sin(Thgrid)
    Zg = embedding_Z(Rgrid, b0)
    for sign in (+1.0, -1.0):
        ax.plot_surface(
            Xg, Yg, sign * Zg, rstride=3, cstride=3, linewidth=0.0,
            color="0.55", alpha=0.45, shade=True, antialiased=True,
        )
    tt = np.linspace(0, 2 * np.pi, 200)
    # Throat ring (minimal circle, areal radius b0).
    ax.plot(b0 * np.cos(tt), b0 * np.sin(tt), np.zeros_like(tt),
            color="k", lw=2.2, zorder=5)
    # J_z spin axis through the throat.
    ax.add_artist(Arrow3D([0, 0], [0, 0], [-2.4, 2.4], mutation_scale=15,
                          lw=2.0, arrowstyle="-|>", color="purple"))
    ax.text(0, 0, 2.7, r"$J_z$", color="purple", ha="center", fontsize=10)
    # Rotation-sense arrows on the rim (CCW).
    for ph in np.linspace(0, 2 * np.pi, 6, endpoint=False):
        px, py = R_zoom * np.cos(ph), R_zoom * np.sin(ph)
        tx, ty = -np.sin(ph), np.cos(ph)
        pz = float(embedding_Z(np.array([R_zoom]), b0)[0])
        ax.add_artist(Arrow3D([px, px + 1.2 * tx], [py, py + 1.2 * ty],
                              [pz, pz], mutation_scale=8, lw=1.1,
                              arrowstyle="-|>", color="darkgreen"))
    ax.text(0, 0, 0.0, f"  throat R={b0:g}", color="k", fontsize=8)
    ax.text2D(0.02, 0.02,
              f"(lumps orbit far out at R={R0:g},\n off the edge of this zoom)",
              transform=ax.transAxes, fontsize=8, color="crimson")
    ax.set_title("A.  Wormhole embedding (3-D), zoomed on the neck\n"
                 rf"$Z(R)=\pm b_0\,\mathrm{{arccosh}}(R/b_0)$, $b_0={b0:g}$",
                 fontsize=10)
    ax.set_xlabel("x [code units]")
    ax.set_ylabel("y [code units]")
    ax.set_zlabel("Z (embedding)")
    ax.set_xlim(-R_zoom, R_zoom)
    ax.set_ylim(-R_zoom, R_zoom)
    ax.set_zlim(-R_zoom, R_zoom)
    try:
        ax.set_box_aspect((1, 1, 1))  # true (equal) aspect -> real neck shape
    except Exception:
        pass
    ax.view_init(elev=18, azim=-60)


def build_panel_cross(ax, b0, R0, amp, prof, r_half, r_tenth):
    """Side-on embedding cross-section: the textbook wormhole profile, with the
    exotic-matter radial footprint overlaid at TRUE scale."""
    R_max = R0 + r_tenth + 1.0
    Rg = np.linspace(b0, R_max, 800)
    Zup = embedding_Z(Rg, b0)
    # Both radial arms (R>0 and its mirror) and both sheets.
    for s_r in (+1.0, -1.0):
        ax.plot(s_r * Rg, Zup, color="0.2", lw=2.0)
        ax.plot(s_r * Rg, -Zup, color="0.2", lw=2.0)
    # Throat neck marker.
    ax.plot([b0, b0, -b0, -b0], [0, 0, 0, 0], "o", color="k", ms=4)
    ax.annotate(rf"throat  $R=b_0={b0:g}$", xy=(b0, 0), xytext=(b0 + 1.2, 1.2),
                arrowprops=dict(arrowstyle="->", color="k"), fontsize=9)

    # Exotic-matter radial footprint: a lump centred at R0 with half-max radius
    # r_half occupies distances-from-origin in [R0 - r_half, R0 + r_half].  Shade
    # that band on the upper sheet to show the matter blankets the neck.
    inner = max(R0 - r_half, 0.0)
    outer = R0 + r_half
    for s_r in (+1.0, -1.0):
        rr = np.linspace(max(inner, b0), outer, 200)
        zz = embedding_Z(rr, b0)
        ax.fill_between(s_r * rr, zz - 0.15, zz + 0.15, color="crimson",
                        alpha=0.30, lw=0)
        # lump centre
        zc = float(embedding_Z(np.array([R0]), b0)[0])
        ax.plot(s_r * R0, zc, "s", color="crimson", ms=7)
    ax.plot([], [], color="crimson", alpha=0.4, lw=6,
            label=rf"exotic Q-ball (half-max reaches $R\approx{inner:.1f}$—{outer:.1f})")
    ax.axvspan(np.nan, np.nan)  # no-op to keep legend tidy

    ax.axhline(0, color="0.7", lw=0.6, ls=":")
    ax.set_aspect("equal")
    ax.set_title("B.  Wormhole cross-section (side view) + exotic-matter reach\n"
                 "true scale: a tiny neck blanketed by fat clouds", fontsize=10)
    ax.set_xlabel(r"areal radius $R$  [code units]  ($\pm$ = two radial arms)")
    ax.set_ylabel("Z (embedding)")
    ax.legend(loc="lower center", fontsize=8, framealpha=0.5)


def build_panel_B(ax, b0, R0, omega_orb, amp, prof, r_half, r_tenth):
    # Equatorial |Phi| from the two overlapping Q-ball clouds (as painted).
    span = max(2.4 * R0, 15.0)
    n = 401
    xs = np.linspace(-span, span, n)
    ys = np.linspace(-span, span, n)
    X, Y = np.meshgrid(xs, ys)
    phi = np.zeros_like(X)
    for cx in (R0, -R0):
        rr = np.hypot(X - cx, Y - 0.0)
        phi += amp * envelope_fn(rr, prof)

    pcm = ax.pcolormesh(X, Y, phi, shading="auto", cmap="magma")
    cb = plt.colorbar(pcm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(r"$|\Phi|$  (exotic scalar modulus, code units)")

    # Throat (to scale: tiny!).
    tt = np.linspace(0, 2 * np.pi, 200)
    ax.plot(b0 * np.cos(tt), b0 * np.sin(tt), color="cyan", lw=1.5)
    ax.text(0, -1.4, "throat\n$R=%.2g$" % b0, color="cyan", ha="center",
            fontsize=8)

    # Orbit circle + lump centres + half-max footprints.
    ax.plot(R0 * np.cos(tt), R0 * np.sin(tt), color="lime", ls="--", lw=1.2,
            label=rf"orbit  $R_0={R0:g}$")
    for cx, col in ((R0, "white"), (-R0, "white")):
        ax.scatter([cx], [0], color=col, s=30, zorder=6)
        ax.add_patch(plt.Circle((cx, 0), r_half, fill=False, ec="white",
                                ls=":", lw=1.0, alpha=0.8))

    # Tangential velocity arrows (v = omega_orb * R0), rotation sense CCW.
    v = omega_orb * R0
    for cx in (R0, -R0):
        sgn = 1.0 if cx > 0 else -1.0  # +x lump -> +y, -x lump -> -y
        ax.annotate("", xy=(cx, sgn * 3.0), xytext=(cx, 0),
                    arrowprops=dict(arrowstyle="-|>", color="yellow", lw=1.6))
    ax.text(R0, 3.4, rf"$v=\omega_{{orb}}R_0={v:g}$", color="yellow",
            ha="center", fontsize=8)

    ax.set_aspect("equal")
    ax.set_title("B.  Equatorial (z=0) exotic-matter map — true proportions\n"
                 rf"2 Q-ball lumps, half-max $R\approx{r_half:.1f}$, "
                 rf"10% $R\approx{r_tenth:.1f}$; $\omega_{{orb}}={omega_orb:g}$",
                 fontsize=10)
    ax.set_xlabel("x  [code units  (G=c=1)]")
    ax.set_ylabel("y  [code units  (G=c=1)]")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.4)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--id-dir", default=DEFAULT_ID_DIR,
                    help="ID directory containing qball_profile.dat")
    ap.add_argument("--out", default=None, help="output PNG path")
    ap.add_argument("--b0", type=float, default=0.5,
                    help="throat areal radius (wormhole_throat_radius)")
    ap.add_argument("--orbit-radius", type=float, default=6.0)
    ap.add_argument("--orbit-omega", type=float, default=0.1)
    ap.add_argument("--amp", type=float, default=0.06,
                    help="per-lump central amplitude (kappa*base)")
    args = ap.parse_args()

    id_dir = Path(args.id_dir)
    prof = load_qball_profile(id_dir)
    if prof is None:
        print(f"[warn] no qball_profile.dat in {id_dir}; using Gaussian fallback")
        r_tab = np.linspace(0, 20, 400)
        env_tab = np.exp(-0.5 * (r_tab / 5.7) ** 2)
    else:
        r_tab, env_tab = prof
        print(f"[ok] loaded Q-ball profile ({len(r_tab)} pts) from {id_dir.name}")

    r_half = profile_radius(r_tab, env_tab, 0.5)
    r_tenth = profile_radius(r_tab, env_tab, 0.1)
    print(f"    half-max R={r_half:.2f}, 10% R={r_tenth:.2f} (code units)")

    fig = plt.figure(figsize=(20, 6.6))
    axA = fig.add_subplot(1, 3, 1, projection="3d")
    axC = fig.add_subplot(1, 3, 2)
    axB = fig.add_subplot(1, 3, 3)
    build_panel_A(axA, args.b0, args.orbit_radius, args.orbit_omega, args.amp,
                  prof, r_half)
    build_panel_cross(axC, args.b0, args.orbit_radius, args.amp, prof,
                      r_half, r_tenth)
    build_panel_B(axB, args.b0, args.orbit_radius, args.orbit_omega, args.amp,
                  prof, r_half, r_tenth)

    fig.suptitle(
        "Spinning supported wormhole (config B: kappa=0.6, 2 lumps, "
        rf"$R_0={args.orbit_radius:g}$, $\omega_{{orb}}={args.orbit_omega:g}$)  "
        "— to scale, code units (G=c=1)",
        fontsize=12, y=0.99,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent / "figures" / "wormhole_embedding_configB.png"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"[saved] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
