#!/usr/bin/env python3
"""Generate the figures for the Bondi-dipole article (PRD style).

Everything is derived from the packed results in results/bondi-dipole-runaway/.
Text is rendered by pdflatex through matplotlib's pgf backend, so figure fonts
are the same Computer Modern as the article body.

Outputs (research/bondi_dipole/figures/):
  fig_chase_frames.pdf   chi-1 snapshots with measured barycentres overlaid
                         (finest uniform-grid campaign cell, dx = 0.25)
  fig_trajectories.pdf   worldlines / drifts+controls / gap vs point mass
  fig_controls.pdf       PM vs its MP mirror vs the still PP and MM pairs
  fig_mirror.pdf         mirror-control overlay
  fig_family.pdf         dressed-star family M(omega)
  fig_velocity.pdf       barycentre velocities vs the point-mass model
  fig_convergence.pdf    resolution ladder: the drift converges, the artefact dies
  fig_weyl.pdf           l=2 Weyl amplitude: growth, controls, velocity lockstep
  fig_waveconv.pdf       r*psi4 on four wave-zone shells at retarded time (boxC)
  fig_constraints.pdf    evolution-time constraint norms (appendix)
"""

import csv
import os
import subprocess

import numpy as np
from PIL import Image

import matplotlib

matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PACK = os.path.join(REPO, "results", "bondi-dipole-runaway")
OUT = os.path.join(HERE, "figures")
SCRATCH = os.environ.get("SCRATCH", "/tmp")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "pgf.rcfonts": False,
        "font.family": "serif",
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 8.5,
        "legend.fontsize": 7.3,
        "xtick.labelsize": 7.8,
        "ytick.labelsize": 7.8,
        "axes.linewidth": 0.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.6,
        "ytick.minor.size": 1.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "lines.linewidth": 1.1,
        "legend.frameon": False,
        "legend.handlelength": 1.7,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)

C_CANON = "#0072B2"  # canonical sector (Okabe-Ito blue)
C_PHANT = "#D55E00"  # phantom sector (Okabe-Ito vermillion)
C_CTRL = "0.62"      # null controls
C_GAP = {8: "#0072B2", 12: "#009E73", 16: "#E69F00"}

SINGLE = 3.375  # PRD column width, inches
DOUBLE = 7.0


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


TRAJ = read_csv(os.path.join(PACK, "analysis", "trajectories.csv"))


def series(rows, cell, ykey, xkey="t"):
    xs, ys = [], []
    for r in rows:
        if r["cell"] == cell and r.get(ykey, "") not in ("", None):
            xs.append(float(r[xkey]))
            ys.append(float(r[ykey]))
    return np.array(xs), np.array(ys)


# -------------------------------------------------- campaign-cell helpers
# The quantitative anchor of the article is the uniform-grid convergence
# campaign, campaign/convA_* (dx = 0.50/0.33/0.25) and campaign/boxC_*.
# Alongside them campaign/ also holds the earlier production runs at
# dx = 0.5, named without a prefix (pair_pm, pair_mp_mirror, single_p,
# single_m, and the *_v2 momentum cells).  Those are used only where the
# convergence campaign has no twin -- the mirrored pair, the single stars,
# and the sector-momentum streams -- and every figure that touches them
# says "dx = 0.5" on its face.

M_CANON, M_PHANT = 0.06395, -0.07696  # dressed ADM masses at omega = 0.550


def point_mass(sep, m_phant=M_PHANT, t_end=60.0, dt=0.01):
    """RK4 point-mass Bondi pair (same model as analysis/separation_scaling.py).

    Returns t, drift_canon, drift_phantom, gap."""
    s = np.array([sep / 2, 0.0, -sep / 2, 0.0])

    def deriv(s):
        x1, v1, x2, v2 = s
        r = x2 - x1
        r3 = max(abs(r), 1e-6) ** 3
        return np.array([v1, m_phant * r / r3, v2, M_CANON * (-r) / r3])

    ts, traj = [0.0], [s.copy()]
    for i in range(int(t_end / dt)):
        k1 = deriv(s); k2 = deriv(s + dt / 2 * k1)
        k3 = deriv(s + dt / 2 * k2); k4 = deriv(s + dt * k3)
        s = s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ts.append((i + 1) * dt); traj.append(s.copy())
    ts, traj = np.array(ts), np.array(traj)
    return ts, traj[:, 0] - sep / 2, traj[:, 2] + sep / 2, traj[:, 0] - traj[:, 2]


# ------------------------------------------------------------------ frame crops
def crop_axes_area(png_path):
    """Cut the plot interior out of a packed movie frame.

    The frames are matplotlib output with full-length spines; the main axes is
    the first pair of long vertical dark lines (the colorbar's pair comes
    further right).  Returns the interior pixels, which span the full
    simulation domain [0, 64]^2.
    """
    im = np.asarray(Image.open(png_path).convert("RGB"))
    gray = im.mean(axis=2)
    dark = gray < 128
    h, w = dark.shape

    def line_groups(counts, thresh):
        idx = np.where(counts > thresh)[0]
        groups, start = [], None
        for a, b in zip(idx, np.append(idx[1:], -99)):
            if start is None:
                start = a
            if b != a + 1:
                groups.append((start, a))
                start = None
        return groups

    vgroups = line_groups(dark.sum(axis=0), 0.55 * h)  # vertical lines
    hgroups = line_groups(dark.sum(axis=1), 0.45 * w)  # horizontal lines
    if not vgroups or not hgroups:
        raise RuntimeError(f"spine detection failed for {png_path}")

    def span(groups):
        # light interior: spines show up as two thin lines; dark interior
        # (e.g. viridis): the whole axes area is one wide block
        if groups[0][1] - groups[0][0] > 50:
            return groups[0][0] + 2, groups[0][1] - 1
        if len(groups) < 2:
            raise RuntimeError(f"spine detection failed for {png_path}")
        return groups[0][1] + 1, groups[1][0]

    left, right = span(vgroups)
    top, bottom = span(hgroups)
    return im[top:bottom, left:right]


def extract_frame(movie, index, dest):
    subprocess.run(
        ["ffmpeg", "-loglevel", "error", "-y", "-i", movie,
         "-vf", f"select=eq(n\\,{index})", "-vsync", "vfr",
         "-frames:v", "1", dest],
        check=True,
    )


# ---------------------------------------------------------- fig: chase snapshots
def read_campaign_barycenters(cell):
    """t, bary_x_canon, bary_x_phantom from a campaign cell's stream."""
    path = os.path.join(PACK, "campaign", cell, "sector_barycenters.dat")
    d = np.loadtxt(path)
    return d[:, 0], d[:, 2], d[:, 7]


def fig_chase_frames():
    # The reference cell: separation 12 on the finest uniform grid
    # (dx = 0.25, twice the production resolution).  Its movie stores one
    # frame every 40 steps = 0.2 time units, so index = 5 t.
    cell = "convA_pm_sep12_n256"
    movies = {
        "act": os.path.join(PACK, "movies", cell,
                            "scalar_activity_proj_z.mp4"),
        "chi": os.path.join(PACK, "movies", cell, "chi_minus_1_z.mp4"),
    }
    frames = [(0, 0.0, "$t=0$"), (200, 40.0, "$t=40$"),
              (250, 50.0, "$t=50$"), (300, 60.0, "$t=60$")]

    tb, bxc, bxp = read_campaign_barycenters(cell)
    fig, axes = plt.subplots(2, 4, figsize=(DOUBLE, 3.95),
                             sharex=True, sharey=True)
    for col, (idx, t, title) in enumerate(frames):
        xc = float(np.interp(t, tb, bxc)) - 32.0
        xp = float(np.interp(t, tb, bxp)) - 32.0

        for row, view in enumerate(("act", "chi")):
            ax = axes[row, col]
            raw = os.path.join(SCRATCH, f"_chase_{view}_{idx}.png")
            extract_frame(movies[view], idx, raw)
            ax.imshow(crop_axes_area(raw), extent=(-32, 32, -32, 32),
                      origin="upper", interpolation="bilinear",
                      rasterized=True)
            ink = "white" if view == "act" else "0.25"
            for x, mk in ((xc, "o"), (xp, "s")):
                ax.plot(x, 0, mk, ms=4.0, markerfacecolor="white",
                        markeredgecolor="k", markeredgewidth=0.7)
                ax.plot([x, x], [-3.2, -13], color=ink, lw=0.5,
                        ls=(0, (1, 2)))
            if abs(xc - xp) > 3:
                ax.annotate("", xy=(xp, -13), xytext=(xc, -13),
                            arrowprops=dict(arrowstyle="<->",
                                            mutation_scale=5, lw=0.6,
                                            color=ink, shrinkA=0, shrinkB=0))
            ax.text((xc + xp) / 2, -15.5, f"${abs(xc - xp):.2f}$",
                    ha="center", va="top", fontsize=7.3, color=ink)
            # displacement from the release positions x = +-6: dashed
            # origin ticks and a thin arrow per sector, so the common +x
            # motion is visible panel by panel
            for x0, x1, yy in ((6.0, xc, 11.5), (-6.0, xp, 17.0)):
                ax.plot([x0, x0], [8.5, 19.5], color=ink, lw=0.45,
                        ls=(0, (2, 2)), alpha=0.9)
                if x1 - x0 > 0.4:
                    ax.annotate("", xy=(x1, yy), xytext=(x0, yy),
                                arrowprops=dict(arrowstyle="-|>",
                                                mutation_scale=5, lw=0.7,
                                                color=ink, shrinkA=0,
                                                shrinkB=0))
            ax.set_xlim(-26, 26)
            ax.set_ylim(-26, 26)

        axes[0, col].set_title(title)
        axes[1, col].set_xticks([-20, -10, 0, 10, 20])
        axes[1, col].set_xlabel("$x - x_c$")

    for row in (0, 1):
        axes[row, 0].set_yticks([-20, -10, 0, 10, 20])
        axes[row, 0].set_ylabel("$y - y_c$")

    # row tags and sector identification, first column only
    axes[0, 0].text(0.035, 0.035, r"$\mathcal{A}$",
                    transform=axes[0, 0].transAxes,
                    ha="left", va="bottom", fontsize=8.5, color="white")
    axes[1, 0].text(0.035, 0.035, r"$\chi-1$", transform=axes[1, 0].transAxes,
                    ha="left", va="bottom", fontsize=8.5, color="0.2")
    axes[0, 0].text(-11.5, 7.5, r"$\Phi_-$", ha="center", fontsize=7.8,
                    color="white")
    axes[0, 0].text(11.5, 7.5, r"$\Phi_+$", ha="center", fontsize=7.8,
                    color="white")
    fig.subplots_adjust(wspace=0.08, hspace=0.10)
    fig.savefig(os.path.join(OUT, "fig_chase_frames.pdf"), dpi=300)
    plt.close(fig)


# ------------------------------------------------------------ fig: trajectories
def fig_trajectories():
    """All curves from the finest uniform-grid campaign cells (dx = 0.25);
    the reference cell is the separation-12 pair."""
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(DOUBLE, 2.35))

    # (a) worldlines of the two barycentres in the reference cell: the chase
    t, xc, xp = read_campaign_barycenters("convA_pm_sep12_n256")
    xc, xp = xc - 32.0, xp - 32.0
    ax.fill_between(t, xp, xc, color="0.93", lw=0)
    ax.plot(t, xc, color="k", ls="-", lw=1.0)
    ax.plot(t, xp, color="k", ls="--", lw=1.0)
    ax.text(3.5, 7.6, r"canonical, $M_+>0$", color="k", fontsize=7.3)
    ax.text(3.5, -8.4, r"phantom, $M_-<0$", color="k", fontsize=7.3)
    for tt, side in ((5.0, "r"), (30.0, "r"), (57.0, "l")):
        a = float(np.interp(tt, t, xc))
        b = float(np.interp(tt, t, xp))
        ax.annotate("", xy=(tt, b), xytext=(tt, a),
                    arrowprops=dict(arrowstyle="<->", lw=0.6, color="0.35",
                                    shrinkA=0, shrinkB=0))
        if side == "r":
            ax.annotate(f"${a - b:.1f}$", xy=(tt + 1.6, (a + b) / 2),
                        fontsize=7.0, color="0.35", va="center", ha="left")
        else:
            ax.annotate(f"${a - b:.1f}$", xy=(tt - 1.6, (a + b) / 2),
                        fontsize=7.0, color="0.35", va="center", ha="right")
    ax.set_xlim(0, 66)
    ax.set_ylim(-9.5, 9.5)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$x - x_c$")
    ax.set_title("(a)", loc="left")

    # (b) drifts: reference d0=12 in black, the contact-inflated d0=8 cell in
    # dash-dotted dark gray for contrast, controls dotted.  Runs are
    # distinguished by LINE STYLE, not by gray level alone.
    for cell, col_, color, ls, lw in [
        ("convA_pm_sep12_n256", "c", "k", "-", 1.1),
        ("convA_pm_sep12_n256", "p", "k", "--", 1.1),
        ("convA_pm_n256", "c", "0.4", "-.", 0.85),
        ("convA_pm_n256", "p", "0.4",
         (0, (4, 1.2, 1, 1.2, 1, 1.2)), 0.85),
    ]:
        tt, c_, p_ = read_campaign_barycenters(cell)
        d = (c_ - c_[0]) if col_ == "c" else (p_ - p_[0])
        bx.plot(tt, d, color=color, ls=ls, lw=lw)
    tt, c_, _ = read_campaign_barycenters("convA_pp_n256")
    bx.plot(tt, c_ - c_[0], color="0.45", ls=(0, (1, 1.5)), lw=0.7)
    tt, _, p_ = read_campaign_barycenters("convA_mm_n256")
    bx.plot(tt, p_ - p_[0], color="0.45", ls=(0, (1, 1.5)), lw=0.7)

    # curve identification at the right margin, clear of the data
    for cell, col_, color, dy, lab in [
        ("convA_pm_n256", "p", "0.35", 0.0, r"$\Phi_-$"),
        ("convA_pm_n256", "c", "0.35", +0.3, r"$\Phi_+$"),
        ("convA_pm_sep12_n256", "p", "k", -0.3, r"$\Phi_-$"),
        ("convA_pm_sep12_n256", "c", "k", 0.0, r"$\Phi_+$"),
    ]:
        tt, c_, p_ = read_campaign_barycenters(cell)
        d = (c_ - c_[0]) if col_ == "c" else (p_ - p_[0])
        # x in axes fraction, y in data units: the labels sit clear of both
        # the curve ends and the right spine whatever the panel width
        bx.text(0.79, float(d[-1]) + dy, lab, color=color, fontsize=7.3,
                va="center", ha="left", transform=bx.get_yaxis_transform())

    handles = [
        Line2D([], [], color="k", ls="-", lw=1.1, label=r"$d_0=12$"),
        Line2D([], [], color="0.4", ls="-.", lw=0.85, label=r"$d_0=8$"),
        Line2D([], [], color="0.45", ls=(0, (1, 1.5)), lw=0.7,
               label="controls"),
    ]
    bx.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.06, 0.98),
              borderaxespad=0.0, fontsize=6.6, handlelength=1.9,
              labelspacing=0.35, borderpad=0.2)
    bx.set_xlim(0, 80)
    bx.set_ylim(-0.8, 10.6)
    bx.set_xticks([0, 20, 40, 60])
    bx.set_xlabel("$t$")
    bx.set_ylabel(r"$\Delta x$")
    bx.set_title("(b)", loc="left")

    # (c) gap: NR closes it, the point-mass model grows it slightly.
    for sep, cell, color, ls, lspm, label in [
        (12, "convA_pm_sep12_n256", "0.15", "-", "--", "$d_0=12$"),
        (8, "convA_pm_n256", "0.4", "-.", (0, (1, 1.5)), "$d_0=8$"),
    ]:
        tt, c_, p_ = read_campaign_barycenters(cell)
        cx.plot(tt, c_ - p_, color=color, ls=ls, lw=1.0,
                label=f"NR, {label}")
        tm, _, _, gm = point_mass(sep)
        cx.plot(tm, gm, color=color, ls=lspm, lw=0.9,
                label=f"model, {label}")
    cx.legend(loc="lower left", borderaxespad=0.4, fontsize=6.6,
              handlelength=1.6, labelspacing=0.35, borderpad=0.2)
    cx.set_xlim(0, 62)
    cx.set_ylim(0, 13.4)
    cx.set_xlabel("$t$")
    cx.set_ylabel("$d$")
    cx.set_title("(c)", loc="left")

    fig.subplots_adjust(wspace=0.38)
    fig.savefig(os.path.join(OUT, "fig_trajectories.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- fig: controls
def fig_controls():
    """The four sign assignments side by side.

    Only the mixed pair moves: PM runs towards $+x$, its mirror image MP
    runs the same distance towards $-x$, and the same-sign pairs PP and MM sit
    still.  (The magnified residual ladder lives in fig_convergence(b).)"""
    fig, ax = plt.subplots(figsize=(SINGLE, 2.55))

    # pair barycentre of each configuration, all at d0 = 8
    t, xc, xp = read_campaign_barycenters("convA_pm_n256")
    ax.plot(t, 0.5 * (xc + xp) - 32.0, color="k", ls="-", lw=1.1,
            label=r"\texttt{PM} $(+,-)$")
    tp, dc = series(TRAJ, "pair_pm", "drift_canon")
    _, dp = series(TRAJ, "pair_pm", "drift_phantom")
    ax.plot(tp, 0.5 * (dc + dp), color="0.35", ls="-.", lw=0.85,
            label=r"\texttt{PM}, $\Delta x=0.5$")
    tm, dcm = series(TRAJ, "pair_mp_mirror", "drift_canon")
    _, dpm = series(TRAJ, "pair_mp_mirror", "drift_phantom")
    ax.plot(tm[::2], 0.5 * (dcm + dpm)[::2], ls="none", marker="o", ms=3.0,
            markerfacecolor="white", markeredgewidth=0.8, markeredgecolor="k",
            label=r"\texttt{MP} mirror, $\Delta x=0.5$")
    t, xc, _ = read_campaign_barycenters("convA_pp_n256")
    ax.plot(t, xc - xc[0], color="0.45", ls=(0, (1, 1.5)), lw=0.8,
            label=r"\texttt{PP} $(+,+)$")
    t, _, xp = read_campaign_barycenters("convA_mm_n256")
    ax.plot(t, xp - xp[0], color="0.45", ls=(0, (3, 1.2)), lw=0.8,
            label=r"\texttt{MM} $(-,-)$")
    ax.axhline(0.0, color="0.8", lw=0.5, zorder=0)
    ax.legend(loc="upper left", borderaxespad=0.4, fontsize=6.6,
              handlelength=1.9, labelspacing=0.35, borderpad=0.2)
    ax.set_xlim(0, 62)
    ax.set_ylim(-8.5, 8.5)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$\Delta \bar{x}$")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_controls.pdf"))
    plt.close(fig)

    # numbers quoted in the article
    print("fig_controls stats (drift at t=60):")
    t, xc, xp = read_campaign_barycenters("convA_pm_n256")
    pm60 = float(np.interp(60.0, t, 0.5 * (xc + xp))) - 32.0
    mir60 = 0.5 * (float(np.interp(60.0, tm, dcm)) + float(np.interp(60.0, tm, dpm)))
    amr60 = 0.5 * (float(np.interp(60.0, tp, dc)) + float(np.interp(60.0, tp, dp)))
    print(f"  PM (dx=0.25) {pm60:+.4f} | PM (dx=0.5) {amr60:+.4f} | "
          f"MP mirror {mir60:+.4f} | mirror asymmetry "
          f"{abs(amr60 + mir60) / abs(amr60) * 100:.1f}%")
    for pre, kind, tag in (("convA_pp_n", "c", "PP"), ("convA_mm_n", "p", "MM")):
        row = []
        for n in (128, 192, 256):
            t, xc, xp = read_campaign_barycenters(f"{pre}{n}")
            y = xc if kind == "c" else xp
            row.append(float(np.interp(60.0, t, y)) - 32.0)
        print(f"  {tag}: {['%+.4f' % v for v in row]} "
              f"(finest = {abs(pm60 / row[-1]):.0f}x smaller than \\texttt{{PM}})")


# ------------------------------------------------------------------ fig: mirror
def fig_mirror():
    """Sector-resolved version of Fig. controls(a): run MP is run PM with the
    two signs exchanged in space, so reflecting its drifts about x = 0 must
    reproduce PM star by star if the motion follows the sign pattern rather
    than the grid."""
    fig, ax = plt.subplots(figsize=(SINGLE, 2.55))
    for key, ls, label in [("drift_canon", "-", r"$\Phi_+$, run \texttt{PM}"),
                           ("drift_phantom", "--",
                            r"$\Phi_-$, run \texttt{PM}")]:
        t, d = series(TRAJ, "pair_pm", key)
        ax.plot(t, d, color="k", ls=ls, lw=1.0, label=label)
    for key, marker, label in [
        ("drift_canon", "o", r"$\Phi_+$, run \texttt{MP} reflected"),
        ("drift_phantom", "s", r"$\Phi_-$, run \texttt{MP} reflected"),
    ]:
        t, d = series(TRAJ, "pair_mp_mirror", key)
        ax.plot(t[::2], -d[::2], ls="none", marker=marker, ms=3.0,
                markerfacecolor="white", markeredgewidth=0.8,
                markeredgecolor="k", label=label)
    ax.legend(loc="upper left", borderaxespad=0.4, fontsize=6.8,
              labelspacing=0.35, borderpad=0.2, handlelength=1.9)
    ax.text(60.0, -0.95, r"swapping the signs in space $\Rightarrow$ "
            "motion reversed", fontsize=6.8, color="0.3", ha="right",
            va="bottom")
    ax.set_xlim(0, 62)
    ax.set_ylim(-1.15, 10)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$\Delta x$")
    fig.savefig(os.path.join(OUT, "fig_mirror.pdf"))
    plt.close(fig)


# ------------------------------------------------------------------ fig: family
def fig_family():
    fam = {"canonical": {}, "phantom": {}}
    with open(os.path.join(PACK, "stars", "star_family.csv")) as f:
        for r in csv.DictReader(f):
            if r["status"] == "ok":
                fam[r["sector"]][float(r["omega_requested"])] = float(r["adm_mass"])

    fig, ax = plt.subplots(figsize=(SINGLE, 2.7))
    for sector, ls, marker, label in [
        ("canonical", "-", "o", "canonical"),
        ("phantom", "--", "s", "phantom"),
    ]:
        om = np.array(sorted(fam[sector]))
        m = np.array([fam[sector][o] for o in om])
        ax.plot(om, m, color="k", ls=ls, marker=marker, ms=2.9, lw=0.9,
                markerfacecolor="white", markeredgewidth=0.8, label=label)
    # branch identification, instead of a legend box in the crowded corner
    ax.text(0.755, 0.0245, "canonical", fontsize=7.3, color="0.1",
            ha="center", va="bottom")
    ax.text(0.755, -0.0245, "phantom", fontsize=7.3, color="0.1",
            ha="center", va="top")

    ax.axhline(0, color="0.5", lw=0.5)
    ax.axvspan(0.500, 0.5265, color="0.90", zorder=0, lw=0)
    ax.text(0.5133, -0.012, "no dressed star", rotation=90, ha="center",
            va="center", fontsize=7.0, color="0.4", zorder=4,
            bbox=dict(facecolor="0.90", edgecolor="none", pad=1.2))

    # the three stars actually used: the reference pair takes one frequency
    # for both sectors, the equal-|M| variant retunes the phantom only
    ax.plot([0.550, 0.550], [0.063951, -0.076960], color="0.35", lw=0.6,
            ls=(0, (2, 2)), zorder=4)
    ax.plot([0.550, 0.550], [0.063951, -0.076960], marker="*", ms=9,
            color="k", ls="none", zorder=5)
    ax.plot([0.56598], [-0.063950], marker="*", ms=9, color="k",
            markerfacecolor="white", markeredgewidth=0.8, ls="none", zorder=5)
    ax.plot([0.550, 0.56598], [-0.063950, -0.063950], color="0.35", lw=0.6,
            ls=(0, (1, 1.5)), zorder=4)

    ax.annotate("reference pair: both sectors\n"
                r"at $\omega=0.550$, $M_\pm=+0.064/{-0.077}$",
                xy=(0.5545, 0.0625), xytext=(0.598, 0.0795), fontsize=7.0,
                color="0.2", ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=2, shrinkB=2))
    ax.annotate(r"equal-$|M|$ variant: $\Phi_-$ retuned" "\n"
                r"to $\omega=0.566$, $M_-=-0.064$",
                xy=(0.5695, -0.0645), xytext=(0.613, -0.0925), fontsize=7.0,
                color="0.2", ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=2, shrinkB=2))
    ax.set_xlim(0.50, 0.86)
    ax.set_ylim(-0.108, 0.092)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$M_{\rm ADM}$")
    fig.savefig(os.path.join(OUT, "fig_family.pdf"))
    plt.close(fig)


# ------------------------------------------------------------ fig: schematic
def fig_schematic():
    """Bondi's sign rules: the three pairings of active mass.

    Grayscale like the rest of the article -- the sector is carried by the
    outline style and the label, never by shade alone.  Each row states what
    happens to the pair's barycentre, which is the quantity the runs measure.
    """
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots(figsize=(SINGLE, 2.62))
    ax.set_aspect("equal")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.75)
    ax.axis("off")

    R = 0.50
    XL, XR = 2.15, 4.85
    Y = {"a": 5.55, "b": 3.55, "c": 1.35}

    ax.text(5.0, 7.45, r"each body is accelerated by the \emph{other's} "
            r"active mass;", ha="center", va="center", fontsize=6.6,
            color="0.15")
    ax.text(5.0, 6.95, r"inertia is $|M|$ and positive in every case",
            ha="center", va="center", fontsize=6.6, color="0.15")

    def body(x, y, sector):
        if sector == "+":
            fc, ec, lab, ls = "0.88", "k", "$+M$", "-"
        else:
            fc, ec, lab, ls = "white", "0.25", "$-M$", (0, (2.6, 1.4))
        ax.add_patch(mpatches.Circle((x, y), R, facecolor=fc, edgecolor=ec,
                                     lw=0.9, linestyle=ls, zorder=3))
        ax.text(x, y, lab, ha="center", va="center", fontsize=6.8, color="k",
                zorder=4)

    def force(x0, x1, y):
        ax.annotate("", xy=(x1, y), xytext=(x0, y), zorder=5,
                    arrowprops=dict(arrowstyle="-|>", color="k", lw=0.9,
                                    shrinkA=0, shrinkB=0, mutation_scale=7))

    def verdict(y, text, sub, bold=False):
        ax.text(6.05, y + 0.28, text, ha="left", va="center", fontsize=7.3,
                fontweight="bold" if bold else "normal")
        ax.text(6.05, y - 0.34, sub, ha="left", va="center", fontsize=5.9,
                color="0.35")

    def letter(y, s):
        ax.text(0.18, y + 0.72, s, ha="left", va="center", fontsize=7.5)

    def caption(y, s, x=3.5):
        ax.text(x, y + 0.72, s, ha="center", va="center", fontsize=6.0,
                color="0.3")

    # (a) + + : each attracted by the other
    y = Y["a"]
    letter(y, "(a)")
    body(XL, y, "+"); body(XR, y, "+")
    force(XL + R + 0.12, XL + R + 0.80, y)
    force(XR - R - 0.12, XR - R - 0.80, y)
    caption(y, "both attracted")
    verdict(y, "fall together", "barycentre fixed")

    # (b) - - : each repelled by the other
    y = Y["b"]
    letter(y, "(b)")
    body(XL, y, "-"); body(XR, y, "-")
    force(XL - R - 0.12, XL - R - 0.80, y)
    force(XR + R + 0.12, XR + R + 0.80, y)
    caption(y, "both repelled")
    verdict(y, "push apart", "barycentre fixed")

    # (c) - + : the Bondi dipole; both forces point the same way
    y = Y["c"]
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.02, y - 0.98), 9.96, 2.12, boxstyle="round,pad=0,rounding_size=0.12",
        facecolor="0.95", edgecolor="0.75", lw=0.5, zorder=0))
    letter(y, "(c)")
    body(XL, y, "-"); body(XR, y, "+")
    force(XL + R + 0.12, XL + R + 0.80, y)
    force(XR + R + 0.12, XR + R + 0.80, y)
    caption(y, "attracted", x=XL + R + 0.46)
    caption(y, "repelled", x=XR + R + 0.46)
    verdict(y, "run away", "barycentre accelerates", bold=True)

    # orientation of the drift axis used throughout the article
    ax.annotate("", xy=(XR + R + 0.95, y - 1.28), xytext=(XL - R, y - 1.28),
                arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.6,
                                shrinkA=0, shrinkB=0, mutation_scale=6))
    ax.text(XR + R + 1.10, y - 1.28, "$x$", ha="left", va="center",
            fontsize=7.0, color="0.45")

    fig.savefig(os.path.join(OUT, "fig_schematic.pdf"))
    plt.close(fig)


# ------------------------------------------------------------ fig: constraints
def read_norms(cell, tree="campaign"):
    """t, composite L2 Ham, composite L2 Mom, composite Linf Ham.

    Columns 12/13/16 of constraint_norms.dat: the whole-hierarchy norms
    (covered coarse cells dropped, outermost 4 coarse cells excluded), which
    the diagnostic code marks as the accuracy figures; cols 2/3 are the
    level-0-only governor input.
    """
    path = os.path.join(PACK, tree, cell, "constraint_norms.dat")
    rows = []
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append((float(v[0]), float(v[11]), float(v[12]),
                         float(v[15])))
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def fig_constraints():
    """Constraint record of the reference cell's full resolution ladder: the
    three grids' L2 curves lie on top of one another (the elliptic-solve
    floor), against the campaign controls run under identical numerics."""
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(SINGLE, 3.6), sharex=True)

    ladder = {n: read_norms(f"convA_pm_sep12_n{n}") for n in (128, 192, 256)}
    tp, hamp, momp, _ = read_norms("convA_pp_n256")
    tm, hamm, momm, _ = read_norms("convA_mm_n256")
    t, ham, mom, linf = ladder[256]

    # (a) Hamiltonian: L2 for the three grids of the reference cell (they
    # coincide -- the initial-data floor), Linf for the finest grid, controls
    # in dotted / tightly dashed gray: style, not shade, carries them.
    ax.plot(t, linf, color="k", ls="-.", lw=0.8)
    ax.plot(tp, hamp, color="0.45", ls=(0, (1, 1.5)), lw=0.7)
    ax.plot(tm, hamm, color="0.45", ls=(0, (3, 1.2)), lw=0.7)
    for n, (tn, hamn, _, _) in ladder.items():
        ax.plot(tn[tn >= 2], hamn[tn >= 2], color="k", ls="-",
                lw=1.0 if n == 256 else 0.55)
    handles = [
        Line2D([], [], color="k", ls="-", lw=1.0,
               label=r"$L_2$, all three grids"),
        Line2D([], [], color="k", ls="-.", lw=0.8,
               label=r"$L_\infty$, $\Delta x=0.25$"),
        Line2D([], [], color="0.45", ls=(0, (1, 1.5)), lw=0.7,
               label=r"\texttt{PP} control"),
        Line2D([], [], color="0.45", ls=(0, (3, 1.2)), lw=0.7,
               label=r"\texttt{MM} control"),
    ]
    ax.legend(handles=handles, loc="upper left", borderaxespad=0.4, ncol=2,
              fontsize=6.6, handlelength=1.9, columnspacing=1.0,
              labelspacing=0.3, borderpad=0.2)
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 300)
    ax.set_ylabel(r"$\mathcal{H}$")
    ax.set_title("(a)", loc="left")

    # (b) momentum: L2 only
    bx.plot(tp, momp, color="0.45", ls=(0, (1, 1.5)), lw=0.7)
    bx.plot(tm, momm, color="0.45", ls=(0, (3, 1.2)), lw=0.7)
    bx.plot(t, mom, color="k", ls="-", lw=1.0)
    bx.annotate("controls", xy=(20.0, 1.5e-4), xytext=(9.0, 6.5e-4),
                fontsize=7.0, color="0.35",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=1, shrinkB=2))
    bx.set_yscale("log")
    bx.set_ylim(1.5e-6, 4e-3)
    bx.set_xlim(0, 62)
    bx.set_xlabel("$t$")
    bx.set_ylabel(r"$\mathcal{M}$")
    bx.set_title("(b)", loc="left")

    fig.subplots_adjust(hspace=0.24)
    fig.savefig(os.path.join(OUT, "fig_constraints.pdf"))
    plt.close(fig)


# ------------------------------------------------------------------ fig: velocity
def central_diff(t, x, width=2.0):
    """Centred difference over +-width time units, whatever the cadence."""
    v = np.full_like(x, np.nan)
    half = max(1, int(round(width / float(np.median(np.diff(t))))))
    v[half:-half] = (x[2 * half:] - x[:-2 * half]) / (t[2 * half:] - t[:-2 * half])
    return v


def fig_velocity():
    fig, ax = plt.subplots(figsize=(SINGLE, 2.45))

    t, xc, xp = read_campaign_barycenters("convA_pm_sep12_n256")
    ax.plot(t, central_diff(t, xc), color="k", ls="-", lw=1.1)
    ax.plot(t, central_diff(t, xp), color="k", ls="--", lw=1.1)

    # matched point-mass model, for scale: dotted (canonical) and dash-dotted
    # (phantom) so it survives grayscale print next to the black NR curves
    tm, dc, dp, _ = point_mass(12)
    ax.plot(tm, np.gradient(dc, tm), color="0.4", ls=(0, (1, 1.5)), lw=0.9)
    ax.plot(tm, np.gradient(dp, tm), color="0.4", ls="-.", lw=0.9)

    ax.text(38.5, 0.078, r"$\Phi_-$", fontsize=8.5)
    ax.text(49.5, 0.056, r"$\Phi_+$", fontsize=8.5)
    ax.annotate("point mass", xy=(20.0, 0.0085), xytext=(11.5, 0.043),
                fontsize=7.0, color="0.35",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=1, shrinkB=2))

    ax.set_xlim(0, 62)
    ax.set_ylim(-0.01, 0.17)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$v_x$")
    fig.savefig(os.path.join(OUT, "fig_velocity.pdf"))
    plt.close(fig)


# ---------------------------------------------------------- fig: convergence
def fig_convergence():
    """The uniform-grid resolution ladder.

    (a) the reference (d0 = 12) cell in time: the mixed-pair drift is
    grid-stable (the three black curves lie on top of one another), while the
    two-positive control's residual drift roughly halves at every refinement
    step (the gray fan).  (b) the end-of-window drift against grid spacing for
    every campaign family: the mixed pairs sit on plateaus, the same-sign
    controls fall towards zero."""
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE, 2.7))

    styles = {256: "-", 192: (0, (4, 1.6)), 128: (0, (1, 1.5))}
    for n in (128, 192, 256):
        t, xc, xp = read_campaign_barycenters(f"convA_pm_sep12_n{n}")
        drift = 0.5 * (xc + xp) - 32.0
        ax.plot(t, np.abs(drift), color="k", ls=styles[n],
                lw=1.1 if n == 256 else 0.85, zorder=5)
    for n in (128, 192, 256):
        t, xc, _ = read_campaign_barycenters(f"convA_pp_n{n}")
        ax.plot(t, np.abs(xc - 32.0), color="0.55", ls=styles[n],
                lw=0.95 if n == 256 else 0.75)

    ax.annotate("mixed pair, $d_0=12$:\ngrid-stable", xy=(46.0, 0.72),
                xytext=(40.0, 7.0), fontsize=7.0, color="0.1", ha="center",
                va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=2, shrinkB=2))
    ax.annotate("control: halves per\nrefinement step", xy=(31.0, 0.030),
                xytext=(31.0, 0.0032), fontsize=7.0, color="0.35",
                ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=2, shrinkB=2))

    handles = [
        Line2D([], [], color="k", ls=styles[128], lw=0.85,
               label=r"$\Delta x=0.50$"),
        Line2D([], [], color="k", ls=styles[192], lw=0.85,
               label=r"$\Delta x=0.33$"),
        Line2D([], [], color="k", ls=styles[256], lw=1.1,
               label=r"$\Delta x=0.25$"),
    ]
    ax.legend(handles=handles, loc="upper left", borderaxespad=0.4)
    ax.set_yscale("log")
    ax.set_xlim(0, 62)
    ax.set_ylim(1.5e-3, 30)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$|\Delta \bar{x}|$")
    ax.set_title("(a)", loc="left")

    # (b) the same test read as a refinement sequence: one point per grid,
    # each family labelled at its own curve so no legend is needed
    dx = np.array([0.50, 1.0 / 3.0, 0.25])
    families = [
        ("convA_pm_n", "c", "k", "-", "o", r"$d_0=8$"),
        ("convA_pm_sep12_n", "c", "k", "-", "D", r"$d_0=12$"),
        ("convA_pm_sep16_n", "c", "k", "-", "^", r"$d_0=16$"),
        ("convA_pp_n", "pp", "0.55", (0, (1, 1.5)), "v", r"\texttt{PP}"),
        ("convA_mm_n", "mm", "0.55", (0, (3, 1.2)), "s", r"\texttt{MM}"),
    ]
    for pre, kind, color, ls, mk, label in families:
        vals = []
        for n in (128, 192, 256):
            t, xc, xp = read_campaign_barycenters(f"{pre}{n}")
            y = {"c": 0.5 * (xc + xp), "pp": xc, "mm": xp}[kind]
            vals.append(abs(float(np.interp(60.0, t, y)) - 32.0))
        bx.plot(dx, vals, color=color, ls=ls, lw=0.9, marker=mk, ms=3.4,
                markerfacecolor="white", markeredgewidth=0.7)
        bx.text(0.238, vals[-1], "~" + label, color=color, fontsize=7.0,
                va="center", ha="left")
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.set_xlim(0.56, 0.135)          # finer grids to the right
    bx.set_ylim(2.2e-3, 60)
    bx.set_xticks([0.5, 1 / 3, 0.25])
    bx.set_xticklabels(["$0.50$", "$0.33$", "$0.25$"])
    bx.minorticks_off()
    bx.text(0.50, 22.0, "pairs: same answer\non every grid", fontsize=7.0,
            color="0.1", ha="left", va="center")
    bx.text(0.52, 0.0040, "controls: shrink with the grid", fontsize=7.0,
            color="0.4", ha="left", va="center")
    bx.set_xlabel(r"$\Delta x$ \quad (grid refined $\rightarrow$)")
    bx.set_ylabel(r"$|\Delta \bar{x}(t=60)|$")
    bx.set_title("(b)", loc="left")

    fig.subplots_adjust(wspace=0.26)
    fig.savefig(os.path.join(OUT, "fig_convergence.pdf"))
    plt.close(fig)

    # numbers quoted in the article
    print("fig_convergence stats (pair drift / control drift):")
    for pre, tag in (("convA_pm_sep12_n", "sep12"), ("convA_pm_n", "sep8"),
                     ("convA_pm_eqm_n", "eqm"), ("convA_pm_sep16_n", "sep16")):
        for tt in (15.0, 30.0, 60.0):
            row = []
            for n in (128, 192, 256):
                t, xc, xp = read_campaign_barycenters(f"{pre}{n}")
                row.append(float(np.interp(tt, t, 0.5 * (xc + xp))) - 32.0)
            spread = (max(row) - min(row)) / abs(row[-1]) * 100
            ctrl = []
            for n in (128, 192, 256):
                t, xc, _ = read_campaign_barycenters(f"convA_pp_n{n}")
                ctrl.append(float(np.interp(tt, t, xc)) - 32.0)
            print(f"  {tag} t={tt:.0f}: {['%+.4f' % v for v in row]} "
                  f"spread {spread:.1f}% | pp {['%+.4f' % v for v in ctrl]}"
                  f" | ratio(fine) {abs(row[-1] / ctrl[-1]):.0f}x")


# ---------------------------------------------------------------- fig: weyl
def read_l2_modes(cell, tree="campaign"):
    """t and the five complex (l=2,m) amplitudes r*psi4 at each shell.

    psi4_mode_l2_all.dat columns: time, then (Re,Im) pairs for m=-2..2 at
    R=8, then the same ten columns at R=16.
    """
    path = os.path.join(PACK, tree, cell, "psi4_mode_l2_all.dat")
    d = np.loadtxt(path)
    modes = {}
    for ri, R in enumerate((8, 16)):
        base = 1 + ri * 10
        for mi, m in enumerate((-2, -1, 0, 1, 2)):
            modes[(R, m)] = d[:, base + 2 * mi] + 1j * d[:, base + 2 * mi + 1]
    return d[:, 0], modes


def l2_amplitude(modes, R):
    return np.sqrt(sum(np.abs(modes[(R, m)]) ** 2 for m in (-2, -1, 0, 1, 2)))


def fig_weyl():
    """The l=2 Weyl signal, anchored on the campaign's finest d0 = 8 cells
    (the strongest signal); the mirror overlay comes from the only mirrored
    run, which exists at the production resolution dx = 0.5."""
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(SINGLE, 3.8), sharex=True)
    t_pm, modes_pm = read_l2_modes("convA_pm_n256")

    # (a) total l=2 amplitude at R=16: mixed cells vs the same-sign controls.
    # Style, not shade, separates the runs (as everywhere in the article).
    A_pm = l2_amplitude(modes_pm, 16)
    ax.plot(t_pm, A_pm, color="k", ls="-", lw=1.1, zorder=5)
    t_eq, modes_eq = read_l2_modes("convA_pm_eqm_n256")
    ax.plot(t_eq, l2_amplitude(modes_eq, 16), color="0.35", ls="-.", lw=0.9,
            zorder=4)
    t_mp, modes_mp = read_l2_modes("pair_mp_mirror", tree="campaign")
    A_mp = l2_amplitude(modes_mp, 16)
    ax.plot(t_mp[::5], A_mp[::5], ls="none", marker="o", ms=2.7,
            markerfacecolor="white", markeredgewidth=0.7,
            markeredgecolor="k", zorder=6)
    for cell, ls in (("convA_pp_n256", (0, (1, 1.5))),
                     ("convA_mm_n256", (0, (3, 1.2)))):
        tc, mc = read_l2_modes(cell)
        ax.plot(tc, l2_amplitude(mc, 16), color="0.45", ls=ls, lw=0.7)

    handles = [
        Line2D([], [], color="k", ls="-", lw=1.1, label=r"\texttt{PM}"),
        Line2D([], [], color="0.35", ls="-.", lw=0.9,
               label=r"\texttt{PM-eq}"),
        Line2D([], [], color="k", lw=0, marker="o", ms=2.7,
               markerfacecolor="white", markeredgewidth=0.7,
               label=r"\texttt{MP} mirror"),
        Line2D([], [], color="0.45", ls=(0, (1, 1.5)), lw=0.7,
               label=r"\texttt{PP}"),
        Line2D([], [], color="0.45", ls=(0, (3, 1.2)), lw=0.7,
               label=r"\texttt{MM}"),
    ]
    # one legend for the whole stack, above the frame: inside either panel it
    # would land on the wandering control curves
    ax.legend(handles=handles, loc="lower left",
              bbox_to_anchor=(0.0, 1.01, 1.0, 0.10), mode="expand", ncol=3,
              borderaxespad=0.0, fontsize=6.8, handlelength=1.6,
              columnspacing=1.0, handletextpad=0.5)
    ax.set_yscale("log")
    ax.set_ylim(3e-6, 3e-2)
    ax.set_ylabel(r"$A_{\ell=2}$")
    ax.text(0.025, 0.94, "(a)", transform=ax.transAxes, ha="left", va="top")

    # (b) the d0=8 amplitude against the rescaled pair speed: the near-zone
    # quadrupole curvature rises in lockstep with the runaway.
    tb, bxc, bxp = read_campaign_barycenters("convA_pm_n256")
    v = central_diff(tb, 0.5 * (bxc + bxp))
    Ai = np.interp(tb, t_pm, A_pm)
    ok = np.isfinite(v) & (tb >= 5) & (tb <= 55)
    kappa = float(np.sum(Ai[ok] * np.abs(v[ok])) / np.sum(v[ok] ** 2))
    corr = float(np.corrcoef(Ai[ok], np.abs(v[ok]))[0, 1])

    bx.plot(t_pm, A_pm, color="k", ls="-", lw=1.1,
            label=r"$A_{\ell=2}$, run \texttt{PM}")
    bx.plot(tb, kappa * np.abs(v), color="0.4", ls="--", lw=0.9,
            label=rf"$\kappa\,|\dot X|$, $\kappa={kappa:.3f}$")
    bx.legend(loc="upper left", bbox_to_anchor=(0.075, 0.99),
              borderaxespad=0.0, fontsize=6.8, handlelength=1.7,
              labelspacing=0.35, borderpad=0.2)
    bx.set_xlim(0, 62)
    bx.set_ylim(-0.0002, 0.0068)
    bx.set_xlabel("$t$")
    bx.set_ylabel(r"$A_{\ell=2}$")
    bx.text(0.025, 0.94, "(b)", transform=bx.transAxes, ha="left", va="top")

    fig.subplots_adjust(hspace=0.14)
    fig.savefig(os.path.join(OUT, "fig_weyl.pdf"))
    plt.close(fig)

    # ---- numbers quoted in the article text -------------------------------
    print("fig_weyl stats:")
    print(f"  corr(A16,|v|) 5<=t<=55 (central diff +-2): {corr:.3f}")
    print(f"  kappa: {kappa:.4f}")
    for cell, (tc, mc) in (("convA_pm_n256", (t_pm, modes_pm)),
                           ("convA_pm_eqm_n256", (t_eq, modes_eq)),
                           ("pair_mp_mirror", (t_mp, modes_mp))):
        A = l2_amplitude(mc, 16)
        floor = np.mean(A[(tc >= 0) & (tc <= 10)])
        late = np.mean(A[(tc >= 55) & (tc <= 60)])
        print(f"  growth {cell}: {late / floor:.0f}x "
              f"({floor:.2e} -> {late:.2e})")
    t_ref, modes_ref = read_l2_modes("convA_pm_sep12_n256")
    A_ref = l2_amplitude(modes_ref, 16)
    print("  growth convA_pm_sep12_n256: "
          f"{np.mean(A_ref[(t_ref >= 55) & (t_ref <= 60)]) / np.mean(A_ref[(t_ref >= 0) & (t_ref <= 10)]):.0f}x")
    # the mirror run exists only at dx=0.5: compare it against its own-
    # resolution reference (campaign/pair_pm) so the test is same-numerics
    t_pmamr, modes_pmamr = read_l2_modes("pair_pm", tree="campaign")
    A_pmamr = l2_amplitude(modes_pmamr, 16)
    Ai_mp = np.interp(t_pmamr, t_mp, A_mp)
    rel = np.max(np.abs(A_pmamr - Ai_mp)) / np.max(A_pmamr)
    print(f"  mirror amplitude reproduction (same dx=0.5): max rel {rel:.3f}")
    sel = (t_pm >= 45) & (t_pm <= 60)
    p = {m: np.mean(np.abs(modes_pm[(16, m)][sel]) ** 2)
         for m in (-2, -1, 0, 1, 2)}
    print(f"  m=+-1 power fraction (45-60): {(p[1] + p[-1]) / sum(p.values()):.4f}")
    A8 = l2_amplitude(modes_pm, 8)
    for tt in (30, 60):
        i = int(np.argmin(np.abs(t_pm - tt)))
        print(f"  A(R=8)/A(R=16) at t={tt}: {A8[i] / A_pm[i]:.1f}")
    # spectral content of the amplitude signal, growth phase
    from numpy.fft import rfft, rfftfreq
    selg = (t_pm >= 15) & (t_pm <= 60)
    y = A_pm[selg] - np.mean(A_pm[selg])
    Y = np.abs(rfft(y * np.hanning(len(y)))) ** 2
    f = rfftfreq(len(y), t_pm[1] - t_pm[0]) * 2.0 * np.pi
    print(f"  A16 spectral fraction below omega=0.4: {Y[f < 0.4].sum() / Y.sum():.4f}")
    print(f"  ... in [0.45,0.65] (field omega): {Y[(f > 0.45) & (f < 0.65)].sum() / Y.sum():.2e}")
    print(f"  ... in [1.0,1.2] (2 omega): {Y[(f > 1.0) & (f < 1.2)].sum() / Y.sum():.2e}")


# ------------------------------------------------------------ fig: wave zone
def read_mode_shells(cell):
    """t, {R: complex r*psi4 l=2 m=0} from a campaign cell's mode stream.

    The header labels each (Re, Im) column pair with its shell radius R=...
    """
    import re

    path = os.path.join(PACK, "campaign", cell, "psi4_mode_l2m0.dat")
    with open(path) as f:
        header = f.readline()
    radii = []
    for r in re.findall(r"R=([0-9.]+)", header):
        if float(r) not in radii:
            radii.append(float(r))
    d = np.loadtxt(path)
    modes = {R: R * (d[:, 1 + 2 * i] + 1j * d[:, 2 + 2 * i])
             for i, R in enumerate(radii)}
    return d[:, 0], modes


def fig_waveconv():
    """One outgoing wave: |r psi4^{2,0}| on the four wave-zone shells of the
    double-size box agrees at fixed retarded time inside u = 10-25; before it
    the shells carry start-up junk, after it the expanding matter arrives.
    Circles: the half-size box at R = 16 -- the amplitude is box-independent."""
    fig, ax = plt.subplots(figsize=(SINGLE, 2.55))

    t, modes = read_mode_shells("boxC_pm_L128_n256")
    styles = {16.0: "-", 24.0: (0, (4, 1.6)), 32.0: "-.", 40.0: (0, (1, 1.5))}
    ax.axvspan(10, 25, color="0.93", zorder=0, lw=0)
    for R, ls in styles.items():
        ax.plot(t - R, np.abs(modes[R]), color="k", ls=ls,
                lw=1.1 if R == 16.0 else 0.85, label=f"$r={R:.0f}$")

    ts, small = read_mode_shells("convA_pm_n128")
    u = ts - 16.0
    keep = (u >= 0) & (u <= 30)
    ax.plot(u[keep][::12], np.abs(small[16.0])[keep][::12], ls="none",
            marker="o", ms=2.7, markerfacecolor="white",
            markeredgewidth=0.7, markeredgecolor="k",
            label="$r=16$, $L=64$ box")

    ax.text(17.5, 0.0295, "one outgoing\nwave", ha="center", va="bottom",
            fontsize=7.0, color="0.25")
    ax.annotate("matter reaches\nthe shells", xy=(41.0, 0.0505),
                xytext=(34.0, 0.0625), fontsize=7.0, color="0.35",
                ha="center", va="top",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=2, shrinkB=3))
    # the legend sits above the frame: inside, it would land on the shaded
    # comparison window
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01, 1.0, 0.14),
              mode="expand", ncol=3, borderaxespad=0.0, fontsize=6.8,
              handlelength=1.7, columnspacing=1.0, handletextpad=0.5)
    ax.set_xlim(0, 46)
    ax.set_ylim(0, 0.065)
    ax.set_xlabel("$u = t - r$")
    ax.set_ylabel(r"$|r\,\psi_4^{2,0}|$")
    fig.savefig(os.path.join(OUT, "fig_waveconv.pdf"))
    plt.close(fig)

    # numbers quoted in the article
    print("fig_waveconv stats:")
    for uu in (10.0, 15.0, 20.0, 25.0):
        a = [float(np.interp(uu + R, t, np.abs(modes[R]))) for R in styles]
        print(f"  u={uu:.0f}: spread {max(a) / min(a):.2f}x")
    for uu in (10.0, 15.0, 20.0, 25.0):
        big = float(np.interp(uu + 16.0, t, np.abs(modes[16.0])))
        sml = float(np.interp(uu + 16.0, ts, np.abs(small[16.0])))
        print(f"  u={uu:.0f}: box L=128 vs 64 at R=16: {abs(big - sml) / big * 100:.2f}%")


if __name__ == "__main__":
    fig_chase_frames()
    fig_trajectories()
    fig_controls()
    fig_mirror()
    fig_schematic()
    fig_family()
    fig_velocity()
    fig_convergence()
    fig_weyl()
    fig_waveconv()
    fig_constraints()
    print("wrote figures to", OUT)
