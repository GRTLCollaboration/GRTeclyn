#!/usr/bin/env python3
"""Generate the figures for the Bondi-dipole article (PRD style).

Everything is derived from the packed results in results/bondi-dipole-runaway/.
Text is rendered by pdflatex through matplotlib's pgf backend, so figure fonts
are the same Computer Modern as the article body.

Outputs (research/bondi_dipole/figures/):
  fig_chase_frames.pdf   chi-1 snapshots with measured barycentres overlaid
  fig_trajectories.pdf   worldlines / drifts+controls / gap vs point mass
  fig_mirror.pdf         mirror-control overlay
  fig_family.pdf         dressed-star family M(omega)
  fig_sepscaling.pdf     separation scaling vs matched point-mass integrations
  fig_velocity.pdf       barycentre velocities vs the point-mass model
  fig_weyl.pdf           l=2 Weyl amplitude: growth, controls, velocity lockstep
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
NEWT = read_csv(os.path.join(PACK, "analysis", "newtonian_reference.csv"))
SEPS = read_csv(os.path.join(PACK, "analysis", "separation_scaling.csv"))


def series(rows, cell, ykey, xkey="t"):
    xs, ys = [], []
    for r in rows:
        if r["cell"] == cell and r.get(ykey, "") not in ("", None):
            xs.append(float(r[xkey]))
            ys.append(float(r[ykey]))
    return np.array(xs), np.array(ys)


def at_time(cell, key, t):
    xs, ys = series(TRAJ, cell, key)
    return float(np.interp(t, xs, ys))


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
def fig_chase_frames():
    movies = {
        "act": os.path.join(PACK, "movies", "pair_pm",
                            "scalar_activity_proj_z.mp4"),
        "chi": os.path.join(PACK, "movies", "pair_pm", "chi_minus_1_z.mp4"),
    }
    frames = [(0, 0.0, "$t=0$"), (101, 40.4, "$t=40.4$"),
              (132, 52.8, "$t=52.8$"), (152, 60.0, "$t=60$")]

    fig, axes = plt.subplots(2, 4, figsize=(DOUBLE, 3.95),
                             sharex=True, sharey=True)
    for col, (idx, t, title) in enumerate(frames):
        xc = at_time("pair_pm", "x_canon", t) - 32.0
        xp = at_time("pair_pm", "x_phantom", t) - 32.0

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
            # displacement from the release positions x = +-4: dashed
            # origin ticks and a thin arrow per sector, so the common +x
            # motion is visible panel by panel
            for x0, x1, yy in ((4.0, xc, 11.5), (-4.0, xp, 17.0)):
                ax.plot([x0, x0], [8.5, 19.5], color=ink, lw=0.45,
                        ls=(0, (2, 2)), alpha=0.9)
                if x1 - x0 > 0.5:
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
    axes[0, 0].text(0.04, 0.96, r"$\mathcal{A}$", transform=axes[0, 0].transAxes,
                    ha="left", va="top", fontsize=8.5, color="white")
    axes[1, 0].text(0.04, 0.96, r"$\chi-1$", transform=axes[1, 0].transAxes,
                    ha="left", va="top", fontsize=8.5, color="0.2")
    axes[0, 0].text(-9.5, 6.5, r"$\Phi_-$", ha="center", fontsize=7.8,
                    color="white")
    axes[0, 0].text(9.5, 6.5, r"$\Phi_+$", ha="center", fontsize=7.8,
                    color="white")
    axes[1, 0].annotate(r"$\Phi_-$\,: $\chi>1$", xy=(-7.5, 5.5),
                        xytext=(-24.5, 16.5), fontsize=7.3, ha="left",
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="0.2",
                                        shrinkA=1, shrinkB=1))
    axes[1, 0].annotate(r"$\Phi_+$\,: $\chi<1$", xy=(6.5, -5.5),
                        xytext=(2.5, -21.5), fontsize=7.3, ha="left",
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="0.2",
                                        shrinkA=1, shrinkB=1))
    fig.subplots_adjust(wspace=0.08, hspace=0.10)
    fig.savefig(os.path.join(OUT, "fig_chase_frames.pdf"), dpi=300)
    plt.close(fig)


# ------------------------------------------------------------ fig: trajectories
def fig_trajectories():
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(DOUBLE, 2.35))

    # (a) worldlines of the two barycentres: the chase
    t, xc = series(TRAJ, "pair_pm", "x_canon")
    _, xp = series(TRAJ, "pair_pm", "x_phantom")
    xc, xp = xc - 32.0, xp - 32.0
    ax.fill_between(t, xp, xc, color="0.93", lw=0)
    ax.plot(t, xc, color="k", ls="-", lw=1.0)
    ax.plot(t, xp, color="k", ls="--", lw=1.0)
    ax.text(4, 6.5, r"canonical, $M_+>0$", color="k", fontsize=7.3)
    ax.text(3, -6.4, r"phantom, $M_-<0$", color="k", fontsize=7.3)
    # third arrow sits at t=55, where the wedge is wide enough that its label
    # can be centered on the arrow without touching either worldline
    for tt, side, dy in ((5.0, "r", 0.0), (40.0, "r", 0.0), (55.0, "l", 0.0)):
        a = float(np.interp(tt, t, xc))
        b = float(np.interp(tt, t, xp))
        ax.annotate("", xy=(tt, b), xytext=(tt, a),
                    arrowprops=dict(arrowstyle="<->", lw=0.6, color="0.35",
                                    shrinkA=0, shrinkB=0))
        if side == "r":
            ax.annotate(f"${a - b:.1f}$", xy=(tt + 1.6, (a + b) / 2 + dy),
                        fontsize=7.0, color="0.35", va="center", ha="left")
        else:
            ax.annotate(f"${a - b:.1f}$", xy=(tt - 1.6, (a + b) / 2 + dy),
                        fontsize=7.0, color="0.35", va="center", ha="right")
    ax.set_xlim(0, 63)
    ax.set_ylim(-7.5, 9.5)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$x - x_c$")
    ax.set_title("(a)", loc="left")

    # (b) drifts, with the null controls; solid/dashed = reference canonical/
    # phantom, dash-dot/dash-dot-dot = equal-|ADM| rerun, dotted = controls.
    # Runs are distinguished by LINE STYLE, not by gray level alone: closely
    # overlapping gray-on-black curves blend together in print.
    for cell, key, color, ls, lw in [
        ("pair_pm", "drift_canon", "k", "-", 1.1),
        ("pair_pm", "drift_phantom", "k", "--", 1.1),
        ("pair_pm_eqm", "drift_canon", "0.35", "-.", 0.9),
        ("pair_pm_eqm", "drift_phantom", "0.35",
         (0, (4, 1.2, 1, 1.2, 1, 1.2)), 0.9),
    ]:
        tt, d = series(TRAJ, cell, key)
        bx.plot(tt, d, color=color, ls=ls, lw=lw)
    for cell, key in [("pair_pp", "drift_canon"), ("pair_mm", "drift_phantom"),
                      ("single_p", "drift_canon"), ("single_m", "drift_phantom")]:
        tt, d = series(TRAJ, cell, key)
        bx.plot(tt, d, color="0.45", ls=(0, (1, 1.5)), lw=0.7)
    bx.axvline(30, ymax=0.82, color="0.45", ls=":", lw=0.7)
    bx.text(44.0, 6.0, r"$\Phi_-$", color="k", fontsize=8.5)
    bx.text(49.5, 1.15, r"$\Phi_+$", color="k", fontsize=8.5)
    bx.text(2.5, 9.15, r"dash--dot: equal-$|M|$ rerun", fontsize=7.0,
            color="0.35", va="top")
    # label the controls by proximity, late, where the flat gray lines are the
    # only curves in the neighbourhood (a leader would have to cross curves)
    bx.text(60.5, 0.45, "controls", fontsize=7.0, color="0.35", ha="right")
    bx.set_xlim(0, 62)
    bx.set_ylim(-0.6, 10)
    bx.set_xlabel("$t$")
    bx.set_ylabel(r"$\Delta x$")
    bx.set_title("(b)", loc="left")

    # (c) gap: NR closes it, the point-mass model does not.  Four curves,
    # four line styles (the equal-|M| pair in dark gray as a secondary cue).
    for cell, color, ls, label in [
        ("pair_pm", "0.15", "-", "NR"),
        ("pair_pm_eqm", "0.35", "-.", r"NR, equal $|M|$"),
    ]:
        tt, s = series(TRAJ, cell, "separation")
        cx.plot(tt, np.abs(s), color=color, ls=ls, label=label)
    for cell, color, ls, label in [
        ("pair_pm", "0.15", "--", "point mass"),
        ("pair_pm_eqm", "0.35", (0, (1, 1.5)), r"point mass, equal $|M|$"),
    ]:
        tt, s = series(NEWT, cell, "sep_pointmass")
        cx.plot(tt, s, color=color, ls=ls, lw=0.9, label=label)
    cx.legend(loc="lower left", borderaxespad=0.4)
    cx.set_xlim(0, 62)
    cx.set_ylim(0, 9.3)
    cx.set_xlabel("$t$")
    cx.set_ylabel("$d$")
    cx.set_title("(c)", loc="left")

    fig.subplots_adjust(wspace=0.38)
    fig.savefig(os.path.join(OUT, "fig_trajectories.pdf"))
    plt.close(fig)


# ------------------------------------------------------------------ fig: mirror
def fig_mirror():
    fig, ax = plt.subplots(figsize=(SINGLE, 2.55))
    for key, ls, label in [("drift_canon", "-", "canonical, reference"),
                           ("drift_phantom", "--", "phantom, reference")]:
        t, d = series(TRAJ, "pair_pm", key)
        ax.plot(t, d, color="k", ls=ls, lw=1.0, label=label)
    for key, marker, label in [
        ("drift_canon", "o", r"canonical, mirror $\times(-1)$"),
        ("drift_phantom", "s", r"phantom, mirror $\times(-1)$"),
    ]:
        t, d = series(TRAJ, "pair_mp_mirror", key)
        ax.plot(t[::2], -d[::2], ls="none", marker=marker, ms=3.0,
                markerfacecolor="white", markeredgewidth=0.8,
                markeredgecolor="k", label=label)
    ax.legend(loc="upper left", borderaxespad=0.4)
    ax.set_xlim(0, 62)
    ax.set_ylim(-0.4, 10)
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

    ax.axhline(0, color="0.5", lw=0.5)
    ax.axvspan(0.500, 0.5265, color="0.90", zorder=0, lw=0)
    ax.text(0.5133, -0.012, "no dressed star", rotation=90, ha="center",
            va="center", fontsize=7.0, color="0.4", zorder=4,
            bbox=dict(facecolor="0.90", edgecolor="none", pad=1.2))

    ax.plot([0.550, 0.56598], [0.063951, -0.063950], marker="*", ms=9,
            color="k", ls="none", zorder=5)
    ax.annotate(r"equal-$|M|$ pair", xy=(0.5705, -0.0655),
                xytext=(0.628, -0.0885), fontsize=7.0, color="0.2",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=1, shrinkB=2))
    ax.legend(loc="upper right", borderaxespad=0.6)
    ax.set_xlim(0.50, 0.82)
    ax.set_ylim(-0.108, 0.092)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$M_{\rm ADM}$")
    fig.savefig(os.path.join(OUT, "fig_family.pdf"))
    plt.close(fig)


# ------------------------------------------------------------ fig: schematic
def fig_schematic():
    """Bondi's sign rules, Forward-style: the three pairings of active mass."""
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots(figsize=(SINGLE, 2.1))
    ax.set_aspect("equal")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.25)
    ax.axis("off")

    R = 0.52
    XL, XR = 2.0, 5.0
    Y = {"a": 5.2, "b": 3.15, "c": 1.05}

    def body(x, y, sector):
        if sector == "+":
            fc, ec, lab, ls = "#E9F2F9", C_CANON, "$+M$", "-"
        else:
            fc, ec, lab, ls = "#FDF0E6", C_PHANT, "$-M$", (0, (3, 1.6))
        ax.add_patch(mpatches.Circle((x, y), R, facecolor=fc, edgecolor=ec,
                                     lw=1.0, linestyle=ls))
        ax.text(x, y, lab, ha="center", va="center", fontsize=6.8, color=ec)

    def force(x0, x1, y):
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="-|>", color="k", lw=0.9,
                                    shrinkA=0, shrinkB=0, mutation_scale=8))

    def verdict(y, text, bold=False):
        ax.text(7.15, y, text, ha="left", va="center", fontsize=7.3,
                fontweight="bold" if bold else "normal")

    def letter(y, s):
        ax.text(0.1, y + 0.75, s, ha="left", va="center", fontsize=7.5)

    # (a) + + : each attracted by the other
    y = Y["a"]
    letter(y, "(a)")
    body(XL, y, "+"); body(XR, y, "+")
    force(XL + R + 0.1, XL + R + 0.78, y)
    force(XR - R - 0.1, XR - R - 0.78, y)
    verdict(y, "fall together")

    # (b) - - : each repelled by the other
    y = Y["b"]
    letter(y, "(b)")
    body(XL, y, "-"); body(XR, y, "-")
    force(XL - R - 0.1, XL - R - 1.0, y)
    force(XR + R + 0.1, XR + R + 1.0, y)
    verdict(y, "push apart")

    # (c) - + : the Bondi dipole; both forces point the same way
    y = Y["c"]
    ax.axhspan(y - 0.95, y + 0.95, color="0.955", zorder=0, lw=0)
    letter(y, "(c)")
    body(XL, y, "-"); body(XR, y, "+")
    force(XL + R + 0.1, XL + R + 1.0, y)
    force(XR + R + 0.1, XR + R + 1.0, y)
    ax.text(XL + R + 0.68, y + 0.62, "attracted", ha="center", va="bottom",
            fontsize=6.0, color="0.25")
    ax.text(XR + R + 0.68, y + 0.62, "repelled", ha="center", va="bottom",
            fontsize=6.0, color="0.25")
    verdict(y, "run away", bold=True)

    fig.savefig(os.path.join(OUT, "fig_schematic.pdf"))
    plt.close(fig)


# ------------------------------------------------------------ fig: sep scaling
def fig_sepscaling():
    fig, ax = plt.subplots(figsize=(SINGLE, 2.7))
    styles = {8: "-", 12: "--", 16: ":"}
    for gap, cell in [(8, "pair_pm"), (12, "pair_pm_sep12"),
                      (16, "pair_pm_sep16")]:
        t, d = series(TRAJ, cell, "drift_phantom")
        keep = (t >= 14) & (t <= 56)
        ax.plot(t[keep], d[keep], color="k", ls=styles[gap], lw=1.0,
                label=f"$d_0={gap}$")
        # tie each matched point-mass value to its own NR curve
        rows = [r for r in SEPS if int(r["separation"]) == gap]
        for r in rows:
            tt = float(r["t"])
            pm = float(r["drift_phantom_pointmass"])
            nr = float(np.interp(tt, t, d))
            ax.plot([tt, tt], [pm, nr], color="0.75", lw=0.5, zorder=1)
        ax.plot([float(r["t"]) for r in rows],
                [float(r["drift_phantom_pointmass"]) for r in rows],
                color="k", ls="none", marker="D", ms=2.8,
                markerfacecolor="white", markeredgewidth=0.7, zorder=3)

    handles = [Line2D([], [], color="k", ls=styles[g], lw=1.0,
                      label=f"$d_0={g}$") for g in (8, 12, 16)]
    handles.append(Line2D([], [], color="k", lw=0, marker="D", ms=2.8,
                          markerfacecolor="white", markeredgewidth=0.7,
                          label="point mass, matched"))
    ax.legend(handles=handles, loc="upper left", borderaxespad=0.5)
    ax.set_xlim(13, 58)
    ax.set_yscale("log")
    ax.set_ylim(0.04, 12)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$\Delta x_{\Phi_-}$")
    fig.savefig(os.path.join(OUT, "fig_sepscaling.pdf"))
    plt.close(fig)


# ------------------------------------------------------------ fig: constraints
def read_norms(cell):
    """t, composite L2 Ham, composite L2 Mom, composite Linf Ham.

    Columns 12/13/16 of constraint_norms.dat: the whole-hierarchy norms
    (covered coarse cells dropped, outermost 4 coarse cells excluded), which
    the diagnostic code marks as the accuracy figures; cols 2/3 are the
    level-0-only governor input.
    """
    path = os.path.join(PACK, "data", cell, "constraint_norms.dat")
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
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(SINGLE, 3.6), sharex=True)

    t, ham, mom, linf = read_norms("pair_pm")
    tp, hamp, momp, _ = read_norms("single_p")
    tm, hamm, momm, _ = read_norms("single_m")

    # (a) Hamiltonian: L2 for all cells, Linf for the reference cell.
    # Controls in dotted / tightly dashed gray: style, not shade, carries them
    ax.plot(t, linf, color="k", ls="-.", lw=0.8)
    ax.plot(tp, hamp, color="0.45", ls=(0, (1, 1.5)), lw=0.7)
    ax.plot(tm, hamm, color="0.45", ls=(0, (3, 1.2)), lw=0.7)
    ax.plot(t, ham, color="k", ls="-", lw=1.0)
    ax.axvline(30, color="0.45", ls=":", lw=0.7)
    ax.text(10.0, 3.2e-2, r"$L_\infty$", fontsize=7.0, color="0.2")
    ax.text(50.0, 2.3e-3, r"$L_2$", fontsize=7.0, color="0.2")
    ax.annotate("single-star controls", xy=(36.5, 2.05e-3),
                xytext=(38.5, 3.1e-4), fontsize=7.0, color="0.35",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=1, shrinkB=2))
    ax.set_yscale("log")
    ax.set_ylim(4e-5, 3)
    ax.set_ylabel(r"$\mathcal{H}$")
    ax.set_title("(a)", loc="left")

    # (b) momentum: L2 only
    bx.plot(tp, momp, color="0.45", ls=(0, (1, 1.5)), lw=0.7)
    bx.plot(tm, momm, color="0.45", ls=(0, (3, 1.2)), lw=0.7)
    bx.plot(t, mom, color="k", ls="-", lw=1.0)
    bx.axvline(30, color="0.45", ls=":", lw=0.7)
    bx.annotate("controls", xy=(28.0, 1.35e-4), xytext=(13.0, 3.6e-5),
                fontsize=7.0, color="0.35",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45",
                                shrinkA=1, shrinkB=2))
    bx.set_yscale("log")
    bx.set_ylim(1.5e-6, 3e-3)
    bx.set_xlim(0, 62)
    bx.set_xlabel("$t$")
    bx.set_ylabel(r"$\mathcal{M}$")
    bx.set_title("(b)", loc="left")

    fig.subplots_adjust(hspace=0.24)
    fig.savefig(os.path.join(OUT, "fig_constraints.pdf"))
    plt.close(fig)


# ------------------------------------------------------------------ fig: velocity
def read_barycenters(cell):
    """t, x_canon, x_phantom from the continuously sampled barycentre stream."""
    path = os.path.join(PACK, "data", cell, "sector_barycenters.dat")
    d = np.loadtxt(path)
    return d[:, 0], d[:, 2], d[:, 7]


def central_diff(t, x, half=5):
    """Centred difference over +-half samples (+-2.0 time units at dt=0.4)."""
    v = np.full_like(x, np.nan)
    v[half:-half] = (x[2 * half:] - x[:-2 * half]) / (t[2 * half:] - t[:-2 * half])
    return v


def fig_velocity():
    fig, ax = plt.subplots(figsize=(SINGLE, 2.45))

    t, xc, xp = read_barycenters("pair_pm")
    ax.plot(t, central_diff(t, xc), color="k", ls="-", lw=1.1)
    ax.plot(t, central_diff(t, xp), color="k", ls="--", lw=1.1)

    # matched point-mass model, for scale: dotted (canonical) and dash-dotted
    # (phantom) so it survives grayscale print next to the black NR curves
    tn, dc = series(NEWT, "pair_pm", "drift_canon_pointmass")
    _, dp = series(NEWT, "pair_pm", "drift_phantom_pointmass")
    ax.plot(tn, np.gradient(dc, tn), color="0.4", ls=(0, (1, 1.5)), lw=0.9)
    ax.plot(tn, np.gradient(dp, tn), color="0.4", ls="-.", lw=0.9)

    ax.axvline(30, color="0.45", ls=":", lw=0.7)
    ax.text(37.5, 0.30, r"$\Phi_-$", fontsize=8.5)
    ax.text(50.5, 0.15, r"$\Phi_+$", fontsize=8.5)
    ax.text(44.0, 0.016, "point mass", fontsize=7.0, color="0.35")

    ax.set_xlim(0, 62)
    ax.set_ylim(-0.02, 0.48)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$v_x$")
    fig.savefig(os.path.join(OUT, "fig_velocity.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- fig: weyl
def read_l2_modes(cell):
    """t and the five complex (l=2,m) amplitudes r*psi4 at each shell.

    psi4_mode_l2_all.dat columns: time, then (Re,Im) pairs for m=-2..2 at
    R=8, then the same ten columns at R=16.
    """
    path = os.path.join(PACK, "data", cell, "psi4_mode_l2_all.dat")
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
    fig, (wx, ax, bx) = plt.subplots(3, 1, figsize=(SINGLE, 5.5), sharex=True)

    # (a) the raw waveform at the outer shell: a monotone, chirpless ramp in
    # the mixed run against the non-secular wander of the same-sign controls
    t_pm, modes_pm = read_l2_modes("pair_pm")
    wx.plot(t_pm, 1e3 * np.real(modes_pm[(16, 0)]), color="k", ls="-", lw=1.1,
            zorder=5, label=r"\texttt{PM}")
    for cell, ls in (("pair_pp", (0, (1, 1.5))), ("pair_mm", (0, (3, 1.2)))):
        tc, mc = read_l2_modes(cell)
        wx.plot(tc, 1e3 * np.real(mc[(16, 0)]), color="0.45", ls=ls, lw=0.7,
                label=rf"\texttt{{{cell[5:].upper()}}}")
    wx.legend(loc="lower left", borderaxespad=0.4)
    wx.set_ylabel(r"$10^3\,r\,\mathrm{Re}\,\psi_4^{2,0}$")
    wx.set_title("(a)", loc="left")

    # (b) total l=2 amplitude at R=16: mixed cells vs the same-sign controls.
    # Style, not shade, separates the runs (as everywhere in the article).
    A_pm = l2_amplitude(modes_pm, 16)
    ax.plot(t_pm, A_pm, color="k", ls="-", lw=1.1, zorder=5)
    t_eq, modes_eq = read_l2_modes("pair_pm_eqm")
    ax.plot(t_eq, l2_amplitude(modes_eq, 16), color="0.35", ls="-.", lw=0.9,
            zorder=4)
    t_mp, modes_mp = read_l2_modes("pair_mp_mirror")
    A_mp = l2_amplitude(modes_mp, 16)
    ax.plot(t_mp[::5], A_mp[::5], ls="none", marker="o", ms=2.7,
            markerfacecolor="white", markeredgewidth=0.7,
            markeredgecolor="k", zorder=6)
    for cell, ls in (("pair_pp", (0, (1, 1.5))), ("pair_mm", (0, (3, 1.2)))):
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
    ax.legend(handles=handles, loc="lower right", borderaxespad=0.4, ncol=2,
              columnspacing=1.0)
    ax.set_yscale("log")
    ax.set_ylim(3e-6, 3e-2)
    ax.set_ylabel(r"$A_{\ell=2}$")
    ax.set_title("(b)", loc="left")

    # (c) the reference-run amplitude against the rescaled pair speed: the
    # near-zone quadrupole curvature rises in lockstep with the runaway.
    tb, xc, xp = read_barycenters("pair_pm")
    v = central_diff(tb, 0.5 * (xc + xp))
    Ai = np.interp(tb, t_pm, A_pm)
    ok = np.isfinite(v) & (tb >= 5) & (tb <= 55)
    kappa = float(np.sum(Ai[ok] * np.abs(v[ok])) / np.sum(v[ok] ** 2))
    corr = float(np.corrcoef(Ai[ok], np.abs(v[ok]))[0, 1])

    bx.plot(t_pm, A_pm, color="k", ls="-", lw=1.1,
            label=r"$A_{\ell=2}$, run \texttt{PM}")
    bx.plot(tb, kappa * np.abs(v), color="0.4", ls="--", lw=0.9,
            label=rf"$\kappa\,|\dot X|$, $\kappa={kappa:.3f}$")
    bx.axvline(30, color="0.45", ls=":", lw=0.7)
    bx.legend(loc="upper left", borderaxespad=0.4)
    bx.set_xlim(0, 62)
    bx.set_ylim(-0.0002, 0.0068)
    bx.set_xlabel("$t$")
    bx.set_ylabel(r"$A_{\ell=2}$")
    bx.set_title("(c)", loc="left")

    fig.subplots_adjust(hspace=0.14)
    fig.savefig(os.path.join(OUT, "fig_weyl.pdf"))
    plt.close(fig)

    # ---- numbers quoted in the article text -------------------------------
    print("fig_weyl stats:")
    print(f"  corr(A16,|v|) 5<=t<=55 (central diff +-2): {corr:.3f}")
    print(f"  kappa: {kappa:.4f}")
    for cell, (tc, mc) in (("pair_pm", (t_pm, modes_pm)),
                           ("pair_pm_eqm", (t_eq, modes_eq)),
                           ("pair_mp_mirror", (t_mp, modes_mp))):
        A = l2_amplitude(mc, 16)
        floor = np.mean(A[(tc >= 0) & (tc <= 10)])
        late = np.mean(A[(tc >= 55) & (tc <= 60)])
        print(f"  growth {cell}: {late / floor:.0f}x "
              f"({floor:.2e} -> {late:.2e})")
    n = min(len(t_pm), len(t_mp))
    rel = np.max(np.abs(A_pm[:n] - A_mp[:n])) / np.max(A_pm[:n])
    print(f"  mirror amplitude reproduction: max rel {rel:.3f}")
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


if __name__ == "__main__":
    fig_chase_frames()
    fig_trajectories()
    fig_mirror()
    fig_schematic()
    fig_family()
    fig_sepscaling()
    fig_velocity()
    fig_weyl()
    fig_constraints()
    print("wrote figures to", OUT)
