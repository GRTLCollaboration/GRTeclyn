import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import yt
from scipy.signal import savgol_filter, welch

from grteclyn_wrapper.core.config import default_sim_data_dir

def _fmt_for_filename(x: float) -> str:
    """Compact float formatting for filenames (avoid '.' and '-')."""
    s = f"{x:g}"
    s = s.replace("-", "m").replace(".", "p")
    return s

def _smooth_psd(freqs: np.ndarray, psd: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    """
    Smooth PSD in log space to reduce jaggedness:
    smooth log10(PSD) vs index (Welch freqs are monotone).
    """
    psd = np.asarray(psd, dtype=float)
    out = psd.copy()

    # Only smooth where PSD is positive and finite
    m = np.isfinite(psd) & (psd > 0)
    if np.sum(m) < 7:
        return out

    y = np.log10(psd[m])
    n = y.size

    # Ensure valid Savitzky-Golay parameters
    w = int(window)
    if w < 5:
        return out
    if w % 2 == 0:
        w += 1
    if w > n:
        w = n if (n % 2 == 1) else (n - 1)
    p = int(polyorder)
    if p < 1:
        return out
    if p >= w:
        p = max(1, w - 2)

    y_s = savgol_filter(y, window_length=w, polyorder=p, mode="interp")
    out[m] = 10 ** y_s
    return out

def get_extraction_points(radius, n_points=24):
    """Generates coordinates and weights for integration on a sphere."""
    theta = np.linspace(0, np.pi, n_points)
    # Important: don't include phi=2pi endpoint (duplicates phi=0)
    phi = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')

    X = radius * np.sin(THETA) * np.cos(PHI)
    Y = radius * np.sin(THETA) * np.sin(PHI)
    Z = radius * np.cos(THETA)

    # Integration weights dOmega
    dtheta = np.pi / (n_points - 1)
    dphi = 2*np.pi / n_points
    Weights = np.sin(THETA) * dtheta * dphi

    return X, Y, Z, THETA, PHI, Weights

def spin_weighted_sph_harm_2_0(theta):
    # Exact (real) spin-weighted spherical harmonic:  _-2Y_{2,0}
    # _-2Y20(theta,phi) = sqrt(15/(32*pi)) * sin^2(theta)
    return np.sqrt(15.0 / (32.0 * np.pi)) * (np.sin(theta) ** 2)

def main():
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _default_data = str(default_sim_data_dir())

    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default=_default_data, help="Directory containing WormholePlt* / plt*")
    parser.add_argument("--radii", type=float, nargs='+', default=[10.0, 20.0, 30.0], help="Extraction radii")
    parser.add_argument("--out", type=str, default=_script_dir, help="Output directory")
    parser.add_argument("--n-points", type=int, default=24, help="Angular resolution N (N×N points)")
    parser.add_argument(
        "--time-axis",
        choices=["simulation", "retarded"],
        default="simulation",
        help="X-axis for plots. 'simulation' uses t, 'retarded' uses (t - R_ext).",
    )
    parser.add_argument("--t-min", type=float, default=0.0, help="Min x value for plot (in chosen time-axis units).")
    parser.add_argument("--t-max", type=float, default=50.0, help="Max x value for plot (in chosen time-axis units).")
    parser.add_argument(
        "--psd-smooth-window",
        type=int,
        default=21,
        help="Savitzky-Golay smoothing window for PSD (odd integer; applied in log10-PSD space).",
    )
    parser.add_argument(
        "--psd-smooth-polyorder",
        type=int,
        default=5,
        help="Savitzky-Golay polynomial order for PSD smoothing.",
    )
    parser.add_argument(
        "--psd-hide-raw",
        action="store_true",
        help="Hide raw (unsmoothed) PSD points; show only smoothed curve.",
    )
    args = parser.parse_args()

    # Find files safely
    raw = glob.glob(os.path.join(args.data, "WormholePlt*")) + glob.glob(os.path.join(args.data, "plt*"))
    files = sorted(f for f in raw if os.path.isdir(f))

    if not files:
        print("No plotfiles found!")
        return

    results = {R: {'t': [], 'psi4': []} for R in args.radii}
    # Distance from extraction center to the (physical) outer boundary.
    D_boundary = None
    print(f"Processing {len(files)} files...")

    # Sanity check: do plotfiles contain Weyl4?
    try:
        ds0 = yt.load(files[0])
        has_re = ("boxlib", "Weyl4_Re") in ds0.field_list
        has_im = ("boxlib", "Weyl4_Im") in ds0.field_list
        if not (has_re and has_im):
            print("ERROR: Plotfiles do not contain ('boxlib','Weyl4_Re/Im').")
            print("Fix: in your simulation params, set:")
            print("  amr.derive_plot_vars = Weyl4   (record name), or use 'ALL'")
            print("Then re-run the simulation to regenerate plotfiles.")
            return
    except Exception as e:
        print(f"Warning: could not validate fields in first plotfile: {e}")

    # Extraction Loop
    for f in files:
        try:
            ds = yt.load(f)
            t = float(ds.current_time)

            # Robust Center finding
            # center = (np.array(ds.domain_left_edge) + np.array(ds.domain_right_edge)) / 2.0

            # FIX: With symmetry (Octant), the center is at the origin [0,0,0], not the middle of the box
            center = np.array([0.0, 0.0, 0.0])

            left = np.array(ds.domain_left_edge)
            right = np.array(ds.domain_right_edge)

            if D_boundary is None:
                D_boundary = float(np.min(right - center))

            for R in args.radii:
                X, Y, Z, THETA, PHI, W = get_extraction_points(R, n_points=args.n_points)
                sample_x = X.flatten() + center[0]
                sample_y = Y.flatten() + center[1]
                sample_z = Z.flatten() + center[2]

                # Filter points inside domain
                in_domain = (
                    (sample_x >= left[0]) & (sample_x <= right[0]) &
                    (sample_y >= left[1]) & (sample_y <= right[1]) &
                    (sample_z >= left[2]) & (sample_z <= right[2])
                )
                valid_idx = np.where(in_domain)[0]

                if len(valid_idx) < 4: continue

                # Manual Point Sampling (complex Psi4 = Weyl4_Re + i Weyl4_Im)
                weyl = np.full(sample_x.shape, np.nan + 1j * np.nan, dtype=np.complex128)

                for i in valid_idx:
                    try:
                        pt = ds.point([sample_x[i], sample_y[i], sample_z[i]])
                        val_re = pt[("boxlib", "Weyl4_Re")]
                        val_im = pt[("boxlib", "Weyl4_Im")]
                        arr_re = np.asarray(val_re)
                        arr_im = np.asarray(val_im)
                        if arr_re.size > 0 and arr_im.size > 0:
                            weyl[i] = float(arr_re.flat[0]) + 1j * float(arr_im.flat[0])
                    except: pass

                # Renormalize weights for symmetry (Octant -> Full Sphere)
                mask = ~np.isnan(weyl.real)
                if mask.sum() < 4: continue

                weyl = np.where(mask, weyl, 0.0).reshape(X.shape)
                W_masked = np.where(mask, W.flatten(), 0.0)
                omega_tot = W_masked.sum()
                if omega_tot == 0: continue

                # Renormalize to 4pi (handles symmetry automatically)
                W_renorm = (W_masked / omega_tot * 4 * np.pi).reshape(X.shape)

                Y_lm = spin_weighted_sph_harm_2_0(THETA)
                mode_amp = np.sum(weyl * np.conj(Y_lm) * W_renorm) * R

                results[R]['t'].append(t)
                results[R]['psi4'].append(mode_amp)

            if len(results[args.radii[0]]['t']) % 5 == 0:
                print(f"Processed t={t:.2f}")

        except Exception as e: print(e)

    # Plotting — LaTeX scientific paper style
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
        "grid.alpha": 0.5,
        "legend.fontsize": 12,
        "figure.figsize": (10, 8),
    })
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Distinct color palettes for Psi4 (blues) vs PSD (greens)
    n_radii = len(args.radii)
    psi4_colors = plt.cm.Blues(np.linspace(0.5, 0.85, n_radii))[::-1]
    psd_colors = plt.cm.Greens(np.linspace(0.5, 0.85, n_radii))[::-1]

    plotted_any = False

    for i, R in enumerate(args.radii):
        t = np.array(results[R]['t'])
        psi4 = np.array(results[R]['psi4'])

        # Allow plotting even for short runs; PSD needs more points.
        if len(t) < 2:
            print(f"Warning: R={R}: not enough samples to plot (n={len(t)}).")
            continue

        # Sort
        idx = np.argsort(t)
        t = t[idx]
        psi4 = psi4[idx]

        # --- FIX: Boundary Reflection Cutoff ---
        if D_boundary is not None:
            t_cutoff = 2 * D_boundary - R
            valid_idx = t < t_cutoff
            if np.sum(valid_idx) < 20:
                print(f"Warning: R={R} is too close to boundary or simulation too short. Valid data points: {np.sum(valid_idx)}")
                # If too short, plot what we have anyway, just warn
            else:
                t = t[valid_idx]
                psi4 = psi4[valid_idx]
                print(f"R={R}: Truncated data at t={t_cutoff:.1f} to remove boundary reflections.")

        # Retarded time
        t_ret = t - R
        x = t if args.time_axis == "simulation" else t_ret

        # Plot Psi4 (Real part of 2,0 mode)
        y = np.real(psi4)
        m = np.ones_like(x, dtype=bool)
        if args.t_min is not None:
            m &= x >= args.t_min
        if args.t_max is not None:
            m &= x <= args.t_max
        ax1.plot(x[m], y[m], color=psi4_colors[i], label=f"R={R}")

        # Spectrum of Psi4: raw as dots, smoothed as line
        dt = np.mean(np.diff(t))
        if dt > 0:
            fs = 1.0/dt
            if len(psi4) >= 8:
                freqs, psd = welch(np.real(psi4), fs, nperseg=min(len(psi4)//2, 128))
                psd_s = _smooth_psd(freqs, psd, window=args.psd_smooth_window, polyorder=args.psd_smooth_polyorder)
                c = psd_colors[i]
                # Plot raw data as dots
                if not args.psd_hide_raw:
                    ax2.semilogy(freqs, psd, "o", color="red", markersize=3.0, markeredgewidth=0, alpha=0.4, zorder=1, label=None)
                # Smoothed data as green line
                ax2.semilogy(freqs, psd_s, "-", color=c, linewidth=1.5, zorder=2, label=f"R={R}")

    radius_suffix = "_".join(f"R{R}" for R in args.radii)
    output_path = os.path.join(args.out, f"psi4_analysis_{radius_suffix}_n{args.n_points}.png")

    ax1.set_xlabel(r"$t$" if args.time_axis == "simulation" else r"$(t - R_{\mathrm{ext}})$")
    ax1.set_ylabel(r"$r\,\mathrm{Re}\,\Psi_4^{20}$")
    ax1.set_title(r"Gravitational wave curvature $\Psi_4$ (mode $\ell=2$, $m=0$)")
    if plotted_any:
        ax1.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax1.grid(True, which="both", ls="--", alpha=0.6)
    ax1.tick_params(axis='both', which='major', direction='in', top=True, right=True)

    ax2.set_xlabel(r"Frequency ($M^{-1}$)")
    ax2.set_ylabel(r"PSD of $\Psi_4$")
    ax2.grid(True, which="major", ls="--", alpha=0.6)
    ax2.tick_params(axis='both', which='major', direction='in', top=True, right=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    main()