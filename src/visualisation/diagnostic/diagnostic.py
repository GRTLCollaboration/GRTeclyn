#!/usr/bin/env python3
"""
Plot streaming collapse diagnostics written by WormholeCollapse.

Expected input: collapse_diagnostics.dat (SmallDataIO ASCII), typically located at:
  <run_output>/data/collapse_diagnostics.dat

Columns (current format):
  time  min_lapse  min_chi  max_abs_K  min_lapse_x  min_lapse_y  min_lapse_z  max_ah_r  min_theta_plus  r_at_min_theta_plus

Optional additional input: areal_radius.dat (from consume_plotfiles.py --areal-radius)
  time  R_areal_min  r_at_R_areal_min

This script produces a publication-style multi-panel figure with:
  - Rows 1-2: standard collapse diagnostics
  - Row 3 (when areal_radius.dat available): R_areal_min, expansion velocity, K-decay lifetime
"""

from __future__ import annotations

import argparse
import string
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

M_SUN_SEC = 4.9255e-6


def _default_run_dir() -> Path:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent.parent
    d2 = (project_root.parent / "data_2gpu").resolve()
    if d2.exists():
        return d2
    return (project_root.parent / "data").resolve()


def _resolve_input_path(data_dir: Path, explicit_input: str | None) -> Path:
    if explicit_input:
        p = Path(explicit_input).expanduser().resolve()
        if not p.exists():
            raise SystemExit(f"Input not found: {p}")
        return p

    candidates = [
        data_dir / "data" / "collapse_diagnostics.dat",
        data_dir / "collapse_diagnostics.dat",
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise SystemExit(
        "Could not find collapse_diagnostics.dat. Tried:\n"
        + "\n".join([f"  - {c}" for c in candidates])
        + "\nProvide it explicitly as a positional argument, or set --data to your run directory."
    )


def _resolve_areal_path(data_dir: Path, explicit: str | None) -> Optional[Path]:
    if explicit:
        p = Path(explicit).expanduser().resolve()
        return p if p.exists() else None

    candidates = [
        data_dir / "small_data" / "areal_radius.dat",
        data_dir / "areal_radius.dat",
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    return None


def load_collapse_diagnostics(path: Path) -> Dict[str, np.ndarray]:
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            rows.append(vals)

    if not rows:
        raise SystemExit(f"No data rows found in {path}")

    arr = np.asarray(rows, dtype=float)
    if arr.shape[1] not in (4, 7, 8, 10, 14):
        raise SystemExit(
            f"Unexpected number of columns in {path}: got {arr.shape[1]}, expected 4, 7, 8, 10, or 14"
        )

    t = arr[:, 0]
    out: Dict[str, np.ndarray] = {
        "t": t,
        "min_lapse": arr[:, 1],
        "min_chi": arr[:, 2],
        "max_abs_K": arr[:, 3],
    }
    if arr.shape[1] >= 7:
        out.update({
            "min_lapse_x": arr[:, 4],
            "min_lapse_y": arr[:, 5],
            "min_lapse_z": arr[:, 6],
        })
    else:
        out.update({
            "min_lapse_x": np.full_like(t, np.nan),
            "min_lapse_y": np.full_like(t, np.nan),
            "min_lapse_z": np.full_like(t, np.nan),
        })

    out["max_ah_r"] = arr[:, 7] if arr.shape[1] >= 8 else np.full_like(t, np.nan)

    if arr.shape[1] >= 10:
        out["min_theta_plus"] = arr[:, 8]
        out["r_at_min_theta_plus"] = arr[:, 9]
    else:
        out["min_theta_plus"] = np.full_like(t, np.nan)
        out["r_at_min_theta_plus"] = np.full_like(t, np.nan)

    if arr.shape[1] >= 14:
        out["min_phi"] = arr[:, 10]
        out["max_phi"] = arr[:, 11]
        out["min_Pi"] = arr[:, 12]
        out["max_Pi"] = arr[:, 13]

    idx = np.argsort(out["t"])
    for k in list(out.keys()):
        out[k] = out[k][idx]
    return out


def load_areal_radius(path: Path) -> Dict[str, np.ndarray]:
    """Parse areal_radius.dat -> dict with keys t, R_areal_min, r_at_min."""
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) >= 3:
                rows.append(vals[:3])

    if not rows:
        raise SystemExit(f"No data rows found in {path}")

    arr = np.asarray(rows, dtype=float)
    idx = np.argsort(arr[:, 0])
    return {
        "t": arr[idx, 0],
        "R_areal_min": arr[idx, 1],
        "r_at_min": arr[idx, 2],
    }


def _compute_expansion_velocity(
    t: np.ndarray, R_areal: np.ndarray, smooth_window: int = 11
) -> np.ndarray:
    """Numerical derivative dR_areal/dt, optionally smoothed."""
    v = np.gradient(R_areal, t)
    n = len(v)
    w = smooth_window
    if w % 2 == 0:
        w += 1
    if n > w >= 5:
        try:
            v = savgol_filter(v, window_length=w, polyorder=3, mode="interp")
        except Exception:
            pass
    return v


def _exp_decay(t: np.ndarray, A: float, tau: float, C: float) -> np.ndarray:
    return A * np.exp(-t / tau) + C


def _extract_bh_mass(
    t: np.ndarray,
    R_areal: np.ndarray,
    max_abs_K: np.ndarray,
    min_lapse: np.ndarray,
    plateau_fraction: float = 0.3,
) -> Optional[Dict[str, float]]:
    """Extract the final BH mass from late-time plateaux.

    Uses three independent diagnostics:
      1) Areal radius plateau -> M_irr = sqrt(A/16pi) with A = 4*pi*R^2
      2) max|K| plateau -> M from trumpet K = 2.517/M (1+log, Gamma-driver)
      3) min(alpha) plateau -> consistency check

    The 1+log trumpet equilibrium values for Schwarzschild are from
    Hannam et al. (2007): R_trumpet ~ 1.312 M, K_trumpet ~ 2.517 / M.
    """
    n = len(t)
    i_start = int((1.0 - plateau_fraction) * n)
    if i_start >= n - 3:
        return None

    R_late = R_areal[i_start:]
    K_late = max_abs_K[i_start:]
    alpha_late = min_lapse[i_start:]

    R_plateau = float(np.median(R_late[np.isfinite(R_late)]))
    K_plateau = float(np.median(K_late[np.isfinite(K_late) & (K_late > 0)]))
    alpha_plateau = float(np.median(alpha_late[np.isfinite(alpha_late) & (alpha_late > 0)]))

    R_TRUMPET_OVER_M = 1.312
    K_TRUMPET_TIMES_M = 2.517

    M_from_R = R_plateau / R_TRUMPET_OVER_M if R_plateau > 0 else 0.0
    M_from_K = K_TRUMPET_TIMES_M / K_plateau if K_plateau > 0 else 0.0

    A_AH = 4.0 * np.pi * R_plateau ** 2
    M_irr = np.sqrt(A_AH / (16.0 * np.pi))

    return {
        "R_plateau": R_plateau,
        "K_plateau": K_plateau,
        "alpha_plateau": alpha_plateau,
        "M_from_R": M_from_R,
        "M_from_K": M_from_K,
        "M_irr": M_irr,
        "A_AH": A_AH,
    }


def _compute_lyapunov_exponent(
    t: np.ndarray,
    signal: np.ndarray,
    t_start: float = 0.0,
    t_end_frac: float = 0.3,
) -> Optional[Tuple[float, np.ndarray, np.ndarray]]:
    """Fit exp growth to early-time signal to extract the Lyapunov exponent.

    Fits log|signal| = lambda*t + const during the early growth phase.
    Returns (lambda, t_fit, fit_line) or None if fit fails.
    """
    t_end = t[0] + t_end_frac * (t[-1] - t[0])
    mask = (t >= t_start) & (t <= t_end)
    if np.sum(mask) < 5:
        return None

    t_sel = t[mask]
    sig_sel = np.abs(signal[mask])
    valid = sig_sel > 0
    if np.sum(valid) < 5:
        return None

    log_sig = np.log(sig_sel[valid])
    t_valid = t_sel[valid]

    finite = np.isfinite(log_sig)
    if np.sum(finite) < 5:
        return None

    coeffs = np.polyfit(t_valid[finite], log_sig[finite], 1)
    lyapunov = float(coeffs[0])
    fit_line = np.polyval(coeffs, t_valid[finite])
    return (lyapunov, t_valid[finite], fit_line)


def _fit_K_lifetime(
    t: np.ndarray, max_abs_K: np.ndarray, fit_start_fraction: float = 0.4
) -> Optional[Tuple[float, float, float, np.ndarray, np.ndarray]]:
    """Fit max|K| late-time tail to A*exp(-t/tau) + C.

    Returns (tau, A, C, fit_t, fit_K) or None if fit fails.
    """
    t_start = t[0] + fit_start_fraction * (t[-1] - t[0])
    mask = t >= t_start
    if np.sum(mask) < 5:
        return None

    t_fit = t[mask]
    K_fit = max_abs_K[mask]

    valid = np.isfinite(K_fit) & (K_fit > 0)
    if np.sum(valid) < 5:
        return None
    t_fit = t_fit[valid]
    K_fit = K_fit[valid]

    try:
        A0 = float(K_fit[0] - K_fit[-1])
        tau0 = float(t_fit[-1] - t_fit[0]) / 3.0
        C0 = float(K_fit[-1])
        popt, _ = curve_fit(
            _exp_decay, t_fit, K_fit,
            p0=[max(A0, 0.1), max(tau0, 0.01), max(C0, 0.0)],
            maxfev=5000,
            bounds=([0, 1e-6, 0], [np.inf, np.inf, np.inf]),
        )
        A, tau, C = popt
        return (float(tau), float(A), float(C), t_fit, _exp_decay(t_fit, *popt))
    except Exception:
        return None


def _apply_scientific_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.linewidth": 1.0,
        "grid.alpha": 0.5,
    })


def plot_collapse_diagnostics(
    data: Dict[str, np.ndarray],
    out_path: Path,
    areal_data: Optional[Dict[str, np.ndarray]] = None,
    fit_lifetime: bool = True,
    mass_msun: float = 30.0,
) -> None:
    _apply_scientific_style()

    t = data["t"]
    has_areal = areal_data is not None and len(areal_data.get("t", [])) > 2
    has_phantom = "min_phi" in data and len(data["min_phi"]) > 0

    base_rows = 2
    if has_areal or fit_lifetime:
        base_rows += 1
    if has_phantom:
        base_rows += 1
    
    n_rows = base_rows

    fig, axes = plt.subplots(n_rows, 3, figsize=(15, 4.0 * n_rows), sharex=True)
    axes = np.asarray(axes)
    
    if len(axes.shape) == 1:
        axes = np.expand_dims(axes, axis=0)

    top_specs: Tuple[Tuple[str, str, str], ...] = (
        ("min_lapse", r"$\min(\alpha)$", r"Minimum lapse: $\alpha$"),
        ("min_chi", r"$\min(\chi)$", r"Minimum conformal factor: $\chi$"),
        ("max_abs_K", r"$\max(|K|)$", r"Maximum curvature: $|K|$"),
    )

    for i, (key, ylabel, title) in enumerate(top_specs):
        ax = axes[0, i]
        y = np.asarray(data[key], dtype=float)
        y_plot = np.where(y > 0, y, np.nan)
        ax.semilogy(t, y_plot, color="black", linewidth=1.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, which="both", ls="--", alpha=0.6)
        ax.tick_params(axis="both", which="major", direction="in")

    bottom_specs: Tuple[Tuple[str, str, str], ...] = (
        ("max_ah_r", r"$r_{\rm AH}$", r"Max trapped surface radius: $\theta_+ \leq 0$"),
        ("min_theta_plus", r"$\min(\theta_+)$", r"Minimum null expansion proxy: $\theta_+$"),
        ("r_at_min_theta_plus", r"$r_{\min\theta_+}$", r"Radius at min $\theta_+$"),
    )

    for i, (key, ylabel, title) in enumerate(bottom_specs):
        ax = axes[1, i]
        y = np.asarray(data[key], dtype=float)
        if key == "min_theta_plus":
            ax.plot(t, y, color="black", linewidth=1.5)
            ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.7)
        else:
            ax.plot(t, y, color="black", linewidth=1.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, which="both", ls="--", alpha=0.6)
        ax.tick_params(axis="both", which="major", direction="in")

    if n_rows >= 3:
        ax_areal = axes[2, 0]
        if has_areal:
            ta = areal_data["t"]
            Ra = areal_data["R_areal_min"]
            ax_areal.plot(ta, Ra, color="black", linewidth=1.5)
            ax_areal.set_ylabel(r"$R_{\mathrm{areal,min}}$")
            ax_areal.set_title(r"Throat areal radius $R_{\mathrm{areal}}$")

            print("\n=== Areal Radius ===")
            print(f"  Initial R_areal_min = {Ra[0]:.6f}")
            print(f"  Final   R_areal_min = {Ra[-1]:.6f}")
        else:
            ax_areal.text(0.5, 0.5, "No areal_radius.dat", transform=ax_areal.transAxes,
                          ha="center", va="center", fontsize=11, color="gray")
            ax_areal.set_title(r"Throat areal radius $R_{\mathrm{areal}}$")
        ax_areal.grid(True, which="both", ls="--", alpha=0.6)
        ax_areal.tick_params(axis="both", which="major", direction="in")

        ax_vel = axes[2, 1]
        if has_areal:
            ta = areal_data["t"]
            Ra = areal_data["R_areal_min"]
            vel = _compute_expansion_velocity(ta, Ra)
            ax_vel.plot(ta, vel, color="black", linewidth=1.5)
            ax_vel.axhline(1.0, color="red", linewidth=0.8, linestyle="--", alpha=0.7, label=r"$c=1$")
            ax_vel.axhline(-1.0, color="red", linewidth=0.8, linestyle="--", alpha=0.7)
            ax_vel.axhline(0.0, color="gray", linewidth=0.5, alpha=0.5)
            ax_vel.set_ylabel(r"$dR_{\mathrm{areal}}/dt\;(c)$")
            ax_vel.set_title("Expansion velocity")
            ax_vel.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=9)

            mean_vel = float(np.mean(vel))
            peak_vel = float(vel[np.argmax(np.abs(vel))])
            print("\n=== Expansion Velocity ===")
            print(f"  Mean velocity: {mean_vel:.4f} c")
            print(f"  Peak velocity: {peak_vel:.4f} c")
            T_phys = mass_msun * M_SUN_SEC
            print(f"  (At M = {mass_msun:g} M_sun, c = 1 code unit = {1.0 / T_phys:.2e} code_length/s)")
        else:
            ax_vel.text(0.5, 0.5, "No areal_radius.dat", transform=ax_vel.transAxes,
                        ha="center", va="center", fontsize=11, color="gray")
            ax_vel.set_title("Expansion velocity")
        ax_vel.grid(True, which="both", ls="--", alpha=0.6)
        ax_vel.tick_params(axis="both", which="major", direction="in")

        ax_life = axes[2, 2]
        K_data = np.asarray(data["max_abs_K"], dtype=float)
        ax_life.semilogy(t, np.where(K_data > 0, K_data, np.nan), color="black", linewidth=1.5, label=r"$\max|K|$")

        if fit_lifetime:
            fit_result = _fit_K_lifetime(t, K_data)
            if fit_result is not None:
                tau, A, C, t_fit, K_fitted = fit_result
                ax_life.semilogy(t_fit, K_fitted, color="red", linewidth=1.5, linestyle="--",
                                 label=rf"Fit: $\tau={tau:.3f}$")
                ax_life.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=9)

                T_phys = mass_msun * M_SUN_SEC
                tau_sec = tau * T_phys
                print(f"\n=== K-Decay Lifetime ===")
                print(f"  tau = {tau:.4f} [code units]")
                print(f"  At M = {mass_msun:g} M_sun: tau = {tau_sec:.4e} s")
                print(f"  Fit: A={A:.4f}, C={C:.4f}")
            else:
                print("\n=== K-Decay Lifetime ===")
                print("  Exponential fit did not converge.")

        ax_life.set_ylabel(r"$\max(|K|)$")
        ax_life.set_title(r"$|K|$ decay and lifetime $\tau$")
        ax_life.grid(True, which="both", ls="--", alpha=0.6)
        ax_life.tick_params(axis="both", which="major", direction="in")

    if has_phantom:
        row_idx = base_rows - 1

        ax_phi = axes[row_idx, 0]
        y_min_phi = np.asarray(data["min_phi"], dtype=float)
        y_max_phi = np.asarray(data["max_phi"], dtype=float)
        ax_phi.plot(t, y_max_phi, color="black", linewidth=1.5, linestyle="-", label=r"$\max(\phi)$")
        ax_phi.plot(t, y_min_phi, color="black", linewidth=1.5, linestyle="--", label=r"$\min(\phi)$")
        ax_phi.set_ylabel(r"$\phi$")
        ax_phi.set_title(r"Scalar field profile $\phi$")
        ax_phi.legend(loc="best", frameon=True, framealpha=0.9, fontsize=9)
        ax_phi.grid(True, which="both", ls="--", alpha=0.6)
        ax_phi.tick_params(axis="both", which="major", direction="in")

        ax_pi = axes[row_idx, 1]
        y_min_pi = np.asarray(data["min_Pi"], dtype=float)
        y_max_pi = np.asarray(data["max_Pi"], dtype=float)
        ax_pi.plot(t, y_max_pi, color="black", linewidth=1.5, linestyle="-", label=r"$\max(\Pi)$")
        ax_pi.plot(t, y_min_pi, color="black", linewidth=1.5, linestyle="--", label=r"$\min(\Pi)$")
        ax_pi.set_ylabel(r"$\Pi$")
        ax_pi.set_title(r"Scalar field momentum $\Pi$")
        ax_pi.legend(loc="best", frameon=True, framealpha=0.9, fontsize=9)
        ax_pi.grid(True, which="both", ls="--", alpha=0.6)
        ax_pi.tick_params(axis="both", which="major", direction="in")

        ax_lyap = axes[row_idx, 2]
        if has_areal:
            ta = areal_data["t"]
            Ra = areal_data["R_areal_min"]
            departure = np.abs(Ra - Ra[0])
            departure = np.where(departure > 1e-15, departure, 1e-15)
            
            ax_lyap.semilogy(ta, departure, color="black", linewidth=1.5, label=r"$|\Delta R_{\mathrm{areal}}|$")
            
            lyap_result = _compute_lyapunov_exponent(ta, departure, t_start=0.0, t_end_frac=0.35)
            if lyap_result is not None:
                lam, t_fit_lyap, fit_line_lyap = lyap_result
                ax_lyap.semilogy(t_fit_lyap, np.exp(fit_line_lyap), color="red", linewidth=1.5, linestyle="--",
                                 label=rf"Fit: $\lambda={lam:.4f}$")
                
            ax_lyap.set_ylabel(r"$|\Delta R_{\mathrm{areal}}|$")
            ax_lyap.set_title(r"Instability growth ($\lambda$)")
            ax_lyap.legend(loc="best", frameon=True, framealpha=0.9, fontsize=9)
        else:
            ax_lyap.text(0.5, 0.5, "No areal_radius.dat", transform=ax_lyap.transAxes,
                         ha="center", va="center", fontsize=11, color="gray")
            ax_lyap.set_title(r"Instability growth ($\lambda$)")

        ax_lyap.grid(True, which="both", ls="--", alpha=0.6)
        ax_lyap.tick_params(axis="both", which="major", direction="in")

    K_data = np.asarray(data["max_abs_K"], dtype=float)
    min_lapse = np.asarray(data["min_lapse"], dtype=float)

    if has_areal:
        bh_info = _extract_bh_mass(
            areal_data["t"], areal_data["R_areal_min"], K_data, min_lapse
        )
        if bh_info is not None:
            T_phys = mass_msun * M_SUN_SEC
            print("\n=== Black-Hole Remnant Characterisation ===")
            print(f"  Late-time R_areal plateau : {bh_info['R_plateau']:.6f}")
            print(f"  Late-time max|K| plateau  : {bh_info['K_plateau']:.4f}")
            print(f"  Late-time min(alpha)      : {bh_info['alpha_plateau']:.6f}")
            M_from_R = bh_info["M_from_R"]
            print(f"  --- Mass estimates (code units) ---")
            print(f"  M (from R_trumpet/1.312)  : {M_from_R:.6f}")
            print(f"  M_irr  = sqrt(A_trumpet/16pi): {bh_info['M_irr']:.6f}")
            print(f"  Note: M_irr uses trumpet area, not AH area; R_trumpet < R_AH = 2M")
            print(f"  --- Trumpet verification (Hannam+ 2007 reference values) ---")
            K_pred = 2.517 / M_from_R if M_from_R > 0 else 0
            alpha_pred_range = "0.02--0.05"
            print(f"  Predicted K_trumpet (at M={M_from_R:.4f}) = {K_pred:.4f}  (observed: {bh_info['K_plateau']:.4f})")
            print(f"  Predicted alpha_min (literature)   ~ {alpha_pred_range}  (observed: {bh_info['alpha_plateau']:.6f})")
            ratio = bh_info['K_plateau'] / K_pred if K_pred > 0 else float('inf')
            print(f"  K ratio (observed/predicted) = {ratio:.3f}")
            print(f"  --- Physical mass ---")
            M_phys = M_from_R * mass_msun
            print(f"  At M_total = {mass_msun:g} M_sun: M_BH = {M_phys:.4f} M_sun")

    if has_areal:
        Ra = areal_data["R_areal_min"]
        ta = areal_data["t"]
        departure = np.abs(Ra - Ra[0])
        departure = np.where(departure > 1e-15, departure, 1e-15)
        lyap_result = _compute_lyapunov_exponent(ta, departure, t_start=0.0, t_end_frac=0.35)
        if lyap_result is not None:
            lam, t_fit_lyap, fit_line_lyap = lyap_result
            print(f"\n=== Instability Growth Rate (Lyapunov Exponent) ===")
            print(f"  lambda = {lam:.4f} M^-1  (from early-time |R_areal(t) - R_areal(0)| growth)")
            if lam > 0:
                print(f"  e-folding time = {1.0/lam:.4f} M")
                T_phys = mass_msun * M_SUN_SEC
                t_efold_phys = (1.0 / lam) * T_phys
                print(f"  At M = {mass_msun:g} M_sun: lambda = {lam / T_phys:.4e} s^-1, t_e-fold = {t_efold_phys:.4e} s")
            else:
                print(f"  (decay, not growth -- collapse too rapid for clear exponential phase)")

    for ax, letter in zip(axes.flatten(), string.ascii_lowercase):
        title = ax.get_title()
        if title:
            ax.set_title(f"({letter}) {title}")

    for i in range(axes.shape[1]):
        axes[n_rows - 1, i].set_xlabel(r"$t$")

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".eps"), dpi=600, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    print(f"\nSaved: {out_path.with_suffix('.png')}, .eps, and .pdf")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    default_data = _default_run_dir()

    parser = argparse.ArgumentParser(
        description="Plot collapse diagnostics with areal radius and K-decay lifetime."
    )
    parser.add_argument(
        "input", nargs="?", default=None,
        help="Path to collapse_diagnostics.dat (optional; otherwise inferred from --data).",
    )
    parser.add_argument(
        "--data", default=str(default_data),
        help="Run output directory (expects data/collapse_diagnostics.dat inside).",
    )
    parser.add_argument(
        "--out", default=str(script_dir),
        help="Output directory for collapse_diagnostics_plot.eps",
    )
    parser.add_argument(
        "--name", default="collapse_diagnostics_plot.eps",
        help="Output filename (default: collapse_diagnostics_plot.eps)",
    )
    parser.add_argument(
        "--areal-radius-file", default=None,
        help="Explicit path to areal_radius.dat (auto-detected from --data by default).",
    )
    parser.add_argument(
        "--no-fit-lifetime", action="store_true",
        help="Disable exponential fit to K decay.",
    )
    parser.add_argument(
        "--mass-msun", type=float, default=30.0,
        help="Total mass in solar masses for physical unit conversion (default: 30).",
    )
    args = parser.parse_args()

    data_dir = Path(args.data).expanduser().resolve()
    in_path = _resolve_input_path(data_dir, args.input)
    data = load_collapse_diagnostics(in_path)

    areal_path = _resolve_areal_path(data_dir, args.areal_radius_file)
    areal_data = None
    if areal_path is not None:
        try:
            areal_data = load_areal_radius(areal_path)
            print(f"Loaded areal radius data from {areal_path}")
        except Exception as e:
            print(f"WARNING: Could not load areal radius data: {e}")

    out_dir = Path(args.out).expanduser().resolve()
    out_path = out_dir / str(args.name)
    plot_collapse_diagnostics(
        data,
        out_path=out_path,
        areal_data=areal_data,
        fit_lifetime=not args.no_fit_lifetime,
        mass_msun=args.mass_msun,
    )


if __name__ == "__main__":
    main()
