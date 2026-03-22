#!/usr/bin/env python3
"""
Plot extracted Psi4 time-series produced by consume_plotfiles.py.

Features:
  - Time-domain waveform r*Re(Psi4)
  - PSD of Psi4 (--plot-psd)
  - Strain PSD conversion h_tilde(f) = Psi4_tilde(f)/(2*pi*f)^2 (--strain)
  - Physical scaling to Hz and LIGO noise overlay (--strain --mass-msun --distance-mpc)
  - Propagation speed measurement across extraction radii (--propagation-speed)

Input format (whitespace separated):
  time  Re(R=r1) Im(R=r1)  Re(R=r2) Im(R=r2)  ...
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter, welch, find_peaks

matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
M_SUN_KG = 1.98892e30
G_SI = 6.67430e-11
C_SI = 2.99792458e8
M_SUN_SEC = G_SI * M_SUN_KG / C_SI**3      # ~4.9255e-6 s
M_SUN_METER = G_SI * M_SUN_KG / C_SI**2     # ~1476.6 m
MPC_METER = 3.08568e22


# ---------------------------------------------------------------------------
# Existing helpers
# ---------------------------------------------------------------------------
def _default_data_dir() -> Path:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent.parent
    d2 = (project_root.parent / "data_2gpu").resolve()
    if d2.exists():
        return d2
    return (project_root.parent / "data").resolve()


def _smooth_psd(psd: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    psd = np.asarray(psd, dtype=float)
    out = psd.copy()
    m = np.isfinite(psd) & (psd > 0)
    if np.sum(m) < 7:
        return out
    y = np.log10(psd[m])
    n = y.size
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


def _parse_header_radii(header_line: str) -> List[float]:
    rs = re.findall(r"R=([0-9.+-eE]+)", header_line)
    radii: List[float] = []
    for x in rs[::2]:
        try:
            radii.append(float(x))
        except ValueError:
            continue
    out: List[float] = []
    for r in radii:
        if not out or abs(out[-1] - r) > 0:
            out.append(r)
    return out


def load_extracted(path: Path) -> Tuple[np.ndarray, List[float], Dict[float, np.ndarray]]:
    header_radii: List[float] = []
    times: List[float] = []
    cols: List[List[float]] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if not header_radii and "R=" in s:
                    header_radii = _parse_header_radii(s)
                continue
            parts = s.split()
            try:
                t = float(parts[0])
            except ValueError:
                continue
            nums = []
            ok = True
            for p in parts[1:]:
                try:
                    nums.append(float(p))
                except ValueError:
                    ok = False
                    break
            if not ok:
                continue
            times.append(t)
            cols.append(nums)

    if not times:
        raise ValueError(f"No data rows found in {path}")

    t_arr = np.asarray(times, dtype=float)
    data = np.asarray(cols, dtype=float)

    if data.shape[1] % 2 != 0:
        raise ValueError(f"Expected Re/Im pairs; got {data.shape[1]} columns in {path}")

    n_r = data.shape[1] // 2
    if header_radii and len(header_radii) != n_r:
        header_radii = []

    radii = header_radii if header_radii else [float(i) for i in range(n_r)]

    out: Dict[float, np.ndarray] = {}
    for i, r in enumerate(radii):
        re_col = data[:, 2 * i]
        im_col = data[:, 2 * i + 1]
        out[r] = re_col + 1j * im_col

    idx = np.argsort(t_arr)
    t_arr = t_arr[idx]
    for r in list(out.keys()):
        out[r] = out[r][idx]

    return t_arr, radii, out


# ---------------------------------------------------------------------------
# Strain conversion
# ---------------------------------------------------------------------------
def _psd_psi4_to_strain(freqs: np.ndarray, psd_psi4: np.ndarray) -> np.ndarray:
    """Convert Psi4 PSD to strain PSD: S_h(f) = S_{Psi4}(f) / (2*pi*f)^4."""
    strain_psd = np.zeros_like(psd_psi4)
    nz = freqs > 0
    strain_psd[nz] = psd_psi4[nz] / (2.0 * np.pi * freqs[nz]) ** 4
    return strain_psd


def _scale_to_physical(
    freqs_code: np.ndarray,
    strain_psd_code: np.ndarray,
    mass_msun: float,
    distance_mpc: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert code-unit strain PSD to physical (Hz, 1/Hz) units.

    In code units (G=c=1, M=1) the time unit is T_M = M * M_sun_sec.
    The radius-scaled quantity r*Psi4 is dimensionless; after dividing by
    (2*pi*f_code)^2 the strain PSD in code units has dimensions M^3.
    Physical scaling:
        f_phys = f_code / T_M
        S_h_phys = S_h_code * T_M * (M * M_sun_meter / D)^2
    """
    T_M = mass_msun * M_SUN_SEC
    D_m = distance_mpc * MPC_METER

    f_phys = freqs_code / T_M
    amp_scale = (mass_msun * M_SUN_METER / D_m) ** 2 * T_M
    S_h_phys = strain_psd_code * amp_scale
    return f_phys, S_h_phys


def _aLIGO_noise_psd(freqs_hz: np.ndarray) -> np.ndarray:
    """Advanced LIGO design sensitivity noise PSD S_n(f) [1/Hz].

    Analytic fit from Ajith et al. (2011) / LIGO-T0900288-v3,
    valid for 10 Hz < f < 5000 Hz.  Outside this band we return +inf.
    """
    S = np.full_like(freqs_hz, np.inf)
    f0 = 215.0  # Hz reference frequency
    valid = (freqs_hz >= 10.0) & (freqs_hz <= 5000.0)
    x = freqs_hz[valid] / f0

    S0 = 1.0e-49  # 1/Hz overall scale (approximate design)
    S[valid] = S0 * (
        x ** (-4.14)
        - 5.0 * x ** (-2)
        + 111.0 * (1.0 - x**2 + 0.5 * x**4) / (1.0 + 0.5 * x**2)
    )
    S[valid] = np.abs(S[valid])
    S[valid] = np.where(S[valid] > 0, S[valid], np.inf)
    return S


def _compute_snr(
    freqs_hz: np.ndarray, strain_psd: np.ndarray, noise_psd: np.ndarray
) -> float:
    """Optimal matched-filter SNR^2 = 4 * int |h(f)|^2 / S_n(f) df."""
    valid = np.isfinite(noise_psd) & (noise_psd > 0) & np.isfinite(strain_psd)
    if np.sum(valid) < 2:
        return 0.0
    integrand = strain_psd[valid] / noise_psd[valid]
    snr_sq = 4.0 * np.trapz(integrand, freqs_hz[valid])
    return float(np.sqrt(max(0.0, snr_sq)))


def _frequency_band_label(f_peak_hz: float) -> str:
    if f_peak_hz < 1.0:
        return "sub-Hz band (LISA / pulsar timing)"
    if f_peak_hz < 10.0:
        return "decihertz band (DECIGO / BBO)"
    if f_peak_hz < 2000.0:
        return "most sensitive band of Advanced LIGO"
    if f_peak_hz < 5000.0:
        return "high-frequency band of Advanced LIGO"
    return "above the LIGO band"


# ---------------------------------------------------------------------------
# Propagation speed analysis
# ---------------------------------------------------------------------------
def _find_peak_times(
    t: np.ndarray,
    series: Dict[float, np.ndarray],
    radii: List[float],
    t_skip_frac: float = 0.05,
) -> Dict[float, List[Tuple[float, float]]]:
    """For each radius find (t_peak, amplitude) of dominant and secondary peaks."""
    result: Dict[float, List[Tuple[float, float]]] = {}
    t_skip = t[0] + t_skip_frac * (t[-1] - t[0])

    for R in radii:
        psi4 = series[R]
        envelope = np.abs(psi4)
        mask = t >= t_skip
        env_masked = envelope.copy()
        env_masked[~mask] = 0.0

        peaks_list: List[Tuple[float, float]] = []

        prominence = 0.1 * np.max(env_masked[mask]) if np.any(mask) else 0.0
        idxs, props = find_peaks(env_masked, prominence=max(prominence, 1e-30))
        if len(idxs) > 0:
            order = np.argsort(-env_masked[idxs])
            for idx in idxs[order[:5]]:
                peaks_list.append((float(t[idx]), float(envelope[idx])))
        else:
            i_max = np.argmax(env_masked)
            peaks_list.append((float(t[i_max]), float(envelope[i_max])))

        result[R] = peaks_list
    return result


def _compute_propagation_speeds(
    radii: List[float], peak_data: Dict[float, List[Tuple[float, float]]]
) -> List[Tuple[float, float, float]]:
    """Return list of (R1, R2, speed) using dominant peak at each radius."""
    speeds = []
    for i in range(len(radii) - 1):
        R1, R2 = radii[i], radii[i + 1]
        t1 = peak_data[R1][0][0]
        t2 = peak_data[R2][0][0]
        dt = t2 - t1
        if abs(dt) > 1e-15:
            v = (R2 - R1) / dt
        else:
            v = np.inf
        speeds.append((R1, R2, v))
    return speeds


def _classify_signal(speed: float, tol: float = 0.15) -> str:
    if not np.isfinite(speed):
        return "instantaneous (numerical artifact)"
    if abs(speed - 1.0) < tol:
        return "consistent with physical GW (v ~ c)"
    if speed > 1.0:
        return "superluminal -- likely CCZ4 constraint-damping mode"
    return "subluminal -- likely gauge/constraint mode"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    default_data = _default_data_dir() / "small_data" / "psi4_mode_l2m0.dat"
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Plot extracted Psi4 time-series with PSD, strain, and propagation speed."
    )
    parser.add_argument(
        "input", nargs="?", default=str(default_data),
        help=f"Input .dat (default: {default_data})",
    )
    parser.add_argument("--out", default=str(script_dir))
    parser.add_argument("--name", default=None)
    parser.add_argument("--radii", type=float, nargs="*", default=None)
    parser.add_argument(
        "--time-axis", choices=["simulation", "retarded"], default="simulation",
    )
    parser.add_argument("--plot-psd", action="store_true")
    parser.add_argument("--t-min", type=float, default=None)
    parser.add_argument("--t-max", type=float, default=None)
    parser.add_argument("--psd-smooth-window", type=int, default=21)
    parser.add_argument("--psd-smooth-polyorder", type=int, default=5)
    parser.add_argument("--psd-hide-raw", action="store_true")

    parser.add_argument("--strain", action="store_true", help="Enable strain PSD conversion and LIGO overlay")
    parser.add_argument("--mass-msun", type=float, default=30.0, help="Total mass in solar masses (default: 30)")
    parser.add_argument("--distance-mpc", type=float, default=10.0, help="Luminosity distance in Mpc (default: 10)")
    parser.add_argument("--propagation-speed", action="store_true", help="Measure propagation speed across extraction radii")

    args = parser.parse_args()

    in_path = Path(args.input).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input not found: {in_path}")

    t, radii_all, series = load_extracted(in_path)

    radii = radii_all
    if args.radii:
        radii = [r for r in args.radii if r in series]
        if not radii:
            raise SystemExit(f"No requested radii found in file. Available: {radii_all}")

    # ------------------------------------------------------------------
    # Determine subplot layout
    # ------------------------------------------------------------------
    n_panels = 1  # always: time-domain waveform
    if args.plot_psd:
        n_panels += 1
    if args.strain:
        n_panels += 2  # code-unit strain PSD + physical h_char with LIGO
    if args.propagation_speed and len(radii) >= 2:
        n_panels += 1

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
        "grid.alpha": 0.5,
        "legend.fontsize": 12,
    })

    fig, axes_arr = plt.subplots(n_panels, 1, figsize=(10, 4.0 * n_panels))
    if n_panels == 1:
        axes_arr = [axes_arr]
    else:
        axes_arr = list(axes_arr)

    ax_idx = 0

    # ------------------------------------------------------------------
    # Panel: time-domain waveform
    # ------------------------------------------------------------------
    ax_wave = axes_arr[ax_idx]; ax_idx += 1
    linestyles = ["-", "--", "-.", ":"]

    dt = np.median(np.diff(t)) if t.size >= 2 else np.nan
    fs = (1.0 / dt) if (np.isfinite(dt) and dt > 0) else None

    stored_psd: Dict[float, Tuple[np.ndarray, np.ndarray]] = {}

    for i, R in enumerate(radii):
        color = "black"
        ls = linestyles[i % len(linestyles)]
        psi4 = series[R]
        x = t if args.time_axis == "simulation" else (t - float(R))
        y = np.real(psi4)

        m = np.ones_like(x, dtype=bool)
        if args.t_min is not None:
            m &= x >= args.t_min
        if args.t_max is not None:
            m &= x <= args.t_max

        ax_wave.plot(x[m], y[m], color=color, linestyle=ls, linewidth=1.2, label=rf"$R={R:g}$")

        if fs is not None and psi4.size >= 8:
            freqs, psd = welch(np.real(psi4), fs, nperseg=min(128, max(8, psi4.size // 2)))
            stored_psd[R] = (freqs, psd)

    ax_wave.set_xlabel(r"$t$" if args.time_axis == "simulation" else r"$t - R_{\mathrm{ext}}$")
    ax_wave.set_ylabel(r"$r\,\mathrm{Re}\!\left(\Psi_4^{2,0}\right)$")
    ax_wave.set_title(r"Radius-scaled waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$")
    ax_wave.legend(loc="upper right", frameon=True, framealpha=0.9)
    ax_wave.grid(True, which="both", ls="--", alpha=0.6)
    ax_wave.tick_params(axis="both", which="major", direction="in", top=True, right=True)

    # ------------------------------------------------------------------
    # Panel: Psi4 PSD (code units)
    # ------------------------------------------------------------------
    if args.plot_psd:
        ax_psd = axes_arr[ax_idx]; ax_idx += 1
        for i, R in enumerate(radii):
            if R not in stored_psd:
                continue
            freqs, psd = stored_psd[R]
            psd_s = _smooth_psd(psd, window=args.psd_smooth_window, polyorder=args.psd_smooth_polyorder)
            ls = linestyles[i % len(linestyles)]
            ax_psd.semilogy(freqs, psd_s, color="black", linestyle=ls, linewidth=1.2, label=rf"$R={R:g}$")
        ax_psd.set_xlabel(r"$f\,(M^{-1})$")
        ax_psd.set_ylabel(r"$\mathrm{PSD}\left[r\,\mathrm{Re}(\Psi_4^{2,0})\right]$")
        ax_psd.set_title(r"Power Spectral Density of $\Psi_4$")
        ax_psd.legend(loc="upper right", frameon=True, framealpha=0.9)
        ax_psd.grid(True, which="major", ls="--", alpha=0.6)
        ax_psd.tick_params(axis="both", which="major", direction="in", top=True, right=True)

    # ------------------------------------------------------------------
    # Strain panels
    # ------------------------------------------------------------------
    if args.strain:
        # Use the outermost radius for the strain calculation (cleanest signal)
        R_strain = radii[-1]
        if R_strain in stored_psd:
            freqs, psd_psi4 = stored_psd[R_strain]
            psd_psi4_s = _smooth_psd(psd_psi4, window=args.psd_smooth_window, polyorder=args.psd_smooth_polyorder)
            strain_psd_code = _psd_psi4_to_strain(freqs, psd_psi4_s)

            # -- Code-unit strain PSD --
            ax_strain_code = axes_arr[ax_idx]; ax_idx += 1
            nz = freqs > 0
            ax_strain_code.loglog(freqs[nz], strain_psd_code[nz], color="black", linewidth=1.2)
            ax_strain_code.set_xlabel(r"$f\,(M^{-1})$")
            ax_strain_code.set_ylabel(r"$S_h(f)\,(M^3)$")
            ax_strain_code.set_title(rf"Strain PSD (code units, $R={R_strain:g}$)")
            ax_strain_code.grid(True, which="major", ls="--", alpha=0.6)
            ax_strain_code.tick_params(axis="both", which="major", direction="in", top=True, right=True)

            # -- Physical characteristic strain with LIGO curve --
            ax_phys = axes_arr[ax_idx]; ax_idx += 1
            f_phys, S_h_phys = _scale_to_physical(
                freqs, strain_psd_code, args.mass_msun, args.distance_mpc,
            )
            valid_phys = (f_phys > 0) & np.isfinite(S_h_phys) & (S_h_phys > 0)
            h_char = np.zeros_like(f_phys)
            h_char[valid_phys] = np.sqrt(f_phys[valid_phys] * S_h_phys[valid_phys])

            ax_phys.loglog(
                f_phys[valid_phys], h_char[valid_phys],
                color="black", linewidth=1.5,
                label=rf"Signal ($M={args.mass_msun:g}\,M_\odot$, $D={args.distance_mpc:g}$ Mpc)",
            )

            f_ligo = np.geomspace(10.0, 5000.0, 500)
            S_n = _aLIGO_noise_psd(f_ligo)
            h_n = np.sqrt(f_ligo * S_n)
            ax_phys.loglog(f_ligo, h_n, color="gray", linewidth=1.0, linestyle="--", label="Advanced LIGO design")

            ax_phys.set_xlabel(r"$f$ (Hz)")
            ax_phys.set_ylabel(r"Characteristic strain $h_{\mathrm{char}}(f)$")
            ax_phys.set_title(r"Strain vs.\ Advanced LIGO sensitivity")
            ax_phys.legend(loc="upper right", frameon=True, framealpha=0.9)
            ax_phys.grid(True, which="major", ls="--", alpha=0.6)
            ax_phys.tick_params(axis="both", which="major", direction="in", top=True, right=True)

            # -- SNR calculation and console output --
            S_n_signal = _aLIGO_noise_psd(f_phys)
            snr = _compute_snr(f_phys, S_h_phys, S_n_signal)

            peak_idx = np.argmax(h_char[valid_phys]) if np.any(valid_phys) else 0
            f_peak_code = freqs[valid_phys][peak_idx] if np.any(valid_phys) else 0.0
            f_peak_hz = f_phys[valid_phys][peak_idx] if np.any(valid_phys) else 0.0
            band = _frequency_band_label(f_peak_hz)

            print("\n=== Strain Analysis ===")
            print(f"  Extraction radius: R = {R_strain:g}")
            print(f"  Mass: M = {args.mass_msun:g} M_sun")
            print(f"  Distance: D = {args.distance_mpc:g} Mpc")
            print(f"  Peak strain frequency (code): f = {f_peak_code:.4g} M^-1")
            print(f"  Peak strain frequency (physical): f = {f_peak_hz:.1f} Hz")
            print(f"  Optimal SNR: {snr:.2f}")
            print(f"  -> The burst peaks at {f_peak_hz:.0f} Hz (M={args.mass_msun:g} M_sun),")
            print(f"     placing it in the {band}.")
            print(f"     SNR = {snr:.2f} at D = {args.distance_mpc:g} Mpc.")
        else:
            print("WARNING: No PSD data available for strain conversion (need fs > 0 and >= 8 samples).")
            ax_idx += 2  # skip placeholder panels

    # ------------------------------------------------------------------
    # Panel: propagation speed
    # ------------------------------------------------------------------
    if args.propagation_speed and len(radii) >= 2:
        ax_prop = axes_arr[ax_idx]; ax_idx += 1
        peak_data = _find_peak_times(t, series, radii)
        speeds = _compute_propagation_speeds(radii, peak_data)

        for i, R in enumerate(radii):
            psi4 = series[R]
            retarded = t - float(R)
            ls = linestyles[i % len(linestyles)]
            ax_prop.plot(
                retarded, np.abs(psi4), color="black", linestyle=ls,
                linewidth=1.0, label=rf"$R={R:g}$",
            )
            if peak_data[R]:
                t_pk, amp_pk = peak_data[R][0]
                ax_prop.plot(t_pk - float(R), amp_pk, "o", color="red", markersize=6)

        y_text = 0.92
        for R1, R2, v in speeds:
            classification = _classify_signal(v)
            label_text = rf"$R={R1:g}\to{R2:g}$: $v={v:.3f}\,c$ ({classification})"
            ax_prop.text(
                0.02, y_text, label_text, transform=ax_prop.transAxes,
                fontsize=10, verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7),
            )
            y_text -= 0.08

        ax_prop.set_xlabel(r"Retarded time $t - R_{\mathrm{ext}}$")
        ax_prop.set_ylabel(r"$|r\,\Psi_4^{2,0}|$")
        ax_prop.set_title("Propagation speed analysis")
        ax_prop.legend(loc="upper right", frameon=True, framealpha=0.9)
        ax_prop.grid(True, which="both", ls="--", alpha=0.6)
        ax_prop.tick_params(axis="both", which="major", direction="in", top=True, right=True)

        print("\n=== Propagation Speed Analysis ===")
        for R1, R2, v in speeds:
            classification = _classify_signal(v)
            print(f"  R={R1:g} -> R={R2:g}: v = {v:.4f} c  [{classification}]")

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_".join([f"R{r:g}" for r in radii])
    out_name = args.name if args.name else f"psi4_extracted_{suffix}.eps"
    out_path = out_dir / out_name
    plt.tight_layout()
    plt.savefig(out_path.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.savefig(out_path.with_suffix(".eps"), dpi=600, bbox_inches="tight")
    plt.savefig(out_path.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    print(f"\nSaved to {out_path.with_suffix('.png')}, .eps, and .pdf")


if __name__ == "__main__":
    main()
