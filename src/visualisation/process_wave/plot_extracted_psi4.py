#!/usr/bin/env python3
"""
Plot extracted Psi4 time-series produced by consume_plotfiles.py.

Features:
  - Time-domain waveform r*Re(Psi4)
  - Tukey-windowed FFT energy spectral density of Psi4 (both polarisations)
  - Strain PSD conversion h_tilde(f) = Psi4_tilde(f)/(2*pi*f)^2 with high-pass cutoff
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
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter, find_peaks
from scipy.signal.windows import tukey

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
# Burst energy spectral density (replaces Welch for short transients)
# ---------------------------------------------------------------------------
def _burst_esd(
    psi4_complex: np.ndarray, fs: float, tukey_alpha: float = 0.25
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the one-sided energy spectral density of r*Psi4 for a burst.

    Uses a Tukey window (to taper edges) and a straight FFT rather than
    Welch's method, which is designed for stationary noise and destroys
    frequency resolution of short transients.  Both the + and x
    polarizations are included: ESD = |FFT(Re)|^2 + |FFT(Im)|^2.
    The result is normalised to a one-sided PSD-like quantity (units of
    amplitude^2 / frequency) so downstream code is unchanged.
    """
    N = len(psi4_complex)
    win = tukey(N, alpha=tukey_alpha)
    dt = 1.0 / fs
    T = N * dt

    re_part = np.real(psi4_complex) * win
    im_part = np.imag(psi4_complex) * win

    re_fft = np.fft.rfft(re_part)
    im_fft = np.fft.rfft(im_part)
    freqs = np.fft.rfftfreq(N, d=dt)

    # |FFT|^2 for both polarisations, normalised to one-sided PSD units
    norm = dt**2 / T
    esd = (np.abs(re_fft) ** 2 + np.abs(im_fft) ** 2) * norm
    esd[1:-1] *= 2.0  # one-sided doubling (exclude DC and Nyquist)

    return freqs, esd


# ---------------------------------------------------------------------------
# Strain conversion
# ---------------------------------------------------------------------------
def _psd_psi4_to_strain(
    freqs: np.ndarray,
    psd_psi4: np.ndarray,
    f_low_frac: float = 0.05,
) -> np.ndarray:
    """Convert Psi4 PSD to strain PSD: S_h(f) = S_{Psi4}(f) / (2*pi*f)^4.

    A 4th-order Butterworth-style high-pass roll-off is applied below
    ``f_low = f_low_frac * f_max`` to suppress the unphysical divergence
    from dividing numerical noise by f^4 as f -> 0.
    """
    strain_psd = np.zeros_like(psd_psi4)
    nz = freqs > 0
    f_max = freqs[nz].max() if np.any(nz) else 1.0
    f_low = f_low_frac * f_max

    omega4 = (2.0 * np.pi * freqs[nz]) ** 4
    strain_psd[nz] = psd_psi4[nz] / omega4

    # Butterworth high-pass suppression (4th order) for f < f_low
    hp = 1.0 / (1.0 + (f_low / freqs[nz]) ** 8)
    strain_psd[nz] *= hp

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
    snr_sq = 4.0 * np.trapezoid(integrand, freqs_hz[valid])
    return float(np.sqrt(max(0.0, snr_sq)))


def _frequency_band_label(f_peak_hz: float) -> str:
    if f_peak_hz < 1.0:
        return "in the sub-Hz band (LISA / pulsar timing)"
    if f_peak_hz < 10.0:
        return "in the decihertz band (DECIGO / BBO)"
    if f_peak_hz < 2000.0:
        return "in the most sensitive band of Advanced LIGO"
    if f_peak_hz < 5000.0:
        return "in the high-frequency band of Advanced LIGO"
    return "above the LIGO band"


# ---------------------------------------------------------------------------
# QNM ringdown fitting
# ---------------------------------------------------------------------------
def _damped_sinusoid(t, A, tau, f0, phi):
    return A * np.exp(-t / tau) * np.sin(2.0 * np.pi * f0 * t + phi)


def _fit_qnm(
    t: np.ndarray,
    psi4_complex: np.ndarray,
    R: float,
    tail_fraction: float = 0.4,
) -> dict | None:
    """Fit A*exp(-t/tau)*sin(2*pi*f*t + phi) to the late-time ringdown.

    Works on the retarded-time Re(r*Psi4) waveform at a single extraction
    radius.  Returns dict with fit parameters, or None if the fit fails.
    """
    t_ret = t - R
    y = np.real(psi4_complex)

    envelope = np.abs(psi4_complex)
    i_peak = np.argmax(envelope)
    if i_peak >= len(t) - 10:
        return None

    t_tail = t_ret[i_peak:]
    y_tail = y[i_peak:]
    t_tail = t_tail - t_tail[0]

    n_start = max(1, int(tail_fraction * len(t_tail)))
    t_fit = t_tail[n_start:]
    y_fit = y_tail[n_start:]
    if len(t_fit) < 10:
        return None

    env_fit = np.abs(y_fit)
    A0 = float(np.max(env_fit)) if np.max(env_fit) > 0 else 1e-6
    tau0 = float(t_fit[-1] - t_fit[0]) / 2.0

    peaks_idx, _ = find_peaks(np.abs(y_fit), prominence=0.1 * A0)
    if len(peaks_idx) >= 2:
        dt_peaks = np.median(np.diff(t_fit[peaks_idx]))
        f0_guess = 1.0 / (2.0 * dt_peaks) if dt_peaks > 0 else 1.0
    else:
        f0_guess = 2.0

    try:
        popt, pcov = curve_fit(
            _damped_sinusoid, t_fit, y_fit,
            p0=[A0, tau0, f0_guess, 0.0],
            bounds=([0, 1e-6, 1e-3, -2 * np.pi], [np.inf, np.inf, np.inf, 2 * np.pi]),
            maxfev=10000,
        )
        A, tau, f_qnm, phi = popt
        perr = np.sqrt(np.diag(pcov))
        return {
            "A": A, "tau": tau, "f_qnm": f_qnm, "phi": phi,
            "A_err": perr[0], "tau_err": perr[1], "f_qnm_err": perr[2],
            "t_fit": t_fit + t_ret[i_peak] + t_tail[n_start],
            "y_fit": _damped_sinusoid(t_fit, *popt),
            "t_ret_offset": t_ret[i_peak] + t_tail[n_start],
        }
    except Exception:
        return None


def _schwarzschild_qnm_l2(M_bh: float) -> Tuple[float, float]:
    """Fundamental l=2 QNM of Schwarzschild: (f_code, tau_code) for mass M_bh."""
    f_code = 0.08896 / M_bh
    tau_code = 1.0 / (0.08896 * 2 * np.pi / (2 * 0.0950 / M_bh))
    omega_R = 2 * np.pi * 0.08896 / M_bh
    omega_I = 0.0950 / M_bh
    return omega_R / (2 * np.pi), 1.0 / omega_I


# ---------------------------------------------------------------------------
# Total radiated energy
# ---------------------------------------------------------------------------
def _compute_radiated_energy(
    t: np.ndarray, psi4_complex: np.ndarray
) -> float:
    """E_rad = (1/16*pi) * int |r*Psi4|^2 dt  (code units, G=c=1)."""
    integrand = np.abs(psi4_complex) ** 2
    return float(np.trapezoid(integrand, t) / (16.0 * np.pi))


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
    """Classify propagation speed (used only for console diagnostics)."""
    if not np.isfinite(speed):
        return "instantaneous (numerical artifact)"
    if abs(speed - 1.0) < tol:
        return "consistent with physical GW (v ~ c)"
    if speed > 1.0:
        return "superluminal -- likely CCZ4 constraint-damping mode"
    return "subluminal -- likely gauge/constraint mode"


# ---------------------------------------------------------------------------
# Shared panel-drawing helpers
# ---------------------------------------------------------------------------
def _draw_waveform(
    ax, t, radii, series, linestyles, time_axis, t_min=None, t_max=None
):
    for i, R in enumerate(radii):
        ls = linestyles[i % len(linestyles)]
        x = t if time_axis == "simulation" else (t - float(R))
        y = np.real(series[R])
        m = np.ones_like(x, dtype=bool)
        if t_min is not None:
            m &= x >= t_min
        if t_max is not None:
            m &= x <= t_max
        ax.plot(x[m], y[m], color="black", linestyle=ls, linewidth=1.0, label=rf"$R={R:g}$")
    ax.set_xlabel(r"$t$" if time_axis == "simulation" else r"$t - R_{\mathrm{ext}}$")
    ax.set_ylabel(r"$r\,\mathrm{Re}\!\left(\Psi_4^{2,0}\right)$")
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=9)
    ax.grid(True, which="both", ls="--", alpha=0.6)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)


def _draw_psd(ax, radii, stored_psd, linestyles, smooth_w, smooth_p):
    for i, R in enumerate(radii):
        if R not in stored_psd:
            continue
        freqs, psd = stored_psd[R]
        psd_s = _smooth_psd(psd, window=smooth_w, polyorder=smooth_p)
        ls = linestyles[i % len(linestyles)]
        ax.semilogy(freqs, psd_s, color="black", linestyle=ls, linewidth=1.0, label=rf"$R={R:g}$")
    ax.set_xlabel(r"$f\,(M^{-1})$")
    ax.set_ylabel(r"$\mathrm{ESD}\left[r\,\Psi_4^{2,0}\right]$")
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=9)
    ax.grid(True, which="major", ls="--", alpha=0.6)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)


def _draw_strain_code(ax, freqs, strain_psd_code, R_strain):
    nz = freqs > 0
    ax.loglog(freqs[nz], strain_psd_code[nz], color="black", linewidth=1.0)
    ax.set_xlabel(r"$f\,(M^{-1})$")
    ax.set_ylabel(r"$S_h(f)\,(M^3)$")
    ax.grid(True, which="major", ls="--", alpha=0.6)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)


def _draw_strain_ligo(ax, freqs, strain_psd_code, mass_msun, distance_mpc, R_ext=None):
    f_phys, S_h_phys = _scale_to_physical(freqs, strain_psd_code, mass_msun, distance_mpc)
    valid_phys = (f_phys > 0) & np.isfinite(S_h_phys) & (S_h_phys > 0)
    h_char = np.zeros_like(f_phys)
    h_char[valid_phys] = np.sqrt(f_phys[valid_phys] * S_h_phys[valid_phys])

    ax.loglog(
        f_phys[valid_phys], h_char[valid_phys],
        color="black", linewidth=1.2,
        label=rf"Signal ($M\!=\!{mass_msun:g}\,M_\odot$, $D\!=\!{distance_mpc:g}$ Mpc)",
    )
    f_ligo = np.geomspace(10.0, 5000.0, 500)
    S_n = _aLIGO_noise_psd(f_ligo)
    h_n = np.sqrt(f_ligo * S_n)
    ax.loglog(f_ligo, h_n, color="gray", linewidth=1.0, linestyle="--", label="Advanced LIGO design")
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$h_{\mathrm{char}}(f)$")
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8)
    ax.grid(True, which="major", ls="--", alpha=0.6)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)

    S_n_signal = _aLIGO_noise_psd(f_phys)
    snr = _compute_snr(f_phys, S_h_phys, S_n_signal)
    peak_idx = np.argmax(h_char[valid_phys]) if np.any(valid_phys) else 0
    f_peak_code = freqs[valid_phys][peak_idx] if np.any(valid_phys) else 0.0
    f_peak_hz = f_phys[valid_phys][peak_idx] if np.any(valid_phys) else 0.0
    band = _frequency_band_label(f_peak_hz)

    snr_threshold = 8.0
    if snr > 0:
        D_horizon = distance_mpc * (snr / snr_threshold)
    else:
        D_horizon = 0.0

    print("\n=== Strain Analysis ===")
    print(f"  Extraction radius: R = {R_ext:g}" if R_ext is not None else "  Extraction radius: (unknown)")
    print(f"  Mass: M = {mass_msun:g} M_sun")
    print(f"  Distance: D = {distance_mpc:g} Mpc")
    print(f"  Peak strain frequency (code): f = {f_peak_code:.4g} M^-1")
    print(f"  Peak strain frequency (physical): f = {f_peak_hz:.1f} Hz")
    print(f"  Optimal SNR: {snr:.2f}")
    print(f"  -> The burst peaks at {f_peak_hz:.0f} Hz,")
    print(f"     placing it {band}.")
    print(f"     SNR = {snr:.2f} at D = {distance_mpc:g} Mpc.")
    if D_horizon > 0:
        print(f"  -> Horizon distance (SNR >= {snr_threshold:.0f}): D_horizon = {D_horizon:.1f} Mpc")
    else:
        print(f"  -> Signal below LIGO noise floor; no detectable horizon distance.")


def _draw_propagation(ax, t, radii, series, linestyles):
    peak_data = _find_peak_times(t, series, radii)
    speeds = _compute_propagation_speeds(radii, peak_data)

    for i, R in enumerate(radii):
        psi4 = series[R]
        retarded = t - float(R)
        ls = linestyles[i % len(linestyles)]
        ax.plot(retarded, np.abs(psi4), color="black", linestyle=ls, linewidth=0.8, label=rf"$R={R:g}$")
        if peak_data[R]:
            t_pk, amp_pk = peak_data[R][0]
            ax.plot(t_pk - float(R), amp_pk, "o", color="red", markersize=5)

    y_text = 0.92
    for R1, R2, v in speeds:
        label_text = rf"$R={R1:g}\!\to\!{R2:g}$: $v={v:.3f}\,c$"
        ax.text(
            0.02, y_text, label_text, transform=ax.transAxes,
            fontsize=8, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="wheat", alpha=0.7),
        )
        y_text -= 0.08

    ax.set_xlabel(r"Retarded time $t - R_{\mathrm{ext}}$")
    ax.set_ylabel(r"$|r\,\Psi_4^{2,0}|$")
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=9)
    ax.grid(True, which="both", ls="--", alpha=0.6)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)

    print("\n=== Propagation Speed Analysis ===")
    for R1, R2, v in speeds:
        classification = _classify_signal(v)
        print(f"  R={R1:g} -> R={R2:g}: v = {v:.4f} c  [{classification}]")


# ---------------------------------------------------------------------------
# Combined 3x2 publication figure
# ---------------------------------------------------------------------------
def _plot_combined(t, radii, series, stored_psd, args, linestyles, fs):
    import string

    fig, axes = plt.subplots(3, 2, figsize=(14, 14))

    # (a) Waveform — simulation time
    _draw_waveform(axes[0, 0], t, radii, series, linestyles, "simulation",
                   t_min=args.t_min, t_max=args.t_max)

    # (b) Waveform — retarded time + QNM overlay
    _draw_waveform(axes[0, 1], t, radii, series, linestyles, "retarded",
                   t_min=args.t_min, t_max=args.t_max)

    R_inner = radii[0]
    qnm_result = _fit_qnm(t, series[R_inner], R_inner)
    if qnm_result is not None:
        t_overlay = qnm_result["t_fit"]
        y_overlay = qnm_result["y_fit"]
        t_ret_start = qnm_result["t_ret_offset"]
        t_plot = np.linspace(t_overlay[0], t_overlay[-1], 300)
        y_plot = _damped_sinusoid(
            t_plot - t_overlay[0],
            qnm_result["A"], qnm_result["tau"],
            qnm_result["f_qnm"], qnm_result["phi"],
        )
        axes[0, 1].plot(
            t_plot, y_plot, color="red", linewidth=1.0, linestyle="--", alpha=0.85,
            label=rf"QNM fit: $f={qnm_result['f_qnm']:.3f}$, $\tau={qnm_result['tau']:.2f}$",
        )
        axes[0, 1].legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8)

    # Radiated energy and QNM diagnostics
    E_rad = _compute_radiated_energy(t, series[R_inner])
    print(f"\n=== Radiated Energy (R={R_inner:g}) ===")
    print(f"  E_rad = {E_rad:.6e} M  (code units)")
    if E_rad > 0:
        print(f"  E_rad / M_initial ~ {E_rad:.4e}")
    if qnm_result is not None:
        f_qnm = qnm_result["f_qnm"]
        tau_qnm = qnm_result["tau"]
        print(f"\n=== QNM Ringdown Fit (R={R_inner:g}) ===")
        print(f"  f_QNM = {f_qnm:.4f} M^-1  (tau = {tau_qnm:.4f} M)")
        print(f"  Quality factor Q = pi * f * tau = {np.pi * f_qnm * tau_qnm:.2f}")
        M_from_f = 0.08896 / f_qnm
        M_from_tau = 0.0950 / (1.0 / tau_qnm) if tau_qnm > 0 else 0
        f_schw, tau_schw = _schwarzschild_qnm_l2(M_from_f)
        print(f"  -> If Schwarzschild l=2 QNM:  M_BH = {M_from_f:.4f} (from f)")
        print(f"     Predicted tau = {1.0 / (0.0950 / M_from_f):.4f} M  (observed: {tau_qnm:.4f})")
        T_M = args.mass_msun * M_SUN_SEC
        f_phys_qnm = f_qnm / T_M
        print(f"  -> At M = {args.mass_msun:g} M_sun: f_QNM = {f_phys_qnm:.0f} Hz")

    # (c) PSD of Ψ4
    _draw_psd(axes[1, 0], radii, stored_psd, linestyles,
              args.psd_smooth_window, args.psd_smooth_polyorder)

    # (d) Propagation speed analysis
    if len(radii) >= 2:
        _draw_propagation(axes[1, 1], t, radii, series, linestyles)
    else:
        axes[1, 1].text(0.5, 0.5, "Need $\\geq 2$ radii", transform=axes[1, 1].transAxes,
                        ha="center", va="center", fontsize=12)

    # (e) Strain PSD (code units) and (f) Strain vs LIGO
    R_strain = radii[-1]
    if R_strain in stored_psd:
        freqs, psd_psi4 = stored_psd[R_strain]
        psd_psi4_s = _smooth_psd(psd_psi4, window=args.psd_smooth_window, polyorder=args.psd_smooth_polyorder)
        strain_psd_code = _psd_psi4_to_strain(freqs, psd_psi4_s)

        _draw_strain_code(axes[2, 0], freqs, strain_psd_code, R_strain)
        _draw_strain_ligo(axes[2, 1], freqs, strain_psd_code, args.mass_msun, args.distance_mpc, R_ext=R_strain)
    else:
        for ax in [axes[2, 0], axes[2, 1]]:
            ax.text(0.5, 0.5, "No PSD data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=12)

    titles = [
        r"Waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$",
        r"Waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$ (retarded) + QNM fit",
        r"Energy Spectral Density of $\Psi_4$",
        "Propagation speed analysis",
        rf"Strain PSD (code units, $R={R_strain:g}$)",
        r"Strain vs.\ Advanced LIGO sensitivity",
    ]
    for ax, letter, title in zip(axes.flatten(), string.ascii_lowercase, titles):
        ax.set_title(f"({letter}) {title}")

    fig.tight_layout()
    _save_figure(fig, args, radii)


# ---------------------------------------------------------------------------
# Legacy stacked (single-column) figure
# ---------------------------------------------------------------------------
def _plot_stacked(t, radii, series, stored_psd, args, linestyles, fs):

    n_panels = 1
    if args.plot_psd:
        n_panels += 1
    if args.strain:
        n_panels += 2
    if args.propagation_speed and len(radii) >= 2:
        n_panels += 1

    fig, axes_arr = plt.subplots(n_panels, 1, figsize=(10, 4.0 * n_panels))
    if n_panels == 1:
        axes_arr = [axes_arr]
    else:
        axes_arr = list(axes_arr)

    ax_idx = 0

    # Waveform
    ax_wave = axes_arr[ax_idx]; ax_idx += 1
    _draw_waveform(ax_wave, t, radii, series, linestyles, args.time_axis,
                   t_min=args.t_min, t_max=args.t_max)
    label = (r"Radius-scaled waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$"
             if args.time_axis == "simulation"
             else r"Radius-scaled waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$ (retarded)")
    ax_wave.set_title(label)

    # PSD
    if args.plot_psd:
        ax_psd = axes_arr[ax_idx]; ax_idx += 1
        _draw_psd(ax_psd, radii, stored_psd, linestyles,
                  args.psd_smooth_window, args.psd_smooth_polyorder)
        ax_psd.set_title(r"Power Spectral Density of $\Psi_4$")

    # Strain
    if args.strain:
        R_strain = radii[-1]
        if R_strain in stored_psd:
            freqs, psd_psi4 = stored_psd[R_strain]
            psd_psi4_s = _smooth_psd(psd_psi4, window=args.psd_smooth_window,
                                     polyorder=args.psd_smooth_polyorder)
            strain_psd_code = _psd_psi4_to_strain(freqs, psd_psi4_s)

            ax_sc = axes_arr[ax_idx]; ax_idx += 1
            _draw_strain_code(ax_sc, freqs, strain_psd_code, R_strain)
            ax_sc.set_title(rf"Strain PSD (code units, $R={R_strain:g}$)")

            ax_ph = axes_arr[ax_idx]; ax_idx += 1
            _draw_strain_ligo(ax_ph, freqs, strain_psd_code, args.mass_msun, args.distance_mpc, R_ext=R_strain)
            ax_ph.set_title(r"Strain vs.\ Advanced LIGO sensitivity")
        else:
            print("WARNING: No PSD data for strain conversion.")
            ax_idx += 2

    # Propagation speed
    if args.propagation_speed and len(radii) >= 2:
        ax_prop = axes_arr[ax_idx]; ax_idx += 1
        _draw_propagation(ax_prop, t, radii, series, linestyles)
        ax_prop.set_title("Propagation speed analysis")

    fig.tight_layout()
    _save_figure(fig, args, radii)


# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------
def _save_figure(fig, args, radii):
    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_".join([f"R{r:g}" for r in radii])
    out_name = args.name if args.name else f"psi4_extracted_{suffix}.eps"
    out_path = out_dir / Path(out_name)
    fig.savefig(out_path.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".eps"), dpi=600, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    print(f"\nSaved to {out_path.with_suffix('.png')}, .eps, and .pdf")
    plt.close(fig)


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
    parser.add_argument(
        "--combined", action="store_true",
        help="Produce a single 3x2 publication figure with all analysis panels."
    )

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
    # Common style
    # ------------------------------------------------------------------
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
        "legend.fontsize": 11,
    })

    linestyles = ["-", "--", "-.", ":"]

    dt = np.median(np.diff(t)) if t.size >= 2 else np.nan
    fs = (1.0 / dt) if (np.isfinite(dt) and dt > 0) else None

    stored_psd: Dict[float, Tuple[np.ndarray, np.ndarray]] = {}
    for i, R in enumerate(radii):
        psi4 = series[R]
        if fs is not None and psi4.size >= 8:
            freqs, psd = _burst_esd(psi4, fs)
            stored_psd[R] = (freqs, psd)

    if args.combined:
        _plot_combined(t, radii, series, stored_psd, args, linestyles, fs)
    else:
        _plot_stacked(t, radii, series, stored_psd, args, linestyles, fs)


if __name__ == "__main__":
    main()
