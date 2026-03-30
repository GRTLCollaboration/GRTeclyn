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

M_SUN_KG = 1.98892e30
G_SI = 6.67430e-11
C_SI = 2.99792458e8
M_SUN_SEC = G_SI * M_SUN_KG / C_SI**3
M_SUN_METER = G_SI * M_SUN_KG / C_SI**2
MPC_METER = 3.08568e22


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

    norm = dt**2 / T
    esd = (np.abs(re_fft) ** 2 + np.abs(im_fft) ** 2) * norm
    esd[1:-1] *= 2.0  # one-sided doubling (exclude DC and Nyquist)

    return freqs, esd


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
    t_fit_raw = t_tail[n_start:]
    y_fit = y_tail[n_start:]
    if len(t_fit_raw) < 10:
        return None

    t_fit = t_fit_raw - t_fit_raw[0]

    env_fit = np.abs(y_fit)
    A0 = float(np.max(env_fit)) if np.max(env_fit) > 0 else 1e-6
    tau0 = float(t_fit[-1] - t_fit[0]) / 2.0

    peaks_idx, _ = find_peaks(np.abs(y_fit), prominence=0.1 * A0)
    if len(peaks_idx) >= 2:
        dt_peaks = np.median(np.diff(t_fit[peaks_idx]))
        f0_guess = 1.0 / (2.0 * dt_peaks) if dt_peaks > 0 else 1.0
    else:
        f0_guess = 2.0

    t_ret_fit_start = t_ret[i_peak] + t_fit_raw[0]

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
            "t_ret_start": t_ret_fit_start,
            "t_ret_end": t_ret_fit_start + t_fit[-1],
        }
    except Exception:
        return None


_QNM_OMEGA_R_DIMLESS = 0.37367
_QNM_OMEGA_I_DIMLESS = 0.08896


def _schwarzschild_qnm_l2(M_bh: float) -> Tuple[float, float]:
    """Fundamental l=2 QNM of Schwarzschild: (f_code, tau_code) for mass M_bh."""
    f_code = _QNM_OMEGA_R_DIMLESS / (2 * np.pi * M_bh)
    tau_code = M_bh / _QNM_OMEGA_I_DIMLESS
    return f_code, tau_code


def _compute_radiated_energy(
    t: np.ndarray, psi4_complex: np.ndarray
) -> float:
    """E_rad = (1/16*pi) * int |r*Psi4|^2 dt  (code units, G=c=1)."""
    integrand = np.abs(psi4_complex) ** 2
    return float(np.trapezoid(integrand, t) / (16.0 * np.pi))


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
    radii: List[float], peak_data: Dict[float, List[Tuple[float, float]]],
    t: np.ndarray = None, series: Dict[float, np.ndarray] = None,
) -> List[Tuple[float, float, float]]:
    """Return list of (R1, R2, speed).

    Uses wavefront tracking: the dominant peak at the innermost radius
    defines a reference retarded time.  At each subsequent radius, the
    peak whose retarded time is closest to that reference is selected,
    ensuring we track the *same* physical wavefront rather than jumping
    to a different (potentially constraint-dominated) feature.
    """
    if not radii or not peak_data:
        return []

    R_ref = radii[0]
    t_ref_sim = peak_data[R_ref][0][0]
    t_ref_ret = t_ref_sim - R_ref

    matched_sim: Dict[float, float] = {R_ref: t_ref_sim}
    for R in radii[1:]:
        best_t_sim = peak_data[R][0][0]
        best_dt_ret = abs((best_t_sim - R) - t_ref_ret)
        for (t_pk, _) in peak_data[R]:
            dt_ret = abs((t_pk - R) - t_ref_ret)
            if dt_ret < best_dt_ret:
                best_dt_ret = dt_ret
                best_t_sim = t_pk
        matched_sim[R] = best_t_sim

    speeds = []
    for i in range(len(radii) - 1):
        R1, R2 = radii[i], radii[i + 1]
        dt = matched_sim[R2] - matched_sim[R1]
        v = (R2 - R1) / dt if abs(dt) > 1e-15 else np.inf
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


def _draw_psd(ax, radii, stored_psd, linestyles, smooth_w, smooth_p,
              pert_sigma=None, f_max: float | None = None):
    for i, R in enumerate(radii):
        if R not in stored_psd:
            continue
        freqs, psd = stored_psd[R]
        psd_s = _smooth_psd(psd, window=smooth_w, polyorder=smooth_p)
        ls = linestyles[i % len(linestyles)]
        if f_max is not None:
            m = freqs <= float(f_max)
            ax.semilogy(freqs[m], psd_s[m], color="black", linestyle=ls, linewidth=1.0, label=rf"$R={R:g}$")
        else:
            ax.semilogy(freqs, psd_s, color="black", linestyle=ls, linewidth=1.0, label=rf"$R={R:g}$")

    if pert_sigma is not None and pert_sigma > 0:
        any_R = next(iter(stored_psd))
        f_arr, psd_arr = stored_psd[any_R]
        nz = f_arr > 0
        psd_peak = np.max(_smooth_psd(psd_arr, window=smooth_w, polyorder=smooth_p)[nz])
        gauss_psd = np.exp(-2.0 * (np.pi * pert_sigma * f_arr[nz]) ** 2)
        gauss_psd *= psd_peak
        visible = gauss_psd > psd_peak * 1e-8
        if np.any(visible):
            ax.semilogy(f_arr[nz][visible], gauss_psd[visible], color="red",
                        linewidth=1.2, linestyle=":",
                        label=rf"Gaussian pert. ($\sigma_K={pert_sigma:g}$)")
        f_char = 1.0 / (np.pi * pert_sigma * np.sqrt(2))
        ax.axvline(f_char, color="red", linewidth=0.8, linestyle="--", alpha=0.6)

    if f_max is not None:
        ax.set_xlim(0.0, float(f_max))

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


def _draw_strain_ligo(ax, freqs, strain_psd_code, mass_msun, distance_mpc, R_ext=None, ligo_quantity: str = "asd"):
    f_phys, S_h_phys = _scale_to_physical(freqs, strain_psd_code, mass_msun, distance_mpc)
    valid_phys = (f_phys > 0) & np.isfinite(S_h_phys) & (S_h_phys > 0)
    ligo_quantity = (ligo_quantity or "asd").lower().strip()
    if ligo_quantity not in {"asd", "hchar"}:
        ligo_quantity = "asd"

    y_sig = np.zeros_like(f_phys)
    if ligo_quantity == "asd":
        y_sig[valid_phys] = np.sqrt(S_h_phys[valid_phys])
        y_label = r"Strain noise, $1/\sqrt{\mathrm{Hz}}$"
    else:
        y_sig[valid_phys] = np.sqrt(f_phys[valid_phys] * S_h_phys[valid_phys])
        y_label = r"$h_{\mathrm{char}}(f)$"

    ax.loglog(
        f_phys[valid_phys], y_sig[valid_phys],
        color="black", linewidth=1.2,
        label=rf"Signal ($M\!=\!{mass_msun:g}\,M_\odot$, $D\!=\!{distance_mpc:g}$ Mpc)",
    )
    f_ligo = np.geomspace(10.0, 5000.0, 500)
    S_n = _aLIGO_noise_psd(f_ligo)
    if ligo_quantity == "asd":
        y_noise = np.sqrt(S_n)
    else:
        y_noise = np.sqrt(f_ligo * S_n)
    ax.loglog(f_ligo, y_noise, color="gray", linewidth=1.0, linestyle="--", label="Advanced LIGO design")
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(y_label)
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8)
    ax.grid(True, which="major", ls="--", alpha=0.6)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)
    ax.set_xlim(10.0, 5000.0)

    S_n_signal = _aLIGO_noise_psd(f_phys)
    snr = _compute_snr(f_phys, S_h_phys, S_n_signal)
    h_char_diag = np.zeros_like(f_phys)
    h_char_diag[valid_phys] = np.sqrt(f_phys[valid_phys] * S_h_phys[valid_phys])
    peak_idx = np.argmax(h_char_diag[valid_phys]) if np.any(valid_phys) else 0
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
    speeds = _compute_propagation_speeds(radii, peak_data, t=t, series=series)

    R_ref = radii[0]
    t_ref_ret = peak_data[R_ref][0][0] - R_ref
    matched_sim: Dict[float, float] = {}
    for R in radii:
        best_t = peak_data[R][0][0]
        best_dt = abs((best_t - R) - t_ref_ret)
        for (t_pk, _) in peak_data[R]:
            dt_ret = abs((t_pk - R) - t_ref_ret)
            if dt_ret < best_dt:
                best_dt = dt_ret
                best_t = t_pk
        matched_sim[R] = best_t

    for i, R in enumerate(radii):
        psi4 = series[R]
        retarded = t - float(R)
        envelope = np.abs(psi4)
        ls = linestyles[i % len(linestyles)]
        ax.plot(retarded, envelope, color="black", linestyle=ls, linewidth=0.8, label=rf"$R={R:g}$")
        t_matched = matched_sim[R]
        idx = np.argmin(np.abs(t - t_matched))
        ax.plot(t_matched - float(R), envelope[idx], "o", color="red", markersize=5)

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


def _draw_spectrogram(ax, t, psi4, R, args, fs):
    if fs is None or len(t) < 64:
        ax.text(0.5, 0.5, "Insufficient data for spectrogram", transform=ax.transAxes, ha="center")
        return

    import scipy.signal as signal
    from matplotlib.ticker import ScalarFormatter
    
    y_re = np.real(psi4)
    dt = 1.0 / fs

    f_max = 3.0
    if args.esd_fmax and float(args.esd_fmax) < f_max:
        f_max = float(args.esd_fmax)

    f_min = 0.1
    n_freq = 260
    freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)
    w = 8.0

    N_orig = len(y_re)
    pad_len = N_orig // 2
    y_padded = np.pad(y_re, (pad_len, pad_len), 'constant')
    N_padded = len(y_padded)
    
    cwtmatr = np.zeros((len(freqs), N_padded), dtype=complex)
    t_wave = np.arange(-N_padded//2, N_padded//2) * dt
    
    for i, f in enumerate(freqs):
        s = w / (2.0 * np.pi * f)
        wavelet = np.exp(2j * np.pi * f * t_wave) * np.exp(-t_wave**2 / (2.0 * s**2))
        wavelet *= (np.pi * s**2)**(-0.25)
        cwtmatr[i, :] = signal.fftconvolve(y_padded, wavelet, mode='same')
        
    amplitude = np.abs(cwtmatr)
    amplitude = amplitude[:, pad_len:pad_len+N_orig]
    
    t_ret = t - R
    
    if amplitude.size > 0 and np.nanmax(amplitude) > 0:
        amp_max = float(np.nanmax(amplitude))
        amp = amplitude / amp_max
        amp = np.nan_to_num(amp, nan=0.0, posinf=0.0, neginf=0.0)

        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad(cmap(0.0))

        amp_m = np.ma.masked_invalid(amp)

        pcm = ax.pcolormesh(
            t_ret,
            freqs,
            amp_m,
            shading="auto",
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            rasterized=True,
            antialiased=False,
        )

        cbar = plt.colorbar(pcm, ax=ax, pad=0.02)
        cbar.set_label("Normalized amplitude", fontsize=10)
        cbar.ax.tick_params(labelsize=9)
    else:
        ax.text(0.5, 0.5, "Spectrogram power is zero", transform=ax.transAxes, ha="center")
        
    ax.set_yscale("log")
    ax.set_ylim(f_min, f_max)
    ax.set_yticks([0.1, 0.2, 0.5, 1.0, 2.0, 3.0])
    ax.yaxis.set_major_formatter(ScalarFormatter())

    if args.t_min is not None or args.t_max is not None:
        if args.t_min is not None:
            ax.set_xlim(left=args.t_min)
        if args.t_max is not None:
            ax.set_xlim(right=args.t_max)
    else:
        x_left = float(np.nanmin(t_ret))
        x_right = min(15.0, float(np.nanmax(t_ret)))
        ax.set_xlim(x_left, x_right)
        
    ax.set_ylabel(r"$f\,(M^{-1})$")
    ax.set_xlabel(r"Retarded time $t - R_{\mathrm{ext}}$")

    ax.set_axisbelow(False)
    ax.grid(True, which="major", color="white", alpha=0.4, linestyle="-", linewidth=0.8)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True)


def _plot_combined(t, radii, series, stored_psd, args, linestyles, fs):
    import string

    fig, axes = plt.subplots(3, 2, figsize=(14, 14))

    _draw_waveform(axes[0, 0], t, radii, series, linestyles, "simulation",
                   t_min=args.t_min, t_max=args.t_max)
    if args.pert_sigma is not None and args.pert_sigma > 0:
        A0 = getattr(args, "pert_A0", None) or 0.0
        A2 = getattr(args, "pert_A2", None) or 0.0
        sigma = args.pert_sigma
        r_prof = np.linspace(0, 5 * sigma, 300)
        K_prof = (A0 + A2) * np.exp(-r_prof ** 2 / sigma ** 2)
        ax_k = axes[0, 0].twinx()
        ax_k.plot(r_prof, K_prof, color="blue", linewidth=1.0, linestyle="-", alpha=0.5,
                  label=r"$K(r,\theta{=}0)$ at $t{=}0$")
        ax_k.set_ylabel(r"$K(r)$", color="blue", fontsize=12)
        ax_k.tick_params(axis="y", colors="blue", labelsize=10)
        ax_k.legend(loc="center right", frameon=True, framealpha=0.8, fontsize=8)

    _draw_waveform(axes[0, 1], t, radii, series, linestyles, "retarded",
                   t_min=args.t_min, t_max=args.t_max)

    R_inner = radii[0]
    qnm_result = _fit_qnm(t, series[R_inner], R_inner)
    if qnm_result is not None:
        t0 = qnm_result["t_ret_start"]
        t1 = qnm_result["t_ret_end"]
        t_plot = np.linspace(t0, t1, 300)
        y_plot = _damped_sinusoid(
            t_plot - t0,
            qnm_result["A"], qnm_result["tau"],
            qnm_result["f_qnm"], qnm_result["phi"],
        )
        axes[0, 1].plot(
            t_plot, y_plot, color="red", linewidth=1.0, linestyle="--", alpha=0.85,
            label=rf"QNM fit: $f={qnm_result['f_qnm']:.3f}$, $\tau={qnm_result['tau']:.2f}$",
        )
        axes[0, 1].legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8)

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
        M_from_f = _QNM_OMEGA_R_DIMLESS / (2 * np.pi * f_qnm)
        M_from_tau = tau_qnm * _QNM_OMEGA_I_DIMLESS if tau_qnm > 0 else 0
        f_schw, tau_schw = _schwarzschild_qnm_l2(M_from_f)
        print(f"  -> If Schwarzschild l=2 QNM:  M_BH = {M_from_f:.4f} (from f),  {M_from_tau:.4f} (from tau)")
        print(f"     Predicted tau = {tau_schw:.4f} M  (observed: {tau_qnm:.4f})")
        T_M = args.mass_msun * M_SUN_SEC
        f_phys_qnm = f_qnm / T_M
        print(f"  -> At M = {args.mass_msun:g} M_sun: f_QNM = {f_phys_qnm:.0f} Hz")
        if args.pert_sigma is not None and args.pert_sigma > 0:
            f_pert = 1.0 / (np.pi * args.pert_sigma * np.sqrt(2))
            ratio = f_qnm / f_pert
            print(f"  --- Initial K perturbation vs QNM comparison ---")
            print(f"  Gaussian perturbation sigma_K = {args.pert_sigma:g}")
            print(f"  Perturbation char. freq  = {f_pert:.4f} M^-1")
            print(f"  Fitted QNM freq          = {f_qnm:.4f} M^-1")
            print(f"  Ratio f_QNM / f_pert     = {ratio:.3f}")
            if 0.7 < ratio < 1.4:
                print(f"  WARNING: f_QNM ~ f_pert -- cannot distinguish QNM from initial perturbation response")
            else:
                print(f"  -> Frequencies differ significantly: signal is NOT dominated by initial perturbation shape")

    _draw_psd(axes[1, 0], radii, stored_psd, linestyles,
              args.psd_smooth_window, args.psd_smooth_polyorder,
              pert_sigma=args.pert_sigma, f_max=args.esd_fmax)

    if len(radii) >= 2:
        _draw_propagation(axes[1, 1], t, radii, series, linestyles)
    else:
        axes[1, 1].text(0.5, 0.5, "Need $\\geq 2$ radii", transform=axes[1, 1].transAxes,
                        ha="center", va="center", fontsize=12)

    R_strain = radii[-1]
    R_spec = args.spectrogram_radius if (getattr(args, "spectrogram_radius", None) is not None) else radii[0]
    _draw_spectrogram(axes[2, 0], t, series[R_spec], R_spec, args, fs)

    if R_strain in stored_psd:
        freqs, psd_psi4 = stored_psd[R_strain]
        psd_psi4_s = _smooth_psd(psd_psi4, window=args.psd_smooth_window, polyorder=args.psd_smooth_polyorder)
        strain_psd_code = _psd_psi4_to_strain(freqs, psd_psi4_s)

        _draw_strain_ligo(
            axes[2, 1],
            freqs,
            strain_psd_code,
            args.mass_msun,
            args.distance_mpc,
            R_ext=R_strain,
            ligo_quantity=args.ligo_quantity,
        )
    else:
        axes[2, 1].text(0.5, 0.5, "No PSD data", transform=axes[2, 1].transAxes,
                        ha="center", va="center", fontsize=12)

    titles = [
        r"Waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$",
        r"Waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$ (retarded) + QNM fit",
        r"Energy Spectral Density of $\Psi_4$",
        "Propagation speed analysis",
        rf"Spectrogram (time-frequency, $R={R_spec:g}$)",
        r"Strain vs.\ Advanced LIGO sensitivity",
    ]
    for ax, letter, title in zip(axes.flatten(), string.ascii_lowercase, titles):
        ax.set_title(f"({letter}) {title}")

    fig.tight_layout()
    _save_figure(fig, args, radii)


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

    ax_wave = axes_arr[ax_idx]; ax_idx += 1
    _draw_waveform(ax_wave, t, radii, series, linestyles, args.time_axis,
                   t_min=args.t_min, t_max=args.t_max)
    label = (r"Radius-scaled waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$"
             if args.time_axis == "simulation"
             else r"Radius-scaled waveform $r\,\mathrm{Re}(\Psi_4^{2,0})$ (retarded)")
    ax_wave.set_title(label)

    if args.plot_psd:
        ax_psd = axes_arr[ax_idx]; ax_idx += 1
        _draw_psd(ax_psd, radii, stored_psd, linestyles,
                  args.psd_smooth_window, args.psd_smooth_polyorder,
                  pert_sigma=args.pert_sigma, f_max=args.esd_fmax)
        ax_psd.set_title(r"Power Spectral Density of $\Psi_4$")

    if args.strain:
        R_strain = radii[-1]
        
        ax_spec = axes_arr[ax_idx]; ax_idx += 1
        R_spec = args.spectrogram_radius if (getattr(args, "spectrogram_radius", None) is not None) else radii[0]
        _draw_spectrogram(ax_spec, t, series[R_spec], R_spec, args, fs)
        ax_spec.set_title(rf"Spectrogram (time-frequency, $R={R_spec:g}$)")

        if R_strain in stored_psd:
            freqs, psd_psi4 = stored_psd[R_strain]
            psd_psi4_s = _smooth_psd(psd_psi4, window=args.psd_smooth_window,
                                     polyorder=args.psd_smooth_polyorder)
            strain_psd_code = _psd_psi4_to_strain(freqs, psd_psi4_s)

            ax_ph = axes_arr[ax_idx]; ax_idx += 1
            _draw_strain_ligo(
                ax_ph,
                freqs,
                strain_psd_code,
                args.mass_msun,
                args.distance_mpc,
                R_ext=R_strain,
                ligo_quantity=args.ligo_quantity,
            )
            ax_ph.set_title(r"Strain vs.\ Advanced LIGO sensitivity")
        else:
            print("WARNING: No PSD data for strain conversion.")
            ax_idx += 1

    if args.propagation_speed and len(radii) >= 2:
        ax_prop = axes_arr[ax_idx]; ax_idx += 1
        _draw_propagation(ax_prop, t, radii, series, linestyles)
        ax_prop.set_title("Propagation speed analysis")

    fig.tight_layout()
    _save_figure(fig, args, radii)


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
    parser.add_argument(
        "--esd-fmax",
        type=float,
        default=None,
        help="Maximum frequency (in code units M^-1) to show on the ESD panel. "
             "Example: --esd-fmax 20 cuts off the flat high-frequency tail.",
    )

    parser.add_argument("--strain", action="store_true", help="Enable strain PSD conversion and LIGO overlay")
    parser.add_argument("--mass-msun", type=float, default=30.0, help="Total mass in solar masses (default: 30)")
    parser.add_argument("--distance-mpc", type=float, default=10.0, help="Luminosity distance in Mpc (default: 10)")
    parser.add_argument(
        "--ligo-quantity",
        choices=["asd", "hchar"],
        default="asd",
        help="Quantity to plot in the LIGO comparison panel: "
             "'asd' plots amplitude spectral density sqrt(S_h) and sqrt(S_n) (paper-style); "
             "'hchar' plots characteristic strain sqrt(f S_h) and sqrt(f S_n).",
    )
    parser.add_argument("--propagation-speed", action="store_true", help="Measure propagation speed across extraction radii")
    parser.add_argument(
        "--pert-sigma", type=float, default=None,
        help="Width sigma_K of the initial Gaussian K perturbation (code units). "
             "When set, overlays the perturbation's spectral content on the ESD panel "
             "to distinguish initial-data artifacts from genuine QNM ringdown.",
    )
    parser.add_argument("--pert-A0", type=float, default=0.0,
                        help="Monopole amplitude A0 of the initial K perturbation")
    parser.add_argument("--pert-A2", type=float, default=0.0,
                        help="Quadrupole amplitude A2 of the initial K perturbation")
    parser.add_argument(
        "--combined", action="store_true",
        help="Produce a single 3x2 publication figure with all analysis panels."
    )
    parser.add_argument(
        "--spectrogram-radius",
        type=float,
        default=None,
        help="Extraction radius R to use for the time-frequency spectrogram panel. "
             "Default: use the smallest available radius (to preserve the retarded-time ringdown tail).",
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
