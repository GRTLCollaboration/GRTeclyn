#!/usr/bin/env python3
"""
Search public LIGO data for a wormhole collapse transient.
This script uses the simulated Psi4 waveform to build a template,
downloads GW190521 data from GWOSC, and performs matched filtering.

Dependencies:
    pip install pycbc gwpy scipy numpy matplotlib
"""

import sys
from pathlib import Path
import argparse
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal.windows import tukey

from pycbc.types import TimeSeries
from pycbc.filter import matched_filter
from pycbc.psd import interpolate as pycbc_interpolate, inverse_spectrum_truncation
from pycbc.catalog import Catalog
from gwpy.timeseries import TimeSeries as GwpyTimeSeries

# Add the visualization path to import load_extracted
project_root = Path(__file__).resolve().parent.parent.parent
_DEFAULT_PSI4_DAT = project_root / "src/SimResults/Run_R0.5_A00.0_A20.02_sigma0.5_perturbed/psi4_mode_l2m0.dat"
sys.path.append(str(project_root))
from src.visualisation.process_wave.plot_extracted_psi4 import load_extracted, M_SUN_SEC, M_SUN_METER, MPC_METER

def _resolve_event_key(cat: Catalog, event: str) -> str:
    # PyCBC catalog keys are often full IDs like "GW190521_074359"
    if event in cat.data:
        return event
    keys = list(cat.data.keys())
    matches = [k for k in keys if (k == event) or k.startswith(event) or (event in k)]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        return sorted(matches)[0]
    raise KeyError(
        f"Event '{event}' not found in catalog. "
        f"Example keys include: {', '.join(keys[:10])}"
    )


def _select_strain_segment(event_dict: dict, ifo: str, sample_rate: int) -> tuple[float, float]:
    """Pick a (start, end) segment from GWOSC metadata in the PyCBC catalog dict."""
    strain_entries = event_dict.get("strain", [])
    if not strain_entries:
        raise KeyError("No 'strain' entries found in event metadata.")

    # Prefer short segments (32s) at the requested sample rate.
    candidates = [
        e for e in strain_entries
        if e.get("detector") == ifo and int(e.get("sampling_rate", -1)) == int(sample_rate)
    ]
    if not candidates:
        dets = sorted({e.get("detector") for e in strain_entries if e.get("detector")})
        rates = sorted({e.get("sampling_rate") for e in strain_entries if e.get("sampling_rate")})
        raise KeyError(f"No strain entries for {ifo} at {sample_rate} Hz. Available detectors={dets}, rates={rates}")

    # Prefer 32 s, else take the smallest available duration.
    candidates_sorted = sorted(candidates, key=lambda e: (abs(int(e.get("duration", 0)) - 32), int(e.get("duration", 0))))
    chosen = candidates_sorted[0]
    start = float(chosen["GPSstart"])
    dur = float(chosen["duration"])
    return start, start + dur


def _fetch_open_data_pycbc(event_dict: dict, ifo: str, sample_rate: int) -> TimeSeries:
    """Fetch GWOSC open data via GWpy and convert to PyCBC TimeSeries."""
    start, end = _select_strain_segment(event_dict, ifo=ifo, sample_rate=sample_rate)
    print(f"Fetching GWOSC open data: {ifo}, [{start}, {end}) at {sample_rate} Hz")
    gw = GwpyTimeSeries.fetch_open_data(ifo, start, end, sample_rate=sample_rate)
    # Convert to PyCBC TimeSeries (epoch in GPS seconds)
    return TimeSeries(gw.value, delta_t=float(gw.dt.value), epoch=float(gw.t0.value))

def _event_gps(event_dict: dict) -> float | None:
    gps = event_dict.get("GPS")
    if gps is None:
        return None
    try:
        return float(gps)
    except Exception:
        return None


def _gwosc_fetch_worker_range(
    ifo: str,
    start: float,
    end: float,
    sample_rate: int,
    host: str,
    cache: bool,
    out_q: "mp.Queue",
) -> None:
    """Module-scope worker for multiprocessing 'spawn' (must be picklable)."""
    try:
        gw = GwpyTimeSeries.fetch_open_data(
            ifo,
            float(start),
            float(end),
            sample_rate=float(sample_rate),
            host=host,
            cache=cache,
        )
        out_q.put(("ok", (np.asarray(gw.value, dtype=np.float64), float(gw.dt.value), float(gw.t0.value))))
    except Exception as e:
        out_q.put(("err", repr(e)))


def _gwosc_fetch_worker(
    event_dict: dict,
    ifo: str,
    sample_rate: int,
    host: str,
    cache: bool,
    out_q: "mp.Queue",
) -> None:
    """Module-scope worker for multiprocessing 'spawn' (must be picklable)."""
    try:
        start, end = _select_strain_segment(event_dict, ifo=ifo, sample_rate=sample_rate)
        gw = GwpyTimeSeries.fetch_open_data(
            ifo,
            start,
            end,
            sample_rate=sample_rate,
            host=host,
            cache=cache,
        )
        out_q.put(("ok", (np.asarray(gw.value, dtype=np.float64), float(gw.dt.value), float(gw.t0.value))))
    except Exception as e:
        out_q.put(("err", repr(e)))


def _fetch_open_data_pycbc_with_timeout(
    event_dict: dict,
    ifo: str,
    sample_rate: int,
    *,
    timeout_s: float = 30.0,
    cache: bool = True,
    host: str = "https://gwosc.org",
    gwosc_request_timeout_s: float | None = 20.0,
) -> TimeSeries:
    """
    Robust GWOSC fetch with a hard wall-clock timeout.

    Why: looping over a whole catalog can "hang" if a single GWOSC request stalls.
    Threads can't reliably stop a stuck network call; a separate process can.
    """

    ctx = mp.get_context("spawn")
    q: mp.Queue = ctx.Queue(maxsize=1)
    p = ctx.Process(
        target=_gwosc_fetch_worker,
        args=(event_dict, ifo, sample_rate, host, cache, q),
        daemon=True,
    )
    p.start()
    p.join(timeout_s)
    if p.is_alive():
        p.terminate()
        p.join(2.0)
        raise TimeoutError(f"Downloading data timed out after {timeout_s:.0f} seconds")

    if q.empty():
        raise RuntimeError("GWOSC fetch failed without returning an error message.")

    status, payload = q.get()
    if status != "ok":
        raise RuntimeError(f"GWOSC fetch failed: {payload}")

    values, dt, epoch = payload
    return TimeSeries(values, delta_t=dt, epoch=epoch)


def _fetch_open_data_range_pycbc_with_timeout(
    ifo: str,
    start: float,
    end: float,
    sample_rate: int,
    *,
    timeout_s: float = 300.0,
    cache: bool = True,
    host: str = "https://gwosc.org",
    gwosc_request_timeout_s: float | None = 30.0,
) -> TimeSeries:
    """Fetch an explicit [start, end) GPS window with a hard wall-clock timeout."""
    ctx = mp.get_context("spawn")
    q: mp.Queue = ctx.Queue(maxsize=1)
    p = ctx.Process(
        target=_gwosc_fetch_worker_range,
        args=(ifo, float(start), float(end), int(sample_rate), host, cache, q),
        daemon=True,
    )
    p.start()
    p.join(float(timeout_s))
    if p.is_alive():
        p.terminate()
        p.join(2.0)
        raise TimeoutError(f"Downloading data timed out after {timeout_s:.0f} seconds")

    if q.empty():
        raise RuntimeError("GWOSC fetch failed without returning an error message.")

    status, payload = q.get()
    if status != "ok":
        raise RuntimeError(f"GWOSC fetch failed: {payload}")

    values, dt, epoch = payload
    return TimeSeries(values, delta_t=dt, epoch=epoch)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Matched-filter search in GWOSC data")
    p.add_argument(
        "--data-path",
        default=str(_DEFAULT_PSI4_DAT),
        help="Path to extracted psi4_mode_l2m0.dat (override with your simulation output)",
    )
    p.add_argument("--mass-msun", type=float, default=1000.0, help="Template mass in solar masses")
    p.add_argument("--distance-mpc", type=float, default=1.0, help="Template distance in Mpc")
    p.add_argument("--catalog", default="GWTC-2", help="PyCBC catalog name (e.g. GWTC-2)")
    p.add_argument(
        "--event",
        default="GW190521",
        help="Event key or prefix (e.g. GW190521 or GW190521_074359)",
    )
    p.add_argument("--ifo", default="H1", help="Detector (H1 or L1)")
    p.add_argument("--sample-rate", type=int, default=4096, help="Template sample rate (Hz)")
    p.add_argument("--low-frequency-cutoff", type=float, default=20.0, help="Matched-filter low-f cutoff (Hz)")
    p.add_argument("--fetch-timeout-s", type=float, default=300.0, help="Hard timeout for each GWOSC fetch (seconds)")
    p.add_argument("--fetch-retries", type=int, default=1, help="Retries per event if GWOSC fetch fails/timeouts")
    p.add_argument("--gwosc-request-timeout-s", type=float, default=30.0, help="Per-request network timeout passed to GWpy/GWOSC (seconds)")
    p.add_argument("--no-cache", action="store_true", help="Disable GWOSC caching (default caches downloads)")
    p.add_argument("--max-events", type=int, default=0, help="For --event all: stop after N events (0 = no limit)")
    p.add_argument(
        "--fetch-strategy",
        choices=["gps", "segment"],
        default="gps",
        help="gps: fetch a fixed window around event GPS (recommended for full-catalog loops); segment: use catalog-provided strain segments",
    )
    p.add_argument("--window-s", type=float, default=16.0, help="Half-width of GPS window (seconds) when --fetch-strategy gps")
    args = p.parse_args(argv)

    # --- 1. Load and scale NR data ---
    data_path = Path(args.data_path)
    if not data_path.exists():
        # Fallback to current directory or fail
        print(f"Warning: {data_path} not found. Make sure the path is correct.")
        return 1

    print(f"Loading Psi4 data from {data_path}...")
    times, radii, data = load_extracted(data_path)

    # Pick the outermost radius available for the best signal, e.g. R=24 if available
    if len(radii) > 0:
        R = max(radii)
        print(f"Using extraction radius R = {R}")
        psi4_complex = data[R]
    else:
        raise ValueError("No extraction radii found.")

    # Multiply by R to get the asymptotic waveform r*Psi4
    r_psi4 = R * psi4_complex

    # Convert Psi4 to strain in frequency domain: h_tilde = - Psi4_tilde / (2 * pi * f)**2
    # First, regularize and window the time domain signal
    dt_code = np.median(np.diff(times))
    N = len(times)
    
    # Tapering to avoid spectral leakage
    win = tukey(N, alpha=0.1)
    r_psi4_windowed = r_psi4 * win

    # Extract the real part (h_+ polarization) before taking the Real-FFT
    psi4_tilde = np.fft.rfft(np.real(r_psi4_windowed))
    freqs = np.fft.rfftfreq(N, d=dt_code)

    # Avoid division by zero at f=0
    omega = 2.0 * np.pi * freqs
    omega[0] = 1e-15 # small number

    # Frequency-domain integration
    h_tilde = -psi4_tilde / (omega**2)

    # High-pass filter in frequency domain to avoid low-frequency amplification artifacts
    # The signal is a short burst, so we can cut off frequencies below some f_min (in code units)
    f_max_code = freqs[-1]
    f_min_code = 0.05 * f_max_code
    
    # 4th-order Butterworth-like high-pass suppression for f < f_min
    hp_mask = 1.0 / (1.0 + (f_min_code / omega)**8)
    h_tilde *= hp_mask

    # Inverse FFT to get strain in time domain (still in code units, scaled by R)
    h_strain_code = np.fft.irfft(h_tilde, n=N)

    # Select real part (one polarization) for the template
    nr_strain = np.real(h_strain_code)

    # --- Scale to physical units ---
    M_wormhole = float(args.mass_msun)  # solar masses
    distance_mpc = float(args.distance_mpc)
    
    T_M = M_wormhole * M_SUN_SEC
    D_m = distance_mpc * MPC_METER
    
    # Physical time
    time_sec = times * T_M
    
    # Physical strain: h_phys = h_code * (M * M_sun_meter / D)
    amp_scale = (M_wormhole * M_SUN_METER) / D_m
    nr_strain_phys = nr_strain * amp_scale

    sample_rate = int(args.sample_rate)
    dt_ligo = 1.0 / sample_rate
    time_interp = np.arange(time_sec[0], time_sec[-1], dt_ligo)

    # The extracted time array can contain duplicate timestamps; de-duplicate
    # before interpolation to a uniform 4096 Hz grid.
    t_unique, unique_idx = np.unique(time_sec, return_index=True)
    h_unique = nr_strain_phys[unique_idx]
    if t_unique.size < 4:
        raise ValueError("Too few unique time samples for interpolation.")

    interpolator = interp1d(
        t_unique,
        h_unique,
        kind="linear",
        fill_value=0.0,
        bounds_error=False,
        assume_sorted=True,
    )
    strain_interp = interpolator(time_interp)

    # Taper in numpy to avoid chasing pycbc version imports
    N_interp = len(strain_interp)
    taper_win = tukey(N_interp, alpha=0.1)
    strain_interp *= taper_win

    template = TimeSeries(strain_interp, delta_t=dt_ligo)

    # --- 2. LIGO strain search ---
    print(f"Loading LIGO Catalog '{args.catalog}'...")
    cat = Catalog(args.catalog)

    if args.event.lower() == "all":
        events_to_search = list(cat.data.keys())
        print(f"Searching all {len(events_to_search)} events in the catalog...")
        if int(args.max_events) > 0:
            events_to_search = events_to_search[: int(args.max_events)]
            print(f"Limiting to first {len(events_to_search)} events due to --max-events.")
    else:
        event_key = _resolve_event_key(cat, args.event)
        if event_key != args.event:
            print(f"Resolved event key: {event_key}")
        events_to_search = [event_key]

    results = []

    for ev_name in events_to_search:
        try:
            print(f"\n--- Searching {ev_name} ({args.ifo}) ---")
            catalog_event = cat.data[ev_name]

            cache = not bool(args.no_cache)
            last_err: Exception | None = None
            data_h1 = None
            for attempt in range(int(args.fetch_retries) + 1):
                try:
                    print(f"Fetching strain (attempt {attempt + 1}/{int(args.fetch_retries) + 1})...")
                    if args.fetch_strategy == "gps":
                        gps = _event_gps(catalog_event)
                        if gps is None:
                            raise KeyError("Event has no GPS field; use --fetch-strategy segment for this catalog.")
                        start = gps - float(args.window_s)
                        end = gps + float(args.window_s)
                    else:
                        start, end = _select_strain_segment(catalog_event, ifo=args.ifo, sample_rate=sample_rate)

                    data_h1 = _fetch_open_data_range_pycbc_with_timeout(
                        args.ifo,
                        start,
                        end,
                        sample_rate,
                        timeout_s=float(args.fetch_timeout_s),
                        cache=cache,
                        gwosc_request_timeout_s=float(args.gwosc_request_timeout_s),
                    )
                    break
                except Exception as e:
                    last_err = e
                    print(f"Fetch failed: {e}")
            if data_h1 is None:
                raise RuntimeError(f"All fetch attempts failed. Last error: {last_err}")

            # Resize template and shift
            current_template = template.copy()
            current_template.resize(len(data_h1))
            current_template = current_template.cyclic_time_shift(current_template.start_time)

            # Calculate PSD
            psd = data_h1.psd(4)
            psd = pycbc_interpolate(psd, data_h1.delta_f)
            psd = inverse_spectrum_truncation(
                psd, int(4 * data_h1.sample_rate), low_frequency_cutoff=15.0
            )

            # Matched filter
            snr = matched_filter(current_template, data_h1, psd=psd, low_frequency_cutoff=float(args.low_frequency_cutoff))
            snr_timeseries = abs(snr)

            # Throw away the first and last 4 seconds of corrupted FFT data
            snr_timeseries = snr_timeseries.crop(4.0, 4.0)

            peak_snr = snr_timeseries.max()
            peak_time = snr_timeseries.sample_times[np.argmax(snr_timeseries)]
            print(f"--> Peak SNR for {ev_name}: {peak_snr:.2f} at GPS time {peak_time:.2f}")
            if args.event.lower() == "all":
                results.append((ev_name, peak_snr, peak_time))
            else:
                results.append((ev_name, peak_snr, peak_time, snr_timeseries))

        except Exception as e:
            print(f"Skipping {ev_name} due to error: {e}")

    # After the loop finishes
    if args.event.lower() == "all":
        results.sort(key=lambda x: x[1], reverse=True)
        print("\n=== TOP 5 WORMHOLE CANDIDATES ===")
        for ev, snr, ptime in results[:5]:
            print(f"{ev}: SNR = {snr:.2f} at GPS {ptime:.2f}")
    else:
        if not results:
            print("No results to plot.")
            return 1
            
        ev_name, peak_snr, peak_time, snr_timeseries = results[0]
        # --- Plot ---
        plt.figure(figsize=(10, 5))
        plt.plot(snr_timeseries.sample_times, snr_timeseries, color="purple")
        plt.axhline(8.0, color="red", linestyle="--", label="Reference threshold (SNR = 8)")
        plt.ylabel("Signal-to-noise ratio (|SNR|)")
        plt.xlabel("GPS time (s)")
        plt.title(
            rf"Search for $1000\,M_\odot$ (D=1 Mpc) Wormhole Template in {ev_name} ({args.ifo})" + "\n" +
            f"Peak SNR = {peak_snr:.2f}"
        )
        
        # Try to use the event GPS time if possible, otherwise use peak time
        try:
            true_event_time = float(cat.data[ev_name].get("GPS"))
        except KeyError:
            true_event_time = peak_time
        except Exception:
            true_event_time = peak_time

        # Zoom in 16 seconds around the actual event
        plt.xlim(true_event_time - 16.0, true_event_time + 16.0)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plot_path = Path(__file__).parent / "search_snr.png"
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {plot_path}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
