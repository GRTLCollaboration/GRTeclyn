#!/usr/bin/env python3
"""
Search public LIGO data for a wormhole collapse transient.
This script uses the simulated Psi4 waveform to build a template,
downloads GW190521 data from GWOSC, and performs matched filtering.

Run from the repository root (``-m`` takes a module name, not ``*.py``)::

    python -m src.search.main

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

project_root = Path(__file__).resolve().parent.parent.parent
_DEFAULT_PSI4_DAT = project_root / "src/SimResults/Run_R0.5_A00.0_A20.02_sigma0.5_perturbed/psi4_mode_l2m0.dat"
sys.path.append(str(project_root))
from src.visualisation.process_wave.plot_extracted_psi4 import load_extracted, M_SUN_SEC, M_SUN_METER, MPC_METER

def _resolve_event_key(cat: Catalog, event: str) -> str:
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

    candidates = [
        e for e in strain_entries
        if e.get("detector") == ifo and int(e.get("sampling_rate", -1)) == int(sample_rate)
    ]
    if not candidates:
        dets = sorted({e.get("detector") for e in strain_entries if e.get("detector")})
        rates = sorted({e.get("sampling_rate") for e in strain_entries if e.get("sampling_rate")})
        raise KeyError(f"No strain entries for {ifo} at {sample_rate} Hz. Available detectors={dets}, rates={rates}")

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
    request_timeout_s: float | None,
    out_q: "mp.Queue",
    tmp_file: str,
) -> None:
    """Module-scope worker for multiprocessing 'spawn' (must be picklable)."""
    try:
        kw = {}
        if request_timeout_s is not None:
            kw["timeout"] = float(request_timeout_s)
        gw = GwpyTimeSeries.fetch_open_data(
            ifo,
            float(start),
            float(end),
            sample_rate=int(sample_rate),
            host=host,
            cache=cache,
            **kw,
        )
        # Save to disk instead of pipe to avoid deadlocks with large arrays
        np.savez(tmp_file, values=np.asarray(gw.value, dtype=np.float64), dt=float(gw.dt.value), epoch=float(gw.t0.value))
        out_q.put("ok")
    except Exception as e:
        out_q.put(repr(e))


def _gwosc_fetch_worker(
    event_dict: dict,
    ifo: str,
    sample_rate: int,
    host: str,
    cache: bool,
    request_timeout_s: float | None,
    out_q: "mp.Queue",
    tmp_file: str,
) -> None:
    """Module-scope worker for multiprocessing 'spawn' (must be picklable)."""
    try:
        start, end = _select_strain_segment(event_dict, ifo=ifo, sample_rate=sample_rate)
        kw = {}
        if request_timeout_s is not None:
            kw["timeout"] = float(request_timeout_s)
        gw = GwpyTimeSeries.fetch_open_data(
            ifo,
            start,
            end,
            sample_rate=int(sample_rate),
            host=host,
            cache=cache,
            **kw,
        )
        # Save to disk instead of pipe to avoid deadlocks with large arrays
        np.savez(tmp_file, values=np.asarray(gw.value, dtype=np.float64), dt=float(gw.dt.value), epoch=float(gw.t0.value))
        out_q.put("ok")
    except Exception as e:
        out_q.put(repr(e))


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
    """GWOSC fetch with a hard wall-clock timeout (separate process)."""
    import queue
    import tempfile
    import os

    # Before spawning process, see if we already have the cached .npz
    # (Optional explicit caching beyond what gwpy does)
    start, end = _select_strain_segment(event_dict, ifo=ifo, sample_rate=sample_rate)
    cache_dir = Path(__file__).parent / "data_cache"
    cache_dir.mkdir(exist_ok=True)
    cache_file = cache_dir / f"{ifo}_{start}_{end}_{sample_rate}.npz"
    
    if cache_file.exists():
        print(f"Loading {ifo} data from local disk cache: {cache_file.name}")
        data = np.load(str(cache_file))
        return TimeSeries(data["values"], delta_t=float(data["dt"]), epoch=float(data["epoch"]))

    ctx = mp.get_context("spawn")
    q: mp.Queue = ctx.Queue(maxsize=1)
    
    # We will tell the worker to save directly to the cache_file
    tmp_path = str(cache_file)
    
    p = ctx.Process(
        target=_gwosc_fetch_worker,
        args=(event_dict, ifo, sample_rate, host, cache, gwosc_request_timeout_s, q, tmp_path),
        daemon=True,
    )
    p.start()
    
    try:
        status = q.get(timeout=float(timeout_s))
    except queue.Empty:
        p.terminate()
        p.join(2.0)
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise TimeoutError(f"Downloading data timed out after {timeout_s:.0f} seconds")

    p.join(2.0)

    if status != "ok":
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise RuntimeError(f"GWOSC fetch failed: {status}")

    data = np.load(tmp_path)
    values = data["values"]
    dt = float(data["dt"])
    epoch = float(data["epoch"])
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
    import queue
    import tempfile
    import os

    cache_dir = Path(__file__).parent / "data_cache"
    cache_dir.mkdir(exist_ok=True)
    cache_file = cache_dir / f"{ifo}_{float(start)}_{float(end)}_{sample_rate}.npz"

    if cache_file.exists():
        print(f"Loading {ifo} data from local disk cache: {cache_file.name}")
        data = np.load(str(cache_file))
        return TimeSeries(data["values"], delta_t=float(data["dt"]), epoch=float(data["epoch"]))

    ctx = mp.get_context("spawn")
    q: mp.Queue = ctx.Queue(maxsize=1)
    
    tmp_path = str(cache_file)

    p = ctx.Process(
        target=_gwosc_fetch_worker_range,
        args=(ifo, float(start), float(end), int(sample_rate), host, cache, gwosc_request_timeout_s, q, tmp_path),
        daemon=True,
    )
    p.start()
    
    try:
        status = q.get(timeout=float(timeout_s))
    except queue.Empty:
        p.terminate()
        p.join(2.0)
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise TimeoutError(f"Downloading data timed out after {timeout_s:.0f} seconds")

    p.join(2.0)

    if status != "ok":
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise RuntimeError(f"GWOSC fetch failed: {status}")

    data = np.load(tmp_path)
    values = data["values"]
    dt = float(data["dt"])
    epoch = float(data["epoch"])
    return TimeSeries(values, delta_t=dt, epoch=epoch)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Matched-filter search in GWOSC data")
    p.add_argument(
        "--data-path",
        default=str(_DEFAULT_PSI4_DAT),
        help="Path to extracted psi4_mode_l2m0.dat (override with your simulation output)",
    )
    p.add_argument("--mass-msun", type=float, default=1000.0, help="Template mass in solar masses")
    p.add_argument("--distance-mpc", type=float, default=0.002, help="Template distance in Mpc (e.g. 0.002 for 2 kpc)")
    p.add_argument("--catalog", default="GWTC-2", help="PyCBC catalog name (e.g. GWTC-2). Used if --event is specified.")
    p.add_argument(
        "--event",
        default="GW190521",
        help="Event key or prefix (e.g. GW190521), or 'all' for catalog search. If provided, overrides --run.",
    )
    p.add_argument("--run", action="store_true", help="Search a continuous block of open data instead of catalog events.")
    p.add_argument("--gps-start", type=float, help="GPS start time for continuous search (required if --run is used)")
    p.add_argument("--gps-end", type=float, help="GPS end time for continuous search (required if --run is used)")
    p.add_argument("--chunk-duration", type=float, default=512.0, help="Duration of chunks to split continuous data into (seconds)")
    p.add_argument("--snr-threshold", type=float, default=8.0, help="SNR threshold to report triggers in continuous search")
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

    data_path = Path(args.data_path)
    if not data_path.exists():
        print(f"Warning: {data_path} not found. Make sure the path is correct.")
        return 1

    print(f"Loading Psi4 data from {data_path}...")
    times, radii, data = load_extracted(data_path)

    if len(radii) > 0:
        R = max(radii)
        print(f"Using extraction radius R = {R}")
        psi4_complex = data[R]
    else:
        raise ValueError("No extraction radii found.")

    r_psi4 = R * psi4_complex

    dt_code = np.median(np.diff(times))
    N = len(times)

    win = tukey(N, alpha=0.1)
    r_psi4_windowed = r_psi4 * win

    psi4_tilde = np.fft.rfft(np.real(r_psi4_windowed))
    freqs = np.fft.rfftfreq(N, d=dt_code)

    omega = 2.0 * np.pi * freqs
    omega[0] = 1e-15

    h_tilde = -psi4_tilde / (omega**2)

    f_max_code = freqs[-1]
    f_min_code = 0.05 * f_max_code

    hp_mask = 1.0 / (1.0 + (f_min_code / omega)**8)
    h_tilde *= hp_mask

    h_strain_code = np.fft.irfft(h_tilde, n=N)

    nr_strain = np.real(h_strain_code)

    M_wormhole = float(args.mass_msun)
    distance_mpc = float(args.distance_mpc)
    
    T_M = M_wormhole * M_SUN_SEC
    D_m = distance_mpc * MPC_METER

    time_sec = times * T_M

    amp_scale = (M_wormhole * M_SUN_METER) / D_m
    nr_strain_phys = nr_strain * amp_scale

    sample_rate = int(args.sample_rate)
    dt_ligo = 1.0 / sample_rate
    time_interp = np.arange(time_sec[0], time_sec[-1], dt_ligo)

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

    N_interp = len(strain_interp)
    taper_win = tukey(N_interp, alpha=0.1)
    strain_interp *= taper_win

    template = TimeSeries(strain_interp, delta_t=dt_ligo)

    if args.run:
        if args.gps_start is None or args.gps_end is None:
            print("Error: --gps-start and --gps-end must be provided when using --run")
            return 1
        
        print(f"\n=== Continuous Search Mode ===")
        print(f"Detector: {args.ifo}")
        print(f"Time: {args.gps_start} to {args.gps_end} ({args.gps_end - args.gps_start} seconds)")
        print(f"Chunk Size: {args.chunk_duration} seconds")
        print(f"SNR Threshold: {args.snr_threshold}\n")
        
        triggers = []
        cache = not bool(args.no_cache)
        chunk_start = float(args.gps_start)
        
        while chunk_start < float(args.gps_end):
            chunk_end = min(chunk_start + float(args.chunk_duration), float(args.gps_end))
            
            # We need to fetch a slightly larger chunk of data to avoid edge effects in the filter
            # PyCBC recommends padding the data by at least the maximum template length + some extra for corruption
            pad = 16.0
            fetch_start = chunk_start - pad
            fetch_end = chunk_end + pad
            
            print(f"Processing chunk: {chunk_start} to {chunk_end} ...")
            
            try:
                data_h1 = _fetch_open_data_range_pycbc_with_timeout(
                    args.ifo,
                    fetch_start,
                    fetch_end,
                    sample_rate,
                    timeout_s=float(args.fetch_timeout_s),
                    cache=cache,
                    gwosc_request_timeout_s=float(args.gwosc_request_timeout_s),
                )
                
                current_template = template.copy()
                current_template.resize(len(data_h1))
                current_template = current_template.cyclic_time_shift(current_template.start_time)

                # Estimate PSD
                psd = data_h1.psd(4)
                psd = pycbc_interpolate(psd, data_h1.delta_f)
                psd = inverse_spectrum_truncation(
                    psd, int(4 * data_h1.sample_rate), low_frequency_cutoff=15.0
                )

                # Filter
                snr = matched_filter(current_template, data_h1, psd=psd, low_frequency_cutoff=float(args.low_frequency_cutoff))
                snr_timeseries = abs(snr)
                
                # Crop the padding off the SNR timeseries to remove edge corruption
                snr_timeseries = snr_timeseries.crop(pad + 4.0, pad + 4.0)

                # Find peaks above threshold using scipy.signal.find_peaks to avoid clustering
                from scipy.signal import find_peaks
                peaks, _ = find_peaks(snr_timeseries, height=float(args.snr_threshold), distance=int(0.5 * sample_rate))
                
                for i in peaks:
                    t = snr_timeseries.sample_times[i]
                    # Ensure we don't double count if a peak is right on the boundary
                    if chunk_start <= t < chunk_end:
                        s = snr_timeseries[i]
                        triggers.append((t, s))
                        print(f"  *** TRIGGER FOUND *** GPS: {t:.2f} | SNR: {s:.2f}")
                            
            except Exception as e:
                print(f"  Skipping chunk due to error: {e}")
                
            chunk_start += float(args.chunk_duration)
            
        print("\n=== CONTINUOUS SEARCH RESULTS ===")
        if not triggers:
            print(f"No triggers found above SNR {args.snr_threshold}.")
        else:
            triggers.sort(key=lambda x: x[1], reverse=True)
            print(f"Found {len(triggers)} triggers above threshold:")
            for t, s in triggers[:20]:
                print(f"  SNR = {s:.2f} at GPS {t:.2f}")
        return 0

    # --- Original Catalog Mode Below ---
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

            current_template = template.copy()
            current_template.resize(len(data_h1))
            current_template = current_template.cyclic_time_shift(current_template.start_time)

            psd = data_h1.psd(4)
            psd = pycbc_interpolate(psd, data_h1.delta_f)
            psd = inverse_spectrum_truncation(
                psd, int(4 * data_h1.sample_rate), low_frequency_cutoff=15.0
            )

            snr = matched_filter(current_template, data_h1, psd=psd, low_frequency_cutoff=float(args.low_frequency_cutoff))
            snr_timeseries = abs(snr)

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
        plt.figure(figsize=(10, 5))
        plt.plot(snr_timeseries.sample_times, snr_timeseries, color="purple")
        plt.axhline(8.0, color="red", linestyle="--", label="Reference threshold (SNR = 8)")
        plt.ylabel("Signal-to-noise ratio (|SNR|)")
        plt.xlabel("GPS time (s)")
        plt.title(
            rf"Search for $1000\,M_\odot$ (D=1 Mpc) Wormhole Template in {ev_name} ({args.ifo})" + "\n" +
            f"Peak SNR = {peak_snr:.2f}"
        )

        try:
            true_event_time = float(cat.data[ev_name].get("GPS"))
        except KeyError:
            true_event_time = peak_time
        except Exception:
            true_event_time = peak_time

        plt.xlim(true_event_time - 16.0, true_event_time + 16.0)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plot_path = Path(__file__).parent / "search_snr.png"
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {plot_path}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
