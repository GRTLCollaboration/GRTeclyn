from __future__ import annotations

import argparse
import json
import os
import shutil
import time
from pathlib import Path

import yt

from grteclyn_wrapper.objective_modes import QD_OBJECTIVE_MODES

from .config import _default_data_dir, _default_frames_out_dir, _frames_auto_zlim_enabled
from .extraction.central import CENTRAL_TIMESERIES_HEADER
from .extraction.confinement import CONFINEMENT_TIMESERIES_HEADER
from .extraction.sector_barycenters import SECTOR_BARYCENTERS_HEADER
from .extraction.ftl import FTL_TIMESERIES_HEADER
from .extraction.shell import _shell_stats_header
from .fields import _canonical_field_name
from .frames.cleanup import (
    _cleanup_embedding_frames,
    _cleanup_existing_frames,
    _cleanup_projection_frames,
)
from .frames.zlim import _lock_frame_zlims_from_plotfile
from .plotfiles import (
    _is_plotfile_ready,
    _iter_plotfile_dirs,
    _parse_plot_index,
    _should_auto_reset,
    _truncate_if_exists,
)
from .state import _append_line, _load_state, _save_state
from .worker import _process_single_plotfile

from grteclyn_wrapper.metrics.aggregation.incremental import IncrementalScoreWriter
from grteclyn_wrapper.metrics.splash_early_term import evaluate_splash_early_term


def _append_radial_block(path: Path, block: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(block)


def _parse_central_row(line: str) -> tuple[float, float | None, float | None, float | None]:
    parts = line.split()
    rho = float(parts[1]) if len(parts) > 1 else None
    lapse = float(parts[2]) if len(parts) > 2 else None
    activity = float(parts[3]) if len(parts) > 3 else None
    return rho, lapse, activity


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract Psi4 (l=2,m=0) from plotfiles and optionally delete them."
    )
    parser.add_argument("--data", default=_default_data_dir(), help="Directory containing WormholePlt*/plt*")
    parser.add_argument("--out", default=None, help="Output directory (default: <data>/small_data)")
    parser.add_argument("--radii", type=float, nargs="+", default=[14.0, 30.0], help="Extraction radii")
    parser.add_argument("--n-points", type=int, default=64, help="Angular resolution N (N×N points)")
    parser.add_argument(
        "--center",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        help="Extraction center (x y z) in code units",
    )
    parser.add_argument("--stable-seconds", type=float, default=30.0, help="Require Header mtime older than this (30s default for NFS)")
    parser.add_argument("--poll-seconds", type=float, default=2.0, help="Polling interval when --watch")
    parser.add_argument("--watch", action="store_true", help="Keep running and process new plotfiles")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-plotfile details (name, time, delete/keep, timings).",
    )
    parser.add_argument(
        "--psi4",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable/disable Psi4 mode extraction to .dat (default: enabled).",
    )
    parser.add_argument(
        "--shell-fields",
        nargs="+",
        default=[],
        help="Extract mean/min/max on spherical shells for these fields (e.g. chi lapse K).",
    )
    parser.add_argument(
        "--frames-fields",
        nargs="+",
        default=[],
        help="Render SlicePlot frames for these fields (e.g. chi K Weyl4_Re).",
    )
    parser.add_argument("--frames-axis", default="z", choices=["x", "y", "z"], help="Slice axis for frames.")
    parser.add_argument("--frames-coord", type=float, default=None, help="Slice coordinate for frames (axis coordinate).")
    parser.add_argument("--frames-zoom", type=float, default=None, help="Zoom width in code units for frames.")
    parser.add_argument("--frames-center", type=float, nargs=3, default=None, help="Center (x y z) for frames.")
    parser.add_argument("--frames-corner", action="store_true", help="Corner mode for symmetry-reduced domains (frames).")
    parser.add_argument(
        "--frames-auto-zlim",
        action="store_true",
        help="Per-frame colorbar from slice percentiles (faint fields only; causes movie flicker).",
    )
    parser.add_argument(
        "--frames-global-zlim",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Lock per-field colorbar limits from the first plotfile (stable movies, visible faint fields).",
    )
    parser.add_argument("--frames-out", default=_default_frames_out_dir(), help="Frames output base dir (default: grteclyn_wrapper/visualisation/visualize).")
    parser.add_argument(
        "--projection-fields",
        nargs="+",
        default=[],
        help="Render 3D line-of-sight projection frames for these fields.",
    )
    parser.add_argument(
        "--projection-axes",
        nargs="+",
        choices=["x", "y", "z"],
        default=[],
        help="Axes to project along for --projection-fields.",
    )
    parser.add_argument(
        "--projection-method",
        choices=["mip", "integrate", "sum"],
        default="mip",
        help="Projection method. mip is best for locating compact blobs.",
    )
    parser.add_argument(
        "--areal-radius",
        action="store_true",
        help="Extract minimum areal radius R_areal = r/sqrt(chi) along x-axis to areal_radius.dat.",
    )
    parser.add_argument(
        "--ftl-timeseries",
        action="store_true",
        help="Per-plotfile FTL features (operational + gated geodesic) to ftl_timeseries.dat.",
    )
    parser.add_argument(
        "--confinement-timeseries",
        action="store_true",
        help="Per-plotfile matter-confinement moments (rms_radius, confined_frac, "
        "true matter barycentre) to confinement.dat -- the trustworthy "
        "'matter dispersed / flew away' detector.",
    )
    parser.add_argument(
        "--confinement-well-width",
        type=float,
        default=1.5,
        help="Lump scale for confinement R_conf = 4*well_width (default 1.5).",
    )
    parser.add_argument(
        "--sector-barycenters",
        action="store_true",
        help="Per-plotfile canonical/phantom sector centroids to "
        "sector_barycenters.dat -- the Bondi-dipole runaway diagnostic "
        "(the aggregate barycentre cancels for a mixed-sign pair).",
    )
    parser.add_argument(
        "--matter-model",
        default="",
        help="Authoritative recipe_matter_model tag for the sector split "
        "(field-sniffing cannot classify runs whose canonical sector is zero).",
    )
    parser.add_argument(
        "--central-timeseries",
        action="store_true",
        help="Per-plotfile central ball rho/lapse/scalar activity to central_timeseries.dat.",
    )
    parser.add_argument(
        "--central-ball",
        action="store_true",
        help="Average central fields over a small sphere instead of a single point.",
    )
    parser.add_argument(
        "--central-ball-radius",
        type=float,
        default=None,
        help="Central ball radius in code units (default: 2*dx_finest).",
    )
    parser.add_argument(
        "--central-radial-profile",
        action="store_true",
        help="Append rho/lapse/activity vs r blocks to central_radial_profile.dat.",
    )
    parser.add_argument(
        "--central-radial-r-max",
        type=float,
        default=6.0,
        help="Outer radius for central radial profile extraction.",
    )
    parser.add_argument(
        "--splash-early-term",
        action="store_true",
        help="Write .stop_sim when splash early-termination predicates fire.",
    )
    parser.add_argument(
        "--stop-sim-path",
        default=None,
        help="Path for early-termination sentinel JSON (default: <data>/.stop_sim).",
    )
    parser.add_argument(
        "--ftl-l",
        type=float,
        default=8.0,
        help="FTL probe half-width / corridor extent (matches collector ftl_L).",
    )
    parser.add_argument(
        "--evolving-geodesic",
        action="store_true",
        help="Set GRTECLYN_EVOLVING_GEODESIC=1 for end-of-run 4D null trace scoring.",
    )
    parser.add_argument(
        "--incremental-score",
        action="store_true",
        help="After each FTL timeseries row, append prefix score to score_timeseries.jsonl.",
    )
    parser.add_argument(
        "--objective-mode",
        default="weighted",
        choices=list(QD_OBJECTIVE_MODES),
        help="Objective mode for incremental scoring (matches final score_episode).",
    )
    parser.add_argument(
        "--target-stop-time",
        type=float,
        default=None,
        help="Configured stop time for survival fraction in incremental scores.",
    )
    parser.add_argument(
        "--score-weight",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Score weight override for incremental scoring (repeatable).",
    )
    parser.add_argument(
        "--embedding",
        action="store_true",
        help="Render 3D embedding-diagram frames (surface of revolution from chi profile).",
    )
    parser.add_argument(
        "--embedding-rmax",
        type=float,
        default=None,
        help="Maximum coordinate radius for embedding diagram (default: full domain).",
    )
    parser.add_argument(
        "--delete",
        action="store_true",
        help="Delete plotfile directory after successful extraction",
    )
    parser.add_argument(
        "--keep-last",
        type=int,
        default=2,
        help="Never delete the newest N plotfiles (safety, default: 2)",
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=1,
        help="Number of parallel worker processes to use",
    )
    parser.add_argument(
        "--keep-existing-frames",
        action="store_true",
        help="Do not delete existing frames at startup.",
    )
    args = parser.parse_args()
    if args.evolving_geodesic:
        os.environ["GRTECLYN_EVOLVING_GEODESIC"] = "1"
    args.metric_stack_cache = bool(args.evolving_geodesic)
    if args.evolving_geodesic:
        from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic_options import (
            metric_stack_n_space_from_env,
        )

        args.metric_stack_n_space = metric_stack_n_space_from_env()
    else:
        args.metric_stack_n_space = 65

    # Reduce yt logging overhead/spam (can be noisy in watch mode).
    try:
        # WARNING by default; INFO is extremely verbose.
        yt.funcs.mylog.setLevel(30 if args.verbose else 40)  # WARNING/ERROR
    except Exception:
        pass

    data_dir = os.path.abspath(args.data)
    out_dir = Path(args.out) if args.out else Path(data_dir) / "small_data"
    state_path = out_dir / "consume_state.json"
    out_path = out_dir / "psi4_mode_l2m0.dat"
    psi4_all_out_path = out_dir / "psi4_mode_l2_all.dat"
    psi4_directional_out_path = out_dir / "psi4_directional.dat"
    areal_out_path = out_dir / "areal_radius.dat"
    shell_out_path = out_dir / "shell_profiles.dat"
    boundary_flux_out_path = out_dir / "boundary_flux.dat"
    ftl_out_path = out_dir / "ftl_timeseries.dat"
    confinement_out_path = out_dir / "confinement.dat"
    sector_barycenters_out_path = out_dir / "sector_barycenters.dat"
    central_out_path = out_dir / "central_timeseries.dat"
    central_radial_out_path = out_dir / "central_radial_profile.dat"
    score_ts_path = out_dir / "score_timeseries.jsonl"
    stop_sim_path = Path(args.stop_sim_path) if args.stop_sim_path else Path(data_dir) / ".stop_sim"
    header = "# time  " + "  ".join([f"Re(R={R:g})  Im(R={R:g})" for R in args.radii])
    psi4_all_header = (
        "# time  "
        + "  ".join(
            [
                f"Re_m{m}(R={R:g})  Im_m{m}(R={R:g})"
                for R in args.radii
                for m in (-2, -1, 0, 1, 2)
            ]
        )
    )
    psi4_directional_header = "# time  P_total  P_z_beam  beam_ratio  beaming_gain  wavezone_std"
    areal_header = "# time  R_areal_min  r_at_R_areal_min"
    shell_header = _shell_stats_header(args.radii, args.shell_fields)

    state = _load_state(state_path)

    # Auto-reset on simulation restart (same output dir reused):
    plot_dirs_now = _iter_plotfile_dirs(data_dir)
    if _should_auto_reset(plot_dirs_now, state):
        print("Detected a likely simulation restart in the same output directory.")
        print(f"Resetting: {out_path} and {state_path}")
        _truncate_if_exists(out_path)
        _truncate_if_exists(psi4_all_out_path)
        _truncate_if_exists(psi4_directional_out_path)
        _truncate_if_exists(areal_out_path)
        if args.shell_fields:
            _truncate_if_exists(shell_out_path)
        if args.ftl_timeseries:
            _truncate_if_exists(ftl_out_path)
        if args.confinement_timeseries:
            _truncate_if_exists(confinement_out_path)
        if args.sector_barycenters:
            _truncate_if_exists(sector_barycenters_out_path)
        if args.central_timeseries:
            _truncate_if_exists(central_out_path)
        if args.central_radial_profile:
            _truncate_if_exists(central_radial_out_path)
        if args.incremental_score:
            _truncate_if_exists(score_ts_path)
        _save_state(state_path, {})
        state = {}

    score_weights = {}
    for pair in args.score_weight:
        if "=" in pair:
            key, value = pair.split("=", 1)
            score_weights[key.strip()] = float(value.strip())

    incremental_writer: IncrementalScoreWriter | None = None
    use_central_incremental = (
        args.incremental_score
        and args.central_timeseries
        and args.objective_mode == "critical_collapse"
    )
    if args.incremental_score and (args.ftl_timeseries or use_central_incremental):
        incremental_writer = IncrementalScoreWriter(
            Path(data_dir),
            objective_mode=args.objective_mode,
            target_stop_time=args.target_stop_time,
            ftl_L=args.ftl_l,
            score_weights=score_weights or None,
            evolving_geodesic_mode=bool(args.evolving_geodesic),
            out_path=score_ts_path,
        )

    peak_rho_state = {"value": 0.0, "last_activity": None}

    def _append_incremental_score(t: float) -> None:
        if incremental_writer is None:
            return
        try:
            record = incremental_writer.append(t)
            if record is not None and args.verbose:
                print(
                    f"[score] t={record['t']:.3g} total={record['score']:.2f} "
                    f"f_geo={record['f_geo'] * 100:.2f}% horizon={record['horizon_penalty']:.2f}",
                    flush=True,
                )
        except Exception as exc:  # noqa: BLE001 - scoring must not stop the consumer
            print(f"WARNING: incremental score at t={t:.6g} failed: {exc}")

    def _handle_central_outputs(res: dict) -> None:
        if res.get("central_line"):
            _append_line(
                central_out_path,
                header=CENTRAL_TIMESERIES_HEADER,
                line=res["central_line"],
            )
            if use_central_incremental:
                _append_incremental_score(float(res["t"]))
            if args.splash_early_term:
                rho, lapse, activity = _parse_central_row(res["central_line"])
                if rho is not None and rho > peak_rho_state["value"]:
                    peak_rho_state["value"] = rho
                decision = evaluate_splash_early_term(
                    t=float(res["t"]),
                    rho=rho,
                    lapse=lapse,
                    activity=activity,
                    peak_rho_so_far=peak_rho_state["value"],
                    previous_activity=peak_rho_state["last_activity"],
                )
                peak_rho_state["last_activity"] = activity
                if decision.should_stop and not stop_sim_path.exists():
                    stop_sim_path.write_text(
                        json.dumps({"reason": decision.reason, "t": float(res["t"])}),
                        encoding="utf-8",
                    )
        if res.get("central_radial_block"):
            _append_radial_block(central_radial_out_path, res["central_radial_block"])

    # If rendering frames, clear existing frames for the requested fields/axis at startup.
    frame_fields_startup = [_canonical_field_name(f) for f in args.frames_fields]
    if frame_fields_startup and not args.keep_existing_frames:
        _cleanup_existing_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            fields=frame_fields_startup,
            axis=args.frames_axis,
            verbose=args.verbose,
        )

    projection_fields_startup = [_canonical_field_name(f) for f in args.projection_fields]
    if projection_fields_startup and args.projection_axes and not args.keep_existing_frames:
        _cleanup_projection_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            fields=projection_fields_startup,
            axes=args.projection_axes,
            verbose=args.verbose,
        )

    if args.embedding and not args.keep_existing_frames:
        _cleanup_embedding_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            verbose=args.verbose,
        )

    def process_once() -> int:
        plot_dirs = _iter_plotfile_dirs(data_dir)
        if not plot_dirs:
            return 0
        processed_count = 0

        # Never delete newest keep_last plotfiles. With keep_last=0, protect none.
        keep_last = max(0, int(args.keep_last))
        protected = set(plot_dirs[-keep_last:]) if keep_last > 0 else set()

        # Reclaim already-processed plotfiles before starting another expensive
        # extraction batch. Otherwise GC is delayed by several minutes while
        # the next multi-GB files are read from NFS.
        if args.delete:
            for p in plot_dirs:
                if p in protected:
                    continue
                key = os.path.basename(p)
                if not state.get(key):
                    continue
                if not _is_plotfile_ready(p, stable_seconds=float(args.stable_seconds)):
                    continue
                try:
                    shutil.rmtree(p)
                    if args.verbose:
                        print(f"[gc] deleted previously-processed {key}")
                except Exception as e:
                    print(f"WARNING: failed to delete {p}: {e}")

        to_process = []
        for p in plot_dirs:
            key = os.path.basename(p)
            if state.get(key):
                continue
            if not _is_plotfile_ready(p, stable_seconds=float(args.stable_seconds)):
                continue
            to_process.append(p)

        if not to_process:
            return 0

        # In watch mode, process only one worker-sized batch per pass. This
        # returns to the cleanup pass promptly instead of queueing the entire
        # backlog while multi-GB plotfiles continue accumulating.
        if args.watch:
            to_process = to_process[: max(1, int(args.jobs))]

        print(
            f"Processing {len(to_process)} plotfile(s) "
            f"(~1–5 min each on large AMR; use -j 1 for reliability)...",
            flush=True,
        )

        args_dict = vars(args)
        args_dict["frames_out"] = os.path.abspath(args.frames_out)
        frame_zlims = state.get("frame_zlims", {})
        use_global_zlim = bool(
            args.frames_global_zlim and not _frames_auto_zlim_enabled(args.frames_auto_zlim)
        )
        if frame_fields_startup and use_global_zlim and not frame_zlims and to_process:
            first = min(
                to_process,
                key=lambda p: _parse_plot_index(os.path.basename(p)) if _parse_plot_index(os.path.basename(p)) is not None else 10**12,
            )
            frame_zlims = _lock_frame_zlims_from_plotfile(first, args_dict)
            state["frame_zlims"] = frame_zlims
            _save_state(state_path, state)
        args_dict["frame_zlims"] = frame_zlims
        args_dict["frames_global_zlim"] = use_global_zlim

        if args.jobs > 1:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor, as_completed

            mp_ctx = mp.get_context("spawn")
            with ProcessPoolExecutor(max_workers=args.jobs, mp_context=mp_ctx) as executor:
                futures = {
                    executor.submit(
                        _process_single_plotfile, p, args_dict, protected, processed_count + i
                    ): p
                    for i, p in enumerate(to_process)
                }
                done = 0
                for f in as_completed(futures):
                    p = futures[f]
                    done += 1
                    key = os.path.basename(p)
                    print(f"[{done}/{len(to_process)}] finished {key}", flush=True)
                    try:
                        res = f.result()
                        if res["success"]:
                            if res["psi4_line"]:
                                _append_line(out_path, header=header, line=res["psi4_line"])
                            if res.get("psi4_all_line"):
                                _append_line(
                                    psi4_all_out_path,
                                    header=psi4_all_header,
                                    line=res["psi4_all_line"],
                                )
                            if res.get("psi4_directional_line"):
                                _append_line(
                                    psi4_directional_out_path,
                                    header=psi4_directional_header,
                                    line=res["psi4_directional_line"],
                                )
                            if res["areal_line"]:
                                _append_line(areal_out_path, header=areal_header, line=res["areal_line"])
                            if res["shell_line"]:
                                _append_line(shell_out_path, header=shell_header, line=res["shell_line"])
                            if res.get("boundary_flux_line"):
                                _append_line(
                                    boundary_flux_out_path,
                                    header="time  net_outward_flux  psi4_boundary_amp",
                                    line=res["boundary_flux_line"],
                                )
                            if res.get("ftl_line"):
                                _append_line(
                                    ftl_out_path,
                                    header=FTL_TIMESERIES_HEADER,
                                    line=res["ftl_line"],
                                )
                                if args.ftl_timeseries:
                                    _append_incremental_score(float(res["t"]))
                            if res.get("confinement_line"):
                                _append_line(
                                    confinement_out_path,
                                    header=CONFINEMENT_TIMESERIES_HEADER,
                                    line=res["confinement_line"],
                                )
                            if res.get("sector_barycenters_line"):
                                _append_line(
                                    sector_barycenters_out_path,
                                    header=SECTOR_BARYCENTERS_HEADER,
                                    line=res["sector_barycenters_line"],
                                )
                            _handle_central_outputs(res)

                            state[res["key"]] = True
                            _save_state(state_path, state)
                            processed_count += 1
                            
                            if args.verbose:
                                print(f"[ok] {res['key']}  t={res['t']:.6g}  {res['status_str']}  ({res['dt_s']:.2f}s)")
                        else:
                            print(f"WARNING: failed to process {p}: {res.get('error', 'Unknown error')}")
                    except Exception as e:
                        print(f"WARNING: worker failed for {p}: {e}")
        else:
            for i, p in enumerate(to_process):
                key = os.path.basename(p)
                print(f"[{i + 1}/{len(to_process)}] loading {key} ...", flush=True)
                res = _process_single_plotfile(p, args_dict, protected, processed_count + i)
                if res["success"]:
                    if res["psi4_line"]:
                        _append_line(out_path, header=header, line=res["psi4_line"])
                    if res.get("psi4_all_line"):
                        _append_line(
                            psi4_all_out_path,
                            header=psi4_all_header,
                            line=res["psi4_all_line"],
                        )
                    if res.get("psi4_directional_line"):
                        _append_line(
                            psi4_directional_out_path,
                            header=psi4_directional_header,
                            line=res["psi4_directional_line"],
                        )
                    if res["areal_line"]:
                        _append_line(areal_out_path, header=areal_header, line=res["areal_line"])
                    if res["shell_line"]:
                        _append_line(shell_out_path, header=shell_header, line=res["shell_line"])
                    if res.get("boundary_flux_line"):
                        _append_line(
                            boundary_flux_out_path,
                            header="time  net_outward_flux  psi4_boundary_amp",
                            line=res["boundary_flux_line"],
                        )
                    if res.get("ftl_line"):
                        _append_line(
                            ftl_out_path,
                            header=FTL_TIMESERIES_HEADER,
                            line=res["ftl_line"],
                        )
                        if args.ftl_timeseries:
                            _append_incremental_score(float(res["t"]))
                    if res.get("confinement_line"):
                        _append_line(
                            confinement_out_path,
                            header=CONFINEMENT_TIMESERIES_HEADER,
                            line=res["confinement_line"],
                        )
                    if res.get("sector_barycenters_line"):
                        _append_line(
                            sector_barycenters_out_path,
                            header=SECTOR_BARYCENTERS_HEADER,
                            line=res["sector_barycenters_line"],
                        )
                    _handle_central_outputs(res)

                    state[res["key"]] = True
                    _save_state(state_path, state)
                    processed_count += 1
                    
                    if args.verbose:
                        print(f"[ok] {res['key']}  t={res['t']:.6g}  {res['status_str']}  ({res['dt_s']:.2f}s)")
                else:
                    print(f"WARNING: failed to process {p}: {res.get('error', 'Unknown error')}")

        return processed_count

    if args.watch:
        print(f"Watching {data_dir} for plotfiles. Writing {out_path}")
        while True:
            n = process_once()
            if n:
                if not args.verbose:
                    print(f"Processed {n} plotfile(s).")
            time.sleep(float(args.poll_seconds))
    else:
        n = process_once()
        print(f"Processed {n} plotfile(s). Output: {out_path}")
