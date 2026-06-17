from __future__ import annotations

import argparse
import os
import shutil
import time
from pathlib import Path

import yt

from .config import _default_data_dir, _default_frames_out_dir, _frames_auto_zlim_enabled
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
    parser.add_argument("--stable-seconds", type=float, default=5.0, help="Require Header mtime older than this")
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
        choices=["weighted", "ftl_first", "robust_ftl", "general_ftl"],
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
    areal_out_path = out_dir / "areal_radius.dat"
    shell_out_path = out_dir / "shell_profiles.dat"
    boundary_flux_out_path = out_dir / "boundary_flux.dat"
    ftl_out_path = out_dir / "ftl_timeseries.dat"
    score_ts_path = out_dir / "score_timeseries.jsonl"
    header = "# time  " + "  ".join([f"Re(R={R:g})  Im(R={R:g})" for R in args.radii])
    areal_header = "# time  R_areal_min  r_at_R_areal_min"
    shell_header = _shell_stats_header(args.radii, args.shell_fields)

    state = _load_state(state_path)

    # Auto-reset on simulation restart (same output dir reused):
    plot_dirs_now = _iter_plotfile_dirs(data_dir)
    if _should_auto_reset(plot_dirs_now, state):
        print("Detected a likely simulation restart in the same output directory.")
        print(f"Resetting: {out_path} and {state_path}")
        _truncate_if_exists(out_path)
        _truncate_if_exists(areal_out_path)
        if args.shell_fields:
            _truncate_if_exists(shell_out_path)
        if args.ftl_timeseries:
            _truncate_if_exists(ftl_out_path)
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
    if args.incremental_score and args.ftl_timeseries:
        incremental_writer = IncrementalScoreWriter(
            Path(data_dir),
            objective_mode=args.objective_mode,
            target_stop_time=args.target_stop_time,
            ftl_L=args.ftl_l,
            score_weights=score_weights or None,
            evolving_geodesic_mode=bool(args.evolving_geodesic),
            out_path=score_ts_path,
        )

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
            from concurrent.futures import ProcessPoolExecutor, as_completed

            with ProcessPoolExecutor(max_workers=args.jobs) as executor:
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
                                _append_incremental_score(float(res["t"]))

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
                        _append_incremental_score(float(res["t"]))

                    state[res["key"]] = True
                    _save_state(state_path, state)
                    processed_count += 1
                    
                    if args.verbose:
                        print(f"[ok] {res['key']}  t={res['t']:.6g}  {res['status_str']}  ({res['dt_s']:.2f}s)")
                else:
                    print(f"WARNING: failed to process {p}: {res.get('error', 'Unknown error')}")

        # Cleanup pass: delete plotfiles that were processed earlier but were
        # inside keep-last at the time. Once they are no longer protected, we
        # can safely delete them.
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
