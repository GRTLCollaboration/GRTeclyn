from __future__ import annotations

import os
import shutil
import time

from .extraction.areal import _extract_areal_radius_min
from .extraction.ftl import _extract_ftl_timeseries_line
from .extraction.psi4 import _extract_mode_amps_l2m0
from .extraction.shell import _extract_shell_field_stats, _format_shell_stats_line
from .fields import _canonical_field_name
from .frames.embedding import _render_embedding_frame
from .frames.projection import _render_projection_frame
from .frames.slice import _render_slice_frame
from .plotfiles import _parse_plot_index


def _process_single_plotfile(p: str, args_dict: dict, protected: set, fallback_frame_idx: int) -> dict:
    import yt
    
    try:
        yt.funcs.mylog.setLevel(30 if args_dict.get("verbose") else 40)
    except Exception:
        pass

    t0 = time.time()
    result = {
        "p": p,
        "key": os.path.basename(p),
        "t": 0.0,
        "psi4_line": None,
        "areal_line": None,
        "shell_line": None,
        "boundary_flux_line": None,
        "ftl_line": None,
        "success": False,
        "deleted": False,
        "status_str": "",
        "dt_s": 0.0,
        "error": None,
    }

    try:
        ds = yt.load(p)
        t = float(ds.current_time)
        result["t"] = t
        key = result["key"]

        if args_dict.get("boundary_flux", True):
            try:
                from grteclyn_wrapper.metrics.probes.boundary import extract_scalar_boundary_flux

                flux = extract_scalar_boundary_flux(p)
                if flux is not None:
                    net, psi4_amp = flux
                    psi4_val = psi4_amp if psi4_amp is not None else 0.0
                    result["boundary_flux_line"] = f"{t:.16e}  {net:.16e}  {psi4_val:.16e}"
            except Exception as exc:
                if args_dict.get("verbose", False):
                    print(f"WARNING: boundary flux extraction failed for {key}: {exc}")

        if args_dict.get("psi4"):
            if ("boxlib", "Weyl4_Re") not in ds.field_list or ("boxlib", "Weyl4_Im") not in ds.field_list:
                raise RuntimeError("Plotfile missing Weyl4_Re/Im. Set: amr.derive_plot_vars = Weyl4 and re-run.")
            amps = _extract_mode_amps_l2m0(
                ds,
                radii=[float(r) for r in args_dict["radii"]],
                n_points=int(args_dict["n_points"]),
                center=args_dict["center"],
            )
            result["psi4_line"] = f"{t:.16e}  " + "  ".join([f"{a.real:.16e}  {a.imag:.16e}" for a in amps])

        shell_fields = list(args_dict.get("shell_fields") or [])
        if shell_fields:
            try:
                stats = _extract_shell_field_stats(
                    ds,
                    radii=[float(r) for r in args_dict["radii"]],
                    n_points=int(args_dict["n_points"]),
                    center=args_dict["center"],
                    fields=shell_fields,
                )
                result["shell_line"] = _format_shell_stats_line(
                    t,
                    stats,
                    [float(r) for r in args_dict["radii"]],
                    shell_fields,
                )
            except Exception as exc:
                if args_dict.get("verbose", False):
                    print(f"WARNING: shell extraction failed for {key}: {exc}")

        frame_fields = [_canonical_field_name(f) for f in args_dict.get("frames_fields", [])]
        if frame_fields:
            idx = _parse_plot_index(key)
            frame_idx = idx if idx is not None else fallback_frame_idx
            for fld in frame_fields:
                _render_slice_frame(
                    ds,
                    field=fld,
                    axis=args_dict["frames_axis"],
                    coord=args_dict.get("frames_coord"),
                    zoom=args_dict.get("frames_zoom"),
                    center_xyz=args_dict.get("frames_center"),
                    corner=args_dict.get("frames_corner"),
                    frames_out_dir=args_dict["frames_out"],
                    frame_idx=int(frame_idx),
                    verbose=args_dict.get("verbose", False),
                    auto_zlim=args_dict.get("frames_auto_zlim"),
                    frame_zlims=args_dict.get("frame_zlims"),
                    use_global_zlim=args_dict.get("frames_global_zlim", True),
                )

        projection_fields = [_canonical_field_name(f) for f in args_dict.get("projection_fields", [])]
        projection_axes = list(args_dict.get("projection_axes", []) or [])
        if projection_fields and projection_axes:
            idx = _parse_plot_index(key)
            frame_idx = idx if idx is not None else fallback_frame_idx
            for fld in projection_fields:
                for axis in projection_axes:
                    _render_projection_frame(
                        ds,
                        field=fld,
                        axis=axis,
                        method=args_dict.get("projection_method", "mip"),
                        zoom=args_dict.get("frames_zoom"),
                        center_xyz=args_dict.get("frames_center"),
                        frames_out_dir=args_dict["frames_out"],
                        frame_idx=int(frame_idx),
                        verbose=args_dict.get("verbose", False),
                    )

        if args_dict.get("areal_radius"):
            if ("boxlib", "chi") in ds.field_list:
                try:
                    R_min, r_min = _extract_areal_radius_min(ds, center=args_dict["center"])
                    result["areal_line"] = f"{t:.16e}  {R_min:.16e}  {r_min:.16e}"
                except Exception as exc:
                    if args_dict.get("verbose", False):
                        print(f"WARNING: areal extraction failed for {key}: {exc}")
            elif args_dict.get("verbose", False):
                print(f"WARNING: plotfile {key} missing chi field; skipping areal radius.")

        if args_dict.get("embedding"):
            if ("boxlib", "chi") in ds.field_list:
                e_idx = _parse_plot_index(key)
                embed_frame_idx = e_idx if e_idx is not None else fallback_frame_idx
                _render_embedding_frame(
                    ds,
                    frames_out_dir=args_dict["frames_out"],
                    frame_idx=int(embed_frame_idx),
                    center=args_dict["center"],
                    r_max=args_dict.get("embedding_rmax"),
                    verbose=args_dict.get("verbose", False),
                )
            else:
                if args_dict.get("verbose", False):
                    print(f"WARNING: plotfile {key} missing chi field; skipping embedding.")

        if args_dict.get("ftl_timeseries"):
            result["ftl_line"] = _extract_ftl_timeseries_line(
                p,
                t=t,
                center=args_dict["center"],
                ftl_L=args_dict.get("ftl_l"),
                verbose=args_dict.get("verbose", False),
            )

        if args_dict.get("metric_stack_cache"):
            try:
                from grteclyn_wrapper.metrics.probes.ftl.metric_stack_cache import (
                    DEFAULT_N_SPACE,
                    append_slice_from_plotfile,
                    metric_stack_dir,
                )

                cache_dir = metric_stack_dir(Path(args_dict["out"]))
                append_slice_from_plotfile(
                    p,
                    cache_dir,
                    t=t,
                    n_space=int(args_dict.get("metric_stack_n_space", DEFAULT_N_SPACE)),
                    half_width=args_dict.get("ftl_l"),
                )
                result["metric_stack_cached"] = True
            except Exception as exc:
                result["metric_stack_cached"] = False
                if args_dict.get("verbose", False):
                    print(f"WARNING: metric stack cache failed for {key}: {exc}")

        result["success"] = bool(
            result["psi4_line"]
            or result["areal_line"]
            or result["shell_line"]
            or result["ftl_line"]
            or frame_fields
            or projection_fields
        )

        if args_dict.get("delete") and (p not in protected):
            shutil.rmtree(p)
            result["deleted"] = True

    except Exception as e:
        result["error"] = str(e)

    result["dt_s"] = time.time() - t0
    status = "deleted" if result["deleted"] else ("kept" if args_dict.get("delete") else "kept(no-delete)")
    result["status_str"] = status
    
    return result
