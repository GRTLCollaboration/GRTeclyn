from __future__ import annotations

import os
import shutil
import time
from pathlib import Path

import numpy as np

from .extraction.areal import _extract_areal_radius_min
from .extraction.central import _extract_central_timeseries_line
from .extraction.confinement import _extract_confinement_line
from .extraction.ftl import _extract_ftl_timeseries_line
from .extraction.psi4 import _extract_mode_amps_l2m0, _extract_mode_amps_l2_all
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
        "psi4_directional_line": None,
        "areal_line": None,
        "shell_line": None,
        "boundary_flux_line": None,
        "ftl_line": None,
        "confinement_line": None,
        "central_line": None,
        "central_radial_block": None,
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
            radii = [float(r) for r in args_dict["radii"]]
            l2m0_amps, directional = _extract_mode_amps_l2_all(
                ds,
                radii=radii,
                n_points=int(args_dict["n_points"]),
                center=args_dict["center"],
            )
            result["psi4_line"] = f"{t:.16e}  " + "  ".join([f"{a.real:.16e}  {a.imag:.16e}" for a in l2m0_amps])
            # Persist all five (l=2,m) amplitudes: for each radius, write
            # Re/Im for m=-2,-1,0,1,2 (from DirectionalPsi4Metrics).
            if directional.modes and directional.modes[0]:
                parts = [f"{t:.16e}"]
                n_radii = len(directional.modes[0])
                for ir in range(n_radii):
                    for m in (-2, -1, 0, 1, 2):
                        c = directional.modes[m][ir]
                        parts.append(f"{c.real:.16e}")
                        parts.append(f"{c.imag:.16e}")
                result["psi4_all_line"] = "  ".join(parts)
            # Aggregate directional metrics across extraction radii (arithmetic mean).
            if directional.p_total:
                p_total = float(np.mean(directional.p_total))
                p_z_beam = float(np.mean(directional.p_z_beam))
                beam_ratio = float(np.mean(directional.beam_ratio))

                # v5: beaming gain — ratio of peak directional power to isotropic.
                # From l=2 modes: dP/dΩ_max ≈ P_total * beaming_gain / (4π).
                # For l=2 m=±2 dominance (Z-beaming): gain ≈ 5/4 * (1+beam_ratio).
                # More accurately: sum |Σ_m C_m Y_m(θ,φ)|^2 over the sphere grid
                # — but that needs the raw mode amplitudes per radius. Use the
                # proxy: gain = (1 + 4*beam_ratio) which maps beam_ratio 0→1 to
                # gain 1→5 (matches l=2 m=±2 angular pattern).
                beaming_gain = 1.0 + 4.0 * beam_ratio

                # v5: 1/r validity check across extraction radii.
                # If the signal is in the wave zone, r*|Ψ₄| should be roughly
                # constant across radii.  Compute relative std.
                wavezone_std = 0.0
                if len(radii) >= 2 and len(directional.p_total) == len(radii):
                    r_psi4 = np.array([r * np.sqrt(p) for r, p in zip(radii, directional.p_total)])
                    r_psi4_mean = float(np.mean(r_psi4))
                    if r_psi4_mean > 0:
                        wavezone_std = float(np.std(r_psi4) / r_psi4_mean)

                result["psi4_directional_line"] = (
                    f"{t:.16e}  {p_total:.16e}  {p_z_beam:.16e}  {beam_ratio:.16e}"
                    f"  {beaming_gain:.16e}  {wavezone_std:.16e}"
                )

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
                # A field requested for frames may be absent from a given matter
                # model's plotfile (e.g. phi_lump1 only exists for the bicomplex
                # model).  Skip it gracefully so one missing field does not abort
                # the remaining frames AND the post-frame FTL/central extraction.
                try:
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
                except Exception as exc:
                    if args_dict.get("verbose", False):
                        print(f"WARNING: frame field {fld!r} skipped for {key}: {exc}")

        projection_fields = [_canonical_field_name(f) for f in args_dict.get("projection_fields", [])]
        projection_axes = list(args_dict.get("projection_axes", []) or [])
        if projection_fields and projection_axes:
            idx = _parse_plot_index(key)
            frame_idx = idx if idx is not None else fallback_frame_idx
            for fld in projection_fields:
                for axis in projection_axes:
                    try:
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
                    except Exception as exc:
                        if args_dict.get("verbose", False):
                            print(
                                f"WARNING: projection field {fld!r} ({axis}) "
                                f"skipped for {key}: {exc}"
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

        if args_dict.get("confinement_timeseries"):
            result["confinement_line"] = _extract_confinement_line(
                p,
                t=t,
                well_width=float(args_dict.get("confinement_well_width", 1.5)),
                verbose=args_dict.get("verbose", False),
            )

        if args_dict.get("central_timeseries"):
            result["central_line"] = _extract_central_timeseries_line(
                p,
                t=t,
                center=args_dict["center"],
                central_ball=bool(args_dict.get("central_ball")),
                central_ball_radius=args_dict.get("central_ball_radius"),
                verbose=args_dict.get("verbose", False),
            )
        if args_dict.get("central_radial_profile"):
            from .extraction.central_radial import _extract_central_radial_profile_block

            result["central_radial_block"] = _extract_central_radial_profile_block(
                p,
                t=t,
                center=args_dict["center"],
                r_max=float(args_dict.get("central_radial_r_max", 6.0)),
                r_min=args_dict.get("central_ball_radius"),
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
            or result.get("confinement_line")
            or result.get("central_line")
            or result.get("central_radial_block")
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
