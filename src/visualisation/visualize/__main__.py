import yt
import numpy as np
import glob
import os
import argparse

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
_DEFAULT_DATA = os.path.abspath(os.path.join(_PROJECT_ROOT, "..", "data"))
_DEFAULT_OUT = _SCRIPT_DIR

def _get_mpi():
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        return comm.Get_rank(), comm.Get_size(), comm
    except ImportError:
        return 0, 1, None

def create_visualizations():
    parser = argparse.ArgumentParser(description="GRTeclyn Visualization Script")
    parser.add_argument("--field", type=str, default="chi")
    parser.add_argument("--axis", type=str, default="z", choices=["x", "y", "z"])
    parser.add_argument("--coord", type=float, default=None, help="Coordinate to slice at")
    parser.add_argument("--animate", action="store_true")
    parser.add_argument("--zoom", type=float, default=None, help="Zoom width in code units (default: full domain)")
    parser.add_argument("--center", type=float, nargs=3, default=None, help="Plot center (x y z). With --corner, this is treated as the symmetry origin (corner).")
    parser.add_argument("--corner", action="store_true", help="Corner mode for symmetry-reduced domains (origin at bottom-left). Requires --zoom for deterministic extents.")
    parser.add_argument("--data", type=str, default=_DEFAULT_DATA)
    parser.add_argument("--out", type=str, default=_DEFAULT_OUT)
    parser.add_argument("--mpi", action="store_true")
    parser.add_argument("--autoscale", action="store_true", help="Ignore the preset color limits and autoscale to the per-frame data range (reveals small deviations).")
    parser.add_argument("--vmin", type=float, default=None, help="Explicit colorbar minimum (overrides preset and --autoscale).")
    parser.add_argument("--vmax", type=float, default=None, help="Explicit colorbar maximum (overrides preset and --autoscale).")

    args = parser.parse_args()

    rank, size, comm = _get_mpi()
    use_mpi = size > 1 or (args.mpi and size > 1)

    configs = {
        "chi":   {"zlim": (0.0, 1.0),  "cmap": "magma", "label": r"Conformal Factor $\chi$"},
        "K":     {"zlim": (-0.1, 0.1), "cmap": "RdBu",  "label": r"Trace of Extrinsic Curvature $K$"},
        "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": r"Z4 Constraint $\Theta$"},
        "lapse": {"zlim": (0.0, 1.0),  "cmap": "viridis", "label": r"Lapse $\alpha$"},
        "Ham":   {"zlim": (-0.1, 0.1), "cmap": "RdBu", "label": r"Hamiltonian Constraint $\mathcal{H}$"},
        "A11":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{11}$"},
        "A12":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{12}$"},
        "A22":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{22}$"},
        "A33":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{33}$"},
        "GW_Plus": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain Proxy $h_+$"},
        "GW_Cross": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain Proxy $h_\times$"},
    }
    cfg = configs.get(args.field, {"zlim": (None, None), "cmap": "viridis", "label": args.field})

    # Match plt* (BinaryBH default) and *Plt* (WormholePlt, etc.)
    # Filter to valid plot roots only (directories with Header, exclude Level_* subdirs)
    raw = set(
        glob.glob(os.path.join(args.data, "plt*")) +
        glob.glob(os.path.join(args.data, "*Plt*"))
    )
    plotfiles = sorted(
        f for f in raw
        if os.path.isdir(f) and os.path.isfile(os.path.join(f, "Header"))
    )
    if not plotfiles:
        print(f"No plotfiles found in {args.data}")
        return

    if not use_mpi or rank == 0:
        print(f"Found {len(plotfiles)} plotfiles")

    output_dir = os.path.join(args.out, f"{args.field}_{args.axis}")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    if not use_mpi or rank == 0:
        for old_frame in glob.glob(os.path.join(frames_dir, "*.png")):
            try:
                os.remove(old_frame)
            except FileNotFoundError:
                pass  # Race with other ranks or already deleted

    if not args.animate:
        plotfiles = [plotfiles[-1]]

    # Define derived fields for Gravitational Wave proxies
    def _gw_plus(field, data):
        # A11 - A22 is roughly h_plus for waves propagating in Z
        return data["boxlib", "A11"] - data["boxlib", "A22"]

    def _gw_cross(field, data):
        # 2 * A12 is roughly h_cross for waves propagating in Z
        return 2.0 * data["boxlib", "A12"]

    if not use_mpi or rank == 0:
        print(f"Target: {args.field} | Axis: {args.axis} | Zoom: {args.zoom}")

    indices = range(rank, len(plotfiles), size) if use_mpi else range(len(plotfiles))
    frame_counter = 0
    for i in indices:
        plt_path = plotfiles[i]
        try:
            ds = yt.load(plt_path)

            # Register derived fields if needed
            if args.field == "GW_Plus":
                ds.add_field(("boxlib", "GW_Plus"), function=_gw_plus, sampling_type="cell", units="")
            if args.field == "GW_Cross":
                ds.add_field(("boxlib", "GW_Cross"), function=_gw_cross, sampling_type="cell", units="")

        except (ValueError, OSError) as e:
            if not use_mpi or rank == 0:
                print(f"Skipping {os.path.basename(plt_path)} (corrupt/incomplete): {e}")
            continue

        # --- UPDATED CENTER LOGIC ---
        # 1. Calculate the physical middle of the box (e.g. 64 for 128 box)
        mid_x = float((ds.domain_right_edge[0] + ds.domain_left_edge[0]) / 2.0)
        mid_y = float((ds.domain_right_edge[1] + ds.domain_left_edge[1]) / 2.0)

        # 2. Physics Center: X,Y are in middle.
        z_center = 0.0

        physics_center = [mid_x, mid_y, z_center]

        # Corner mode: treat --center as the symmetry origin (lower-left corner)
        # and set the plot center to origin + zoom/2 in the slice plane.
        if args.corner and args.zoom is not None:
            slice_plane_val = 0.0 if args.coord is None else float(args.coord)
            origin = np.array(args.center, dtype=float) if args.center is not None else np.array(ds.domain_left_edge, dtype=float)
            zoom = float(args.zoom)
            if args.axis == "z":
                physics_center = [origin[0] + zoom / 2.0, origin[1] + zoom / 2.0, slice_plane_val]
            elif args.axis == "y":
                physics_center = [origin[0] + zoom / 2.0, slice_plane_val, origin[2] + zoom / 2.0]
            elif args.axis == "x":
                physics_center = [slice_plane_val, origin[1] + zoom / 2.0, origin[2] + zoom / 2.0]
        elif args.center is not None:
            physics_center = list(args.center)

        # 3. Adjust based on slicing axis (--coord override)
        if args.axis == 'z':
            if args.coord is not None: physics_center[2] = args.coord
        elif args.axis == 'y':
            if args.coord is not None: physics_center[1] = args.coord
        elif args.axis == 'x':
            if args.coord is not None: physics_center[0] = args.coord

        plot_center = ds.arr(physics_center, 'code_length')

        slc = yt.SlicePlot(ds, args.axis, ('boxlib', args.field), center=plot_center)
        # Use physical (native) dataset coordinates on axes (no auto-centering to [-L/2, L/2]).
        slc.set_origin("native")
        if args.zoom is not None:
            slc.set_width((args.zoom, "code_length"))
        slc.set_log(('boxlib', args.field), False)

        # Calculate and print min/max for the current field
        if args.zoom is not None:
            frb_width = args.zoom
        else:
            frb_width = float(ds.domain_width[0])

        frb = slc.data_source.to_frb((frb_width, "code_length"), 512)
        field_data = frb[('boxlib', args.field)]
        min_val = float(np.nanmin(field_data))
        max_val = float(np.nanmax(field_data))
        if not use_mpi or rank == 0:
            print(f"Frame {frame_counter}: {args.field} min={min_val:.2e}, max={max_val:.2e}")

        if args.vmin is not None or args.vmax is not None:
            lo = args.vmin if args.vmin is not None else min_val
            hi = args.vmax if args.vmax is not None else max_val
            slc.set_zlim(('boxlib', args.field), lo, hi)
        elif args.autoscale:
            # Pad a hair so a near-constant field still renders a usable range.
            if max_val - min_val < 1.0e-12:
                pad = max(abs(max_val), 1.0e-6) * 1.0e-3 + 1.0e-12
                slc.set_zlim(('boxlib', args.field), min_val - pad, max_val + pad)
            else:
                slc.set_zlim(('boxlib', args.field), min_val, max_val)
        elif cfg["zlim"][0] is not None:
            slc.set_zlim(('boxlib', args.field), cfg["zlim"][0], cfg["zlim"][1])
        slc.set_cmap(('boxlib', args.field), cfg["cmap"])

        coord_val = physics_center[{'x':0,'y':1,'z':2}[args.axis]]
        slc.annotate_title(r"%s $\quad t=%.2f \quad %s=%.1f$" % (cfg["label"], float(ds.current_time), args.axis, coord_val))

        # Force pure LaTeX axis labels (no "(code length)" units in label text).
        axis_map = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}
        xlabel_name, ylabel_name = axis_map[args.axis]
        slc.set_xlabel(r"$%s$" % xlabel_name)
        slc.set_ylabel(r"$%s$" % ylabel_name)
        slc.set_colorbar_label(("boxlib", args.field), cfg["label"])

        # Remove the duplicated "0" tick label at the symmetry origin corner.
        # In corner mode the plot lower-left is (0,0); hide the y-axis 0 label
        # so only one "0" is shown in the corner.
        if args.corner:
            slc.render()
            plot = slc.plots[('boxlib', args.field)]
            ax = plot.axes
            ymin = float(ax.get_ylim()[0])
            for tick, lbl in zip(ax.get_yticks(), ax.get_yticklabels()):
                if abs(float(tick) - ymin) < 1.0e-9:
                    lbl.set_visible(False)

        frame_idx = frame_counter if not use_mpi else i
        frame_name = f"frame_{args.axis}_{frame_idx:04d}.png"
        slc.save(os.path.join(frames_dir, frame_name))
        frame_counter += 1

    if use_mpi:
        comm.Barrier()

    if args.animate and len(plotfiles) > 1 and (not use_mpi or rank == 0):
        movie_name = f"movie_{args.field}_{args.axis}.mp4"
        movie_path = os.path.join(output_dir, movie_name)
        print(f"Stitching animation: {movie_name}")
        cmd = (
            f"ffmpeg -y -framerate 10 -i {frames_dir}/frame_{args.axis}_%04d.png "
            f'-vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p '
            f"{movie_path}"
        )
        os.system(cmd)

if __name__ == "__main__":
    create_visualizations()
