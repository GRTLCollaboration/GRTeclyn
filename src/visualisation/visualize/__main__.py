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
    parser.add_argument("--center", type=float, nargs=3, default=None, help="Center of the plot (x y z)")
    parser.add_argument("--corner", action="store_true", help="Place the symmetry origin (0,0) at the bottom-left corner")
    parser.add_argument("--data", type=str, default=_DEFAULT_DATA)
    parser.add_argument("--out", type=str, default=_DEFAULT_OUT)
    parser.add_argument("--mpi", action="store_true")

    args = parser.parse_args()

    rank, size, comm = _get_mpi()
    # Auto-detect MPI when launched via mpirun (size > 1), even without --mpi
    use_mpi = size > 1 or (args.mpi and size > 1)

    configs = {
        "chi":   {"zlim": (0.0, 1.0),  "cmap": "magma", "label": "Conformal Factor (chi)"},
        "K":     {"zlim": (-0.1, 0.1), "cmap": "RdBu",  "label": "Trace of Extrinsic Curvature (K)"},
        "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": "Z4 Constraint (Theta)"},
        "lapse": {"zlim": (0.0, 1.0),  "cmap": "viridis", "label": "Lapse (alpha)"},
        "Ham":   {"zlim": (-0.1, 0.1), "cmap": "RdBu", "label": "Hamiltonian Constraint"},
        "A11":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": "Extrinsic Curvature A11 (hxx)"},
        "A12":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": "Extrinsic Curvature A12 (hxy)"},
        "A22":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": "Extrinsic Curvature A22 (hyy)"},
        "A33":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": "Extrinsic Curvature A33 (hzz)"},
        "GW_Plus": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": "GW Strain Proxy (+): A11 - A22"},
        "GW_Cross": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": "GW Strain Proxy (x): 2 * A12"},
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

        if args.corner and args.zoom is not None:
            slice_plane_val = 0.0
            if args.coord is not None:
                slice_plane_val = args.coord

            if args.axis == 'z':
                physics_center = [args.zoom / 2.0, args.zoom / 2.0, slice_plane_val]
            elif args.axis == 'y':
                physics_center = [args.zoom / 2.0, slice_plane_val, args.zoom / 2.0]
            elif args.axis == 'x':
                physics_center = [slice_plane_val, args.zoom / 2.0, args.zoom / 2.0]

        elif args.center is not None:
            physics_center = args.center

        # 3. Adjust based on slicing axis (--coord override)
        if args.axis == 'z':
            if args.coord is not None: physics_center[2] = args.coord
        elif args.axis == 'y':
            if args.coord is not None: physics_center[1] = args.coord
        elif args.axis == 'x':
            if args.coord is not None: physics_center[0] = args.coord

        plot_center = ds.arr(physics_center, 'code_length')
        # ----------------------------

        slc = yt.SlicePlot(ds, args.axis, ('boxlib', args.field), center=plot_center)
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

        if cfg["zlim"][0] is not None:
            slc.set_zlim(('boxlib', args.field), cfg["zlim"][0], cfg["zlim"][1])
        slc.set_cmap(('boxlib', args.field), cfg["cmap"])

        coord_val = physics_center[{'x':0,'y':1,'z':2}[args.axis]]
        slc.annotate_title(f"{cfg['label']} | T={ds.current_time:.2f} | {args.axis}={coord_val:.1f}")

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
