import yt
import glob
import os
import sys
import argparse

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(os.path.dirname(_SCRIPT_DIR))
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
    parser.add_argument("--zoom", type=float, default=60.0)
    parser.add_argument("--data", type=str, default=_DEFAULT_DATA)
    parser.add_argument("--out", type=str, default=_DEFAULT_OUT)
    parser.add_argument("--mpi", action="store_true")
    
    args = parser.parse_args()

    rank, size, comm = _get_mpi()
    # Auto-detect MPI when launched via mpirun (size > 1), even without --mpi
    use_mpi = size > 1 or (args.mpi and size > 1)

    configs = {
        "chi":   {"zlim": (0.1, 1.0),  "cmap": "magma", "label": "Conformal Factor (chi)"},
        "K":     {"zlim": (-0.01, 0.01), "cmap": "RdBu",  "label": "Trace of Extrinsic Curvature (K)"},
        "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": "Z4 Constraint (Theta)"},
        "lapse": {"zlim": (0.1, 1.0),  "cmap": "viridis", "label": "Lapse (alpha)"},
        "Weyl4_Re": {"zlim": (-0.001, 0.001), "cmap": "RdBu", "label": "Gravitational Waves (Re Psi4)"}
    }
    cfg = configs.get(args.field, {"zlim": (None, None), "cmap": "viridis", "label": args.field})

    plotfiles = sorted(glob.glob(os.path.join(args.data, "plt*")))
    if not plotfiles:
        print(f"No plotfiles found in {args.data}")
        return

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

    if not use_mpi or rank == 0:
        print(f"Target: {args.field} | Axis: {args.axis} | Zoom: {args.zoom}")

    indices = range(rank, len(plotfiles), size) if use_mpi else range(len(plotfiles))
    for i in indices:
        plt_path = plotfiles[i]
        ds = yt.load(plt_path)

        # --- UPDATED CENTER LOGIC ---
        # 1. Calculate the physical middle of the box (e.g. 64 for 128 box)
        mid_x = float((ds.domain_right_edge[0] + ds.domain_left_edge[0]) / 2.0)
        mid_y = float((ds.domain_right_edge[1] + ds.domain_left_edge[1]) / 2.0)
        
        # 2. Physics Center: X,Y are in middle.
        # For Z: If we are slicing a side view (axis x or y), we must center
        # the camera on the DATA (z > 0), not the mirror (z=0).
        # Reflective symmetry at z=0 means data only exists for z >= 0.
        # Shifting z_center by zoom/2 shows z=0..zoom instead of z=-zoom/2..zoom/2.
        z_center = 0.0
        if args.axis in ['x', 'y']:
            z_center = args.zoom / 2.0
        
        physics_center = [mid_x, mid_y, z_center]

        # 3. Adjust based on slicing axis (--coord override)
        # If we slice Z, we look at Z=coord (or 0). Center is X,Y middle.
        # If we slice Y, we look at Y=coord. Center is X-middle, Z-0.
        if args.axis == 'z':
            if args.coord is not None: physics_center[2] = args.coord
        elif args.axis == 'y':
            if args.coord is not None: physics_center[1] = args.coord
        elif args.axis == 'x':
            if args.coord is not None: physics_center[0] = args.coord

        plot_center = ds.arr(physics_center, 'code_length')
        # ----------------------------

        slc = yt.SlicePlot(ds, args.axis, ('boxlib', args.field), center=plot_center)
        slc.set_width((args.zoom, "code_length"))
        slc.set_log(('boxlib', args.field), False)
        
        if cfg["zlim"][0] is not None:
            slc.set_zlim(('boxlib', args.field), cfg["zlim"][0], cfg["zlim"][1])
        slc.set_cmap(('boxlib', args.field), cfg["cmap"])
        
        coord_val = physics_center[{'x':0,'y':1,'z':2}[args.axis]]
        slc.annotate_title(f"{cfg['label']} | T={ds.current_time:.2f} | {args.axis}={coord_val:.1f}")

        frame_name = f"frame_{args.axis}_{i:04d}.png"
        slc.save(os.path.join(frames_dir, frame_name))

    if use_mpi: comm.Barrier()

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