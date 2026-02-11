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
    """Return (rank, size, comm) or (0, 1, None) if not using MPI."""
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        return comm.Get_rank(), comm.Get_size(), comm
    except ImportError:
        return 0, 1, None

def create_visualizations():
    # 1. Setup Command Line Arguments
    parser = argparse.ArgumentParser(description="GRTeclyn Visualization Script")
    parser.add_argument("--field", type=str, default="chi", help="Field to plot (chi, K, Theta, lapse, etc.)")
    parser.add_argument("--animate", action="store_true", help="Create an MP4 animation")
    parser.add_argument("--zoom", type=float, default=60.0, help="Width of the plot in code units")
    parser.add_argument("--data", type=str, default=_DEFAULT_DATA, help="Path to plt folders")
    parser.add_argument("--out", type=str, default=_DEFAULT_OUT, help="Output directory")
    parser.add_argument("--mpi", action="store_true", help="Use MPI to parallelize across multiple cores (run with: mpirun -np N python visualize.py ...)")
    
    args = parser.parse_args()

    rank, size, comm = _get_mpi()
    use_mpi = args.mpi and size > 1
    if args.mpi and not use_mpi:
        if size == 1 and comm is None:
            print("WARNING: --mpi given but mpi4py not installed. Run: pip install mpi4py")
        elif size == 1:
            print("WARNING: --mpi given but only 1 MPI rank. Run with: mpirun -np N python visualize.py ... --mpi")

    # 2. Field-specific Configurations (Colormaps and Scales)
    configs = {
        "chi":   {"zlim": (0.1, 1.0),  "cmap": "magma", "label": "Conformal Factor (chi)"},
        "K":     {"zlim": (-0.01, 0.01), "cmap": "RdBu",  "label": "Trace of Extrinsic Curvature (K)"},
        "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": "Z4 Constraint (Theta)"},
        "lapse": {"zlim": (0.1, 1.0),  "cmap": "viridis", "label": "Lapse (alpha)"}
    }

    # Use custom config if available, otherwise use defaults
    cfg = configs.get(args.field, {"zlim": (None, None), "cmap": "viridis", "label": args.field})

    # 3. Find files
    plotfiles = sorted(glob.glob(os.path.join(args.data, "plt*")))
    if not plotfiles:
        print(f"No plotfiles found in {args.data}")
        return

    output_dir = os.path.join(args.out, args.field)
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    # Delete old frames from previous runs (rank 0 only)
    if not use_mpi or rank == 0:
        for old_frame in glob.glob(os.path.join(frames_dir, "frame_*.png")):
            os.remove(old_frame)

    if not args.animate:
        plotfiles = [plotfiles[-1]]

    if not use_mpi or rank == 0:
        print(f"Target Field: {args.field}")
        print(f"Processing {len(plotfiles)} files" + (f" on {size} MPI ranks" if use_mpi else "..."))

    # Field check (rank 0, or all if serial)
    if not use_mpi or rank == 0:
        ds0 = yt.load(plotfiles[0])
        available = [f[1] for f in ds0.field_list]
        if args.field not in available:
            print(f"Error: Field '{args.field}' not found. Available: {available}")
            if use_mpi:
                comm.Abort(1)
            return

    if use_mpi:
        comm.Barrier()

    # 4. Processing Loop (each rank processes its subset of indices)
    indices = range(rank, len(plotfiles), size) if use_mpi else range(len(plotfiles))
    for i in indices:
        plt_path = plotfiles[i]
        ds = yt.load(plt_path)

        # Plane slice at z=0 (mirror)
        plot_center = ds.arr([64, 64, 0], 'code_length')
        slc = yt.SlicePlot(ds, 'z', ('boxlib', args.field), center=plot_center)
        slc.set_width((args.zoom, "code_length"))
        
        # Force linear scale (avoids SymLogNorm failing when field is all zeros, e.g. K at t=0)
        slc.set_log(('boxlib', args.field), False)
        
        # Apply limits and colormap
        if cfg["zlim"][0] is not None:
            slc.set_zlim(('boxlib', args.field), cfg["zlim"][0], cfg["zlim"][1])
        slc.set_cmap(('boxlib', args.field), cfg["cmap"])
        
        # Add a title
        slc.annotate_title(f"Field: {cfg['label']} | Time: {ds.current_time:.2f}")

        frame_path = os.path.join(frames_dir, f"frame_{i:04d}.png")
        slc.save(frame_path)
        
        # Print progress
        ad = ds.all_data()
        stats = ad.quantities.extrema(('boxlib', args.field))
        rank_tag = f" [rank {rank}]" if use_mpi else ""
        print(f"Frame {i:04d} | Min: {float(stats[0]):.4e} | Max: {float(stats[1]):.4e}{rank_tag}")

    if use_mpi:
        comm.Barrier()

    # 5. Create Animation (rank 0 only)
    if args.animate and len(plotfiles) > 1 and (not use_mpi or rank == 0):
        movie_path = os.path.join(output_dir, f"movie_{args.field}.mp4")
        print(f"Stitching {args.field} animation...")
        cmd = (
            f"ffmpeg -y -framerate 10 -i {frames_dir}/frame_%04d.png "
            f'-vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p '
            f"{movie_path}"
        )
        if os.system(cmd) == 0:
            print(f"SUCCESS: {movie_path}")
        else:
            print("ERROR: ffmpeg failed. Ensure ffmpeg is installed and frames were fully generated.")

if __name__ == "__main__":
    create_visualizations()