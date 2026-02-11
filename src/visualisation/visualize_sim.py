import yt
import glob
import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(os.path.dirname(_SCRIPT_DIR))
_DEFAULT_DATA = os.path.abspath(os.path.join(_PROJECT_ROOT, "..", "data"))
_DEFAULT_OUT = os.path.join(_SCRIPT_DIR, "binary_bh")

def create_visualizations(data_dir, output_dir, animate=False):
    # Find all plt directories
    plotfiles = sorted(glob.glob(os.path.join(data_dir, "plt*")))
    if not plotfiles:
        print(f"No plotfiles found in {data_dir}")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    field_to_plot = ('boxlib', 'chi')
    zoom_width = 60.0 

    # Prepare frames directory
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    if not animate:
        plotfiles = [plotfiles[-1]]

    print(f"Processing {len(plotfiles)} frames...")

    for i, plt_path in enumerate(plotfiles):
        ds = yt.load(plt_path)
        
        # In a 128x128x64 box (reflective symmetry), the center of the orbit is:
        # x=64, y=64, z=0
        plot_center = ds.arr([64, 64, 0], 'code_length')
        
        # Create the slice exactly at the z-mirror (where BHs live)
        slc = yt.SlicePlot(ds, 'z', field_to_plot, center=plot_center)
        slc.set_width((zoom_width, "code_length"))
        
        # Color bar: Chi starts at 1.0. As the BHs form, it drops.
        # We set 0.5 to 1.0 to see the "dents" even if the simulation is early.
        slc.set_zlim(field_to_plot, 0.5, 1.0)
        slc.set_cmap(field_to_plot, "magma")
        
        # Debug info: This tells you if the black holes are actually being resolved
        ad = ds.all_data()
        min_chi = float(ad.quantities.extrema(field_to_plot)[0])
        print(f"Frame {i:04d} ({os.path.basename(plt_path)}): Min Chi = {min_chi:.4f}")

        frame_path = os.path.join(frames_dir, f"frame_{i:04d}.png")
        slc.save(frame_path)

    # --- THIS PART CREATES THE MP4 ---
    if animate and len(plotfiles) > 1:
        print("Stitching frames into animation...")
        movie_path = os.path.join(output_dir, "binary_bh_movie.mp4")
        
        # Command to combine PNGs into a high-quality MP4
        cmd = (
            f"ffmpeg -y -framerate 10 -i {frames_dir}/frame_%04d.png "
            f'-vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" ' # Required for H.264
            f"-c:v libx264 -pix_fmt yuv420p {movie_path}"
        )
        
        ret = os.system(cmd)
        if ret == 0:
            print(f"\nSUCCESS! Movie saved to: {movie_path}")
        else:
            print("\nError: ffmpeg failed. Check if ffmpeg is installed.")

if __name__ == "__main__":
    create_visualizations(_DEFAULT_DATA, _DEFAULT_OUT, animate=True)