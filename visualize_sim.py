import yt
import glob
import os
import sys
import matplotlib.pyplot as plt

def create_visualizations(data_dir, output_dir, animate=False):
    # Find all plt directories
    plotfiles = sorted(glob.glob(os.path.join(data_dir, "plt*")))
    
    if not plotfiles:
        print(f"No plotfiles found in {data_dir}")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Function to get field
    def get_field(ds):
        field_to_plot = ('boxlib', 'chi')
        if field_to_plot not in ds.derived_field_list:
            # Fallback to first field
            if ds.field_list:
                return ds.field_list[0]
            return None
        return field_to_plot

    if not animate:
        # Just visualize the latest
        latest_plt = plotfiles[-1]
        print(f"Visualizing latest plotfile: {latest_plt}")
        ds = yt.load(latest_plt)
        field = get_field(ds)
        if field:
            slc = yt.SlicePlot(ds, 'z', field)
            output_path = os.path.join(output_dir, f"slice_{os.path.basename(latest_plt)}.png")
            slc.save(output_path)
            print(f"Saved plot to {output_path}")
    else:
        # Create an animation
        print(f"Creating animation from {len(plotfiles)} plotfiles...")
        frames_dir = os.path.join(output_dir, "frames")
        if not os.path.exists(frames_dir):
            os.makedirs(frames_dir)
            
        for i, plt_path in enumerate(plotfiles):
            print(f"Processing frame {i+1}/{len(plotfiles)}: {os.path.basename(plt_path)}")
            ds = yt.load(plt_path)
            field = get_field(ds)
            if field:
                slc = yt.SlicePlot(ds, 'z', field)
                # Keep fixed width for animation (yt wants (number, unit) not unyt_quantity)
                w = float((ds.domain_right_edge[0] - ds.domain_left_edge[0]).value)
                slc.set_width((w, "code_length"))
                
                frame_path = os.path.join(frames_dir, f"frame_{i:04d}.png")
                slc.save(frame_path)
        
        # Combine frames using ffmpeg (libx264 needs width/height divisible by 2)
        print("Combining frames into animation...")
        movie_path = os.path.join(output_dir, "simulation_movie.mp4")
        cmd = (
            f"ffmpeg -y -framerate 10 -i {frames_dir}/frame_%04d.png "
            f'-vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p '
            f"{movie_path}"
        )
        ret = os.system(cmd)
        if ret == 0:
            print(f"Animation saved to {movie_path}")
        else:
            print("ffmpeg failed. Frames are in", frames_dir)

if __name__ == "__main__":
    # Paths
    data_path = "/home/jovyan/nachevsky/test/simulation/data"
    output_path = "/home/jovyan/nachevsky/test/simulation/GRTeclyn/Examples/BinaryBH/visualisation"
    
    # Check if user wants animation (processes all plt* and creates simulation_movie.mp4)
    do_animate = "--animate" in sys.argv
    if do_animate:
        print("Animation mode: will process all plotfiles and create simulation_movie.mp4")
    else:
        print("Single-frame mode. Use --animate to create an animation from all plotfiles.")
    
    create_visualizations(data_path, output_path, animate=do_animate)
