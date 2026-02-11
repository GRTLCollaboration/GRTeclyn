import yt
import glob
import os

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(os.path.dirname(_SCRIPT_DIR))
_DEFAULT_DATA = os.path.abspath(os.path.join(_PROJECT_ROOT, "..", "data"))
_DEFAULT_OUT = os.path.join(_SCRIPT_DIR, "waves")

def plot_waves(data_dir, output_dir):
    # Find all plotfiles, but skip plt00000 because it usually has no wave data
    plotfiles = sorted(glob.glob(os.path.join(data_dir, "plt*")))
    if len(plotfiles) > 1:
        plotfiles = plotfiles[1:] 
        
    os.makedirs(output_dir, exist_ok=True)
    
    # Standard name for GRTeclyn/AMReX derived fields
    field = ('boxlib', 'Weyl4_Re') 
    
    for i, plt_path in enumerate(plotfiles):
        ds = yt.load(plt_path)
        
        # Check if the field actually exists in this file
        available_fields = [f[1] for f in ds.field_list]
        if 'Weyl4_Re' not in available_fields:
            print(f"!!! Warning: Weyl4_Re not found in {os.path.basename(plt_path)}")
            print(f"Available fields are: {available_fields}")
            continue

        # Slice at z=0 (reflective symmetry mirror)
        slc = yt.SlicePlot(ds, 'z', field, center=[64, 64, 0])
        slc.set_width((120, "code_length"))
        
        # Waves are small (10^-3 to 10^-5). Tighten limits to see them.
        slc.set_zlim(field, -0.0005, 0.0005)
        slc.set_cmap(field, "RdBu") # Red/Blue is best for waves
        
        frame_path = os.path.join(output_dir, f"wave_frame_{i:04d}.png")
        slc.save(frame_path)
        print(f"Successfully processed {os.path.basename(plt_path)}")

if __name__ == "__main__":
    plot_waves(_DEFAULT_DATA, _DEFAULT_OUT)