import argparse
import os
import numpy as np
from PIL import Image

def erase_box(img_arr, x0, x1, y0, y1):
    img_arr = img_arr.copy()
    img_arr[y0:y1, x0:x1] = 255
    return img_arr

def create_panel(frame_dir, output_path, frame_indices):
    """
    Creates a 1xN panel of simulation frames, removing overlapping x-axis labels
    and keeping only one colorbar on the right side.
    """
    frames = []
    for f_idx in frame_indices:
        # Check both frame_z_0000.png and frame_0000.png style
        path_z = os.path.join(frame_dir, f"frame_z_{f_idx:04d}.png")
        path_std = os.path.join(frame_dir, f"frame_{f_idx:04d}.png")
        
        if os.path.exists(path_z):
            frames.append(np.array(Image.open(path_z).convert('RGB')))
        elif os.path.exists(path_std):
            frames.append(np.array(Image.open(path_std).convert('RGB')))
        else:
            raise FileNotFoundError(f"Could not find frame {f_idx} in {frame_dir}")

    # Layout bounds determined empirically
    Y_top = 10 
    Y_bottom = 1354
    X_plot_left = 103
    X_plot_right = 1303
    X_cbar_right = 1502

    # Erase overlapping x-axis labels ('0' and '40' roughly)
    # The labels can be slightly wider, so we use a wider bounding box:
    # Left label (around x=0): 70 to 140
    # Right label (around x=40): 1250 to 1340
    # Y range for labels: 1250 to 1315
    
    # Erase overlapping x-axis labels ('0' and '40' roughly)
    # The labels can be slightly wider, so we use a wider bounding box:
    # Left label (around x=0): 70 to 140
    # Right label (around x=40): 1250 to 1340
    # Y range for labels: 1256 to 1315
    
    crops = []
    for i, frame in enumerate(frames):
        f = frame.copy()
        
        # Make the greyish background completely white
        # VisIt sometimes uses 246, but edges might be slightly off.
        # So we threshold anything above 240 as white.
        bg_mask = (f[:, :, 0] > 240) & (f[:, :, 1] > 240) & (f[:, :, 2] > 240)
        f[bg_mask] = [255, 255, 255]
        
        if i == 0:
            # First panel: erase right label, keep left label and y-axis
            f = erase_box(f, 1250, 1340, 1256, 1315)
            crops.append(Image.fromarray(f).crop((15, Y_top, X_plot_right, Y_bottom)))
        elif i == len(frames) - 1:
            # Last panel: erase left label, keep right label and colorbar
            f = erase_box(f, 70, 140, 1256, 1315)
            crops.append(Image.fromarray(f).crop((X_plot_left, Y_top, X_cbar_right, Y_bottom)))
        else:
            # Middle panels: erase both left and right labels
            f = erase_box(f, 70, 140, 1256, 1315)
            f = erase_box(f, 1250, 1340, 1256, 1315)
            crops.append(Image.fromarray(f).crop((X_plot_left, Y_top, X_plot_right, Y_bottom)))

    total_width = sum(crop.width for crop in crops)
    max_height = max(crop.height for crop in crops)

    panel = Image.new('RGB', (total_width, max_height), (255, 255, 255))
    
    current_x = 0
    for crop in crops:
        panel.paste(crop, (current_x, 0))
        current_x += crop.width

    # Get the base path without extension
    import pathlib
    out_path = pathlib.Path(output_path)
    
    panel.save(out_path.with_suffix('.png'), dpi=(600, 600))
    panel.save(out_path.with_suffix('.pdf'), dpi=(600, 600), resolution=600)
    print(f"Saved panel plot to {out_path.with_suffix('.png')} and .pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a panel plot from simulation frames.")
    parser.add_argument("--frame_dir", type=str, required=True, help="Directory containing the frames (e.g., K_z/frames)")
    parser.add_argument("--out", type=str, required=True, help="Output image path")
    parser.add_argument("--frames", type=int, nargs='+', default=[0, 20, 40, 60], help="Frame indices to plot")
    
    args = parser.parse_args()
    
    create_panel(args.frame_dir, args.out, args.frames)
