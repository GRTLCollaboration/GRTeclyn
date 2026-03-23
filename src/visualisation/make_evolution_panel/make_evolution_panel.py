import argparse
import os
from pathlib import Path

import numpy as np
from PIL import Image

# Default output directory: src/visualisation/plots (sibling of this package)
PLOTS_DIR = Path(__file__).resolve().parent.parent / "plots"

def erase_box(img_arr, x0, x1, y0, y1):
    img_arr = img_arr.copy()
    img_arr[y0:y1, x0:x1] = 255
    return img_arr


def to_grayscale_rgb(img_rgb: np.ndarray) -> np.ndarray:
    """
    Convert an RGB uint8 image to grayscale *but keep 3 channels* (RGB)
    so downstream cropping/pasting logic remains unchanged.
    """
    img = img_rgb.astype(np.float32)
    # ITU-R BT.601 luma transform (standard for sRGB-like imagery)
    y = 0.299 * img[:, :, 0] + 0.587 * img[:, :, 1] + 0.114 * img[:, :, 2]
    y8 = np.clip(y, 0, 255).astype(np.uint8)
    return np.stack([y8, y8, y8], axis=-1)


def draw_grid_lines(img_arr: np.ndarray, x_ticks, y_ticks):
    """Draw a dashed grid on the image array at given tick coordinates."""
    img = img_arr.copy()
    h, w = img.shape[:2]
    # Very light gray, dashed pattern
    dash_len = 10
    space_len = 10
    color = np.array([180, 180, 180], dtype=img.dtype)

    # Draw horizontal grid lines (for each y_tick)
    for y in y_ticks:
        if 0 <= y < h:
            for x in range(0, w):
                # only draw inside the plot area, say after the y-axis
                if (x % (dash_len + space_len)) < dash_len:
                    # Blend with existing color (alpha = 0.5)
                    img[y, x] = img[y, x] * 0.5 + color * 0.5

    # Draw vertical grid lines (for each x_tick)
    for x in x_ticks:
        if 0 <= x < w:
            for y in range(0, h):
                if (y % (dash_len + space_len)) < dash_len:
                    img[y, x] = img[y, x] * 0.5 + color * 0.5

    return img

def create_panel(frame_dir, output_path, frame_indices, grayscale: bool = False, grid: bool = False):
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
            frames.append(np.array(Image.open(path_z).convert("RGB")))
        elif os.path.exists(path_std):
            frames.append(np.array(Image.open(path_std).convert("RGB")))
        else:
            raise FileNotFoundError(f"Could not find frame {f_idx} in {frame_dir}")

    # Layout bounds determined empirically
    Y_top = 10 
    Y_bottom = 1354
    X_plot_left = 130  # Increased from 103 to 130 to completely crop out the y-axis numbers
    X_plot_right = 1303
    X_cbar_right = 1502

    # Y-axis ticks are around x=130-140, they are small black horizontal lines.
    # X-axis ticks are around y=1250-1260, they are small black vertical lines.
    
    crops = []
    for i, frame in enumerate(frames):
        f = frame.copy()
        
        # --- Grid Extraction ---
        if grid:
            # We want FEWER horizontal lines (larger spacing on y-axis)
            # and MORE vertical lines (tighter spacing on x-axis).
            # The tick spacing in pixels is approx 30px (y) and 150px (x).
            # Let's use a uniform 100px grid anchored at the bottom-left origin.
            
            plot_y_min = 40
            plot_y_max = 1255
            plot_x_min = 135
            plot_x_max = 1300
            
            # The bottom-left origin of the plot
            origin_x = 135
            origin_y = 1255
            
            # Generate ticks with 100px spacing
            step = 100
            
            # Y ticks (horizontal lines) - fewer of them (spacing = 100)
            y_ticks = np.arange(origin_y, plot_y_min, -step)
            
            # X ticks (vertical lines) - more of them (spacing = 100)
            x_ticks = np.arange(origin_x, plot_x_max, step)
            
            dash_len = 8
            space_len = 6
            grid_col = np.array([200, 200, 200], dtype=f.dtype)
            
            # Horizontal lines
            for yt in y_ticks:
                if plot_y_min < yt < plot_y_max:
                    for xt in range(plot_x_min, plot_x_max):
                        if (xt % (dash_len + space_len)) < dash_len:
                            # Avoid over-drawing the black data lines. Only draw on white background.
                            if np.all(f[yt, xt] > 230):
                                f[yt, xt] = grid_col

            # Vertical lines
            for xt in x_ticks:
                if plot_x_min < xt < plot_x_max:
                    for yt in range(plot_y_min, plot_y_max):
                        if (yt % (dash_len + space_len)) < dash_len:
                            if np.all(f[yt, xt] > 230):
                                f[yt, xt] = grid_col
        # --- End Grid ---

        if grayscale:
            f = to_grayscale_rgb(f)
        
        # Make the greyish background completely white
        # VisIt sometimes uses 246, but edges might be slightly off.
        # So we threshold anything above 240 as white.
        bg_mask = (f[:, :, 0] > 240) & (f[:, :, 1] > 240) & (f[:, :, 2] > 240)
        f[bg_mask] = [255, 255, 255]
        
        # Erase overlapping x-axis labels ('0' and '40' roughly)
        # The labels can be slightly wider, so we use a wider bounding box:
        # Left label (around x=0): 70 to 160
        # Right label (around x=40): 1250 to 1340
        # Y range for labels: 1256 to 1315
        if i == 0:
            # First panel: erase right label, keep left label and y-axis
            f = erase_box(f, 1250, 1340, 1256, 1315)
            crops.append(Image.fromarray(f).crop((15, Y_top, X_plot_right, Y_bottom)))
        elif i == len(frames) - 1:
            # Last panel: erase left label, keep right label and colorbar
            f = erase_box(f, 70, 160, 1256, 1315)
            crops.append(Image.fromarray(f).crop((X_plot_left, Y_top, X_cbar_right, Y_bottom)))
        else:
            # Middle panels: erase both left and right labels
            f = erase_box(f, 70, 160, 1256, 1315)
            f = erase_box(f, 1250, 1340, 1256, 1315)
            crops.append(Image.fromarray(f).crop((X_plot_left, Y_top, X_plot_right, Y_bottom)))

    total_width = sum(crop.width for crop in crops)
    max_height = max(crop.height for crop in crops)

    panel = Image.new('RGB', (total_width, max_height), (255, 255, 255))
    
    current_x = 0
    for crop in crops:
        panel.paste(crop, (current_x, 0))
        current_x += crop.width

    out_path = Path(output_path)
    
    panel.save(out_path.with_suffix('.png'), dpi=(600, 600))
    panel.save(out_path.with_suffix('.pdf'), dpi=(600, 600), resolution=600)
    print(f"Saved panel plot to {out_path.with_suffix('.png')} and .pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a panel plot from simulation frames.")
    parser.add_argument("--frame_dir", type=str, required=True, help="Directory containing the frames (e.g., K_z/frames)")
    parser.add_argument(
        "--out",
        type=str,
        default="evolution_panel",
        help=f"Output base path (stem; .png/.pdf added). Relative paths are saved under {PLOTS_DIR}.",
    )
    parser.add_argument("--frames", type=int, nargs='+', default=[0, 20, 40, 60], help="Frame indices to plot")
    parser.add_argument(
        "--grayscale",
        action="store_true",
        help="Enable grayscale conversion (default is to keep original colors).",
    )
    parser.add_argument(
        "--no-grid",
        action="store_true",
        help="Disable drawing the dashed grid.",
    )
    
    args = parser.parse_args()

    out = Path(args.out)
    output_path = out if out.is_absolute() else PLOTS_DIR / out
    output_path.parent.mkdir(parents=True, exist_ok=True)

    create_panel(args.frame_dir, str(output_path), args.frames, grayscale=args.grayscale, grid=not args.no_grid)
