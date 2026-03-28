import argparse
import os
from pathlib import Path

import numpy as np
from PIL import Image

# Default output directory: src/visualisation/plots (sibling of this package)
PLOTS_DIR = Path(__file__).resolve().parent.parent / "plots"

KZ_MODE = "k_z"
EMBEDDING_MODE = "embedding"


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


def _rgb_to_hsv(rgb01: np.ndarray) -> np.ndarray:
    """Vectorized RGB->HSV. Input/output are float arrays in [0,1]."""
    r, g, b = rgb01[..., 0], rgb01[..., 1], rgb01[..., 2]
    cmax = np.maximum(np.maximum(r, g), b)
    cmin = np.minimum(np.minimum(r, g), b)
    delta = cmax - cmin

    h = np.zeros_like(cmax)
    nz = delta > 1e-12
    # Hue
    mask = nz & (cmax == r)
    h[mask] = ((g[mask] - b[mask]) / delta[mask]) % 6.0
    mask = nz & (cmax == g)
    h[mask] = ((b[mask] - r[mask]) / delta[mask]) + 2.0
    mask = nz & (cmax == b)
    h[mask] = ((r[mask] - g[mask]) / delta[mask]) + 4.0
    h = (h / 6.0) % 1.0

    # Saturation
    s = np.zeros_like(cmax)
    s[cmax > 1e-12] = delta[cmax > 1e-12] / cmax[cmax > 1e-12]

    v = cmax
    return np.stack([h, s, v], axis=-1)


def _hsv_to_rgb(hsv01: np.ndarray) -> np.ndarray:
    """Vectorized HSV->RGB. Input/output are float arrays in [0,1]."""
    h, s, v = hsv01[..., 0], hsv01[..., 1], hsv01[..., 2]
    h6 = (h % 1.0) * 6.0
    i = np.floor(h6).astype(int)
    f = h6 - i
    p = v * (1.0 - s)
    q = v * (1.0 - f * s)
    t = v * (1.0 - (1.0 - f) * s)

    i_mod = i % 6
    r = np.empty_like(v)
    g = np.empty_like(v)
    b = np.empty_like(v)

    m0 = i_mod == 0
    r[m0], g[m0], b[m0] = v[m0], t[m0], p[m0]
    m1 = i_mod == 1
    r[m1], g[m1], b[m1] = q[m1], v[m1], p[m1]
    m2 = i_mod == 2
    r[m2], g[m2], b[m2] = p[m2], v[m2], t[m2]
    m3 = i_mod == 3
    r[m3], g[m3], b[m3] = p[m3], q[m3], v[m3]
    m4 = i_mod == 4
    r[m4], g[m4], b[m4] = t[m4], p[m4], v[m4]
    m5 = i_mod == 5
    r[m5], g[m5], b[m5] = v[m5], p[m5], q[m5]

    return np.stack([r, g, b], axis=-1)


def _apply_embedding_color_tweak(
    f_rgb: np.ndarray,
    *,
    bg_mask: np.ndarray,
    hue_shift_deg: float = 0.0,
    sat_scale: float = 1.0,
    val_scale: float = 1.0,
    green_scale: float = 1.0,
    contrast: float = 1.0,
    gamma: float = 1.0,
    sat_min: float = 0.12,
    neutral_delta_max: float = 0.05,
    green_dom_min_255: float = 8.0,
) -> np.ndarray:
    """
    Post-process already-rendered embedding PNGs (RGB) to improve aesthetics
    without regenerating frames. Operates only on non-background pixels.
    """
    if (
        abs(hue_shift_deg) < 1e-9
        and abs(sat_scale - 1.0) < 1e-9
        and abs(val_scale - 1.0) < 1e-9
        and abs(green_scale - 1.0) < 1e-9
    ):
        return f_rgb

    out = f_rgb.copy()
    fg = ~bg_mask

    rgb01 = (out.astype(np.float32) / 255.0)
    hsv = _rgb_to_hsv(rgb01)
    sat = hsv[..., 1]
    delta = np.max(rgb01, axis=2) - np.min(rgb01, axis=2)
    # Protect near-neutral pixels (grid, axes, text) from getting tinted.
    non_neutral = (sat >= float(sat_min)) & (delta >= float(neutral_delta_max))
    apply_mask = fg & non_neutral

    # HSV adjustments, applied only to sufficiently-saturated pixels.
    if abs(hue_shift_deg) > 1e-9 or abs(sat_scale - 1.0) > 1e-9 or abs(val_scale - 1.0) > 1e-9:
        hsv2 = hsv.copy()
        if abs(hue_shift_deg) > 1e-9:
            hsv2[..., 0] = (hsv2[..., 0] + (hue_shift_deg / 360.0)) % 1.0
        if abs(sat_scale - 1.0) > 1e-9:
            hsv2[..., 1] = np.clip(hsv2[..., 1] * float(sat_scale), 0.0, 1.0)
        if abs(val_scale - 1.0) > 1e-9:
            hsv2[..., 2] = np.clip(hsv2[..., 2] * float(val_scale), 0.0, 1.0)
        rgb_new = _hsv_to_rgb(hsv2)
        rgb8 = np.clip(np.round(rgb_new * 255.0), 0, 255).astype(np.uint8)
        out[apply_mask] = rgb8[apply_mask]

    # Green-channel attenuation, applied selectively to green-dominant pixels.
    if abs(green_scale - 1.0) > 1e-9:
        r = out[:, :, 0].astype(np.float32)
        g = out[:, :, 1].astype(np.float32)
        b = out[:, :, 2].astype(np.float32)
        green_dom = (g - np.maximum(r, b)) >= float(green_dom_min_255)
        m = apply_mask & green_dom
        g[m] = np.clip(g[m] * float(green_scale), 0, 255)
        out[:, :, 1] = g.astype(np.uint8)

    # Gentle RGB-space contrast / gamma, applied only to the foreground.
    if abs(contrast - 1.0) > 1e-9 or abs(gamma - 1.0) > 1e-9:
        rgb01 = (out.astype(np.float32) / 255.0)
        if abs(contrast - 1.0) > 1e-9:
            rgb01 = (rgb01 - 0.5) * float(contrast) + 0.5
        rgb01 = np.clip(rgb01, 0.0, 1.0)
        if abs(gamma - 1.0) > 1e-9 and gamma > 1e-6:
            rgb01 = rgb01 ** (1.0 / float(gamma))
        rgb8 = np.clip(np.round(rgb01 * 255.0), 0, 255).astype(np.uint8)
        out[fg] = rgb8[fg]

    return out


def _find_title_gap_end_y(img_rgb: np.ndarray, *, white_thr: int = 245, search_max_y: int = 260) -> int:
    """
    For "embedding" frames: find the end of the big whitespace gap between the title
    and the 3D plot area. Returns a y index suitable as a crop start.
    """
    h = img_rgb.shape[0]
    max_y = min(h, search_max_y)
    row_nonwhite = np.sum(np.any(img_rgb[:max_y] < white_thr, axis=2), axis=1)
    zero = row_nonwhite == 0

    best_start, best_end = 0, 0
    cur_start = None
    for i, is_zero in enumerate(zero):
        if is_zero and cur_start is None:
            cur_start = i
        elif (not is_zero) and cur_start is not None:
            if (i - cur_start) > (best_end - best_start):
                best_start, best_end = cur_start, i
            cur_start = None
    if cur_start is not None:
        i = len(zero)
        if (i - cur_start) > (best_end - best_start):
            best_start, best_end = cur_start, i

    # If we didn't find a meaningful gap, don't crop the top.
    if (best_end - best_start) < 10:
        return 0
    return int(best_end)


def _union_content_bbox(
    images_rgb: list[np.ndarray],
    *,
    white_thr: int = 245,
    ignore_top: int = 0,
    pad: int = 6,
) -> tuple[int, int, int, int]:
    """
    Return a (left, top, right, bottom) bbox that encloses non-white content
    across all images, with padding and optional ignore_top.
    """
    if not images_rgb:
        return (0, 0, 0, 0)

    h, w = images_rgb[0].shape[:2]
    x0, y0, x1, y1 = w, h, 0, 0

    for img in images_rgb:
        roi = img[ignore_top:, :, :]
        mask = np.any(roi < white_thr, axis=2)
        ys, xs = np.where(mask)
        if xs.size == 0:
            continue
        x0 = min(x0, int(xs.min()))
        x1 = max(x1, int(xs.max() + 1))
        y0 = min(y0, int(ys.min() + ignore_top))
        y1 = max(y1, int(ys.max() + 1 + ignore_top))

    if x1 <= x0 or y1 <= y0:
        return (0, 0, w, h)

    x0 = max(0, x0 - pad)
    y0 = max(0, y0 - pad)
    x1 = min(w, x1 + pad)
    y1 = min(h, y1 + pad)
    return (x0, y0, x1, y1)

def create_panel(
    frame_dir,
    output_path,
    frame_indices,
    grayscale: bool = False,
    grid: bool = False,
    mode: str = KZ_MODE,
    keep_title: bool = False,
    embedding_hue_shift_deg: float = 0.0,
    embedding_sat_scale: float = 1.0,
    embedding_val_scale: float = 1.0,
    embedding_green_scale: float = 1.0,
    embedding_contrast: float = 1.0,
    embedding_gamma: float = 1.0,
    embedding_filter: str = "none",
):
    """
    Creates a 1xN panel of simulation frames, removing overlapping x-axis labels
    and keeping only one colorbar on the right side.
    """
    if mode not in {KZ_MODE, EMBEDDING_MODE}:
        raise ValueError(f"Unknown mode '{mode}'. Expected '{KZ_MODE}' or '{EMBEDDING_MODE}'.")

    frames: list[np.ndarray] = []
    for f_idx in frame_indices:
        path_z = os.path.join(frame_dir, f"frame_z_{f_idx:04d}.png")
        path_std = os.path.join(frame_dir, f"frame_{f_idx:04d}.png")

        # Prefer different filename conventions depending on the mode.
        candidates = [path_z, path_std] if mode == KZ_MODE else [path_std, path_z]
        chosen = next((p for p in candidates if os.path.exists(p)), None)
        if chosen is None:
            raise FileNotFoundError(f"Could not find frame {f_idx} in {frame_dir}")
        frames.append(np.array(Image.open(chosen).convert("RGB")))

    # Layout bounds determined empirically
    if mode == KZ_MODE:
        Y_top = 10
        Y_bottom = 1354
        X_plot_left = 130  # Increased from 103 to 130 to completely crop out the y-axis numbers
        X_plot_right = 1303
        X_cbar_right = 1502

    # Y-axis ticks are around x=130-140, they are small black horizontal lines.
    # X-axis ticks are around y=1250-1260, they are small black vertical lines.
    
    processed_frames: list[np.ndarray] = []
    for frame in frames:
        f = frame.copy()

        # --- Grid Extraction (K_z only) ---
        if mode == KZ_MODE and grid:
            plot_y_min = 40
            plot_y_max = 1255
            plot_x_min = 135
            plot_x_max = 1300

            origin_x = 135
            origin_y = 1255
            step = 100

            y_ticks = np.arange(origin_y, plot_y_min, -step)
            x_ticks = np.arange(origin_x, plot_x_max, step)

            dash_len = 8
            space_len = 6
            grid_col = np.array([200, 200, 200], dtype=f.dtype)

            for yt in y_ticks:
                if plot_y_min < yt < plot_y_max:
                    for xt in range(plot_x_min, plot_x_max):
                        if (xt % (dash_len + space_len)) < dash_len:
                            if np.all(f[yt, xt] > 230):
                                f[yt, xt] = grid_col

            for xt in x_ticks:
                if plot_x_min < xt < plot_x_max:
                    for yt in range(plot_y_min, plot_y_max):
                        if (yt % (dash_len + space_len)) < dash_len:
                            if np.all(f[yt, xt] > 230):
                                f[yt, xt] = grid_col
        # --- End Grid ---

        if grayscale:
            f = to_grayscale_rgb(f)

        bg_mask = (f[:, :, 0] > 240) & (f[:, :, 1] > 240) & (f[:, :, 2] > 240)
        f[bg_mask] = [255, 255, 255]

        if mode == EMBEDDING_MODE and not grayscale:
            # Presets (can still be overridden by explicit numeric flags)
            if embedding_filter == "less-green":
                embedding_green_scale = embedding_green_scale if embedding_green_scale != 1.0 else 0.85
                embedding_sat_scale = embedding_sat_scale if embedding_sat_scale != 1.0 else 0.95
                embedding_val_scale = embedding_val_scale if embedding_val_scale != 1.0 else 1.02
                embedding_contrast = embedding_contrast if embedding_contrast != 1.0 else 1.03
            elif embedding_filter == "warm":
                embedding_hue_shift_deg = embedding_hue_shift_deg if embedding_hue_shift_deg != 0.0 else -12.0
                embedding_sat_scale = embedding_sat_scale if embedding_sat_scale != 1.0 else 0.95
                embedding_val_scale = embedding_val_scale if embedding_val_scale != 1.0 else 1.02
                embedding_contrast = embedding_contrast if embedding_contrast != 1.0 else 1.03
            elif embedding_filter == "cool":
                embedding_hue_shift_deg = embedding_hue_shift_deg if embedding_hue_shift_deg != 0.0 else 12.0
                embedding_sat_scale = embedding_sat_scale if embedding_sat_scale != 1.0 else 0.95
                embedding_val_scale = embedding_val_scale if embedding_val_scale != 1.0 else 1.02
                embedding_contrast = embedding_contrast if embedding_contrast != 1.0 else 1.03

            f = _apply_embedding_color_tweak(
                f,
                bg_mask=bg_mask,
                hue_shift_deg=embedding_hue_shift_deg,
                sat_scale=embedding_sat_scale,
                val_scale=embedding_val_scale,
                green_scale=embedding_green_scale,
                contrast=embedding_contrast,
                gamma=embedding_gamma,
            )
        processed_frames.append(f)

    crops: list[Image.Image] = []
    if mode == KZ_MODE:
        for i, f in enumerate(processed_frames):
            # Erase overlapping x-axis labels ('0' and '40' roughly)
            if i == 0:
                f = erase_box(f, 1250, 1340, 1256, 1315)
                crops.append(Image.fromarray(f).crop((15, Y_top, X_plot_right, Y_bottom)))
            elif i == len(processed_frames) - 1:
                f = erase_box(f, 70, 160, 1256, 1315)
                crops.append(Image.fromarray(f).crop((X_plot_left, Y_top, X_cbar_right, Y_bottom)))
            else:
                f = erase_box(f, 70, 160, 1256, 1315)
                f = erase_box(f, 1250, 1340, 1256, 1315)
                crops.append(Image.fromarray(f).crop((X_plot_left, Y_top, X_plot_right, Y_bottom)))
    else:
        ignore_top = 0
        if not keep_title and processed_frames:
            ignore_top = _find_title_gap_end_y(processed_frames[0])
        bbox = _union_content_bbox(processed_frames, ignore_top=ignore_top, white_thr=245, pad=6)
        crops = [Image.fromarray(f).crop(bbox) for f in processed_frames]

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
        "--mode",
        type=str,
        choices=[KZ_MODE, EMBEDDING_MODE],
        default=KZ_MODE,
        help=f"Frame type/layout. '{KZ_MODE}' keeps the current K_z crop logic; '{EMBEDDING_MODE}' crops embedding frames.",
    )
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
    parser.add_argument(
        "--keep-title",
        action="store_true",
        help="(Embedding mode) Keep the title area at the top instead of cropping it out.",
    )
    parser.add_argument(
        "--embedding-filter",
        type=str,
        choices=["none", "less-green", "warm", "cool"],
        default="none",
        help="(Embedding mode) Simple color-grading preset applied while stitching.",
    )
    parser.add_argument(
        "--embedding-hue-shift-deg",
        type=float,
        default=0.0,
        help="(Embedding mode) Hue shift in degrees applied while stitching (no frame regeneration).",
    )
    parser.add_argument(
        "--embedding-sat-scale",
        type=float,
        default=1.0,
        help="(Embedding mode) Saturation multiplier applied while stitching (e.g. 0.9).",
    )
    parser.add_argument(
        "--embedding-val-scale",
        type=float,
        default=1.0,
        help="(Embedding mode) Value/brightness multiplier applied while stitching (e.g. 1.05).",
    )
    parser.add_argument(
        "--embedding-green-scale",
        type=float,
        default=1.0,
        help="(Embedding mode) Multiply the green channel to reduce green dominance (e.g. 0.75).",
    )
    parser.add_argument(
        "--embedding-contrast",
        type=float,
        default=1.0,
        help="(Embedding mode) RGB contrast multiplier (e.g. 1.05).",
    )
    parser.add_argument(
        "--embedding-gamma",
        type=float,
        default=1.0,
        help="(Embedding mode) Gamma correction (e.g. 1.05).",
    )
    
    args = parser.parse_args()

    out = Path(args.out)
    output_path = out if out.is_absolute() else PLOTS_DIR / out
    output_path.parent.mkdir(parents=True, exist_ok=True)

    create_panel(
        args.frame_dir,
        str(output_path),
        args.frames,
        grayscale=args.grayscale,
        grid=not args.no_grid,
        mode=args.mode,
        keep_title=args.keep_title,
        embedding_filter=args.embedding_filter,
        embedding_hue_shift_deg=args.embedding_hue_shift_deg,
        embedding_sat_scale=args.embedding_sat_scale,
        embedding_val_scale=args.embedding_val_scale,
        embedding_green_scale=args.embedding_green_scale,
        embedding_contrast=args.embedding_contrast,
        embedding_gamma=args.embedding_gamma,
    )
