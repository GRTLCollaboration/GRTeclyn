"""Orchestrate geometry-atlas elite visualisation."""

from __future__ import annotations

import json
from pathlib import Path

from .elites import EliteRecord, load_top_elites
from .panels import plot_elite_panels, plot_gallery, render_elite


def default_output_dir(run_dir: Path) -> Path:
    """``<run_dir>/figures/geometry_atlas_elites``."""
    return Path(run_dir) / "figures" / "geometry_atlas_elites"


def _render_defaults(run_dir: Path) -> tuple[int | None, float | None]:
    """Pull (n, L) from campaign metadata when present."""
    meta_path = run_dir / "metadata.json"
    if not meta_path.exists():
        return None, None
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        render = (meta.get("config") or {}).get("render") or {}
        n = int(render["n"]) if "n" in render else None
        L = float(render["L"]) if "L" in render else None
        return n, L
    except Exception:  # noqa: BLE001
        return None, None


def visualise_top_elites(
    run_dir: Path | str,
    *,
    top_n: int = 5,
    out_dir: Path | str | None = None,
    n: int | None = None,
    L: float | None = None,
    quiver_stride: int = 4,
    gallery: bool = True,
) -> list[Path]:
    """Load top elites, write per-elite panels (+ optional gallery).

    Returns the list of written PNG paths.
    """
    run_dir = Path(run_dir)
    out_dir = Path(out_dir) if out_dir is not None else default_output_dir(run_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("rank*_cell*_fgeo*.png", "gallery_top*.png"):
        for previous_plot in out_dir.glob(pattern):
            previous_plot.unlink()

    meta_n, meta_L = _render_defaults(run_dir)
    if n is None:
        n = meta_n
    if L is None:
        L = meta_L

    elites = load_top_elites(run_dir, top_n=top_n)
    if not elites:
        raise FileNotFoundError(f"No elites with genomes found under {run_dir}")

    written: list[Path] = []
    for elite in elites:
        rendered, cfg = render_elite(elite, n=n, L=L)
        path = out_dir / (
            f"rank{elite.rank:02d}_cell{elite.cell[0]}_{elite.cell[1]}"
            f"_fgeo{elite.f_geo:.3f}.png"
        )
        written.append(
            plot_elite_panels(
                elite, rendered, cfg, out_path=path, quiver_stride=quiver_stride
            )
        )

    if gallery:
        gpath = out_dir / f"gallery_top{len(elites)}.png"
        written.append(plot_gallery(elites, out_path=gpath, n=n, L=L))

    _write_manifest(out_dir / "manifest.txt", elites, written)
    return written


def _write_manifest(
    path: Path, elites: list[EliteRecord], written: list[Path]
) -> None:
    lines = ["# geometry-atlas elite visualisation", ""]
    for e in elites:
        modes = ", ".join(e.active_modes())
        lines.append(
            f"#{e.rank} f_geo={e.f_geo:.4f} cell={list(e.cell)} "
            f"family={e.family_label} shift={e.shift_fraction:.3f} "
            f"E-={e.integral_negative_rho:.3f} modes=[{modes}]"
        )
    lines.append("")
    lines.append("# written files")
    for p in written:
        lines.append(str(p))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


__all__ = ["default_output_dir", "visualise_top_elites"]
