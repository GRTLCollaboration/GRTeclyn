#!/usr/bin/env python3
"""Copy one rendered frame every dt from a run's frames/ tree.

A finished cell holds roughly 250 frames per field across ~19 fields.  All of
them would be hundreds of megabytes of PNG in a repository that only needs
enough stills to show a reader what the movie shows, so this keeps one frame
per dt bucket and leaves the rest -- and the movies stitched from them -- in
the gitignored run tree.

Frame index -> time comes from the run's slice cache, which records the time of
every frame it wrote and shares its file numbering with the PNGs.  A run
launched without ``--frames-cache-slices`` has no such record, and is skipped
rather than guessed at.

    thin_frames.py <run>/frames <dest> [--dt 10]
"""

from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import sys

INDEX_RE = re.compile(r"(\d+)\.(?:npz|png)$")


def frame_times(cache_root: pathlib.Path) -> dict[int, float]:
    """index -> simulation time, read from any one cached series.

    Every series a run caches is written on the same step cadence, so the first
    non-empty one answers for all of them.
    """
    if not cache_root.is_dir():
        return {}
    import numpy as np

    for series in sorted(cache_root.iterdir()):
        if not series.is_dir():
            continue
        times: dict[int, float] = {}
        for f in sorted(series.glob("*.npz")):
            m = INDEX_RE.search(f.name)
            if m:
                times[int(m.group(1))] = float(np.load(f, allow_pickle=True)["time"])
        if times:
            return times
    return {}


def pick_indices(times: dict[int, float], dt: float) -> dict[int, float]:
    """The frame nearest each multiple of dt, from t = 0 up to the last frame."""
    picked: dict[int, float] = {}
    t_max = max(times.values())
    k = 0
    while k * dt <= t_max + 1e-9:
        target = k * dt
        idx = min(times, key=lambda i: abs(times[i] - target))
        picked[idx] = times[idx]
        k += 1
    return picked


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("frames_dir")
    ap.add_argument("dest")
    ap.add_argument("--dt", type=float, default=10.0)
    args = ap.parse_args(argv)

    src, dst = pathlib.Path(args.frames_dir), pathlib.Path(args.dest)
    times = frame_times(src / "_slice_cache")
    if not times:
        print(f"[thin-frames] {src}: no slice cache -- skipping", file=sys.stderr)
        return 0

    wanted = pick_indices(times, args.dt)
    kept, fields = 0, 0
    for series in sorted(src.iterdir()):
        if not series.is_dir() or series.name == "_slice_cache":
            continue
        # Rendered PNGs live one level down for cached series and at the top
        # for the projection series, which are drawn straight from plotfiles.
        pngs = sorted((series / "frames").glob("*.png")) or sorted(series.glob("*.png"))
        if not pngs:
            continue
        by_index = {}
        for f in pngs:
            m = INDEX_RE.search(f.name)
            if m:
                by_index[int(m.group(1))] = f
        if not by_index:
            continue
        out_dir = dst / series.name
        out_dir.mkdir(parents=True, exist_ok=True)
        fields += 1
        for idx in sorted(wanted):
            # A projection series may be rendered on a coarser cadence than the
            # cache, so fall back to its nearest available frame.
            pick = by_index.get(idx)
            if pick is None:
                pick = by_index[min(by_index, key=lambda i: abs(i - idx))]
            shutil.copy2(pick, out_dir / pick.name)
            kept += 1

    if fields:
        kept_times = ", ".join(f"{t:g}" for t in sorted(wanted.values()))
        (dst / "FRAMES.md").write_text(
            "# Frames\n\n"
            f"One frame every dt = {args.dt:g} from each field's series, chosen by the\n"
            "time recorded in the run's slice cache.  The full series (~250 frames per\n"
            "field) and the movies stitched from them stay in the gitignored run tree;\n"
            "these stills exist so the repository shows what the movie shows.\n\n"
            f"times kept: {kept_times}\n"
        )
    print(f"[thin-frames] {kept} stills across {fields} fields (dt = {args.dt:g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
