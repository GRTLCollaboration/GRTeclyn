"""One fixed colour scale per field, measured over a whole run.

WHY THIS EXISTS.  Frames are rendered as plotfiles arrive, so the renderer
cannot know how large a field will grow later.  The two answers it had were both
unusable for a paper movie:

* rescale every frame from its own percentiles (``per_frame_zlim``) -- the
  colourbar and its tick labels change in every single frame, so the movie
  visibly jitters and the same colour means a different number in each frame;
* lock the scale from the first plotfile -- static, but ``chi`` starts nearly
  flat and develops its well later, so everything interesting is off the end of
  the scale by the middle of the run.

Scanning fixes both.  Given the plotfiles of a finished run, this measures the
range each frame *would* have chosen and takes the envelope: the tightest single
scale that clips nothing the per-frame scale would have shown.  Every frame is
then rendered against that one scale, so the bar never moves and colours are
comparable across time.

It needs the plotfiles, so it is for rendering after a run, not during one.
Streaming runs delete plotfiles as they are consumed; keep them (``--keep-last``)
if the movies are going in a paper.

Cost is one slice extraction per plotfile per field -- roughly a rendering pass
without the drawing.  ``stride`` subsamples when that is too slow; the envelope
is then taken over the sampled subset, which is safe as long as the sampling is
dense enough to catch the extremes, so it is reported.
"""

from __future__ import annotations

import json
import os
from typing import Sequence

from .zlim import _SCANNED_KEY, _lock_frame_zlims_from_plotfile

#: Written next to the frames so a movie can be traced back to the scale it used.
SCAN_RECORD_NAME = "zlim_scan.json"


def scan_series_zlims(
    plotfiles: Sequence[str],
    args_dict: dict,
    *,
    stride: int = 1,
    record_dir: str | None = None,
    verbose: bool = False,
) -> dict[str, list[float]]:
    """Envelope of the per-plotfile colour limits over ``plotfiles``.

    Returns a ``frame_zlims`` dict carrying the scanned-series sentinel, ready
    to hand to the renderer.  Returns an empty dict if nothing could be read, so
    the caller falls back to its previous behaviour rather than rendering
    against a bogus scale.
    """
    ordered = sorted(plotfiles)
    if stride > 1:
        # Always keep the last plotfile: the extremes usually live at late times,
        # and dropping it is the one sampling error that would clip the movie.
        sampled = ordered[::stride]
        if ordered and ordered[-1] not in sampled:
            sampled.append(ordered[-1])
    else:
        sampled = ordered
    if not sampled:
        return {}

    envelope: dict[str, list[float]] = {}
    scanned = 0
    for path in sampled:
        try:
            one = _lock_frame_zlims_from_plotfile(
                path, args_dict, include_per_frame_fields=True
            )
        except Exception as exc:
            if verbose:
                print(f"WARNING: zlim-scan skipped {os.path.basename(path)}: {exc}")
            continue
        scanned += 1
        for field, (lo, hi) in one.items():
            cur = envelope.get(field)
            if cur is None:
                envelope[field] = [float(lo), float(hi)]
            else:
                cur[0] = min(cur[0], float(lo))
                cur[1] = max(cur[1], float(hi))

    if not envelope:
        return {}

    print(
        f"[zlim-scan] {len(envelope)} field(s) over {scanned}/{len(ordered)} plotfile(s)"
        + (f" (stride {stride})" if stride > 1 else "")
    )
    if verbose:
        for field in sorted(envelope):
            lo, hi = envelope[field]
            print(f"[zlim-scan]   {field}: {lo:.6g} .. {hi:.6g}")

    if record_dir:
        try:
            os.makedirs(record_dir, exist_ok=True)
            with open(os.path.join(record_dir, SCAN_RECORD_NAME), "w", encoding="utf-8") as fh:
                json.dump(
                    {
                        "zlims": envelope,
                        "plotfiles_scanned": scanned,
                        "plotfiles_total": len(ordered),
                        "stride": stride,
                        "first": os.path.basename(sampled[0]),
                        "last": os.path.basename(sampled[-1]),
                    },
                    fh,
                    indent=2,
                    sort_keys=True,
                )
        except OSError as exc:
            print(f"WARNING: could not write {SCAN_RECORD_NAME}: {exc}")

    envelope[_SCANNED_KEY] = True  # type: ignore[assignment]
    return envelope
