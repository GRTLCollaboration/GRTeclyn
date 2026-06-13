from __future__ import annotations

import os


# Time-resolved FTL feature stream.  One row per plotfile, written alongside the
# other per-plotfile diagnostics in small_data/.  The FTL signal peaks mid-run and
# diffuses, so scoring only the final frame is half-blind; this stream lets the
# scorer/aggregator see the whole f_geo(t) / f_op(t) curve and average over it.
FTL_TIMESERIES_HEADER = (
    "# time  f_op  f_geo  geo_trustworthy  max_local_speed  "
    "superluminal_fraction  max_shift  structure_coherence  reachable  "
    "n_rays  n_reached  max_h_rel_drift"
)


def _extract_ftl_timeseries_line(
    p: str,
    *,
    t: float,
    center,
    ftl_L: float | None,
    verbose: bool = False,
) -> str | None:
    """Compute per-plotfile FTL features (general + gated geodesic) -> dat row.

    The cheap operational-FTL probe runs on every plotfile; the expensive
    gauge-invariant geodesic probe only fires when a coordinate channel is
    present (same gate the final-frame collector uses), keeping per-frame cost
    bounded.  Returns None if the general probe fails entirely (no row written).
    """
    from grteclyn_wrapper.metrics.probes.ftl.general import (
        compute_general_ftl_from_plotfile,
    )

    try:
        gen = compute_general_ftl_from_plotfile(p, n=97, L=ftl_L)
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(f"WARNING: ftl-timeseries general probe failed for {os.path.basename(p)}: {exc}")
        return None
    if gen is None:
        return None

    f_op = float(gen.f_op)
    max_speed = float(gen.max_local_speed)
    superlum = float(gen.superluminal_fraction)
    max_shift = float(gen.max_shift)
    coherence = gen.structure_coherence
    reachable = 1 if gen.reachable else 0

    # Gated gauge-invariant geodesic probe: only worth its cost when the cheap
    # probe already shows a coordinate channel (f_op above floor or a
    # superluminal cell).  Mirrors collector.read_episode_metrics gating.
    f_geo = 0.0
    geo_trustworthy = 0
    n_rays = 0
    n_reached = 0
    max_h_rel = 0.0
    if f_op > 1.0e-3 or max_speed > 1.0:
        from grteclyn_wrapper.metrics.probes.ftl.geodesic import (
            compute_geodesic_ftl_from_plotfile,
        )

        try:
            geo = compute_geodesic_ftl_from_plotfile(p, n=65, half_width=ftl_L)
        except Exception as exc:  # noqa: BLE001
            geo = None
            if verbose:
                print(f"WARNING: ftl-timeseries geodesic probe failed for {os.path.basename(p)}: {exc}")
        if geo is not None:
            f_geo = float(geo.f_geo)
            n_rays = int(geo.n_rays)
            n_reached = int(geo.n_reached)
            max_h_rel = float(geo.max_h_rel_drift)
            geo_trustworthy = 1 if (
                geo.h_quality_ok and geo.n_rays > 0 and geo.n_reached == geo.n_rays
            ) else 0

    coherence_str = f"{float(coherence):.16e}" if coherence is not None else "nan"
    return (
        f"{t:.16e}  {f_op:.16e}  {f_geo:.16e}  {geo_trustworthy:d}  "
        f"{max_speed:.16e}  {superlum:.16e}  {max_shift:.16e}  {coherence_str}  "
        f"{reachable:d}  {n_rays:d}  {n_reached:d}  {max_h_rel:.16e}"
    )
