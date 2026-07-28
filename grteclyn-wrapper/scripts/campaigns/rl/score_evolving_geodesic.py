#!/usr/bin/env python
"""Post-hoc 4D evolving null-geodesic pass over finished runs.

WHY THIS EXISTS.  ``consume_plotfiles --evolving-geodesic`` does NOT run the 4D
trace.  Its only effects are to set ``GRTECLYN_EVOLVING_GEODESIC=1`` in the
consumer's own environment and to switch on the ``small_data/metric_stack``
cache at the right spatial resolution.  The trace itself lives in the metrics
aggregation layer (``metrics.aggregation.collector``), which a campaign queue
script never invokes -- so ``f_geo_evol`` / ``f_geo_evol_ok`` (cols 13/14 of
``ftl_timeseries.dat``) stay at the hardcoded ``0.0  0`` placeholders that
``consume_plotfiles/extraction/ftl.py`` writes on every row, and no
``evolving_geodesic.json`` is ever produced.

This script closes that gap without needing plotfiles: the trace reads the
cached metric slices, which live on NFS and survive plotfile pruning.  Run it
any time after a campaign finishes.

    grteclyn-wrapper/.venv/bin/python \
        grteclyn-wrapper/scripts/campaigns/rl/score_evolving_geodesic.py \
        runs/pump_ladder_m0/lad_m0_tp*

For each run directory it writes ``small_data/evolving_geodesic.json`` and
patches cols 13/14 of ``small_data/ftl_timeseries.dat``.

NOTE ON RESOLUTION.  The cache's spatial resolution is fixed at write time by
``GRTECLYN_METRIC_STACK_N_SPACE`` (default 33 in search mode, 65 in hq mode).
This script cannot raise it after the fact -- the plotfiles are gone.  It
reports the cached resolution so the number is quotable with its caveat.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import replace
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "grteclyn-wrapper" / "src"))

# Keep every cache/temp write inside our own scratch, never ~/.local or ~/.cache.
_SCRATCH = Path(os.environ.get("GRTECLYN_SCRATCH", "/tmp/grteclyn_scratch")) / "_cache"
for _var, _sub in (
    ("XDG_CACHE_HOME", ""),
    ("UV_CACHE_DIR", "uv"),
    ("MPLCONFIGDIR", "mpl"),
    ("TMPDIR", "tmp"),
    ("PYTHONPYCACHEPREFIX", "pyc"),
):
    _p = _SCRATCH / _sub if _sub else _SCRATCH
    _p.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault(_var, str(_p))

from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (  # noqa: E402
    compute_evolving_geodesic_ftl_from_metric_stack_cache,
    evolving_report_trustworthy,
    patch_ftl_timeseries_evolving,
    write_evolving_geodesic_json,
)
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic_options import (  # noqa: E402
    evolving_geodesic_options_from_env,
)
from grteclyn_wrapper.metrics.probes.ftl.metric_stack_cache import (  # noqa: E402
    cache_fidelity,
    list_slice_files,
    metric_stack_dir,
    slice_time,
)

MIN_SLICES = 3
FIDELITY_TOL = 1.5


def true_min_chi_at(run_dir: Path) -> dict[float, float]:
    """time -> min_chi as reported by the simulation (collapse_diagnostics col 3)."""
    path = run_dir / "data" / "collapse_diagnostics.dat"
    out: dict[float, float] = {}
    if not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) > 2:
            try:
                out[round(float(parts[0]), 3)] = float(parts[2])
            except ValueError:
                continue
    return out


def frozen_peak_from_ftl_timeseries(run_dir: Path) -> float | None:
    """Peak trustworthy frozen ``f_geo`` the consumer measured per plotfile.

    The consumer computed these from full-fidelity plotfiles during the
    campaign; recomputing them from the float32 cache costs two full-grid
    inversions plus gradients per slice (the dominant cost of the whole pass at
    257^3) and can only produce a number of LOWER provenance.  Reuse the column.
    """
    path = run_dir / "small_data" / "ftl_timeseries.dat"
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].lstrip().startswith("#"):
        return None
    cols = lines[0].lstrip("#").split()
    try:
        i_geo = cols.index("f_geo")
        i_ok = cols.index("geo_trustworthy")
    except ValueError:
        return None
    peak: float | None = None
    for line in lines[1:]:
        if line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) <= max(i_geo, i_ok):
            continue
        try:
            ok = float(parts[i_ok])
            val = float(parts[i_geo])
        except ValueError:
            continue
        if ok >= 0.5 and (peak is None or val > peak):
            peak = val
    return peak


def cached_n_space(cache_dir: Path) -> int | None:
    """Spatial resolution actually baked into the cache (cannot be changed now)."""
    import numpy as np

    files = list_slice_files(cache_dir)
    if not files:
        return None
    with np.load(files[0]) as slab:
        return int(slab["g"].shape[0])


def score_one(
    run_dir: Path,
    *,
    ftl_L: float,
    dry_run: bool,
    force: bool,
    max_time: float | None = None,
) -> int:
    t_start = time.monotonic()
    small = run_dir / "small_data"
    cache_dir = metric_stack_dir(small)
    all_files = list_slice_files(cache_dir)
    n_cached = len(all_files)
    if max_time is not None:
        n_used = sum(1 for p in all_files if slice_time(p) <= max_time + 1.0e-9)
        print(
            f"  {run_dir.name}: truncating stack to t<={max_time:g} "
            f"({n_used}/{n_cached} slices kept)",
            flush=True,
        )
    else:
        n_used = n_cached
    if n_used < MIN_SLICES:
        print(f"  {run_dir.name}: SKIP -- {n_used} slices (< {MIN_SLICES})")
        return 1

    n_space = cached_n_space(cache_dir)

    # A uniform resample coarser than the finest AMR level erases sharp
    # features silently.  Refuse to report a number the cache cannot support.
    bad = cache_fidelity(
        cache_dir, true_min_chi_at(run_dir), tol=FIDELITY_TOL, max_time=max_time
    )
    if bad:
        worst_t, true, cached, ratio = bad[0]
        print(
            f"  {run_dir.name}: UNFAITHFUL CACHE -- {len(bad)}/{n_used} slices "
            f"under-resolved; worst t={worst_t:.2f} true min_chi={true:.3e} "
            f"cached={cached:.3e} ({ratio:.0f}x too shallow) at {n_space}^3"
        )
        for t_bad, tr, ca, ra in bad[:5]:
            print(f"      t={t_bad:6.2f}  true={tr:.3e}  cached={ca:.3e}  {ra:6.1f}x")
        if not force:
            print(
                "      REFUSING to report f_geo_evol. Rays would not feel the "
                "well -> spurious shortcut. Re-run with a finer "
                "GRTECLYN_METRIC_STACK_N_SPACE, truncate the stack before the "
                "offending slices with --max-time, or pass --force to override."
            )
            return 1
    print(
        f"  {run_dir.name}: fidelity PASS ({n_used} slices @ {n_space}^3) "
        f"[{time.monotonic() - t_start:.0f}s]",
        flush=True,
    )

    # The consumer already measured per-plotfile frozen f_geo at plotfile
    # fidelity; skip the cache-side recompute (dominant cost at 257^3).
    frozen = frozen_peak_from_ftl_timeseries(run_dir)
    opts = replace(evolving_geodesic_options_from_env(), compute_frozen_peak=False)

    report = compute_evolving_geodesic_ftl_from_metric_stack_cache(
        cache_dir,
        options=opts,
        max_time=max_time,
        frozen_peak_override=frozen,
    )
    if report is None:
        print(f"  {run_dir.name}: FAILED -- {n_used} slices but no report")
        return 1

    extra_notes = [
        "f_geo_frozen_peak from ftl_timeseries.dat (consumer, plotfile fidelity)"
        if frozen is not None
        else "f_geo_frozen_peak unavailable (no trustworthy ftl_timeseries rows)"
    ]
    if max_time is not None:
        extra_notes.append(
            f"stack truncated to t<={max_time:g} ({n_used}/{n_cached} slices; "
            "strict no-frozen-tail guard makes truncation conservative)"
        )
    report = replace(report, notes=report.notes + tuple(extra_notes))

    ok = evolving_report_trustworthy(report)
    # Persist BEFORE printing: a broken stdout pipe (dead parent) must not cost
    # the result of an hour-long trace.
    if not dry_run:
        write_evolving_geodesic_json(small / "evolving_geodesic.json", report)
        patch_ftl_timeseries_evolving(
            small / "ftl_timeseries.dat",
            f_geo_evol=float(report.f_geo),
            f_geo_evol_ok=ok,
        )
    print(
        f"  {run_dir.name}: f_geo_evol={report.f_geo:.4e} ok={int(ok)} "
        f"rays={report.n_reached}/{report.n_rays} captured={report.n_captured} "
        f"h_ok={int(report.h_quality_ok)} t_emit={report.t_emit:.2f} "
        f"slices={n_used}@{n_space}^3 [{time.monotonic() - t_start:.0f}s]"
    )
    for note in report.notes:
        print(f"      note: {note}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dirs", type=Path, nargs="+", help="finished run directories")
    ap.add_argument("--ftl-l", type=float, default=8.0, help="match the run's ftl_L")
    ap.add_argument(
        "--dry-run", action="store_true", help="compute and print, write nothing"
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="report even when the cache cannot represent the geometry (unsafe)",
    )
    ap.add_argument(
        "--max-time",
        type=float,
        default=None,
        help=(
            "drop cached slices after this simulation time (e.g. exclude "
            "late under-resolved slices; conservative under the strict "
            "no-frozen-tail guard)"
        ),
    )
    args = ap.parse_args()

    mode = os.environ.get("GRTECLYN_EVOLVING_GEODESIC_MODE", "search")
    print(f"evolving-geodesic mode={mode}  ftl_L={args.ftl_l}")
    failures = 0
    for run_dir in args.run_dirs:
        if not (run_dir / "small_data").is_dir():
            print(f"  {run_dir}: SKIP -- no small_data/")
            failures += 1
            continue
        failures += score_one(
            run_dir,
            ftl_L=args.ftl_l,
            dry_run=args.dry_run,
            force=args.force,
            max_time=args.max_time,
        )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
