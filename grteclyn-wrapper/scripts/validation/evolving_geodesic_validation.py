#!/usr/bin/env python3
"""Compare frozen-slice vs 4D evolving null-geodesic FTL on analytic and candidate runs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from grteclyn_wrapper.metrics.probes import warpfactory as wf
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    compute_evolving_geodesic_ftl_from_analytic,
)
from grteclyn_wrapper.metrics.probes.ftl.general import find_recent_plotfiles
from grteclyn_wrapper.metrics.probes.ftl.geodesic import compute_geodesic_ftl_from_plotfile
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    compute_evolving_geodesic_ftl_from_plotfiles,
)


def _print_row(label: str, frozen: float | None, evolving: float | None) -> None:
    frozen_s = f"{frozen:.4f}" if frozen is not None else "n/a"
    evol_s = f"{evolving:.4f}" if evolving is not None else "n/a"
    if frozen is not None and evolving is not None:
        delta = evolving - frozen
        delta_s = f"{delta:+.4f}"
    else:
        delta_s = "n/a"
    print(f"{label:<28} {frozen_s:>12} {evol_s:>12} {delta_s:>12}")


def _alcubierre_case(*, velocity: float, n_space: int, half_width: float) -> tuple[float, float, bool]:
    g, spacing = wf.alcubierre_metric(
        velocity=velocity,
        n_space=n_space,
        half_width=half_width,
        dt=0.2,
    )
    evo = compute_evolving_geodesic_ftl_from_analytic(g, spacing, n_rays=5)
    frozen_peak = evo.f_geo_frozen_peak
    return frozen_peak or 0.0, evo.f_geo, evo.h_quality_ok


def _candidate_case(eval_dir: Path, *, half_width: float | None) -> tuple[float | None, float | None]:
    plotfiles = find_recent_plotfiles(eval_dir, count=5)
    if len(plotfiles) < 3:
        print(f"WARNING: {eval_dir} has <3 plotfiles; skipping candidate evolving trace.")
        return None, None

    frozen_peak = 0.0
    for path in plotfiles:
        rep = compute_geodesic_ftl_from_plotfile(path, n=65, half_width=half_width)
        if rep is not None and rep.h_quality_ok and rep.f_geo > frozen_peak:
            frozen_peak = rep.f_geo

    evo = compute_evolving_geodesic_ftl_from_plotfiles(
        [str(p) for p in plotfiles],
        n_space=65,
        half_width=half_width,
    )
    if evo is None:
        return frozen_peak or None, None
    return frozen_peak or None, evo.f_geo


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--v-s",
        type=float,
        default=2.0,
        help="Alcubierre bubble speed for the analytic positive control.",
    )
    parser.add_argument(
        "--candidate-eval",
        type=Path,
        default=None,
        help="Optional episode directory (e.g. HQ promote eval 177) for comparison.",
    )
    parser.add_argument("--n-space", type=int, default=65)
    parser.add_argument("--half-width", type=float, default=6.0)
    parser.add_argument("--ftl-l", type=float, default=None)
    args = parser.parse_args()

    print(f"{'Case':<28} {'frozen peak':>12} {'evolving':>12} {'delta':>12}")
    print("-" * 68)

    frozen_a, evo_a, h_ok = _alcubierre_case(
        velocity=args.v_s,
        n_space=args.n_space,
        half_width=args.half_width,
    )
    _print_row("Alcubierre (moving bubble)", frozen_a, evo_a)

    candidate_ok = True
    if args.candidate_eval is not None:
        eval_dir = args.candidate_eval.expanduser().resolve()
        frozen_c, evo_c = _candidate_case(eval_dir, half_width=args.ftl_l)
        _print_row(f"Candidate ({eval_dir.name})", frozen_c, evo_c)

    print()
    if not h_ok or evo_a <= 0.1:
        print(
            f"FAIL: analytic control (f_geo={evo_a:.4f}, h_quality_ok={h_ok}); "
            "expected f_geo > 0.1 with reliable integration."
        )
        return 1
    print("PASS: analytic Alcubierre evolving trace registered a real shortcut.")
    if args.candidate_eval is not None and not candidate_ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
