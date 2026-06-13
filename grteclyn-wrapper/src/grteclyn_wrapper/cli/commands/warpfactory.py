"""Warp Factory analytic metric command handler."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...core.episode import write_json
from ...metrics import dataclass_to_dict
from ...metrics.probes import warpfactory as wf


def run_warpfactory_command(args: argparse.Namespace) -> int:
    if getattr(args, "convergence", False):
        result = wf.convergence_order(
            velocity=args.velocity, half_width=args.half_width, dt=args.dt
        )
        print(json.dumps(result, indent=2))
        return 0

    if args.metric == "minkowski":
        g, spacing = wf.minkowski_metric(
            half_width=args.half_width, n_space=args.n_space, dt=args.dt
        )
    elif args.metric == "alcubierre":
        g, spacing = wf.alcubierre_metric(
            velocity=args.velocity,
            bubble_radius=args.bubble_radius,
            sigma=args.sigma,
            half_width=args.half_width,
            n_space=args.n_space,
            dt=args.dt,
        )
    else:  # pragma: no cover
        raise ValueError(f"unknown metric {args.metric!r}")

    report = wf.evaluate_four_metric(
        g,
        spacing,
        n_directions=args.n_directions,
        n_speeds=args.n_speeds,
        max_speed=args.max_speed,
    )
    payload = dataclass_to_dict(report)
    if args.output:
        write_json(Path(args.output).expanduser().resolve(), payload)
    print(json.dumps(payload, indent=2))
    return 0
