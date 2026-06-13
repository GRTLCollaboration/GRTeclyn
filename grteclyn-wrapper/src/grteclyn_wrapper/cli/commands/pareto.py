"""Pareto front extraction command handler."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...core.episode import write_json
from ...search.pareto import front_to_dict, load_trajectory_points, pareto_front


def run_pareto_command(args: argparse.Namespace) -> int:
    points = load_trajectory_points(Path(args.trajectory).expanduser().resolve())
    front = pareto_front(points)
    payload = front_to_dict(front)
    if args.output:
        write_json(Path(args.output).expanduser().resolve(), payload)
    print(json.dumps(payload, indent=2))
    return 0
