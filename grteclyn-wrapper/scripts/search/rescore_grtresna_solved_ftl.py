#!/usr/bin/env python3
"""Offline re-score GRTresna search evals with solved-geometry operational FTL."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from grteclyn_wrapper.metrics.ftl_solved_geometry import (
    compute_solved_geometry_ftl,
)
from grteclyn_wrapper.search.solved_ftl_gate import (
    SolvedFtlGateConfig,
    solved_ftl_has_signal,
    solved_geometry_ftl_to_dict,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "campaign_dir",
        type=Path,
        help="e.g. runs/grtresna_search/optimize_20260602T181607Z",
    )
    parser.add_argument("--L", type=float, default=8.0, help="integration half-width")
    parser.add_argument("--evals", nargs="*", help="specific eval dirs (default: all)")
    parser.add_argument("--solved-ftl-f-op-floor", type=float, default=1.0e-4)
    parser.add_argument("--solved-ftl-near-luminal-speed-floor", type=float, default=0.99)
    parser.add_argument("--solved-ftl-superluminal-speed-floor", type=float, default=1.01)
    parser.add_argument("--solved-ftl-superluminal-fraction-floor", type=float, default=0.02)
    parser.add_argument("--solved-ftl-max-physical-coord-speed", type=float, default=8.0)
    parser.add_argument("--solved-ftl-max-physical-f-op", type=float, default=0.85)
    parser.add_argument("--solved-ftl-rejection-speed-target", type=float, default=1.01)
    args = parser.parse_args()
    gate_config = SolvedFtlGateConfig(
        f_op_floor=args.solved_ftl_f_op_floor,
        near_luminal_speed_floor=args.solved_ftl_near_luminal_speed_floor,
        superluminal_speed_floor=args.solved_ftl_superluminal_speed_floor,
        superluminal_fraction_floor=args.solved_ftl_superluminal_fraction_floor,
        max_physical_coord_speed=args.solved_ftl_max_physical_coord_speed,
        max_physical_f_op=args.solved_ftl_max_physical_f_op,
        rejection_speed_target=args.solved_ftl_rejection_speed_target,
    )

    campaign = args.campaign_dir.expanduser().resolve()
    if args.evals:
        eval_dirs = [campaign / e for e in args.evals]
    else:
        eval_dirs = sorted(campaign.glob("eval_*"))

    rows: list[dict] = []
    total = len(eval_dirs)
    for n, ep in enumerate(eval_dirs, start=1):
        gi = ep / "initial_data.gridinit"
        if not gi.is_file():
            rows.append({"episode": ep.name, "status": "no_gridinit"})
            print(f"[{n}/{total}] {ep.name}: no_gridinit", flush=True)
            continue
        solved = compute_solved_geometry_ftl(gi, L=args.L)
        if solved is None:
            rows.append({"episode": ep.name, "status": "compute_failed"})
            print(f"[{n}/{total}] {ep.name}: compute_failed", flush=True)
            continue
        summary = solved_geometry_ftl_to_dict(solved, config=gate_config)
        summary["episode"] = ep.name
        summary["would_gate"] = not solved_ftl_has_signal(solved, config=gate_config)
        rows.append(summary)
        print(
            f"[{n}/{total}] {ep.name}: F_op={summary['f_op']:.4f} "
            f"max_c={summary['max_local_speed']:.3f} "
            f"mech={summary['mechanism_descriptor']:.2f} "
            f"gate={summary['would_gate']}",
            flush=True,
        )

    out_path = campaign / "solved_ftl_rescore.json"
    out_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")

    scored = [r for r in rows if "f_op" in r]
    gated = sum(1 for r in scored if r["would_gate"])
    print(f"\nWrote {out_path} ({len(rows)} evals)")
    print(
        f"Summary: {len(scored)} scored, {gated} would-gate (no solved FTL), "
        f"{len(scored) - gated} pass the gate"
    )


if __name__ == "__main__":
    main()
