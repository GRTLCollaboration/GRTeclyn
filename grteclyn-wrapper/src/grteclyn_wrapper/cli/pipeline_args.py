"""CLI arguments for the pipelined evaluation scheduler."""

from __future__ import annotations

import argparse
import os


def add_pipeline_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--gpu-slots-per-device",
        type=int,
        default=1,
        help="Concurrent evolution slots per physical GPU (default 1 until benchmark).",
    )
    parser.add_argument(
        "--gpu-memory-gb-per-job",
        type=float,
        default=None,
        help="Optional VRAM ceiling per evolution job for sizing hints.",
    )
    parser.add_argument(
        "--cluster-cpu-fraction",
        type=float,
        default=None,
        help="Max fraction of node CPUs for this pipeline (default from CLUSTER_CPU_FRACTION or 0.30).",
    )
    parser.add_argument(
        "--pipeline-cpu-share",
        type=float,
        default=None,
        help="Fraction of cluster CPU budget for this process (MODE=par uses 0.333).",
    )
    parser.add_argument(
        "--max-concurrent-grtresna",
        type=int,
        default=0,
        help="Override GRTresna concurrency (0 = auto from CPU budget).",
    )
    parser.add_argument(
        "--grtresna-reserve-cores",
        type=int,
        default=None,
        help="Cores reserved for Python driver overhead within CPU budget.",
    )
    parser.add_argument(
        "--remove-partial",
        action="store_true",
        help="On resume, delete partial eval dirs without score.json (default: keep).",
    )
    parser.add_argument(
        "--no-pipeline",
        action="store_true",
        help="Disable pipelined evaluation (legacy batch+join mode).",
    )


def pipeline_settings_from_args(args: argparse.Namespace) -> dict:
    pipeline_share = getattr(args, "pipeline_cpu_share", None)
    if pipeline_share is None:
        raw = os.environ.get("PIPELINE_CPU_SHARE", "").strip()
        pipeline_share = float(raw) if raw else 1.0

    cluster_fraction = getattr(args, "cluster_cpu_fraction", None)
    if cluster_fraction is None:
        raw = os.environ.get("CLUSTER_CPU_FRACTION", "").strip()
        cluster_fraction = float(raw) if raw else 0.30

    reserve_cores = getattr(args, "grtresna_reserve_cores", None)
    if reserve_cores is None:
        reserve_cores = 2 if pipeline_share < 1.0 else 4

    return {
        "slots_per_gpu": getattr(args, "gpu_slots_per_device", 1),
        "memory_gb_per_job": getattr(args, "gpu_memory_gb_per_job", None),
        "cluster_cpu_fraction": cluster_fraction,
        "pipeline_share": pipeline_share,
        "max_concurrent_grtresna": getattr(args, "max_concurrent_grtresna", 0),
        "reserve_cores": reserve_cores,
        "remove_partial": getattr(args, "remove_partial", False),
        "use_pipeline": not getattr(args, "no_pipeline", False),
    }
