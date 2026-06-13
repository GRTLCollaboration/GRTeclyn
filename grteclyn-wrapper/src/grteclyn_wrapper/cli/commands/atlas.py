"""Atlas command handler."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...core.config import resolve_example, resolve_executable
from ...search.atlas import run_atlas


def run_atlas_command(args: argparse.Namespace, base_overrides: dict) -> int:
    example = resolve_example(args.example)
    executable = None
    if not args.dry_run:
        executable = resolve_executable(
            args.executable,
            example=example,
            mpi_ranks=args.mpi_ranks,
            comp=args.comp,
            cuda=not args.no_cuda,
            debug=args.debug,
        )

    template = Path(args.template).expanduser().resolve() if args.template else example.template
    paths, records, summary = run_atlas(
        runs_dir=Path(args.runs_dir).expanduser().resolve(),
        executable=executable,
        count=args.count,
        seed=args.seed,
        base_overrides=base_overrides,
        template=template,
        example=example,
        name=args.name,
        dry_run=args.dry_run,
        stop_on_failure=args.stop_on_failure,
        cuda_devices=args.cuda_devices,
        check_params=not args.skip_check_params,
        constrained=getattr(args, "constrained", False),
        phantom=getattr(args, "phantom", False),
        preflight=getattr(args, "preflight", False),
    )
    print(json.dumps({"atlas": str(paths.root), "records": len(records), "summary": summary}, indent=2))
    return 0
