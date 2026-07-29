"""GRTeclyn post-load constraint gate for projected initial data."""

from __future__ import annotations

import json
import logging
import os
import shutil
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

_log = logging.getLogger(__name__)

from ..core.config import DEFAULT_RADIAL_RECIPE_TEMPLATE, resolve_example, resolve_executable
from ..core.episode import create_episode, update_metadata
from ..core.params import write_params
from ..core.runner import run_episode
from ..core.scratch import purge_plotfile_scratch
from ..metrics.diagnostics.constraints import read_constraint_metrics
from ..metrics.types.diagnostics import ConstraintMetrics

DEFAULT_MAX_HAM_L2 = 1.0e-2
DEFAULT_MAX_MOM_L2 = 1.0e-2
_HEAVY_ARTIFACT_PREFIXES = ("RadialRecipePlt", "RadialRecipeChk")

# NFS read-after-write retry parameters
_NFS_MAX_RETRIES = 5
_NFS_RETRY_DELAY_S = 0.5


def _invalidate_nfs_cache(path: Path) -> None:
    """Force NFS to re-validate the directory entry cache for *path*.

    On NFS, a negative-lookup or stale attribute cache can cause
    ``path.exists()`` to return ``False`` (or ``open()`` to see a
    truncated file) even though the server already has the data.
    Listing the parent directory forces the NFS client to issue a
    READDIR/LOOKUP, refreshing the dentry cache.
    """
    parent = path.parent
    try:
        os.listdir(parent)
    except OSError:
        pass


def _cleanup_heavy_artifacts(episode_dir: Path) -> None:
    """Drop AMReX plot/checkpoint trees once constraint metrics are extracted."""
    if not episode_dir.is_dir():
        return
    for child in episode_dir.iterdir():
        if not child.is_dir():
            continue
        if any(child.name.startswith(prefix) for prefix in _HEAVY_ARTIFACT_PREFIXES):
            shutil.rmtree(child, ignore_errors=True)
    # The gate runs no consumer, so nothing it writes is ever in the ledger and
    # the constraint numbers it wanted are already out of the .dat file.  Force.
    purge_plotfile_scratch(episode_dir, force=True)


@dataclass(frozen=True)
class PostLoadGateConfig:
    max_hamiltonian_l2: float = DEFAULT_MAX_HAM_L2
    max_momentum_l2: float = DEFAULT_MAX_MOM_L2
    stop_time: float = 0.01
    check_params: bool = False
    # Whole-hierarchy thresholds (cols 12 and 16 of constraint_norms.dat).
    # Deliberately unset: no run has yet produced those columns, so any number
    # here would be invented rather than measured.  While they are None the gate
    # still *reports* the composite norms in `notes`, which is what a later
    # calibration needs.  Set them once real gate runs have accumulated.
    max_hamiltonian_l2_amr: float | None = None
    max_hamiltonian_linf_amr: float | None = None


@dataclass(frozen=True)
class PostLoadGateResult:
    passed: bool
    max_hamiltonian_l2: float | None
    max_momentum_l2: float | None
    final_hamiltonian_l2: float | None
    final_momentum_l2: float | None
    episode_path: str | None
    reason: str | None
    notes: tuple[str, ...]
    max_hamiltonian_l2_amr: float | None = None
    max_hamiltonian_linf_amr: float | None = None
    max_hamiltonian_l2_amr_ref: float | None = None
    has_refined_region: bool = False


def evaluate_constraint_gate(
    constraint_metrics_path: str | Path,
    *,
    config: PostLoadGateConfig | None = None,
) -> PostLoadGateResult:
    """Check an existing constraint_norms.dat without launching GRTeclyn."""
    cfg = config or PostLoadGateConfig()
    metrics = read_constraint_metrics(Path(constraint_metrics_path))
    if metrics is None:
        return PostLoadGateResult(
            passed=False,
            max_hamiltonian_l2=None,
            max_momentum_l2=None,
            final_hamiltonian_l2=None,
            final_momentum_l2=None,
            episode_path=None,
            reason="constraint_norms.dat missing or empty",
            notes=(),
        )

    passed = (
        metrics.max_hamiltonian_l2 is not None
        and metrics.max_momentum_l2 is not None
        and metrics.max_hamiltonian_l2 <= cfg.max_hamiltonian_l2
        and metrics.max_momentum_l2 <= cfg.max_momentum_l2
    )
    reason = None if passed else "post-load constraints exceed threshold"

    composite_passed, notes = _check_composite_norms(metrics, cfg)
    if not composite_passed:
        passed = False
        reason = "post-load composite (whole-hierarchy) constraints exceed threshold"

    return PostLoadGateResult(
        passed=passed,
        max_hamiltonian_l2=metrics.max_hamiltonian_l2,
        max_momentum_l2=metrics.max_momentum_l2,
        final_hamiltonian_l2=metrics.final_hamiltonian_l2,
        final_momentum_l2=metrics.final_momentum_l2,
        episode_path=str(Path(constraint_metrics_path).parent.parent),
        reason=reason,
        notes=notes,
        max_hamiltonian_l2_amr=metrics.max_hamiltonian_l2_amr,
        max_hamiltonian_linf_amr=metrics.max_hamiltonian_linf_amr,
        max_hamiltonian_l2_amr_ref=metrics.max_hamiltonian_l2_amr_ref,
        has_refined_region=metrics.has_refined_region,
    )


def _check_composite_norms(
    metrics: ConstraintMetrics,
    cfg: PostLoadGateConfig,
) -> tuple[bool, tuple[str, ...]]:
    """Gate on the whole-hierarchy norms, and say what could not be judged.

    The pass/fail above uses cols 2-3, which are level-0 only, include the cells
    hidden under the refinement, and are averaged over a domain that is almost
    entirely vacuum.  Cols 12/16/17 are the composite: every level, covered
    cells dropped, plus an undiluted peak and a refined-region-only norm.
    Whenever this run produced them we report them, and enforce whichever
    thresholds the caller has actually calibrated.
    """
    notes: list[str] = []

    if metrics.max_hamiltonian_l2_amr is None:
        notes.append(
            "constraint_norms.dat predates the whole-hierarchy columns; "
            "pass/fail rests on the level-0 domain mean alone"
        )
        return True, tuple(notes)

    notes.append(
        "composite L2(Ham)={:.6e} Linf(Ham)={:.6e}".format(
            metrics.max_hamiltonian_l2_amr,
            metrics.max_hamiltonian_linf_amr or 0.0,
        )
    )

    if metrics.has_refined_region:
        notes.append(
            "refined-region L2(Ham)={:.6e}".format(
                metrics.max_hamiltonian_l2_amr_ref
            )
        )
    else:
        # The gate run is one coarse step long and the chi tagger does not fire
        # on freshly projected data, so in practice no level 1 is ever built.
        # Col 17 is then written as exactly 0.0, which reads as flawless and
        # means unmeasured.  Never let that count as evidence of anything.
        notes.append(
            "no refinement was built: the refined-region norm is unmeasured, "
            "not zero, and the composite columns describe a single-level grid"
        )

    ok = True
    if cfg.max_hamiltonian_l2_amr is not None:
        ok = ok and metrics.max_hamiltonian_l2_amr <= cfg.max_hamiltonian_l2_amr
    if (
        cfg.max_hamiltonian_linf_amr is not None
        and metrics.max_hamiltonian_linf_amr is not None
    ):
        ok = ok and metrics.max_hamiltonian_linf_amr <= cfg.max_hamiltonian_linf_amr
    if cfg.max_hamiltonian_l2_amr is None and cfg.max_hamiltonian_linf_amr is None:
        notes.append(
            "no composite threshold configured; reporting only "
            "(see PostLoadGateConfig)"
        )
    return ok, tuple(notes)


def run_postload_gate(
    gridinit_path: str | Path,
    *,
    out_dir: str | Path,
    config: PostLoadGateConfig | None = None,
    cuda_devices: str | None = None,
    dry_run: bool = False,
    overrides: Mapping[str, Any] | None = None,
) -> PostLoadGateResult:
    """Load a .gridinit in RadialRecipe and verify GRTeclyn constraints."""
    cfg = config or PostLoadGateConfig()
    gridinit_path = Path(gridinit_path).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    example = resolve_example("RadialRecipe")
    template = DEFAULT_RADIAL_RECIPE_TEMPLATE
    episode = create_episode(
        out_dir,
        name="postload_gate",
        metadata={"mode": "postload_gate", "gridinit": str(gridinit_path)},
    )

    # Promotion/search overrides carry stop_time=50 etc.; gate-specific keys must
    # win so this stays a short constraint check, not a full evolution.
    gate_overrides: dict[str, Any] = dict(overrides) if overrides else {}
    gate_overrides.update(
        {
            "recipe_initial_data_file": str(gridinit_path),
            "stop_time": cfg.stop_time,
            "plot_interval": 1000,
            "checkpoint_interval": -1,
        }
    )
    write_params(
        template,
        episode.params_path,
        episode_dir=episode.path,
        example=example,
        overrides=gate_overrides,
    )

    if dry_run:
        update_metadata(episode, {"dry_run": True})
        return PostLoadGateResult(
            passed=True,
            max_hamiltonian_l2=None,
            max_momentum_l2=None,
            final_hamiltonian_l2=None,
            final_momentum_l2=None,
            episode_path=str(episode.path),
            reason=None,
            notes=("dry_run: GRTeclyn not launched",),
        )

    executable = resolve_executable(example=example)
    result = run_episode(
        episode,
        executable,
        check_params=cfg.check_params,
        cuda_devices=cuda_devices,
        consume_plotfiles=False,
    )
    update_metadata(episode, {"simulation_exit_code": result.returncode})

    constraint_path = episode.path / "data" / "constraint_norms.dat"

    # NFS read-after-write mitigation: the simulation subprocess has
    # exited, but NFS attribute / dentry caches may still report the
    # file as missing or return a stale (truncated) view.  Retry with
    # cache invalidation before giving up.
    gate: PostLoadGateResult | None = None
    for attempt in range(_NFS_MAX_RETRIES + 1):
        _invalidate_nfs_cache(constraint_path)
        gate = evaluate_constraint_gate(constraint_path, config=cfg)
        if gate.reason != "constraint_norms.dat missing or empty":
            break
        if attempt < _NFS_MAX_RETRIES:
            _log.debug(
                "constraint_norms.dat not yet visible (attempt %d/%d), "
                "retrying in %.1fs …",
                attempt + 1,
                _NFS_MAX_RETRIES,
                _NFS_RETRY_DELAY_S,
            )
            time.sleep(_NFS_RETRY_DELAY_S)
    assert gate is not None  # loop always executes at least once

    _cleanup_heavy_artifacts(episode.path)
    return PostLoadGateResult(
        passed=gate.passed and result.returncode == 0,
        max_hamiltonian_l2=gate.max_hamiltonian_l2,
        max_momentum_l2=gate.max_momentum_l2,
        final_hamiltonian_l2=gate.final_hamiltonian_l2,
        final_momentum_l2=gate.final_momentum_l2,
        episode_path=str(episode.path),
        reason=gate.reason if gate.passed else gate.reason or f"exit_code={result.returncode}",
        notes=gate.notes,
        max_hamiltonian_l2_amr=gate.max_hamiltonian_l2_amr,
        max_hamiltonian_linf_amr=gate.max_hamiltonian_linf_amr,
        max_hamiltonian_l2_amr_ref=gate.max_hamiltonian_l2_amr_ref,
        has_refined_region=gate.has_refined_region,
    )


def write_gate_result(result: PostLoadGateResult, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(asdict(result), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
