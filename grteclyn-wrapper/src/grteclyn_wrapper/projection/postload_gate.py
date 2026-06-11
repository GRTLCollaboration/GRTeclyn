"""GRTeclyn post-load constraint gate for projected initial data."""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from ..core.config import DEFAULT_RADIAL_RECIPE_TEMPLATE, resolve_example, resolve_executable
from ..core.episode import create_episode, update_metadata
from ..core.params import write_params
from ..core.runner import run_episode
from ..metrics.diagnostics.constraints import read_constraint_metrics

DEFAULT_MAX_HAM_L2 = 1.0e-2
DEFAULT_MAX_MOM_L2 = 1.0e-2
_HEAVY_ARTIFACT_PREFIXES = ("RadialRecipePlt", "RadialRecipeChk")


def _cleanup_heavy_artifacts(episode_dir: Path) -> None:
    """Drop AMReX plot/checkpoint trees once constraint metrics are extracted."""
    if not episode_dir.is_dir():
        return
    for child in episode_dir.iterdir():
        if not child.is_dir():
            continue
        if any(child.name.startswith(prefix) for prefix in _HEAVY_ARTIFACT_PREFIXES):
            shutil.rmtree(child, ignore_errors=True)


@dataclass(frozen=True)
class PostLoadGateConfig:
    max_hamiltonian_l2: float = DEFAULT_MAX_HAM_L2
    max_momentum_l2: float = DEFAULT_MAX_MOM_L2
    stop_time: float = 0.01
    check_params: bool = False


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
    return PostLoadGateResult(
        passed=passed,
        max_hamiltonian_l2=metrics.max_hamiltonian_l2,
        max_momentum_l2=metrics.max_momentum_l2,
        final_hamiltonian_l2=metrics.final_hamiltonian_l2,
        final_momentum_l2=metrics.final_momentum_l2,
        episode_path=str(Path(constraint_metrics_path).parent.parent),
        reason=reason,
        notes=(),
    )


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
    gate = evaluate_constraint_gate(constraint_path, config=cfg)
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
    )


def write_gate_result(result: PostLoadGateResult, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(asdict(result), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
