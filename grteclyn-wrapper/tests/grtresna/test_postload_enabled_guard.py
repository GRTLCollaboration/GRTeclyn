"""Tests that postload gate respects config.postload_enabled."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.projection.postload_gate import PostLoadGateConfig
from grteclyn_wrapper.search.grtresna_evaluation_gates import (
    GRTresnaPreEvolutionGateConfig,
    apply_grtresna_pre_evolution_gates,
)


@pytest.fixture
def gate_context(tmp_path: Path):
    episode = create_episode(tmp_path, name="gate_episode")
    gridinit = tmp_path / "initial_data.gridinit"
    gridinit.write_bytes(b"placeholder")
    return episode, gridinit


def test_postload_gate_skipped_when_disabled(gate_context) -> None:
    episode, gridinit = gate_context
    with patch(
        "grteclyn_wrapper.search.grtresna_evaluation_gates.reject_postload_gate",
    ) as mock_postload:
        rejection = apply_grtresna_pre_evolution_gates(
            episode=episode,
            convergence={"ham_pct": 1.0, "mom_pct": 1.0},
            gridinit_path=gridinit,
            gte_overrides={"recipe_initial_data_file": str(gridinit)},
            cuda_devices="0",
            config=GRTresnaPreEvolutionGateConfig(
                postload_enabled=False,
                postload_config=PostLoadGateConfig(),
            ),
        )
    assert rejection is None
    mock_postload.assert_not_called()


def test_postload_gate_runs_when_enabled(gate_context) -> None:
    episode, gridinit = gate_context
    with patch(
        "grteclyn_wrapper.search.grtresna_evaluation_gates.reject_postload_gate",
        return_value=None,
    ) as mock_postload:
        rejection = apply_grtresna_pre_evolution_gates(
            episode=episode,
            convergence={"ham_pct": 1.0, "mom_pct": 1.0},
            gridinit_path=gridinit,
            gte_overrides={"recipe_initial_data_file": str(gridinit)},
            cuda_devices="0",
            config=GRTresnaPreEvolutionGateConfig(
                postload_enabled=True,
                postload_config=PostLoadGateConfig(),
            ),
            dry_run=True,
        )
    assert rejection is None
    mock_postload.assert_called_once()
