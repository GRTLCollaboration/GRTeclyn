"""Integration checks for the shared GRTresna pre-evolution gate path."""

from __future__ import annotations

from pathlib import Path

import pytest

from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.grtresna.matter_wiring import evolution_overrides_from_config
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.projection.postload_gate import PostLoadGateConfig
from grteclyn_wrapper.search.grtresna_evaluation_gates import (
    GRTresnaPreEvolutionGateConfig,
    apply_grtresna_pre_evolution_gates,
)


def test_apply_pre_evolution_gates_accepts_good_convergence_and_dry_postload(
    tmp_path: Path,
) -> None:
    episode = create_episode(tmp_path, name="gate_episode")
    gridinit = tmp_path / "initial_data.gridinit"
    gridinit.write_bytes(b"placeholder")

    cfg = GRTresnaConfig(
        lumps=[{"amp": 0.1, "width": 4.0, "center": (0.0, 0.0, 0.0), "exotic": 0}],
        scalar_mass=0.1,
    )
    gte_overrides = {
        "recipe_initial_data_file": str(gridinit),
        **evolution_overrides_from_config(cfg),
    }

    rejection = apply_grtresna_pre_evolution_gates(
        episode=episode,
        convergence={"ham_pct": 1.0, "mom_pct": 1.0},
        gridinit_path=gridinit,
        gte_overrides=gte_overrides,
        cuda_devices=None,
        config=GRTresnaPreEvolutionGateConfig(
            postload_enabled=True,
            postload_config=PostLoadGateConfig(),
        ),
        dry_run=True,
    )
    assert rejection is None
