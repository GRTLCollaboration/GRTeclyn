"""Tests for the GRTresna exotic-matter wrapper integration."""

from __future__ import annotations

import json
import random
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from grteclyn_wrapper.grtresna.io import _reflect_half_z_to_full
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, write_grtresna_params
from grteclyn_wrapper.core.config import resolve_example
from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.core.plot_consumer import build_consume_command
from grteclyn_wrapper.core.params import write_params
from grteclyn_wrapper.metrics.episode_metrics import EpisodeMetrics
from grteclyn_wrapper.metrics.ftl_general import GeneralFtlReport
from grteclyn_wrapper.metrics.score import score_episode
from grteclyn_wrapper.search.optimize import (
    _force_exotic_template,
    _grtresna_convergence_rejection_reason,
    _grtresna_rejection_fitness,
    _load_warm_start_vectors,
    build_search_space,
    build_grtresna_config,
)


def test_exotic_indexed_lump_enables_amr_safe_maximal_slicing() -> None:
    base = GRTresnaConfig(
        dphi=0.0,
        dpi=0.0,
        regrid_radius=0.0,
        bh1_bare_mass=0.0,
        bh1_spin=(0.0, 0.0, 0.0),
    )
    cfg = build_grtresna_config(
        {
            "grtresna_lump0_amp": 0.2,
            "grtresna_lump0_width": 7.0,
            "grtresna_lump0_center_x": -5.0,
            "grtresna_lump0_velocity_x": 0.3,
            "grtresna_lump0_mode": 1.8,
            "grtresna_lump0_exotic": 1.0,
            "grtresna_lump1_amp": 0.15,
            "grtresna_lump1_exotic": 0.0,
        },
        base,
    )

    assert cfg.maximal_slicing is True
    assert cfg.psi_relaxation == 0.6
    assert cfg.psi_floor == 0.1
    assert cfg.maximal_jacobian_cap == 25.0
    assert cfg.coefficient_average_type == "arithmetic"
    assert cfg.dphi == 0.0
    assert cfg.dpi == 0.0
    assert len(cfg.lumps) == 2
    assert cfg.lumps[0]["exotic"] == 1
    assert cfg.lumps[0]["mode"] == 2
    assert cfg.lumps[0]["center"] == (-5.0, 0.0, 0.0)
    assert cfg.lumps[0]["velocity"] == (0.3, 0.0, 0.0)
    assert cfg.lumps[1]["exotic"] == 0


def test_canonical_lump_keeps_standard_solver_path() -> None:
    cfg = build_grtresna_config(
        {
            "grtresna_lump0_amp": 0.1,
            "grtresna_lump0_velocity_x": 0.2,
            "grtresna_lump0_exotic": 0.0,
        },
        GRTresnaConfig(coefficient_average_type="harmonic"),
    )

    assert cfg.maximal_slicing is False
    assert cfg.psi_relaxation == 1.0
    assert cfg.psi_floor == -1.0
    assert cfg.maximal_jacobian_cap == -1.0
    assert cfg.coefficient_average_type == "harmonic"


def test_write_grtresna_params_renders_exotic_amr_knobs() -> None:
    cfg = GRTresnaConfig(
        max_level=3,
        refine_threshold=0.3,
        regrid_radius=0.0,
        coefficient_average_type="arithmetic",
        dphi=0.0,
        dpi=0.0,
        maximal_slicing=True,
        psi_relaxation=0.5,
        psi_floor=0.1,
        maximal_jacobian_cap=20.0,
        lumps=[
            {
                "amp": 0.2,
                "width": 8.0,
                "center": (0.0, 0.0, 0.0),
                "velocity": (0.25, 0.0, 0.0),
                "omega": 0.0,
                "mode": 0,
                "exotic": 1,
            }
        ],
    )

    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "params.txt"
        write_grtresna_params(cfg, path)
        text = path.read_text(encoding="utf-8")

    assert "max_level = 3" in text
    assert "refine_threshold = 0.3" in text
    assert "regrid_radius = 0.0" in text
    assert "coefficient_average_type = arithmetic" in text
    assert "dphi = 0.0" in text
    assert "dpi = 0.0" in text
    assert "num_lumps = 1" in text
    assert "lump0_exotic = 1" in text
    assert "maximal_slicing = 1" in text
    assert "psi_relaxation = 0.5" in text
    assert "psi_floor = 0.1" in text
    assert "maximal_jacobian_cap = 20.0" in text


def test_rotating_wormhole_example_is_registered() -> None:
    example = resolve_example("RotatingWormholeCollapse")

    assert example.template.name == "params_rotating.txt"
    assert example.check_prefix == "RotatingWormholeChk"
    assert example.plot_prefix == "RotatingWormholePlt"


def test_rotating_exotic_config_renders_no_kick_momentum_source() -> None:
    cfg = GRTresnaConfig(
        bh1_bare_mass=0.0,
        dphi=0.0,
        dpi=0.0,
        scalar_mass=0.0,
        lump_amp=0.08,
        lump_width=4.0,
        lump_omega=0.1,
        lump_mode=2,
        lump_exotic=1,
        maximal_slicing=True,
        psi_relaxation=0.6,
        psi_floor=0.1,
        maximal_jacobian_cap=25.0,
        coefficient_average_type="arithmetic",
    )

    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "params.txt"
        write_grtresna_params(cfg, path)
        text = path.read_text(encoding="utf-8")

    assert "dphi = 0.0" in text
    assert "dpi = 0.0" in text
    assert "lump_omega = 0.1" in text
    assert "lump_mode = 2" in text
    assert "lump_exotic = 1" in text
    assert "maximal_slicing = 1" in text


def test_rotating_wormhole_consumer_uses_psi4_profile() -> None:
    example = resolve_example("RotatingWormholeCollapse")
    with TemporaryDirectory() as tmp:
        episode = create_episode(Path(tmp), name="rotating-consumer-test")
        write_params(
            example.template,
            episode.params_path,
            episode_dir=episode.path,
            example=example,
            overrides={"recipe_initial_data_file": "/tmp/initial_data.gridinit"},
        )

        command = build_consume_command(
            episode,
            profile="wormhole",
            radii=(12.0, 16.0, 20.0, 24.0),
        )

    assert "--psi4" in command
    assert "--embedding" in command
    assert "--frames-fields" in command


def test_half_z_reflection_places_positive_z_data_in_upper_half() -> None:
    # The loader maps GRTeclyn's z >= 0 half-domain onto the upper half of the
    # full reflected .gridinit. If Chombo's stored z >= 0 slab is left in the
    # lower half, the z=0 slice samples the far-boundary tail instead of matter.
    data = np.zeros((4, 1, 1, 2), dtype=float)
    data[:2, 0, 0, 0] = [10.0, 20.0]  # even component, positive z
    data[:2, 0, 0, 1] = [1.0, 2.0]    # z-odd component, positive z

    _reflect_half_z_to_full(data, ["chi", "shift3"])

    assert data[:, 0, 0, 0].tolist() == [20.0, 10.0, 10.0, 20.0]
    assert data[:, 0, 0, 1].tolist() == [-2.0, -1.0, 1.0, 2.0]


def test_grtresna_search_space_defaults_to_compact_five_lump_basis() -> None:
    dims = build_search_space(grtresna=True)

    assert len(dims) == 55
    assert len([d for d in dims if d.param_key.endswith("_exotic")]) == 5
    width_dims = [d for d in dims if d.param_key.endswith("_width")]
    center_x_dims = [d for d in dims if d.param_key.endswith("_center_x")]
    center_z_dims = [d for d in dims if d.param_key.endswith("_center_z")]
    assert all(d.upper <= 4.5 for d in width_dims)
    assert all(d.lower >= -6.0 and d.upper <= 6.0 for d in center_x_dims)
    assert all(d.lower >= -3.0 and d.upper <= 3.0 for d in center_z_dims)


def test_forced_exotic_template_sets_discrete_exotic_lumps() -> None:
    dims = build_search_space(grtresna=True, grtresna_lumps=5)
    vec = [d.center for d in dims]

    injected = _force_exotic_template(vec, dims, random.Random(7), template_index=1)
    overrides = {dim.param_key: injected[i] for i, dim in enumerate(dims)}

    rounded = [
        round(overrides[f"grtresna_lump{k}_exotic"])
        for k in range(5)
    ]
    assert any(rounded)
    assert all(abs(overrides[f"grtresna_lump{k}_center_x"]) <= 6.0 for k in range(5))
    assert all(1.5 <= overrides[f"grtresna_lump{k}_width"] <= 4.5 for k in range(5))


def test_load_warm_start_vectors_uses_top_trajectory_scores() -> None:
    dims = build_search_space(grtresna=True, grtresna_lumps=1)
    with TemporaryDirectory() as tmp:
        traj = Path(tmp) / "trajectory.jsonl"
        low = {"grtresna_lump0_amp": 0.05, "grtresna_lump0_width": 4.0}
        high = {"grtresna_lump0_amp": 0.2, "grtresna_lump0_width": 6.0}
        traj.write_text(
            "\n".join([
                json.dumps({"score": 1.0, "overrides": low}),
                json.dumps({"score": 9.0, "overrides": high}),
            ])
            + "\n",
            encoding="utf-8",
        )

        vectors = _load_warm_start_vectors([traj], dims, top_k=1)

    values = {dim.param_key: vectors[0][i] for i, dim in enumerate(dims)}
    assert values["grtresna_lump0_amp"] == 0.2
    assert values["grtresna_lump0_width"] == 4.5


def test_ftl_first_score_prioritizes_ftl_precursor_over_health() -> None:
    def metrics_with_precursor(value: float) -> EpisodeMetrics:
        report = GeneralFtlReport(
            f_op=0.0,
            t_min=None,
            t_flat=1.0,
            max_local_speed=1.0 + value,
            superluminal_fraction=value,
            path_offaxis=False,
            reachable=True,
            notes=(),
        )
        return EpisodeMetrics(
            collapse=None,
            constraints=None,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="test",
            general_ftl=report,
            general_ftl_evolved=report,
        )

    no_ftl = score_episode(metrics_with_precursor(0.0), objective_mode="ftl_first")
    some_ftl = score_episode(metrics_with_precursor(0.03), objective_mode="ftl_first")

    assert some_ftl.total > no_ftl.total
    assert some_ftl.components["ftl_precursor"] > 0.0


def test_bad_grtresna_convergence_is_rejected_before_evolution() -> None:
    assert _grtresna_convergence_rejection_reason(None) is not None
    assert _grtresna_convergence_rejection_reason({
        "iteration": 29,
        "ham_pct": float("nan"),
        "mom_pct": 1.0,
    }) is not None
    assert _grtresna_convergence_rejection_reason({
        "iteration": 29,
        "ham_pct": 100.0,
        "mom_pct": 79.0,
    }) is not None
    assert _grtresna_convergence_rejection_reason({
        "iteration": 29,
        "ham_pct": 0.02,
        "mom_pct": 0.01,
    }) is None

    mildly_bad = _grtresna_rejection_fitness({
        "iteration": 29,
        "ham_pct": 8.0,
        "mom_pct": 7.0,
    })
    very_bad = _grtresna_rejection_fitness({
        "iteration": 29,
        "ham_pct": 100.0,
        "mom_pct": 79.0,
    })
    nonfinite = _grtresna_rejection_fitness({
        "iteration": 29,
        "ham_pct": float("nan"),
        "mom_pct": 1.0,
    })
    assert mildly_bad < very_bad < nonfinite
