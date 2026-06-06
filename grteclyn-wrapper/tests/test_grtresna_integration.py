"""Tests for the GRTresna exotic-matter wrapper integration."""

from __future__ import annotations

import json
import importlib.util
import random
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from grteclyn_wrapper.grtresna.io import _reflect_half_z_to_full
from grteclyn_wrapper.grtresna.domain import GRTresnaDomainConfig
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, write_grtresna_params
from grteclyn_wrapper.core.config import resolve_example
from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.core.plot_consumer import build_consume_command
from grteclyn_wrapper.core.params import write_params
from grteclyn_wrapper.metrics.episode_metrics import (
    ComovingMetrics,
    CurvatureInvariantMetrics,
    EpisodeMetrics,
)
from grteclyn_wrapper.metrics.ftl_general import GeneralFtlReport
from grteclyn_wrapper.metrics.null_geodesic import GeodesicFtlReport
from grteclyn_wrapper.metrics.score import score_episode
from grteclyn_wrapper.search.optimize import (
    _force_exotic_template,
    _grtresna_convergence_rejection_reason,
    _grtresna_rejection_fitness,
    _load_warm_start_vectors,
    build_search_space,
    build_grtresna_config,
    GRTRESNA_REJECTION_BASE_FITNESS,
    GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
)

_ROTATING_ID_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "wormhole"
    / "make_rotating_wormhole_id.py"
)
_ROTATING_ID_SPEC = importlib.util.spec_from_file_location(
    "make_rotating_wormhole_id", _ROTATING_ID_SCRIPT
)
assert _ROTATING_ID_SPEC is not None and _ROTATING_ID_SPEC.loader is not None
_ROTATING_ID_MODULE = importlib.util.module_from_spec(_ROTATING_ID_SPEC)
_ROTATING_ID_SPEC.loader.exec_module(_ROTATING_ID_MODULE)
build_config = _ROTATING_ID_MODULE.build_config


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


def test_grtresna_domain_policy_half_z_defaults() -> None:
    domain = GRTresnaDomainConfig(full_z=False, l_full=80.0, n_full=96, grtresna_nx=48)
    cfg = domain.apply_to_solver(GRTresnaConfig())

    assert cfg.N == (48, 64, 24)
    assert cfg.lo_boundary == (0, 0, 1)
    assert cfg.target_center == (40.0, 40.0, 0.0)
    assert domain.evolution_overrides() == {
        "L_full": 80.0,
        "N_full": 96,
        "center": "40 40 0",
        "lo_boundary": "1 1 2",
    }


def test_grtresna_domain_policy_full_z_is_configurable() -> None:
    domain = GRTresnaDomainConfig(
        full_z=True,
        l_full=96.0,
        n_full=128,
        grtresna_l=192.0,
        grtresna_nx=80,
        grtresna_ny=72,
        grtresna_nz=64,
        gridinit_nx=96,
        gridinit_ny=96,
        gridinit_nz=128,
    )
    cfg = domain.apply_to_solver(GRTresnaConfig())

    assert cfg.N == (80, 72, 64)
    assert cfg.L == 192.0
    assert cfg.lo_boundary == (0, 0, 0)
    assert cfg.target_center == (48.0, 48.0, 48.0)
    assert (cfg.gridinit_nx, cfg.gridinit_ny, cfg.gridinit_nz) == (96, 96, 128)
    assert domain.evolution_overrides() == {
        "L_full": 96.0,
        "N_full": 128,
        "center": "48 48 48",
        "lo_boundary": "1 1 1",
    }


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


def test_rotating_wormhole_id_full_z_config() -> None:
    class Args:
        ranks = 2
        iterations = 30
        timeout = 1800
        nx = 64
        ny = 64
        nz = 64
        length = 128.0
        max_level = 2
        block_factor = 16
        max_grid_size = 16
        refine_threshold = 0.5
        regrid_radius = 0.0
        bh_mass = 0.0
        scalar_mass = 0.1
        amp = 0.2
        width = 8.0
        velocity_x = 0.0
        velocity_y = 0.0
        velocity_z = 0.0
        mode = 2
        psi_relaxation = 0.6
        psi_floor = 0.1
        maximal_jacobian_cap = 25.0
        gridinit_nx = 128
        gridinit_ny = 128
        gridinit_nz = 128
        full_z = True
        target_center_x = 32.0
        target_center_y = 32.0
        target_center_z = None
        keep_hdf5 = False

    cfg = build_config(Args, omega=0.05)

    assert cfg.N == (64, 64, 64)
    assert cfg.lo_boundary == (0, 0, 0)
    assert cfg.target_center == (32.0, 32.0, 64.0)
    assert cfg.lumps[0]["omega"] == 0.05


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
    frames_out_idx = command.index("--frames-out")
    assert command[frames_out_idx + 1] == str(episode.frames_dir)


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


def test_sub_luminal_precursor_keeps_curvature_gradient() -> None:
    report = GeneralFtlReport(
        f_op=0.0,
        t_min=None,
        t_flat=1.0,
        max_local_speed=0.95,
        superluminal_fraction=0.0,
        path_offaxis=False,
        reachable=True,
        notes=(),
    )

    def metrics_with_curvature(max_l2: float) -> EpisodeMetrics:
        return EpisodeMetrics(
            collapse=None,
            constraints=None,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="test",
            curvature=CurvatureInvariantMetrics(
                final_time=2.0,
                max_abs_ricci_scalar=0.0,
                max_ricci_tensor_sq=0.0,
                max_kij_sq=0.0,
                max_l2_ricci_scalar=max_l2,
            ),
            general_ftl_evolved=report,
        )

    weak = score_episode(metrics_with_curvature(1.0), objective_mode="ftl_first")
    stronger = score_episode(metrics_with_curvature(4.0), objective_mode="ftl_first")

    assert 0.0 < weak.components["curvature_activity"] < stronger.components["curvature_activity"] < 1.0
    assert 0.0 < weak.components["ftl_precursor"] < stronger.components["ftl_precursor"]


def test_ftl_first_rewards_shift_drive_before_shortcut_exists() -> None:
    def metrics_with_shift(max_shift: float) -> EpisodeMetrics:
        report = GeneralFtlReport(
            f_op=0.0,
            t_min=None,
            t_flat=1.0,
            max_local_speed=0.96,
            superluminal_fraction=0.0,
            path_offaxis=False,
            reachable=True,
            notes=(),
            max_shift=max_shift,
        )
        return EpisodeMetrics(
            collapse=None,
            constraints=None,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="test",
            general_ftl_evolved=report,
        )

    weak = score_episode(metrics_with_shift(0.03), objective_mode="ftl_first")
    stronger = score_episode(metrics_with_shift(0.30), objective_mode="ftl_first")

    assert 0.0 < weak.components["shift_drive"] < stronger.components["shift_drive"]
    assert stronger.total > weak.total


def _metrics_with_general_ftl(
    *,
    t_min: float | None,
    t_flat: float,
    max_local_speed: float,
    superluminal_fraction: float,
    max_shift: float,
    max_l2_ricci: float = 2.0,
    f_op: float = 0.0,
    stationary: bool = False,
) -> EpisodeMetrics:
    report = GeneralFtlReport(
        f_op=f_op,
        t_min=t_min,
        t_flat=t_flat,
        max_local_speed=max_local_speed,
        superluminal_fraction=superluminal_fraction,
        path_offaxis=False,
        reachable=t_min is not None,
        notes=(),
        max_shift=max_shift,
    )
    return EpisodeMetrics(
        collapse=None,
        constraints=None,
        stability=None,
        ftl=None,
        termination_reason="test",
        comoving=ComovingMetrics(
            beta_mean=0.0,
            delta_comoving=None,
            score=None,
            stationary=True,
        ) if stationary else None,
        curvature=CurvatureInvariantMetrics(
            final_time=2.0,
            max_abs_ricci_scalar=0.0,
            max_ricci_tensor_sq=0.0,
            max_kij_sq=0.0,
            max_l2_ricci_scalar=max_l2_ricci,
        ),
        general_ftl_evolved=report,
    )


def test_channel_progress_zero_on_flat_space() -> None:
    flat = _metrics_with_general_ftl(
        t_min=15.75,
        t_flat=15.75,
        max_local_speed=1.0,
        superluminal_fraction=0.0,
        max_shift=0.0,
        max_l2_ricci=0.0,
    )
    score = score_episode(flat, objective_mode="ftl_first")
    assert score.components["channel_progress"] == 0.0


def test_channel_progress_prefers_coupled_candidates_over_single_mechanism() -> None:
    t_flat = 15.75
    precursor_only = _metrics_with_general_ftl(
        t_min=15.87,
        t_flat=t_flat,
        max_local_speed=1.043,
        superluminal_fraction=0.56,
        max_shift=0.02,
    )
    shift_only = _metrics_with_general_ftl(
        t_min=21.87,
        t_flat=t_flat,
        max_local_speed=0.933,
        superluminal_fraction=0.0,
        max_shift=0.147,
    )
    coupled = _metrics_with_general_ftl(
        t_min=15.80,
        t_flat=t_flat,
        max_local_speed=1.02,
        superluminal_fraction=0.20,
        max_shift=0.10,
    )

    prec_score = score_episode(precursor_only, objective_mode="ftl_first")
    shift_score = score_episode(shift_only, objective_mode="ftl_first")
    coupled_score = score_episode(coupled, objective_mode="ftl_first")

    assert coupled_score.components["channel_progress"] > prec_score.components["channel_progress"]
    assert coupled_score.components["channel_progress"] > shift_score.components["channel_progress"]
    assert coupled_score.total > shift_score.total


def test_channel_progress_ranks_eval128_style_below_coupled_but_above_shift_only() -> None:
    t_flat = 15.75
    eval128_style = _metrics_with_general_ftl(
        t_min=15.87016792210433,
        t_flat=t_flat,
        max_local_speed=1.0432600860953292,
        superluminal_fraction=0.5585601826813292,
        max_shift=0.02041905875701828,
    )
    eval57_style = _metrics_with_general_ftl(
        t_min=21.869045287604937,
        t_flat=t_flat,
        max_local_speed=0.932994638314321,
        superluminal_fraction=0.0,
        max_shift=0.14652009373502575,
    )

    s128 = score_episode(eval128_style, objective_mode="ftl_first")
    s57 = score_episode(eval57_style, objective_mode="ftl_first")

    assert s128.components["channel_progress"] > s57.components["channel_progress"]
    assert s128.components["ftl_precursor"] > s57.components["ftl_precursor"]
    assert s57.components["shift_drive"] > s128.components["shift_drive"]


def test_weak_stationary_shortcut_is_not_scored_as_operational_ftl() -> None:
    weak_lens = _metrics_with_general_ftl(
        t_min=15.55,
        t_flat=15.75,
        max_local_speed=1.196,
        superluminal_fraction=1.0,
        max_shift=0.038,
        f_op=0.0127,
        stationary=True,
    )

    score = score_episode(weak_lens, objective_mode="ftl_first")

    assert score.components["operational_ftl"] == 0.0
    assert score.components["channel_progress"] > 0.0
    assert score.components["stationary_artifact_penalty"] < 0.0
    assert score.total < 1000.0


def test_strong_evolved_shortcut_gets_operational_ftl_reward() -> None:
    strong = _metrics_with_general_ftl(
        t_min=14.8,
        t_flat=15.75,
        max_local_speed=1.25,
        superluminal_fraction=0.2,
        max_shift=0.12,
        f_op=0.060,
    )

    score = score_episode(strong, objective_mode="ftl_first")

    assert score.components["operational_ftl"] > 0.0
    assert score.components["stationary_artifact_penalty"] == 0.0
    # Coordinate FTL alone (no gauge-invariant geodesic confirmation) is
    # deliberately demoted: the dominant FTL budget now sits on
    # operational_ftl_geodesic, which is unset here.  A strong coordinate
    # shortcut still earns a large positive score, just below the old
    # coordinate-only weighting.
    assert score.total > 400.0


def _geodesic_report(f_geo: float) -> GeodesicFtlReport:
    return GeodesicFtlReport(
        f_geo=f_geo,
        t_min=15.0 if f_geo > 0 else 15.75,
        t_flat=15.75,
        n_rays=5,
        n_reached=5,
        max_h_drift=1.0e-8,
        h_quality_ok=True,
    )


def test_geodesic_confirmation_dominates_coordinate_shortcut() -> None:
    base = _metrics_with_general_ftl(
        t_min=14.8,
        t_flat=15.75,
        max_local_speed=1.25,
        superluminal_fraction=0.2,
        max_shift=0.12,
        f_op=0.060,
    )
    coord_only = score_episode(base, objective_mode="ftl_first")

    confirmed = EpisodeMetrics(
        collapse=base.collapse,
        constraints=base.constraints,
        stability=base.stability,
        comoving=base.comoving,
        ftl=base.ftl,
        termination_reason=base.termination_reason,
        curvature=base.curvature,
        general_ftl_evolved=base.general_ftl_evolved,
        geodesic_ftl=_geodesic_report(0.05),
    )
    geodesic_score = score_episode(confirmed, objective_mode="ftl_first")

    assert geodesic_score.components["operational_ftl_geodesic"] > 0.0
    assert geodesic_score.total > coord_only.total + 500.0


def test_geodesic_zero_flags_gauge_artifact() -> None:
    artifact = EpisodeMetrics(
        collapse=None,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="test",
        general_ftl_evolved=GeneralFtlReport(
            f_op=0.060,
            t_min=14.8,
            t_flat=15.75,
            max_local_speed=1.25,
            superluminal_fraction=0.2,
            path_offaxis=False,
            reachable=True,
            notes=(),
            max_shift=0.12,
        ),
        geodesic_ftl=_geodesic_report(0.0),
    )
    score = score_episode(artifact, objective_mode="weighted")
    assert score.components["operational_ftl_geodesic"] == 0.0
    assert any("gauge artifact" in note for note in score.notes)


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
    crash_fitness = GRTRESNA_REJECTION_BASE_FITNESS + GRTRESNA_REJECTION_MAX_EXTRA_FITNESS
    assert crash_fitness >= very_bad
    nonfinite = _grtresna_rejection_fitness({
        "iteration": 29,
        "ham_pct": float("nan"),
        "mom_pct": 1.0,
    })
    assert mildly_bad < very_bad < nonfinite
