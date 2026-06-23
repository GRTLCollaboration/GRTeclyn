from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

gym = pytest.importorskip("gymnasium")

from grteclyn_wrapper.rl.env import (  # noqa: E402
    GAUGE_ACTION_DIM,
    GLOBAL_OBS_DIM,
    LUMP_ACTION_STRIDE,
    LUMP_OBS_STRIDE,
    SpacetimeFtlEnv,
    SpacetimeFtlEnvConfig,
)


class _StubIO:
    """Minimal ObservationSource + ActionSink stub (no socket / subprocess)."""

    def recv_obs(self) -> np.ndarray:
        return np.zeros(1, dtype=np.float64)

    def send_action(self, action: np.ndarray) -> None:  # pragma: no cover
        pass


def _make_env(num_lumps: int) -> SpacetimeFtlEnv:
    cfg = SpacetimeFtlEnvConfig(
        episode_path=Path("/tmp/rl_env_dims_test"),
        executable=Path("/bin/true"),
        params_file=Path("/tmp/params.txt"),
        num_lumps=num_lumps,
    )
    stub = _StubIO()
    return SpacetimeFtlEnv(cfg, observation_source=stub, action_sink=stub)


@pytest.mark.parametrize("num_lumps", [1, 2, 5])
def test_spaces_scale_with_num_lumps(num_lumps: int) -> None:
    env = _make_env(num_lumps)
    assert env.action_space.shape == (LUMP_ACTION_STRIDE * num_lumps + GAUGE_ACTION_DIM,)
    assert env.observation_space.shape == (GLOBAL_OBS_DIM + LUMP_OBS_STRIDE * num_lumps,)


def test_map_action_is_passthrough_clip() -> None:
    env = _make_env(2)
    raw = np.array([2.0, -3.0, 0.5, 0.1, -0.2, 0.3, 1.5, -1.5], dtype=np.float64)
    mapped = env._map_action(raw)
    assert mapped.shape == raw.shape
    assert mapped.min() >= -1.0 and mapped.max() <= 1.0
    # values already in range pass through unchanged
    assert mapped[2] == pytest.approx(0.5)
