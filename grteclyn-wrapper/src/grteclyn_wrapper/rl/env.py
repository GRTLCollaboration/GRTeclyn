from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gymnasium as gym
import numpy as np

from grteclyn_wrapper.rl.audit import compute_audit_penalty, wait_consumer_drain, wait_for_frame_record
from grteclyn_wrapper.rl.params import read_rl_num_lumps
from grteclyn_wrapper.rl.process import SubprocessEpisodeLauncher
from grteclyn_wrapper.rl.protocols import ActionSink, ObservationSource
from grteclyn_wrapper.rl.reward import TaxManState, compute_dense_reward, evaluate_fences
from grteclyn_wrapper.rl.zmq_client import ZmqObservationSource


@dataclass
class SpacetimeFtlEnvConfig:
    episode_path: Path
    executable: Path
    params_file: Path
    zmq_port: int = 5555
    stop_time: float = 30.0
    objective_mode: str = "general_ftl"
    use_taxman_audit: bool = True
    gpu_id: int | None = None
    num_lumps: int | None = None  # None => read rl_num_lumps from params_file


# Per-lump observation stride (must match RL_LUMP_OBS_STRIDE in
# Source/GRTeclynCore/RL/RLLumpState.hpp): x, y, z, size, mass, peak,
# min_lapse, min_chi.
LUMP_OBS_STRIDE = 8
# Global observation scalars: min_chi, min_lapse, max_abs_K, L2_Ham, L2_Mom, t.
GLOBAL_OBS_DIM = 6
# Per-lump actions: amplitude, frequency, phase.  Plus 2 gauge actions.
LUMP_ACTION_STRIDE = 3
GAUGE_ACTION_DIM = 2


class SpacetimeFtlEnv(gym.Env):
    metadata = {"render_modes": []}

    def __init__(
        self,
        config: SpacetimeFtlEnvConfig,
        *,
        observation_source: ObservationSource | None = None,
        action_sink: ActionSink | None = None,
    ) -> None:
        super().__init__()
        self._config = config
        self._taxman = TaxManState()
        self._episode_return = 0.0
        self._terminated = False

        if config.num_lumps is None:
            n = read_rl_num_lumps(config.params_file)
        else:
            n = max(1, int(config.num_lumps))
        self._num_lumps = n
        action_dim = LUMP_ACTION_STRIDE * n + GAUGE_ACTION_DIM
        obs_dim = GLOBAL_OBS_DIM + LUMP_OBS_STRIDE * n
        self.action_space = gym.spaces.Box(
            low=-1.0, high=1.0, shape=(action_dim,), dtype=np.float64
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(obs_dim,), dtype=np.float64
        )

        self._zmq = observation_source or ZmqObservationSource(port=config.zmq_port)
        self._action_sink: ActionSink = action_sink or self._zmq
        self._launcher = SubprocessEpisodeLauncher(
            executable=config.executable,
            params_file=config.params_file,
            work_dir=config.episode_path,
            gpu_id=config.gpu_id,
        )

    def _map_action(self, action: np.ndarray) -> np.ndarray:
        # C++ (RLActionApplier) is the single mapping site: it converts these
        # raw-normalized actions to physical bounds and owns the safety caps.
        # Here we only clip to the policy's declared [-1, 1] range.
        return np.clip(np.asarray(action, dtype=np.float64), -1.0, 1.0)

    def reset(self, *, seed: int | None = None, options: dict[str, Any] | None = None):
        super().reset(seed=seed)
        self._taxman = TaxManState()
        self._episode_return = 0.0
        self._terminated = False
        self._launcher.start()
        obs = self._zmq.recv_obs()
        expected = GLOBAL_OBS_DIM + LUMP_OBS_STRIDE * self._num_lumps
        if obs.size != expected:
            raise RuntimeError(
                f"observation size mismatch: got {obs.size}, expected {expected} "
                f"(num_lumps={self._num_lumps}; check params rl_num_lumps)"
            )
        return obs, {}

    def step(self, action: np.ndarray):
        if self._terminated:
            raise RuntimeError("step() called on terminated episode")

        physical = self._map_action(action)
        self._action_sink.send_action(physical)
        obs = self._zmq.recv_obs()
        min_lapse = float(obs[1])
        l2_ham = float(obs[3])
        sim_time = float(obs[5])

        record = wait_for_frame_record(self._config.episode_path, sim_time)
        ftl_geo = float(record.get("f_geo", 0.0)) if record else 0.0
        wec_violation_fraction = record.get("wec_violation_fraction") if record else None
        if wec_violation_fraction is not None:
            wec_violation_fraction = float(wec_violation_fraction)

        reward = compute_dense_reward(
            ftl_geo=ftl_geo,
            l2_ham=l2_ham,
            min_lapse=min_lapse,
            state=self._taxman,
        )
        fence = evaluate_fences(
            min_lapse=min_lapse,
            l2_ham=l2_ham,
            wec_violation_fraction=wec_violation_fraction,
            horizon_detected=min_lapse < 0.05,
            state=self._taxman,
        )
        reward += fence.penalty
        self._episode_return += reward

        terminated = fence.terminated or sim_time >= self._config.stop_time
        truncated = False

        if terminated and self._config.use_taxman_audit:
            wait_consumer_drain(self._config.episode_path)
            audit = compute_audit_penalty(
                self._config.episode_path,
                sum(self._taxman.dense_history),
                objective_mode=self._config.objective_mode,
                target_stop_time=self._config.stop_time,
            )
            reward += audit
            self._episode_return += audit
            self._launcher.stop()

        if terminated:
            self._terminated = True

        return obs, reward, terminated, truncated, {"sim_time": sim_time, "fence": fence.reason}

    def close(self) -> None:
        self._launcher.stop()
        if hasattr(self._zmq, "close"):
            self._zmq.close()
