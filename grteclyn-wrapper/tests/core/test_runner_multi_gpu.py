from pathlib import Path

import pytest

from grteclyn_wrapper.core.config import ExecutableConfig, resolve_example
from grteclyn_wrapper.core.runner import build_command


def _executable(ranks: int) -> ExecutableConfig:
    return ExecutableConfig(
        path=Path("/tmp/main3d.gnu.MPI.CUDA.ex"),
        example=resolve_example("RadialRecipe"),
        mpi_ranks=ranks,
    )


def test_build_command_binds_one_explicit_gpu_per_mpi_rank() -> None:
    command = build_command(
        _executable(2),
        Path("/tmp/params.txt"),
        cuda_devices="4,5",
    )

    assert command[:4] == ["mpirun", "-n", "2", "bash"]
    wrapper = command[-1]
    assert "GRTECLYN_GPU_IDS" in wrapper
    assert "OMPI_COMM_WORLD_LOCAL_RANK" in wrapper
    assert "CUDA_VISIBLE_DEVICES" in wrapper
    assert "/tmp/main3d.gnu.MPI.CUDA.ex /tmp/params.txt" in wrapper


def test_build_command_rejects_rank_gpu_count_mismatch() -> None:
    with pytest.raises(ValueError, match="MPI ranks"):
        build_command(
            _executable(2),
            Path("/tmp/params.txt"),
            cuda_devices="4",
        )
