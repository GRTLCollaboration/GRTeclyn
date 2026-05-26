from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORTED_WORMHOLE_DIR = REPO_ROOT / "Examples" / "SupportedWormholeCollapse"
DEFAULT_TEMPLATE = SUPPORTED_WORMHOLE_DIR / "params_2gpu.txt"


@dataclass(frozen=True)
class ExecutableConfig:
    """Resolve the SupportedWormhole executable and launch mode."""

    path: Path
    mpi_ranks: int = 1

    @property
    def uses_mpi(self) -> bool:
        return self.mpi_ranks > 1


def default_executable_name(*, mpi_ranks: int = 1, comp: str = "gnu", cuda: bool = True, debug: bool = False) -> str:
    parts = ["main3d", comp]
    if mpi_ranks > 1:
        parts.append("MPI")
    if debug:
        parts.append("DEBUG")
    if cuda:
        parts.append("CUDA")
    parts.append("ex")
    return ".".join(parts)


def resolve_executable(
    executable: str | Path | None = None,
    *,
    mpi_ranks: int = 1,
    comp: str = "gnu",
    cuda: bool = True,
    debug: bool = False,
) -> ExecutableConfig:
    """Return an executable config without requiring the binary to exist yet."""

    if executable is None:
        path = SUPPORTED_WORMHOLE_DIR / default_executable_name(
            mpi_ranks=mpi_ranks,
            comp=comp,
            cuda=cuda,
            debug=debug,
        )
    else:
        path = Path(executable).expanduser()
        if not path.is_absolute():
            path = (REPO_ROOT / path).resolve()

    return ExecutableConfig(path=path.resolve(), mpi_ranks=int(mpi_ranks))


def default_runs_dir() -> Path:
    return REPO_ROOT / "runs"
