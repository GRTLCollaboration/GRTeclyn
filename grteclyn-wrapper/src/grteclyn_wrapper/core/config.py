from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


WRAPPER_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class ExampleConfig:
    """GRTeclyn example directory and default params template."""

    name: str
    dir: Path
    template: Path
    check_prefix: str
    plot_prefix: str


def _looks_like_grteclyn(path: Path) -> bool:
    return (
        (path / "Examples" / "SupportedWormholeCollapse" / "params_2gpu.txt").exists()
        or (path / "Examples" / "RadialRecipe" / "params.txt").exists()
    )


def resolve_repo_root(explicit: str | Path | None = None) -> Path:
    """Resolve the GRTeclyn checkout root for a standalone wrapper install."""

    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(Path(explicit).expanduser())
    if os.environ.get("GRTECLYN_ROOT"):
        candidates.append(Path(os.environ["GRTECLYN_ROOT"]).expanduser())
    candidates.extend([
        Path.cwd(),
        WRAPPER_ROOT.parent,
        WRAPPER_ROOT.parent.parent,
    ])

    for candidate in candidates:
        resolved = candidate.resolve()
        if _looks_like_grteclyn(resolved):
            return resolved

    raise FileNotFoundError(
        "Could not find a GRTeclyn checkout. Set GRTECLYN_ROOT=/path/to/GRTeclyn "
        "or run from the GRTeclyn repository root."
    )


REPO_ROOT = resolve_repo_root()
SUPPORTED_WORMHOLE_DIR = REPO_ROOT / "Examples" / "SupportedWormholeCollapse"
RADIAL_RECIPE_DIR = REPO_ROOT / "Examples" / "RadialRecipe"
DEFAULT_TEMPLATE = SUPPORTED_WORMHOLE_DIR / "params_2gpu.txt"
DEFAULT_RADIAL_RECIPE_TEMPLATE = RADIAL_RECIPE_DIR / "params.txt"

EXAMPLES: dict[str, ExampleConfig] = {
    "SupportedWormholeCollapse": ExampleConfig(
        name="SupportedWormholeCollapse",
        dir=SUPPORTED_WORMHOLE_DIR,
        template=DEFAULT_TEMPLATE,
        check_prefix="SupportedWormholeChk",
        plot_prefix="SupportedWormholePlt",
    ),
    "RadialRecipe": ExampleConfig(
        name="RadialRecipe",
        dir=RADIAL_RECIPE_DIR,
        template=DEFAULT_RADIAL_RECIPE_TEMPLATE,
        check_prefix="RadialRecipeChk",
        plot_prefix="RadialRecipePlt",
    ),
}


def resolve_example(name: str = "SupportedWormholeCollapse") -> ExampleConfig:
    try:
        return EXAMPLES[name]
    except KeyError as exc:
        known = ", ".join(sorted(EXAMPLES))
        raise ValueError(f"Unknown example {name!r}. Known examples: {known}") from exc


@dataclass(frozen=True)
class ExecutableConfig:
    """Resolve a GRTeclyn example executable and launch mode."""

    path: Path
    example: ExampleConfig
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
    example: str | ExampleConfig = "SupportedWormholeCollapse",
    mpi_ranks: int = 1,
    comp: str = "gnu",
    cuda: bool = True,
    debug: bool = False,
) -> ExecutableConfig:
    """Return an executable config without requiring the binary to exist yet."""

    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)

    if executable is None:
        path = example_cfg.dir / default_executable_name(
            mpi_ranks=mpi_ranks,
            comp=comp,
            cuda=cuda,
            debug=debug,
        )
    else:
        path = Path(executable).expanduser()
        if not path.is_absolute():
            path = (REPO_ROOT / path).resolve()

    return ExecutableConfig(path=path.resolve(), example=example_cfg, mpi_ranks=int(mpi_ranks))


def default_runs_dir() -> Path:
    return REPO_ROOT / "runs"
