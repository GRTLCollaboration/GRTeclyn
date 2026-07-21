"""Plug-and-play site layout for GRTeclyn / GRTresna / Chombo / OpenMPI.

Single place to configure machine paths (no identity in git):

1. Process environment variables (highest priority)
2. Gitignored ``grteclyn-wrapper/.env`` (see ``.env.example``)
3. Legacy fallback ``site.local.env`` if present
4. Portable discovery relative to the GRTeclyn checkout

Consumers depend on :class:`SiteLayout` (or module facades), never on hardcoded
absolute paths. On a new machine: ``cp .env.example .env`` and edit.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, MutableMapping, Protocol, Sequence, runtime_checkable

WRAPPER_ROOT = Path(__file__).resolve().parents[3]
DOTENV_PATH = WRAPPER_ROOT / ".env"
_LEGACY_SITE_LOCAL = WRAPPER_ROOT / "site.local.env"
_ENV_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}|\$([A-Za-z_][A-Za-z0-9_]*)")


@runtime_checkable
class SiteLayout(Protocol):
    """Dependency-inverted view of the simulation site tree."""

    def grteclyn_root(self) -> Path: ...

    def sim_root(self) -> Path: ...

    def grtresna_root(self) -> Path: ...

    def chombo_home(self) -> Path: ...

    def openmpi_root(self) -> Path: ...

    def grtresna_env(self) -> Path | None: ...


def parse_env_file(text: str) -> dict[str, str]:
    """Parse ``KEY=VALUE`` / ``export KEY=VALUE`` lines (comments ignored)."""

    out: dict[str, str] = {}
    for line in text.splitlines():
        raw = line.strip()
        if not raw or raw.startswith("#"):
            continue
        if raw.startswith("export "):
            raw = raw[len("export ") :].strip()
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        key = key.strip()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
            value = value[1:-1]
        if key:
            out[key] = value
    return out


def expand_env_value(value: str, lookup: Mapping[str, str]) -> str:
    """Expand ``$VAR`` / ``${VAR}`` using *lookup* (and leave unknowns unchanged)."""

    def _replace(match: re.Match[str]) -> str:
        name = match.group(1) or match.group(2)
        return lookup.get(name, match.group(0))

    previous = None
    current = value
    # Allow nested refs a few times (SIM_ROOT → GRTECLYN_ROOT).
    for _ in range(8):
        if current == previous:
            break
        previous = current
        current = _ENV_VAR_PATTERN.sub(_replace, current)
    return os.path.expanduser(current)


@dataclass
class EnvOverlay:
    """Load gitignored dotenv overlay(s) into an environment mapping.

    Existing keys in the target mapping are never overwritten (plug-in, not clobber).
    """

    paths: Sequence[Path] = field(
        default_factory=lambda: (DOTENV_PATH, _LEGACY_SITE_LOCAL)
    )
    _applied: bool = False

    def apply(
        self,
        env: MutableMapping[str, str] | None = None,
        *,
        force: bool = False,
    ) -> Path | None:
        target: MutableMapping[str, str] = os.environ if env is None else env
        if self._applied and not force:
            for path in self.paths:
                if path.is_file():
                    return path
            return None
        self._applied = True

        loaded_from: Path | None = None
        pending: dict[str, str] = {}
        for path in self.paths:
            if not path.is_file():
                continue
            pending.update(parse_env_file(path.read_text(encoding="utf-8")))
            loaded_from = path
            break  # first existing overlay wins (.env preferred over legacy)

        if not pending:
            return None

        lookup: dict[str, str] = {**dict(target), **pending}
        for key, raw in pending.items():
            if key in target:
                continue
            target[key] = expand_env_value(raw, lookup)
            lookup[key] = target[key]
        return loaded_from


_DEFAULT_OVERLAY = EnvOverlay()


def _env_path(name: str, env: Mapping[str, str]) -> Path | None:
    raw = env.get(name)
    if not raw:
        return None
    return Path(expand_env_value(raw, env)).resolve()


@dataclass(frozen=True)
class DefaultSiteLayout:
    """Portable layout: dotenv/env → sibling discovery under ``sim_root``."""

    env: Mapping[str, str] | None = None
    overlay: EnvOverlay | None = None

    def _environ(self) -> Mapping[str, str]:
        if self.overlay is not None:
            self.overlay.apply()
        elif self.env is None:
            _DEFAULT_OVERLAY.apply()
        return os.environ if self.env is None else self.env

    def grteclyn_root(self) -> Path:
        environ = self._environ()
        explicit = _env_path("GRTECLYN_ROOT", environ)
        if explicit is not None:
            return explicit
        from .config import resolve_repo_root

        return resolve_repo_root()

    def sim_root(self) -> Path:
        environ = self._environ()
        explicit = _env_path("SIM_ROOT", environ)
        if explicit is not None:
            return explicit
        return self.grteclyn_root().parent

    def grtresna_root(self) -> Path:
        environ = self._environ()
        explicit = _env_path("GRTRESNA_ROOT", environ)
        if explicit is not None:
            return explicit
        return self.sim_root() / "GRTresna"

    def chombo_home(self) -> Path:
        environ = self._environ()
        explicit = _env_path("CHOMBO_HOME", environ)
        if explicit is not None:
            return explicit
        return self.sim_root() / "Chombo" / "lib"

    def openmpi_root(self) -> Path:
        environ = self._environ()
        explicit = _env_path("OPENMPI_ROOT", environ)
        if explicit is not None:
            return explicit
        # Prefer OPENMPI_ROOT from .env; otherwise a conventional sibling prefix.
        return self.sim_root() / "local" / "openmpi"

    def grtresna_env(self) -> Path | None:
        """Return ``GRTRESNA_ENV`` if set (via process env or ``.env``). No host guesses."""

        environ = self._environ()
        explicit = _env_path("GRTRESNA_ENV", environ)
        if explicit is not None:
            return explicit if explicit.is_dir() else None
        return None

    def grtresna_env_bin(self) -> Path | None:
        env = self.grtresna_env()
        if env is None:
            return None
        bin_dir = env / "bin"
        return bin_dir if bin_dir.is_dir() else None

    def grtresna_env_lib(self) -> Path | None:
        env = self.grtresna_env()
        if env is None:
            return None
        lib_dir = env / "lib"
        return lib_dir if lib_dir.is_dir() else None


_DEFAULT_LAYOUT = DefaultSiteLayout()


def get_site_layout(layout: SiteLayout | None = None) -> SiteLayout:
    """Return the active layout (injectable for tests / alternate sites)."""

    return _DEFAULT_LAYOUT if layout is None else layout


def load_dotenv(*, force: bool = False) -> Path | None:
    """Apply gitignored ``.env`` (or legacy ``site.local.env``) into ``os.environ``."""

    return _DEFAULT_OVERLAY.apply(force=force)


# Back-compat alias
def load_site_local_env(*, force: bool = False) -> Path | None:
    return load_dotenv(force=force)


def grteclyn_root(layout: SiteLayout | None = None) -> Path:
    return get_site_layout(layout).grteclyn_root()


def sim_root(layout: SiteLayout | None = None) -> Path:
    return get_site_layout(layout).sim_root()


def grtresna_root(layout: SiteLayout | None = None) -> Path:
    return get_site_layout(layout).grtresna_root()


def chombo_home(layout: SiteLayout | None = None) -> Path:
    return get_site_layout(layout).chombo_home()


def openmpi_root(layout: SiteLayout | None = None) -> Path:
    return get_site_layout(layout).openmpi_root()


def grtresna_env(layout: SiteLayout | None = None) -> Path | None:
    return get_site_layout(layout).grtresna_env()


def grtresna_env_bin(layout: SiteLayout | None = None) -> Path | None:
    active = get_site_layout(layout)
    if isinstance(active, DefaultSiteLayout):
        return active.grtresna_env_bin()
    env = active.grtresna_env()
    if env is None:
        return None
    bin_dir = env / "bin"
    return bin_dir if bin_dir.is_dir() else None


def grtresna_env_lib(layout: SiteLayout | None = None) -> Path | None:
    active = get_site_layout(layout)
    if isinstance(active, DefaultSiteLayout):
        return active.grtresna_env_lib()
    env = active.grtresna_env()
    if env is None:
        return None
    lib_dir = env / "lib"
    return lib_dir if lib_dir.is_dir() else None
