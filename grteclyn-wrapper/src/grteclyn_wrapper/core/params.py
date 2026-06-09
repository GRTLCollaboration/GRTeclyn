from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from .config import ExampleConfig, resolve_example


_ASSIGNMENT_RE = re.compile(r"^(?P<prefix>\s*(?P<key>[A-Za-z0-9_.]+)\s*=\s*)(?P<value>.*?)(?P<comment>\s+#.*)?$")

# AMReX ParmParse expects space-separated tokens, not one quoted string.
_PARM_PARSE_TOKEN_LIST_KEYS = frozenset(
    {"amr.plot_vars", "amr.derive_plot_vars"},
)


def format_param_value(key: str, value: object) -> str:
    """Format a params.txt assignment value, honoring ParmParse list keys."""
    if key in _PARM_PARSE_TOKEN_LIST_KEYS:
        if isinstance(value, str):
            return value.strip().strip('"')
        if isinstance(value, (list, tuple)):
            return " ".join(str(item) for item in value)
    return format_value(value)


def format_value(value: object) -> str:
    if isinstance(value, Path):
        return f'"{value.expanduser().resolve()}"'
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.startswith('"') and stripped.endswith('"'):
            return stripped
        # A space-separated list of plain numbers is an array literal (e.g. a
        # boundary spec "1 1 1" or a center "32 32 32") and must NOT be quoted,
        # otherwise AMReX ParmParse fails to parse it.
        tokens = stripped.split()
        if len(tokens) > 1 and "/" not in stripped:
            def _is_num(tok: str) -> bool:
                try:
                    float(tok)
                    return True
                except ValueError:
                    return False
            if all(_is_num(tok) for tok in tokens):
                return stripped
        if any(ch.isspace() for ch in stripped) or "/" in stripped:
            return f'"{stripped}"'
        return stripped
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (list, tuple)):
        return " ".join(format_value(item).strip('"') for item in value)
    return str(value)


@dataclass
class ParamsTemplate:
    lines: list[str]

    @classmethod
    def load(cls, path: Path) -> "ParamsTemplate":
        return cls(path.read_text(encoding="utf-8").splitlines())

    def render(self, overrides: Mapping[str, object]) -> str:
        remaining = {
            str(key): format_param_value(str(key), value) for key, value in overrides.items()
        }
        rendered: list[str] = []

        for line in self.lines:
            match = _ASSIGNMENT_RE.match(line)
            if not match:
                rendered.append(line)
                continue

            key = match.group("key")
            if key not in remaining:
                rendered.append(line)
                continue

            comment = match.group("comment") or ""
            rendered.append(f"{match.group('prefix')}{remaining.pop(key)}{comment}")

        if remaining:
            rendered.append("")
            rendered.append("# Wrapper overrides not present in the source template")
            for key in sorted(remaining):
                rendered.append(f"{key} = {remaining[key]}")

        return "\n".join(rendered) + "\n"


def episode_path_overrides(episode_dir: Path, example: ExampleConfig | str = "SupportedWormholeCollapse") -> dict[str, object]:
    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    episode_dir = episode_dir.expanduser().resolve()
    return {
        "output_path": episode_dir,
        "amr.check_file": episode_dir / example_cfg.check_prefix,
        "amr.plot_file": episode_dir / example_cfg.plot_prefix,
    }


def write_params(
    template_path: Path,
    output_path: Path,
    *,
    episode_dir: Path,
    example: ExampleConfig | str = "SupportedWormholeCollapse",
    overrides: Mapping[str, object] | None = None,
) -> Path:
    merged: dict[str, object] = episode_path_overrides(episode_dir, example=example)
    if overrides:
        merged.update(overrides)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(ParamsTemplate.load(template_path).render(merged), encoding="utf-8")
    return output_path
