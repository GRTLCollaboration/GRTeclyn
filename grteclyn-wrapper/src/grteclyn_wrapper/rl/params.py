from __future__ import annotations

from pathlib import Path


def _strip_param_value(value: str) -> str:
    value = value.split("#", 1)[0].strip()
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return value


def read_int_param(params_path: Path, key: str, default: int) -> int:
    try:
        for line in params_path.read_text(encoding="utf-8").splitlines():
            if "=" not in line:
                continue
            lhs, rhs = line.split("=", 1)
            if lhs.strip() != key:
                continue
            raw = _strip_param_value(rhs).split()
            if raw:
                return int(float(raw[0]))
    except (FileNotFoundError, ValueError, OSError):
        pass
    return default


def read_rl_num_lumps(params_path: Path, *, default: int = 1) -> int:
    """Read ``rl_num_lumps`` from a GRTeclyn ``params.txt`` (minimum 1)."""
    return max(1, read_int_param(params_path, "rl_num_lumps", default))
