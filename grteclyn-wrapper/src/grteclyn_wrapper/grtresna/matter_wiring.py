"""Map a GRTresna solve configuration to GRTeclyn RadialRecipe matter params."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from typing import TYPE_CHECKING

from .lump_fields import MAX_INDEPENDENT_LUMPS, lump_sign

if TYPE_CHECKING:
    from .solver import GRTresnaConfig

GRTRESNA_INDEPENDENT_MATTER_MODEL = "grtresna_independent_scalars"
GRTRESNA_COMPLEX_SCALAR_MODEL = "grtresna_complex_scalar"

# Full metric in plotfiles for evolved FTL/EC scoring plus per-lump scalars.
_BASE_PLOT_VARS = (
    "chi h11 h12 h13 h22 h23 h33 K lapse shift1 shift2 shift3 phi Pi"
)


def plot_vars_for_complex_scalar() -> tuple[str, ...]:
    """amr.plot_vars for canonical complex scalar (boson star) evolution."""
    base = tuple(_BASE_PLOT_VARS.split())
    return (*base, "phi2", "Pi2")


def plot_vars_for_independent_scalars(num_fields: int) -> tuple[str, ...]:
    """amr.plot_vars including per-lump scalar channels for frame extraction."""
    base = tuple(_BASE_PLOT_VARS.split())
    n = max(0, min(int(num_fields), MAX_INDEPENDENT_LUMPS))
    lump_names: list[str] = []
    for k in range(n):
        lump_names.extend([f"phi_lump{k}", f"Pi_lump{k}"])
    return (*base, *lump_names) if lump_names else base


def evolution_overrides_from_complex_scalar(
    mass: float,
    lam: float = 0.0,
    sign: float = 1.0,
) -> dict[str, Any]:
    """GRTeclyn params for grtresna_complex_scalar matter model."""
    overrides: dict[str, Any] = {
        "recipe_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
        "recipe_scalar_mass": float(mass),
        "recipe_scalar_lambda": float(lam),
        "calculate_constraint_norms": 1,
        "amr.plot_vars": plot_vars_for_complex_scalar(),
    }
    if sign != 1.0:
        overrides["recipe_scalar_sign"] = float(sign)
    return overrides


@dataclass(frozen=True)
class GRTresnaMatterMetadata:
    """Serialized matter layout for a GRTresna-projected episode."""

    matter_model: str
    num_scalar_fields: int
    scalar_field_signs: tuple[int, ...]
    scalar_mass: float
    scalar_lambda: float
    lump_count: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_config(cls, cfg: GRTresnaConfig) -> GRTresnaMatterMetadata:  # noqa: F821
        lumps = list(cfg.lumps)
        signs = tuple(lump_sign(lump) for lump in lumps[:MAX_INDEPENDENT_LUMPS])
        return cls(
            matter_model=GRTRESNA_INDEPENDENT_MATTER_MODEL,
            num_scalar_fields=len(signs),
            scalar_field_signs=signs,
            scalar_mass=float(cfg.scalar_mass),
            scalar_lambda=float(cfg.scalar_lambda),
            lump_count=len(lumps),
        )


def evolution_overrides_from_config(cfg: GRTresnaConfig) -> dict[str, Any]:  # noqa: F821
    """GRTeclyn params that select the matched matter model."""
    if getattr(cfg, "matter_model", "") == GRTRESNA_COMPLEX_SCALAR_MODEL:
        return evolution_overrides_from_complex_scalar(
            mass=float(cfg.scalar_mass),
            lam=float(cfg.scalar_lambda),
            sign=float(getattr(cfg, "scalar_sign", 1)),
        )

    meta = GRTresnaMatterMetadata.from_config(cfg)
    if meta.num_scalar_fields == 0:
        return {"calculate_constraint_norms": 1}

    return {
        "recipe_matter_model": meta.matter_model,
        "recipe_num_scalar_fields": meta.num_scalar_fields,
        "recipe_scalar_field_signs": " ".join(str(s) for s in meta.scalar_field_signs),
        "recipe_scalar_mass": meta.scalar_mass,
        "recipe_scalar_lambda": meta.scalar_lambda,
        "calculate_constraint_norms": 1,
        "amr.plot_vars": plot_vars_for_independent_scalars(meta.num_scalar_fields),
    }


def write_matter_metadata(path: str | Path, cfg: GRTresnaConfig) -> Path:  # noqa: F821
    path = Path(path)
    payload = GRTresnaMatterMetadata.from_config(cfg).to_dict()
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_matter_metadata(path: str | Path) -> GRTresnaMatterMetadata:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    signs = tuple(int(v) for v in payload.get("scalar_field_signs", ()))
    return GRTresnaMatterMetadata(
        matter_model=str(payload.get("matter_model", GRTRESNA_INDEPENDENT_MATTER_MODEL)),
        num_scalar_fields=int(payload.get("num_scalar_fields", len(signs))),
        scalar_field_signs=signs,
        scalar_mass=float(payload.get("scalar_mass", 0.0)),
        scalar_lambda=float(payload.get("scalar_lambda", 0.0)),
        lump_count=int(payload.get("lump_count", len(signs))),
    )


def merge_evolution_overrides(
    base: Mapping[str, Any] | None,
    cfg: GRTresnaConfig,  # noqa: F821
) -> dict[str, Any]:
    merged = dict(base or {})
    merged.update(evolution_overrides_from_config(cfg))
    return merged
