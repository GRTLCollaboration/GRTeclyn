"""GRTresna matter-model registry and search → solver projection.

Single source of truth for which Chombo example and GRTeclyn matter path
correspond to a given initial-data configuration.  Keeps the constraint-solve
T_ab identical to the evolution T_ab (MapElites hard consistency rule).

Matter is selected independently of geometry ansatz (shell / ring / free):

  sector   — ``scalar`` (independent real lumps) or ``boson_star`` (complex U(1))
  coupling — ``canonical`` (positive energy) or ``exotic`` (phantom / −rho)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

GRTRESNA_INDEPENDENT_MATTER_MODEL = "grtresna_independent_scalars"
GRTRESNA_COMPLEX_SCALAR_MODEL = "grtresna_complex_scalar"

GRTRESNA_EXAMPLE_SCALAR_FIELD_BH = "ScalarFieldBH"
GRTRESNA_EXAMPLE_BOSON_STAR_BH = "BosonStarBH"

MATTER_SECTOR_SCALAR = "scalar"
MATTER_SECTOR_BOSON_STAR = "boson_star"

MATTER_COUPLING_CANONICAL = "canonical"
MATTER_COUPLING_EXOTIC = "exotic"

_SECTOR_ALIASES: dict[str, str] = {
    "scalar": MATTER_SECTOR_SCALAR,
    "scalars": MATTER_SECTOR_SCALAR,
    "independent_scalars": MATTER_SECTOR_SCALAR,
    "grtresna_independent_scalars": MATTER_SECTOR_SCALAR,
    "boson_star": MATTER_SECTOR_BOSON_STAR,
    "boson": MATTER_SECTOR_BOSON_STAR,
    "complex_scalar": MATTER_SECTOR_BOSON_STAR,
    "grtresna_complex_scalar": MATTER_SECTOR_BOSON_STAR,
}

_COUPLING_ALIASES: dict[str, str] = {
    "canonical": MATTER_COUPLING_CANONICAL,
    "real": MATTER_COUPLING_CANONICAL,
    "positive": MATTER_COUPLING_CANONICAL,
    "exotic": MATTER_COUPLING_EXOTIC,
    "phantom": MATTER_COUPLING_EXOTIC,
    "negative": MATTER_COUPLING_EXOTIC,
}


@dataclass(frozen=True)
class MatterModelSpec:
    """Immutable descriptor for one supported matter sector."""

    model_id: str
    example: str


MATTER_MODELS: dict[str, MatterModelSpec] = {
    GRTRESNA_INDEPENDENT_MATTER_MODEL: MatterModelSpec(
        GRTRESNA_INDEPENDENT_MATTER_MODEL,
        GRTRESNA_EXAMPLE_SCALAR_FIELD_BH,
    ),
    GRTRESNA_COMPLEX_SCALAR_MODEL: MatterModelSpec(
        GRTRESNA_COMPLEX_SCALAR_MODEL,
        GRTRESNA_EXAMPLE_BOSON_STAR_BH,
    ),
}


@dataclass(frozen=True)
class MatterSelection:
    """User-facing matter choice: sector × coupling."""

    sector: str
    coupling: str

    @property
    def model_id(self) -> str:
        if self.sector == MATTER_SECTOR_BOSON_STAR:
            return GRTRESNA_COMPLEX_SCALAR_MODEL
        return GRTRESNA_INDEPENDENT_MATTER_MODEL

    @property
    def is_boson_star(self) -> bool:
        return self.sector == MATTER_SECTOR_BOSON_STAR

    @property
    def is_scalar(self) -> bool:
        return self.sector == MATTER_SECTOR_SCALAR

    @property
    def is_exotic(self) -> bool:
        return self.coupling == MATTER_COUPLING_EXOTIC

    @property
    def scalar_sign(self) -> int:
        return -1 if self.is_exotic else 1


def normalize_matter_sector(raw: str) -> str:
    key = str(raw).strip().lower()
    try:
        return _SECTOR_ALIASES[key]
    except KeyError as exc:
        raise ValueError(
            f"unknown matter sector: {raw!r} "
            f"(choices: {MATTER_SECTOR_SCALAR}, {MATTER_SECTOR_BOSON_STAR})"
        ) from exc


def normalize_matter_coupling(raw: str) -> str:
    key = str(raw).strip().lower()
    try:
        return _COUPLING_ALIASES[key]
    except KeyError as exc:
        raise ValueError(
            f"unknown matter coupling: {raw!r} "
            f"(choices: {MATTER_COUPLING_CANONICAL}, {MATTER_COUPLING_EXOTIC})"
        ) from exc


def resolve_matter_selection(sector: str, coupling: str) -> MatterSelection:
    """Build a validated matter selection from CLI / campaign knobs."""
    return MatterSelection(
        sector=normalize_matter_sector(sector),
        coupling=normalize_matter_coupling(coupling),
    )


def matter_selection_base_overrides(selection: MatterSelection) -> dict[str, Any]:
    """Base overrides stamped on every eval for replay / HQ handoff."""
    overrides: dict[str, Any] = {
        "grtresna_matter_model": selection.model_id,
        "grtresna_matter_sector": selection.sector,
        "grtresna_matter_coupling": selection.coupling,
    }
    if selection.is_boson_star and selection.is_exotic:
        overrides["grtresna_scalar_sign"] = -1.0
    if selection.is_scalar and selection.is_exotic:
        overrides["grtresna_force_exotic_lumps"] = 1.0
    return overrides


def matter_model_for_ansatz(ansatz: str) -> str:
    """Backward-compatible helper: ``boson_star`` ansatz → complex scalar model.

    Prefer :func:`resolve_matter_selection` with ``--grtresna-matter-sector``.
    """
    if str(ansatz).strip().lower() == MATTER_SECTOR_BOSON_STAR:
        return GRTRESNA_COMPLEX_SCALAR_MODEL
    return GRTRESNA_INDEPENDENT_MATTER_MODEL


def is_complex_scalar_model(matter_model: str) -> bool:
    return matter_model == GRTRESNA_COMPLEX_SCALAR_MODEL


def is_boson_star_sector(sector: str) -> bool:
    return normalize_matter_sector(sector) == MATTER_SECTOR_BOSON_STAR


def resolve_grtresna_example(
    matter_model: str,
    *,
    has_lumps: bool = False,
) -> str:
    """Pick the GRTresna binary example directory for a matter configuration."""
    if is_complex_scalar_model(matter_model):
        return GRTRESNA_EXAMPLE_BOSON_STAR_BH
    if has_lumps or matter_model == GRTRESNA_INDEPENDENT_MATTER_MODEL:
        return GRTRESNA_EXAMPLE_SCALAR_FIELD_BH
    return GRTRESNA_EXAMPLE_SCALAR_FIELD_BH


def finalize_solver_config(cfg: Any) -> None:
    """Align ``cfg.example`` with ``cfg.matter_model`` before launching GRTresna."""
    matter_model = getattr(cfg, "matter_model", "") or ""
    if not matter_model and getattr(cfg, "lumps", None):
        matter_model = GRTRESNA_INDEPENDENT_MATTER_MODEL
    cfg.example = resolve_grtresna_example(
        matter_model,
        has_lumps=bool(getattr(cfg, "lumps", None)),
    )


def _get_float(
    overrides: Mapping[str, Any],
    key: str,
    default: float,
) -> float:
    try:
        return float(overrides.get(key, default))
    except (TypeError, ValueError):
        return default


def _coupling_from_overrides(overrides: Mapping[str, Any]) -> str:
    return normalize_matter_coupling(
        str(overrides.get("grtresna_matter_coupling", MATTER_COUPLING_CANONICAL))
    )


def apply_scalar_exotic_coupling(
    cfg: Any,
    overrides: Mapping[str, Any],
    *,
    enable_exotic_safe_solver: Callable[[], None] | None = None,
) -> None:
    """Force every independent-scalar lump exotic when coupling is exotic."""
    force = int(round(float(overrides.get("grtresna_force_exotic_lumps", 0)))) >= 1
    if _coupling_from_overrides(overrides) != MATTER_COUPLING_EXOTIC and not force:
        return
    if getattr(cfg, "lumps", None):
        for lump in cfg.lumps:
            lump["exotic"] = 1
    elif float(getattr(cfg, "lump_amp", 0.0)) > 0.0 or "grtresna_lump_amp" in overrides:
        cfg.lump_exotic = 1
    if enable_exotic_safe_solver is not None and (
        any(int(l.get("exotic", 0)) for l in getattr(cfg, "lumps", []) or [])
        or int(getattr(cfg, "lump_exotic", 0))
    ):
        enable_exotic_safe_solver()


def finalize_independent_scalar_config(
    cfg: Any,
    overrides: Mapping[str, Any],
    *,
    enable_exotic_safe_solver: Callable[[], None] | None = None,
) -> Any:
    """Tag cfg as independent scalars and apply campaign-level exotic coupling."""
    cfg.matter_model = GRTRESNA_INDEPENDENT_MATTER_MODEL
    apply_scalar_exotic_coupling(
        cfg,
        overrides,
        enable_exotic_safe_solver=enable_exotic_safe_solver,
    )
    return cfg


def apply_boson_star_overrides(
    cfg: Any,
    overrides: Mapping[str, Any],
    *,
    enable_exotic_safe_solver: Callable[[], None] | None = None,
) -> bool:
    """Configure ``cfg`` for the complex-scalar (mini boson star) sector.

    Returns True when boson-star overrides were applied (search path active).
    """
    matter_model = str(overrides.get("grtresna_matter_model", "")).strip()
    sector_raw = str(overrides.get("grtresna_matter_sector", "")).strip().lower()
    has_bs_keys = any(str(k).startswith("grtresna_bs_") for k in overrides)
    is_boson = (
        is_complex_scalar_model(matter_model)
        or matter_model in ("boson_star", "complex_scalar")
        or has_bs_keys
    )
    if not is_boson and sector_raw:
        try:
            is_boson = normalize_matter_sector(sector_raw) == MATTER_SECTOR_BOSON_STAR
        except ValueError:
            is_boson = False
    if not is_boson:
        return False

    cfg.matter_model = GRTRESNA_COMPLEX_SCALAR_MODEL
    cfg.example = GRTRESNA_EXAMPLE_BOSON_STAR_BH
    cfg.lumps = []

    cfg.scalar_mass = _get_float(overrides, "grtresna_scalar_mass", cfg.scalar_mass)
    cfg.scalar_lambda = _get_float(overrides, "grtresna_scalar_lambda", cfg.scalar_lambda)
    cfg.bs_phi_c = _get_float(overrides, "grtresna_bs_phi_c", cfg.bs_phi_c)
    cfg.bs_profile_width = _get_float(
        overrides, "grtresna_bs_profile_width", cfg.bs_profile_width
    )
    cfg.bs_omega = _get_float(overrides, "grtresna_bs_omega", cfg.bs_omega)
    cfg.shift_seed = _get_float(overrides, "grtresna_shift_seed", cfg.shift_seed)

    if _coupling_from_overrides(overrides) == MATTER_COUPLING_EXOTIC:
        cfg.scalar_sign = -1
    elif "grtresna_scalar_sign" in overrides:
        cfg.scalar_sign = int(round(float(overrides["grtresna_scalar_sign"])))
    else:
        cfg.scalar_sign = int(getattr(cfg, "scalar_sign", 1) or 1)

    if cfg.scalar_sign < 0 and enable_exotic_safe_solver is not None:
        enable_exotic_safe_solver()

    return True
