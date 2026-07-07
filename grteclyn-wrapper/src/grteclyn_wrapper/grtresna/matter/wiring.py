"""Map a GRTresna solve configuration to GRTeclyn RadialRecipe matter params."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from typing import TYPE_CHECKING

from ..fields.lump import MAX_INDEPENDENT_LUMPS, lump_sign
from .models import (
    GRTRESNA_BICOMPLEX_SCALAR_MODEL,
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    GRTRESNA_INDEPENDENT_MATTER_MODEL,
    is_bicomplex_scalar_model,
    is_complex_scalar_model,
)

if TYPE_CHECKING:
    from ..solver import GRTresnaConfig

# Re-export for backward compatibility.
__all__ = [
    "EVOLUTION_MATTER_KEYS",
    "GRTRESNA_COMPLEX_SCALAR_MODEL",
    "GRTRESNA_INDEPENDENT_MATTER_MODEL",
    "GRTresnaMatterMetadata",
    "evolution_overrides_from_complex_scalar",
    "evolution_overrides_from_config",
    "evolution_overrides_from_metadata",
    "merge_evolution_overrides",
    "plot_vars_for_complex_scalar",
    "plot_vars_for_bicomplex_scalar",
    "plot_vars_for_independent_scalars",
    "read_matter_metadata",
    "write_matter_metadata",
]

# Canonical GRTeclyn recipe keys required for a valid gridinit replay.
# Single source of truth — replay_eval and promotion scripts import this tuple.
EVOLUTION_MATTER_KEYS: tuple[str, ...] = (
    "recipe_matter_model",
    "recipe_num_scalar_fields",
    "recipe_scalar_field_signs",
    "recipe_scalar_mass",
    "recipe_scalar_lambda",
    "recipe_scalar_mu",
    "recipe_scalar_sign",
    "recipe_exotic_matter",
    "recipe_support_strength",
    "amr.plot_vars",
    "trajectory_pump_frequency",
    "rl_pump_kp",
    "rl_pump_kd",
    "rl_pump_target_profile",
    "rl_pump_target_width",
)

# Full metric in plotfiles for evolved FTL/EC scoring plus matter channels.
_BASE_PLOT_VARS = (
    "chi h11 h12 h13 h22 h23 h33 K lapse shift1 shift2 shift3 phi Pi"
)


def plot_vars_for_complex_scalar() -> tuple[str, ...]:
    """amr.plot_vars for canonical complex scalar (boson star) evolution.

    GRTeclyn stores the imaginary components in the first lump slots
    (``phi_lump0``, ``Pi_lump0``) — see ``StateVariables.hpp`` — not ``phi2``/``Pi2``.
    """
    base = tuple(_BASE_PLOT_VARS.split())
    return (*base, "phi_lump0", "Pi_lump0")


def plot_vars_for_bicomplex_scalar() -> tuple[str, ...]:
    """amr.plot_vars for the two-complex-field (canonical + phantom) model.

    The bicomplex plotfile lays out 8 matter components (see ``StateVariables.hpp``):
    ``phi``/``Pi`` (canonical Phi+ Re/Pi), ``phi_lump0``/``Pi_lump0`` (canonical
    Phi+ Im), ``phi_lump1``/``Pi_lump1`` (phantom Phi- Re/Pi), and
    ``phi_lump2``/``Pi_lump2`` (phantom Phi- Im).  All four lump slots are emitted
    so both the canonical and phantom fields can be visualised.
    """
    base = tuple(_BASE_PLOT_VARS.split())
    return (
        *base,
        "phi_lump0", "Pi_lump0",
        "phi_lump1", "Pi_lump1",
        "phi_lump2", "Pi_lump2",
    )


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
    bs_omega: float = 0.0,
    mu: float = 0.0,
    selfgrav: bool = False,
    phi_c: float = 0.0,
) -> dict[str, Any]:
    """GRTeclyn params for grtresna_complex_scalar matter model.

    When ``bs_omega`` is nonzero the trajectory pump frequency is set so that
    the pump drives Pi and Pi2 coherently at the boson star's internal
    oscillation frequency (U(1) charge injection).

    ``mu`` is the sextic self-interaction coupling (V = 1/2 m^2|Phi|^2 -
    lambda/4 |Phi|^4 + mu/6 |Phi|^6); with lambda>0, mu>0 the lump is a stable
    Q-ball that binds itself rather than dispersing.

    ``selfgrav`` (with the seed central amplitude ``phi_c``) selects the
    self-gravitating boson star: the pump transport target is fitted to the ODE
    star's shape rather than a tail-matched sech.
    """
    overrides: dict[str, Any] = {
        "recipe_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
        "recipe_scalar_mass": float(mass),
        "recipe_scalar_lambda": float(lam),
        "recipe_scalar_mu": float(mu),
        "calculate_constraint_norms": 1,
        "amr.plot_vars": plot_vars_for_complex_scalar(),
    }
    if sign != 1.0:
        overrides["recipe_scalar_sign"] = float(sign)
    overrides.update(
        _pump_controller_overrides(
            bs_omega,
            mass=float(mass),
            selfgrav=selfgrav,
            lam=float(lam),
            mu=float(mu),
            phi_c=float(phi_c),
        )
    )
    return overrides


# Closed-loop PD "trap" controller gains for the matter pump.  Stiff tracking:
# the pump drives each lump toward its moving target soliton (error-proportional,
# self-limiting) so the soliton is confined/transported along the trajectory
# instead of being dispersed by open-loop forcing.  kd ~ 2*sqrt(kp) => ~critical
# damping; both well within the explicit-stepping stability limit.
_PUMP_CONTROLLER_KP = 12.0
_PUMP_CONTROLLER_KD = 7.0


def _pump_controller_overrides(
    bs_omega: float,
    mass: float = 0.0,
    *,
    selfgrav: bool = False,
    lam: float = 0.0,
    mu: float = 0.0,
    phi_c: float = 0.0,
) -> dict[str, Any]:
    """Evolution params enabling the closed-loop PD trap pump controller.

    ``trajectory_pump_frequency`` carries the boson-star phase velocity (bs_omega)
    so the controller's target Pi* rotates coherently with the field's own U(1)
    phase.  ``rl_pump_kp/kd`` switch the pump from the legacy open-loop source to
    the feedback controller (kp > 0).  The trap target is the BOUND sech lump
    (``rl_pump_target_profile=2``); its width must match the initial-data profile
    so the controller drives the field toward the state it was seeded in.

      * Sech / Q-ball seeds: ``rl_pump_target_width = 1/sqrt(m^2-omega^2)`` (the
        bound-state tail length the sech seed was built with).
      * Self-gravitating star: the seed is the ODE profile, whose core is more
        compact than a tail-matched sech.  We fit the sech width to that profile
        (half-maximum match) so the transport target has the star's actual shape
        instead of a broader, off-equilibrium sech.
    """
    from ..profiles.boson_star import PROFILE_SECH_BOUND, bound_width

    overrides: dict[str, Any] = {
        "trajectory_pump_frequency": float(bs_omega),
        "rl_pump_kp": _PUMP_CONTROLLER_KP,
        "rl_pump_kd": _PUMP_CONTROLLER_KD,
    }
    if mass > 0.0:
        overrides["rl_pump_target_profile"] = int(PROFILE_SECH_BOUND)
        if selfgrav:
            from ..profiles.boson_star_ode import (
                default_selfgrav_phi_c,
                selfgrav_pump_target_width,
            )

            target_width = selfgrav_pump_target_width(mass, lam, mu, phi_c)
            # The trap must aim at the star's OWN central amplitude, not the
            # transport knob well_depth.  Aiming at a much smaller amplitude
            # makes the PD controller fight gravity and collapse the star.
            target_amp = phi_c if phi_c > 0.0 else default_selfgrav_phi_c(lam, mu)
            overrides["rl_pump_target_amp"] = float(target_amp)
        else:
            target_width = bound_width(mass, bs_omega)
        overrides["rl_pump_target_width"] = float(target_width)
    return overrides


def evolution_overrides_from_bicomplex_scalar(
    mass: float,
    lam: float,
    field_signs: tuple[int, ...],
    bs_omega: float = 0.0,
    mu: float = 0.0,
    selfgrav: bool = False,
    phi_c: float = 0.0,
) -> dict[str, Any]:
    """GRTeclyn params for the two-complex-field (canonical + phantom) model.

    ``field_signs`` is the per-lump sign list (+1 canonical, -1 exotic), in the
    trajectory lump order.  It routes each pump spotlight to the canonical
    (Phi+) or phantom (Phi-) field via ``recipe_scalar_field_signs``.

    ``bs_omega`` drives the closed-loop PD trap pump controller (target Pi*
    rotates at the boson phase velocity).
    """
    signs = field_signs[:MAX_INDEPENDENT_LUMPS]
    overrides = {
        "recipe_matter_model": GRTRESNA_BICOMPLEX_SCALAR_MODEL,
        "recipe_scalar_mass": float(mass),
        "recipe_scalar_lambda": float(lam),
        "recipe_scalar_mu": float(mu),
        "recipe_num_scalar_fields": len(signs),
        "recipe_scalar_field_signs": " ".join(str(int(s)) for s in signs),
        "calculate_constraint_norms": 1,
        "amr.plot_vars": plot_vars_for_bicomplex_scalar(),
    }
    overrides.update(
        _pump_controller_overrides(
            bs_omega,
            mass=float(mass),
            selfgrav=selfgrav,
            lam=float(lam),
            mu=float(mu),
            phi_c=float(phi_c),
        )
    )
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
    scalar_mu: float = 0.0
    bs_phi_c: float = 0.0
    bs_profile_width: float = 0.0
    bs_omega: float = 0.0
    scalar_sign: int = 1
    bs_selfgrav: bool = False
    # Lump centres as offsets from the throat (== RadialRecipe centred-frame
    # positions), used to seed the RL lump tracker / pump spotlights.
    lump_centers: tuple[tuple[float, float, float], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_config(cls, cfg: GRTresnaConfig) -> GRTresnaMatterMetadata:  # noqa: F821
        selfgrav = _has_selfgrav_lump(getattr(cfg, "lumps", ()))
        if is_bicomplex_scalar_model(getattr(cfg, "matter_model", "")):
            lumps = list(getattr(cfg, "lumps", ()) or ())
            signs = tuple(lump_sign(lump) for lump in lumps[:MAX_INDEPENDENT_LUMPS])
            return cls(
                matter_model=GRTRESNA_BICOMPLEX_SCALAR_MODEL,
                num_scalar_fields=len(signs),
                scalar_field_signs=signs,
                scalar_mass=float(cfg.scalar_mass),
                scalar_lambda=float(cfg.scalar_lambda),
                scalar_mu=float(getattr(cfg, "scalar_mu", 0.0)),
                lump_count=len(lumps),
                bs_phi_c=float(getattr(cfg, "bs_phi_c", 0.0)),
                bs_profile_width=float(getattr(cfg, "bs_profile_width", 0.0)),
                bs_omega=float(getattr(cfg, "bs_omega", 0.0)),
                scalar_sign=int(getattr(cfg, "scalar_sign", 1)),
                bs_selfgrav=selfgrav,
                lump_centers=_lump_centers(lumps),
            )

        if is_complex_scalar_model(getattr(cfg, "matter_model", "")):
            return cls(
                matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
                num_scalar_fields=0,
                scalar_field_signs=(),
                scalar_mass=float(cfg.scalar_mass),
                scalar_lambda=float(cfg.scalar_lambda),
                scalar_mu=float(getattr(cfg, "scalar_mu", 0.0)),
                lump_count=0,
                bs_phi_c=float(cfg.bs_phi_c),
                bs_profile_width=float(cfg.bs_profile_width),
                bs_omega=float(cfg.bs_omega),
                scalar_sign=int(getattr(cfg, "scalar_sign", 1)),
                bs_selfgrav=selfgrav,
                lump_centers=_lump_centers(getattr(cfg, "lumps", ())),
            )

        lumps = list(cfg.lumps)
        signs = tuple(lump_sign(lump) for lump in lumps[:MAX_INDEPENDENT_LUMPS])
        return cls(
            matter_model=GRTRESNA_INDEPENDENT_MATTER_MODEL,
            num_scalar_fields=len(signs),
            scalar_field_signs=signs,
            scalar_mass=float(cfg.scalar_mass),
            scalar_lambda=float(cfg.scalar_lambda),
            scalar_mu=float(getattr(cfg, "scalar_mu", 0.0)),
            lump_count=len(lumps),
            lump_centers=_lump_centers(lumps),
        )

    def to_evolution_overrides(self) -> dict[str, Any]:
        """GRTeclyn recipe params for GPU replay from serialized matter metadata."""
        if is_bicomplex_scalar_model(self.matter_model):
            return evolution_overrides_from_bicomplex_scalar(
                mass=self.scalar_mass,
                lam=self.scalar_lambda,
                field_signs=self.scalar_field_signs,
                bs_omega=self.bs_omega,
                mu=self.scalar_mu,
                selfgrav=self.bs_selfgrav,
                phi_c=self.bs_phi_c,
            )
        if is_complex_scalar_model(self.matter_model):
            return evolution_overrides_from_complex_scalar(
                mass=self.scalar_mass,
                lam=self.scalar_lambda,
                sign=float(self.scalar_sign),
                bs_omega=self.bs_omega,
                mu=self.scalar_mu,
                selfgrav=self.bs_selfgrav,
                phi_c=self.bs_phi_c,
            )
        if self.num_scalar_fields == 0:
            return {"calculate_constraint_norms": 1}
        return {
            "recipe_matter_model": self.matter_model,
            "recipe_num_scalar_fields": self.num_scalar_fields,
            "recipe_scalar_field_signs": " ".join(str(s) for s in self.scalar_field_signs),
            "recipe_scalar_mass": self.scalar_mass,
            "recipe_scalar_lambda": self.scalar_lambda,
            "calculate_constraint_norms": 1,
            "amr.plot_vars": plot_vars_for_independent_scalars(self.num_scalar_fields),
        }


def evolution_overrides_from_metadata(meta: GRTresnaMatterMetadata) -> dict[str, Any]:
    """Rebuild evolution overrides from ``initial_data.matter.json``."""
    return meta.to_evolution_overrides()


def _has_selfgrav_lump(lumps: Any) -> bool:
    """True when any lump uses the self-gravitating boson-star seed (profile 4)."""
    from ..profiles.boson_star import PROFILE_SELFGRAV_BOUND

    for lump in lumps or ():
        try:
            if int(lump.get("profile", 0)) == PROFILE_SELFGRAV_BOUND:
                return True
        except (AttributeError, TypeError, ValueError):
            continue
    return False


def evolution_overrides_from_config(cfg: GRTresnaConfig) -> dict[str, Any]:  # noqa: F821
    """GRTeclyn params that select the matched matter model."""
    selfgrav = _has_selfgrav_lump(getattr(cfg, "lumps", ()))
    if is_bicomplex_scalar_model(getattr(cfg, "matter_model", "")):
        lumps = list(getattr(cfg, "lumps", ()) or ())
        signs = tuple(lump_sign(lump) for lump in lumps[:MAX_INDEPENDENT_LUMPS])
        return evolution_overrides_from_bicomplex_scalar(
            mass=float(cfg.scalar_mass),
            lam=float(cfg.scalar_lambda),
            field_signs=signs,
            bs_omega=float(getattr(cfg, "bs_omega", 0.0)),
            mu=float(getattr(cfg, "scalar_mu", 0.0)),
            selfgrav=selfgrav,
            phi_c=float(getattr(cfg, "bs_phi_c", 0.0)),
        )

    if is_complex_scalar_model(getattr(cfg, "matter_model", "")):
        return evolution_overrides_from_complex_scalar(
            mass=float(cfg.scalar_mass),
            lam=float(cfg.scalar_lambda),
            sign=float(getattr(cfg, "scalar_sign", 1)),
            bs_omega=float(getattr(cfg, "bs_omega", 0.0)),
            mu=float(getattr(cfg, "scalar_mu", 0.0)),
            selfgrav=selfgrav,
            phi_c=float(getattr(cfg, "bs_phi_c", 0.0)),
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
        scalar_mu=float(payload.get("scalar_mu", 0.0)),
        lump_count=int(payload.get("lump_count", len(signs))),
        bs_phi_c=float(payload.get("bs_phi_c", 0.0)),
        bs_profile_width=float(payload.get("bs_profile_width", 0.0)),
        bs_omega=float(payload.get("bs_omega", 0.0)),
        scalar_sign=int(payload.get("scalar_sign", 1)),
        bs_selfgrav=bool(payload.get("bs_selfgrav", False)),
        lump_centers=tuple(
            (float(c[0]), float(c[1]), float(c[2]))
            for c in payload.get("lump_centers", ())
            if c is not None and len(c) >= 3
        ),
    )


def _lump_centers(lumps: Any) -> tuple[tuple[float, float, float], ...]:
    """Extract lump centres (throat-relative offsets) from lump specs."""
    out: list[tuple[float, float, float]] = []
    for lump in lumps or ():
        try:
            c = lump.get("center", (0.0, 0.0, 0.0))
            out.append((float(c[0]), float(c[1]), float(c[2])))
        except (AttributeError, TypeError, IndexError, ValueError):
            continue
    return tuple(out)


def merge_evolution_overrides(
    base: Mapping[str, Any] | None,
    cfg: GRTresnaConfig,  # noqa: F821
) -> dict[str, Any]:
    merged = dict(base or {})
    merged.update(evolution_overrides_from_config(cfg))
    return merged
