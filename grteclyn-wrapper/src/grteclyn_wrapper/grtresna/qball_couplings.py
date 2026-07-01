"""Q-ball self-interaction couplings — single source of truth for search/replay."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from .boson_star_profile import bound_width


@dataclass(frozen=True)
class QBallCouplings:
    """Flat-space Q-ball potential couplings.

    Potential: V = ½ m² |Φ|² − ¼ λ |Φ|⁴ + ⅙ μ |Φ|⁶

    Design relations (flat space):
      ρ* = 3λ / 4μ,  |Φ|_core ≈ √ρ*
      ω_min² = m² − 3λ² / (16μ)
    """

    mass: float
    lam: float
    mu: float
    omega: float
    phi_core: float = 0.15
    """Target lump paint amplitude (well_depth cap), not necessarily √ρ*."""

    @classmethod
    def free_field(cls, *, mass: float = 1.0, omega: float = 0.4) -> QBallCouplings:
        """Massive Klein–Gordon field with no self-interaction."""
        return cls(mass=mass, lam=0.0, mu=0.0, omega=omega, phi_core=0.0)

    @classmethod
    def standard(cls) -> QBallCouplings:
        """Stage-A test couplings (λ=160, μ=5333) used in eval-122 Q-ball runs."""
        return cls(mass=1.0, lam=160.0, mu=5333.0, omega=0.4)

    @classmethod
    def stiff(cls) -> QBallCouplings:
        """NextSteps item 2: 4× λ, μ chosen to preserve ω_min (λ²/μ fixed).

        Equilibrium core √(3λ/4μ) is ~0.075 (half of ``standard()``).
        Use ``with_equilibrium_paint()`` when seeding to avoid super-equilibrium radiation.
        """
        return cls(mass=1.0, lam=640.0, mu=85333.0, omega=0.4)

    @classmethod
    def compact(cls) -> QBallCouplings:
        """Thick-wall couplings for a *localized* seed that fits the evolution box.

        ``standard()``/``stiff()`` use ω=0.4, deep in the thin-wall regime where the
        equilibrium Q-ball is a flat-top condensate of radius ~20 code units — far
        larger than a multi-lump box, and (before the radial-ODE shoot was fixed)
        the source of the non-localized seed documented in NextSteps.md §3.  Raising
        ω into the thick-wall regime (ω≈0.8·m) collapses the flat top so the soliton
        decays to ~0 within ~7 code units while keeping the sextic stabiliser and
        the same λ, μ.  Use this for ``--qball-ode-profile`` seeding.
        """
        return cls(mass=1.0, lam=640.0, mu=85333.0, omega=0.8)

    def with_equilibrium_paint(self) -> QBallCouplings:
        """Return a copy whose ``phi_core`` matches flat-space √(3λ/4μ)."""
        return QBallCouplings(
            mass=self.mass,
            lam=self.lam,
            mu=self.mu,
            omega=self.omega,
            phi_core=self.core_amplitude if self.lam > 0.0 else self.phi_core,
        )

    def cap_well_depth(self, well_depth: float, *, legacy_cap: float = 0.15) -> float:
        """Cap trajectory well_depth for on-attractor seeding."""
        depth = min(legacy_cap, max(0.0, float(well_depth)))
        if self.lam > 0.0 and self.mu > 0.0:
            return min(depth, self.core_amplitude)
        return depth

    @property
    def omega_min(self) -> float:
        """Lower edge of the Q-ball existence band (ω > ω_min required)."""
        if self.lam <= 0.0 or self.mu <= 0.0:
            return 0.0
        omega_min_sq = self.mass * self.mass - 3.0 * self.lam * self.lam / (16.0 * self.mu)
        return math.sqrt(max(0.0, omega_min_sq))

    @property
    def core_amplitude(self) -> float:
        """Equilibrium core amplitude √ρ* from λ, μ (0 when μ=0)."""
        if self.mu <= 0.0 or self.lam <= 0.0:
            return 0.0
        return math.sqrt(3.0 * self.lam / (4.0 * self.mu))

    @property
    def bound_state_width(self) -> float:
        """Tail decay scale R = 1 / √(m² − ω²)."""
        return bound_width(self.mass, self.omega)

    def scaled_preserving_band(self, lam_factor: float) -> QBallCouplings:
        """Scale λ by *lam_factor* and μ by λ² so ω_min is unchanged."""
        if lam_factor <= 0.0:
            raise ValueError("lam_factor must be positive")
        f = float(lam_factor)
        return QBallCouplings(
            mass=self.mass,
            lam=self.lam * f,
            mu=self.mu * f * f,
            omega=self.omega,
            phi_core=self.phi_core,
        )

    def scaled(self, factor: float) -> QBallCouplings:
        """Scale λ and μ together — preserves equilibrium core amplitude."""
        if factor <= 0.0:
            raise ValueError("scale factor must be positive")
        return QBallCouplings(
            mass=self.mass,
            lam=self.lam * factor,
            mu=self.mu * factor,
            omega=self.omega,
            phi_core=self.phi_core,
        )

    def validate(self) -> None:
        """Fail fast if couplings are inconsistent with a bound Q-ball."""
        if self.mass <= 0.0:
            raise ValueError("mass must be positive")
        if self.lam < 0.0 or self.mu < 0.0:
            raise ValueError("lam and mu must be non-negative")
        if self.lam > 0.0 and self.mu <= 0.0:
            raise ValueError("mu > 0 required when lam > 0 (sextic stabiliser)")
        if self.omega >= self.mass:
            raise ValueError(
                f"omega ({self.omega}) must be < mass ({self.mass}) for confinement"
            )
        if self.lam > 0.0 and self.mu > 0.0:
            if self.omega <= self.omega_min:
                raise ValueError(
                    f"omega ({self.omega}) must exceed omega_min ({self.omega_min:.4f})"
                )

    def seed_overrides(
        self,
        *,
        equilibrium_amplitude: bool = False,
        ode_profile: bool = False,
    ) -> dict[str, float | int]:
        """GRTresna search / replay keys plus optional dispersion-fix flags."""
        out: dict[str, float | int] = dict(self.to_search_overrides())
        if equilibrium_amplitude:
            out["grtresna_qball_equilibrium_amplitude"] = 1
        if ode_profile:
            out["grtresna_qball_ode_profile"] = 1
        return out

    def to_search_overrides(self) -> dict[str, float]:
        """GRTresna search / replay override keys for this Q-ball sector."""
        self.validate()
        return {
            "grtresna_scalar_mass": self.mass,
            "grtresna_scalar_lambda": self.lam,
            "grtresna_scalar_mu": self.mu,
            "grtresna_bs_omega": self.omega,
        }

    def merge_into(self, overrides: dict[str, Any]) -> dict[str, Any]:
        """Return a copy of *overrides* with Q-ball physics keys applied."""
        merged = dict(overrides)
        merged.update(self.to_search_overrides())
        return merged
