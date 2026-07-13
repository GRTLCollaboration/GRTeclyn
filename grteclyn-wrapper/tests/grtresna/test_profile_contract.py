"""Matter-profile CONTRACT rail: cross-code t=0 consistency test.

These tests lock down the invariants that, when silently violated, produce a
"plausible-looking wrong answer" (dispersal / t=0 constraint blow-up) rather than
an error.  They demonstrate the rail both (a) PASSES a correctly-painted gridinit
and (b) FAILS the exact bug classes seen in practice: an empty field (stale C++
binary), a momentum sign flip (exotic mismatch), and a wrong U(1) frequency.

Self-contained: builds a synthetic gridinit from the Python reference painter, so
no GPU / GRTresna solve is required.  An optional test also checks the real torus
gridinit when the artifact is present.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.io.gridinit import read_gridinit, write_gridinit
from grteclyn_wrapper.grtresna.matter.profile_contract import (
    PROFILE_KIND_TO_INT,
    MatterProfileSpec,
    check_gridinit_matches_spec,
    reference_paint,
)

_MATTER = ("phi", "phi2", "Pi", "Pi2")
_COMPS = ["chi", "lapse", *_MATTER]


def _gaussian_spec(**kw) -> MatterProfileSpec:
    base = dict(
        kind="gaussian", amp=0.1, omega=0.25, mass=0.5, m_az=1,
        exotic=True, winding=True, width=4.0, center=(16.0, 16.0, 16.0),
    )
    base.update(kw)
    return MatterProfileSpec(**base)


def _paint_synthetic_gridinit(
    tmp: Path, spec: MatterProfileSpec, *, n: int = 32, L: float = 32.0,
    mutate=None,
) -> Path:
    """Paint a synthetic gridinit from the reference painter (optionally mutated)."""
    dx = L / n
    xs = (np.arange(n) + 0.5) * dx
    Z, Y, X = np.meshgrid(xs, xs, xs, indexing="ij")
    ref = reference_paint(spec, X, Y, Z, alpha=1.0)
    data = np.zeros((n, n, n, len(_COMPS)), dtype=np.float64)
    data[..., _COMPS.index("chi")] = 1.0
    data[..., _COMPS.index("lapse")] = 1.0
    for name in _MATTER:
        data[..., _COMPS.index(name)] = ref[name]
    if mutate is not None:
        mutate(data)
    out = tmp / "synthetic.gridinit"
    write_gridinit(data, list(_COMPS), np.full(3, dx), np.zeros(3), out)
    return out


# --------------------------------------------------------------------------
# The single-source-of-truth keys
# --------------------------------------------------------------------------

def test_lump_keys_emit_correct_profile_int() -> None:
    assert PROFILE_KIND_TO_INT["torus"] == 4
    keys = _gaussian_spec(kind="torus", profile_path="/x/q.dat").grtresna_lump_keys(0)
    assert keys["lump0_profile"] == 4
    assert keys["lump0_winding"] == 1
    assert keys["lump0_exotic"] == 1
    assert keys["lump0_mode"] == 1
    assert keys["lump0_profile_path"] == "/x/q.dat"


def test_reference_paint_obeys_momentum_and_winding_rules() -> None:
    spec = _gaussian_spec(omega=0.3)
    x = np.array([18.0, 16.0]); y = np.array([16.0, 18.0]); z = np.array([16.0, 16.0])
    out = reference_paint(spec, x, y, z, alpha=1.0)
    # Pi2 = -(omega/alpha) phi1  and  Pi1 = +(omega/alpha) phi2  (alpha=1)
    assert np.allclose(out["Pi2"], -spec.omega * out["phi"])
    assert np.allclose(out["Pi"], spec.omega * out["phi2"])
    # winding: |Phi|^2 = phi1^2 + phi2^2 = modulus^2 (axisymmetric)
    mod2 = out["phi"] ** 2 + out["phi2"] ** 2
    assert np.all(mod2 >= 0.0)


# --------------------------------------------------------------------------
# The rail PASSES a correctly-painted field
# --------------------------------------------------------------------------

def test_contract_passes_for_correct_gridinit(tmp_path) -> None:
    spec = _gaussian_spec()
    gi = _paint_synthetic_gridinit(tmp_path, spec)
    rep = check_gridinit_matches_spec(gi, spec, tol=1e-6)
    assert rep.passed, str(rep)
    for name in _MATTER:
        assert rep.rel_err[name] < 1e-6


# --------------------------------------------------------------------------
# The rail CATCHES the real bug classes
# --------------------------------------------------------------------------

def test_contract_catches_empty_field_stale_binary(tmp_path) -> None:
    """Stale binary painted f==0 but 'converged' -> must FAIL, not pass silently."""
    spec = _gaussian_spec()

    def zero_matter(data):
        for name in _MATTER:
            data[..., _COMPS.index(name)] = 0.0

    gi = _paint_synthetic_gridinit(tmp_path, spec, mutate=zero_matter)
    rep = check_gridinit_matches_spec(gi, spec, tol=5e-2)
    assert not rep.passed
    assert rep.rel_err["phi"] > 0.9  # ~fully wrong


def test_contract_catches_momentum_sign_flip(tmp_path) -> None:
    """Exotic/sign mismatch flips Pi -> must FAIL on the momentum channels."""
    spec = _gaussian_spec()

    def flip_pi(data):
        data[..., _COMPS.index("Pi")] *= -1.0
        data[..., _COMPS.index("Pi2")] *= -1.0

    gi = _paint_synthetic_gridinit(tmp_path, spec, mutate=flip_pi)
    rep = check_gridinit_matches_spec(gi, spec, tol=5e-2)
    assert not rep.passed
    assert rep.rel_err["Pi2"] > 1.0  # sign flip -> ~2x peak error
    # the field channels themselves are unaffected by the momentum flip
    assert rep.rel_err["phi"] < 1e-6


def test_contract_catches_wrong_frequency(tmp_path) -> None:
    """A gridinit painted at omega=0.5 checked against a spec at omega=0.25."""
    painted = _gaussian_spec(omega=0.5)
    gi = _paint_synthetic_gridinit(tmp_path, painted)
    checked = _gaussian_spec(omega=0.25)  # wrong frequency in the spec
    rep = check_gridinit_matches_spec(gi, checked, tol=5e-2)
    assert not rep.passed
    # phi unaffected by omega; only Pi = omega*phi differs by factor 2
    assert rep.rel_err["phi"] < 1e-6
    assert rep.rel_err["Pi2"] > 0.4


# --------------------------------------------------------------------------
# Optional: the real torus gridinit (skipped if the artifact is absent)
# --------------------------------------------------------------------------

_TORUS_ID = Path(
    "/home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/rotating_torus_id/"
    "torus_m1_om0p250_kappa1p00_dx0p5_L64_lam170_mu614450_exotic"
)


@pytest.mark.skipif(
    not (_TORUS_ID / "initial_data.gridinit").is_file(),
    reason="real torus gridinit artifact not present",
)
def test_real_torus_gridinit_matches_solver() -> None:
    from grteclyn_wrapper.grtresna.profiles.qball_torus import read_torus_profile

    tor = read_torus_profile(_TORUS_ID / "qball_torus.dat")
    spec = MatterProfileSpec(
        kind="torus", amp=tor.f_max, omega=tor.omega, mass=tor.mass,
        lam=tor.lam, mu6=tor.mu, m_az=tor.m_az, exotic=True,
        center=(32.0, 32.0, 32.0), f_max=tor.f_max,
    )
    rep = check_gridinit_matches_spec(
        _TORUS_ID / "initial_data.gridinit", spec, tol=5e-2, torus=tor
    )
    assert rep.passed, str(rep)
