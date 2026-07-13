"""Tests for the rotating-wormhole kappa (amplitude-reduction) ID family.

Two layers, matching the repo's convention of always-on pure-function tests plus
fixture-gated tests over real solved artifacts:

1. Driver logic (always runs): the ``solve_kappa_family`` script's param-rewriting
   -- amplitude scaling ``f -> kappa*f`` and output-path normalisation -- and the
   L2 guard that its gridinit-export grid equals the evolution's level-0 dx.

2. Winding gridinit properties (skipped when no solved fixture is present): the
   solved ``.gridinit`` must carry the *genuine phase winding* fix (both complex
   channels + two-channel momentum populated, |Phi|^2 axisymmetric, K=0 maximal
   slicing, dx == evolution level-0 dx).
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DRIVER_PATH = REPO_ROOT / "scripts" / "wormhole" / "id" / "solve_kappa_family.py"
EVO_PARAMS = (
    REPO_ROOT.parent
    / "Examples"
    / "RotatingWormholeCollapse"
    / "params_rotating_grtresna_exotic.txt"
)
_ID_ROOT = REPO_ROOT.parent / "runs" / "rotating_wormhole_id"
# Prefer the high-res dx=0.5 (N=128) family; fall back to the dx=1.0 (N=64) one.
KAPPA_GRIDINIT = next(
    (
        p
        for p in (
            _ID_ROOT / "rotwh_omega_p0p05_m1_kappa_1p00_dx0p5" / "initial_data.gridinit",
            _ID_ROOT / "rotwh_omega_p0p05_m1_kappa_1p00" / "initial_data.gridinit",
        )
        if p.is_file()
    ),
    _ID_ROOT / "rotwh_omega_p0p05_m1_kappa_1p00_dx0p5" / "initial_data.gridinit",
)


def _load_driver():
    spec = importlib.util.spec_from_file_location("solve_kappa_family", DRIVER_PATH)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BASE_PARAMS_TEXT = """\
output_path = Outputs_rotwh/
N = 64 64 64
L = 64
num_lumps = 1
lump0_amp = 0.1
lump0_omega = 0.05
lump0_mode = 1
lump0_exotic = 1
lump0_winding = 1
bh1_bare_mass = 0.25
"""


# --------------------------------------------------------------------------- #
# 1. Driver logic (always runs)                                               #
# --------------------------------------------------------------------------- #
def test_driver_reads_base_amplitude() -> None:
    drv = _load_driver()
    assert drv._read_base_amp(BASE_PARAMS_TEXT) == pytest.approx(0.1)


@pytest.mark.parametrize("kappa", [1.0, 0.9, 0.7, 0.5])
def test_driver_scales_amplitude_by_kappa(tmp_path: Path, kappa: float) -> None:
    """f -> kappa*f: only lump0_amp changes; other lump knobs are untouched."""
    drv = _load_driver()
    base_amp = drv._read_base_amp(BASE_PARAMS_TEXT)
    dst = tmp_path / "params.txt"
    drv._write_scaled_params(dst, BASE_PARAMS_TEXT, kappa * base_amp)
    text = dst.read_text()

    m = re.search(r"^\s*lump0_amp\s*=\s*([0-9eE.+-]+)", text, re.MULTILINE)
    assert m is not None
    assert float(m.group(1)) == pytest.approx(kappa * base_amp)

    # Winding/exotic/throat structure must survive the rewrite unchanged.
    assert "lump0_winding = 1" in text
    assert "lump0_exotic = 1" in text
    assert "lump0_mode = 1" in text
    assert "bh1_bare_mass = 0.25" in text


def test_driver_normalises_output_path(tmp_path: Path) -> None:
    """The base file writes to Outputs_rotwh/; the driver must force Outputs/ so
    it always finds Outputs/InitialDataFinal.3d.hdf5 afterwards."""
    drv = _load_driver()
    dst = tmp_path / "params.txt"
    drv._write_scaled_params(dst, BASE_PARAMS_TEXT, 0.05)
    text = dst.read_text()
    assert re.search(r"^\s*output_path\s*=\s*Outputs/\s*$", text, re.MULTILINE)
    assert "Outputs_rotwh" not in text


def test_driver_parse_convergence(tmp_path: Path) -> None:
    drv = _load_driver()
    err = tmp_path / "Ham_and_Mom_errors.txt"
    err.write_text(
        "NL_iteration       Ham error [%]       Mom error [%]\n"
        "0                   1.000000e+02        1.000000e+02\n"
        "1                   4.035989e+01        1.500282e+01\n"
        "7                   7.620815e-01        6.308495e-01\n"
    )
    ham, mom = drv._parse_convergence(err)
    assert ham == pytest.approx(0.7620815)
    assert mom == pytest.approx(0.6308495)
    assert drv._parse_convergence(tmp_path / "missing.txt") is None


@pytest.mark.skipif(not EVO_PARAMS.is_file(), reason="evolution params not present")
def test_driver_grid_matches_evolution_level0_dx() -> None:
    """L2 lesson: the gridinit must be emitted at the exact evolution level-0 dx.
    The driver's EVO_L/EVO_N/center must equal the RotatingWormholeCollapse
    evolution params so dx_gridinit == dx_evolution (no interpolation kinks)."""
    drv = _load_driver()
    text = EVO_PARAMS.read_text()

    def _f(key: str) -> float:
        m = re.search(rf"^\s*{key}\s*=\s*([0-9eE.+-]+)", text, re.MULTILINE)
        assert m is not None, f"{key} not found in {EVO_PARAMS}"
        return float(m.group(1))

    evo_L = _f("L")
    evo_N = _f("N1")
    assert drv.EVO_L == pytest.approx(evo_L)
    assert drv.EVO_N == int(evo_N)
    # dx match is the actual L2 invariant.
    assert drv.EVO_L / drv.EVO_N == pytest.approx(evo_L / evo_N)

    cm = re.search(r"^\s*center\s*=\s*([-0-9.eE]+)\s+([-0-9.eE]+)\s+([-0-9.eE]+)",
                   text, re.MULTILINE)
    assert cm is not None
    evo_center = tuple(float(cm.group(i)) for i in (1, 2, 3))
    assert drv.EVO_CENTER == pytest.approx(evo_center)


# --------------------------------------------------------------------------- #
# 2. Winding gridinit properties (fixture-gated)                              #
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not KAPPA_GRIDINIT.is_file(),
                    reason="no solved kappa=1.0 gridinit fixture; run solve_kappa_family.sh")
def test_solved_gridinit_is_winding_and_constraint_clean() -> None:
    from grteclyn_wrapper.grtresna.io import read_gridinit

    g = read_gridinit(KAPPA_GRIDINIT)
    comp = {n: g.data[..., i] for i, n in enumerate(g.comp_names)}

    # Resolution-agnostic: isotropic dx = L/N with L=64, and target_center
    # (32,32,0) with L=64 -> origin (0,0,-32). Works for both the dx=1.0 (N=64)
    # and dx=0.5 (N=128) families.
    n = g.data.shape[0]
    dx_expected = 64.0 / n
    np.testing.assert_allclose(g.dx_xyz, [dx_expected] * 3)
    np.testing.assert_allclose(g.origin, [0.0, 0.0, -32.0])

    for name in ("chi", "phi", "phi2", "Pi", "Pi2", "K"):
        assert name in comp, f"missing {name} in gridinit"
        assert np.isfinite(comp[name]).all(), f"{name} has non-finite entries"

    # Throat: a genuine wormhole conformal factor departs from flat (chi != 1).
    assert comp["chi"].min() < 0.99
    assert comp["chi"].max() > 1.0

    # Maximal slicing: K == 0 everywhere.
    assert np.max(np.abs(comp["K"])) == pytest.approx(0.0, abs=1e-12)

    # GENUINE PHASE WINDING (the B1-B3 fix): both complex channels are excited
    # with comparable amplitude. The failed real-modulation ansatz (L3) has
    # phi2 == 0 identically, which this rejects.
    amp_phi = np.max(np.abs(comp["phi"]))
    amp_phi2 = np.max(np.abs(comp["phi2"]))
    assert amp_phi > 1e-6 and amp_phi2 > 1e-6
    assert amp_phi2 == pytest.approx(amp_phi, rel=0.2)

    # TWO-CHANNEL MOMENTUM (constraint-clean spin): Pi and Pi2 both populated,
    # scaled ~ omega * phi (omega = 0.05).
    amp_pi = np.max(np.abs(comp["Pi"]))
    amp_pi2 = np.max(np.abs(comp["Pi2"]))
    assert amp_pi > 1e-8 and amp_pi2 > 1e-8
    assert amp_pi == pytest.approx(0.05 * amp_phi, rel=0.35)


@pytest.mark.skipif(not KAPPA_GRIDINIT.is_file(),
                    reason="no solved kappa=1.0 gridinit fixture; run solve_kappa_family.sh")
def test_solved_gridinit_modulus_is_axisymmetric() -> None:
    """|Phi|^2 = phi^2 + phi2^2 for a phase-winding ansatz f e^{i m phi_az} is
    independent of azimuth -- no four-lobe pattern (L3). Sample a ring at fixed
    radius in the throat plane and require a small azimuthal coefficient of
    variation."""
    from grteclyn_wrapper.grtresna.io import read_gridinit

    g = read_gridinit(KAPPA_GRIDINIT)
    comp = {n: g.data[..., i] for i, n in enumerate(g.comp_names)}
    mod2 = comp["phi"] ** 2 + comp["phi2"] ** 2  # (nz, ny, nx)

    nz, ny, nx = mod2.shape
    dx = g.dx_xyz
    origin = g.origin
    # Physical coordinates of cell centres.
    xs = origin[0] + (np.arange(nx) + 0.5) * dx[0]
    ys = origin[1] + (np.arange(ny) + 0.5) * dx[1]
    zs = origin[2] + (np.arange(nz) + 0.5) * dx[2]
    # Throat centre is the evolution center (32, 32, 0).
    cx, cy, cz = 32.0, 32.0, 0.0
    kz = int(np.argmin(np.abs(zs - cz)))  # z-plane through the throat

    plane = mod2[kz]  # (ny, nx)
    XX, YY = np.meshgrid(xs - cx, ys - cy)  # (ny, nx)
    R = np.sqrt(XX ** 2 + YY ** 2)

    # Ring at a mid radius where the field has support.
    r_target = 6.0
    ring = np.abs(R - r_target) < 1.0
    vals = plane[ring]
    vals = vals[vals > 0]
    assert vals.size >= 8, "not enough ring samples to assess axisymmetry"

    cov = float(np.std(vals) / np.mean(vals))
    # A single-real cos(m phi) modulation gives O(1) azimuthal variation; the
    # winding modulus should be far more uniform.
    assert cov < 0.35, f"|Phi|^2 not axisymmetric on the throat ring (cov={cov:.2f})"


# --------------------------------------------------------------------------- #
# 4. Q-ball (Rung 1a) initial-data path                                       #
# --------------------------------------------------------------------------- #
QBALL_BASE_PARAMS_TEXT = """\
output_path = Outputs_rotwh/
N = 64 64 64
L = 64
scalar_mass = 0.0
scalar_lambda = 0.0
scalar_mu = 0.0
num_lumps = 1
lump0_amp = 0.1
lump0_omega = 0.05
lump0_mode = 1
lump0_exotic = 1
lump0_winding = 1
lump0_profile = 0
bh1_bare_mass = 0.25
"""


def _load_driver_with_env(monkeypatch, **env):
    for k, v in env.items():
        monkeypatch.setenv(k, str(v))
    return _load_driver()


def test_qball_path_sets_couplings_and_profile(tmp_path, monkeypatch) -> None:
    """LAMBDA>0 switches the throat lump to the solved Q-ball eigenstate:
    scalar_lambda/scalar_mu are written, lump0_profile becomes 3, and a
    qball_profile.dat + qball_profile_path are emitted (Rung 1a wiring)."""
    drv = _load_driver_with_env(
        monkeypatch, MASS=0.5, LAMBDA=1.0, MU6=0.5, LUMP_OMEGA=0.05, RES_N=64
    )
    dst = tmp_path / "params.txt"
    drv._write_scaled_params(dst, QBALL_BASE_PARAMS_TEXT, 0.1)
    text = dst.read_text()

    assert re.search(r"^\s*scalar_lambda\s*=\s*1\b", text, re.MULTILINE)
    assert re.search(r"^\s*scalar_mu\s*=\s*0\.5\b", text, re.MULTILINE)
    assert re.search(r"^\s*lump0_profile\s*=\s*3\s*$", text, re.MULTILINE)
    m = re.search(r"^\s*qball_profile_path\s*=\s*(\S+)", text, re.MULTILINE)
    assert m is not None
    prof = Path(m.group(1))
    assert prof.is_file()
    # Tabulated profile is normalisable (phi0(0) > 0) and localized.
    lines = [ln for ln in prof.read_text().splitlines() if not ln.startswith("#")]
    r0, phi0 = (float(x) for x in lines[0].split())
    assert r0 == pytest.approx(0.0)
    assert phi0 > 0.0


def test_winding_frequency_and_mode_drive_params_and_tag(tmp_path, monkeypatch) -> None:
    drv = _load_driver_with_env(
        monkeypatch,
        MASS=0.5,
        LAMBDA=120.0,
        MU6=7200.0,
        LUMP_OMEGA=0.4,
        LUMP_MODE=2,
        RES_N=64,
    )
    dst = tmp_path / "params.txt"
    drv._write_scaled_params(dst, QBALL_BASE_PARAMS_TEXT, 0.1)
    text = dst.read_text()

    assert re.search(r"^\s*lump0_omega\s*=\s*0\.4\s*$", text, re.MULTILINE)
    assert re.search(r"^\s*lump0_mode\s*=\s*2\s*$", text, re.MULTILINE)
    assert drv._run_tag(1.0) == (
        "rotwh_omega_p0p40_m2_kappa_1p00_dx1_mass0p5_"
        "qball_lam120_mu67200"
    )


def test_massless_default_leaves_profile_gaussian(tmp_path, monkeypatch) -> None:
    """Backward compat: with LAMBDA unset the lump stays the analytic Gaussian
    (profile 0) and no Q-ball couplings/profile are injected."""
    drv = _load_driver_with_env(monkeypatch, MASS=0.0, RES_N=64)
    # Clear any Q-ball env from other tests in the same process.
    monkeypatch.delenv("LAMBDA", raising=False)
    monkeypatch.delenv("MU6", raising=False)
    drv = _load_driver()
    dst = tmp_path / "params.txt"
    drv._write_scaled_params(dst, QBALL_BASE_PARAMS_TEXT, 0.1)
    text = dst.read_text()
    assert re.search(r"^\s*lump0_profile\s*=\s*0\s*$", text, re.MULTILINE)
    assert "qball_profile_path" not in text
