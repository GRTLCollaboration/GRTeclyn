"""Tests for the spherical-harmonic (SH) matter ansatz."""

import math

from grteclyn_wrapper.__main__ import build_parser
from grteclyn_wrapper.cli.grtresna_context import build_grtresna_search_context
from grteclyn_wrapper.grtresna.sh_fields import (
    cartesian_to_spherical,
    eval_sh_modulation,
    real_sph_harm,
    sh_ell_m,
    sh_flat_index,
)
from grteclyn_wrapper.search.optimize import (
    SH_DEFAULT_ELL_MAX,
    build_grtresna_config,
    build_search_space,
    grtresna_sh_search_space,
)


# ---------------------------------------------------------------------------
# Spherical harmonic math tests
# ---------------------------------------------------------------------------


def test_sh_flat_index_roundtrip() -> None:
    """flat index -> (ell, m) -> flat index is identity."""
    for ell in range(6):
        for m in range(-ell, ell + 1):
            idx = sh_flat_index(ell, m)
            assert sh_ell_m(idx) == (ell, m), f"roundtrip failed for ({ell}, {m})"


def test_y00_is_constant() -> None:
    """Y_{0,0} = 1 / (2 sqrt(pi)) everywhere."""
    expected = 1.0 / (2.0 * math.sqrt(math.pi))
    for theta in (0.0, 0.5, 1.0, math.pi):
        for phi in (0.0, 1.0, 3.0, 2 * math.pi):
            val = real_sph_harm(0, 0, theta, phi)
            assert abs(val - expected) < 1e-12, f"Y00({theta},{phi}) = {val}"


def test_y10_is_cos_theta() -> None:
    """Y_{1,0} proportional to cos(theta)."""
    norm = math.sqrt(3.0 / (4.0 * math.pi))
    for theta in (0.0, 0.3, math.pi / 2, 2.5, math.pi):
        expected = norm * math.cos(theta)
        val = real_sph_harm(1, 0, theta, 0.0)
        assert abs(val - expected) < 1e-10, f"Y10({theta}) = {val} != {expected}"


def test_y11_is_sin_theta_cos_phi() -> None:
    """Y_{1,1} proportional to sin(theta) cos(phi)."""
    norm = math.sqrt(3.0 / (4.0 * math.pi))
    for theta in (0.3, math.pi / 2, 2.0):
        for phi in (0.0, 1.0, math.pi):
            expected = -norm * math.sin(theta) * math.cos(phi)
            val = real_sph_harm(1, 1, theta, phi)
            assert abs(val - expected) < 1e-10, (
                f"Y11({theta},{phi}) = {val} != {expected}"
            )


def test_modulation_is_one_when_all_coeffs_zero() -> None:
    """With zero modulation coefficients, factor is exactly 1.0."""
    coeffs = [0.0] * 25
    for theta in (0.0, 1.0, math.pi):
        for phi in (0.0, 1.5, 3.0):
            mod = eval_sh_modulation(coeffs, theta, phi, ell_max=4)
            assert abs(mod - 1.0) < 1e-14


def test_modulation_is_clamped_nonnegative() -> None:
    """Large negative coefficients cannot produce negative modulation."""
    coeffs = [0.0] * 25
    coeffs[1] = -100.0  # huge negative Y_{1,-1}
    for theta in (0.0, 0.5, 1.0, math.pi):
        mod = eval_sh_modulation(coeffs, theta, 0.0, ell_max=4)
        assert mod >= 0.0


# ---------------------------------------------------------------------------
# Search space tests
# ---------------------------------------------------------------------------


def test_sh_search_space_dimensionality() -> None:
    """SH search space has (ell_max+1)^2 - 1 modulation + global dims."""
    space = grtresna_sh_search_space(ell_max=4)
    keys = {d.param_key for d in space}

    # 24 modulation coefficients + 1 base amp = 25 SH dims.
    sh_keys = {k for k in keys if k.startswith("grtresna_sh_c")}
    assert len(sh_keys) == 24

    # Required global dims.
    assert "grtresna_sh_amp" in keys
    assert "grtresna_sh_radius" in keys
    assert "grtresna_sh_width" in keys
    assert "grtresna_sh_toroidal_velocity" in keys
    assert "grtresna_sh_omega" in keys
    assert "grtresna_sh_exotic_fraction" in keys
    assert "grtresna_shift_seed" in keys
    assert "grtresna_scalar_mass" in keys
    assert "grtresna_sh_static" in keys


def test_sh_search_space_ell3_is_smaller() -> None:
    """ell_max=3 gives fewer dims than ell_max=4."""
    s3 = grtresna_sh_search_space(ell_max=3)
    s4 = grtresna_sh_search_space(ell_max=4)
    sh3 = {k.param_key for k in s3 if k.param_key.startswith("grtresna_sh_c")}
    sh4 = {k.param_key for k in s4 if k.param_key.startswith("grtresna_sh_c")}
    assert len(sh3) == 15  # (3+1)^2 - 1
    assert len(sh4) == 24  # (4+1)^2 - 1
    assert len(s3) < len(s4)


def test_sh_search_space_wired_into_build() -> None:
    """build_search_space dispatches 'sh' ansatz correctly."""
    space = build_search_space(grtresna=True, grtresna_ansatz="sh")
    keys = {d.param_key for d in space}
    assert "grtresna_sh_amp" in keys
    assert "grtresna_sh_c1" in keys
    # Should NOT have shell keys.
    assert "grtresna_shell_amp" not in keys


# ---------------------------------------------------------------------------
# Config expansion tests
# ---------------------------------------------------------------------------


def _sh_base_overrides(num_lumps: int = 8, **extra: float) -> dict:
    """Minimal SH overrides for config expansion."""
    ov = {
        "grtresna_sh_lumps": num_lumps,
        "grtresna_sh_amp": 0.13,
        "grtresna_sh_width": 2.4,
        "grtresna_sh_radius": 4.0,
        "grtresna_sh_thickness": 0.0,
        "grtresna_sh_toroidal_velocity": 0.4,
        "grtresna_sh_poloidal_velocity": 0.0,
        "grtresna_sh_radial_velocity": 0.0,
        "grtresna_sh_omega": 0.0,
        "grtresna_sh_exotic_fraction": 0.0,
        "grtresna_sh_exotic_phase": 0.0,
        "grtresna_sh_static": 0.0,
        "grtresna_scalar_mass": 0.6,
        "grtresna_scalar_lambda": 0.0,
        "grtresna_shift_seed": 0.0,
    }
    ov.update(extra)
    return ov


def test_sh_expands_to_lumps_on_sphere() -> None:
    """SH ansatz produces N lumps on a spherical shell."""
    overrides = _sh_base_overrides(num_lumps=8)
    cfg = build_grtresna_config(overrides)

    assert len(cfg.lumps) == 8

    # All lumps should sit on the requested radius.
    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        r = math.sqrt(cx * cx + cy * cy + cz * cz)
        assert abs(r - 4.0) < 1e-6, f"lump at r={r}, expected 4.0"

    # Lumps should cover the full sphere (spread on all 3 axes).
    for axis in range(3):
        spread = max(abs(lump["center"][axis]) for lump in cfg.lumps)
        assert spread > 0.5, f"lumps collapse along axis {axis}"


def test_sh_uniform_gives_equal_amplitudes() -> None:
    """With zero modulation coefficients, all lumps get equal amp."""
    overrides = _sh_base_overrides(num_lumps=6)
    cfg = build_grtresna_config(overrides)
    amps = [lump["amp"] for lump in cfg.lumps]
    # All should be base amp (0.13) modulated by Y00*0 + 1 = 1.0.
    for a in amps:
        assert abs(a - 0.13) < 1e-10, f"amp={a}, expected 0.13"


def test_sh_modulation_creates_per_lump_variation() -> None:
    """Non-zero SH coefficients produce different amplitudes per lump."""
    overrides = _sh_base_overrides(num_lumps=8)
    # Add a strong ℓ=2, m=0 modulation (axial quadrupole).
    idx_20 = sh_flat_index(2, 0)  # idx = 6
    overrides[f"grtresna_sh_c{idx_20}"] = 0.6

    cfg = build_grtresna_config(overrides)
    amps = [lump["amp"] for lump in cfg.lumps]

    # Not all equal anymore.
    assert max(amps) - min(amps) > 0.01, (
        f"expected per-lump variation, got amps={amps}"
    )


def test_sh_higher_modes_give_finer_structure() -> None:
    """Higher ell modes create more angular variation."""
    # ℓ=1 (dipole) variation.
    ov1 = _sh_base_overrides(num_lumps=12)
    idx_10 = sh_flat_index(1, 0)
    ov1[f"grtresna_sh_c{idx_10}"] = 0.5

    # ℓ=4 (hexadecapole) variation.
    ov4 = _sh_base_overrides(num_lumps=12)
    idx_40 = sh_flat_index(4, 0)
    ov4[f"grtresna_sh_c{idx_40}"] = 0.5

    cfg1 = build_grtresna_config(ov1)
    cfg4 = build_grtresna_config(ov4)

    amps1 = [l["amp"] for l in cfg1.lumps]
    amps4 = [l["amp"] for l in cfg4.lumps]

    # Both should have variation.
    assert max(amps1) - min(amps1) > 0.01
    assert max(amps4) - min(amps4) > 0.01

    # ℓ=4 should have more sign changes (more local extrema) than ℓ=1.
    # Count sign changes in the deviation from mean.
    def sign_changes(vals: list[float]) -> int:
        mean = sum(vals) / len(vals)
        devs = [v - mean for v in vals]
        return sum(1 for i in range(1, len(devs)) if devs[i] * devs[i - 1] < 0)

    sc1 = sign_changes(amps1)
    sc4 = sign_changes(amps4)
    assert sc4 > sc1, (
        f"ℓ=4 should have more sign changes ({sc4}) than ℓ=1 ({sc1})"
    )


def test_sh_static_zeroes_velocities() -> None:
    """sh_static=1 zeroes all lump velocities and omega."""
    overrides = _sh_base_overrides(
        grtresna_sh_toroidal_velocity=0.5,
        grtresna_sh_omega=0.2,
        grtresna_sh_static=1.0,
    )
    cfg = build_grtresna_config(overrides)

    for lump in cfg.lumps:
        vx, vy, vz = lump["velocity"]
        assert vx == 0.0 and vy == 0.0 and vz == 0.0, (
            f"velocity not zeroed: {lump['velocity']}"
        )
        assert lump["omega"] == 0.0, f"omega not zeroed: {lump['omega']}"


def test_sh_dynamics_on_has_nonzero_velocities() -> None:
    """sh_static=0 preserves velocities."""
    overrides = _sh_base_overrides(
        grtresna_sh_toroidal_velocity=0.5,
        grtresna_sh_omega=0.2,
        grtresna_sh_static=0.0,
    )
    cfg = build_grtresna_config(overrides)

    has_velocity = False
    for lump in cfg.lumps:
        vx, vy, vz = lump["velocity"]
        if abs(vx) + abs(vy) + abs(vz) > 0.01:
            has_velocity = True
    assert has_velocity, "no lumps have non-zero velocity with sh_static=0"


def test_sh_exotic_fraction_assigns_exotic_flags() -> None:
    """exotic_fraction > 0 marks some lumps as exotic."""
    overrides = _sh_base_overrides(
        num_lumps=6,
        grtresna_sh_exotic_fraction=0.5,
    )
    cfg = build_grtresna_config(overrides)
    exotic_count = sum(l["exotic"] for l in cfg.lumps)
    assert exotic_count == 3, f"expected 3 exotic lumps, got {exotic_count}"


# ---------------------------------------------------------------------------
# CLI / context wiring tests
# ---------------------------------------------------------------------------


def test_cli_accepts_sh_ansatz() -> None:
    """CLI parser accepts --grtresna-ansatz sh."""
    parser = build_parser()
    args = parser.parse_args([
        "--consume-plotfiles",
        "qd",
        "--grtresna",
        "--grtresna-ansatz",
        "sh",
    ])
    assert args.grtresna_ansatz == "sh"


def test_context_injects_sh_lumps() -> None:
    """grtresna_context sets grtresna_sh_lumps in base_overrides."""
    parser = build_parser()
    args = parser.parse_args([
        "--consume-plotfiles",
        "qd",
        "--grtresna",
        "--grtresna-ansatz",
        "sh",
        "--grtresna-lumps",
        "10",
    ])
    ctx = build_grtresna_search_context(args, {})
    assert ctx.base_overrides.get("grtresna_sh_lumps") == 10
    keys = {d.param_key for d in ctx.search_space}
    assert "grtresna_sh_amp" in keys
    assert "grtresna_shell_amp" not in keys


# ---------------------------------------------------------------------------
# Serialization test
# ---------------------------------------------------------------------------


def test_sh_lump_lines_include_per_lump_data() -> None:
    """Config from SH produces lump lines with per-lump amp variation."""
    from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, _lump_lines

    overrides = _sh_base_overrides(num_lumps=5)
    idx_20 = sh_flat_index(2, 0)
    overrides[f"grtresna_sh_c{idx_20}"] = 0.5

    cfg = build_grtresna_config(overrides, GRTresnaConfig())
    text = "\n".join(_lump_lines(cfg))

    assert "num_lumps = 5" in text
    assert "lump0_amp" in text
    assert "lump4_amp" in text
    assert "lump0_center" in text
    assert "lump0_velocity" in text
