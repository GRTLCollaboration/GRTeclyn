from __future__ import annotations

from grteclyn_wrapper.rl.seed import load_elite_overrides, pump_geometry_from_overrides


def test_pump_geometry_defaults() -> None:
    radius, width = pump_geometry_from_overrides({})
    assert radius == 5.0
    assert width == 1.5


def test_pump_geometry_from_overrides() -> None:
    radius, width = pump_geometry_from_overrides(
        {"recipe_basis_radius_max": 7.0, "recipe_basis_width": 2.0}
    )
    assert radius == 7.0
    assert width == 2.0
