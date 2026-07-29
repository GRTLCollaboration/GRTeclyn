"""A shell that searches no velocity must not be given one behind your back.

The boson-shell space builder strips the four matter-current dimensions for the
static families, but ``grtresna_shell_static`` stayed in the space with its
centre at 0 = moving.  Any candidate landing on 0 then picked up a hard-coded
0.35c toroidal current that no dimension declared and no trajectory record
mentioned.  Half a campaign, advertised as static, was flowing.
"""

from __future__ import annotations

from grteclyn_wrapper.search.optimize.config import (
    SHELL_CURRENT_DEFAULTS,
    _expand_shell_lumps_from_overrides,
    resolved_shell_currents,
    unrecorded_shell_currents,
)
from grteclyn_wrapper.search.optimize.spaces import (
    grtresna_boson_shell_search_space,
    grtresna_shell_search_space,
)

_TOROIDAL = "grtresna_shell_toroidal_velocity"
_VELOCITY_KEYS = tuple(SHELL_CURRENT_DEFAULTS)


def _keys(dims) -> set[str]:
    return {d.param_key for d in dims}


def test_the_toroidal_default_is_still_a_real_current() -> None:
    """Guards the premise: 0.35 is not a typo someone already zeroed."""
    assert SHELL_CURRENT_DEFAULTS[_TOROIDAL] == 0.35


def test_a_space_without_velocities_does_not_search_the_static_switch() -> None:
    dims = grtresna_boson_shell_search_space(static=True)
    keys = _keys(dims)
    assert not (keys & set(_VELOCITY_KEYS))
    assert "grtresna_shell_static" not in keys


def test_a_space_with_velocities_keeps_the_static_switch() -> None:
    dims = grtresna_boson_shell_search_space(static=False, matched_bounds=True)
    keys = _keys(dims)
    assert _TOROIDAL in keys
    assert "grtresna_shell_static" in keys


def test_the_scalar_shell_space_is_untouched() -> None:
    """The unmodified family still searches every current and the switch."""
    keys = _keys(grtresna_shell_search_space())
    assert set(_VELOCITY_KEYS) <= keys
    assert "grtresna_shell_static" in keys


def test_no_searched_velocity_means_no_velocity() -> None:
    """The exact broken configuration: static bit says 'moving', space says nothing."""
    overrides = {"grtresna_shell_static": 0.0, "grtresna_shell_amp": 0.1}
    resolved = resolved_shell_currents(overrides)
    assert resolved[_TOROIDAL] == 0.0
    assert all(resolved[key] == 0.0 for key in _VELOCITY_KEYS)


def test_a_searched_velocity_is_honoured() -> None:
    overrides = {"grtresna_shell_static": 0.0, _TOROIDAL: 0.2}
    assert resolved_shell_currents(overrides)[_TOROIDAL] == 0.2


def test_the_static_switch_still_stops_a_searched_velocity() -> None:
    overrides = {"grtresna_shell_static": 1.0, _TOROIDAL: 0.2}
    assert resolved_shell_currents(overrides)[_TOROIDAL] == 0.0


def test_an_undeclared_current_is_reported() -> None:
    """A pinned velocity leaves the record; a defaulted one used not to."""
    overrides = {_TOROIDAL: 0.35, "grtresna_shell_static": 0.0}
    assert unrecorded_shell_currents(overrides, frozenset()) == {}

    searched = frozenset({_TOROIDAL, "grtresna_shell_static"})
    assert unrecorded_shell_currents(overrides, searched) == {}


def test_the_lumps_themselves_come_out_still() -> None:
    """End of the chain: the actual lump velocities, not just the knobs.

    This is the assertion that fails on the old decoder -- it hands every lump a
    0.35c toroidal velocity for these overrides.
    """
    overrides = {
        "grtresna_shell_static": 0.0,
        "grtresna_shell_amp": 0.1,
        "grtresna_shell_width": 3.0,
        "grtresna_shell_lumps": 6.0,
    }

    def get_float(key: str, default: float) -> float:
        return float(overrides.get(key, default))

    lumps = _expand_shell_lumps_from_overrides(overrides, get_float)
    assert lumps
    for lump in lumps:
        assert all(abs(v) == 0.0 for v in lump["velocity"])


def test_the_recorded_currents_match_what_the_decoder_uses() -> None:
    """Whatever the record claims has to be what the lumps actually get."""
    overrides = {_TOROIDAL: 0.35, "grtresna_shell_poloidal_velocity": -0.1}
    resolved = resolved_shell_currents(overrides)
    assert resolved[_TOROIDAL] == 0.35
    assert resolved["grtresna_shell_poloidal_velocity"] == -0.1
    # Both were supplied, so neither is hidden.
    assert unrecorded_shell_currents(overrides, frozenset(overrides)) == {}
