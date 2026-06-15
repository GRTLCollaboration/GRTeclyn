from grteclyn_wrapper.search.optimize import (
    build_grtresna_config,
    build_search_space,
)


def test_ring_ansatz_search_space_is_low_dimensional() -> None:
    space = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="ring")

    assert len(space) == 14
    assert {dim.param_key for dim in space} >= {
        "grtresna_ring_amp",
        "grtresna_ring_radius",
        "grtresna_ring_tangential_velocity",
        "grtresna_ring_exotic_fraction",
    }


def test_ring_ansatz_expands_to_lumps() -> None:
    overrides = {
        "grtresna_ring_lumps": 5,
        "grtresna_ring_amp": 0.15,
        "grtresna_ring_width": 3.0,
        "grtresna_ring_radius": 4.0,
        "grtresna_ring_phase": 0.0,
        "grtresna_ring_tangential_velocity": 0.4,
        "grtresna_ring_radial_velocity": 0.0,
        "grtresna_ring_vertical_velocity": 0.0,
        "grtresna_ring_exotic_fraction": 0.4,
        "grtresna_ring_exotic_phase": 0.0,
        "grtresna_ring_mode": 1.0,
    }

    cfg = build_grtresna_config(overrides)

    assert len(cfg.lumps) == 5
    assert sum(lump["exotic"] for lump in cfg.lumps) == 2
    assert cfg.lumps[0]["center"][0] == 4.0
    assert cfg.lumps[0]["velocity"][1] == 0.4
    assert cfg.maximal_slicing
