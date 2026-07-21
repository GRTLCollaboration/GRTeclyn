"""Site layout is plug-and-play: .env overlay + injectable SiteLayout."""

from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.core.site_paths import (
    DefaultSiteLayout,
    EnvOverlay,
    expand_env_value,
    parse_env_file,
)


def test_parse_env_file_supports_export_and_quotes() -> None:
    text = """
# comment
export FOO=/tmp/a
BAR='/tmp/b'
BAZ="/tmp/c"
"""
    assert parse_env_file(text) == {
        "FOO": "/tmp/a",
        "BAR": "/tmp/b",
        "BAZ": "/tmp/c",
    }


def test_expand_env_value_resolves_nested_refs() -> None:
    lookup = {"SIM_ROOT": "/sim", "HOME": "/home/user"}
    assert expand_env_value("${SIM_ROOT}/GRTeclyn", lookup) == "/sim/GRTeclyn"
    assert expand_env_value("$HOME/envs/grtresna", lookup) == "/home/user/envs/grtresna"


def test_env_overlay_sets_missing_keys_only(tmp_path: Path) -> None:
    env_file = tmp_path / ".env"
    env_file.write_text(
        "SIM_ROOT=/overlay/sim\nGRTECLYN_ROOT=${SIM_ROOT}/GRTeclyn\nGRTRESNA_ENV=/already/ignored\n",
        encoding="utf-8",
    )
    target = {"GRTRESNA_ENV": "/already/set"}
    overlay = EnvOverlay(paths=(env_file,))
    assert overlay.apply(target) == env_file
    assert target["SIM_ROOT"] == "/overlay/sim"
    assert target["GRTECLYN_ROOT"] == "/overlay/sim/GRTeclyn"
    assert target["GRTRESNA_ENV"] == "/already/set"


def test_default_layout_honors_injected_env(tmp_path: Path) -> None:
    sim = tmp_path / "sim"
    grteclyn = sim / "GRTeclyn"
    grtresna = sim / "GRTresna"
    chombo = sim / "Chombo" / "lib"
    openmpi = sim / "local" / "openmpi"
    for path in (grteclyn, grtresna, chombo, openmpi):
        path.mkdir(parents=True)

    env = {
        "GRTECLYN_ROOT": str(grteclyn),
        "SIM_ROOT": str(sim),
        "GRTRESNA_ENV": str(tmp_path / "missing_env"),
    }
    layout = DefaultSiteLayout(
        env=env,
        overlay=EnvOverlay(paths=(tmp_path / "missing.env",)),
    )
    assert layout.grteclyn_root() == grteclyn.resolve()
    assert layout.sim_root() == sim.resolve()
    assert layout.grtresna_root() == grtresna.resolve()
    assert layout.chombo_home() == chombo.resolve()
    assert layout.openmpi_root() == openmpi.resolve()
    assert layout.grtresna_env() is None  # path does not exist


def test_grteclyn_root_facade_finds_checkout() -> None:
    from grteclyn_wrapper.core.site_paths import grteclyn_root

    root = grteclyn_root()
    assert (root / "Examples" / "RadialRecipe").is_dir()
    assert (root / "grteclyn-wrapper").is_dir()
