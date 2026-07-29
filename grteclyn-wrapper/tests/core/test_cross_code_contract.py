"""GRTresna's config vs GRTeclyn's config vs the .gridinit between them.

WHY THIS FILE EXISTS
--------------------
Every run is configured twice, in two files, using two different names for the
same physics.  ``grtresna/params.txt`` says ``scalar_mu``; ``params.txt`` says
``recipe_scalar_mu``.  ``lump3_exotic = 1`` upstream means
``recipe_scalar_field_signs[3] = -1`` downstream.  Until now nothing compared
them.

A disagreement does not raise: GRTresna converges, GRTeclyn runs, and the output
is a spacetime whose initial data solves a different Lagrangian than the one
evolving it.  The constraint violation is baked in at t=0 and looks like
physics.

The ``.gridinit`` handoff has the mirror-image problem.  GRTeclyn's loader
copies component ``c`` of the file into state variable ``c`` and never reads the
header's ``component_names``.  The Python analysis code indexes the very same
file *by* name.  Both work only while the file's order happens to match
GRTeclyn's slot order, and nothing enforces that.

Synthetic tests below break each contract deliberately; the live-tree tests hold
the real artifacts to it.  See ``research/neuralspacetime/DebugPreGPU.md``.
"""

from __future__ import annotations

import re
import struct
from pathlib import Path

import pytest

from grteclyn_wrapper.core.cross_code import (
    STATE_NAME_ALIASES,
    audit_run_pair,
    parse_state_variable_names,
    read_gridinit_header,
)

REPO = Path(__file__).resolve().parents[3]
CCZ4_HEADER = REPO / "Source" / "CCZ4" / "CCZ4StateVariables.hpp"
EXAMPLE_HEADER = REPO / "Examples" / "RadialRecipe" / "StateVariables.hpp"
BOSON_STAR_PY = (
    REPO
    / "grteclyn-wrapper"
    / "src"
    / "grteclyn_wrapper"
    / "grtresna"
    / "fields"
    / "boson_star.py"
)

_HAVE_CXX = CCZ4_HEADER.exists() and EXAMPLE_HEADER.exists()

# A minimal but faithful component list: the 25 CCZ4 vars followed by the eight
# bicomplex matter components, exactly as the real writer emits them.
GOOD_COMPONENTS = [
    "chi", "h11", "h12", "h13", "h22", "h23", "h33", "K",
    "A11", "A12", "A13", "A22", "A23", "A33", "Theta",
    "Gamma1", "Gamma2", "Gamma3", "lapse", "shift1", "shift2", "shift3",
    "B1", "B2", "B3",
    "phi", "Pi", "phi2", "Pi2",
    "phi_lump1", "Pi_lump1", "phi_lump2", "Pi_lump2",
]


def _write_run(
    tmp_path: Path,
    *,
    tresna: str,
    teclyn: str,
    components: list[str] | None = None,
    declared_ncomp: int | None = None,
) -> Path:
    run = tmp_path / "run"
    (run / "grtresna").mkdir(parents=True)
    (run / "grtresna" / "params.txt").write_text(tresna)
    (run / "params.txt").write_text(teclyn)
    if components is not None:
        _write_gridinit(run / "initial_data.gridinit", components, declared_ncomp)
    return run


def _write_gridinit(
    path: Path, components: list[str], declared_ncomp: int | None = None
) -> None:
    ncomp = declared_ncomp if declared_ncomp is not None else len(components)
    header = (
        "GRTECLYN_GRID_INIT_V2\n"
        f"num_components {ncomp}\n"
        f"component_names {' '.join(components)}\n"
        "nx_ny_nz 2 2 2\n"
        "dx 5.0e-01 5.0e-01 5.0e-01\n"
        "origin 0.0e+00 0.0e+00 0.0e+00\n"
        "END_HEADER\n"
    )
    payload = struct.pack(f"<{8 * len(components)}d", *([0.0] * 8 * len(components)))
    path.write_bytes(header.encode() + payload)


BASE_TRESNA = (
    "num_lumps = 2\n"
    "scalar_mass = 1.0\n"
    "scalar_lambda = 640.0\n"
    "scalar_mu = 85333.0\n"
    "lump0_exotic = 0\n"
    "lump1_exotic = 1\n"
)
BASE_TECLYN = (
    "recipe_scalar_mass = 1.0\n"
    "recipe_scalar_lambda = 640.0\n"
    "recipe_scalar_mu = 85333.0\n"
    "recipe_scalar_field_signs = 1 -1\n"
)


# ---------------------------------------------------------------------------
# Synthetic: each contract, broken on purpose
# ---------------------------------------------------------------------------


def test_agreeing_configs_produce_no_findings(tmp_path: Path) -> None:
    """The control. Without this the tests below prove only that something errors."""
    run = _write_run(tmp_path, tresna=BASE_TRESNA, teclyn=BASE_TECLYN)
    assert not audit_run_pair(run).findings


@pytest.mark.parametrize(
    "constant, tresna_key, teclyn_key",
    [
        ("mass", "scalar_mass", "recipe_scalar_mass"),
        ("lambda", "scalar_lambda", "recipe_scalar_lambda"),
        ("mu", "scalar_mu", "recipe_scalar_mu"),
    ],
)
def test_catches_a_constant_that_differs_between_the_codes(
    tmp_path: Path, constant: str, tresna_key: str, teclyn_key: str
) -> None:
    """Solve with one potential, evolve with another.

    This is Debug.md #14 promoted to the cross-code level: there, a diagnostic
    read mu = 0 while evolution used 85333.  Here the *solver* would use a
    different value from the *evolution*, so the initial data would not satisfy
    the constraints of the system integrating it.
    """
    teclyn = "\n".join(
        f"{teclyn_key} = 999.0" if line.startswith(teclyn_key) else line
        for line in BASE_TECLYN.splitlines()
    )
    run = _write_run(tmp_path, tresna=BASE_TRESNA, teclyn=teclyn)

    blocking = audit_run_pair(run).blocking()
    assert [f.kind for f in blocking] == ["cross-code-constant-mismatch"]
    assert tresna_key in blocking[0].key


def test_catches_a_flipped_field_sign(tmp_path: Path) -> None:
    """``lump1_exotic = 1`` means sign -1; a +1 downstream flips the kinetic term.

    A ghost field solved as a ghost and evolved as ordinary matter (or the
    reverse) is not a small numerical difference -- it is a different theory.
    """
    teclyn = BASE_TECLYN.replace(
        "recipe_scalar_field_signs = 1 -1", "recipe_scalar_field_signs = 1 1"
    )
    run = _write_run(tmp_path, tresna=BASE_TRESNA, teclyn=teclyn)

    blocking = audit_run_pair(run).blocking()
    assert [f.kind for f in blocking] == ["field-sign-mismatch"]
    assert "lump1_exotic" in blocking[0].key


def test_catches_the_truncated_sign_array(tmp_path: Path) -> None:
    """Bug #19 seen from the config side: five lumps upstream, one sign downstream.

    The C++-side check in ``test_param_contract`` catches the reader; this
    catches the emitted pair, which is what a preflight can inspect before a
    campaign burns GPU hours.
    """
    tresna = BASE_TRESNA.replace("num_lumps = 2", "num_lumps = 5")
    teclyn = BASE_TECLYN.replace(
        "recipe_scalar_field_signs = 1 -1", "recipe_scalar_field_signs = 1"
    )
    run = _write_run(tmp_path, tresna=tresna, teclyn=teclyn)

    kinds = [f.kind for f in audit_run_pair(run).blocking()]
    assert "field-signs-truncated" in kinds


@pytest.mark.skipif(not _HAVE_CXX, reason="C++ headers not present")
def test_catches_a_reordered_gridinit(tmp_path: Path) -> None:
    """Swap two components and the loader evolves each as the other.

    GRTeclyn copies by position and never reads component_names, so a reorder
    produces no error at all -- just a spacetime where the lapse is stored in
    the shift.
    """
    swapped = list(GOOD_COMPONENTS)
    swapped[18], swapped[19] = swapped[19], swapped[18]  # lapse <-> shift1
    run = _write_run(
        tmp_path, tresna=BASE_TRESNA, teclyn=BASE_TECLYN, components=swapped
    )

    blocking = audit_run_pair(
        run, ccz4_header=CCZ4_HEADER, example_header=EXAMPLE_HEADER
    ).blocking()
    assert [f.kind for f in blocking] == ["gridinit-slot-mismatch"] * 2


@pytest.mark.skipif(not _HAVE_CXX, reason="C++ headers not present")
def test_accepts_the_documented_slot_aliases(tmp_path: Path) -> None:
    """``phi2``/``Pi2`` and ``phi_lump0``/``Pi_lump0`` are the same two slots.

    StateVariables.hpp:12-13 says so explicitly: the single-complex-scalar model
    reuses the first lump slots.  The writer emits the enum names, the state
    vector holds the lump names, and both are right.  A checker that flagged
    this would be wrong on every file in the tree.
    """
    run = _write_run(
        tmp_path,
        tresna=BASE_TRESNA,
        teclyn=BASE_TECLYN,
        components=GOOD_COMPONENTS,
    )
    report = audit_run_pair(
        run, ccz4_header=CCZ4_HEADER, example_header=EXAMPLE_HEADER
    )
    assert not report.findings, report.format()


def test_catches_a_component_count_that_lies(tmp_path: Path) -> None:
    """``num_components`` sets the payload stride; a wrong one shears the data.

    The loader trusts the count, not the name list.  Off by one and every
    component after the first cell is read from the wrong offset.
    """
    run = _write_run(
        tmp_path,
        tresna=BASE_TRESNA,
        teclyn=BASE_TECLYN,
        components=GOOD_COMPONENTS,
        declared_ncomp=len(GOOD_COMPONENTS) + 1,
    )
    blocking = audit_run_pair(run).blocking()
    assert [f.kind for f in blocking] == ["gridinit-component-count-mismatch"]


def test_header_reader_does_not_load_the_payload(tmp_path: Path) -> None:
    """Guard the guard: these files reach gigabytes.

    A validator that has to read the whole array to check its header is one
    nobody will run before a campaign, which makes it worth nothing.
    """
    path = tmp_path / "big.gridinit"
    _write_gridinit(path, GOOD_COMPONENTS)
    header = read_gridinit_header(path)

    assert header["component_names"] == GOOD_COMPONENTS
    assert header["nx_ny_nz"] == [2, 2, 2]
    assert header["num_components"] == [len(GOOD_COMPONENTS)]


# ---------------------------------------------------------------------------
# Live tree
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAVE_CXX, reason="C++ headers not present")
def test_python_slot_constants_match_the_cxx_enum() -> None:
    """``boson_star.py`` hand-copies GRTeclyn's state layout. Hold it to the copy.

    Its own comment says "must match StateVariables.hpp" -- and a comment is not
    a check.  These indices decide which component the boson-star painter writes
    into; if the C++ enum gains a variable, the Python silently paints the wrong
    slots and the initial data is scrambled with no error anywhere.
    """
    names = parse_state_variable_names(CCZ4_HEADER, EXAMPLE_HEADER)
    source = BOSON_STAR_PY.read_text(encoding="utf-8")

    consts: dict[str, int] = {}
    for line in source.splitlines():
        if line.startswith("_NUM_CCZ4_VARS = "):
            consts["_NUM_CCZ4_VARS"] = int(line.split("=")[1])
        elif line.startswith("_C_") and "_NUM_CCZ4_VARS +" in line:
            name = line.split("=")[0].strip()
            consts[name] = consts["_NUM_CCZ4_VARS"] + int(line.rsplit("+", 1)[1])

    n_ccz4 = names.index("phi")
    assert consts["_NUM_CCZ4_VARS"] == n_ccz4, (
        f"boson_star.py hardcodes _NUM_CCZ4_VARS = {consts['_NUM_CCZ4_VARS']} "
        f"but the C++ headers declare {n_ccz4} CCZ4 variables."
    )

    # Slot -> the name the C++ actually gives it, resolving the documented alias.
    canonical = {v: k for k, v in STATE_NAME_ALIASES.items()}
    for const, expected in [
        ("_C_PHI_P", "phi"),
        ("_C_PI_P", "Pi"),
        ("_C_PHI2_P", "phi2"),
        ("_C_PI2_P", "Pi2"),
        ("_C_PHI_M", "phi_lump1"),
        ("_C_PI_M", "Pi_lump1"),
        ("_C_PHI2_M", "phi_lump2"),
        ("_C_PI2_M", "Pi_lump2"),
    ]:
        slot = consts[const]
        actual = names[slot]
        assert canonical.get(actual, actual) == expected, (
            f"{const} points at slot {slot}, which the C++ calls "
            f"{actual!r}, not {expected!r}."
        )


LOADER_HEADER = REPO / "Examples" / "RadialRecipe" / "ExternalGridInitialData.hpp"


@pytest.mark.skipif(not LOADER_HEADER.exists(), reason="C++ loader not present")
def test_the_loader_resolves_components_by_name() -> None:
    """The gridinit contract must be nominal on both sides, not nominal on one.

    The loader used to copy file component ``c`` into state variable ``c`` and
    never read ``component_names`` at all, while the Python analysis code indexed
    the same file by name.  They agreed only because nobody had reordered a
    column yet.  This test fails if the positional shortcut ever comes back.
    """
    source = LOADER_HEADER.read_text(encoding="utf-8")

    assert "StateVariables::names[v] == name" in source, (
        "ExternalGridInitialData no longer resolves components through "
        "StateVariables::names; the contract is positional again."
    )
    assert "cell(i, j, k, v) = static_cast<amrex::Real>(val);" in source, (
        "the interpolation writes to the file's column index rather than to the "
        "state slot that column was resolved to."
    )
    assert "Refusing to" in source, (
        "an unrecognised component name must abort, not be guessed or skipped."
    )


@pytest.mark.skipif(not LOADER_HEADER.exists(), reason="C++ loader not present")
def test_the_two_alias_tables_agree() -> None:
    """The C++ loader and this package must resolve the same aliases.

    Both sides special-case ``phi2``/``Pi2`` -> ``phi_lump0``/``Pi_lump0``
    (StateVariables.hpp:12-13).  Two hand-maintained copies of one table is how
    the original split happened; this holds them together.
    """
    source = LOADER_HEADER.read_text(encoding="utf-8")
    cxx_aliases = dict(
        re.findall(
            r'if\s*\(name\s*==\s*"([A-Za-z0-9_]+)"\)\s*\n\s*return\s+"([A-Za-z0-9_]+)";',
            source,
        )
    )
    assert cxx_aliases == STATE_NAME_ALIASES, (
        f"C++ loader resolves {cxx_aliases}, this package resolves "
        f"{STATE_NAME_ALIASES}. A slot the two disagree about is loaded into "
        "one variable and read back as another."
    )


def _run_pairs() -> list[Path]:
    runs = REPO / "runs"
    if not runs.exists():
        return []
    return [p.parent.parent for p in sorted(runs.rglob("grtresna/params.txt"))]


@pytest.mark.skipif(not _HAVE_CXX, reason="C++ headers not present")
@pytest.mark.parametrize("run_dir", _run_pairs(), ids=lambda p: p.name)
def test_real_runs_are_cross_code_consistent(run_dir: Path) -> None:
    """Every solved-and-evolved run in the tree must agree with itself."""
    report = audit_run_pair(
        run_dir, ccz4_header=CCZ4_HEADER, example_header=EXAMPLE_HEADER
    )
    assert not report.blocking(), report.format()


@pytest.mark.skipif(not _HAVE_CXX, reason="C++ headers not present")
def test_state_variable_parse_is_not_silently_empty() -> None:
    """Guard the guard: an empty name list would make every check above pass."""
    names = parse_state_variable_names(CCZ4_HEADER, EXAMPLE_HEADER)
    assert len(names) > 30
    assert names[0] == "chi"
    assert "lapse" in names and "phi" in names
