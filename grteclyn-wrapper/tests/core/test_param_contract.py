"""The params.txt boundary contract: writer vs reader.

WHY THIS FILE EXISTS
--------------------
``params.txt`` is the only interface between this wrapper and the two C++ codes.
It is untyped and unvalidated: a key nobody reads, or a five-token vector read
into a ``std::string``, produces no error -- it produces a plausible wrong
answer.  Both sides are individually self-consistent, so neither side's own unit
tests can see the mismatch.  Only a cross-check can.

The tests below come in two halves:

1. **Synthetic** -- a hand-built C++ file and params.txt reproducing each known
   bug shape.  These pin the checker itself.  Every one of them is a bug that
   actually shipped, so if the checker cannot catch it here it is worthless in
   the tree.
2. **Live tree** -- run the checker against the real GRTeclyn/GRTresna sources
   and real emitted params.txt files.  These are the regression guard.

Debug.md #19 (recipe_scalar_field_signs), #14 (recipe_scalar_mu = 0), and the
retracted ``trajectory_lump<k>_exotic`` finding are the three shapes.
See ``research/neuralspacetime/DebugPreGPU.md``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from grteclyn_wrapper.core.param_contract import (
    ARITY_ARRAY,
    ARITY_SCALAR,
    audit_params_file,
    parse_params_file,
    scan_cxx,
)

REPO = Path(__file__).resolve().parents[3]
GRTECLYN_SRC = REPO / "Source"
GRTECLYN_EX = REPO / "Examples"
GRTRESNA = REPO.parent / "GRTresna"
WRAPPER_SRC = REPO / "grteclyn-wrapper" / "src"


def _cxx_roots() -> list[Path]:
    return [
        GRTECLYN_SRC,
        GRTECLYN_EX,
        GRTRESNA / "Source",
        GRTRESNA / "Examples",
    ]


# ---------------------------------------------------------------------------
# Synthetic: each known bug shape must be caught
# ---------------------------------------------------------------------------


def test_catches_vector_read_into_a_string(tmp_path: Path) -> None:
    """Bug #19 verbatim: five signs written, one std::string reading them.

    ParmParse's scalar query returns token 0 and drops the rest, so
    ``1 -1 -1 1 -1`` became ``[1, 0, 0, 0, 0]``.  Every bicomplex campaign ran
    with zero pump force on the phantom sector.
    """
    cxx = tmp_path / "Broken.cpp"
    cxx.write_text(
        "void read(GRParmParse &pp) {\n"
        "    std::string recipe_scalar_field_signs_str;\n"
        '    pp.load("recipe_scalar_field_signs", '
        "recipe_scalar_field_signs_str);\n"
        "}\n"
    )
    params = tmp_path / "params.txt"
    params.write_text("recipe_scalar_field_signs = 1 -1 -1 1 -1\n")

    report = audit_params_file(params, [cxx], require_physics_keys=False)
    blocking = report.blocking()

    assert [f.key for f in blocking] == ["recipe_scalar_field_signs"]
    assert blocking[0].kind == "vector-read-as-scalar"
    assert "5 tokens" in blocking[0].detail


def test_accepts_a_vector_read_into_a_container(tmp_path: Path) -> None:
    """The counterpart: the same key read correctly must NOT be flagged.

    ``GRParmParse::load`` is overloaded on the destination type, so the function
    name is identical in the broken and fixed versions -- only the variable's
    type differs.  A checker that keys off the function name reports both, and a
    checker that cries wolf gets muted.
    """
    cxx = tmp_path / "Fixed.cpp"
    cxx.write_text(
        "void read(GRParmParse &pp) {\n"
        "    std::array<double, 5> recipe_scalar_field_signs;\n"
        '    pp.load("recipe_scalar_field_signs", recipe_scalar_field_signs);\n'
        "}\n"
    )
    params = tmp_path / "params.txt"
    params.write_text("recipe_scalar_field_signs = 1 -1 -1 1 -1\n")

    report = audit_params_file(params, [cxx], require_physics_keys=False)
    assert not report.blocking(), report.format()


def test_catches_a_missing_physics_constant(tmp_path: Path) -> None:
    """Bug #14: mu read with a silent default of 0 while evolution used 85333.

    The wrapper stopping emitting the key is indistinguishable, at the C++ side,
    from the wrapper asking for zero.  That is what fabricated the
    energy-condition violations.
    """
    cxx = tmp_path / "Matter.cpp"
    cxx.write_text(
        "void read(GRParmParse &pp) {\n"
        "    double recipe_scalar_mu;\n"
        '    pp.load("recipe_scalar_mu", recipe_scalar_mu, 0.0);\n'
        "}\n"
    )
    params = tmp_path / "params.txt"
    params.write_text("recipe_scalar_mass = 0.5\n")  # mu conspicuously absent

    report = audit_params_file(params, [cxx])
    blocking = report.blocking()

    assert [f.key for f in blocking] == ["recipe_scalar_mu"]
    assert blocking[0].kind == "missing-physics-constant"


def test_catches_a_dead_key(tmp_path: Path) -> None:
    """A knob written into every params.txt and read by nobody.

    Harmless to the physics, corrosive to the provenance record: the config
    claims a lever exists, and anyone reading it later believes the campaign
    swept that lever.
    """
    cxx = tmp_path / "Real.cpp"
    cxx.write_text(
        'void read(GRParmParse &pp) { double x; pp.load("real_key", x, 0.0); }\n'
    )
    params = tmp_path / "params.txt"
    params.write_text("real_key = 1.0\nphantom_knob = 7\n")

    report = audit_params_file(params, [cxx], require_physics_keys=False)
    dead = [f for f in report.findings if f.kind == "dead-key"]

    assert [f.key for f in dead] == ["phantom_knob"]


def test_a_python_consumed_key_is_not_dead(tmp_path: Path) -> None:
    """``trajectory_retrograde_only`` is read by the wrapper, not by C++.

    Emitted for provenance and consumed by the Python clamp.  Only a scan that
    looks at both sides can tell it apart from a genuinely dead key.
    """
    cxx = tmp_path / "Empty.cpp"
    cxx.write_text("// no reads here\n")
    py = tmp_path / "clamp.py"
    py.write_text('enabled = bool(int(overrides.get("trajectory_retrograde_only", 0)))\n')
    params = tmp_path / "params.txt"
    params.write_text("trajectory_retrograde_only = 1\n")

    report = audit_params_file(
        params, [cxx], python_roots=[py], require_physics_keys=False
    )
    assert not [f for f in report.findings if f.kind == "dead-key"], report.format()


def test_unresolvable_arity_warns_but_does_not_accuse(tmp_path: Path) -> None:
    """An undeclared target is a hole in the checker, not evidence of a bug.

    Reporting it as CRITICAL would be an accusation the scan cannot support.
    It is still reported, because a silently skipped check is worse than a
    noisy one -- just not at a severity that blocks a build.
    """
    cxx = tmp_path / "Opaque.cpp"
    cxx.write_text('void read(GRParmParse &pp) { pp.load("mystery", some_alias); }\n')
    params = tmp_path / "params.txt"
    params.write_text("mystery = 1 2 3\n")

    report = audit_params_file(params, [cxx], require_physics_keys=False)

    assert not report.blocking()
    assert [f.kind for f in report.findings] == ["unverified-arity"]


# ---------------------------------------------------------------------------
# Static-scan blind spots that once produced false dead-key reports
# ---------------------------------------------------------------------------


def test_catches_indexed_keys_past_the_slot_cap(tmp_path: Path) -> None:
    """Ask for eight lumps, get five, hear nothing.

    ``trajectory_num_lumps`` is clamped to GRTRESNA_MAX_INDEPENDENT_SCALARS on
    load and the per-lump read loop runs to the clamped count, so
    ``trajectory_lump5_*`` and beyond are emitted and read by nobody.  The
    dynamic-prefix exemption used to hide this: "assume a reader exists" is the
    same mistake as a silent default -- it turns an unknown into a reassurance.
    """
    cxx = tmp_path / "Layout.hpp"
    cxx.write_text("static constexpr int GRTRESNA_MAX_INDEPENDENT_SCALARS = 5;\n")
    params = tmp_path / "params.txt"
    params.write_text(
        "trajectory_num_lumps = 8\n"
        "trajectory_lump4_R0 = 5.0\n"
        "trajectory_lump5_R0 = 5.4\n"
        "trajectory_lump7_omega_rot = -0.3\n"
    )

    report = audit_params_file(params, [cxx], require_physics_keys=False)
    keys = [f.key for f in report.blocking()]

    assert keys == ["trajectory_lump5_R0", "trajectory_lump7_omega_rot"]
    assert "trajectory_lump4_R0" not in keys  # index 4 is within the cap


def test_grtresna_lump_keys_are_not_capped(tmp_path: Path) -> None:
    """GRTresna's own ``lump<k>_*`` honours num_lumps -- only GRTeclyn clamps.

    The eight-lump initial data is genuinely solved with eight lumps.  Flagging
    those keys too would blame the wrong side of the boundary.
    """
    cxx = tmp_path / "Layout.hpp"
    cxx.write_text("static constexpr int GRTRESNA_MAX_INDEPENDENT_SCALARS = 5;\n")
    params = tmp_path / "params.txt"
    params.write_text("num_lumps = 8\nlump7_amp = 0.075\n")

    report = audit_params_file(params, [cxx], require_physics_keys=False)
    assert not report.blocking(), report.format()


def test_helper_built_keys_are_not_reported_dead(tmp_path: Path) -> None:
    """``load_coeff_array(pp, "recipe_chi_coeff", coeffs)`` builds key_0..key_N.

    Twenty-four live keys were reported dead before the scan understood that a
    literal handed to a helper alongside the ParmParse object may be a prefix.
    """
    cxx = tmp_path / "Recipe.hpp"
    cxx.write_text(
        "void load_coeff_array(GRParmParse &pp, const char *prefix,\n"
        "                      std::array<double, 8> &coeffs) {\n"
        "    for (int n = 0; n < num_bases; ++n) {\n"
        "        std::ostringstream key;\n"
        '        key << prefix << "_" << n;\n'
        "        pp.load(key.str().c_str(), coeffs[n], 0.0);\n"
        "    }\n"
        "}\n"
        "void read(GRParmParse &pp) {\n"
        '    load_coeff_array(pp, "recipe_chi_coeff", recipe_params.chi_coeffs);\n'
        "}\n"
    )
    params = tmp_path / "params.txt"
    params.write_text("recipe_chi_coeff_0 = 1.0\nrecipe_chi_coeff_1 = 2.0\n")

    report = audit_params_file(params, [cxx], require_physics_keys=False)
    assert not [f for f in report.findings if f.kind == "dead-key"], report.format()


def test_constructor_style_declarations_resolve(tmp_path: Path) -> None:
    """``std::vector<int> v(2 * n);`` is a declaration, not a call.

    SimulationParametersBase.hpp:128 declares the extraction-mode buffer this
    way; missing it left a real array read looking unverifiable.
    """
    cxx = tmp_path / "Extraction.hpp"
    cxx.write_text(
        "void read(GRParmParse &pp) {\n"
        "    std::vector<int> extraction_modes_vect(2 * num_modes);\n"
        '    pp.load("modes", extraction_modes_vect, 2 * num_modes);\n'
        "}\n"
    )
    index = scan_cxx([cxx])
    assert index.reads["modes"][0].arity == ARITY_ARRAY


def test_indexed_target_is_scalar(tmp_path: Path) -> None:
    """``pp.load("k", arr[0])`` absorbs exactly one token, container or not."""
    cxx = tmp_path / "Indexed.cpp"
    cxx.write_text(
        "void read(GRParmParse &pp) {\n"
        "    std::array<double, 3> coeffs;\n"
        '    pp.load("k", coeffs[0], 0.0);\n'
        "}\n"
    )
    index = scan_cxx([cxx])
    assert index.reads["k"][0].arity == ARITY_SCALAR


# ---------------------------------------------------------------------------
# Live tree: the actual regression guard
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not GRTECLYN_SRC.exists(), reason="C++ sources not present")
def test_live_sources_still_read_the_sign_array_as_an_array() -> None:
    """The #19 fix must stay fixed.

    Pinned as its own test rather than folded into the whole-file audit, so a
    regression here names itself instead of arriving as one line in a report.
    """
    index = scan_cxx([GRTECLYN_SRC, GRTECLYN_EX])
    sites = index.reads.get("recipe_scalar_field_signs", [])

    assert sites, "no reader found for recipe_scalar_field_signs at all"
    assert any(s.is_array_read for s in sites), (
        "every reader of recipe_scalar_field_signs is scalar again -- "
        f"{[(Path(s.file).name, s.line, s.fn, s.var, s.arity) for s in sites]}"
    )


# The eight-lump trajectory campaign is a known-defective artifact: it emits 24
# keys for lumps 5-7 that GRTeclyn cannot reach (DebugPreGPU.md PG-9).  It is
# excluded from the sweep below and pinned by its own test instead, so the
# defect stays visible as a fact about that campaign rather than as recurring
# noise that eventually gets the whole sweep muted.
KNOWN_INERT_CAMPAIGN = "qball_traj_bicomplex_8lump_v1"


def _emitted_params_files(limit: int = 12) -> list[Path]:
    runs = REPO / "runs"
    if not runs.exists():
        return []
    found = [
        p
        for p in sorted(runs.rglob("params.txt"))
        if KNOWN_INERT_CAMPAIGN not in p.parts
    ]
    return found[:limit]


@pytest.mark.skipif(not GRTECLYN_SRC.exists(), reason="C++ sources not present")
@pytest.mark.parametrize("params", _emitted_params_files(), ids=lambda p: p.parent.name)
def test_emitted_params_files_have_no_blocking_findings(params: Path) -> None:
    """No real params.txt may carry a CRITICAL finding.

    This is the check that would have caught #19 and #14 the day they were
    introduced, on the first campaign that emitted a params.txt.
    """
    report = audit_params_file(
        params, _cxx_roots(), python_roots=[WRAPPER_SRC], require_physics_keys=False
    )
    assert not report.blocking(), report.format()


@pytest.mark.skipif(not GRTECLYN_SRC.exists(), reason="C++ sources not present")
def test_eight_lump_campaign_dimensions_are_inert() -> None:
    """Pin PG-9 on the campaign that has it.

    The QD search declared 8 lumps x 8 params and GRTeclyn evolved 5 lumps, so
    24 emitted keys -- and every search dimension wired to them -- did nothing.
    If this ever stops failing, either the cap was raised or the campaign was
    re-emitted; both are worth noticing.
    """
    hits = sorted(
        (REPO / "runs").rglob(f"{KNOWN_INERT_CAMPAIGN}/*/params.txt")
    )
    if not hits:
        pytest.skip("eight-lump campaign artifacts not present")

    report = audit_params_file(
        hits[0], _cxx_roots(), python_roots=[WRAPPER_SRC], require_physics_keys=False
    )
    inert = [f for f in report.blocking() if f.kind == "indexed-key-beyond-cap"]

    assert len(inert) == 24, report.format()
    assert {f.key.split("_")[1] for f in inert} == {"lump5", "lump6", "lump7"}


@pytest.mark.skipif(not GRTECLYN_SRC.exists(), reason="C++ sources not present")
def test_an_over_cap_campaign_cannot_be_launched() -> None:
    """PG-9's fix: refuse at search-space construction, not after the GPU time.

    The 8-lump campaign got as far as writing 200 evaluations' worth of archive
    before anyone noticed 21 of its dimensions were inert. Both trajectory space
    builders now read the cap out of the C++ header and refuse.
    """
    from grteclyn_wrapper.core.param_contract import scalar_slot_cap
    from grteclyn_wrapper.search.optimize.spaces import (
        grtresna_trajectory_boson_search_space,
        grtresna_trajectory_search_space,
    )

    cap = scalar_slot_cap()
    assert cap is not None, "the cap must be read from the C++, never guessed"

    for build in (grtresna_trajectory_search_space,
                  grtresna_trajectory_boson_search_space):
        build(cap)  # at the cap is fine
        with pytest.raises(ValueError, match="can drive at most"):
            build(cap + 1)


@pytest.mark.skipif(not GRTECLYN_SRC.exists(), reason="C++ sources not present")
def test_the_lump_count_clamp_is_gone_from_the_cxx() -> None:
    """The evolution must refuse an over-cap count, not truncate it quietly.

    A silent clamp is what made PG-9 invisible: the config said 8, the run did 5,
    and nothing anywhere said so.
    """
    src = (REPO / "Examples" / "RadialRecipe" /
           "SimulationParameters.hpp").read_text(encoding="utf-8")

    assert "n_traj = GRTRESNA_MAX_INDEPENDENT_SCALARS;" not in src, (
        "trajectory_num_lumps is being clamped to the cap again instead of "
        "aborting; over-cap lumps would silently stop being driven."
    )
    assert "rl_num_lumps = GRTRESNA_MAX_INDEPENDENT_SCALARS;" not in src, (
        "rl_num_lumps is being clamped to the cap again instead of aborting."
    )
    assert src.count("amrex::Abort(") >= 2


@pytest.mark.skipif(not GRTECLYN_SRC.exists(), reason="C++ sources not present")
def test_params_parsing_is_not_silently_empty() -> None:
    """Guard the guard: an audit over zero keys passes vacuously.

    Every check above is 'no finding of kind X'.  If the parser regressed to
    returning nothing, all of them would go green while checking nothing.
    """
    files = _emitted_params_files(limit=1)
    if not files:
        pytest.skip("no emitted params.txt in runs/")
    assert len(parse_params_file(files[0])) > 20
    assert len(scan_cxx(_cxx_roots()).reads) > 50
