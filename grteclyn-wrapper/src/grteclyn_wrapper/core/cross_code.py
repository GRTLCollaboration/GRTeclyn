"""Cross-code consistency: what goes to GRTresna vs what goes to GRTeclyn.

WHAT THIS IS FOR
----------------
A run is configured twice.  The wrapper writes ``grtresna/params.txt`` for the
elliptic solve and ``params.txt`` for the evolution, and the two files use
**different names for the same physics**::

    GRTresna              GRTeclyn
    --------              --------
    scalar_mass      <->  recipe_scalar_mass
    scalar_lambda    <->  recipe_scalar_lambda
    scalar_mu        <->  recipe_scalar_mu
    lump<k>_exotic   <->  recipe_scalar_field_signs[k]     (sign = 1 - 2*exotic)

Nothing checks that they agree.  If they ever disagree, the constraint solve
produces initial data for one Lagrangian and the evolution integrates it under
another: the data is no longer a solution of the system evolving it, and the
constraint violation is there at t=0 by construction.  Nothing in either code
would report an error -- GRTresna would converge, GRTeclyn would run, and the
resulting spacetime would be a plausible-looking artifact.

The third file in the handoff is the ``.gridinit``, and it carries its own trap.
Its header lists ``component_names``, but GRTeclyn's loader never reads them:
``ExternalGridInitialData::compute`` copies component ``c`` of the file into
state variable ``c``, for ``c < min(num_components, NUM_VARS)``.  The contract
is **positional**.  The Python side, meanwhile, indexes the same file **by
name** (``comp_names.index("lapse")`` and friends).  Two codes, two different
notions of what identifies a component, no cross-check.  Reorder or insert a
state variable and every field lands one slot off, silently.

See ``research/neuralspacetime/DebugPreGPU.md``.
"""

from __future__ import annotations

import re
from pathlib import Path

from .param_contract import ContractReport, Finding, parse_params_file

__all__ = [
    "PHYSICS_CORRESPONDENCE",
    "STATE_NAME_ALIASES",
    "parse_state_variable_names",
    "read_gridinit_header",
    "audit_run_pair",
]

#: (GRTresna key, GRTeclyn key) pairs that must carry the same number.
PHYSICS_CORRESPONDENCE: tuple[tuple[str, str], ...] = (
    ("scalar_mass", "recipe_scalar_mass"),
    ("scalar_lambda", "recipe_scalar_lambda"),
    ("scalar_mu", "recipe_scalar_mu"),
    ("G_Newton", "G_Newton"),
)

#: GRTeclyn declares two names for the same slot.  ``StateVariables.hpp:12-13``:
#: the single-complex-scalar model reuses the first lump slots, so ``c_phi2``
#: and ``phi_lump0`` are the same component.  The wrapper writes the enum name,
#: ``StateVariables::names`` holds the lump name.  Both are correct.
STATE_NAME_ALIASES: dict[str, str] = {
    "phi2": "phi_lump0",
    "Pi2": "Pi_lump0",
}

_CCZ4_NAMES_RE = re.compile(
    r"namespace\s+CCZ4StateVariables.*?names\s*=\s*\{(?P<body>.*?)\}", re.DOTALL
)
_ADDITIONAL_NAMES_RE = re.compile(
    r"additional_names\s*=\s*\{(?P<body>.*?)\}", re.DOTALL
)
_QUOTED = re.compile(r'"([^"]+)"')


def parse_state_variable_names(
    ccz4_header: Path, example_header: Path
) -> list[str]:
    """The authoritative component order, straight from the C++ headers.

    Parsed rather than duplicated in Python on purpose: a hand-kept copy is the
    thing that goes stale.  ``grtresna/fields/boson_star.py`` currently holds
    such a copy, marked only by the comment "must match StateVariables.hpp".
    """
    ccz4 = _CCZ4_NAMES_RE.search(
        ccz4_header.read_text(encoding="utf-8", errors="replace")
    )
    extra = _ADDITIONAL_NAMES_RE.search(
        example_header.read_text(encoding="utf-8", errors="replace")
    )
    if not ccz4 or not extra:
        raise ValueError("could not parse state variable names from C++ headers")
    return _QUOTED.findall(ccz4.group("body")) + _QUOTED.findall(
        extra.group("body")
    )


def read_gridinit_header(path: Path) -> dict[str, object]:
    """Parse a ``.gridinit`` header without touching the binary payload.

    These files run to gigabytes; a validator that has to load one in order to
    check its header would never be run.
    """
    out: dict[str, object] = {}
    with path.open("rb") as fin:
        for _ in range(64):  # headers are ~8 lines; bail out rather than scan
            raw = fin.readline()
            if not raw:
                break
            line = raw.decode("utf-8", errors="replace").strip()
            if line == "END_HEADER":
                break
            parts = line.split()
            if not parts:
                continue
            key, vals = parts[0], parts[1:]
            if key == "component_names":
                out[key] = vals
            elif key in {"nx_ny_nz", "num_components"}:
                out[key] = [int(v) for v in vals]
            elif key in {"dx", "origin"}:
                out[key] = [float(v) for v in vals]
    return out


def _as_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    try:
        return float(raw.strip().strip('"'))
    except ValueError:
        return None


def _check_physics_constants(
    tresna: dict[str, str], teclyn: dict[str, str]
) -> list[Finding]:
    findings: list[Finding] = []
    for t_key, e_key in PHYSICS_CORRESPONDENCE:
        t_val, e_val = _as_float(tresna.get(t_key)), _as_float(teclyn.get(e_key))
        if t_val is None or e_val is None:
            continue  # not every model configures every constant
        if t_val != e_val:
            findings.append(
                Finding(
                    "CRITICAL",
                    "cross-code-constant-mismatch",
                    f"{t_key}/{e_key}",
                    f"GRTresna solved with {t_key} = {t_val:g} but GRTeclyn "
                    f"evolves with {e_key} = {e_val:g}. The initial data is not "
                    "a solution of the system integrating it; the constraint "
                    "violation is present at t=0 by construction.",
                )
            )
    return findings


def _check_field_signs(
    tresna: dict[str, str], teclyn: dict[str, str]
) -> list[Finding]:
    """``lump<k>_exotic`` and ``recipe_scalar_field_signs[k]`` must agree.

    GRTresna flags a ghost lump with ``exotic = 1``; GRTeclyn carries the same
    information as a sign of -1.  The relation is ``sign = 1 - 2*exotic``.  Bug
    #19 lived in exactly this handoff.
    """
    raw_signs = teclyn.get("recipe_scalar_field_signs")
    if raw_signs is None:
        return []
    try:
        signs = [float(tok) for tok in raw_signs.split()]
    except ValueError:
        return [
            Finding(
                "MAJOR",
                "unparsable-field-signs",
                "recipe_scalar_field_signs",
                f"could not parse {raw_signs!r} as a list of signs.",
            )
        ]

    findings: list[Finding] = []
    if len(signs) == 1 and _as_float(tresna.get("num_lumps")) not in (None, 1.0):
        # The bug-#19 signature seen from the other side: five lumps upstream,
        # one sign token downstream.
        findings.append(
            Finding(
                "CRITICAL",
                "field-signs-truncated",
                "recipe_scalar_field_signs",
                f"GRTresna configures {tresna.get('num_lumps')} lumps but only "
                "one sign token survived into the GRTeclyn config. This is the "
                "shape of Debug.md #19, where the sign array was parsed into a "
                "std::string and every token after the first was dropped.",
            )
        )

    for k, sign in enumerate(signs):
        exotic = _as_float(tresna.get(f"lump{k}_exotic"))
        if exotic is None:
            continue
        expected = 1.0 - 2.0 * round(exotic)
        if sign != expected:
            findings.append(
                Finding(
                    "CRITICAL",
                    "field-sign-mismatch",
                    f"lump{k}_exotic/recipe_scalar_field_signs[{k}]",
                    f"GRTresna painted lump {k} with exotic = {exotic:g} "
                    f"(sign {expected:+g}) but GRTeclyn evolves it with sign "
                    f"{sign:+g}. The kinetic term flips between the solve and "
                    "the evolution.",
                )
            )
    return findings


def _check_gridinit_layout(
    gridinit: Path, cxx_names: list[str] | None
) -> list[Finding]:
    header = read_gridinit_header(gridinit)
    names = header.get("component_names")
    ncomp_field = header.get("num_components")
    findings: list[Finding] = []

    if not isinstance(names, list) or not names:
        return [
            Finding(
                "MAJOR",
                "gridinit-header-unreadable",
                gridinit.name,
                "no component_names in the header; the positional contract "
                "cannot be checked at all.",
            )
        ]

    ncomp = ncomp_field[0] if isinstance(ncomp_field, list) and ncomp_field else None
    if ncomp is not None and ncomp != len(names):
        findings.append(
            Finding(
                "CRITICAL",
                "gridinit-component-count-mismatch",
                gridinit.name,
                f"header declares num_components = {ncomp} but lists "
                f"{len(names)} component names. GRTeclyn trusts the count and "
                "reads the payload with that stride, so every component after "
                "the first row would be misaligned.",
            )
        )

    if cxx_names is None:
        return findings

    for i, name in enumerate(names):
        canonical = STATE_NAME_ALIASES.get(name, name)
        if i >= len(cxx_names):
            findings.append(
                Finding(
                    "MAJOR",
                    "gridinit-longer-than-state-vector",
                    gridinit.name,
                    f"component {i} ({name!r}) has no slot in GRTeclyn's state "
                    f"vector ({len(cxx_names)} vars). The loader truncates at "
                    "min(num_components, NUM_VARS) without saying so.",
                )
            )
            break
        if canonical != cxx_names[i]:
            findings.append(
                Finding(
                    "CRITICAL",
                    "gridinit-slot-mismatch",
                    gridinit.name,
                    f"component {i} is {name!r} in the file but slot {i} is "
                    f"{cxx_names[i]!r} in GRTeclyn. The loader matches by "
                    "position and never reads component_names, so this field "
                    "would be evolved as the wrong variable, silently.",
                )
            )
    return findings


def audit_run_pair(
    run_dir: Path,
    *,
    ccz4_header: Path | None = None,
    example_header: Path | None = None,
) -> ContractReport:
    """Cross-check one run's GRTresna config, GRTeclyn config, and gridinit.

    ``run_dir`` is the evolution directory: it holds ``params.txt``, the
    ``grtresna/`` subdirectory with the solver's own params.txt, and the
    ``.gridinit`` handed between them.
    """
    teclyn_path = run_dir / "params.txt"
    tresna_path = run_dir / "grtresna" / "params.txt"
    report = ContractReport(params_file=str(run_dir))

    if not teclyn_path.exists() or not tresna_path.exists():
        report.findings.append(
            Finding(
                "MINOR",
                "incomplete-run-pair",
                run_dir.name,
                "missing params.txt on one side; nothing to cross-check.",
            )
        )
        return report

    teclyn = {e.key: e.raw_value for e in parse_params_file(teclyn_path)}
    tresna = {e.key: e.raw_value for e in parse_params_file(tresna_path)}
    report.n_emitted = len(teclyn) + len(tresna)

    report.findings.extend(_check_physics_constants(tresna, teclyn))
    report.findings.extend(_check_field_signs(tresna, teclyn))

    cxx_names = None
    if ccz4_header and example_header:
        try:
            cxx_names = parse_state_variable_names(ccz4_header, example_header)
        except (OSError, ValueError):
            cxx_names = None

    for gridinit in sorted(run_dir.glob("*.gridinit")):
        report.findings.extend(_check_gridinit_layout(gridinit, cxx_names))

    return report
