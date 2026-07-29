"""Static cross-check of the params.txt boundary against its C++ readers.

WHAT THIS IS FOR
----------------
``params.txt`` is the only interface between the Python wrapper and the two C++
codes (GRTresna's elliptic solver, GRTeclyn's evolution).  It is untyped: a
misspelled key, a vector written where a scalar is read, or a key nobody reads
produces **no error at all** -- it produces a plausible-looking wrong answer.
Three real bugs of exactly this class are on record:

* ``recipe_scalar_field_signs`` -- five tokens written, read into a single
  ``std::string``, so only token 0 parsed.  ``[1,-1,-1,1,-1]`` became
  ``[1,0,0,0,0]`` and the phantom sector received zero pump force in every
  bicomplex campaign.  (Debug.md bug #19)
* ``recipe_scalar_mu`` -- read with a silent default of 0.0 by a diagnostic
  while the evolution used 85333, fabricating energy-condition violations that
  did not exist.  (Debug.md bug #14)
* ``trajectory_lump<k>_exotic`` -- written into every params.txt, read by
  nothing, making the provenance record claim a knob was live when it was not.
  (Debug.md, abandoned/retracted)

None of the three was catchable by a unit test on either side alone, because
each side is individually self-consistent.  They are only visible by comparing
the writer against the reader.  This module does that comparison statically:
parse the C++ readers, parse a real emitted params.txt, and diff.

USAGE
-----
As a test (see ``tests/core/test_param_contract.py``), and as a preflight check
before a campaign::

    from grteclyn_wrapper.core.param_contract import audit_params_file
    report = audit_params_file(Path("params.txt"), cxx_roots=[...])
    assert not report.blocking(), report.format()

DELIBERATE LIMITATION
---------------------
This is a *static* check.  It proves a key is read somewhere with a compatible
arity; it cannot prove the reader interprets the value correctly (units, sign
conventions, normalization).  For that, the runtime cross-code comparison in
``grtresna/matter/profile_contract.py`` is the right instrument -- it paints the
Python reference field and compares against what GRTresna actually produced.
The two are complementary: this one catches plumbing, that one catches physics.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

__all__ = [
    "CxxRead",
    "EmittedKey",
    "Finding",
    "ContractReport",
    "parse_cxx_reads",
    "parse_params_file",
    "audit_params_file",
    "scalar_slot_cap",
]

# pp.load("key", var, default)      -> scalar or std::array, with a default
# pp.get("key", var)                -> scalar, required
# pp.query("key", var)              -> scalar, optional
# pp.getarr("key", vec, start, n)   -> vector, required
# pp.queryarr("key", vec, start, n) -> vector, optional
# pp.countval("key")                -> arity probe; presence means the author
#                                      thought about multi-token input
_READ_RE = re.compile(
    r"""\bpp\w*\.(?P<fn>load|get|query|getarr|queryarr|countval)\s*\(\s*
        (?P<key>"(?P<lit>[A-Za-z0-9_.]+)"|[A-Za-z_]\w*)
        (?:\s*,\s*(?P<var>[A-Za-z_][\w:.\[\]]*))?""",
    re.VERBOSE,
)

_ARRAY_FNS = {"getarr", "queryarr", "countval"}
_PARMPARSE_FNS = {"load", "get", "query", "getarr", "queryarr", "countval"}

# Helper functions that take the ParmParse object and a key literal and do the
# reading themselves:
#
#   load_coeff_array(pp, "recipe_chi_coeff", recipe_params.chi_coeffs)
#   StateVariablesParmParse::load_vars_to_vector(pp, "nonzero_asymptotic_vars", v)
#
# Without this the scan calls 25 live keys dead.  The literal may be a whole key
# or a *prefix* the helper suffixes with an index -- we cannot tell statically,
# so it is registered as both and the audit accepts either match.
_HELPER_READ_RE = re.compile(
    r"""\b(?P<fn>[A-Za-z_][\w:]*)\s*\(\s*
        (?P<pp>[A-Za-z_]\w*)\s*,\s*
        "(?P<lit>[A-Za-z0-9_.]+)"
        (?:\s*,\s*(?P<var>[A-Za-z_][\w:.\[\]]*))?""",
    re.VERBOSE,
)

# CRITICAL: arity is decided by the TYPE OF THE TARGET VARIABLE, not by the
# function name.  ``GRParmParse::load`` is overloaded -- it forwards to the
# array ``get``/``getarr`` when handed a std::array or std::vector, and to the
# scalar ``get`` otherwise (Source/IO/GRParmParse.hpp:27-72).  A checker that
# keys off the function name alone reports every `pp.load("center", center)` as
# a truncation bug, which is wrong and trains people to ignore it.
#
# Bug #19 is precisely a *type* error: `recipe_scalar_field_signs` was declared
# std::string, so the overload resolved to the scalar branch and ParmParse
# returned token 0.  So the type is the discriminator, and resolving it is the
# whole job of this checker.
_ARRAY_TYPE_RE = re.compile(
    r"\b(?:std::array|std::vector|amrex::Vector|Vector|std::pair)\s*<"
)
_SCALAR_TYPE_RE = re.compile(
    r"\b(?:double|float|int|bool|unsigned|size_t|std::string|string|"
    r"amrex::Real|Real)\b"
)

ARITY_ARRAY = "array"
ARITY_SCALAR = "scalar"
ARITY_UNKNOWN = "unknown"

# A params.txt assignment: `key = v1 v2 v3` or `key = "string"`.
_ASSIGN_RE = re.compile(r"^\s*([A-Za-z_][\w.]*)\s*=\s*(.*?)\s*$")

# Keys whose value is legitimately a free-text string that may contain spaces;
# token-counting them would produce false positives.
_STRING_VALUED = {
    "recipe_matter_model",
    "recipe_initial_data_file",
    "amr.plot_file",
    "amr.check_file",
    "output_path",
}


@dataclass(frozen=True)
class CxxRead:
    """One ParmParse read site found in the C++ sources."""

    key: str
    fn: str
    file: str
    line: int
    var: str | None = None
    arity: str = ARITY_UNKNOWN

    @property
    def is_array_read(self) -> bool:
        """True when this site can consume more than one token.

        Either the function is inherently array-valued, or the destination
        variable is a container type.
        """
        return self.fn in _ARRAY_FNS or self.arity == ARITY_ARRAY

    @property
    def is_proven_scalar(self) -> bool:
        """True only when the destination type was resolved AND is scalar.

        Unresolved types are never reported as bugs -- an unverifiable site is
        a gap in the checker, not evidence against the code.
        """
        return self.fn not in _ARRAY_FNS and self.arity == ARITY_SCALAR

    @property
    def has_silent_default(self) -> bool:
        """``pp.load``/``query`` supply a default; a missing key is not an error."""
        return self.fn in {"load", "query", "queryarr"}


@dataclass(frozen=True)
class EmittedKey:
    """One assignment found in an emitted params.txt."""

    key: str
    raw_value: str
    line: int

    @property
    def n_tokens(self) -> int:
        if self.key in _STRING_VALUED or self.raw_value.startswith('"'):
            return 1
        return len(self.raw_value.split())

    @property
    def is_vector(self) -> bool:
        return self.n_tokens > 1


@dataclass(frozen=True)
class Finding:
    severity: str  # "CRITICAL" | "MAJOR" | "MINOR"
    kind: str
    key: str
    detail: str

    def __str__(self) -> str:
        return f"[{self.severity}] {self.kind}: {self.key} -- {self.detail}"


@dataclass
class ContractReport:
    params_file: str
    findings: list[Finding] = field(default_factory=list)
    n_emitted: int = 0
    n_cxx_keys: int = 0

    def blocking(self) -> list[Finding]:
        """Findings that must fail a build, not merely warn."""
        return [f for f in self.findings if f.severity == "CRITICAL"]

    def by_severity(self, severity: str) -> list[Finding]:
        return [f for f in self.findings if f.severity == severity]

    def format(self) -> str:
        if not self.findings:
            return (
                f"params contract OK: {self.n_emitted} keys emitted, "
                f"{self.n_cxx_keys} keys read by C++, no mismatches."
            )
        lines = [
            f"params contract: {len(self.findings)} finding(s) for "
            f"{self.params_file}",
        ]
        for sev in ("CRITICAL", "MAJOR", "MINOR"):
            group = self.by_severity(sev)
            if group:
                lines.append(f"  -- {sev} ({len(group)}) --")
                lines.extend(f"    {f}" for f in group)
        return "\n".join(lines)


# A C++ declaration of `name`: capture whatever type token precedes it.  One
# level of template nesting is enough for `std::vector<std::array<double, 3>>`
# and everything else these codes actually declare.
_TYPE_FRAGMENT = (
    r"(?:const\s+|static\s+|mutable\s+)*"
    r"[A-Za-z_][\w:]*"
    r"(?:\s*<[^<>;]*(?:<[^<>;]*>)?[^<>;]*>)?"
)


# `(` is in the lookahead for constructor-style declarations --
# `std::vector<int> extraction_modes_vect(2 * n);`.  It also matches ordinary
# calls, but those never classify as a container or scalar type, so they fall
# out harmlessly.
_DECL_RE = re.compile(
    rf"(?P<type>{_TYPE_FRAGMENT})\s*[&*]?\s+(?P<name>[A-Za-z_]\w*)\s*(?=[;=,()\[{{])"
)


def _declarations(text: str) -> dict[str, str]:
    """Map every declared name in `text` to ARITY_ARRAY or ARITY_SCALAR.

    One pass per file, not one pass per (name, file) pair -- the corpus is
    thousands of files and the naive form takes minutes.  A name declared as a
    container anywhere wins over a scalar declaration elsewhere: the question
    is whether *this* target can hold more than one token, and a false "array"
    only costs us a missed warning, while a false "scalar" costs a false
    accusation.
    """
    out: dict[str, str] = {}
    for m in _DECL_RE.finditer(text):
        type_str = m.group("type")
        if _ARRAY_TYPE_RE.search(type_str):
            out[m.group("name")] = ARITY_ARRAY
        elif _SCALAR_TYPE_RE.search(type_str):
            out.setdefault(m.group("name"), ARITY_SCALAR)
    return out


def _resolve_arity(
    var: str | None, own_decls: dict[str, str], pool: dict[str, str]
) -> str:
    """Decide how many tokens a read site can absorb, from its target's type.

    ``GRParmParse::load`` is overloaded on the destination type
    (Source/IO/GRParmParse.hpp:27-72), so the function name says nothing.  Look
    at the variable: same file first, then the whole scanned corpus, because
    members are declared in a header and read in a .cpp.
    """
    if var is None:
        return ARITY_UNKNOWN
    if var.endswith("]"):
        return ARITY_SCALAR  # `pp.load("k", arr[0])` absorbs exactly one token
    name = var.split(".")[-1].split("::")[-1]
    if not name:
        return ARITY_UNKNOWN
    return own_decls.get(name) or pool.get(name, ARITY_UNKNOWN)


@dataclass
class CxxIndex:
    """Everything the static scan learned about how C++ reads params.txt."""

    reads: dict[str, list[CxxRead]] = field(default_factory=dict)
    #: literal -> site, where the helper may suffix the literal to build keys
    prefixes: dict[str, CxxRead] = field(default_factory=dict)

    def sites_for(self, key: str) -> list[CxxRead]:
        return self.reads.get(key, [])

    def has_prefix_reader(self, key: str) -> bool:
        return any(key.startswith(p + "_") for p in self.prefixes)


def scan_cxx(roots: list[Path]) -> CxxIndex:
    """Scan C++ trees for ParmParse read sites, keyed by param name.

    Two passes: collect the sources, then resolve each read site's destination
    type.  The second pass needs the first to have finished, because a read in
    ``FooLevel.cpp`` usually targets a member declared in ``FooParams.hpp``.

    Non-literal keys (built from a ``std::ostringstream``, as the per-lump
    ``lump<k>_*`` readers are) cannot be resolved statically.  They are skipped
    here and handled by ``_DYNAMIC_KEY_PREFIXES`` in the audit, so that a
    dynamically-constructed reader does not masquerade as a dead key.
    """
    texts: dict[Path, str] = {}
    for root in roots:
        if not root.exists():
            continue
        paths = (
            [root]
            if root.is_file()
            else [
                p
                for ext in ("*.hpp", "*.cpp", "*.H", "*.cpp.in", "*.impl.hpp")
                for p in root.rglob(ext)
            ]
        )
        for path in paths:
            try:
                texts[path] = path.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue

    per_file: dict[Path, dict[str, str]] = {
        path: _declarations(text) for path, text in texts.items()
    }
    pool: dict[str, str] = {}
    for decls in per_file.values():
        for name, arity in decls.items():
            if arity == ARITY_ARRAY or name not in pool:
                pool[name] = arity

    index = CxxIndex()
    for path, text in texts.items():
        for m in _READ_RE.finditer(text):
            literal = m.group("lit")
            if literal is None:  # key built at runtime; not resolvable here
                continue
            var = m.group("var")
            index.reads.setdefault(literal, []).append(
                CxxRead(
                    key=literal,
                    fn=m.group("fn"),
                    file=str(path),
                    line=text.count("\n", 0, m.start()) + 1,
                    var=var,
                    arity=_resolve_arity(var, per_file[path], pool),
                )
            )

        for m in _HELPER_READ_RE.finditer(text):
            if m.group("fn") in _PARMPARSE_FNS:
                continue  # already handled above, and pp.get("k", v) is not this
            if "pp" not in m.group("pp").lower():
                continue  # first argument is not a ParmParse; unrelated call
            literal = m.group("lit")
            var = m.group("var")
            site = CxxRead(
                key=literal,
                fn=m.group("fn"),
                file=str(path),
                line=text.count("\n", 0, m.start()) + 1,
                var=var,
                arity=_resolve_arity(var, per_file[path], pool),
            )
            index.reads.setdefault(literal, []).append(site)
            # The helper may append "_0", "_1", ... to build the real keys.
            # Recorded with arity unknown: the read that matters happens inside
            # the helper, on a variable we are not looking at.
            index.prefixes.setdefault(
                literal, CxxRead(literal, site.fn, site.file, site.line, var)
            )
    return index


def parse_cxx_reads(roots: list[Path]) -> dict[str, list[CxxRead]]:
    """Exact-key view of :func:`scan_cxx`, kept for direct/legacy callers."""
    return scan_cxx(roots).reads


def parse_params_file(path: Path) -> list[EmittedKey]:
    """Parse an emitted params.txt into (key, value, line) records."""
    out: list[EmittedKey] = []
    for i, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.split("#", 1)[0]
        if not line.strip():
            continue
        m = _ASSIGN_RE.match(line)
        if m:
            out.append(EmittedKey(m.group(1), m.group(2), i))
    return out


# Keys assembled at runtime in C++ (ostringstream) that the static scan cannot
# see.  A key matching one of these prefixes is assumed to have a reader.
_DYNAMIC_KEY_PREFIXES: tuple[str, ...] = (
    "lump",
    "trajectory_lump",
    "rl_lump",
    "bh1_",
    "bh2_",
)

# ...but "assumed to have a reader" is only true up to the loop bound.  GRTeclyn
# reads `trajectory_lump<k>_*` for k < trajectory_params.num_lumps, and that
# count is silently clamped to GRTRESNA_MAX_INDEPENDENT_SCALARS on load
# (Examples/RadialRecipe/SimulationParameters.hpp:126).  Ask for eight lumps and
# you get five, with no warning: the keys for the rest are emitted, are wired to
# live search dimensions, and are read by nothing.
#
# This exemption is why the check missed it the first time.  A blanket "assume a
# reader exists" is the same mistake as a silent default -- it turns an unknown
# into a reassuring answer.
_INDEXED_KEY_RE = re.compile(r"^(?P<family>[a-z_]+?)(?P<idx>\d+)_")
_CAP_CONSTANT_RE = re.compile(
    r"\b(?:static\s+)?(?:constexpr|const)\s+int\s+"
    r"GRTRESNA_MAX_INDEPENDENT_SCALARS\s*=\s*(\d+)"
)
#: families whose index is capped by GRTRESNA_MAX_INDEPENDENT_SCALARS.
#: GRTresna's own bare ``lump<k>_*`` is NOT capped -- it honours num_lumps.
_CAPPED_KEY_FAMILIES: tuple[str, ...] = ("trajectory_lump", "rl_lump")


def scalar_slot_cap(roots: Sequence[Path] | None = None) -> int | None:
    """The number of independent scalar lumps GRTeclyn can actually drive.

    Read out of the C++ header rather than duplicated here, so raising the cap
    is a one-line change that the wrapper picks up automatically.  Returns
    ``None`` when the header is not reachable (installed wrapper, no source
    tree), in which case callers must not invent a value.
    """
    return _scalar_slot_cap(list(roots) if roots is not None else _default_cxx_roots())


def _default_cxx_roots() -> list[Path]:
    """Best-effort location of the GRTeclyn source tree from an installed wrapper."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "Source" / "Matter" / "GRTresnaScalarLayout.hpp"
        if candidate.exists():
            return [candidate]
    return []


def _scalar_slot_cap(roots: list[Path]) -> int | None:
    """Read GRTRESNA_MAX_INDEPENDENT_SCALARS from the C++ headers."""
    for root in roots:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else root.rglob("*.hpp")
        for path in paths:
            try:
                m = _CAP_CONSTANT_RE.search(
                    path.read_text(encoding="utf-8", errors="replace")
                )
            except OSError:
                continue
            if m:
                return int(m.group(1))
    return None

# Namespaces owned by AMReX/Chombo infrastructure rather than our physics.
_INFRA_PREFIXES: tuple[str, ...] = (
    "amr.",
    "amrex.",
    "geometry.",
    "particles.",
    "fabarray.",
    "eb2.",
    "blueprint.",
    "vismf.",
    "DistributionMapping.",
)

# Physics constants whose silent default is dangerous: if the wrapper ever stops
# emitting one, the C++ default takes over and the run is quietly wrong.  This
# is the bug-#14 (mu = 0) shape.  Emitting these is mandatory.
REQUIRED_PHYSICS_KEYS: tuple[str, ...] = (
    "recipe_scalar_mass",
    "recipe_scalar_lambda",
    "recipe_scalar_mu",
    "recipe_num_scalar_fields",
)


# Some emitted keys are consumed by the *wrapper*, not by C++:
# `trajectory_r_min` and `trajectory_retrograde_only` configure the Python
# clamps in ``search/optimize/candidates.py`` and are written to params.txt as
# provenance.  They are live knobs; only a scan that looks at both sides can
# tell them apart from a genuinely dead one like ``hdf5_subpath``.
_PY_KEY_RE = re.compile(r"""['"]([A-Za-z_][\w.]{3,})['"]""")


def _python_referenced_keys(roots: list[Path]) -> set[str]:
    """Every string literal in the Python sources that could name a param."""
    found: set[str] = set()
    for root in roots:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else root.rglob("*.py")
        for path in paths:
            try:
                text = path.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            found.update(m.group(1) for m in _PY_KEY_RE.finditer(text))
    return found


def _is_infra(key: str) -> bool:
    return key.startswith(_INFRA_PREFIXES)


def _is_dynamic(key: str) -> bool:
    return key.startswith(_DYNAMIC_KEY_PREFIXES)


def _index_beyond_cap(key: str, cap: int | None) -> tuple[str, int] | None:
    """Return (family, index) if `key` names a lump slot GRTeclyn cannot reach."""
    if cap is None:
        return None
    m = _INDEXED_KEY_RE.match(key)
    if not m:
        return None
    family = m.group("family")
    if family not in _CAPPED_KEY_FAMILIES:
        return None
    idx = int(m.group("idx"))
    return (family, idx) if idx >= cap else None


def audit_params_file(
    params_path: Path,
    cxx_roots: list[Path],
    *,
    python_roots: list[Path] | None = None,
    require_physics_keys: bool = True,
) -> ContractReport:
    """Cross-check one emitted params.txt against the C++ readers.

    Checks, in descending severity:

    CRITICAL  vector written, scalar read -- the bug-#19 shape.  Silently drops
              every token after the first.
    CRITICAL  a required physics constant is not emitted, so the C++ silent
              default (typically 0.0) governs -- the bug-#14 shape.
    MAJOR     key emitted that no C++ reader consumes -- a dead knob.  If it is
              wired to a search dimension, that dimension does nothing.
    MINOR     key read only via a silent default and never emitted -- benign
              today, a trap the moment the default stops being right.
    """
    emitted = parse_params_file(params_path)
    index = scan_cxx(cxx_roots)
    reads = index.reads
    py_keys = _python_referenced_keys(python_roots or [])
    slot_cap = _scalar_slot_cap(cxx_roots)
    report = ContractReport(
        params_file=str(params_path),
        n_emitted=len(emitted),
        n_cxx_keys=len(reads),
    )

    for item in emitted:
        if _is_infra(item.key):
            continue

        sites = index.sites_for(item.key)

        beyond = _index_beyond_cap(item.key, slot_cap)
        if beyond is not None:
            family, idx = beyond
            report.findings.append(
                Finding(
                    "CRITICAL",
                    "indexed-key-beyond-cap",
                    item.key,
                    f"line {item.line}: {family}{idx}_* is emitted, but GRTeclyn "
                    f"reads only indices 0..{slot_cap - 1} "
                    f"(GRTRESNA_MAX_INDEPENDENT_SCALARS = {slot_cap}, silently "
                    "clamped at SimulationParameters.hpp). The config asks for "
                    f"{idx + 1}+ lumps and gets {slot_cap}, with no diagnostic. "
                    "Any search dimension mapped here is inert.",
                )
            )
            continue

        if not sites:
            if (
                not _is_dynamic(item.key)
                and not index.has_prefix_reader(item.key)
                and item.key not in py_keys
            ):
                report.findings.append(
                    Finding(
                        "MAJOR",
                        "dead-key",
                        item.key,
                        f"line {item.line}: emitted but read by nothing -- no C++ "
                        "reader and no Python consumer. If a search dimension "
                        "maps here, it is inert.",
                    )
                )
            continue

        if item.is_vector and not any(s.is_array_read for s in sites):
            where = ", ".join(
                f"{Path(s.file).name}:{s.line} (pp.{s.fn} -> "
                f"{s.var or '?'}: {s.arity})"
                for s in sites
            )
            if all(s.is_proven_scalar for s in sites):
                report.findings.append(
                    Finding(
                        "CRITICAL",
                        "vector-read-as-scalar",
                        item.key,
                        f"line {item.line}: {item.n_tokens} tokens emitted "
                        f"({item.raw_value!r}) but every reader targets a scalar "
                        f"variable -- {where}. Only token 0 will be parsed. This "
                        "is the recipe_scalar_field_signs bug shape; declare the "
                        "target as a container or use countval+getarr.",
                    )
                )
            else:
                # We could not read the destination's type.  Say so; do not
                # accuse.  A checker that cries wolf gets muted, and then it
                # catches nothing at all.
                report.findings.append(
                    Finding(
                        "MINOR",
                        "unverified-arity",
                        item.key,
                        f"line {item.line}: {item.n_tokens} tokens emitted but "
                        f"the reader's destination type could not be resolved -- "
                        f"{where}. Confirm by hand that all tokens are consumed.",
                    )
                )

    if require_physics_keys:
        emitted_keys = {e.key for e in emitted}
        for key in REQUIRED_PHYSICS_KEYS:
            if key in emitted_keys or key not in reads:
                continue
            defaults = [s for s in reads[key] if s.has_silent_default]
            if defaults:
                s = defaults[0]
                report.findings.append(
                    Finding(
                        "CRITICAL",
                        "missing-physics-constant",
                        key,
                        f"not emitted, but read at {Path(s.file).name}:{s.line} "
                        f"with a silent default (pp.{s.fn}). The C++ default "
                        "governs and the run is quietly wrong -- this is the "
                        "recipe_scalar_mu = 0 bug shape.",
                    )
                )

    return report
