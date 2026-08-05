"""Single source of truth for ``--objective-mode`` whitelists.

Every parser that validates an objective mode MUST import these tuples instead
of holding its own literal copy.  The mode lists were previously copy-pasted in
four places (reproduce/optimize/qd in ``cli/parser.py`` plus the plotfile
consumer's own CLI); adding ``f_geo_max`` to only three of them made every
consumer launch die at argparse time -- and the pipeline swallowed that crash,
silently degrading the evolving-geodesic trace to a too-short plotfile window
where every f_geo reads 0.  A whitelist that can drift is a physics bug, not a
style problem.

Deliberately dependency-free so the CLI and the consumer subprocess can both
import it without dragging in numpy/scoring at argparse time.
"""

from __future__ import annotations

# Modes accepted by every entry point (reproduce, optimize/CMA-ES, qd, and the
# plotfile consumer's incremental scorer).  Must stay in sync with the dispatch
# in ``metrics/score/objectives.py::compute_total``.
OBJECTIVE_MODES: tuple[str, ...] = (
    "weighted",
    "ftl_first",
    "robust_ftl",
    "general_ftl",
    "f_geo_max",
    "f_geo_depth",
    "critical_collapse",
    "gw_beam",
    "spacetime_shear",
)

# Modes only meaningful for the MAP-Elites qd search (need qd-side state such
# as a target motif).  The consumer still accepts them: it merely tags
# incremental scores and must never crash on a mode the campaign runs.
QD_ONLY_MODES: tuple[str, ...] = ("geometry_first",)

QD_OBJECTIVE_MODES: tuple[str, ...] = OBJECTIVE_MODES + QD_ONLY_MODES
