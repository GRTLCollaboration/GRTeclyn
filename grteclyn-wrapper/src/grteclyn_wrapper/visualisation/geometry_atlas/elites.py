"""Load top geometry-atlas elites from a campaign directory."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ...search.geometry_atlas.ansatz import ANSATZE, ansatz_offset
from ...search.geometry_atlas.genome import PARAMS_PER_CENTER, GeometryGenome


@dataclass(frozen=True)
class EliteRecord:
    """One archive elite with genome + scalar diagnostics."""

    rank: int
    cell: tuple[int, int]
    eval_id: int | None
    f_geo: float
    score: float
    shift_fraction: float
    log_exotic_energy: float
    integral_negative_rho: float
    min_rho: float
    h_quality_ok: bool | None
    genome: GeometryGenome
    diagnostics: dict[str, Any]
    source_path: Path

    @property
    def family_label(self) -> str:
        sf = self.shift_fraction
        if sf >= 0.9:
            return "shift-tube"
        if sf >= 0.5:
            return "hybrid"
        if sf >= 0.2:
            return "curvature-leaning"
        return "curvature-dominated"

    def active_modes(self, *, amp_tol: float = 0.05) -> list[str]:
        """Human-readable list of active analytic topologies."""
        coeffs = self.genome.coeffs
        n_c = self.genome.config.n_centers
        base = n_c * PARAMS_PER_CENTER
        labels: list[str] = []
        for ansatz in ANSATZE:
            off = base + ansatz_offset(ansatz.name)
            block = coeffs[off : off + ansatz.n_params]
            if ansatz.name == "alcubierre":
                if abs(float(block[0])) > amp_tol:
                    labels.append(
                        f"alc v={block[0]:.2f} R={block[1]:.1f} σ={block[2]:.1f}"
                    )
            elif abs(float(block[0])) > amp_tol:
                labels.append(f"{ansatz.name} s={block[0]:.2f}")
        return labels or ["RBF-only"]


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _genome_from_payload(payload: dict) -> GeometryGenome | None:
    gdata = payload.get("genome")
    if not isinstance(gdata, dict):
        return None
    try:
        return GeometryGenome.from_dict(gdata)
    except Exception:  # noqa: BLE001 — skip unreadable elites
        return None


def _elite_from_eval_payload(
    payload: dict,
    *,
    rank: int,
    source: Path,
    cell: tuple[int, int] | None = None,
) -> EliteRecord | None:
    genome = _genome_from_payload(payload)
    if genome is None:
        return None
    cell_raw = cell or tuple(payload.get("cell", (0, 0)))
    desc = payload.get("descriptors") or [0.0, 0.0]
    # Prefer descriptor_details when present (elite JSON).
    details = payload.get("descriptor_details") or {}
    shift = float(details.get("shift_fraction", desc[0] if desc else 0.0))
    log_e = float(details.get("log_exotic_energy", desc[1] if len(desc) > 1 else 0.0))
    diag = dict(payload.get("diagnostics") or {})
    return EliteRecord(
        rank=rank,
        cell=(int(cell_raw[0]), int(cell_raw[1])),
        eval_id=payload.get("eval_id"),
        f_geo=float(payload.get("f_geo", 0.0)),
        score=float(payload.get("score", 0.0)),
        shift_fraction=shift,
        log_exotic_energy=log_e,
        integral_negative_rho=float(payload.get("integral_negative_rho", 0.0)),
        min_rho=float(payload.get("min_rho", diag.get("min_rho", 0.0))),
        h_quality_ok=payload.get("h_quality_ok"),
        genome=genome,
        diagnostics=diag,
        source_path=source,
    )


def load_top_elites(run_dir: Path | str, *, top_n: int = 5) -> list[EliteRecord]:
    """Return the top-``n`` archive elites ranked by frozen ``f_geo``.

    Prefers ``elites/cell_i_j.json`` (has genome); falls back to the matching
    ``evals/eval_XXXXXX.json`` referenced by the archive.
    """
    run_dir = Path(run_dir)
    archive = _load_json(run_dir / "archive.json")
    cells = archive.get("cells") or {}
    items = list(cells.values()) if isinstance(cells, dict) else list(cells)

    ranked = sorted(
        items,
        key=lambda e: float((e.get("objectives") or e).get("f_geo", 0.0)),
        reverse=True,
    )

    out: list[EliteRecord] = []
    for rank, elite in enumerate(ranked[:top_n], start=1):
        cell = tuple(elite.get("cell", (0, 0)))
        cell_t = (int(cell[0]), int(cell[1]))
        obj = elite.get("objectives") or {}
        details = elite.get("descriptor_details") or {}

        # 1) Prefer persisted elite JSON.
        elite_json = run_dir / "elites" / f"cell_{cell_t[0]}_{cell_t[1]}.json"
        payload: dict | None = None
        source = elite_json
        if elite_json.exists():
            payload = _load_json(elite_json)
            # Fill missing f_geo from archive objectives.
            payload.setdefault("f_geo", obj.get("f_geo", 0.0))
            payload.setdefault("score", elite.get("score", 0.0))
            payload.setdefault("cell", list(cell_t))
            payload.setdefault("descriptor_details", details)
            payload.setdefault(
                "integral_negative_rho", obj.get("integral_negative_rho", 0.0)
            )
            payload.setdefault("min_rho", obj.get("min_rho", 0.0))
            if "eval_id" not in payload and "params" in elite:
                payload["eval_id"] = int(elite["params"].get("eval_id", -1))

        # 2) Fallback: eval JSON.
        if payload is None or payload.get("genome") is None:
            eval_id = int((elite.get("params") or {}).get("eval_id", -1))
            eval_json = run_dir / "evals" / f"eval_{eval_id:06d}.json"
            if eval_json.exists():
                payload = _load_json(eval_json)
                source = eval_json
                payload.setdefault("descriptor_details", details)

        if payload is None:
            continue
        rec = _elite_from_eval_payload(payload, rank=rank, source=source, cell=cell_t)
        if rec is not None:
            out.append(rec)
    return out


__all__ = ["EliteRecord", "load_top_elites"]
