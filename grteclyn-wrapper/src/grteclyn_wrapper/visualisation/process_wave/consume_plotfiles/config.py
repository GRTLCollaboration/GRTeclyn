from __future__ import annotations

import os
from typing import Dict

from ....core.config import VISUALISATION_DIR, default_sim_data_dir

def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


# Fast frame defaults (~9x speedup on 256^3 AMR HQ runs; override to restore quality).
# env: GRTECLYN_FRAMES_DPI, GRTECLYN_FRAMES_BUFF_CAP, GRTECLYN_FRAMES_SAMPLES_PER_CELL
_FRAME_DPI = _env_int("GRTECLYN_FRAMES_DPI", 90)
_FRAME_SAMPLES_PER_CELL = _env_int("GRTECLYN_FRAMES_SAMPLES_PER_CELL", 1)
_FRAME_BUFF_CAP = _env_int("GRTECLYN_FRAMES_BUFF_CAP", 512)


def _default_data_dir() -> str:
    return str(default_sim_data_dir())


def _default_frames_out_dir() -> str:
    return str(VISUALISATION_DIR / "visualize")


def _frames_auto_zlim_enabled(explicit: bool | None = None) -> bool:
    """Per-frame percentile scaling (colorbar jumps). Default: fixed limits for movies."""
    if explicit is not None:
        return explicit
    return os.environ.get("GRTECLYN_FRAMES_AUTO_ZLIM", "").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


# Fixed color limits keep movie colorbars stable across time steps.
# Override any field via GRTECLYN_FRAMES_ZLIM_<FIELD>=lo,hi (e.g. shift1=-0.05,0.05).
_FIELD_FRAME_CONFIGS: Dict[str, dict] = {
    # chi develops a deep, localised conformal well (min_chi can reach ~0.4, and
    # deeper near strong fields) while the far field stays ~1.0.  The window is
    # widened to (0.1, 1.05) so the well is resolved without per-frame rescaling:
    # a fixed scale keeps the movie colorbar stable (no "bouncing").  Re-enable
    # per-frame scaling explicitly via GRTECLYN_FRAMES_AUTO_ZLIM if needed.
    # per_frame_zlim: local percentile contrast each step (global t=0 lock washes out
    # late-time wells on chi; fixed ±0.9 on chi-1 hides O(0.01) features).
    "chi": {
        "zlim": (0.1, 1.05),
        "cmap": "magma",
        "label": r"Conformal Factor $\chi$",
        "per_frame_zlim": True,
    },
    "chi_minus_1": {
        "zlim": (-0.9, 0.9),
        "cmap": "RdBu",
        "label": r"$\chi - 1$",
        "per_frame_zlim": True,
    },
    "phi": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\phi$"},
    "Pi": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\Pi$"},
    "phi_lump0": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\phi_2$"},
    "Pi_lump0": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\Pi_2$"},
    # Bicomplex phantom field (Phi-) channels: phi_lump1/Pi_lump1 = Re/Pi,
    # phi_lump2/Pi_lump2 = Im components (see StateVariables.hpp).
    "phi_lump1": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\phi^-_{\mathrm{Re}}$"},
    "Pi_lump1": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\Pi^-_{\mathrm{Re}}$"},
    "phi_lump2": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\phi^-_{\mathrm{Im}}$"},
    "Pi_lump2": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\Pi^-_{\mathrm{Im}}$"},
    "phi_lump_sum": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\sum_k\phi_k$"},
    "Pi_lump_sum": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\sum_k\Pi_k$"},
    "scalar_activity": {"zlim": (0.0, 0.20), "cmap": "viridis", "label": r"$\sum_k\sqrt{\phi_k^2+\Pi_k^2}$"},
    "lump_activity": {"zlim": (0.0, 0.20), "cmap": "viridis", "label": r"$\sum_k\sqrt{\phi_k^2+\Pi_k^2}$"},
    "local_speed": {"zlim": (0.90, 1.30), "cmap": "magma", "label": r"Local Coordinate Speed"},
    "shift1": {
        "zlim": (-0.05, 0.05),
        "cmap": "RdBu",
        "label": r"shift1",
        "per_frame_zlim": True,
    },
    "shift2": {
        "zlim": (-0.05, 0.05),
        "cmap": "RdBu",
        "label": r"shift2",
        "per_frame_zlim": True,
    },
    "shift3": {
        "zlim": (-0.05, 0.05),
        "cmap": "RdBu",
        "label": r"shift3",
        "per_frame_zlim": True,
    },
    "rho_req": {"zlim": (-3.0e-3, 3.0e-3), "cmap": "RdBu", "label": r"$\rho_{\mathrm{req}}$"},
    "K": {"zlim": (-5.0e-4, 5.0e-4), "cmap": "RdBu", "label": r"Trace of Extrinsic Curvature $K$"},
    "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": r"Z4 Constraint $\Theta$"},
    "lapse": {"zlim": (0.995, 1.005), "cmap": "viridis", "label": r"Lapse $\alpha$"},
    "Ham": {"zlim": (-0.1, 0.1), "cmap": "RdBu", "label": r"Hamiltonian Constraint $\mathcal{H}$"},
    "A11": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{11}$"},
    "A12": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{12}$"},
    "A22": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{22}$"},
    "A33": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{33}$"},
    "GW_Plus": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain $h_+$"},
    "GW_Cross": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain $h_\times$"},
    "Weyl4_Re": {"zlim": (-1.0e-4, 1.0e-4), "cmap": "PiYG", "label": r"$\mathrm{Re}(\Psi_4)$"},
    "Weyl4_Im": {"zlim": (-1.0e-4, 1.0e-4), "cmap": "PiYG", "label": r"$\mathrm{Im}(\Psi_4)$"},
    "Weyl4_Mag": {"zlim": (0.0, 1.0e-4), "cmap": "inferno", "label": r"$|\Psi_4|$"},
}


def _field_frame_config(field: str) -> dict:
    cfg = dict(_FIELD_FRAME_CONFIGS.get(field, {"zlim": (None, None), "cmap": "viridis", "label": field}))
    env_key = f"GRTECLYN_FRAMES_ZLIM_{field.upper()}"
    override = os.environ.get(env_key, "").strip()
    if override:
        parts = [p.strip() for p in override.split(",")]
        if len(parts) == 2:
            cfg["zlim"] = (float(parts[0]), float(parts[1]))
    return cfg
