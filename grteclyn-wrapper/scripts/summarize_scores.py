import json
import sys
from pathlib import Path

from grteclyn_wrapper.metrics import read_episode_metrics
from grteclyn_wrapper.metrics import score_episode

runs = Path(sys.argv[1])
target = float(sys.argv[2]) if len(sys.argv) > 2 else 16.0

eps = sorted(p.parent for p in runs.glob("*/metadata.json") if "_dry_" not in p.parent.name)
rows = []
for ep in eps:
    meta = json.loads((ep / "metadata.json").read_text())
    ov = meta.get("overrides", {}) or meta.get("params_overrides", {}) or {}
    exotic = ov.get("recipe_exotic_matter", meta.get("recipe_exotic_matter", "?"))
    try:
        m = read_episode_metrics(ep, ftl_L=8.0)
        s = score_episode(m, target_stop_time=target)
        c = s.components
        ec = m.energy_conditions
        cv = m.curvature
        g0 = m.general_ftl
        ge = m.general_ftl_evolved
        rows.append({
            "name": ep.name.split("_gpu_")[0],
            "exotic": exotic,
            "score": round(s.total, 2),
            "survival": round(c.get("survival", 0.0), 3),
            "stability": round(c.get("stability", 0.0), 3),
            "s_comov": round(c.get("comoving_stability", 0.0), 3),
            "F_log": round(c.get("ftl_shortcut", 0.0), 3),
            "op_ftl": round(c.get("operational_ftl", 0.0), 3),
            "curv_act": round(c.get("curvature_activity", 0.0), 3),
            "min_nec": None if not ec else (None if ec.min_nec is None else round(ec.min_nec, 4)),
            "intNECviol": None if not ec else (None if ec.max_integral_nec_violation is None else round(ec.max_integral_nec_violation, 3)),
            "maxRicci": None if not cv else (None if cv.max_abs_ricci_scalar is None else round(cv.max_abs_ricci_scalar, 3)),
            "L2Ricci": None if not cv else (None if cv.max_l2_ricci_scalar is None else round(cv.max_l2_ricci_scalar, 3)),
            "ftl0_fop": None if not g0 else round(g0.f_op, 4),
            "ftlE_fop": None if not ge else round(ge.f_op, 4),
            "gate": round(c.get("nontriviality_gate", 1.0), 3),
            "eff_exotic": round(c.get("effective_exoticity", 0.0), 4),
            "eff_nec": None if not m.effective_ec else (None if m.effective_ec.nec_min is None else round(m.effective_ec.nec_min, 6)),
        })
    except Exception as exc:
        rows.append({"name": ep.name, "error": repr(exc)})

print(json.dumps(rows, indent=2))
print()
hdr = ["name", "score", "gate", "survival", "stability", "F_log", "op_ftl", "curv_act", "eff_exotic", "eff_nec", "min_nec", "ftlE_fop"]
print(" ".join(f"{h:>11}" for h in hdr))
for r in rows:
    if "error" in r:
        print(f"{r['name']:>11}  ERROR {r['error']}")
        continue
    print(" ".join(f"{str(r.get(h)):>11}" for h in hdr))
