# Lab Journal

## 2026-07-03 — `gw_beam_qd100_v4` complete: eval 61 / 88 analysis

Campaign finished **100/100**. Hard gates held: **77/100** collapse modes @ ~**−116** (`mult=0`); **22** healthy survivors; **5** archive elites. Best score **eval 88 @ 3.09**; best beaming **eval 61 @ 2.82** (beam_ratio **~30%**).

Run dir: `runs/grtresna_qd/gw_beam_qd100_v4/`

### Campaign headline

| | eval 88 (best score) | eval 61 (best beam) |
|--|----------------------|---------------------|
| Score | **3.09** | 2.82 |
| Tier | constructed | constructed |
| `gw_health_multiplier` | 1.0 | 1.0 |
| max ‖Ham‖₂ | 0.08 | 0.14 |
| mean Ψ₄ power | 4.5×10⁻⁴ | **6.4×10⁻⁴** |
| beam_ratio | 14% | **~30%** (to 40% late) |
| Scorer bias | high survival (1.0) | better GW, slightly noisier grid |

**Critical verdict:** Neither run is a strong GW emitter. Both produce a **weak steady hum** (P ~ 10⁻⁵–10⁻⁴ after t≈5), not a merger chirp or beamed burst. The t=0 Ψ₄ spike is **initial/near-zone transient**, not radiative physics. Eval **61** is the better physics candidate for directional GW; eval **88** wins total score because survival/stability bonuses dominate over `gw_beam_quality`.

### Eval 61 — lump dynamics (best beaming)

**Geometry:** five Q-balls in a **compact ring** R ≈ 2.4–4.1; all tangential speeds **~0.27–0.30** (near trajectory speed cap).

| Lump | R₀ | v_t | v_rad |
|------|-----|-----|-------|
| 0 | 4.06 | 0.30 | −0.01 |
| 1 | 2.48 | 0.30 | +0.03 |
| 2 | 2.73 | 0.30 | −0.01 |
| 3 | 2.38 | 0.27 | −0.03 |
| 4 | 3.69 | 0.28 | **+0.10** (outward) |

**Global:** breathing A ≈ 1.55 (T ≈ 4.7); Z bob amp ≈ 0.93 (T ≈ 4.5). Slow retrograde orbits (T_orb 52–85) — azimuthal drift is small over t = 24.

**Ψ₄ at R = 20:** t = 0 spike P ≈ 3.6×10⁻³; steady P ≈ 3–4×10⁻⁵ with **beam_ratio ≈ 30%** flat from t ≈ 5–24. **Constraints:** Ham rises slowly 0.06 → 0.14; Mom ~10⁻³–10⁻²; **no spike**.

**Mechanism (critical):** ~30% Z-beaming comes from a **coherent, fast, compact multi-lump cluster** + breathing quadrupole — not a clean radiative-zone binary. Extraction at R = 20 still sees **near-zone** motion from R ~ 2–4 lumps.

**Plots:** `eval_000061/plots/constraints_plot.png`, `eval_000061/plots/psi4_analysis_r20.png`

### Eval 88 — lump dynamics (best score)

**Geometry:** four **inner** lumps R ≈ 1.7–1.9 + one **outer** lump R = 6.5 (v_rad = +0.065 outward). Tangential speeds **slower** (0.08–0.30) than eval 61.

**Global:** stronger breathing A ≈ 1.77 (T ≈ 4.9); larger Z bob (amp 1.23). Orbits very slow (T_orb 67–144).

**Ψ₄ at R = 20:** t = 0 spike; steady P ≈ 2–3×10⁻⁵; beam_ratio **~14%** only. max Ham ≈ 0.08 — very stable, weaker directional signal.

**Why higher score:** `survival = 1.0`, lower instability penalty — scorer rewards numerical cleanliness over beaming.

### Implications

1. v4 gates work; no repeat of v3 **336** exploit.
2. Search has **not** found a strong GW laser — only stable breathing multi-lump configs with 10⁻⁵–10⁻⁴ power.
3. Total score **underweights beam_ratio** relative to survival — consider rebalancing or CMA-ES warm-start from **eval 61**.
4. Next: longer `stop_time`, R = 24 extraction, tighter `v_rad` / breathing bounds, inspect eval 61 lump trajectories visually.

---

## 2026-07-03 — `gw_beam_qd100_v3` → v4: collapse-mode reward hacking

### The exploit (classic ML reward hacking in NR)

Asked MAP-Elites to maximize GW power. The optimizer did **not** build a gravitational-wave laser — it built a **numerical bomb**.

> Ψ₄ uses second derivatives of the metric. Crash the Hamiltonian constraint → grid fills with high-frequency numerical noise → second-derivative operators report “infinite wave power.”

**Mechanism:** 5 pumped Q-balls, strong breathing + inward `v_rad` → lapse pinches → **Ham/Mom cliff at t≈4.5** → sim keeps running (`exit_code=0`) → Ψ₄ burst scored as GW. Eval 51 hit **336** vs trustworthy eval 7 @ **3.4**.

**Why additive penalties failed:** `−113` spike penalty vs `+197` from uncapped **mean** Ψ₄. You cannot fight unbounded exploits with additive terms.

### Partial fixes (v3, insufficient)

| Fix | Effect |
|-----|--------|
| Pre-spike **peak** cap only | Peak → 0; **mean** still post-collapse → eval 51 @ **81** |
| `constraint_spike_penalty ×150` | Additive; still dominated by quality term |
| Tier = rejected | Ignored by archive — rejected runs still filled cells |

### Permanent closure (v4 hard gates)

**Fix A — Ψ₄ time-series truncation** (`read_psi4_metrics`)
- At `t_spike`, truncate **peak, mean, and final** to `[0, t_spike]` only
- If `t_spike < 10` (wave hasn't reached R=12 cleanly), GW metrics invalid

**Fix B — Archive admission** (`qd_search/driver.py`)
- `tier >= CONSTRUCTED` required for `archive.insert()`
- Rejected collapse runs stay in trajectory / near-miss pool, never occupy descriptor cells

**Fix C — Multiplicative health gate** (`gw_beam_total`)
- `gw_health_multiplier = 0` if `max(‖Ham‖₂) > 100` or early spike (`t_spike < 10`)
- `score = (1000×quality + 100×peak + health bonuses) × multiplier + penalties`
- Collapse eval 51: **336 → −116** (GW terms zeroed via `mult=0`; remainder is additive penalties)

### v3 rescore (hard gates applied 2026-07-03)

- Stopped v3 at 100/100 evals; rescored 19 rows with surviving eval dirs
- Collapse evals 42/51/70: **−116**, `gw_health_multiplier=0`
- Archive rebuilt: **8 elites**, best **3.37** (eval 7 trajectory row; eval 91 @ 2.70 best on disk)

### Physics notes for v4 campaign

1. **R = 20** extraction (was R = 12 in v3) — farther wave zone, early-spike gate at t ≈ 18.
2. **Inward `v_rad`** with canonical matter → central collapse. Clamped `v_rad > −0.1`.

### Next run

v4 **complete** (see post-run analysis above). v2/v3 run dirs deleted.

Plots: `eval_000051/plots/`, `eval_000042/plots/` (`constraints_plot.png`, `psi4_analysis_r12.png`).
