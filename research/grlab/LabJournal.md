# Lab Journal

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
2. **Inward `v_rad`** with canonical (attractive) matter → instant central collapse. Consider clamping `v_rad > −0.1` in search space.

### Next run

Launch **`gw_beam_qd100_v4`** with hard gates + rescore/restart workflow. Stop v3; rescore with `scripts/search/rescore_gw_beam_campaign.py` before resume if scorer changes mid-campaign.

Plots: `eval_000051/plots/`, `eval_000042/plots/` (`constraints_plot.png`, `psi4_analysis_r12.png`).
