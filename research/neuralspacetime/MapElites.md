## MAP-Elites FTL Discovery Status - 2026-06-09

### What changed

The previous QD/HQ runs did not really score evolved FTL because plotfiles were deleted before evolved/geodesic FTL and effective energy conditions were computed. The new pipeline keeps the last 3 plotfiles, adds an evolved-FTL persistence window, adds a `speed_horizon` descriptor (`max_local_speed` vs `min_theta_plus`), and fixes scalar frame extraction to use evolved `phi_lump*` / `Pi_lump*` channels instead of frozen aggregate `phi` / `Pi`.

### Confirmed HQ survivor: `000114`

Source: `runs/grtresna_promote/revalidate_keep3_l128n256t30_qd_eval000114`

Settings: `L=128`, `N=256`, `dx=0.5`, `t=30`, `consumer_keep_last=3`.

Result:

- Evolved operational FTL: `f_op = 0.07847`
- Geodesic FTL: `f_geo = 0.06936`
- Sustained window: `f_op_min = 0.07857`, `f_op_median = 0.07879`
- Evolved max local speed: `1.11575`
- Evolved max shift: `0.03299`
- Horizon proxy: `min_theta_plus = 0.00529`, `max_horizon_radius = 0.0`
- Final constraints: Ham `3.06e-4`, Mom `2.97e-5`

Interpretation: this is the first strong HQ candidate from the current campaign that passes the working FTL discovery gate: evolved `f_op > 0.03`, geodesic `f_geo > 0.001`, no trapped-surface proxy, and persistence across the retained plotfile window.

### Matter content of `000114`

Matter model: `grtresna_independent_scalars`

Scalar fields: 5 independent lumps, scalar mass `0.1`.

Exotic matter: yes. Signs are `[-1, -1, -1, -1, +1]`, so 4 of 5 scalar fields are phantom/exotic and 1 is canonical. Matter/energy-condition diagnostics still show exotic support:

- `integral_negative_rho = 0.7093`
- `min_rho_required = -1.58e-3`
- matter NEC min `-1.79e-3`
- effective evolved NEC min `-2.66e-4`

Shell ansatz parameters:

- shell amplitude `0.1977`
- radius `3.566`
- width `3.0`
- thickness `0.306`
- toroidal velocity `0.326`
- poloidal velocity `0.607`
- radial velocity `-0.123`
- omega `0.00479`
- mode `2`
- exotic fraction `0.769`

GRTresna convergence: iteration `6`, Ham `0.996%`, Mom `0.069%`; postload gate passed.

### Is it a new FTL geometry?

It is a new candidate geometry from the MAP-Elites/GRTresna search path, and it survives the current HQ evolved/geodesic FTL validation at `L=128,N=256,t=30`.

It is not yet a final physical discovery. Remaining caveats:

- It uses exotic/phantom scalar matter.
- Geodesic diagnostic reports null-constraint drift above the strict quality flag threshold (`max H = 4.83e-4`), though all 5 rays reached and `f_geo` is strongly positive.
- Needs repeat validation at higher resolution / longer time and ideally reduced exotic matter.

### Current new QD campaign

Path: `runs/grtresna_qd/ftl_discovery/qd_ftl_discovery_20260609T162553Z`

Purpose: retargeted FTL-first MAP-Elites search with retention, persistence scoring, `speed_horizon` descriptors, and `grtresna_shift_seed`.

Current notable candidates:

- `eval_000025`: highest score so far (`1797`), strong FTL metrics but trapped-surface penalty active, so not the safe target.
- `eval_000024`: horizon-safe channel candidate at low resolution (`f_op = 0.0438`, `f_geo = 0.0230`, sustained `f_op_min = 0.0394`, `min_theta_plus = 0.0126`). It is currently being promoted to HQ as `runs/grtresna_promote/l128n256t30_ftl_discovery_qd_eval000024`.
- `eval_000001`: seeded survivor with high score (`1784`) and no horizon penalty in the low-res QD run.

### Next restart point

Use `000114` as the confirmed HQ anchor. Continue watching/promoting `eval_000024` and future horizon-safe channel elites. The next science target is to reduce exotic matter while retaining evolved/geodesic FTL and to repeat the HQ validation at higher resolution/longer time.
