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

### HQ promotion result: `ftl_discovery/eval_000024`

Source: `runs/grtresna_promote/l128n256t30_ftl_discovery_qd_eval000024`

Settings: `L=128`, `N=256`, `dx=0.5`, `t=30`, `consumer_keep_last=3`.

Result:

- Total score: `1283.36`
- Evolved operational FTL: `f_op = 0.04773`
- Sustained window: `f_op_min = 0.04774`, `f_op_median = 0.04801`, `n = 3`
- Geodesic FTL: `f_geo = 0.04230` (`4/5` rays reached; quality flag fails because `max H = 5.22e-4`)
- Evolved max local speed: `1.07070`
- Evolved max shift: `0.01497`
- Horizon/trapped proxy: `min_theta_plus = -0.02666`, `max_horizon_radius = 113.48`, `horizon_penalty = -1.0`
- Final constraints: Ham `1.57e-4`, Mom `1.34e-5`
- Exotic support: WEC min `-2.48e-2`, NEC min `-4.21e-3`, effective NEC min `-1.31e-4`

Interpretation: `eval_000024` is a real HQ evolved/geodesic FTL channel and the FTL signal is persistent over the retained final plotfiles. It is not the safe target from this campaign because the trapped-surface proxy is active. In short: good operational channel, bad horizon safety.

Movies were generated for all 15 frame sequences under `runs/grtresna_promote/l128n256t30_ftl_discovery_qd_eval000024/frames`. The most useful first-pass movies are `lump_activity_z/movie_lump_activity_z.mp4`, `phi_lump_sum_z/movie_phi_lump_sum_z.mp4`, `Pi_lump_sum_z/movie_Pi_lump_sum_z.mp4`, `local_speed_z/movie_local_speed_z.mp4`, `rho_req_z/movie_rho_req_z.mp4`, and the `lump_activity_proj_*` projection movies.

Comparison with the previous confirmed HQ survivor `000114`:

- `000114` is still the better safe anchor: it has stronger FTL (`f_op = 0.07847`, `f_geo = 0.06936`) and remains horizon-safe (`min_theta_plus = 0.00529`, no trapped proxy).
- `000024` is cleaner on final constraints than `000114` (Ham `1.57e-4` vs `3.06e-4`, Mom `1.34e-5` vs `2.97e-5`), but its negative `min_theta_plus` and `horizon_penalty = -1.0` make it physically less useful.
- Both require exotic matter. `000024` has clearly negative WEC/NEC diagnostics and should be treated as a trapped high-FTL mechanism to learn from, not as the current discovery target.

### HQ promotion results: `ftl_discovery/eval_000055` and `eval_000025`

Sources:

- `runs/grtresna_promote/l128n256t30_ftl_discovery_qd_eval000055`
- `runs/grtresna_promote/l128n256t30_ftl_discovery_qd_eval000025`

Settings for both: `L=128`, `N=256`, `dx=0.5`, `t=30`, `consumer_keep_last=3`.

`eval_000055` result:

- Total score: `1731.33`
- Evolved operational FTL: `f_op = 0.06554`
- Sustained window: `f_op_min = 0.06569`, `f_op_median = 0.06595`, `n = 3`
- Geodesic FTL: `f_geo = 0.05071` (`5/5` rays reached; quality flag fails because `max H = 3.48e-4`)
- Evolved max local speed: `1.12182`
- Evolved max shift: `0.00994`
- Horizon/trapped proxy: `min_theta_plus = -0.05234`, `max_horizon_radius = 118.70`, `horizon_penalty = -1.0`
- Final constraints: Ham `1.34e-4`, Mom `1.08e-5`
- Exotic support: WEC min `-1.03e-2`, NEC min `-1.75e-3`, effective NEC min `-1.46e-4`

Interpretation: `eval_000055` is the strongest new HQ promotion from the retargeted QD campaign by score and has a real geodesic/evolved FTL signal. It is not safe: the trapped proxy is stronger than in `000024`, and the stability component is very poor (`0.0083`), so this is a useful high-FTL trapped mechanism, not a discovery anchor.

`eval_000025` result:

- Total score: `141.25`
- Evolved operational FTL: `f_op = 0.05466`
- Sustained window: `f_op_min = 0.05461`, `f_op_median = 0.05510`, `n = 3`
- Geodesic FTL: `f_geo = 0.0` (`0/5` rays reached)
- Evolved max local speed: `1.17399`
- Evolved max shift: `0.01978`
- Horizon/trapped proxy: `min_theta_plus = -0.05051`, `max_horizon_radius = 116.24`, `horizon_penalty = -1.0`
- Final constraints: Ham `2.90e-4`, Mom `2.96e-5`
- Exotic support: WEC min `-1.88e-2`, NEC min `-3.21e-3`, effective NEC min `-2.69e-4`

Interpretation: `eval_000025` should be downgraded. It has a persistent coordinate superluminal channel, but the geodesic diagnostic finds no successful shortcut (`f_geo = 0`, `0/5` rays reached). Combined with the trapped proxy, large exotic penalty, and very poor stability (`0.0084`), this looks like a gauge/trapped artifact rather than a robust FTL geometry.

Movies were generated for all 15 frame sequences in both HQ folders. First-pass review should focus on `local_speed_z`, `shift1_z`, `chi_minus_1_z`, `rho_req_z`, `lump_activity_z`, and the `lump_activity_proj_*` movies.

HQ comparison table. Scores are the objective totals from `score.json`; they are useful for ranking search pressure, but the critical scientific ranking must also require geodesic FTL and horizon safety.

<table>
  <colgroup>
    <col style="width: 10%;">
    <col style="width: 7%;">
    <col style="width: 7%;">
    <col style="width: 7%;">
    <col style="width: 7%;">
    <col style="width: 8%;">
    <col style="width: 8%;">
    <col style="width: 12%;">
    <col style="width: 6%;">
    <col style="width: 40%; min-width: 420px;">
  </colgroup>
  <thead>
    <tr>
      <th>HQ run</th>
      <th align="right">Score</th>
      <th align="right">Evolved <code>f_op</code></th>
      <th align="right">Geodesic <code>f_geo</code></th>
      <th align="right">Rays reached</th>
      <th align="right"><code>min_theta_plus</code></th>
      <th align="right">Horizon penalty</th>
      <th>Final Ham/Mom L2</th>
      <th align="right">Stability</th>
      <th>Critical status</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>000114</code></td>
      <td align="right"><code>2170.53</code></td>
      <td align="right"><code>0.07847</code></td>
      <td align="right"><code>0.06936</code></td>
      <td align="right"><code>5/5</code></td>
      <td align="right"><code>0.00529</code></td>
      <td align="right"><code>0.0</code></td>
      <td><code>3.06e-4</code> / <code>2.97e-5</code></td>
      <td align="right"><code>0.788</code></td>
      <td>Best current safe anchor: strong evolved/geodesic FTL and no trapped proxy.</td>
    </tr>
    <tr>
      <td><code>eval_000055</code></td>
      <td align="right"><code>1731.33</code></td>
      <td align="right"><code>0.06554</code></td>
      <td align="right"><code>0.05071</code></td>
      <td align="right"><code>5/5</code></td>
      <td align="right"><code>-0.05234</code></td>
      <td align="right"><code>-1.0</code></td>
      <td><code>1.34e-4</code> / <code>1.08e-5</code></td>
      <td align="right"><code>0.008</code></td>
      <td>Strong trapped mechanism: real geodesic/evolved FTL, but fails horizon safety and stability.</td>
    </tr>
    <tr>
      <td><code>eval_000024</code></td>
      <td align="right"><code>1283.36</code></td>
      <td align="right"><code>0.04773</code></td>
      <td align="right"><code>0.04230</code></td>
      <td align="right"><code>4/5</code></td>
      <td align="right"><code>-0.02666</code></td>
      <td align="right"><code>-1.0</code></td>
      <td><code>1.57e-4</code> / <code>1.34e-5</code></td>
      <td align="right"><code>0.439</code></td>
      <td>Secondary trapped mechanism: real FTL signal, cleaner constraints, not safe.</td>
    </tr>
    <tr>
      <td><code>000314</code></td>
      <td align="right"><code>816.36</code></td>
      <td align="right"><code>0.02718</code></td>
      <td align="right"><code>0.02646</code></td>
      <td align="right"><code>5/5</code></td>
      <td align="right"><code>0.00904</code></td>
      <td align="right"><code>0.0</code></td>
      <td><code>6.66e-5</code> / <code>6.23e-6</code></td>
      <td align="right"><code>0.860</code></td>
      <td>Horizon-safe but below the working evolved <code>f_op &gt; 0.03</code> gate.</td>
    </tr>
    <tr>
      <td><code>000358</code></td>
      <td align="right"><code>299.90</code></td>
      <td align="right"><code>0.03108</code></td>
      <td align="right"><code>0.0</code></td>
      <td align="right"><code>0/5</code></td>
      <td align="right"><code>0.00904</code></td>
      <td align="right"><code>0.0</code></td>
      <td><code>1.21e-4</code> / <code>9.57e-6</code></td>
      <td align="right"><code>0.841</code></td>
      <td>Horizon-safe coordinate channel, but geodesic check fails.</td>
    </tr>
    <tr>
      <td><code>eval_000025</code></td>
      <td align="right"><code>141.25</code></td>
      <td align="right"><code>0.05466</code></td>
      <td align="right"><code>0.0</code></td>
      <td align="right"><code>0/5</code></td>
      <td align="right"><code>-0.05051</code></td>
      <td align="right"><code>-1.0</code></td>
      <td><code>2.90e-4</code> / <code>2.96e-5</code></td>
      <td align="right"><code>0.008</code></td>
      <td>Reject as discovery target: coordinate FTL only, trapped and unstable.</td>
    </tr>
  </tbody>
</table>

Critical ranking after these promotions:

- Best current safe anchor: `000114`. It is the only HQ result here with strong evolved/geodesic FTL and no trapped-surface proxy.
- Best trapped high-FTL mechanism: `eval_000055`. It has strong `f_op` and `f_geo`, but fails horizon safety.
- Secondary trapped mechanism: `eval_000024`. It has real geodesic/evolved FTL, cleaner constraints, but weaker FTL than `000055` and still fails horizon safety.
- Reject as discovery target: `eval_000025`. The evolved coordinate FTL does not survive the geodesic cross-check.

Current notable candidates:

- `eval_000055`: promoted to HQ and confirmed as a strong trapped geodesic/evolved FTL mechanism (`f_op = 0.06554`, `f_geo = 0.05071`, `horizon_penalty = -1.0`).
- `eval_000025`: promoted to HQ and downgraded; evolved coordinate FTL is present, but geodesic FTL fails (`f_geo = 0`, `0/5` rays reached).
- `eval_000024`: promoted to HQ and confirmed as a persistent evolved/geodesic FTL channel, but no longer horizon-safe at HQ (`min_theta_plus = -0.02666`, `horizon_penalty = -1.0`).
- `eval_000001`: seeded survivor with high score (`1784`) and no horizon penalty in the low-res QD run.

### Next restart point

Use `000114` as the confirmed HQ safe anchor. Treat `eval_000055` and `eval_000024` as trapped high-FTL mechanism references, not safe targets. Do not prioritize `eval_000025` except as a cautionary example of coordinate FTL that fails the geodesic check. The next science target is to find/promote horizon-safe channel elites, reduce exotic matter, and repeat the HQ validation at higher resolution/longer time.
