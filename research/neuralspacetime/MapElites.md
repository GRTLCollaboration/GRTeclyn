## MAP-Elites FTL Discovery Status

Status: **reset**. The previous QD/HQ campaign results were discarded.

### Why the reset

The apparent-horizon / expansion (`theta_plus`) diagnostic in the GRTeclyn
example levels measured the radius and radial direction from the coordinate
origin (the domain corner) instead of the physics `grid_center`. With the
domain at `[0, L]` and `grid_center = (L/2, L/2, L/2)`, the regularizing
`2*sqrt(chi)/r` term was evaluated at `r ~ |grid_center|` instead of `r ~ 0`,
so it collapsed near the center and drove `theta_plus` spuriously negative.
This produced false trapped-surface detections (`max_horizon_radius ~ |grid_center|`,
`min_theta_plus < 0` located at `r ~ |grid_center|`) and a `-1.0` horizon
penalty that vetoed otherwise-viable candidates and inflated the stability
violation. The whole "horizon-safe vs trapped" ranking from that campaign was
therefore unreliable, so all of its runs were deleted.

### Fix applied

`theta_plus` is now measured about `grid_center` in:

- `Examples/RadialRecipe/RadialRecipeLevel.cpp`
- `Examples/RotatingWormholeCollapse/SupportedWormholeLevel.cpp`
- `Examples/SupportedWormholeCollapse/SupportedWormholeLevel.cpp`

The `RadialRecipe` binary was rebuilt with the fix.

### Next

A fresh QD FTL-discovery campaign starts from scratch (no resume, no seed
trajectory) on the corrected binary. Results will be recorded here.
