# Gauge Classes

GRTeclyn implements the gauge-independent evolution equations and the gauge conditions in separate classes. `CCZ4RHS` or `CCZ4RHSWithMatter` updates variables related to the induced metric and the extrinsic curvature, while a gauge class calculates the lapse and shift right-hand sides. In contrast to GRChombo, this lets each example select its gauge without making the gauge a template parameter of the CCZ4 or matter RHS class.

The two moving puncture gauge classes currently used by the examples are:

- `MovingPunctureGauge`, which evolves the standard Gamma driver auxiliary field $B^i$ and is used by the Binary Black Hole example, and
- `IntegratedMovingPunctureGauge`, which eliminates the evolving Gamma driver auxiliary field and is used by the Scalar Field example.

Both classes are templated over the derivative class and constructed with the grid spacing. For example:

```cpp
CCZ4RHS<FourthOrderDerivatives> ccz4_rhs(dx);
MovingPunctureGauge<FourthOrderDerivatives> gauge(dx);
```

## Standard moving puncture gauge

`MovingPunctureGauge` implements

$$
\begin{aligned}
\partial_t \alpha
&= a_\alpha \beta^j \partial_j \alpha
   - c_\alpha \alpha^p (K - 2\Theta), \\
\partial_t \beta^i
&= a_\beta \beta^j \partial_j \beta^i + F B^i, \\
\partial_t B^i
&= a_\beta \beta^j \partial_j B^i
   - a_\beta \beta^j \partial_j \hat{\Gamma}^i
   + (\partial_t\hat{\Gamma}^i)_{\rm CCZ4}
   - \eta(r) B^i.
\end{aligned}
$$

Here $B^i$ is an evolved auxiliary field. The term $(\partial_t\hat{\Gamma}^i)_{\rm CCZ4}$ in the final equation must already contain every vacuum and matter contribution.

The following runtime parameters set the corresponding coefficients in the equations above:

| Parameter | Symbol | Default |
| --- | --- | --- |
| `gauge.lapse_advec_coeff` | $a_\alpha$ | `1.0` |
| `gauge.lapse_coeff` | $c_\alpha$ | `2.0` |
| `gauge.lapse_power` | $p$ | `1.0` |
| `gauge.shift_advec_coeff` | $a_\beta$ | `0.0` |
| `gauge.shift_Gamma_coeff` | $F$ | `0.75` |
| `gauge.eta` | $\eta$ | `1.0` |

### Far-field eta cutoff

`MovingPunctureGauge` can optionally make the Gamma driver damping parameter spatially varying according to [Eq. (11) of Schnetter (2010)](https://arxiv.org/abs/1003.0859):

$$
\eta(r)=\eta_*\left[c_f\frac{R^2}{r^2+R^2}+(1-c_f)\right],
$$

where $\eta_*$ is set by `gauge.eta`, $r=|\boldsymbol{x}-\boldsymbol{x}_{\rm center}|$ with $\boldsymbol{x}_{\rm center}$ set by `geometry.center`, and $R$ is given by `gauge.eta_cutoff_radius`. Here, $c_f$ is `1` when `gauge.enable_eta_cutoff = true` and `0` otherwise. The cutoff is disabled by default, in which case $\eta(r)=\eta_*$ everywhere. When enabled, $\eta(r)$ is approximately $\eta_*$ near the centre and decays as $R^2/r^2$ far outside the scale $R$. `gauge.eta_cutoff_radius` defaults to `500.0` and must be positive.

The eta cutoff is implemented in the standard `MovingPunctureGauge` and is also applied in `IntegratedMovingPunctureGauge`.

## Integrated moving puncture gauge

The `IntegratedMovingPunctureGauge` class implements

$$
\begin{aligned}
\partial_t \alpha
&= a_\alpha \beta^j \partial_j \alpha
   - c_\alpha \alpha^p (K - 2\Theta), \\
\partial_t \beta^i
&= a_\beta \beta^j \partial_j \beta^i
   + F\hat{\Gamma}^i - \eta\beta^i - B_0^i.
\end{aligned}
$$

The second equation follows from the observation that $C^i\equiv B^i-\left(\hat{\Gamma}^i-\tfrac{\eta}{F}\beta^i\right)$ is conserved along the advective flow. Solving for $B^i$ and substituting it into the original shift condition gives the integrated Gamma driver above, with $B_0^i=-F C^i$ as the integration constant.

GRTeclyn fixes this constant by requiring the initial advective time derivative of the shift to vanish. The `IntegratedMovingPunctureGauge` class therefore stores

$$
B_0^i\equiv F\hat{\Gamma}^i(t_0)-\eta\beta^i(t_0)
$$

in `vars.B`. This explains the `-vars.B(i)` term in the integrated shift equation.

To initialise `vars.B`, `IntegratedMovingPunctureGauge` must be called after the conformal connection functions have been calculated. The Scalar Field example does this in `initData()` by calling `GammaCalculator` and then the gauge object's function-call operator:

```cpp
GammaCalculator<> gamma_calculator(dx);
IntegratedMovingPunctureGauge<FourthOrderDerivatives> gauge(dx);

gamma_calculator(ix, iy, iz, state);
gauge(ix, iy, iz, state);
```

`GammaCalculator` must run first so that the value of $\hat{\Gamma}^i(t_0)$ used in $B_0^i$ is consistent with the initial conformal metric. The gauge initialisation then evaluates $B_0^i=F\hat{\Gamma}^i(t_0)-\eta\beta^i(t_0)$ and stores it in `vars.B` before time evolution starts.

## Composing the right hand side

A vacuum example should calculate the right hand side in the following order:

1. call `ccz4_rhs.compute_chi_and_h_ij(...)`;
2. call `ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(...)`;
3. call `gauge.calculate_gauge_rhs(...)`; and
4. call `ccz4_rhs.apply_dissipation(...)`.

The first two calculations should remain in separate `ParallelFor` launches to control GPU register use. The gauge and dissipation calls may share a third launch. See [Performance optimisation](performance_optimisation.md) for more detail.

For an evolution with matter, add the Einstein sources before calculating the gauge RHS, because the standard Gamma driver uses the fully accumulated right hand side of the conformal connection functions. The complete order is:

1. calculate the two vacuum CCZ4 parts;
2. add the Einstein sources;
3. calculate the gauge RHS;
4. add the matter RHS; and
5. apply dissipation.

See [`ScalarFieldLevel::specificEvalRHS()`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/ScalarField/ScalarFieldLevel.cpp) for an example of this ordering and [Matter classes](matter_classes.md) for guidance on adding a new matter type.
