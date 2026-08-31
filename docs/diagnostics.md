# Diagnostics

Diagnostics are quantities calculated from the current state, rather than variables which need to be evolved. Examples include the Hamiltonian constraint, the energy density and the Weyl scalar $\Psi_4$. AMReX calls these *derived variables*.

GRTeclyn does not keep a separate diagnostic state on the grid. A diagnostic is registered during setup and calculated from the evolution state only when something requests it. For example, a diagnostic selected with `amr.derive_plot_vars` is calculated when a plot file is written. A derived variable may also be requested directly by an extraction routine. This avoids storing and updating values which may be needed only occasionally.

## Using existing diagnostics

Reusable diagnostics provide a `set_up()` function. Call this from the example level's `variableSetUp()` after `state_variable_set_up()`. For example, the Scalar Field example registers the constraints and energy density with

```cpp
ScalarFieldConstraints::set_up(state_index);
ScalarFieldEnergyDensity::set_up(state_index);
```

The names registered by these classes can then be added to `amr.derive_plot_vars` in the parameter file. Setting

```text
amr.derive_plot_vars = ALL
```

outputs every registered diagnostic, although selecting only the variables needed will save calculation and output time.

## Adding a diagnostic

The [Klein Gordon derived variables](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/KleinGordon/DerivedVariables.hpp) provide a compact example, while [`Constraints`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/CCZ4/Constraints.hpp) shows how a reusable diagnostic can be packaged as a class.

To add your own diagnostic:

1. Define the function which calculates it, for example `calc_energy_density`. Take a look at how this is done in `DerivedVariables.hpp` and `DerivedVariables.impl.hpp` in the Klein Gordon example. The function is given the current state and fills the new diagnostic variable.

2. Add the diagnostic to AMReX's list of derived variables, `derive_lst`. This happens in the example level's `variableSetUp()` function. Give the diagnostic a name and pass in the calculation function from step 1. You must also tell AMReX which state variables the calculation needs; the existing examples show the usual setup.

3. Select the diagnostic in `amr.derive_plot_vars`. When `writePlotFile()` is called, AMReX checks this list, calculates the selected diagnostics and writes them to the plot file. Diagnostics which have not been selected are not calculated.

For a diagnostic that will be reused, follow the pattern used by `Constraints`, where the calculation and setup are grouped together in a class. Example-specific diagnostics can instead be registered directly in the example level, as in [`KleinGordonLevel::variableSetUp()`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/KleinGordon/KleinGordonLevel.cpp).

No call needs to be added to the RHS calculation or every time step. Plot-file output requests the selected diagnostics automatically. The particle interpolator can also request a registered diagnostic for an extraction; see [Extraction](extraction.md) for an example. In both cases it is calculated from the state at the time it is needed.
