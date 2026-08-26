# Matter classes

GRTeclyn separates the evolution of the spacetime variables from the details of a particular matter model. The `CCZ4RHSWithMatter` class inherits the vacuum CCZ4 calculations from `CCZ4RHS`, then adds the matter sources to the spacetime equations and calculates the right hand side (RHS) for the matter variables themselves.

The matter-aware RHS code is templated over the matter class and the derivative class. This means that the compiler creates a version of the RHS for the particular types that an example selects. There is no runtime lookup: inside the generic code, an object of type `matter_t` can be used as if its exact class had been written there directly.

For example, if a template is declared schematically as

```cpp
template <class matter_t, class deriv_t = FourthOrderDerivatives>
class MatterEvolution
{
    matter_t m_matter;
    deriv_t m_deriv;
};
```

then `matter_t` and `deriv_t` are placeholders. Supplying a scalar field class and `FourthOrderDerivatives` makes those the concrete member types. A default such as `FourthOrderDerivatives` means that the second type need only be supplied when a different derivative operator is wanted. Since the type is needed at compile time, template implementations generally live in header files; in GRTeclyn the longer definitions are often placed in an `.impl.hpp` file which is included at the end of the main `.hpp` file.

The [ScalarField class](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/Matter/ScalarField.hpp) and its [implementation](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/Matter/ScalarField.impl.hpp) provide the main example to follow when adding a matter type.

## Adding a new matter type

First decide which new variables need to be evolved. Add their indices, names, boundary parities and asymptotic values to the example's `StateVariables.hpp`. It is also useful to make a small variables class equivalent to [`ScalarFieldVars`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/Matter/ScalarFieldVars.hpp). This gives readable accessors for the CCZ4 and matter variables at one grid cell.

Next create header and implementation files equivalent to `ScalarField.hpp` and `ScalarField.impl.hpp`. The matter class should provide:

- a `Vars` alias for its variables class;

- a default constructor which reads any runtime parameters needed by the matter model (this is required when the class is used by the [diagnostic callbacks](diagnostics.md));

- `compute_emtensor(...)`, which returns an `emtensor_t` containing the physical energy density $\rho$, momentum density $j_i$, spatial stress $S_{ij}$ and its trace $S$;

- `compute_einstein_sources(...)`, which returns the corresponding `einstein_sources_t` used in the Einstein and gauge RHS equations; and

- `add_matter_rhs(...)`, which adds the evolution equations for the matter variables to the RHS array.

The distinction between the two source functions is important. `compute_emtensor(...)` returns the 3+1 decomposition of the physical stress-energy tensor $T_{ab}$. These are the quantities normally wanted for diagnostics, for example when outputting the energy density. `compute_einstein_sources(...)` returns the quantities that appear in the Einstein equations, including the gravitational coupling. For minimally coupled matter this will usually call `compute_emtensor(...)` and multiply each component by $8\pi G$. Keeping the two functions separate prevents a diagnostic intended to output $\rho$ from accidentally outputting $8\pi G\rho$ instead.

Finally, in the example level class:

1. Select the concrete matter type, as `ScalarFieldLevel` does with `ScalarField<Potential>`.
2. Construct the matter-aware CCZ4 RHS with that matter type and the required derivative type.
3. During `specificEvalRHS()`, call the inherited vacuum CCZ4 calculations, add the Einstein sources, add the matter RHS, and then apply dissipation. See [`ScalarFieldLevel::specificEvalRHS()`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/ScalarField/ScalarFieldLevel.cpp) for the order of these operations.
4. Register any matter diagnostics, such as the energy density or constraints with matter, in `variableSetUp()`.

If the new model has a potential or another interchangeable piece of physics, that may itself be templated in the same way. The Scalar Field example keeps its model-specific potential in [`Examples/ScalarField/Potential.hpp`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/ScalarField/Potential.hpp), while the reusable scalar field evolution remains in `Source/Matter`.
