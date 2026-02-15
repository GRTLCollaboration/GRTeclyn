# A Testing with the Klein Gordon example

The repository includes an example which solves the Klein Gordon equation in flat spacetime (i.e. without GR, or any non trivial metric background). It casts the equation as first order in time, and therefore only contains two variables - the field `phi` and its conjugate momentum `Pi`.

This can be useful if you:

- want a simpler example from which to learn the structure of the code and how it interacts with AMReX, before diving into the full CCZ4 equations
- want to work on optimizating the code (sometimes the register pressure may be too high on the PVCs for an analysis using the BinaryBH example)