# Parameters guide

GRTeclyn uses the AMReX [`ParmParse`](https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse) system for runtime parameters. At start-up, AMReX reads the parameter file and any command-line overrides into a global table. The GRTeclyn wrapper `GRParmParse` provides access to this table: an unscoped parser can look up a full name, while a scoped parser automatically adds its prefix. For example,

```cpp
GRParmParse evolution_pp("evolution");
evolution_pp.get("sigma", sigma);
```

reads `evolution.sigma`. Required values are normally read with `get`, while `queryAdd` is used to supply a default and add the resolved value to the table.

## Where parameters are loaded

Parameters normally live in the class that uses them, often in a nested `params_t` structure. That class reads its values from `GRParmParse` when they are needed. Keeping the values close to the relevant code makes ownership and defaults easier to identify and avoids passing unrelated parameters throughout the code.

!!! note "For GRChombo users"

    GRTeclyn does not construct and pass around a single `sim_params` object. `SimulationParameters.hpp` coordinates parameter validation, but it is not a container holding every parameter in the simulation.

Each example's `SimulationParameters.hpp` calls the shared checks in `BaseParameterChecker` and the checks supplied by its problem-specific classes. These checks run during AMReX initialisation, after the parameter table has been populated but before the classes load the values for use in the simulation. A new parameter-owning class may provide a static `check_params()` function and register it here when early validation is useful. This can check allowed values and relationships between parameters, and can issue warnings or stop with a clear error message.

## Common scopes

Parameter names use dot-separated scopes. The main scopes include:

| Scope | Typical contents |
| --- | --- |
| `amr` | Cell counts, refinement levels, regridding, grid-box sizes and output intervals |
| `geometry` and `boundary` | Domain size and centre, periodicity and physical boundary conditions |
| `evolution` | Time step, stopping conditions, derivative order, dissipation and runtime checks |
| `ccz4` and `gauge` | Formulation and moving-puncture gauge parameters |
| `tagging` | Refinement thresholds and tagging controls |
| `puncture_tagging` and `extraction_tagging` | Separation and size of refinement levels |
| `weyl_extraction` and `puncture_tracking` | Gravitational-wave extraction and black-hole tracking |
| Example-specific scopes | For example `bh1` and `bh2` |

This is only a guide to the organisation. The example `params_*.txt` files document the available parameters, defaults and suggested values for each executable, and we have tried to add useful comments where the naming is unclear.

Note that for parameters that are common with AMReX, we keep the (sometimes ugly) AMReX naming, to avoid any confusion or duplication. Again the comments in the params files are there to help.

### Physical boundary conditions

For each non-periodic direction, `boundary.hi_condition` and `boundary.lo_condition` select one of the following options:

- `FIRST_ORDER_EXTRAPOLATION_BC` uses AMReX's `foextrap` condition, which copies the current value in the nearest interior cell into the exterior ghost cells. This is constant extrapolation in space, not a boundary value held fixed in time.
- `SOMMERFELD_BC` applies an outgoing-radiation condition to the RHS in the outer valid cells. Its asymptotic values come from the state-variable definitions and may be changed through `set_vars_asymptotic_values()`.
- `REFLECTIVE_BC` fills ghost cells using the parity assigned to each state variable.

Use `UNSET_BC` for periodic directions, where physical boundary conditions do not apply. A condition supplied for a periodic direction is ignored.

## TwoPunctures parameters

When GRTeclyn is built with TwoPunctures support, its initial-data parameters use the `two_punctures` scope. They include the bare or target masses, the plus and minus puncture offsets, momenta and spins, the initial lapse, spectral resolution and nonlinear solver settings. See the BinaryBH parameter files and the [TwoPunctures guide](two_punctures.md) for some helpful tips.

## Checking a parameter file without running

To validate a parameter file and exit before constructing the grid or evolving, place the GRTeclyn option after AMReX's `--` command-line separator:

```bash
./BinaryBH3d.gnu.ex params_test.txt -- -check_params
```

For an MPI build, use the same arguments after the executable, for example:

```bash
mpiexec -n 4 ./BinaryBH3d.gnu.MPI.ex params_test.txt -- -check_params
```

The `--` and `-check_params` are separate arguments. The checker loads the parameter table, applies registered defaults, reports warnings or errors, writes the resolved table to `parameters_and_version.txt`, and then exits without starting the simulation.

The checker prints a short summary stating whether any warnings were found. Warnings do not prevent a simulation from running, but you should understand why they can safely be ignored before continuing. The checks catch selected parameter errors but cannot guarantee that a simulation will run successfully, so inspect `parameters_and_version.txt` to confirm that the resolved values match your intentions.
