# GRTeclyn Agents Guide

## Overview

GRTeclyn is a numerical relativity code for solving Einstein's equations using finite difference methods and adaptive mesh refinement methods. It is built on top of the [AMReX](https://github.com/AMReX-Codes/amrex) library. AMReX supplies the supporting functions for mesh refinement, tagging, load balancing, the data containers, parallelization and GPU offloading.

## Purpose & Personas

- **GRTeclyn developers** - You are an expert in GPU programming with AMReX. You follow a predictable, documented workflow that respects the coding conventions defined in the documentation in `docs/`.
- **GRTeclyn users** - You run, configure and post-process simulations but do not modify the code apart from the parameter files or to select compile time options. You will provide accurate build/test/docs guidance sourced from the repo.
- **Downstream module developers** — You are responsible for downstream repositories e.g. `GRFolres`, `GRDzhadzha`, that are built *on top of* GRTeclyn and pinned directly against its class/method names (`GRAmrLevel`, `StateVariables` enums, the `matter_t` interface `CCZ4RHSWithMatter` expects). As with GRTeclyn developers, but also be mindful of breaking GRTeclyn changes.


## Documentation

Refer to the documentation as an authoritative source. The documentation is written using `mkdocs` and is kept in the `docs` directory.
It is built automatically on pushes to `main` by a GitHub action. For local builds, `mkdocs` can be installed using `pip install -r docs/requirements.txt`
but prefer making a virtual environment (or using `pipx` if available) to keep the user's base env clean.

## Code structure

The code is structured as follows:

- Source: This is where the main code is located. It is not specific to a particular physical scenario. Nothing is built until the file has been included into an example; it is almost library like in effect. See the subdirectory map below.
- Tests: These are the unit tests run via [doctest](https://github.com/doctest/doctest). There is one unit test per directory, except for the `Common` directory which contains definitions used across tests. Create a new directory for testing new features here with the directory name ending in Tests or Test and add a `TEST_CASE` entry to `Tests/TestCases.hpp`. Refer to `Building and Running` in the documentation for more information.
- Examples: These are specific physics examples using the underlying code. BinaryBH runs a binary black hole merger. KerrBH evolves a single rotating black hole. ScalarField runs an oscillaton. KleinGordon solves the Klein-Gordon equation without GR. Refer to `Building and Running` for more information.
- Tools: Contains doxygen documentation (not used), the default top-level Makefiles required when compiling the examples (including the dependencies required by TwoPunctures), the tool to print out quotes following a compilation (StoicQuotes) and the copyright info of the repo.
- external: Contains vendored external dependencies such as `doctest`.
- docs: Documentation lives here.
- .github: Contains GitHub specific configuration files such as GitHub actions workflows responsible for running CI/CD on GitHub. Also contains comparison data for CI/CD.
- .gitlab: The repository is also mirrored to a GitLab instance in order to run the GitLab pipeline for the unit and regression tests using GPUs on the A100 partition on the CSD3 HPC system. GitLab specific files reside here.

### Source Sub-directories
 - `GRTeclynCore` — the boundary between GRTeclyn and AMReX. The derived classes of
    the AMReX Amr and AmrLevel classes, the hooks an example level overrides, and
    parameter and setup infrastructure. Nothing here should be specific to an
    evolution system or a spacetime.
 - `Grids` — everything that knows about the mesh: the tensor and dimension
    types the pointwise code is written in, a cell's coordinates, finite
    difference stencils that reach into neighbouring cells, and boundary and
    parity handling.
 - `Maths` — mathematical operations that are pure functions of local values,
    with no knowledge of the mesh or of AMReX. If it needs a neighbouring cell
    it belongs in `Grids` instead.
 - `CCZ4` — the vacuum evolution system: the equations themselves, the gauge
    choices, the constraints, and diagnostics derived from the metric.
 - `Matter` — the matter-coupled counterparts of the `CCZ4` pieces, and the
    matter interface that a matter class must satisfy. Kept separate so that
    vacuum examples e.g. the binary black hole example need not include it.
 - `BlackHoles` — physics specific to black hole spacetimes: initial data
    and puncture tracking. It lives here rather than in an example directory because
    more than one example could use it.
 - `Tagging` — refinement criteria, kept together because they are
    interchangeable policies that an example selects between.
 - `ParticleInterpolator` — This provides the infrastructure to interpolate state
    and derived quantities at arbitrary points within the domain. It is built
    upon the AMReX particle infrastructure. There are also classes built on top
    of the ParticleInterpolator capabilities to interpolate quantities over user
    specified surfaces e.g. spheres and lines.
 - `IO` — parameter parsing and the small text outputs GRTeclyn writes itself,
    as distinct from the plotfiles and checkpoints AMReX writes.

  See `docs/code_structure.md` for the full discussion of the
  AMReX/GRTeclyn/example hierarchy.


## Building and Running

### Tests

The tests are compiled into a single binary but can be run individually using the `doctest` backend by passing the flag `-dt-tc="<name of test>"`.
All `doctest` command line argument flags are prefaced with `-dt-` e.g. `-dt-help` to list all options and `-dt-ltc` to list the names of all test cases.
To make and run the unit tests:

```
make run <AMReX build arguments>
```
from the top level `GRTeclyn` directory.

Notes:
- You may need to add the build argument `COMP`, to specify a toolchain. The default is `gnu` if nothing is specified. An appropriate choice depends on system you are running on and what modules have been loaded. Check the documentation in `docs/example_configs.md` to see if there are recommended setups for the user's environment. This is an AMReX setting.
- There are three GPU build options: `USE_CUDA`, `USE_HIP` and `USE_SYCL`. Set only one of these to be true per build as appropriate to the target system's GPU architecture. It may not be necessary to set `COMP` at the same time; check the recommended setups for the user's environment in the documentation. By default these are set to `FALSE`. This is an AMReX setting.
- Some tests need to be built against HDF5 using the build argument `USE_HDF5=TRUE` as `h5diff` is used to compare with trusted outputs either from GRChombo or a tagged version of GRTeclyn. If `USE_HDF5=FALSE` is set (the default option) then the tests requiring HDF5 will be skipped. If both `USE_MPI=TRUE` and `USE_HDF5=TRUE` are set, the HDF5 library will need to have been built with parallel HDF5 enabled.
- The CI running on GitHub Actions has a test matrix for the values of `DEBUG`, `USE_MPI` and `USE_OMP`. Any of these values can be set to either (`TRUE` or `FALSE`). `COMP` can also be varied between `gnu` and `llvm` but this will depend on what compilers you can access. It is not necessary to test all combinations locally as these will be tested in CI on push to the remote if the conditions are met in `.github/workflows/check-changes.yml`.
- If on a login node for a HPC cluster, do not attempt to run the tests directly. Prefer to submit a job to the job scheduler to run the tests on a compute node. Some modules may need to be loaded depending on the target partition. Check `docs/example_configs.md` to see if there are modules and build configs for the system that you are on. If in doubt, ask the user directly.

#### Adding a Unit Test

A unit test should be added to accompany each new addition to the code that covers its typical use case, preferably with a known result from e.g. analytics, GRChombo, or an external calculation, for comparison. Initial data for the CCZ4 and BSSN variables, the scalar field, $phi$, and its first derivative, $Pi$, can be found in the `Tests/Common` directory. It may be necessary to setup an isolated simulation with no RHS solve to use the `AmrLevel` object as in e.g. `PunctureTrackerTest` and `SphericalExtractionTest`. New tests can be added by creating a new subdirectory in the `Tests` directory with the word `Test` appended to the end of the class/function name and adding a new `TEST_CASE` to `Tests/TestCases.hpp` in which the test function is called. The test function should live in a `.cpp` source file that calls `amrex::Initialize` at the beginning and `amrex::Finalize` at the end and include the `doctest` headers, `doctest.h` and `doctestCLIArgs.hpp`. The header file should only have the test function declaration and header guard. For example a class called `MyClass` should have a subdirectory inside `Tests` called `MyClassTest` and call the test using `TEST_CASE("MyClass"){run_my_class_test();}` from `Tests/TestCases.hpp`. `run_my_class_test` needs to be declared in `Tests/MyClassTest/MyClassTest.hpp` and defined in `Tests/MyClassTest/MyClassTest.cpp`.

### Examples

Each subdirectory in the `Examples` directory contains one fully self-contained example that can each be built and run separately. Use `make` to build with optional flags for GPU builds as necessary (refer to the documentation `docs/building_gpus.md` for specific flags). The examples can be run using `mpirun` or `mpiexec`, e.g. `mpirun -np 2 ./BinaryBH3d.<build options>.ex ./params_test.txt` if built with MPI (enabled by default). Otherwise run by calling the binary directly. Adjust the number of ranks to account for locally available resources.

Notes:
- Take note of the first two points from the Tests section as they are still valid AMReX options when building the examples.
- Check the requirements in the documentation. By default `USE_MPI=TRUE` so you either need access to an MPI library or set this to false manually.
- When adding a new example to solve a specific physics problem, prefer to create a child class for the `AmrLevel` object named as `SpecificPhysicsLevel` with `GRAmrLevel` as the parent class.
- The BinaryBH, KerrBH and Klein-Gordon examples also have a regression test with `params_test.txt` as the parameter file for each directory. Test it manually by building `fcompare.<build config>.ex` in the `${AMREX_HOME}/Tools/Plotfile` directory and run it against the outputs stored in `.github/workflows/data/`. For example: `${AMREX_HOME}/Tools/Plotfile/fcompare.gnu.ex --abs_tol 1e-10 --rel_tol 1e-10 ../../.github/workflows/data/BinaryBH/plt00008_compare plots/plt00008`.
- `params_test.txt` should be small enough to be run on a local machine without MPI. Prefer this when running with limited resources e.g. on a laptop.

#### Adding a Regression Test
Adding a new example also requires a regression test to cover that example. The tests involve running a small simulation for a few timesteps (currently hardcoded to 8) and comparing the output with known data stored in `.github/workflows/data`. The regression tests are run remotely on the GitHub server and on Nvidia A100 GPUs using a GitLab instance hosted on CSD3. All materials related to the GitHub regression testing is stored in `.github/workflows` and the GitLab regression test is contained in `.gitlab-ci.yml` with additional helper scripts/files in `.gitlab`. For the GitHub action, common actions such as checking out the repositories and installing dependencies are kept in `.github/workflows/regression-tests-setup/action.yml`. These steps should be reused as much as possible to avoid duplication. The new example should be added to the matrix of test cases in `.github/workflows/regression-tests.yml`, which uses the example template, `.github/workflows/regression-tests-examples-template.yml` to run the code and compare the output with an earlier trusted result using the `fcompare` tool supplied by AMReX. It will also test particle positions (if any) and whether the code can restart from a checkpoint cleanly. For the GitLab pipeline, you will need to add steps to build the example in `.gitlab-ci.yml` and add another line to the Slurm script, `.gitlab/run-examples.sh` to run the new example. In both instances, the test data are stored in `.github/workflows/data`. New data may be added here for comparison with `fcompare`.

## Recommended workflow

When setting up new users to the code or on a new system/HPC cluster:

- Install [pre-commit](https://pre-commit.com/) according to the docs.
- Always format code with clang-format. Follow the rules in .clang-format.
- Git clone the `development` branch of AMReX if not already present on user's system.
- Git clone GRTeclyn if not already present on user's system. Prefer the GRTeclyn directory to live next to the AMReX directory not as a sub-directory or vice versa.
- Fetch the latest version before making a new branch for changes, preferably from `main` if appropriate.
- Refer to the coding conventions in the documentation.
- Linters are used to enforce good coding practices. There is a C++ linter contained in `.github/workflows/lint.yml` that is triggered by pushes and pull requests to the GitHub repo. Refer to the workflow for the version number and the `.clang-tidy` file for the rules that we bypass. A commonly flagged issue by the linter is `bugprone-easily-swappable-parameters` which is acceptable to bypass. The unit test `.cpp` file needs `// NOLINTNEXTLINE(bugprone-casting-through-void)` before `amrex::Initialize` for Open MPI builds. There are also linters that are pre-commit hooks e.g. `yamllint` and `shfmt` for `yaml` and shell scripts respectively.
- Prefer to manually run the unit and regression tests after each commit. Use git diff to see whether a unit test run is appropriate, e.g. run the puncture tracker test if changes to PunctureTracker.impl.hpp, PunctureTracker.hpp, or any of the files in Tests/PunctureTrackerTest are made. If the changes involve the entire example pipeline or are not covered by the unit test, run the regression test. Note that passing the unit and regression tests does not mean that the code is bug free as the tests do not have 100% coverage.


## Rules for pull requests

- Never push directly to `main`. This is protected by branch protection rules also.
- New feature branches should be named `feature/short_description_of_feature`.
- Branches containing code improvements should be named `enhancement/short_description_of_enhancement`.
- Branches for simulations, e.g. to record config files for later reproduction or sharing, should be named `simulation/short_description_of_simulation`.
- Bugfixes should be named `bugfix/github_issue_number_or_name`. If no relevant issue number exists, create an issue.
- Try to use a descriptive title written in the imperative mood.
- Require that all the CI/CD tests have passed, except for draft pull requests.
- Draft PRs can be opened at any time, but may invite comments from other users!
- Before opening a pull request, fetch from the remote to ensure that the branch is based off an up-to-date `main`. Rebase if necessary.
- Small commits e.g. changes for the linter, should be squashed or fixed up.
- Remind the user to update the documentation with the pull request. Always update the documentation if changes made by the pull request also update information in the documentation. Always update the documentation if changes made by the pull request are breaking changes.
- Do not add or track changes to the `Make.local-pre` file. Prefer not to make changes to the top level `GNUMakefile`, example specific changes can be made in the `GNUMakefile` in that example directory.
- Prefer to keep pull requests atomic. If the pull request addresses more than one issue, prefer to open a separate pull request for each.

## Rules for Git commit messages

- Try to use a descriptive title with less than 50 characters written in the imperative mood.
- The body of the commit message should describe the changes made and the reasoning behind them. Preferably the text should be wrapped at 72 characters.

## Rules for opening a new issue

- Check if the branch is up to date with `main`.
- Check if there is an existing open issue that covers the same topic.
- Write a detailed report on how to reproduce the issue, including the build environment and flags used, and quote the GitHub SHA.


## GPU Lambda Safety

Identical caveat to AMReX's own `AGENTS.md`: AMREX_GPU_DEVICE/ParallelFor lambdas run on the GPU, so never capture host-only pointers (e.g., Geometry::CellSize(), Geometry::ProbLo()). Take the device-safe views first—e.g., auto const dx = geom.CellSizeArray();, auto const problo = geom.ProbLoArray();—and pass those by value into the lambda.