# GRTeclyn Tests

The tests in this repository are implemented using the [doctest
framework](https://github.com/doctest/doctest). If you want to just build and
run the tests skip to the section on [Building and running the
tests](#building-and-running-the-tests).

## doctest

### Source

The doctest code is contained within a single header: [doctest.h](../External/doctest/doctest.h).

### Basic doctest syntax

A simple doctest test case is created by adding code such as

```cpp
TEST_CASE("<test case name>")
{
    <test code>

    CHECK(<condition that is true if test passes and false if test fail>);    
}
```

In practice, due to device kernel linking issues with
[HIP](https://github.com/GRTLCollaboration/GRTeclyn/issues/48) and
[SYCL](https://github.com/GRTLCollaboration/GRTeclyn/issues/46) we have instead
created each test case as function declared in a header and defined in a cpp
file:

```cpp
void run_my_test_case()
{
    <test code>

    CHECK(<condition that is true if test passes and false if test fail>);    
}
```

The header is included in [TestCases.hpp](./TestCases.hpp) and this function is
called inside a `TEST_CASE` there:

```cpp
TEST_CASE("<test case name>")
{
   run_my_test_case();
}
```


Instead of a `CHECK()` clause which will report failure in the output and return
value but not abort the application (thereby allowing later tests to run), one
can instead use a `REQUIRE()` clause which will immediately abort the
application if the relevant test fails (i.e. if the condition is false).
Unless, there is a good reason, one should default to using `CHECK()` clauses.

A common scenario is that one wishes to compare a floating point value produced
by the code to a known correct value. In this case, one can test equality with
the a `doctest::Approx` (see [this docs
page](https://github.com/doctest/doctest/blob/ae7a13539fb71f270b87eb2e874fbac80bc8dda2/doc/markdown/assertions.md#floating-point-comparisons))
object, for example
```cpp
double tol = 1e-10;
CHECK(computed_value == doctest::Approx(correct_value).epsilon(tol));
```

### doctest documentation

The doctest documentation can be found in the doctest repository
[here](https://bit.ly/doctest-docs).

## Building and running the tests

### Prerequisites

You will need the usual AMReX prerequisites:

* Git
* GNU Make >= 3.81
* Python >= 2.7
* A Unix-like environment with perl and sed commands
* C compiler with C99 support
* C++ compiler with C++17 support

Additionally, you will need to make sure you have cloned the [AMReX
repository](https://github.com/AMReX-Codes/amrex) locally. The makefile assumes
you have cloned it into the same directory you cloned this repository (i.e.
[../../amrex](../../amrex/) relative to where this `README.md` file is located)
but you can instead clone it elsewhere and set the `AMREX_HOME` environment
variable appropriately:
```bash 
export AMREX_HOME=/path/to/amrex
```

### Building

The tests are all contained in one AMReX application with its `main()` defined
in [`Tests.cpp`](./Tests.cpp). To build this application, simply run
```
make -j <number of build jobs>
```
Like any other AMReX GNU Make application, one can pass [AMReX configuration
options](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#id1)
such as `USE_OMP=TRUE` on the command line e.g.
```
make -j 4 USE_OMP=TRUE
```
or modify these directly in the [GNUmakefile](./GNUmakefile). Alternatively if
one wishes to use the same configuration options for all AMReX
applications, one can set configuration options in
[`${AMREX_HOME}/Tools/GNUMake/Make.local-pre`](../../amrex/Tools/GNUMake/Make.local-pre)
(this file will need to be created if it doesn't exist).

### Running

One can run all the tests by simply executing the created executable (e.g.
`Tests3d.gnu.ex`) without any arguments. To run a specific test, pass its name
as an argument to the `-dt-tc` flag:
```
./Tests3d.gnu.ex -dt-tc="CCZ4 RHS"
```
and multiple tests can be run by passing their names in a comma separated list
e.g.
```
./Tests3d.gnu.ex -dt-tc="CCZ4 RHS,Derivative*"
```

doctest provides a plethora of command line flags to customize the output which
can be found by passing the `-dt-h` flag e.g.
```
./Tests3d.gnu.ex -dt-h
```
Some particularly useful ones include
```
-dt-s,   --dt-success              include successful tests in output
-dt-aa,  --dt-abort-after=<int>    abort after <int> failed assertions
-dt-ltc, --dt-list-test-cases      list all/matching test cases
```
All doctest flags are prefixed with `-dt-` or `--dt-`. Any unprefixed arguments
may be passed as command line arguments to test cases so that one can pass
arguments to AMReX e.g.
```bash
./Tests3d.gnu.ex amrex.verbose=1 -dt-tc="CCZ4 RHS"
```

## Adding a new test

See the [CCZ4 RHS test](./CCZ4RHSTest/) as an example. Here is an outline of the
basic steps to adding a new test to the [Tests application](./Tests.cpp).

1. Create a new test directory with an appropriate name (e.g. `NewGRTeclynTest`).
   Make sure that the directory name ends with `Test` or `Tests` (so that Make
   can find it).
2. In that directory, create a cpp file with the appropriate name (e.g.
   `NewGRTeclynTest.cpp`). In that file make sure you include the doctest
   and base AMReX headers (if you are using any AMReX classes).
   ```cpp
   // Doctest header
   #include "doctest.h"

   // Provides doctest::cli_args
   #include "doctestCLIArgs.hpp"
   
   // AMReX includes
   #include "AMReX.H"
   ```
3. In the cpp file, create a new `run_my_new_grteclyn_test()` as described
   [above](#basic-doctest-syntax). Make sure to initialize and finalize AMReX if
   using any AMReX data structures (otherwise no memory will be allocated):
   ```cpp
    void run_my_new_grteclyn_test()
    {
        // doctest::cli_args stores the non -dt- command line args
        int amrex_argc    = doctest::cli_args.argc();
        char **amrex_argv = doctest::cli_args.argv();
        // MPI_COMM_WORLD defined in AMReX_ccse-mpi.H even when compiling 
        // without MPI
        amrex::Initialize(amrex_argc, amrex_argv, true, MPI_COMM_WORLD);
        {

            bool test_passes = true;

            // This code should set test_passes = false if some failure happens
            <test code>

            CHECK(test_passes);
        }
        amrex::Finalize();
    }
    ```
4. Create a hpp header file with the same filename stem (e.g.
   `NewGRTeclynTest.hpp`).  In that file, declare the test function:
   ```cpp
   void run_my_new_grteclyn_test();
   ```
   Don't forget to include a header guard!
5. Create a `Make.package` file in the same directory with the following
   content: 
   ```makefile
   GRTECLYN_CEXE_headers += NewGRTeclynTest.hpp <any other headers>

   GRTECLYN_CEXE_sources += NewGRTeclynTest.cpp <any other cpp files>
   ```
6. Include the header file in [TestCases.hpp](./TestCases.hpp), create the test
   case and call the test function:
   ```cpp
   #include "NewGRTeclynTest.hpp"

   ...

   TEST_CASE("New GRTeclyn")
   {
      run_my_new_grteclyn_test();
   }
   ```
7. In the GNUmakefile, one may need to uncomment some lines of the form
   ```makefile
   # include $(AMREX_HOME)/Src/<component>/Make.package
   ```
   depending on which AMReX subcomponents are required.
8. Build and run the test as described [above](#building-and-running-the-tests).