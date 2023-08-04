# GRTeclyn Tests

The tests in this repository are implemented using the [Catch2
framework](https://github.com/catchorg/Catch2). If you want to just build and
run the tests make sure you have set up the git submodule correctly as described 
[below](#git-submodule) and then skip to the section on [Building and running
the tests](#building-and-running-the-tests).

## Catch2

### git submodule

The Catch2 code is contained within a [git
submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in the
[root/External/Catch2](../External/Catch2/) directory. If this repository was
cloned without using the `--recurse-submodule` flag, then you will need to first
initialize git submodule using the command

```bash
git submodule init
```

from somewhere within the repository. Then, in order to obtain the Catch2 code
at the relevant commit, run the command

```bash
git submodule update
```

The above command can also be used to update the checked out version of Catch2
to the version required by the current commit of this code. This is necessary
when the version of Catch2 is updated in the upstream repository (i.e. on
GitHub), and you fetch/pull this into your local clone.

### Basic Catch2 syntax

A simple Catch2 test case is created by adding code such as

```cpp
TEST_CAST("<test case name>")
{
    <test code>

    CHECK(<condition that is true if test passes and false if test fail>);    
}
```

and then making sure that this code is included in/linked with the main
[Tests.cpp](Tests.cpp) file.  Instead of a `CHECK()` clause which will report
failure in the output and return value but not abort the application (thereby
allowing later tests to run), one can instead use a `REQUIRE()` clause which
will immediately abort the application if the relevant test fails (i.e. if the
condition is false). Unless, there is a good reason, one should default to using
`CHECK()` clauses.

A common scenario is that one wishes to compare a floating point value produced
by the code to a known correct value. In this case instead of using the
`CHECK()` or `REQUIRE()` macros, one can instead use [Catch2's floating point
matchers](https://github.com/catchorg/Catch2/blob/v3.3.2/docs/comparing-floating-point-numbers.md)
in a `CHECK_THAT()` or `REQUIRE_THAT()` clause, such as

```cpp
double tol = 1e-10;
CHECK_THAT(computed_value, Catch::Matchers::WithinAbs(correct_value, tol));
```

### Catch2 documentation

The Catch2 documentation can be found in the Catch2 repository
[here](https://github.com/catchorg/Catch2/blob/v3.3.2/docs/Readme.md).

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

Finally, you will need to make sure the Catch2 git submodule is set up correctly
as described [above](#git-submodule).

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
as an argument e.g.
```
./Tests3d.gnu.ex "CCZ4 RHS"
```
and multiple tests can be run by passing their names in a comma separated list
e.g.
```
./Tests3d.gnu.ex "CCZ4 RHS,Derivatives"
```

Catch2 provides a plethora of command line flags to customize the output which
can be found by passing the `-h` flag e.g.
```
./Tests3d.gnu.ex -h
```
Some particularly useful ones include
```
-s, --success        include successful tests in output
-a, --abort          abort at first failure
--list-tests         list all/matching test cases                                          
```


## Adding a new test

See the [CCZ4 RHS test](./CCZ4RHSTest/) as an example. Here is an outline of the
basic steps to adding a new test to the [Tests application](./Tests.cpp).

1. Create a new test directory with an appropriate name (e.g. `NewGRTeclynTest`).
   Make sure that the directory name ends with `Test` or `Tests` (so that Make
   can find it).
2. In that directory, create a cpp file with the appropriate name (e.g.
   `NewGRTeclynTest.cpp`). In that file make sure you include the Catch2
   and base AMReX headers (if you are using any AMReX classes).
   ```cpp
   // Catch2 header
   #include "catch_amalgamated.hpp"
   
   // AMReX includes
   #include "AMReX.H"
   ```
3. In the cpp file, create a new `TEST_CASE()` as described
   [above](#basic-catch2-syntax). Make sure to initialize and finalize AMReX if
   using any AMReX data structures (otherwise no memory will be allocated):
   ```cpp
    TEST_CASE("New GRTeclyn")
    {
        // MPI_COMM_WORLD defined in AMReX_ccse-mpi.H even when compiling 
        // without MPI
        amrex::Initialize(MPI_COMM_WORLD)

        bool test_passes = true;

        // This code should set test_passes = false if some failure happens
        <test code>

        CHECK(test_passes);
        
        amrex::Finalize();
    }
    ```
4. Create a `Make.package` file in the same directory with the following
   content: 
   ```makefile
   # Can omit the following line if you don't have any extra headers 
   GRTECLYN_CEXE_headers += <any headers you need for your test>

   GRTECLYN_CEXE_sources += NewGRTeclynTest.cpp <any other cpp files>
   ```
5. In the GNUmakefile, one may need to uncomment some lines of the form
   ```makefile
   # include $(AMREX_HOME)/Src/<component>/Make.package
   ```
   depending on which AMReX components are required.
6. Build and run the test as described [above](#building-and-running-the-tests).