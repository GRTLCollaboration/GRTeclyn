DEBUG = FALSE

USE_MPI = FALSE
USE_OMP = FALSE

USE_HDF5 = FALSE

TINY_PROFILE = FALSE

TEST = FALSE
USE_ASSERTION = FALSE
BOUND_CHECK = FALSE

USE_CUDA  = FALSE
USE_HIP   = FALSE
USE_SYCL = FALSE

COMP = gnu

DIM = 3

BL_NO_FORT = TRUE
AMREX_NO_PROBINIT = TRUE

CXXSTD = c++17

EBASE = Tests

AMREX_HOME ?= $(realpath ../../amrex)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
# include $(AMREX_HOME)/Src/Boundary/Make.package
# include $(AMREX_HOME)/Src/AmrCore/Make.package
# include $(AMREX_HOME)/Src/Amr/Make.package

GRAMREX_HOME = $(realpath ..)
GRAMREX_SOURCE = $(GRAMREX_HOME)/Source

src_dirs := $(GRAMREX_SOURCE)/utils \
            $(GRAMREX_SOURCE)/CCZ4  \
            $(GRAMREX_SOURCE)/BoxUtils  \
            $(GRAMREX_SOURCE)/simd

GRAMREX_TESTS_HOME = $(GRAMREX_HOME)/Tests

# test_dirs = $(GRAMREX_TESTS_HOME)/Common \
#             $(GRAMREX_TESTS_HOME)/CCZ4RHSTest \
#             $(GRAMREX_TESTS_HOME)/DerivativeUnitTests \
#             $(GRAMREX_TESTS_HOME)/PositiveChiAndAlphaUnitTest \
#             $(GRAMREX_TESTS_HOME)/CCZ4GeometryUnitTests \
#             $(GRAMREX_TESTS_HOME)/CoordinateTransformationsTest \
#             $(GRAMREX_TESTS_HOME)/SphericalHarmonicTest

# Find all directories that end in "Test" or "Tests"
# Note that hidden directories (i.e. those starting with a '.') are not included
test_dirs_match := $(wildcard $(GRAMREX_TESTS_HOME)/*Test/ $(GRAMREX_TESTS_HOME)/*Tests/)

# Strip the trailing path separator (probably not necessary)
test_dirs := $(patsubst %/, %, $(test_dirs_match))

# Finally add in the Common directory
test_dirs += $(GRAMREX_TESTS_HOME)/Common

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
test_pack := $(foreach dir, $(test_dirs), $(dir)/Make.package)
include $(src_pack)
include $(test_pack)

catch2_dir := $(GRAMREX_HOME)/External/Catch2/extras
catch2_sources := catch_amalgamated.cpp
catch2_headers := catch_amalgamated.hpp

CEXE_sources += $(GRAMREX_CEXE_sources) $(catch2_sources)
CEXE_headers += $(GRAMREX_CEXE_headers) $(catch2_headers)
VPATH_LOCATIONS += $(src_dirs) $(test_dirs) $(catch2_dir)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs) $(test_dirs)
SYSTEM_INCLUDE_LOCATIONS += $(catch2_dir)

# We want to define our own main rather than using Catch2's main
DEFINES += -DCATCH_AMALGAMATED_CUSTOM_MAIN

include $(GRAMREX_HOME)/Tools/GNUMake/Make.rules