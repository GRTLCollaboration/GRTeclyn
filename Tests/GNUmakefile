GRTECLYN_HOME = $(realpath ..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

# We don't want to use MPI by default for the tests
USE_MPI = FALSE

EBASE = Tests

AMREX_HOME ?= $(realpath ../../amrex)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

src_dirs := $(GRTECLYN_SOURCE)/utils \
            $(GRTECLYN_SOURCE)/CCZ4  \
            $(GRTECLYN_SOURCE)/BoxUtils  \
            $(GRTECLYN_SOURCE)/simd

GRTECLYN_TESTS_HOME = $(GRTECLYN_HOME)/Tests

# test_dirs = $(GRTECLYN_TESTS_HOME)/Common \
#             $(GRTECLYN_TESTS_HOME)/CCZ4RHSTest \
#             $(GRTECLYN_TESTS_HOME)/DerivativeUnitTests \
#             $(GRTECLYN_TESTS_HOME)/PositiveChiAndAlphaUnitTest \
#             $(GRTECLYN_TESTS_HOME)/CCZ4GeometryUnitTests \
#             $(GRTECLYN_TESTS_HOME)/CoordinateTransformationsTest \
#             $(GRTECLYN_TESTS_HOME)/SphericalHarmonicTest

# Find all directories that end in "Test" or "Tests"
# Note that hidden directories (i.e. those starting with a '.') are not included
test_dirs_match := $(wildcard $(GRTECLYN_TESTS_HOME)/*Test/ $(GRTECLYN_TESTS_HOME)/*Tests/)

# Strip the trailing path separator (probably not necessary)
test_dirs := $(patsubst %/, %, $(test_dirs_match))

# Finally add in the Common directory
test_dirs += $(GRTECLYN_TESTS_HOME)/Common

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
test_pack := $(foreach dir, $(test_dirs), $(dir)/Make.package)
include $(src_pack)
include $(test_pack)

doctest_dir := $(GRTECLYN_HOME)/External/doctest
doctest_headers := doctest.h

CEXE_sources += $(GRTECLYN_CEXE_sources)
CEXE_headers += $(GRTECLYN_CEXE_headers) $(doctest_headers)
VPATH_LOCATIONS += $(src_dirs) $(test_dirs)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs) $(test_dirs)
SYSTEM_INCLUDE_LOCATIONS += $(doctest_dir)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_MPI), TRUE)
LAUNCH = mpiexec -n 2 ./
else
LAUNCH = ./
endif

DOCTEST_OPTIONS := -dt-d=yes

run: $(executable)
	$(info Running tests for configuration: $(optionsSuffix))
	$(SILENT) $(LAUNCH)$(executable) $(DOCTEST_OPTIONS)

