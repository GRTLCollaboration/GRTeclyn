GRTECLYN_HOME = $(realpath ..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

# We don't want to use MPI by default for the tests
USE_MPI = FALSE
USE_PARTICLES = TRUE

EBASE = Tests

AMREX_HOME ?= $(realpath ../../amrex)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

GRTECLYN_TESTS_HOME = $(GRTECLYN_HOME)/Tests

src_dirs := $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/BlackHoles \
            $(GRTECLYN_SOURCE)/CCZ4 \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/Matter \
            $(GRTECLYN_SOURCE)/Tagging

# Find all directories that end in "Test" or "Tests"
# Note that hidden directories (i.e. those starting with a '.') are not included
test_dirs_match := $(wildcard $(GRTECLYN_TESTS_HOME)/*Test/ $(GRTECLYN_TESTS_HOME)/*Tests/)

# Strip the trailing path separator (probably not necessary)
test_dirs := $(patsubst %/, %, $(test_dirs_match))

# Finally add in the Common directory
test_dirs += $(GRTECLYN_TESTS_HOME)/Common

search_dirs := . $(src_dirs) $(test_dirs)

doctest_dir := $(GRTECLYN_HOME)/external/doctest
doctest_headers := doctest.h

# Find all source and header files
CEXE_sources += $(notdir $(foreach dir, $(search_dirs), $(wildcard $(dir)/*.cpp)))
CEXE_headers += $(notdir $(foreach dir, $(search_dirs), $(wildcard $(dir)/*.hpp)))
CEXE_headers += $(notdir $(foreach dir, $(search_dirs), $(wildcard $(dir)/*.impl.hpp)))

# Add doctest header
CEXE_headers += $(doctest_headers)

VPATH_LOCATIONS += $(src_dirs) $(test_dirs)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs) $(test_dirs)
SYSTEM_INCLUDE_LOCATIONS += $(doctest_dir)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_MPI), TRUE)
  ifeq ($(origin SLURM_NTASKS),undefined)
    LAUNCH = mpiexec -n 2 ./
  else
    LAUNCH = mpiexec -n ${SLURM_NTASKS} ./
  endif
else
LAUNCH = ./
endif

DOCTEST_OPTIONS := -dt-d=yes

run: $(executable)
	$(info Running tests for configuration: $(optionsSuffix))
	$(SILENT) $(LAUNCH)$(executable) $(DOCTEST_OPTIONS)
