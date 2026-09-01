GRTECLYN_HOME = $(realpath ../..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

EBASE = BinaryBH

# We use particles for the puncture tracking
USE_PARTICLES = TRUE

AMREX_HOME ?= $(realpath ../../../amrex)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

GRTECLYN_HOME = $(realpath ../..)
GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

# Name all directories for relevant source and header files
# Don't forget the dot! (i.e., example specific files are in current directory)
src_dirs := $(GRTECLYN_SOURCE)/BlackHoles \
            $(GRTECLYN_SOURCE)/CCZ4 \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/Tagging \
            $(realpath .)

# Find all source and header files in the src_dirs
CEXE_sources += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.cpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.hpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.impl.hpp)))

VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(src_dirs)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules
