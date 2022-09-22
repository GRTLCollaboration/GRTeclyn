DEBUG = FALSE

USE_MPI = TRUE
USE_OMP = FALSE

USE_HDF5 = FALSE

TINY_PROFILE = FALSE

TEST = FALSE
USE_ASSERTION = FALSE
BOUND_CHECK = FALSE

USE_CUDA  = FALSE
USE_HIP   = FALSE
USE_DPCPP = FALSE

COMP = gnu

DIM = 3

BL_NO_FORT = TRUE
AMREX_NO_PROBINIT = TRUE

CXXSTD = c++17

AMREX_HOME ?= ../../../amrex

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

GRAMREX_SOURCE = ../../Source

src_dirs := $(GRAMREX_SOURCE)/utils \
            $(GRAMREX_SOURCE)/simd  \
            $(GRAMREX_SOURCE)/CCZ4  \
            $(GRAMREX_SOURCE)/BoxUtils  \
            $(GRAMREX_SOURCE)/GRAMReXCore  \
            $(GRAMREX_SOURCE)/AMRInterpolator  \
            $(GRAMREX_SOURCE)/TaggingCriteria  \
            $(GRAMREX_SOURCE)/InitialConditions/BlackHoles \
            $(GRAMREX_SOURCE)/BlackHoles

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
include $(src_pack)
VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(src_dirs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
