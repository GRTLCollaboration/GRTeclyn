DEBUG = FALSE

USE_MPI = TRUE
USE_OMP = FALSE

COMP = gnu

DIM = 3

BL_NO_FORT = TRUE

AMREX_HOME ?= ../../../amrex

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

GRAMREX_SOURCE = ../../Source

src_dirs := $(GRAMREX_SOURCE)/utils \
            $(GRAMREX_SOURCE)/simd  \
            $(GRAMREX_SOURCE)/CCZ4  \
            $(GRAMREX_SOURCE)/BoxUtils  \
            $(GRAMREX_SOURCE)/GRChomboCore  \
            $(GRAMREX_SOURCE)/AMRInterpolator  \
            $(GRAMREX_SOURCE)/TaggingCriteria  \
            $(GRAMREX_SOURCE)/InitialConditions/BlackHoles \
            $(GRAMREX_SOURCE)/BlackHoles

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
include $(src_pack)
VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(src_dirs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
