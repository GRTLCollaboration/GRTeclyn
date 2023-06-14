GRAMREX_HOME = $(realpath ../..)

include $(GRAMREX_HOME)/Tools/GNUMake/Make.defaults

AMREX_HOME ?= $(realpath ../../../amrex)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

GRAMREX_HOME = $(realpath ../..)
GRAMREX_SOURCE = $(GRAMREX_HOME)/Source

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
CEXE_sources += $(GRAMREX_CEXE_sources)
CEXE_headers += $(GRAMREX_CEXE_headers)
VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs)

include $(GRAMREX_HOME)/Tools/GNUMake/Make.rules