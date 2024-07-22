GRTECLYN_HOME ?= $(realpath ../..)
include $(GRTECLYN_HOME)/Tools/GNUMake/Make.base

EBASE = BinaryBH

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

src_dirs := . \
            $(GRTECLYN_SOURCE)/utils \
            $(GRTECLYN_SOURCE)/simd  \
            $(GRTECLYN_SOURCE)/CCZ4  \
            $(GRTECLYN_SOURCE)/BoxUtils  \
            $(GRTECLYN_SOURCE)/GRTeclynCore  \
            $(GRTECLYN_SOURCE)/AMRInterpolator  \
            $(GRTECLYN_SOURCE)/TaggingCriteria  \
            $(GRTECLYN_SOURCE)/InitialConditions/BlackHoles \
            $(GRTECLYN_SOURCE)/BlackHoles

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
include $(src_pack)
CEXE_sources += $(GRTECLYN_CEXE_sources)
CEXE_headers += $(GRTECLYN_CEXE_headers)
VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules