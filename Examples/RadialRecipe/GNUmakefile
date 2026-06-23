GRTECLYN_HOME = $(realpath ../..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

USE_PARTICLES = FALSE
USE_RL ?= FALSE

AMREX_HOME ?= $(realpath ../../../amrex)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

GRTECLYN_HOME = $(realpath ../..)
GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

src_dirs := $(GRTECLYN_SOURCE)/AMRInterpolator \
            $(GRTECLYN_SOURCE)/CCZ4 \
            $(GRTECLYN_SOURCE)/Matter \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/Tagging

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
include $(src_pack)
CEXE_sources += $(GRTECLYN_CEXE_sources)
CEXE_headers += $(GRTECLYN_CEXE_headers)
VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs)

ifeq ($(USE_RL),TRUE)
  DEFINES += -DUSE_RL
  ZMQ_PREFIX ?= $(GRTECLYN_HOME)/local/zeromq
  INCLUDE_LOCATIONS += $(GRTECLYN_SOURCE)/GRTeclynCore/RL $(ZMQ_PREFIX)/include
  # -L must go in LDFLAGS here: Make.defs (already included above) consumes
  # LIBRARY_LOCATIONS into LDFLAGS before this block runs, so adding it there
  # would be too late.  LDFLAGS is expanded at the link rule, so this works.
  LDFLAGS += -L$(ZMQ_PREFIX)/lib -Wl,-rpath,$(ZMQ_PREFIX)/lib
  LIBRARIES += -lzmq
endif

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules
