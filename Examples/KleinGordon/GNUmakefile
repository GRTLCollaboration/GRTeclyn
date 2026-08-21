GRTECLYN_HOME = $(realpath ../..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

EBASE = KleinGordon

AMREX_HOME ?= $(realpath ../../../amrex)


include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

src_dirs := $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/Tagging

# Find all source and header files in the src_dirs
CEXE_sources += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.cpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.hpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.impl.hpp)))

# Example specific source files
CEXE_sources += KleinGordonLevel.cpp \
                Main_KleinGordon.cpp

# Example specific header files
CEXE_headers += DerivedVariables.hpp \
                DerivedVariables.impl.hpp \
                KleinGordonLevel.hpp \
                KleinGordonRHS.hpp \
                SimulationParameters.hpp \
                SineGordon.hpp \
                StateVariables.hpp \
                Wave.hpp

VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules
