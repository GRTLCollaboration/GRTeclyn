GRTECLYN_HOME = $(realpath ../..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

AMREX_HOME ?= $(realpath ../../../amrex)

DEBUG = FALSE

DIM = 3

USE_MPI = TRUE
USE_OMP = FALSE

USE_CUDA  = FALSE
USE_HIP   = FALSE
USE_SYCL  = FALSE

#No Fortran
BL_NO_FORT = TRUE
# No Fortran probin file
AMREX_NO_PROBINIT = TRUE

CXXSTD=c++17

COMP=intel-llvm

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

src_dirs := $(GRTECLYN_SOURCE)/utils \
            $(GRTECLYN_SOURCE)/simd  \
            $(GRTECLYN_SOURCE)/GRTeclynCore  






# Use a common temporary build directory for all GRTeclyn applications
TMP_BUILD_DIR ?= $(GRTECLYN_HOME)/tmp_build_dir


