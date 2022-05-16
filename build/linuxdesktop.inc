# BUILD INCLUDE build/linuxdesktop.inc
# Authors: P.D. Boom
# ------------------------------------------------------------------------------
# Include file for standard linux desktop workstation compilation using gfortran
# and OpenMPI. Tested using Tested using Ubuntu 20.04 and gfortran 10
# ------------------------------------------------------------------------------
#export PETSC_DIR=$(PARAFEMGEMS_HOME)/lib_paragems/src/libraries/PETSC
#export PETSC_ARCH=arch-linux2-c-debug
#include ${PETSC_DIR}/lib/petsc/conf/variables

# ------------------------------------------------------------------------------
# Generic settings for compilation
# ------------------------------------------------------------------------------
FC=mpif90
LD=
INSTALL_CMD=mv

# ------------------------------------------------------------------------------
# FLAGS for MODULE/LIBRARY compilation
# ------------------------------------------------------------------------------
# These flags are directly used by the module makefiles
# ------------------------------------------------------------------------------
#LIB_REL_FFLAGS=-O2 -cpp $(PETSC_FC_INCLUDES) -fallow-argument-mismatch
#LIB_DBG_FFLAGS=-g -fbounds-check -cpp $(PETSC_FC_INCLUDES) -fallow-argument-mismatch
LIB_REL_FFLAGS=-O2 -cpp -fallow-argument-mismatch
LIB_DBG_FFLAGS=-g -fbounds-check -cpp -fallow-argument-mismatch

# ------------------------------------------------------------------------------
# DEFAULT FLAGS for all program groups
# ------------------------------------------------------------------------------
# These are NOT directly used in any makefiles; specific program flags are
# assigned these defaults below and may be overridden.
# ------------------------------------------------------------------------------
# generic darcy group FFLAGS and LDFLAGS
DEFAULT_GEMS_REL_FFLAGS=-O2 -fallow-argument-mismatch
DEFAULT_GEMS_DBG_FFLAGS=-g -fbounds-check -fallow-argument-mismatch
DEFAULT_GEMS_REL_LDFLAGS=-lblas -llapack $(PETSC_LIB)
DEFAULT_GEMS_DBG_LDFLAGS=-lblas -llapack $(PETSC_LIB)

# ------------------------------------------------------------------------------
# PROGRAM FLAGS for all programs.
# ------------------------------------------------------------------------------
GEMS_REL_FFLAGS=$(DEFAULT_darcy_REL_FFLAGS)
GEMS_REL_LDFLAGS=$(DEFAULT_darcy_REL_LDFLAGS)
GEMS_DBG_FFLAGS=$(DEFAULT_darcy_DBG_FFLAGS)
GEMS_DBG_LDFLAGS=$(DEFAULT_darcy_DBG_LDFLAGS)
