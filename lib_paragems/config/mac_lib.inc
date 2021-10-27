# BUILD INCLUDE config/mac.inc
# Authors: P.D. Boom
# ------------------------------------------------------------------------------
# Include file for standard MacOS workstation compilation using gfortran
# and OpenMPI. Tested using MacOS 10.14 and gfortran 10
# ------------------------------------------------------------------------------
ifeq ($(FOR_paraFEM),0)
  export PETSC_DIR=$(PARAGEMS_HOME)/src/libraries/PETSC
  export PETSC_ARCH=arch-darwin-c-debug
  include ${PETSC_DIR}/lib/petsc/conf/variables
endif

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
ifeq ($(FOR_paraFEM),0)
  LIB_REL_FFLAGS=-O2 -cpp -fallow-argument-mismatch $(PETSC_FC_INCLUDES)
  LIB_DBG_FFLAGS=-g -fbounds-check -cpp -fallow-argument-mismatch $(PETSC_FC_INCLUDES)
else
  LIB_REL_FFLAGS=-O2 -cpp -fallow-argument-mismatch -DFORparaFEM
  LIB_DBG_FFLAGS=-g -fbounds-check -cpp -fallow-argument-mismatch -DFORparaFEM
endif

# ------------------------------------------------------------------------------
# DEFAULT FLAGS for all program groups
# ------------------------------------------------------------------------------
# These are NOT directly used in any makefiles; specific program flags are
# assigned these defaults below and may be overridden.
# ------------------------------------------------------------------------------
# generic darcy group FFLAGS and LDFLAGS
DEFAULT_darcy_REL_FFLAGS=-O2 -fallow-argument-mismatch
DEFAULT_darcy_DBG_FFLAGS=-g -fbounds-check -fallow-argument-mismatch
DEFAULT_darcy_REL_LDFLAGS=-lblas -llapack $(PETSC_LIB)
DEFAULT_darcy_DBG_LDFLAGS=-lblas -llapack $(PETSC_LIB)

# ------------------------------------------------------------------------------
# PROGRAM FLAGS for all programs.
# ------------------------------------------------------------------------------
darcy_REL_FFLAGS=$(DEFAULT_darcy_REL_FFLAGS)
darcy_REL_LDFLAGS=$(DEFAULT_darcy_REL_LDFLAGS)
darcy_DBG_FFLAGS=$(DEFAULT_darcy_DBG_FFLAGS)
darcy_DBG_LDFLAGS=$(DEFAULT_darcy_DBG_LDFLAGS)

# ------------------------------------------------------------------------------
# TESTING FLAGS for all tests.
# ------------------------------------------------------------------------------
paragems_unit_REL_FFLAGS=$(DEFAULT_darcy_REL_FFLAGS)
paragems_unit_REL_LDFLAGS=$(DEFAULT_darcy_REL_LDFLAGS)
paragems_unit_DBG_FFLAGS=$(DEFAULT_darcy_DBG_FFLAGS)
paragems_unit_DBG_LDFLAGS=$(DEFAULT_darcy_DBG_LDFLAGS)
