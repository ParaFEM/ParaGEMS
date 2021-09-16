# BUILD INCLUDE config/csf.inc
# Authors: P.D. Boom
# ------------------------------------------------------------------------------
# Include file for Computational Shared Facility at the University of Manchester
# ------------------------------------------------------------------------------
export PETSC_DIR=$(PARAGEMS_HOME)/src/libraries/PETSC
export PETSC_ARCH=arch-linux2-c-debug
include ${PETSC_DIR}/lib/petsc/conf/variables

# ------------------------------------------------------------------------------
# Generic settings for compilation
# ------------------------------------------------------------------------------
FC=mpifort
LD=
INSTALL_CMD=mv

# ------------------------------------------------------------------------------
# FLAGS for MODULE/LIBRARY compilation
# ------------------------------------------------------------------------------
# These flags are directly used by the module makefiles
# ------------------------------------------------------------------------------
LIB_REL_FFLAGS=-O2 -cpp $(PETSC_FC_INCLUDES)
LIB_DBG_FFLAGS=-g -fbounds-check -cpp $(PETSC_FC_INCLUDES)

# ------------------------------------------------------------------------------
# DEFAULT FLAGS for all program groups
# ------------------------------------------------------------------------------
# These are NOT directly used in any makefiles; specific program flags are
# assigned these defaults below and may be overridden.
# ------------------------------------------------------------------------------
# generic darcy group FFLAGS and LDFLAGS
DEFAULT_darcy_REL_FFLAGS=-O2 -cpp $(PETSC_FC_INCLUDES)
DEFAULT_darcy_DBG_FFLAGS=-g -fbounds-check -cpp $(PETSC_FC_INCLUDES)
DEFAULT_darcy_REL_LDFLAGS=$(PETSC_LIB)
DEFAULT_darcy_DBG_LDFLAGS=$(PETSC_LIB)

# ------------------------------------------------------------------------------
# DEFAULT FLAGS for all test groups
# ------------------------------------------------------------------------------
# These are NOT directly used in any makefiles; specific test flags are
# assigned these defaults below and may be overridden.
# ------------------------------------------------------------------------------

# generic paragems_unit group FFLAGS and LDFLAGS
DEFAULT_paragems_unit_REL_FFLAGS=-O2 -cpp $(PETSC_FC_INCLUDES)
DEFAULT_paragems_unit_DBG_FFLAGS=-g -fbounds-check -cpp $(PETSC_FC_INCLUDES)
DEFAULT_paragems_unit_REL_LDFLAGS=$(PETSC_LIB)
DEFAULT_paragems_unit_DBG_LDFLAGS=$(PETSC_LIB)

# ------------------------------------------------------------------------------
# PROGRAM FLAGS for all programs.
# ------------------------------------------------------------------------------
# All use their program group defaults defined above unless overridden below
# ------------------------------------------------------------------------------
darcy_REL_FFLAGS=$(DEFAULT_darcy_REL_FFLAGS)
darcy_REL_LDFLAGS=$(DEFAULT_darcy_REL_LDFLAGS)
darcy_DBG_FFLAGS=$(DEFAULT_darcy_DBG_FFLAGS)
darcy_DBG_LDFLAGS=$(DEFAULT_darcy_DBG_LDFLAGS)

# ------------------------------------------------------------------------------
# TEST FLAGS for all tests.
# ------------------------------------------------------------------------------
# All use their test group defaults defined above unless overridden below
# ------------------------------------------------------------------------------
paragems_unit_REL_FFLAGS=$(DEFAULT_paragems_unit_REL_FFLAGS)
paragems_unit_REL_LDFLAGS=$(DEFAULT_paragems_unit_REL_LDFLAGS)
paragems_unit_DBG_FFLAGS=$(DEFAULT_paragems_unit_DBG_FFLAGS)
paragems_unit_DBG_LDFLAGS=$(DEFAULT_paragems_unit_DBG_LDFLAGS)
