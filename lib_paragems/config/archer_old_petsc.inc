# BUILD INCLUDE build/mac.inc
# Authors: P.D. Boom
# ------------------------------------------------------------------------------
# Include file for standard MacOS workstation compilation using gfortran
# and OpenMPI. Tested using MacOS 10.14 and gfortran 6.3
# ------------------------------------------------------------------------------

export PARAGEMS_LIB=/work/ecseaa12/ecseaa12/pboom/codes/paragems_old_petsc/lib
export PARAGEMS_INC=/work/ecseaa12/ecseaa12/pboom/codes/paragems_old_petsc/include
export PARAGEMS_BIN=/work/ecseaa12/ecseaa12/pboom/codes/paragems_old_petsc/bin

export PETSC_DIR=/work/ecseaa12/ecseaa12/pboom/
export PETSC_ARCH=""
include ${PETSC_DIR}/lib/petsc/conf/variables

# ------------------------------------------------------------------------------
# Generic settings for compilation
# ------------------------------------------------------------------------------
FC=ftn
LD=
INSTALL_CMD=mv

# ------------------------------------------------------------------------------
# FLAGS for MODULE/LIBRARY compilation
# ------------------------------------------------------------------------------
# Form is LIB_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
# These flags are directly used by the module makefiles
# ------------------------------------------------------------------------------
LIB_REL_FFLAGS=-O2 -e Z $(PETSC_FC_INCLUDES)
LIB_DBG_FFLAGS=-g -fbounds-check -e Z $(PETSC_FC_INCLUDES)

# ------------------------------------------------------------------------------
# DEFAULT FLAGS for all program groups
# ------------------------------------------------------------------------------
# Form is DEFAULT_(darcy)_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
#   where applicable.
# These are NOT directly used in any makefiles; specific program flags are
# assigned these defaults below and may be overridden.
# ------------------------------------------------------------------------------
# generic darcy group FFLAGS and LDFLAGS
DEFAULT_darcy_REL_FFLAGS=-O2 -e Z $(PETSC_FC_INCLUDES)
DEFAULT_darcy_DBG_FFLAGS=-g -fbounds-check -e Z $(PETSC_FC_INCLUDES)
DEFAULT_darcy_REL_LDFLAGS=$(PETSC_LIB)
DEFAULT_darcy_DBG_LDFLAGS=$(PETSC_LIB)

# ------------------------------------------------------------------------------
# DEFAULT FLAGS for all test groups
# ------------------------------------------------------------------------------
# Form is DEFAULT_(paragems_unit)_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
#   where applicable.
# These are NOT directly used in any makefiles; specific test flags are
# assigned these defaults below and may be overridden.
# ------------------------------------------------------------------------------

# generic paragems_unit group FFLAGS and LDFLAGS
DEFAULT_paragems_unit_REL_FFLAGS=-O2 -e Z $(PETSC_FC_INCLUDES)
DEFAULT_paragems_unit_DBG_FFLAGS=-g -fbounds-check -e Z $(PETSC_FC_INCLUDES)
DEFAULT_paragems_unit_REL_LDFLAGS=$(PETSC_LIB)
DEFAULT_paragems_unit_DBG_LDFLAGS=$(PETSC_LIB)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SPECIFIC FLAGS for PROGRAMS and TOOLS
# By default these inherit the defaults defined above.
# NOTE: If a specific program needs different flags then modify below.
# NOTE: It may be useful to create a new include file for a modified
#       architecture.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

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

