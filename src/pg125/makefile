# MAKEFILE src/programs/5th_ed/pg125/makefile
# Author: Pieter Boom (adapted from p125 by Louise Lever)
# -----------------------------------------------------------------------
# Compiles pg125 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build pg125 executable as RELEASE edition
#   debug: Build pg125 executable as DEBUG edition
#   install: Copy pg125 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pg125 executable
#   relink: Rebuild pg125 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../build/parafem/mk_defs.inc
include ../../build/$(MACHINE).inc

MODULE_FEM	= ../../include/parafem/mpi
LIB_FEM	= ../../lib/parafem
MODULE_GEMS	= ../../include/paragems
LIB_GEMS	= ../../lib/paragems
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building PG125 RELEASE"
	$(MAKE) pg125 \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG125_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB_FEM) -lParaFEM_mpi.$(VERSION) -L$(LIB_GEMS) -lParaGEMS.2022.02.28 $(PG125_MPI_REL_LDFLAGS) -lblas -llapack"
	@echo "Done pg125 RELEASE"

debug:
	@echo
	@echo "Building pg125 DEBUG"
	$(MAKE) pg125 \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG125_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB_FEM) -lParaFEM_mpi_D.$(VERSION) -L$(LIB_GEMS) -lParaGEMS_D.2022.02.28  $(PG125_MPI_DBG_LDFLAGS) -lblas -llapack"
	@echo "Done pg125 DEBUG"

pg125: pg125.o
	$(FC) pg125.o -o pg125 $(LDFLAGS)

pg125.o: pg125.f90
	$(FC) pg125.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pg125 ../../bin/pg125

relink: pg125

install:
	cp pg125 ../../bin
