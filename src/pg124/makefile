# MAKEFILE src/pg124/makefile
# Author: Pieter Boom (adapted from p124 by Louise Lever)
# -----------------------------------------------------------------------
# Compiles pg124 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build pg124 executable as RELEASE edition
#   debug: Build pg124 executable as DEBUG edition
#   install: Copy pg124 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pg124 executable
#   relink: Rebuild pg124 from object files and libraries
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
	@echo "Building PG124 RELEASE"
	$(MAKE) pg124 \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG124_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB_FEM) -lParaFEM_mpi.$(VERSION) -L$(LIB_GEMS) -lParaGEMS.2022.02.28 $(PG124_MPI_REL_LDFLAGS) -lblas -llapack"
	@echo "Done pg124 RELEASE"

debug:
	@echo
	@echo "Building pg124 DEBUG"
	$(MAKE) pg124 \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG124_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB_FEM) -lParaFEM_mpi_D.$(VERSION) -L$(LIB_GEMS) -lParaGEMS_D.2022.02.28  $(PG124_MPI_DBG_LDFLAGS) -lblas -llapack"
	@echo "Done pg124 DEBUG"

pg124: pg124.o
	$(FC) pg124.o -o pg124 $(LDFLAGS)

pg124.o: pg124.f90
	$(FC) pg124.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pg124 ../../bin/pg124

relink: pg124

install:
	cp pg124 ../../bin
