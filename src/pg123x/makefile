# MAKEFILE src/pg123x/makefile
# Author: Pieter Boom (adapted from p123 by Louise Lever)
# -----------------------------------------------------------------------
# Compiles pg123x program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build pg123x executable as RELEASE edition
#   debug: Build pg123x executable as DEBUG edition
#   install: Copy pg123x executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pg123x executable
#   relink: Rebuild pg123x from object files and libraries
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
	@echo "Building PG123x RELEASE"
	$(MAKE) pg123x \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG123_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB_FEM) -lParaFEM_mpi.$(VERSION) -L$(LIB_GEMS) -lParaGEMS.2022.02.28 $(PG123_MPI_REL_LDFLAGS) -lblas -llapack"
	@echo "Done pg123x RELEASE"

debug:
	@echo
	@echo "Building PG123x DEBUG"
	$(MAKE) pg123x \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG123_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB_FEM) -lParaFEM_mpi_D.$(VERSION) -L$(LIB_GEMS) -lParaGEMS_D.2022.02.28  $(PG123_MPI_DBG_LDFLAGS) -lblas -llapack"
	@echo "Done pg123x DEBUG"

pg123x: pg123x.o
	$(FC) pg123x.o -o pg123x $(LDFLAGS)

pg123x.o: pg123x.f90
	$(FC) pg123x.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pg123x ../../bin/pg123x

relink: pg123x

install:
	cp pg123x ../../bin
