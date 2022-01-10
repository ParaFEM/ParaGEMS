# MAKEFILE src/pg121/makefile
# Author: Pieter Boom (adapted from pg121 by Louise Lever)
# -----------------------------------------------------------------------
# Compiles pg121 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build pg121 executable as RELEASE edition
#   debug: Build pg121 executable as DEBUG edition
#   install: Copy pg121 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pg121 executable
#   relink: Rebuild pg121 from object files and libraries
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
	@echo "Building PG121 RELEASE"
	$(MAKE) pg121 \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG121_MPI_REL_FFLAGS)" \
	LDFLAGS="-lblas -llapack -L$(LIB_FEM) -lParaFEM_mpi.$(VERSION) -L$(LIB_GEMS) -lParaGEMS.2021.03.16 $(PG121_MPI_REL_LDFLAGS)"
	@echo "Done pg121 RELEASE"

debug:
	@echo
	@echo "Building pg121 DEBUG"
	$(MAKE) pg121 \
	FFLAGS="-c -I$(MODULE_FEM) -I$(MODULE_GEMS) $(PG121_MPI_DBG_FFLAGS)" \
	LDFLAGS="-lblas -llapack -L$(LIB_FEM) -lParaFEM_mpi_D.$(VERSION) -L$(LIB_GEMS) -lParaGEMS_D.2021.03.16  $(PG121_MPI_DBG_LDFLAGS)"
	@echo "Done pg121 DEBUG"

pg121: pg121.o
	$(FC) pg121.o -o pg121 $(LDFLAGS)

pg121.o: pg121.f90
	$(FC) pg121.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pg121 ../../bin/pg121

relink: pg121

install:
	cp pg121 ../../bin

test:
	./test.sh pg121
