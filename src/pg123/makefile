# MAKEFILE src/pg123/makefile
# Author: Pieter Boom (adapted from p123 by Louise Lever)
# -----------------------------------------------------------------------
# Compiles pg123 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build pg123 executable as RELEASE edition
#   debug: Build pg123 executable as DEBUG edition
#   install: Copy pg123 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pg123 executable
#   relink: Rebuild pg123 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../build/parafem/mk_defs.inc
include ../../build/parafem/$(MACHINE).inc

MODULE	= ../../include/parafem/mpi
LIB	= ../../lib
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building PG123 RELEASE"
	$(MAKE) pg123 \
	FFLAGS="-c -I$(MODULE) $(PG123_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(PG123_MPI_REL_LDFLAGS)"
	@echo "Done pg123 RELEASE"

debug:
	@echo
	@echo "Building pg123 DEBUG"
	$(MAKE) pg123 \
	FFLAGS="-c -I$(MODULE) $(PG123_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(PG123_MPI_DBG_LDFLAGS)"
	@echo "Done pg123 DEBUG"

pg123: pg123.o
	$(FC) pg123.o -o pg123 $(LDFLAGS)

pg123.o: pg123.f90
	$(FC) pg123.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pg123 ../../bin/pg123

relink: pg123

install:
	cp pg123 ../../bin
