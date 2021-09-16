# MAKEFILE src/tests/unit/makefile
# Author: Pieter D Boom
# Adapted from parafem (Louise M. Lever)
# -----------------------------------------------------------------------
# Compiles paragems_unit program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaGEMS.$(VERSION).a
#    or with RELEASE STATIC library libParaGEMS_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build paragems_unit executable as RELEASE edition
#   debug: Build paragems_unit executable as DEBUG edition
#   install: Copy paragems_unit executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local paragems_unit executable
#   relink: Rebuild paragems_unit from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include $(PARAGEMS_HOME)/config/mk_defs.inc
include $(PARAGEMS_HOME)/config/$(MACHINE).inc

MODULE	= $(PARAGEMS_INC)
LIB	= $(PARAGEMS_LIB)

all:	release

release:
	@echo
	@echo "Building paragems_unit RELEASE"
	$(MAKE) paragems_unit \
	FFLAGS="-c -I$(MODULE) $(paragems_unit_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaGEMS.$(VERSION) $(paragems_unit_REL_LDFLAGS)"
	@echo "Done paragems_unit RELEASE"

debug:
	@echo
	@echo "Building paragems_unit DEBUG"
	$(MAKE) paragems_unit \
	FFLAGS="-c -I$(MODULE) $(paragems_unit_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaGEMS_D.$(VERSION) $(paragems_unit_DBG_LDFLAGS)"
	@echo "Done paragems_unit DEBUG"

paragems_unit: paragems_unit.o
	$(FC) paragems_unit.o -o paragems_unit $(LDFLAGS)

paragems_unit.o: FORCE
	$(FC) paragems_unit.f90 $(FFLAGS)

FORCE:

clean:
	rm -f *.o

execlean:
	rm -f paragems_unit $(PARAGEMS_BIN)/paragems_unit

relink: paragems_unit

install:
	cp paragems_unit $(PARAGEMS_HOME)/bin
