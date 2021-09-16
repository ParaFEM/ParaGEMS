# MAKEFILE src/programs/darcy_flow/makefile
# Author: Pieter D Boom
# Adapted from parafem (Louise M. Lever)
# -----------------------------------------------------------------------
# Compiles darcy_flow program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaGEMS.$(VERSION).a
#    or with RELEASE STATIC library libParaGEMS_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build darcy executable as RELEASE edition
#   debug: Build darcy executable as DEBUG edition
#   install: Copy darcy executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local darcy executable
#   relink: Rebuild darcy from object files and libraries
# -----------------------------------------------------------------------

include $(PARAGEMS_HOME)/config/mk_defs.inc
include $(PARAGEMS_HOME)/config/$(MACHINE).inc

MODULE	= $(PARAGEMS_INC)
LIB	= $(PARAGEMS_LIB)

TARGETS = ohm_2f\
	darcy_2f #\
	darcy_crkp_2f \
	darcy_2f_time \
	ohm_2f \
	ohm_2f_time

all:	release

release:
	@echo
	@echo "Building darcy RELEASE"
	$(MAKE) darcy_prgms \
	FFLAGS="-c -I$(MODULE) $(darcy_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaGEMS.$(VERSION) $(darcy_REL_LDFLAGS)"
	@echo "Done darcy RELEASE"

debug:
	@echo
	@echo "Building darcy DEBUG"
	$(MAKE) darcy_prgms \
	FFLAGS="-c -I$(MODULE) $(darcy_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaGEMS_D.$(VERSION) $(darcy_DBG_LDFLAGS)"
	@echo "Done darcy DEBUG"

darcy_prgms: $(TARGETS)

%: %.o
	$(FC) $< -o $@ $(LDFLAGS)

%.o: %.f90 FORCE
	$(FC) $< $(FFLAGS)

FORCE:

clean:
	rm -f *.o

execlean:
	rm -f $(TARGETS) $(PARAGEMS_BIN)/$(TARGETS)

relink: darcy_prgms

install:
	cp $(TARGETS) $(PARAGEMS_BIN)
