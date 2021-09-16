# Makefile for Utilities
#
# Pieter Boom

FOBJS = simpleLoadBalance simpleLoadBalance_2
FC = gfortran

all: $(FOBJS) install

$(FOBJS): % : %.f90
	$(FC) $< -o $@

clean::
	rm -rf $(FOBJS) *~

install:
	cp $(FOBJS) ../bin
