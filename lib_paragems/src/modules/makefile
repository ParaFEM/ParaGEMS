# MAKEFILE src/modules/makefile
# Author: Pieter D. Boom
# Adapted from parafem (Louise M. Lever)
# -----------------------------------------------------------------------
# Options:
# (default) : Compile and build library
# install: Copy modules into global include directory
#          Copy library into global lib directory
# -----------------------------------------------------------------------

include $(PARAGEMS_HOME)/config/mk_defs.inc
include $(PARAGEMS_HOME)/config/$(MACHINE).inc

FOBJS =	\
	common/common_mod.o \
	mpi/mpi_mod.o \
	time_marching/time_marching_mod.o \
	math/math_mod.o \
	dec/dec_mod.o \
	dec/geometry_mod.o \
	dec/forman_mod.o \
	solvers/solver_mod.o	\
	io/io_mod.o \
	mpi/partition_mod.o \
	phy_diffusion/diffusion_mod.o	\
	testing/testing_mod.o	\
	common/test_common_mod.o \
	dec/test_dec_mod.o \
	io/test_io_mod.o \
	math/test_math_mod.o \
	mpi/test_mpi_mod.o \
	phy_diffusion/test_diffusion_mod.o \
	solvers/test_solver_mod.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE"
	$(MAKE) FFLAGS="$(LIB_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE"

debug:
	@echo
	@echo "Building MODULE DEBUG"
	$(MAKE) FFLAGS="$(LIB_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG"

lib-release:	$(FOBJS)
	ar -r libParaGEMS.$(VERSION).a $(FOBJS)

lib-debug:	$(FOBJS)
	ar -r libParaGEMS_D.$(VERSION).a $(FOBJS)

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

clean:
	rm -f */*.o *.mod *.a

execlean:
	rm -f $(PARAGEMS_INC)/*.mod
	rm -f $(PARAGEMS_LIB)/libParaGEMS*.a

install:
	cp *.mod $(PARAGEMS_INC)
	cp *.a $(PARAGEMS_LIB)
