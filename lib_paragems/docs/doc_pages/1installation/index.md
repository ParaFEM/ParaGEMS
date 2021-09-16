title: Installation
author: Pieter Boom
date: 2021/01/26

Prerequisites:

- Compilers: ParaGEMS has been compiled and tested using using both Intel (17.0 and 18.0) and GCC (8, 9 and 10) compilers. It may be possible to compile with other compilers, but this is has not be tested.

  NOTE: on Apple Macs you CANNOT use the built-in GCC compiler - you MUST install a complete set of GCC compilers from, for example, Homebrew. Make sure you update your PATH environment variable to point to this new installation

- MPI: ParaGEMS has been compiled and tested with OpenMPI 4.0.1

- Python: ParaGEMS has been compiled and tested using PETSc 3.12 and python 2.7. It has also been compiled and tested with PETSc 3.14 and python 3.

- BLAS and LAPACK: An installation of BLAS and LAPACK is required for compiling both PETSC and ParaGEMS.

- PETSc: ParaGEMS has been compiled and tested using PETSc 3.12 (recommended). It has also been compiled and tested with PETSc 3.14. PETSC must be compiled with the Hypre package (--download-hypre).

Installation:

- Enter the ParaGEMS directory, which is referred to as $(PARAGEMS_HOME)

- Explore the machine configuration files in config/ directory. These can be adapted for your own system. Pay particular attention to the PATH described for your PETSc installation directory. Save your configuration as [machine].inc

- Compilation is executed for your [machine] with the command:
  ./make_paragemms MACHINE=[machine]

- If successful, the compiled ParaGEMS executables can then be found in the $(PARAGEMS_HOME)/bin directory
