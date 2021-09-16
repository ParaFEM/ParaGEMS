Description
===========

ParaGEMS is a parallel math library for discrete exterior calculus, along with
miniApps for various applications in science and engineering.

Authors: Dr Pieter Boom (pieter.boom@manchester.ac.uk);
Dr Lee Margetts (lee.margetts@manchester.ac.uk)

Updated repository: https://bitbucket.org/pieterboom/paragems/

File/Directory Structure
========================

bin/
----
Directory for compiled executables created at time of compilation

config/
-------
Directory for system specific configuration files (includes template)

docs/
-----
Automatic documentation of ParaGEMS  

- doc_ford:  documentation built with FORD  
- doc_robo:  documentation built with ROBODoc  
- make-docs: script to recompile documentation  

examples/
---------
Directory with sample test cases  

- ICOS:        steady Darcy flow through icosahedron with 20 cells  
- CUBE:        steady Darcy flow through unit cube with 20981 cells  
- CUBE_SCrack: steady Darcy flow through unit cube with stochastic cracking  
- CUBE_DCrack: steady Darcy flow through unit cube with deterministic cracking  

include/
--------
Directory of compiled module files (.mod) created at time of compilation

lib/
----
Directory of compiled library files (.a) created at time of compilation

license/
--------
Directory with software license

src/
----
Directory with source code, including any optional libraries, modules, programs
(mini apps), and code templates.  

- libraries/:             directory for optional libraries, such as BLAS, LAPACK
    & PETSc  
- modules/common:         global variable definitions  
- modules/dec:            geometry and Discrete Exterior Calculus (DEC)
    subroutines  
- modules/io:             I/O subroutines  
- modules/math:           general math subroutines  
- modules/mpi:            MPI subroutines  
- modules/phy_darcy_flow: Darcy flow specific subroutines  
- modules/solvers:        Integration of PETSc solvers (currently only KSP
    solvers)  
- programs/darcy_flow/:   MiniApps for steady Darcy flow (with cracking)  
- templates:              Coding standard - templates  

utils/
------
Software tools created for use with ParaGEMS

make-paragems
-------------
Script to compile ParaGEMS code (see installation instruction below)


Installation
============

Prerequisites:
--------------
- Compilers: ParaGEMS has been compiled and tested using using both Intel (17.0
    and 18.0) and GCC (8, 9 and 10) compilers. It may be possible to compile
    with other compilers, but this is has not be tested. NOTE: on Apple Macs you
    CANNOT use the built-in GCC compiler and MPI library - you MUST install an
    independent set of compilers and libraries from, for example, Homebrew or
    source. Make sure you update your PATH environment variable to point to
    these new installations  
- MPI: ParaGEMS has been compiled and tested with OpenMPI 4.0.1  
- BLAS and LAPACK: An installation of BLAS and LAPACK is required for compiling
    both PETSC and ParaGEMS.  
- PETSc: ParaGEMS has been compiled and tested using PETSc 3.12 (recommended).
    It has also been compiled and tested with PETSc 3.14. PETSC must be compiled
    with the Hypre package (--download-hypre).  
- Python: ParaGEMS has been compiled and tested using PETSc 3.12 and python 2.7.
    It has also been compiled and tested with PETSc 3.14 and python 3.  

Installation of ParaGEMS:
-------------------------
- Enter the ParaGEMS directory, which is referred to as $(PARAGEMS_HOME)  
- Explore the machine configuration files in config/ directory. These can be
    adapted for your own system. Pay particular attention to the PATH described
    for your PETSc installation directory. Save the configuration as
    [machine].inc  
- Compilation is executed for your [machine] with the command:  
    \>> ./make_paragemms MACHINE=[machine]
- If successful, the compiled ParaGEMS executables can then be found in the
    $(PARAGEMS_HOME)/bin directory  
- If using MacOS or Linux you can add the $(PARAGEMS_HOME)/bin directory to your
    PATH environment variable:  
    \>> export PATH=$(PARAGEMS_HOME)/bin:$PATH  

Installation of utilities:
--------------------------
- Enter utils/ directory
- Compile utilities with command:  
    \>> make

Execution
=========

Test case: steady Darcy flow through icosahedron in parallel with 2 processes  

Directory: $(PARAGEMS_HOME)/examples/test_darcy/ICOS  

Files:  

- input.param    (input): input file  
- icos.1.node    (mesh): node location of mesh  
- icos.1.edge    (mesh): line definition of mesh (node indices)  
- icos.1.face    (mesh): face definition of mesh (node indices) including
    boundary conditions  
- icos.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- paragems.log\* (output): solver log file  
- solution.m     (output): solution vector in MATLAB format (order associated
    with parallel partitioning)  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order
    associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters
    in MATLAB format (order associated with parallel partitioning)  
- solution.vtk   (output): solution output file (can be opened with ParaView)  

Execution:  

\>> cd $(PARAGEMS_HOME)/examples/test_darcy/ICOS  
\>> mpirun -np 2 $(PARAGEMS_HOME)/bin/darcy_2f input.param  

> **or** if bin is in PATH environment variable:  
> \>> mpirun -np 2 darcy_2f input.param

----

Test case: steady Darcy flow through unit cube in parallel with 4 processes
(optional use of simple load balancing tool)  

Directory: $(PARAGEMS_HOME)/examples/test_darcy/CUBE  

Files:  

- input.param    (input): input file  
- cube.poly      (mesh): mesh definition for use with TetGen software  
- cube.1.node    (mesh): node location of mesh  
- cube.1.edge    (mesh): line definition of mesh (node indices)  
- cube.1.face    (mesh): face definition of mesh (node indices) with boundary
    conditions  
- cube.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- cube.1.\*.old  (mesh): mesh files before load balancing  
- paragems.log\* (output): solver log file  
- solution.m     (output): solution vector in MATLAB format (order associated
    with parallel partitioning)  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order
  associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters
    in MATLAB format (order associated with parallel partitioning)  
- solution.vtk   (output): solution output file (can be opened with ParaView)  

Execution:  

\>> cd $(PARAGEMS_HOME)/examples/test_darcy/CUBE  
\>> mpirun -np 4 darcy_2f input.param  
\>> simpleloadbalance  
    cube.1  
\>> mpirun -np 4 darcy_2f input.param  

----

Test case: steady Darcy flow through unit cube with deterministic cracking in
parallel with 4 processes (10 cracks)

Directory: $(PARAGEMS_HOME)/examples/test_darcy/CUBE_DCrack

Files:  

-  input.param   (input): input file  
- cube.poly      (mesh): mesh definition for use with TetGen software  
- cube.1.node    (mesh): node location of mesh  
- cube.1.edge    (mesh): line definition of mesh (node indices)  
- cube.1.face    (mesh): face definition of mesh (node indices) with boundary
    conditions  
- cube.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- cube.1.\*.old  (mesh): mesh files before load balancing  
- paragems.log\* (output): solver log file  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order
    associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters
    in MATLAB format (order associated with parallel partitioning)  
- unsteady\*.vtk (output): solution output file (can be opened with ParaView)  
- unsteady\*.m   (output): solution vector in MATLAB format (order associated
    with parallel partitioning)  
- unsteady\*press.m (output): interpolated pressure solution vector in MATLAB
    format (order associated with parallel partitioning) - face values scaled by
    area  
- unsteady.log   (output): log of cracking process  

Execution:  

\>> cd $(PARAGEMS_HOME)/examples/test_darcy/CUBE_DCrack  
\>> simpleloadbalance  
    cube.1  
\>> mpirun -np 4 darcy_crkp_2f input.param  

----

Test case: steady Darcy flow through unit cube with stochastic cracking in
parallel with 4 processes (100 random cracks in 10 steps)

Directory: $(PARAGEMS_HOME)/examples/test_darcy/CUBE_SCrack

Files:  

-  input.param   (input): input file  
- cube.poly      (mesh): mesh definition for use with TetGen software  
- cube.1.node    (mesh): node location of mesh  
- cube.1.edge    (mesh): line definition of mesh (node indices)  
- cube.1.face    (mesh): face definition of mesh (node indices) with boundary
    conditions  
- cube.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- cube.1.\*.old  (mesh): mesh files before load balancing  
- paragems.log\* (output): solver log file  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order
    associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters
    in MATLAB format (order associated with parallel partitioning)  
- unsteady\*.vtk (output): solution output file (can be opened with ParaView)  
- unsteady\*.m   (output): solution vector in MATLAB format (order associated
    with parallel partitioning)  
- unsteady\*press.m (output): interpolated pressure solution vector in MATLAB
    format (order associated with parallel partitioning) - face values scaled by
    area  
- unsteady.log   (output): log of cracking process  

Execution:  

\>> cd $(PARAGEMS_HOME)/examples/test_darcy/CUBE_SCrack  
\>> simpleloadbalance  
    cube.1  
\>> mpirun -np 4 darcy_crkp_2f input.param
