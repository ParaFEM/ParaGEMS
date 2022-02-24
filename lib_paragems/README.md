Description
===========

ParaGEMS is a parallel math library for discrete exterior calculus, along with miniApps for various applications in science and engineering. The geometric and topological subroutines/functions can also be compiled as a separate library to be integrated with the ParaFEM finite-element library.

Authors: Dr Pieter Boom (pieter.boom@manchester.ac.uk); Dr Lee Margetts (lee.margetts@manchester.ac.uk)

Updated repository: https://bitbucket.org/pieterboom/paragems/; https://github.com/ParaFEM/ParaGEMS

Installation
============

The following section describes the installation of the ParaGEMS library and associated prerequisites on a Linux machine, with notes for MacOS and Windows Subsystem for Linux (WSL).

Prerequisites
-------------

The software has been installed and tested on both Linux (Xubuntu 20.04 LTS, CentOS 7.9.2009) and MacOS (High Sierra, Catalina, and Big Sur). Other operating systems are NOT supported.

ParaGEMS has been compiled and tested using using the following compilers and libraries (Other compilers and libraries are NOT supported):

 - Intel ifort (17.0 and 18.0), GCC gfortran (8, 9 and 10), Cray Compiling Environment (CCE), and AMD compiler environment (AOCC)
 - OpenMPI 4.0.1, MPICH (8.1.4)
 - PETSc (3.12, 3.14) compiled with MPI and including the Hypre package (not needed for ParaFEM integration)

Other software needed:

 - Git (see section: Getting the software)
 - BLAS and LAPACK, which can be installed on Ubunutu using the command:

    sudo apt-get install libblas-dev liblapack-dev


Optional software:

 - Triangle for 2D mesh generation (https://www.cs.cmu.edu/~quake/triangle.html)
 - TetGen for 3D mesh generation (https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1)

Compilers and MPI on MacOS
--------------------------

On MacOS you CANNOT use the pre-installed GCC compilers and MPI library; you MUST install an independent set of compilers and MPI library.

> Example installing GCC 9 from Homebrew (see https://formulae.brew.sh/formula/gcc):

    brew install gcc@9

> OpenMPI can likewise be installed from Homebrew. For installing OpenMPI from source see section 7 from: https://www.open-mpi.org/faq/?category=osx

Make sure you update your $PATH environment variable to point to the newly installed compilers and libraries

> Example for GCC 9 installed via Homebrew:

    export PATH="/opt/local/Cellar/gcc@9/9.3.0_1/bin:$PATH"

This command can be appended to the end of your ~/.zshrc file to execute it every time a new terminal window is opened.

Windows Subsystem for Linux (WSL) Ubuntu 20.04
----------------------------------------------

Known issues:

- One user reported having to first install valgrind and build-essential (c++) before compiling PETSc
- OpenMPI has difficulty running on >1 core on WSL due to a kernel setting. The fix is to type the following into the command line (https://github.com/Microsoft/WSL/issues/3397):


    echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope

>This command can also be added to your .bashrc or .profile

Getting the software
--------------------

Make sure you have Git (https://git-scm.com) installed on your machine. This can often be done through an appstore like Snap or via the command line, for example on Ubuntu:

    sudo apt-get install git

or on MacOS through Homebrew:

    brew install git

For ParaGEMS standalone (no ParaFEM integration): the ParaGEMS software repository can be cloned to your machine via ssh (requires a Bitbucket account and an SSH key):

    git clone git@bitbucket.org:pieterboom/paragems.git

or via https:

    git clone https://pieterboom@bitbucket.org/pieterboom/paragems.git

For ParaFEM integration: the ParaGEMS software repository should be cloned into the root directory of an existing ParaFEM installation via ssh (requires a git hub account and an SSH key):

    cd [Root ParaFEM directory]
    ls
        README.md          demo-installer-win fruit              
        fruit_3.4.3        parafem            parafem-viewer
    git clone git@github.com:ParaFEM/ParaGEMS.git

or via https:

    cd [Root ParaFEM directory]
    ls
        README.md          demo-installer-win fruit              
        fruit_3.4.3        parafem            parafem-viewer
    git clone https://github.com/ParaFEM/ParaGEMS.git


PETSc
-----

ParaGEMS standalone (no ParaFEM integration) requires PETSc, specifically with MPI and the Hypre package installed. For ParaFEM integration this is not needed.

>If you already have such an installation, you will need to update the path to the PETSc library files in the linuxdesktop.inc file discussed below. You can also skip the rest of this section

To install PETSc, first navigate to the src/libraries/PETSC subdirectory and uncompress the tarball included:

    cd ~/paragems/src/libraries/PETSC
    tar -xvf petsc-3.14.0.tar.gz

>This assumes you cloned the ParaGEMS library to your home directory. If not you will need to update the path in the 'cd' command

Now is a good time to ensure that the configuration script will use the correct compilers and MPI library (especially for Mac users). This can be tested by examining the output of commands like:

    which gcc
    which mpif90

If the output does not match the path to the desired compilers and MPI library, it is recommended that you update your PATH environment variable

>example for compilers assuming the path to their location is: /opt/local/Cellar/gcc@9/9.3.0_1/bin

    export PATH=/opt/local/Cellar/gcc@9/9.3.0_1/bin:$PATH

Once you are satisfied that the output of the 'which' commands are correct, navigate to the uncompressed petsc directory and run the configure script (this may take a while):

    cd petsc-3.14.0
    ./configure --download-hypre

Once the configuration script is complete, you can compile PETSc using the final lines of the configuration output or the configure.log file (this may take a while). The necessary command will look something like:

>make PETSC_DIR=~/paragems/src/libraries/PETSC/petsc-3.14.0 PETSC_ARCH=arch-darwin-c-debug all

Once PETSc is compiled, you can run the provided checks to verify the installation has worked properly using MPI and Hypre. The require command can be found on the last few lines of the make output or make.log file. The necessary command will look something like:

>make PETSC_DIR=~/paragems/src/libraries/PETSC/petsc-3.14.0 PETSC_ARCH=arch-darwin-c-debug check

Compiling ParaGEMS standalone (no ParaFEM integration)
------------------------------------------------------

Navigate to the root ParaGEMS directory:

    cd ~/paragems

>This assumes the root ParaGEMS directory is in your home directory. If not you will need to update the path displayed above

In the config/ subdirectory there are include files for a few different systems, including for a standard linux desktop: linuxdesktop.inc.

>Existing include files (.inc) can be adapted for the particulars of your own system, for example if you are using a previously installed version of PETSc.

ParaGEMS can now be compiled with the command:

    ./make-paragems MACHINE=linuxdesktop

If successful, the compiled ParaGEMS executables can then be found in the bin/ subdirectory of the ParaGEMS directory.

Compiling ParaGEMS for ParaFEM integration
------------------------------------------------------

Navigate to the root ParaGEMS directory:

    cd ~/parafem/paragems/lib_paragems

>This assumes the root ParaFEM directory is in your home directory. If not you will need to update the path displayed above

In the config/ subdirectory there are include files for a few different systems, including for a standard linux desktop: linuxdesktop_lib.inc.

>Existing include files (.inc) can be adapted for the particulars of your own system.

ParaGEMS can now be compiled with the command:

    ./make-paragems MACHINE=linuxdesktop_lib --for-parafem

If successful, the compiled ParaGEMS library and modules can then be found in the include/ and lib/ subdirectories of the ParaGEMS directory.

Installation of utilities:
--------------------------
Enter utils/ subdirectory and compile utilities with make:  

    cd ~/paragems/utils
    make

>This assumes you cloned the ParaGEMS library to your home directory. If not you will need to update the path in the 'cd' command

Setting the PATH environment variable for ParaGEMS (no ParaFEM integration)
---------------------------------------------------------------------------

To make the ParaGEMS and utility executables available while in any directory and every time you open a new terminal window, you can add the paragems/bin directory to your PATH environment variable by adding:

    export PATH=~/paragems/bin:$PATH

>This assumes the root ParaGEMS directory is in your home directory. If not you will need to update the path displayed above

to the end of your ~/.bashrc file (for bash shell) and sourcing it:

    source ~/.bashrc

ParaGEMS File/Directory Structure
=================================

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
Directory with source code, including any optional libraries, modules, programs (mini apps), and code templates.  

- libraries/:             directory for optional libraries, such as BLAS, LAPACK & PETSc  
- modules/common:         global variable definitions  
- modules/dec:            geometry and Discrete Exterior Calculus (DEC) subroutines  
- modules/io:             I/O subroutines  
- modules/math:           general math subroutines  
- modules/mpi:            MPI subroutines  
- modules/phy_darcy_flow: Darcy flow specific subroutines  
- modules/solvers:        Integration of PETSc solvers (currently only KSP solvers)  
- programs/darcy_flow/:   MiniApps for steady Darcy flow (with cracking)  
- templates:              Coding standard - templates  

utils/
------
Software tools created for use with ParaGEMS

make-paragems
-------------
Script to compile ParaGEMS code (see installation instruction below)

Execution (no ParaFEM integration)
==================================

>In this section [PARAGEMS_HOME] is used as a placeholder for your root ParaGEMS directory. If following on from the installation instructions where it is assumed that ParaGEMS is cloned into your home directory, this would be ~/paragems

Test case: steady Darcy flow through icosahedron in parallel with 2 processes  

Directory: [PARAGEMS_HOME]/examples/ICOS  

Files:  

- input.param    (input): input file  
- icos.1.node    (mesh): node location of mesh  
- icos.1.edge    (mesh): line definition of mesh (node indices)  
- icos.1.face    (mesh): face definition of mesh (node indices) including boundary conditions  
- icos.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- paragems.log\* (output): solver log file  
- solution.m     (output): solution vector in MATLAB format (order associated with parallel partitioning)  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters in MATLAB format (order associated with parallel partitioning)  
- solution.vtk   (output): solution output file (can be opened with ParaView)  

Execution:  

    cd [PARAGEMS_HOME]/examples/ICOS  
    mpirun -np 2 [PARAGEMS_HOME]/bin/darcy_2f input.param  

> **or** if [PARAGEMS_HOME]/bin/ is in your PATH environment variable:  

    mpirun -np 2 darcy_2f input.param

----

Test case: steady Darcy flow through unit cube in parallel with 4 processes (optional use of simple load balancing tool)  

Directory: [PARAGEMS_HOME]/examples/CUBE  

Files:  

- input.param    (input): input file  
- cube.poly      (mesh): mesh definition for use with TetGen software  
- cube.1.node    (mesh): node location of mesh  
- cube.1.edge    (mesh): line definition of mesh (node indices)  
- cube.1.face    (mesh): face definition of mesh (node indices) with boundary conditions  
- cube.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- cube.1.\*.old  (mesh): mesh files before load balancing  
- paragems.log\* (output): solver log file  
- solution.m     (output): solution vector in MATLAB format (order associated with parallel partitioning)  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters in MATLAB format (order associated with parallel partitioning)  
- solution.vtk   (output): solution output file (can be opened with ParaView)  

Execution:  

    cd [PARAGEMS_HOME]/examples/CUBE  
    mpirun -np 4 darcy_2f input.param  
    simpleloadbalance  
      cube.1  
    mpirun -np 4 darcy_2f input.param  
> Assuming [PARAGEMS_HOME]/bin/ is in your PATH environment variable  

----

Test case: steady Darcy flow through unit cube with deterministic cracking in parallel with 4 processes (10 cracks)

Directory: [PARAGEMS_HOME]/examples/CUBE_DCrack

Files:  

- input.param   (input): input file  
- cube.poly      (mesh): mesh definition for use with TetGen software  
- cube.1.node    (mesh): node location of mesh  
- cube.1.edge    (mesh): line definition of mesh (node indices)  
- cube.1.face    (mesh): face definition of mesh (node indices) with boundary conditions  
- cube.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- cube.1.\*.old  (mesh): mesh files before load balancing  
- paragems.log\* (output): solver log file  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters in MATLAB format (order associated with parallel partitioning)  
- unsteady\*.vtk (output): solution output file (can be opened with ParaView)  
- unsteady\*.m   (output): solution vector in MATLAB format (order associated with parallel partitioning)  
- unsteady\*press.m (output): interpolated pressure solution vector in MATLAB format (order associated with parallel partitioning) - face values scaled by area  
- unsteady.log   (output): log of cracking process  

Execution:  

    cd [PARAGEMS_HOME]/examples/CUBE_DCrack  
    simpleloadbalance  
      cube.1  
    mpirun -np 4 darcy_crkp_2f input.param  
> Assuming [PARAGEMS_HOME]/bin/ is in your PATH environment variable  

----

Test case: steady Darcy flow through unit cube with stochastic cracking in parallel with 4 processes (100 random cracks in 10 steps)

Directory: [PARAGEMS_HOME]/examples/CUBE_SCrack

Files:  

-  input.param   (input): input file  
- cube.poly      (mesh): mesh definition for use with TetGen software  
- cube.1.node    (mesh): node location of mesh  
- cube.1.edge    (mesh): line definition of mesh (node indices)  
- cube.1.face    (mesh): face definition of mesh (node indices) with boundary conditions  
- cube.1.ele     (mesh): element (cell) definition of mesh (node indices)  
- cube.1.\*.old  (mesh): mesh files before load balancing  
- paragems.log\* (output): solver log file  
- solution_vol.m (output): volumes of geometric entities in MATLAB format (order associated with parallel partitioning)  
- solution_xyz.m (output): x,y,z locations of geometric entities' circumcenters in MATLAB format (order associated with parallel partitioning)  
- unsteady\*.vtk (output): solution output file (can be opened with ParaView)  
- unsteady\*.m   (output): solution vector in MATLAB format (order associated with parallel partitioning)  
- unsteady\*press.m (output): interpolated pressure solution vector in MATLAB format (order associated with parallel partitioning) - face values scaled by area  
- unsteady.log   (output): log of cracking process  

Execution:  

    cd [PARAGEMS_HOME]/examples/CUBE_SCrack  
    simpleloadbalance  
      cube.1  
    mpirun -np 4 darcy_crkp_2f input.param
> Assuming [PARAGEMS_HOME]/bin/ is in your PATH environment variable  
