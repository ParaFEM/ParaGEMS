Description
===========

Integration of ParaGEMS and ParaFEM libraries, along with example miniApps. From this point forward, the integrated libraries will be referred to as ParaFEM+GEMS

Authors: Dr Pieter Boom (pieter.boom@manchester.ac.uk); Dr Lee Margetts (lee.margetts@manchester.ac.uk)

Updated repository: https://github.com/ParaFEM/ParaGEMS

Installation
============

The following section describes the installation of the ParaFEM+GEMS on a Linux machine, with notes for MacOS and Windows Subsystem for Linux (WSL).

Prerequisites
-------------

A standard compilation of ParaFEM is required, along with a "--for-parafem" compilation of the ParaGEMS library found in the lib_paragems subdirectory. It is expected that this software is compiled in a subdirectory of the main ParaFEM installation:

    cd  [Root ParaFEM directory]
    ls
        README.md          demo-installer-win fruit              
        fruit_3.4.3        parafem            parafem-viewer        
        paragems

Installation instructions for these libraries is found in their respective subdirectories.

The software has been installed and tested on both Linux (Xubuntu 20.04 LTS, CentOS 7.9.2009) and MacOS (High Sierra, Catalina, and Big Sur). Other operating systems are NOT supported.

ParaFEM+GEMS has been compiled and tested using using the following compilers and libraries (Other compilers and libraries are NOT supported):

 - Intel ifort (17.0 and 18.0), GCC gfortran (8, 9 and 10), Cray Compiling Environment (CCE), and AMD compiler environment (AOCC)
 - OpenMPI 4.0.1, MPICH (8.1.4)

Other software needed:

 - Git (see section: Getting the software)
 - BLAS and LAPACK, which can be installed on Ubunutu using the command:

    sudo apt-get install libblas-dev liblapack-dev


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

First clone the ParaFEM library via ssh (requires a git hub account and an SSH key):

    git clone git@github.com:ParaFEM/ParaFEM.git

or via https:

    git clone https://github.com/ParaFEM/ParaFEM.git

Next, enter the ParaFEM directory and clone the ParaGEMS integration via ssh (requires a git hub account and an SSH key):

    cd parafem
    ls
        README.md          demo-installer-win fruit              
        fruit_3.4.3        parafem            parafem-viewer
    git clone git@github.com:ParaFEM/ParaGEMS.git

or via https:

    cd parafem
    ls
        README.md          demo-installer-win fruit              
        fruit_3.4.3        parafem            parafem-viewer
    git clone https://github.com/ParaFEM/ParaGEMS.git



Compiling ParaFEM+GEMS
----------------------

Navigate to the root ParaFEM directory:

    cd ~/parafem/parafem

>This assumes the root ParaFEM directory is in your home directory. If not you will need to update the path displayed above

In the build/ subdirectory there are include files for a few different systems, including for a standard linux desktop: linuxdesktop.inc.

>Existing include files (.inc) can be adapted for the particulars of your own system.

ParaFEM can now be compiled with the command:

    ./make-parafem MACHINE=linuxdesktop

Next, navigate the paragems library directory:

    cd ~/parafem/paragems/lib_paragems

>This assumes the root ParaFEM directory is in your home directory. If not you will need to update the path displayed above

In the config/ subdirectory there are include files for a few different systems, including for a standard linux desktop: linuxdesktop_lib.inc.

>Existing include files (.inc) can be adapted for the particulars of your own system.

ParaGEMS can now be compiled with the command:

    ./make-paragems MACHINE=linuxdesktop_lib --for-parafem

Finally, return to the ParaGEMS integration directory:

    cd ~/parafem/paragems/

>This assumes the root ParaFEM directory is in your home directory. If not you will need to update the path displayed above

In the build/ subdirectory there are include files for a few different systems, including for a standard linux desktop: linuxdesktop.inc.

>Existing include files (.inc) can be adapted for the particulars of your own system.

ParaFEM+GEMS can now be compiled with the command:

    ./make-para-fem-gems MACHINE=linuxdesktop

Setting the PATH environment variable for ParaFEM+GEMS
------------------------------------------------------

To make the ParaGEMS and utility executables available while in any directory and every time you open a new terminal window, you can add the parafem/paragems/bin directory to your PATH environment variable by adding:

    export PATH=~/parafem/paragems/bin:$PATH

>This assumes the root ParaGEMS directory is in your home directory. If not you will need to update the path displayed above

to the end of your ~/.bashrc file (for bash shell) and sourcing it:

    source ~/.bashrc

ParaGEMS File/Directory Structure
=================================

bin/
----
Directory for compiled executables created at time of compilation

build/
------
Directory for system specific configuration files (includes template)

examples/
---------
Directory with sample test cases  

- pg121:        linear elasticity  
- pg123:        steady Laplace equation  
- pg123x:       steady Laplace equation with emerging discontinuities  
- pg124:        implicit time dependent hear flow  
- pg125:        explicit time dependent hear flow  

include/
--------
Directory of compiled module files (.mod) created at time of compilation

lib/
----
Directory of compiled library files (.a) created at time of compilation

lib_paragems/
-------------
Directory of paragems library

license/
--------
Directory with software license

src/
----
Directory with source code (mini apps):

- pg121/        linear elasticity  
- pg123/        steady Laplace equation  
- pg123x/       steady Laplace equation with emerging discontinuities  
- pg124/        implicit time dependent hear flow  
- pg125/        explicit time dependent hear flow  


make-para-fem-gems
------------------
Script to compile ParaFEM+GEMS code

Execution (ParaFEM+GEMS integration)
==================================

>In this section [PARAFEMGEMS_HOME] is used as a placeholder for your root ParaFEM+GEMS directory. If following on from the installation instructions where it is assumed that ParaGEMS is cloned into your home directory, this would be ~/parafem/paragems

Test case: Program 12.1 three dimensional analysis of an elastic solid using 8-node brick elements, DEC, preconditioned conjugate gradient solver; diagonal preconditioner diag_precon; parallel version loaded_nodes only (note that ParaFEM's p121 uses 20-node bricks)

Directory: [PARAFEMGEMS_HOME]/examples/pg121

Subdirectories and Files:  

- book/             example from textbook
- demo/             basic example file
- mg/               additional mesh files
- demo/readme.txt   description of simulation input and output files

Execution:  

    cd [PARAFEMGEMS_HOME]/examples/pg121/demo  
    mpirun -np 2 [PARAFEMGEMS_HOME]/bin/pg121 pg121_demo  

> **or** if [PARAFEMGEMS_HOME]/bin/ is in your PATH environment variable:  

    mpirun -np 2 pg121 pg121_demo

----

Test case: program p12.3 three dimensional analysis of Laplace's equation using 8-node bricks, DEC, preconditioned conjugate gradient solver diagonal preconditioner; parallel; externally generated model

Directory: [PARAFEMGEMS_HOME]/examples/pg123

Subdirectories and Files:  

- book/             example from textbook
- demo/             basic example file
- mg/               additional mesh files
- book/readme.txt   description of simulation input and output files

Execution:  

    cd [PARAFEMGEMS_HOME]/examples/pg123/demo  
    mpirun -np 2 [PARAFEMGEMS_HOME]/bin/pg123 pg123_demo

> **or** if [PARAFEMGEMS_HOME]/bin/ is in your PATH environment variable:  

    mpirun -np 2 pg123 pg123_demo

----

Test case: program p12.3 three dimensional analysis of Laplace's equation using 8-node bricks, DEC, preconditioned conjugate gradient solver diagonal preconditioner; parallel; externally generated model; evolution of faces with zero diffusivity

Directory: [PARAFEMGEMS_HOME]/examples/pg123x

Subdirectories and Files:  

- book/             example from textbook
- demo/             basic example file
- mg/               additional mesh files
- book/readme.txt   description of simulation input and output files

Execution:  

    cd [PARAFEMGEMS_HOME]/examples/pg123x/demo  
    mpirun -np 2 [PARAFEMGEMS_HOME]/bin/pg123x pg123_demo

> **or** if [PARAFEMGEMS_HOME]/bin/ is in your PATH environment variable:  

    mpirun -np 2 pg123x pg123_demo

----

Test case: program 12.4 three dimensional transient analysis of heat conduction equation using 8-node hexahedral elements; DEC, parallel pcg version implicit; integration in time using 'theta' method

Directory: [PARAFEMGEMS_HOME]/examples/pg124

Subdirectories and Files:  

- book/             example from textbook
- demo/             basic example file
- mg/               additional mesh files
- demo/readme.txt   description of simulation input and output files

Execution:  

    cd [PARAFEMGEMS_HOME]/examples/pg124/demo  
    mpirun -np 2 [PARAFEMGEMS_HOME]/bin/pg124 pg124_demo  

> **or** if [PARAFEMGEMS_HOME]/bin/ is in your PATH environment variable:  

    mpirun -np 2 pg124 pg124_demo

----

Test case: Program 12.5 conduction equation on a 3-d box volume using 8-node hexahedral elements and a simple explicit algorithm: DEC, parallel version write on processor it at freedom nres

Directory: [PARAFEMGEMS_HOME]/examples/pg125

Subdirectories and Files:  

- book/             example from textbook
- demo/             basic example file
- mg/               additional mesh files
- demo/readme.txt   description of simulation input and output files

Execution:  

    cd [PARAFEMGEMS_HOME]/examples/pg125/demo  
    mpirun -np 2 [PARAFEMGEMS_HOME]/bin/pg125 pg125_demo  

> **or** if [PARAFEMGEMS_HOME]/bin/ is in your PATH environment variable:  

    mpirun -np 2 pg125 pg125_demo
