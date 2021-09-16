
  *** READ ONLY PROGRAM ***

  PROGRAM: pg125.f90

  pg125 performs the 3D explicit analysis of the transient conduction equation.
  The implementation is based on discrete exterior calculus through the ParaGEMS
  library.

    Usage: pg125 <job_name>

  ERRATUM

  In the book, the argument NUMVAR in the subroutine SCATTER_NODES (line 110)
  is given the value NDIM(=3). This is incorrect and has been replaced with
  NODOF(=1).

  AUTHOR

  pieter.boom@manchester.ac.uk
  lee.margetts@manchester.ac.uk
