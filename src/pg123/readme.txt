
  *** READ ONLY PROGRAM ***

  PROGRAM: pg123.f90

  pg123 performs the three dimensional analysis of Laplace's equation using
  8-node bricks. The implementation is based on discrete exterior calculus
  through the ParaGEMS library.

    Usage: pg123 <job_name>

  ERRATUM

  In the book, the argument NUMVAR in the subroutine SCATTER_NODES (line 176)
  is given the value NDIM(=3). This is incorrect and has been replaced with
  NODOF(=1).

  AUTHOR

  pieter.boom@manchester.ac.uk
  lee.margetts@manchester.ac.uk
