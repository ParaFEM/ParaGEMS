
  *** READ ONLY PROGRAM ***

  PROGRAM: pg124.f90

  pg124 performs the 3D analysis of the transient heat conduction equation. The
  implementation is based on discrete exterior calculus through the ParaGEMS
  library.

    Usage: pg124 <job_name>

  ERRATUM

  On some platforms, the original program fails because loaded_freedoms_pp
  has not been set to zero. The program has been updated on line 132:

    IF(loaded_freedoms==0) loaded_freedoms_pp=0

  AUTHOR

  pieter.boom@manchester.ac.uk
  lee.margetts@manchester.ac.uk
