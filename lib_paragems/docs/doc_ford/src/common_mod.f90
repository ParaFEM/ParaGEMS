!
!===============================================================================
!-- Common Module
!> Module for global variables, included library, and common code
!===============================================================================
!/****/h* modules|comon/common_mod
!* SYNOPSIS
MODULE common_mod
!* PURPOSE
!*   Module for global variables, included library, and common code
!* INCLUDES
!*   Name                  Purpose
!*   petsc.h                 - petsc library (includes mpi execution
!*                             environment library)
!*   global_precision.inc    - precision of real numbers
!*   global_vars.inc         - general global variables and parameters
!*   global_mpi.inc          - global mpi variables
!*   global_io.inc           - global IO variables
!*   global_mesh.inc         - global mesh variables
!*   global_element.inc      - global primary structures
!* CONTAINS
!*   Subroutine              Purpose
!*   clean_up()              - deallocated global structures/variables
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Module for global variables, included library, and common code
!===============================================================================

  !-- implicit none, and include files for variables and libraries --
#include <petsc/finclude/petsc.h>
  USE petsc

  IMPLICIT NONE

#include "global_precision.inc"
#include "global_vars.inc"
#include "global_mpi.inc"
#include "global_io.inc"
#include "global_mesh.inc"
#include "global_element.inc"
#include "global_petsc.inc"

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* common_mod/clean_up
!* SYNOPSIS
  SUBROUTINE clean_up()
!* PURPOSE
!*   Deallocate global structures/variables
!* SIDE EFFECTS
!*   Global structure (lcl_complex) deallocated
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Deallocate global structures/variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- deallocate primary elements --
    IF (ALLOCATED(lcl_complex)) DEALLOCATE(lcl_complex)

    RETURN

  END SUBROUTINE
! clean_up
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
END MODULE
! common_mod
!===============================================================================
!
