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
!*   Name                     Purpose
!*   petsc.h                  petsc library (includes mpi execution env library)
!*   global_precision.inc     precision of real numbers
!*   global_vars.inc          general global variables and parameters
!*   global_mpi.inc           global mpi variables
!*   global_io.inc            global IO variables
!*   global_mesh.inc          global mesh variables
!*   global_element.inc       global primary structures
!*   global_petsc.inc         global petsc variables
!* CONTAINS
!*   Subroutine               Purpose
!*   deallocate_lcl_complex   deallocated global structures/variables
!*   chkerr                   check for errors, print message and exit gracefully
!*
!*   Function                 Purpose
!*   chkerr_log               check for errors, print message and return logical
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
!/****/s* common_mod/deallocate_lcl_complex
!* SYNOPSIS
  SUBROUTINE deallocate_lcl_complex()
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
    IF (ALLOCATED(lcl_complex)) DEALLOCATE(lcl_complex); RETURN

  END SUBROUTINE
! deallocate_lcl_complex
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* common_mod/chkerr
!* SYNOPSIS
  SUBROUTINE chkerr(ierr,str,io_unit)
!* PURPOSE
!*   check for errors, print message if error, and exit gracefully
!* INPUTS
!*   Name                    Description
!*   ierr                    error integer
!*   str                     error message string to be written to file
!*   io_unit                 io unit for outpute file (assumed open)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> check for errors, print message if error, and exit gracefully
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)           :: ierr     !> error integer
    INTEGER, INTENT(IN)           :: io_unit  !> io unit
    CHARACTER(LEN=*), INTENT(IN)  :: str      !> error string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- check for errors and exit graacefully if any found --
    IF (ierr /= 0) THEN; WRITE(io_unit,'(A,I5)') 'Error :: '//TRIM(str)//&
      ' : errcode = ',ierr; FLUSH(io_unit); STOP 'ParaGEMS: shutdown'; END IF
    RETURN

  END SUBROUTINE
! chkerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f* common_mod/chkerr_log
!* SYNOPSIS
  FUNCTION chkerr_log(ierr,str,io_unit)
!* PURPOSE
!*   check for errors, print message if error, and return logical
!*   Name                    Description
!*   ierr                    error integer
!*   str                     error message string to be written to file
!*   io_unit                 io unit for outpute file (assumed open)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> check for errors and return logical
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments and result --
    LOGICAL                       :: chkerr_log   !> error logical (return)
    INTEGER, INTENT(IN)           :: ierr         !> error integer
    INTEGER, INTENT(IN)           :: io_unit      !> io unit
    CHARACTER(LEN=*), INTENT(IN)  :: str          !> error string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- check for errors and exit graacefully if any found --
    IF (ierr /= 0) THEN; WRITE(io_unit,'(A,I5)') 'Error :: '//TRIM(str)//&
      ' : errcode = ',ierr; chkerr_log = .TRUE.; ELSE; chkerr_log = .FALSE.
    END IF;  RETURN

  END FUNCTION
! chkerr_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
END MODULE
! common_mod
!===============================================================================
!
