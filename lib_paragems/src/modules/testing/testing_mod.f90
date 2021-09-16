!
!===============================================================================
!-- Testing Module
!> Functions and subroutines for testing
!===============================================================================
!/****/h* tests|testing_mod
!* SYNOPSIS
MODULE testing_mod
!* PURPOSE
!*   Functions and subroutines for testing
!* CONTAINS
!*   Subroutine              Purpose
!*   common_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> Functions and subroutines for testing
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE common_mod
  USE mpi_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/test/testing_mod|chkerr_unit
!* SYNOPSIS
  SUBROUTINE chkerr_unit(cnt,err,msg)
!* PURPOSE
!*   error output for unit test
!* INPUTS
!*   test counter
!* OUTPUTS
!*   test counter
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> error output for unit test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT)  :: cnt  !> test counter
    LOGICAL, INTENT(IN)     :: err  !> error logical
    CHARACTER(LEN=*)        :: msg  !> test msg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- output result
    IF (err) THEN;  cnt=cnt+1
      WRITE(*,'(I5,A)') cnt,' : XXXXXX  FAILED : '//TRIM(msg)
    ELSE;  WRITE(*,'(I5,A)')cnt,' : v/      PASSED : '//TRIM(msg);  END IF

    RETURN

  END SUBROUTINE
! testing_mod|chkerr_unit
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/test/testing_mod|chkerr_unit_MPI
!* SYNOPSIS
  SUBROUTINE chkerr_unit_MPI(cnt,err,msg)
!* PURPOSE
!*   error output for unit test
!* INPUTS
!*   test counter
!* OUTPUTS
!*   test counter
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> error output for unit test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT)  :: cnt  !> test counter
    LOGICAL, INTENT(IN)     :: err  !> error logical
    CHARACTER(LEN=*)        :: msg  !> test msg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    !-- output result
    IF (rank==root) THEN
      IF (err) THEN;  cnt=cnt+1
        WRITE(*,'(I5,A)') cnt,' : XXXXXX  FAILED : '//TRIM(msg)
      ELSE;  WRITE(*,'(I5,A)')cnt,' : v/      PASSED : '//TRIM(msg);  END IF
    END IF

    RETURN

  END SUBROUTINE
! testing_mod|chkerr_unit_MPI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! testing
!===============================================================================
!
