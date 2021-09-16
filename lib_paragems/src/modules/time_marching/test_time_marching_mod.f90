!
!===============================================================================
!-- time_marching Module Tests
!> Tests for time_marching_mod module
!===============================================================================
!/****/h* modules|time_marching/test_time_marching_mod
!* SYNOPSIS
MODULE test_time_marching_mod
!* PURPOSE
!*   Tests for time_marching_mod module
!* INCLUDES
!*   time_marching_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   time_marching_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for time_marching_mod module
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE time_marching_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/time_marching/test_time_marching_mod|time_marching_test
!* SYNOPSIS
  SUBROUTINE time_marching_tests(cnt)
!* PURPOSE
!*   tests time_marching_mod
!* INPUTS
!*   test counter
!* OUTPUTS
!*   test counter
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> tests time_marching_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT) :: cnt             !> test counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !---------------------------------------------------------------------
    ! 1. precision
    !---------------------------------------------------------------------
    cnt = cnt + 1
    IF (1.d0 - 9.999999999999999d-1 > 0.d0) THEN
      WRITE(*,'(A,I5,A)')'time_marching_test : ',cnt,&
      ' : v/      PASSED : precision accurate to 16 decimal places'
    ELSE
      WRITE(*,'(A,I5,A)')'time_marching_test : ',cnt,&
      ' : XXXXXX  FAILED : precision NOT accurate to 16 decimal places'
    END IF

    !---------------------------------------------------------------------
    ! DO NOTHING
    !------------
    ! 2. global variables
    ! 3. global mpi
    ! 4. global io
    ! 5. global mesh
    ! 6. primary element
    ! 7. variable clean up
    !---------------------------------------------------------------------

  END SUBROUTINE
! test_time_marching_mod|time_marching_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! time_marching_test
!===============================================================================
!
