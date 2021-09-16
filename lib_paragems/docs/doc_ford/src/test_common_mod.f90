!
!===============================================================================
!-- Common Module Tests
!> Tests for common_mod module, primairly global variables
!===============================================================================
!/****/h* modules|common/test_common_mod
!* SYNOPSIS
MODULE test_common_mod
!* PURPOSE
!*   Tests for common_mod module, primairly global variables
!* INCLUDES
!*   common_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   test_vars(cnt)          - deallocated global structures/variables
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for common_mod module, primairly global variables
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE common_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/common/common_test|test_vars
!* SYNOPSIS
  SUBROUTINE test_vars(cnt)
!* PURPOSE
!*   test global structures/variables
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
!> deallocate global structures/variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TO DO:
! - global variables?
! - primary element data structure
! - clean up
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
      WRITE(*,'(A,I5,A)')'common_test : ',cnt,&
      ' : v/      PASSED : precision accurate to 16 decimal places'
    ELSE
      WRITE(*,'(A,I5,A)')'common_test : ',cnt,&
      ' : XXXXXX  FAILED : precision NOT accurate to 16 decimal places'
    END IF

    !---------------------------------------------------------------------
    ! 2. global variables
    ! 3. global mpi
    ! 4. global io
    ! 5. global mesh
    ! 6. primary element
    ! 7. variable clean up
    !---------------------------------------------------------------------

  END SUBROUTINE
! common_test|test_vars
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! common_test
!===============================================================================
!
