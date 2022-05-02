!
!===============================================================================
!-- elasticity Module Tests
!> Tests for elasticity_mod module
!===============================================================================
!/****/h* modules|elasticity/test_elasticity_mod
!* SYNOPSIS
MODULE test_elasticity_mod
!* PURPOSE
!*   Tests for elasticity_mod module
!* INCLUDES
!*   elasticity_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   elasticity_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for elasticity_mod module
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE elasticity_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/elasticity/test_elasticity_mod|elasticity_test
!* SYNOPSIS
  SUBROUTINE elasticity_tests(cnt)
!* PURPOSE
!*   tests elasticity_mod
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
!> tests elasticity_mod
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
    !*   Subroutine              Purpose
    !---------------------------------------------------------------------

  END SUBROUTINE
! test_elasticity_mod|elasticity_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! elasticity_test
!===============================================================================
!
