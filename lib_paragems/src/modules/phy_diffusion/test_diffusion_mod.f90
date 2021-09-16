!
!===============================================================================
!-- diffusion Module Tests
!> Tests for diffusion_mod module
!===============================================================================
!/****/h* modules|diffusion/test_diffusion_mod
!* SYNOPSIS
MODULE test_diffusion_mod
!* PURPOSE
!*   Tests for diffusion_mod module
!* INCLUDES
!*   diffusion_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   diffusion_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for diffusion_mod module
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE diffusion_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/diffusion/test_diffusion_mod|diffusion_test
!* SYNOPSIS
  SUBROUTINE diffusion_tests(cnt)
!* PURPOSE
!*   tests diffusion_mod
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
!> tests diffusion_mod
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
    !*   read_input_diffusion        read inputs, and check parameters (phys/io/solver)
    !*   check_param_diffusion       check parameters for diffusion simulations
    !*   initialise_diffusion
    !*   finalise_diffusion
    !*   get_RHS_diffusion
    !*   get_LHS_diffusion
    !*   identify_crack
    !*   identify_crack2
    !*   identify_crack3
    !*   exchange_bndry_cond
    !*   initialise_diffusion2
    !*   finalise_diffusion2
    !*   get_RHS_diffusion2
    !*   get_LHS_diffusion2
    !*   identify_crack4
    !*   identify_crack5
    !---------------------------------------------------------------------

  END SUBROUTINE
! test_diffusion_mod|diffusion_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! diffusion_test
!===============================================================================
!
