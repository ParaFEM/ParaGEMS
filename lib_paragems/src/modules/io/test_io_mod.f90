!
!===============================================================================
!-- io Module Tests
!> Tests for io_mod module
!===============================================================================
!/****/h* modules|io/test_io_mod
!* SYNOPSIS
MODULE test_io_mod
!* PURPOSE
!*   Tests for io_mod module
!* INCLUDES
!*   io_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   io_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for io_mod module
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE io_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/io/test_io_mod|io_test
!* SYNOPSIS
  SUBROUTINE io_tests(cnt)
!* PURPOSE
!*   tests io_mod
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
!> tests io_mod
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
    !*   read_prml_elms          read local node indices of highest order primal elements
    !*   read_glb_indx
    !*   read_nodes_prml         read the nodal locations of the primal mesh
    !*   root_open_file_read
    !*   root_open_file_write
    !*   write_solution_D0S
    !*   write_solution_D0S2
    !*   write_unsteady_D0S
    !*   write_unsteady_D0S2
    !*   read_crack_faces
    !*   read_bndry_cond
    !*   read_bndry_cond2
    !*   init_write_MATLAB
    !*   clean_PETSc_output
    !*   write_centers_MATLAB2
    !*   write_prml_volumes_MATLAB2
    !*   write_solution_MATLAB2
    !*   write_pressure_MATLAB2
    !*   write_centers_MATLAB
    !*   write_prml_volumes_MATLAB
    !*   write_solution_MATLAB
    !*   write_pressure_MATLAB
    !---------------------------------------------------------------------

  END SUBROUTINE
! test_io_mod|io_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! io_test
!===============================================================================
!
