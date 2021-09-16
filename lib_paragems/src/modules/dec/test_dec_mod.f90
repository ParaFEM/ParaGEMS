!
!===============================================================================
!-- Geometry, DEC and Forman Module Tests
!> Tests for geometry_mod, dec_mod, forman_mod
!===============================================================================
! /****m*/ /src/modules/dec/test_dec_mod
MODULE test_dec_mod
!* PURPOSE
!*   Tests for geometry_mod, dec_mod, forman_mod
!* INCLUDES
!*   geometry_mod
!*   dec_mod
!*   forman_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   geometry_tests(cnt)
!*   dec_tests(cnt)
!*   forman_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for geometry_mod, dec_mod, forman_mod
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE geometry_mod
  USE dec_mod
  USE forman_mod
  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/dec/test_dec_mod|geometry_test
  SUBROUTINE geometry_tests(cnt)
!* PURPOSE
!*   tests geometry_mod
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
!> tests geometry_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT) :: cnt             ! test counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------------------------------------------------------------
    !*   initialise_geo          initialise geometric quantities for the given mesh
    !*   get_lcl_node_indx       create map between global and local node indices
    !*   calc_circumcenters      compute circumcenter of given elements
    !*   calc_prml_sgnd_vlm      compute signed volume of primal elements
    !*   calc_prml_unsgnd_vlm    compute unsigned volume of primal elements
    !*   calc_dual_vlm           compute volume for dual elements of all geometric order
    !*   calc_dual_vlm_i         recursively add to dual volume calculation
    !*   add_points_i            recursively adds points for dual volume calculation
    !*   calc_unsgnd_vlm         compute unsigned volume for given set of points
    !*   exchange_dual_vlm       exchange external dual volumes between adjacent processes
    !*   exchange_dual_dir       exchange external dual edge directions between adjacent processes
    !*   calc_prml_unsgnd_vlm    compute unsigned volume of primal elements
    !*   calc_prml_dir           compute the unit direction of primal edges
    !*   calc_dual_dir           compute the unit direction of dual edges
    !---------------------------------------------------------------------

  END SUBROUTINE
! geometry_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/dec/test_dec_mod|dec_test
  SUBROUTINE dec_tests(cnt)
!* PURPOSE
!*   tests dec_mod
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
!> tests dec_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT) :: cnt             ! test counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------------------------------------------------------------
    !*   calc_orientation        sort nodal indices and compute ±1 orientation
    !*   calc_bndry_cobndry      recursively compute element (co-)boundaries
    !*   build_bndry_work_array  build boundary data working array
    !*   count_bndry_cobndry     count [co-]boundaries: internal, external, surface
    !*   set_bndry_cobndry       set [co-]boundaries: internal, external, surface
    !*   calc_hodge_star         compute hodge star and it's inverse from primal and dual volumes
    !---------------------------------------------------------------------

  END SUBROUTINE
! dec_mod_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/dec/test_dec_mod|forman_test
  SUBROUTINE forman_tests(cnt)
!* PURPOSE
!*   tests forman_mod
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
!> tests forman_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT) :: cnt             ! test counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------------------------------------------------------------
    !*
    !---------------------------------------------------------------------

  END SUBROUTINE
! forman_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! dec_test
!===============================================================================
!
