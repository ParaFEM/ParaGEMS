!
!===============================================================================
!-- ParaGEMS Unit testing
!> Driver routine for ParaGEMS' unit tests
!===============================================================================
!/****p* tests/paragems_unit
!* SYNOPSIS
PROGRAM paragems_unit
!* PURPOSE
!*   Driver routine for ParaGEMS' unit tests
!* ASSUMPTION
!*
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!* SIDE EFFECTS
!*
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Driver routine for ParaGEMS' unit tests
!===============================================================================

  !-----------------------------------------------------------------------------
  ! modules and implicit none
  !-----------------------------------------------------------------------------
  USE common_mod;       USE test_common_mod

  IMPLICIT NONE

  !---------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------
  INTEGER :: cnt=0              ! test counter

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-----------------------------------------------------------------------
  !  1. common_test
  !-----------------------------------------------------------------------
  CALL test_vars(cnt)

END PROGRAM
! paragems_unit
!===============================================================================
!
