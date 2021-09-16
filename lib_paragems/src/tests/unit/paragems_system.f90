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
  USE common_mod;         USE test_common_mod
  !USE geometric_mod
  !USE forman_mod
  USE dec_mod;            USE test_dec_mod
  USE io_mod;             USE test_io_mod
  USE math_mod;           USE test_math_mod
  USE mpi_mod;            USE test_mpi_mod
  USE partition_mod;      USE test_partition_mod
  USE darcy_mod;          USE test_darcy_mod
  USE solver_mod;         USE test_solver_mod


  IMPLICIT NONE

  !---------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------
  INTEGER :: cnt=0              ! test counter

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  CALL common_test(cnt)
  CALL geometric_test(cnt)
  CALL dec_test(cnt)
  CALL forman_test(cnt)
  CALL io_test(cnt)
  CALL math_test(cnt)
  CALL mpi_test(cnt)
  CALL partition_test(cnt)
  CALL darcy_test(cnt)
  CALL solver_test(cnt)

END PROGRAM
! paragems_unit
!===============================================================================
!
