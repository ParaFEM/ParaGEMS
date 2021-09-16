!
!===============================================================================
!-- Darcy flow - two-field formulation
!> Miniapp to solve two-field Darcy flow in parallel using PETSc KSP
!===============================================================================
!/****p* programs|darcy_flow/darcy
!* SYNOPSIS
PROGRAM darcy
!* PURPOSE
!*   Miniapp to solve two-field Darcy flow in parallel using PETSc KSP
!* ASSUMPTION
!*   Mesh file exists (TetGen tetrahedral mesh)
!*   Input file (default: input.param)
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 ???
!*   partition_mod
!*   darcy_mod
!*   solver_mod
!* SIDE EFFECTS
!*   Solution to Darcy flow computed
!*   Solution written to file
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Miniapp to solve two-field Darcy flow in parallel using PETSc KSP
!===============================================================================

  !-- modules and implicit none --
  USE common_mod
  USE mpi_mod
  USE partition_mod
  USE darcy_mod
  USE solver_mod
  IMPLICIT NONE

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- start MPI processes and write miniapp details --
  CALL start_mpi()
  CALL rootwrite_log('ParaGEMS: Solving two field Darcy flow in parallel using PETSc KSP')

  !-- read input --
  CALL read_input_darcy()

  !-- setup partitioning and connectivity --
  CALL parallel_setup()

  !-- initialise geometric information --
  CALL initialise_geo()

  !-- initialise solution --
  CALL start_petsc()
  CALL initialise_darcy()

  !-- get necessary matrices and vectors (A x = b) --
  CALL get_LHS_darcy()
  CALL get_RHS_darcy()

  !-- solve linear system --
  CALL solve_KSP()

  !-- write solution (darcy) --
  CALL syncwrite_log('> Whitney interpolcation()');  CALL calc_whitney_C2_BC()
  CALL syncwrite_log('> write_solution()');  CALL write_solution_D0S('pressure','velocity')

  !-- clean up allocated variables and end MPI --
  CALL syncwrite_log('ParaGEMS: cleaning up Darcy flow solution')
  CALL clean_up()
  CALL finalise_darcy()
  CALL end_petsc()
  CALL end_mpi()

END PROGRAM
! darcy
!===============================================================================
!
