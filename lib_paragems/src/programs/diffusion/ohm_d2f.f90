!
!===============================================================================
!-- Ohm flow - two-field formulation
!> Miniapp to solve two-field Ohm flow in parallel using PETSc KSP
!===============================================================================
!/****p* programs|darcy_flow/ohm_d2f
!* SYNOPSIS
PROGRAM ohm_d2f
!* PURPOSE
!*   Miniapp to solve two-field Ohm flow in parallel using PETSc KSP
!* ASSUMPTION
!*   Mesh file exists (TetGen tetrahedral mesh)
!*   Input file (default: input.param)
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 general mpi routines (start, end, syncwrite,etc)
!*   partition_mod           parallel partitioning
!*   darcy_mod               Ohm flow specific routines
!*   solver_mod              Solver routines (PETSc)
!* SIDE EFFECTS
!*   Solution to Ohm flow computed
!*   Solution written to file
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Miniapp to solve two-field Ohm flow in parallel using PETSc KSP
!===============================================================================

  !-- modules and implicit none --
  USE common_mod
  USE mpi_mod
  USE geometry_mod
  USE partition_mod
  USE darcy_mod
  USE solver_mod

  IMPLICIT NONE

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- start MPI processes and write miniapp details --
  CALL start_petsc_mpi()
  CALL rootwrite_log('ParaGEMS: Solving two field Ohm flow in parallel '//&
    'using PETSc KSP')

  !-- read input --
  CALL read_input_darcy()

  !-- setup partitioning and connectivity --
  !-- includes: DEC operators (topology)
  CALL setup_partitioning();  CALL calc_bndry_cobndry();  CALL setup_connectivity()

  !-- initialise geometry --
  !-- includes: DEC operators (geometry)
  CALL read_nodes_prml();  CALL exchange_prml_nodes()
  CALL initialise_geo();  CALL initialise_hodge_star()

  !-- initialise --
  CALL initialise_darcy_d();  CALL set_IC_darcy_d();  CALL setup_LHS_darcy_d()
  CALL setup_RHS_darcy_d();  CALL setup_IC_BC_darcy_d();  CALL get_RHS_darcy_d()

  !-- solve linear system --
  CALL solve_KSP_SC()

  !-- MATLAB output of solution and pressures
  CALL update_solution();  CALL extract_sol_KSP_d()
  CALL write_solution_MATLAB();  CALL write_pressure_MATLAB_d()

  !-- MATLAB output of centers and volumes
  CALL write_centers_MATLAB_d();  CALL write_dual_volumes_MATLAB()

  !-- write solution (ohm) --
  CALL syncwrite_log('> Whitney interpolcation()');  CALL calc_whitney_C2_BC()
  CALL syncwrite_log('> write_solution()')
  CALL write_solution_D0S2('pressure','velocity')

  !-- clean up allocated variables and end MPI --
  CALL syncwrite_log('ParaGEMS: cleaning up Ohm flow solution')
  CALL deallocate_lcl_complex();  CALL finalise_darcy();  CALL end_petsc_mpi()

END PROGRAM
! ohm
!===============================================================================
!
