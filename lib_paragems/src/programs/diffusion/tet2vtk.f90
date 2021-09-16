!
!===============================================================================
!-- Ohm flow - two-field formulation
!> Miniapp to solve two-field Ohm flow in parallel using PETSc KSP
!===============================================================================
!/****p* programs|darcy_flow/tet2vtk_2f
!* SYNOPSIS
PROGRAM tet2vtk_2f
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
  CALL initialise_darcy();
  !CALL set_IC_darcy();  CALL setup_LHS_darcy()
  ! CALL setup_RHS_darcy();  CALL setup_IC_BC_darcy();  CALL get_RHS_darcy()
  !
  ! !-- solve linear system --
  ! CALL solve_KSP_SC()
  !
  ! !-- MATLAB output of solution and pressures
  ! CALL update_solution();  CALL extract_sol_KSP()
  ! CALL write_solution_MATLAB();  CALL write_pressure_MATLAB()
  !
  ! !-- MATLAB output of centers and volumes
  ! CALL write_centers_MATLAB();  CALL write_prml_volumes_MATLAB()

  !-- write solution (tet2vtk) --
  ! CALL syncwrite_log('> Whitney interpolcation()');  CALL calc_whitney_C2_BC()
  CALL syncwrite_log('> write_solution()')
  CALL write_solution_D0S2('pressure','velocity')

  !-- clean up allocated variables and end MPI --
  CALL syncwrite_log('ParaGEMS: cleaning up Ohm flow solution')
  CALL deallocate_lcl_complex();
  CALL end_petsc_mpi()

END PROGRAM
! tet2vtk
!===============================================================================
!
