!
!===============================================================================
!-- Ohm flow - two-field formulation
!> Miniapp to solve two-field Ohm flow in parallel using PETSc KSP
!===============================================================================
!/****p* programs|diffusion/ohm_2f
!* SYNOPSIS
PROGRAM ohm_2f
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
!*   diffusion_mod           Ohm flow specific routines
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
  USE diffusion_mod
  USE solver_mod

  IMPLICIT NONE

  !-- local variables --
  INTEGER :: k  !> simplicial order

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- start MPI processes and write miniapp details --
  CALL start_petsc_mpi(); CALL rootwrite('ParaGEMS: Solving two field '//&
    'Ohm''s law in parallel using PETSc KSP',log_unit)

  !-- read input --
  CALL rootwrite('> read_input_diffusion()',log_unit)
  CALL read_input_diffusion()

  !-- setup partitioning and connectivity - DEC operators (topology) --
  CALL syncwrite('> setup_partitioning()',log_unit)
  CALL partition_dual(); CALL read_prml_elms()
  CALL syncwrite_time(io_unit=log_unit)

  CALL syncwrite('> calc_bndry_cobndry()',log_unit)
  DO k=dim_cmplx+1,2,-1; CALL calc_bndry_cobndry(k); END DO
  CALL syncwrite_time(io_unit=log_unit)

  CALL syncwrite('> setup_connectivity()',log_unit)
  CALL get_glb_indx_dual_vlm(); CALL get_connectivity()
  CALL get_glb_indx(dim_cmplx); CALL syncwrite_time(io_unit=log_unit)

  !-- initialise geometry - DEC operators (geometry) --
  CALL syncwrite('> read_nodes_prml()',log_unit)
  CALL read_nodes_prml(); CALL exchange_prml_nodes()
  CALL syncwrite_time(io_unit=log_unit)
  CALL syncwrite('> initialise_geo()',log_unit)
  CALL initialise_geo(dim_cmplx-1,dim_cmplx+1); CALL calc_hodge_star(dim_cmplx)
  CALL syncwrite_time(io_unit=log_unit)

  !-- initialise  simulation --
  CALL syncwrite('> initialise_diffusion() and set_IC_diffusion()',log_unit)
  CALL initialise_diffusion(); CALL set_IC_diffusion()
  CALL syncwrite_time(io_unit=log_unit)
  CALL syncwrite('> setup LHS, RHS, IC, BC and get_RHS_diffusion()',log_unit)
  CALL setup_LHS_diffusion(); CALL setup_RHS_diffusion()
  CALL setup_IC_BC_diffusion(); CALL get_RHS_diffusion()
  CALL syncwrite_time(io_unit=log_unit)

  !-- solve linear system --
  CALL syncwrite('> solve_KSP_SC()',log_unit)
  CALL solve_KSP_SC(); CALL syncwrite_time(io_unit=log_unit)

  !-- MATLAB output of solution and pressures
  CALL syncwrite('> update_solution() and extract_sol_KSP()',log_unit)
  CALL update_solution(); CALL extract_sol_KSP()
  CALL syncwrite_time(io_unit=log_unit)

  !-- MATLAB output --
  CALL syncwrite('> MATLAB output',log_unit)
  CALL write_centers_MATLAB(); CALL write_prml_volumes_MATLAB()
  CALL write_solution_MATLAB(); CALL write_pressure_MATLAB()
  CALL syncwrite_time(io_unit=log_unit)

  !-- write solution (ohm) --
  CALL syncwrite('> Whitney interpolcation()',log_unit)
  CALL calc_whitney_C2_BC(); CALL syncwrite_time(io_unit=log_unit)
  CALL syncwrite('> write_solution()',log_unit)
  CALL write_solution_D0S2('charge','current')
  CALL syncwrite_time(io_unit=log_unit)

  !-- clean up allocated variables and end MPI --
  CALL syncwrite('ParaGEMS: cleaning up Ohm flow solution',log_unit)
  CALL deallocate_lcl_complex(); CALL finalise_diffusion(); CALL end_petsc_mpi()

END PROGRAM
! ohm
!===============================================================================
!
