!
!===============================================================================
!-- Darcy flow - two-field formulation with cracking
!> Miniapp to solve two-field Darcy flow in parallel using PETSc KSP with cracking
!===============================================================================
!/****p* programs|darcy_flow/darcy_crkp_2f
!* SYNOPSIS
PROGRAM darcy_crkp_2f
!* PURPOSE
!*   Miniapp to solve two-field Darcy flow in parallel using PETSc KSP with cracking
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
!> Miniapp to solve two-field Darcy flow in parallel using PETSc KSP with cracking
!===============================================================================

  !-- modules and implicit none --
  USE common_mod
  USE mpi_mod
  USE partition_mod
  USE darcy_mod
  USE solver_mod

  IMPLICIT NONE

  INTEGER :: crck_id
  LOGICAL :: exit_cond
  CHARACTER(LEN=slen) :: msg,fname

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- start MPI processes and write miniapp details --
  CALL start_mpi()
  CALL rootwrite_log('ParaGEMS: Solving two field Darcy flow in parallel '// &
    'using PETSc KSP')

  !-- read input --
  CALL read_input_darcy()

  !-- setup partitioning and connectivity --
  CALL parallel_setup()

  !-- initialise geometric information --
  CALL initialise_geo()

  !-- initialise solution --
  CALL start_petsc();   CALL initialise_darcy()

  !-- get necessary matrices and vectors (A x = b) --
  CALL get_LHS_darcy();   CALL get_RHS_darcy()

  !-- solve linear system --
  CALL solve_KSP()

  !-- MATLAB output of solution and pressures
  IF (sol_output) THEN
    !-- extract solution data from PETSc arrays --
    CALL extract_sol_KSP()

    !-- write MATLAB unsteady solution (darcy) --
    CALL write_solution_MATLAB2(0)
    CALL write_pressure_MATLAB2(0)

    !-- MATLAB output of centers and volumes
    CALL write_centers_MATLAB2();   CALL write_prml_volumes_MATLAB2()

    !-- write vtk solution (darcy) --
    CALL syncwrite_log('> Whitney interpolcation()');  CALL calc_whitney_C2_BC()
    CALL syncwrite_log('> write_solution()')
    CALL write_unsteady_D0S('pressure','velocity',0)
  END IF

  !--
  CALL open_unsteady_log()

  !--
  DO crck_id=1,max_crcks
    write(msg,*) 'Crack propagation iteration: ', crck_id
    CALL syncwrite_log(msg)

    !-- identify the faces with maximum velocity and crack --
    IF (crck_type==0) THEN
      write(msg,*) '>>> no cracking';   CALL syncwrite_log(msg);   EXIT
    ELSEIF (crck_type==1) THEN
      !CALL identify_crack5(crck_id,exit_cond)
      CALL identify_crack3(crck_id,exit_cond)
    ELSEIF (crck_type==2) THEN
      !CALL identify_crack4(crck_id,exit_cond)
      CALL identify_crack(exit_cond)
    END IF
    IF (exit_cond) EXIT

    !-- solve linear system --
    CALL solve_KSP()

    !-- MATLAB output of solution and pressures
    IF (sol_output .and. MOD(crck_id,output_frqcy)==0) THEN
      !-- extract solution data from PETSc arrays --
      CALL extract_sol_KSP()

      !-- write MATLAB unsteady solution (darcy) --
      CALL write_solution_MATLAB2(crck_id)
      CALL write_pressure_MATLAB2(crck_id)

      !-- write vtk unsteady solution (darcy) --
      CALL syncwrite_log('> Whitney interpolcation()')
      CALL calc_whitney_C2_BC()
      CALL syncwrite_log('> write_solution()')
      CALL write_unsteady_D0S('pressure','velocity',crck_id)
    END IF
  END DO

  !-- clean up allocated variables and end MPI --
  CALL syncwrite_log('ParaGEMS: cleaning up Darcy flow solution')
  CALL finalise_darcy();   CALL close_unsteady_log()
  CALL clean_up()
  CALL end_petsc();   CALL end_mpi()

END PROGRAM
! darcy
!===============================================================================
!
