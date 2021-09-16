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
  USE geometry_mod
  USE partition_mod
  USE darcy_mod
  USE solver_mod

  IMPLICIT NONE

  INTEGER :: crck_id, i
  LOGICAL :: exit_cond
  CHARACTER(LEN=slen) :: msg,fname

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- start MPI processes and write miniapp details --
  CALL start_petsc_mpi()
  CALL rootwrite_log('ParaGEMS: Solving two field Darcy flow in parallel '//&
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
  CALL initialise_darcy();  CALL set_IC_darcy();  CALL setup_LHS_darcy()
  CALL setup_RHS_darcy();  CALL setup_IC_BC_darcy();  CALL get_RHS_darcy()

  CALL write_mat_MATLAB(Asub(1),'A1','solution_A1.m')
  CALL write_mat_MATLAB(Asub(2),'A2','solution_A2.m')
  CALL write_mat_MATLAB(Asub(3),'A3','solution_A3.m')
  CALL write_mat_MATLAB(Asub(4),'A4','solution_A4.m')
  CALL write_vec_MATLAB(r,'r','solution_r.m')
  CALL write_vec_MATLAB(b,'b','solution_b.m')
  CALL write_vec_MATLAB(q,'qp','solution_qp.m')

  !-- solve linear system --
  CALL solve_KSP_SC();  CALL update_solution();  CALL extract_sol_KSP()

  CALL write_vec_MATLAB(dq,'dq','solution_dq.m')
  CALL write_vec_MATLAB(q,'q','solution_q.m')

  !-- MATLAB output of solution and pressures
  IF (sol_output) THEN
    !-- write MATLAB unsteady solution (darcy) --
    CALL write_solution_MATLAB(0);   CALL write_pressure_MATLAB(0)

    !-- MATLAB output of centers and volumes
    CALL write_centers_MATLAB();   CALL write_prml_volumes_MATLAB()

    !-- write vtk solution (darcy) --
    CALL syncwrite_log('> Whitney interpolcation()');  CALL calc_whitney_C2_BC()
    CALL syncwrite_log('> write_solution()')
    CALL write_unsteady_D0S2('pressure','velocity',0)
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
      DO i = 1, crcks_pstep
        CALL identify_random_crack(crck_id,exit_cond)
      END DO
    ELSEIF (crck_type==2) THEN
      CALL identify_threshold_crack(crck_id,exit_cond)
    ELSEIF (crck_type==3) THEN
      CALL identify_maxFD_crack(crck_id,exit_cond)
    ELSEIF (crck_type==4) THEN
      CALL identify_minFD_crack(crck_id,exit_cond)
    END IF
    IF (exit_cond) EXIT

    !-- solve linear system --
    CALL get_RHS_darcy()
    CALL write_mat_MATLAB(Asub(1),'A12','solution_A12.m')
    CALL write_mat_MATLAB(Asub(2),'A22','solution_A22.m')
    CALL write_mat_MATLAB(Asub(3),'A32','solution_A32.m')
    CALL write_mat_MATLAB(Asub(4),'A42','solution_A42.m')
    CALL write_vec_MATLAB(q,'q1','solution_q.m')
    CALL write_vec_MATLAB(r,'r2','solution_r2.m')
    CALL solve_KSP_SC();  CALL update_solution();  CALL extract_sol_KSP()

    CALL write_vec_MATLAB(q,'q2','solution_q2.m')


    !-- MATLAB output of solution and pressures
    IF (sol_output .and. MOD(crck_id,output_frqcy)==0) THEN
      !-- write MATLAB unsteady solution (darcy) --
      CALL write_solution_MATLAB(crck_id)
      CALL write_pressure_MATLAB(crck_id)

      !-- write vtk unsteady solution (darcy) --
      CALL syncwrite_log('> Whitney interpolcation()')
      CALL calc_whitney_C2_BC()
      CALL syncwrite_log('> write_solution()')
      CALL write_unsteady_D0S2('pressure','velocity',crck_id)
    END IF
  END DO

  !-- clean up allocated variables and end MPI --
  CALL syncwrite_log('ParaGEMS: cleaning up Darcy flow solution')
  CALL close_unsteady_log()
  CALL deallocate_lcl_complex();  CALL finalise_darcy();   CALL end_petsc_mpi()

END PROGRAM
! darcy
!===============================================================================
!
