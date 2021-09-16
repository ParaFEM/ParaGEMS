!
!===============================================================================
!-- Ohm flow - two-field formulation
!> Miniapp to solve two-field Ohm flow in parallel using PETSc KSP
!===============================================================================
!/****p* programs|darcy_flow/ohm_2f
!* SYNOPSIS
PROGRAM ohm_2f_time
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
  USE time_marching_mod
  USE mpi_mod
  USE geometry_mod
  USE partition_mod
  USE darcy_mod
  USE solver_mod

  IMPLICIT NONE

  !-- local variables --
  CHARACTER(len=slen) :: msg

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- start MPI processes and write miniapp details --
  CALL start_petsc_mpi()
  CALL rootwrite_log('ParaGEMS: Solving two field Ohm flow in parallel '//&
    'using PETSc KSP')

  !-- read input --
  time_dependent = .TRUE.
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
  CALL initialise_time_marching(); CALL set_IC_darcy();
  CALL setup_LHS_darcy(); CALL setup_RHS_darcy();  CALL setup_IC_BC_darcy2();

  !-- write initial data to file --
  IF (sol_output) THEN
    CALL write_centers_MATLAB();  CALL write_prml_volumes_MATLAB()
    CALL write_solution_MATLAB(0);   CALL write_pressure_MATLAB(0)
    CALL write_unsteady_D0S2('pressure','velocity',0)
  END IF

  !-- linear multistep startup --
  !TODO: CALL multistep_startup()

  !-- loop through time steps --
  IF (time_marching%s_initial==2) THEN
    CALL get_RHS_darcy();
    CALL VecCopy(r,r_stage(time_marching%s),petsc_ier)
    CALL VecAssemblyBegin(r_stage(time_marching%s),petsc_ier);
    CALL VecAssemblyEnd(r_stage(time_marching%s),petsc_ier)
  END IF

  DO step = initial_step, final_step
    !-- loop through time stages --
    IF (time_marching%s_initial==2) THEN
      CALL VecCopy(r_stage(time_marching%s),r_stage(1),petsc_ier)
      CALL VecAssemblyBegin(r_stage(1),petsc_ier)
      CALL VecAssemblyEnd(r_stage(1),petsc_ier)
    END IF

    DO stage = time_marching%s_initial,time_marching%s
      !-- update simulation time --
      time = time + dt*time_marching%dc(stage)

      WRITE(msg,'(A,I5,A,I5,A,F8.5)') &
        '           >>>>>>>>>>>>>>>>>>>>> step: ',step,&
        ';         stage: ',stage,';         time: ', time
      CALL syncwrite_log(msg)

      !-- extrapolate solution --
      !TODO: CALL extrapolate_time_stage()

      !-- get RHS for darcy and time marching --
      CALL update_LHS_unsteady();  CALL get_RHS_darcy();  CALL get_RHS_DIMRK()

      !-- solve linear system --
      CALL solve_KSP_SC()

      !-- update solution --
      CALL update_solution()
      CALL extract_sol_KSP();  CALL get_RHS_darcy()
      CALL VecCopy(r,r_stage(stage),petsc_ier)
      CALL VecAssemblyBegin(r_stage(stage),petsc_ier)
      CALL VecAssemblyEnd(r_stage(stage),petsc_ier)

    END DO
    !-- end time stages --

    !-- update simulation time and solution --
    time = time + dt*time_marching%dc(stage)
    WRITE(msg,'(A,I5,A,F10.5)') &
      '           >>>>>>>>>>>>>>>>>>>>> step complete: ',step,&
      ';         time: ', time
    IF (rank==root)  WRITE(*,'(A,I5,A,F10.5,I5)') &
      '           >>>>>>>>>>>>>>>>>>>>> step complete: ',step,&
      ';         time: ', time
    CALL update_time_step()
    CALL extract_sol_KSP()

    !-- MATLAB output of solution and pressures
    IF (sol_output .and. MOD(step,output_frqcy)==0) THEN
      !-- write MATLAB unsteady solution (darcy) --
      CALL write_solution_MATLAB(step);  CALL write_pressure_MATLAB(step)

      !-- write vtk unsteady solution (darcy) --
      CALL syncwrite_log('> Whitney interpolcation()')
      CALL calc_whitney_C2_BC()
      CALL syncwrite_log('> write_solution()')
      CALL write_unsteady_D0S2('pressure','velocity',step)
    END IF
  END DO
  !-- end time steps --

  !-- clean up allocated variables and end MPI --
  CALL syncwrite_log('ParaGEMS: cleaning up Ohm flow solution')
  CALL deallocate_lcl_complex();  CALL finalise_darcy();  CALL end_petsc_mpi()

END PROGRAM
! ohm
!===============================================================================
!
