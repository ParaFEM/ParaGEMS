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
  USE testing_mod
  USE common_mod;         USE test_common_mod
  ! USE geometry_mod;       USE forman_mod
  ! USE dec_mod;            USE test_dec_mod
  ! USE io_mod;             USE test_io_mod
  ! USE math_mod;           USE test_math_mod
  USE mpi_mod;            USE partition_mod;      USE test_mpi_mod
  ! USE darcy_mod;          USE test_darcy_mod
  ! USE solver_mod;         USE test_solver_mod


  IMPLICIT NONE

  !---------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------
  INTEGER             :: cnt=0    !> test counter
  CHARACTER(LEN=slen) :: fname    !> file name
  LOGICAL             :: fexists  !> logical for file existence

!===============================================================================
! MAIN EXECUTION
!===============================================================================

  !---------------------------------------------------------------------
  ! Initial MPI checks
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! start_mpi
  !---------------------------------------------------------------------
  CALL start_mpi(); CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  IF (rank==root) THEN
    WRITE(*,*) '------------------------------------------------------------'
    WRITE(*,*) ' TESTING :: some initial MPI routines from mpi_mod :'
    CALL chkerr_unit(cnt,.FALSE.,'mpi_tests - '//&
      'start_mpi (did not fail)')
  END IF

  !---------------------------------------------------------------------
  ! open_log
  !---------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  IF (rank == root) THEN !-- delete old log files first
    fname = TRIM(log_prefix)//'.log'; INQUIRE(FILE=fname, EXIST=fexists)
    IF (fexists) THEN; OPEN(log_unit,FILE=fname,STATUS="OLD")
      CLOSE(log_unit,STATUS="DELETE"); END IF
    fname = TRIM(log_prefix) // '.log.1'; INQUIRE(FILE=fname, EXIST=fexists)
    IF (fexists) THEN; OPEN(log_unit,FILE=fname,STATUS="OLD")
      CLOSE(log_unit,STATUS="DELETE"); END IF
  END IF
  CALL open_log(); CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  IF (rank==root) CALL chkerr_unit(cnt,.FALSE.,'mpi_tests - '//&
      'open_log (did not fail)')

  !---------------------------------------------------------------------
  ! open_unsteady_log
  !---------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL open_unsteady_log(); CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  IF (rank==root) CALL chkerr_unit(cnt,.FALSE.,'mpi_tests - '//&
    'open_unsteady_log (did not fail)')

  !---------------------------------------------------------------------
  CALL common_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL geometry_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL dec_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL forman_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL io_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL math_tests(cnt)
  !---------------------------------------------------------------------
  CALL mpi_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL partition_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL darcy_tests(cnt)
  !---------------------------------------------------------------------
  ! CALL solver_tests(cnt)
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! final MPI checks
  !---------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  IF (rank==root) THEN
    WRITE(*,*) '------------------------------------------------------------'
    WRITE(*,*) ' TESTING :: some final MPI routines from mpi_mod :'
  END IF

  !---------------------------------------------------------------------
  ! close unsteady log file (already tested)
  !---------------------------------------------------------------------
  CALL close_unsteady_log()

  !---------------------------------------------------------------------
  ! end_mpi
  !---------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  IF (rank==root) THEN
    CALL end_mpi(); CALL chkerr_unit(cnt,.FALSE.,'mpi_tests - '//&
      'end_mpi (did not fail)')
  ELSE; CALL end_mpi(); ENDIF

END PROGRAM
! paragems_unit
!===============================================================================
!
