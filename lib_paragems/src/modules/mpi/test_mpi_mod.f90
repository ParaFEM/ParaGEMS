!
!===============================================================================
!-- mpi Module Tests
!> Tests for mpi_mod module
!===============================================================================
!/****/h* modules|mpi/test_mpi_mod
!* SYNOPSIS
MODULE test_mpi_mod
!* PURPOSE
!*   Tests for mpi_mod module
!* INCLUDES
!*   mpi_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   mpi_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for mpi_mod module
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE testing_mod
  USE mpi_mod
  USE partition_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/mpi/test_mpi_mod|mpi_test
!* SYNOPSIS
  SUBROUTINE mpi_tests(cnt)
!* PURPOSE
!*   tests mpi_mod
!* INPUTS
!*   test counter
!* OUTPUTS
!*   test counter
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> tests mpi_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- use statements and implicit none --
    USE common_mod

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT) :: cnt             !> test counter

    !---------------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------------
    CHARACTER(LEN=128)  :: fname='test_output_mpi.dat'  !> test filename
    CHARACTER(LEN=128)  :: str,str2                     !> test strings
    INTEGER             :: io_unit = 20                 !> io_unit
    INTEGER             :: i                            !> counter
    LOGICAL             :: test_log                     !> test logical

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !---------------------------------------------------------------------
    ! Error checks
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN
      WRITE(*,*) '------------------------------------------------------------'
      WRITE(*,*) ' TESTING :: mpi_mod :'; END IF

    !---------------------------------------------------------------------
    ! chkerrMPI
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN
      OPEN(io_unit,FILE=fname,STATUS='UNKNOWN',ACTION='READWRITE',IOSTAT=ier)
      IF (ier/=0) THEN; WRITE(*,*) '... error opening test file '//&
        TRIM(fname)//' for write'; STOP; END IF; END IF
    !-- negative value -- will call MPI_FINALIZE
    !-- positive value -- will call MPI_FINALIZE
    !-- zero value
    CALL chkerrMPI(0,'no error',io_unit)
    IF (rank==root) THEN; BACKSPACE(io_unit); READ(io_unit,*,IOSTAT=ier) str
      CALL chkerr_unit(cnt,ier/=-1,'mpi_tests - chkerrMPI zero val '//&
      '(file output)'); CLOSE(io_unit,STATUS='delete'); END IF

    !---------------------------------------------------------------------
    ! chkerrMPI_BCAST
    !---------------------------------------------------------------------
    !-- TRUE - will call MPI_FINALIZE
    !-- FALSE
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN
      test_log = .FALSE.; ELSE; test_log = .TRUE.; END IF
    CALL chkerrMPI_BCAST(test_log); CALL chkerr_unit_MPI(cnt,.FALSE.,&
      'mpi_tests - chkerrMPI_BCAST FALSE (did not exit)')

    !---------------------------------------------------------------------
    ! elm2proc
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN
      !-- negative value - no built in checks for this
      WRITE(*,*) 'WARNING :: mpi_tests - elm2proc values less than or equal to zero'
      WRITE(*,*) '                     * no built in checks, but this should not occur'

      WRITE(*,*) 'WARNING :: mpi_tests - elm2proc num_pelm_pp(1) == [0,1]'
      WRITE(*,*) '                     * no built in checks, but this should not occur'
      WRITE(*,*) '                       # elements must be > # processes'

      !-- no extra_pelm
      ALLOCATE(num_pelm_pp(1)); extra_pelm = 0; num_pelm_pp(1) = 10
      !-- proc 1
      CALL chkerr_unit(cnt,elm2proc(1)/=0, 'mpi_tests - elm2proc proc 1 LB - extra_pelm=0');
      CALL chkerr_unit(cnt,elm2proc(10)/=0,'mpi_tests - elm2proc proc 1 UB - extra_pelm=0');
      !-- proc 2
      CALL chkerr_unit(cnt,elm2proc(11)/=1,'mpi_tests - elm2proc proc 2 LB - extra_pelm=0');
      CALL chkerr_unit(cnt,elm2proc(20)/=1,'mpi_tests - elm2proc proc 2 UB - extra_pelm=0');
      !-- proc n
      CALL chkerr_unit(cnt,elm2proc(21)/=2,'mpi_tests - elm2proc proc n LB - extra_pelm=0');

      !-- no extra_pelm
      extra_pelm = 2; num_pelm_pp(1) = 11
      !-- proc 1
      CALL chkerr_unit(cnt,elm2proc(1)/=0, 'mpi_tests - elm2proc proc 1 LB - extra_pelm=2');
      CALL chkerr_unit(cnt,elm2proc(11)/=0,'mpi_tests - elm2proc proc 1 UB - extra_pelm=2');
      !-- proc 2
      CALL chkerr_unit(cnt,elm2proc(12)/=1,'mpi_tests - elm2proc proc 2 LB - extra_pelm=2');
      CALL chkerr_unit(cnt,elm2proc(22)/=1,'mpi_tests - elm2proc proc 2 UB - extra_pelm=2');
      !-- proc 2
      CALL chkerr_unit(cnt,elm2proc(23)/=2,'mpi_tests - elm2proc proc 3 LB - extra_pelm=2');
      CALL chkerr_unit(cnt,elm2proc(32)/=2,'mpi_tests - elm2proc proc 3 UB - extra_pelm=2');
      !-- proc n
      CALL chkerr_unit(cnt,elm2proc(33)/=3,'mpi_tests - elm2proc proc n LB - extra_pelm=2');

      !-- clean up
      DEALLOCATE(num_pelm_pp); END IF

    !---------------------------------------------------------------------
    ! close_log
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    CALL close_log(); CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    CALL chkerr_unit_MPI(cnt,.FALSE.,'mpi_tests - close_log (did not fail)')

    !---------------------------------------------------------------------
    ! close_unsteady_log
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    CALL close_unsteady_log(); CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    CALL chkerr_unit_MPI(cnt,.FALSE.,'mpi_tests - close_unsteady_log '//&
        '(did not fail)')

    !---------------------------------------------------------------------
    ! open_log
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); CALL open_log(); IF (rank==root) THEN
      !-- test log naming
      INQUIRE(FILE=TRIM(log_prefix)//'.log.1',EXIST=test_log);
      CALL chkerr_unit(cnt,.NOT.test_log,'mpi_tests - open_log naming')

      !-- check open_log header
      OPEN(io_unit,FILE=TRIM(log_prefix) // '.log',STATUS='OLD',&
        ACTION='READ',IOSTAT=ier)
      READ(io_unit,'(A)') str; test_log = (TRIM(str)/=&
        '================================================')
      READ(io_unit,'(A)') str; test_log = test_log .OR. (TRIM(str)/=&
        'ParaGEMS: Parallel GEometric MechanicS')
      READ(io_unit,'(A)') str; test_log = test_log .OR. (TRIM(str)/=&
        '------------------------------------------------')
      READ(io_unit,'(A)') str; test_log = test_log .OR. (TRIM(str)/=&
        'Developed by: Pieter Boom')
      READ(io_unit,'(A)') str; test_log = test_log .OR. (TRIM(str)/=&
        'University of Manchester')
      READ(io_unit,'(A)') str; test_log = test_log .OR. (TRIM(str)/=&
        '================================================')
      WRITE(str2,'(A,I7,A)') 'ParaGEMS: running: ',num_procs,' processes'
      READ(io_unit,'(A,I7,A)') str; test_log = test_log .OR. &
        (TRIM(str)/=TRIM(str2)); CLOSE(io_unit)
      CALL chkerr_unit(cnt,test_log,'mpi_tests - open_log header')
    END IF

    !---------------------------------------------------------------------
    ! open_unsteady_log
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN;
      !-- check open_unsteady_log header
      OPEN(io_unit,FILE=TRIM(unstdy_prefix)//'.log',STATUS='OLD',&
        ACTION='READ',IOSTAT=ier)
      WRITE(str2,'(A,I7,A)') 'ParaGEMS: running: ', num_procs ,' processes'
      READ(io_unit,'(A,I7,A)') str; CALL chkerr_unit(cnt,TRIM(str)/=&
        'iter  max_flx indx   area   area^(3/2)   csum_area   csum_area^(3/2)',&
        'mpi_tests - open_unsteady_log header'); CLOSE(io_unit)
    END IF; CALL MPI_BARRIER(MPI_COMM_WORLD,ier); CALL open_unsteady_log()

    !---------------------------------------------------------------------
    ! test log file write functions
    !---------------------------------------------------------------------
    !-- preparing (tested further down)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN
      OPEN(io_unit,FILE=fname,STATUS='UNKNOWN',ACTION='READWRITE',IOSTAT=ier)
      IF (ier/=0) THEN; WRITE(*,*) '... error opening test file '//&
        TRIM(fname)//' for write'; STOP; END IF; END IF

    !-- rootwrite
    CALL rootwrite('test_mpi_mod - rootwrite',io_unit); IF (rank==root) THEN;
      BACKSPACE(io_unit); READ(io_unit,'(A)',IOSTAT=ier) str
      CALL chkerr_unit(cnt,TRIM(str)/='test_mpi_mod - rootwrite',&
      'mpi_tests - rootwrite'); END IF

    !-- syncwrite
    CALL syncwrite('test_mpi_mod - syncwrite',io_unit); IF (rank==root) THEN;
      BACKSPACE(io_unit); READ(io_unit,'(A)',IOSTAT=ier) str
      CALL chkerr_unit(cnt,TRIM(str)/='test_mpi_mod - syncwrite',&
      'mpi_tests - syncwrite'); END IF

    !-- syncwrite_time
    CALL syncwrite_time('test_mpi_mod - syncwrite_time',io_unit)
    IF (rank==root) THEN
      BACKSPACE(io_unit); BACKSPACE(io_unit); BACKSPACE(io_unit);
      READ(io_unit,'(A)') str; test_log = (TRIM(str)==&
        'test_mpi_mod - syncwrite_time')
      READ(io_unit,'(A)') str; test_log = test_log .AND. &
        (str(1:22)=='   - TIME for process:')
      READ(io_unit,'(A)') str; test_log = test_log .AND. &
        (str(1:25)=='   - TIME from MPI start:')
      CALL chkerr_unit(cnt,.NOT. test_log,'mpi_tests - syncwrite_time')
      CLOSE(io_unit,STATUS='delete'); END IF

    !---------------------------------------------------------------------
    ! DO NOTHING ???
    !------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN
      WRITE(*,*) ' NO TESTS FOR :: mpi_mod :'
      WRITE(*,*) ' - chkerrMPI - positive values (will call MPI_FINALIZE)'
      WRITE(*,*) ' - chkerrMPI - negative values (will call MPI_FINALIZE)'
      WRITE(*,*) ' - chkerrMPI_BCAST - TRUE (will call MPI_FINALIZE)'; END IF

  END SUBROUTINE
! test_mpi_mod|mpi_tests
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/mpi/test_mpi_mod|partition_test
!* SYNOPSIS
  SUBROUTINE partition_tests(cnt)
!* PURPOSE
!*   tests partition_mod
!* INPUTS
!*   test counter
!* OUTPUTS
!*   test counter
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> tests partition_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT) :: cnt             !> test counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------------------------------------------------------------
    !*   parallel_setup          sets up partitioning, read mesh and nodal
    !*                           locations, and sets connectivity
    !*   partition_dual          partitions dual volumes across available processes
    !*   get_lcl_prml_elms       ??? get local node indices of highest order primal elements
    !*   get_connectivity
    !*   get_glb_indx_dual_vlm
    !*   get_glb_indx
    !*   exchange_prml_nodes
    !*   exchange_glb_indx
    !---------------------------------------------------------------------

  END SUBROUTINE
! test_mpi_mod|partition_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! mpi_test
!===============================================================================
!
