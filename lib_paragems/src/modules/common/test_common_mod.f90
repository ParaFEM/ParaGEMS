!
!===============================================================================
!-- Common Module Tests
!> Tests for common_mod module
!===============================================================================
!/****/h* modules|common/test_common_mod
!* SYNOPSIS
MODULE test_common_mod
!* PURPOSE
!*   Tests for common_mod module
!* INCLUDES
!*   common_mod
!* CONTAINS
!*   Subroutine              Purpose
!*   common_tests(cnt)
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Tests for common_mod module
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE testing_mod
  USE common_mod
  USE mpi_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/common/test_common_mod|common_test
!* SYNOPSIS
  SUBROUTINE common_tests(cnt)
!* PURPOSE
!*   tests common_mod
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
!> tests common_mod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------------
    INTEGER, INTENT(INOUT)  :: cnt  !> test counter

    !---------------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------------
    CHARACTER(LEN=128)  :: fname='test_output_common.dat' !> test filename
    CHARACTER(LEN=128)  :: str                            !> test string
    INTEGER             :: io_unit = 20                   !> io_unit
    LOGICAL             :: test_log                       !> test logical

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !---------------------------------------------------------------------
    ! Only run checks on root (no MPI to test in Functions and Subroutines)
    !---------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); IF (rank==root) THEN

      !---------------------------------------------------------------------
      ! Error checks
      !---------------------------------------------------------------------
      WRITE(*,*) '------------------------------------------------------------'
      WRITE(*,*) ' TESTING :: Common_mod :'

      !---------------------------------------------------------------------
      !-- precision
      !---------------------------------------------------------------------
      CALL chkerr_unit(cnt,1.d0 - 9.999999999999999d-1 <= 0.d0,&
        'common_tests - accuracy to 16 decimal places')

      !---------------------------------------------------------------------
      !-- chkerr
      !---------------------------------------------------------------------
      OPEN(io_unit,FILE=fname,STATUS='UNKNOWN',ACTION='READWRITE',IOSTAT=ier)
      IF (ier/=0) THEN; WRITE(*,*) '... error opening test file '//&
        TRIM(fname)//' for write'; STOP; END IF
      !-- negative value -- will stop \ _ / can't test here, but chkerr_log is
      !-- positive value -- will stop /   \ a similar function tested below
      !-- zero value
      CALL chkerr(0,'no error',io_unit)
      BACKSPACE(io_unit); READ(io_unit,*,IOSTAT=ier) str
      CALL chkerr_unit(cnt,ier/=-1,'common_tests - chkerr zero val '//&
        '(file output)'); BACKSPACE(io_unit)

      !---------------------------------------------------------------------
      !-- chkerr_log
      !---------------------------------------------------------------------
      !-- negative value
      test_log = chkerr_log(-1,'negative error',io_unit)
      CALL chkerr_unit(cnt,.NOT. test_log,'common_tests - chkerr_log - val')
      BACKSPACE(io_unit); READ(io_unit,*) str
      CALL chkerr_unit(cnt,TRIM(str)=='Error :: negative error : errcode '//&
        '=    -1','common_tests - chkerr_log - val (file output)')

      !-- zero value
      test_log = chkerr_log(0,'no error',io_unit)
      CALL chkerr_unit(cnt,test_log,'common_tests - chkerr_log zero val')

      !-- positive value
      test_log = chkerr_log(1,'positive error',io_unit)
      CALL chkerr_unit(cnt,.NOT. test_log,'common_tests - chkerr_log + val')
      BACKSPACE(io_unit); READ(io_unit,*) str
      CALL chkerr_unit(cnt,TRIM(str)=='Error :: negative error : errcode '//&
        '=     1','common_tests - chkerr_log + val (file output)')
      CLOSE(io_unit,STATUS='delete')

      !---------------------------------------------------------------------
      ! DO NOTHING
      !------------
      WRITE(*,*) ' NO TESTS FOR :: common_mod :'
      WRITE(*,*) ' - deallocate_lcl_complex'
      WRITE(*,*) ' - chkerr - positive values (will call STOP)'
      WRITE(*,*) ' - chkerr - negative values (will call STOP)'
      WRITE(*,*) ' - global element'
      WRITE(*,*) ' - global io'
      WRITE(*,*) ' - global mesh'
      WRITE(*,*) ' - global mpi'
      WRITE(*,*) ' - global petsc'
      WRITE(*,*) ' - global vars'
      !---------------------------------------------------------------------
    END IF

  END SUBROUTINE
! test_common_mod|common_test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! common_test
!===============================================================================
!
