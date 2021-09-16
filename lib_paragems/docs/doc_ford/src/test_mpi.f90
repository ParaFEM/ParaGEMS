!
!===============================================================================
! /****m*/ test_mpi
MODULE test_mpi
!
! PURPOSE:    Module contains tests for mpi_mod module
!
! CONTAINS:   Subroutine            Purpose
!             t_start_mpi()         - start MPI execution environment & get details
!             t_end_mpi()           - end MPI execution environment & stop ParaGEMS
!             t_parallel_setup()    - sets up partitioning, read mesh and nodal
!                                   locations, and sets connectivity
!             t_partition()         - partitions work across available processes
!             t_get_connectivity()  - gets process connectivity
!
! UPDATES:  created (PDB) :: 2019/08/21
!
!===============================================================================

  !-----------------------------------------------------------------------------
  ! use statements and implicit none
  !-----------------------------------------------------------------------------
  USE common_mod
  USE mpi_mod
  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/mpi/modules/mpi_test|t_start_mpi
  SUBROUTINE t_start_mpi()
!
! PURPOSE:  Start MPI execution environment & get details
!
! TESTS:    Number  Purpose
!
! UPDATES:  created (PDB) :: 2019/08/21
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------------------------------------------------------------
    !  1. start MPI
    !---------------------------------------------------------------------
    CALL MPI_INIT(ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_INIT: errcode = ',ier
      STOP "ParaGEMS: shutdown"
    END IF

    !---------------------------------------------------------------------
    !  2. get local rank
    !---------------------------------------------------------------------
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_RANK: errcode = ',ier
      CALL end_mpi()
    END IF

    !---------------------------------------------------------------------
    !  3. get the global number of processes
    !---------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_SIZE: errcode = ',ier
      CALL end_mpi()
    END IF

    RETURN

  END SUBROUTINE
! mpi_test|t_start_mpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/modules/mpi/mpi_test|t_end_mpi
  SUBROUTINE t_end_mpi()
!
! PURPOSE:  End MPI execution environment & stop ParaGEMS
!
! TESTS:
!
! UPDATES:  created (PDB) :: 2019/08/21
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------------------------------------------------------------
    !  1. end MPI execution environment and stop ParaGEMS
    !---------------------------------------------------------------------
    CALL MPI_FINALIZE(ier)
    STOP "ParaGEMS: shutdown: the program terminated successfully"

    RETURN

  END SUBROUTINE
! mpi_test|t_end_mpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/modules/common/common_test|t_parallel_setup
  SUBROUTINE t_parallel_setup()
!
! PURPOSE:  sets up partitioning, read mesh and nodal locations, and sets
!           connectivity
!
! TESTS:
!
! UPDATES:  created (PDB) :: 2019/08/21
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    RETURN

  END SUBROUTINE
! common_test|t_parallel_setup
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! /****s*/ /src/modules/common/common_test|t_partition
  SUBROUTINE t_partition()
!
! PURPOSE:
!
! TESTS:
!
! UPDATES:  created (PDB) :: 2019/08/21
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !---------------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------------
    CHARACTER(LEN=slen) :: fname          ! file name
    INTEGER             :: extras         ! modulus of the number of volumes
                                          ! and available processes
    INTEGER             :: int_buffer(2)  ! integer buffer for MPI comms
    INTEGER             :: junk           ! junk IO variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    RETURN

  END SUBROUTINE
! common_test|t_partition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! mpi_test
!===============================================================================
!
