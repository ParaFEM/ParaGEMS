!
!===============================================================================
!-- MPI Module
!> Module contains routines for controlling MPI execution environment
!===============================================================================
!/****/h* modules|mpi/mpi_mod
!* SYNOPSIS
MODULE mpi_mod
!* PURPOSE
!*   Module contains routines for controlling MPI execution environment
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!* CONTAINS
!*   Subroutine              Purpose
!*   start_mpi()             start MPI execution environment & get details
!*   end_mpi()               end MPI execution environment & stop ParaGEMS
!*   elm2proc()              convert partitioning element to process num (rank)
!*   open_log                Open log file from root process & write header
!*   close_log               Close log file from root process
!*   rootwrite_log           Write unsynchronised message to log from root
!*   write_log               Write unsynchronised messages to log from all ranks
!*   syncwrite_log           Write synchronised message to log from root
!*   syncwrite_log_mpidata   Write a synchronised message to log with num ranks
!*   syncwrite_log_time      Write a synchronised message to log with timings
!*   open_unsteady_log       Open log file for unsteady simulations from root
!*                           and write header
!*   close_unsteady_log      Close log file for unsteady simulations from root
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Module contains routines for controlling MPI execution environment
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/start_mpi
!* SYNOPSIS
  SUBROUTINE start_mpi()
!* PURPOSE
!*   Start MPI execution environment & get details
!* SIDE EFFECTS
!*   - MPI execution environment started, rank identified, and number of
!*     processes evaluated
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Start MPI execution environment & get details
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- start MPI --
    CALL MPI_INIT(ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_INIT: errcode = ',ier
      STOP 'ParaGEMS: shutdown'
    END IF

    !-- get local rank --
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_RANK: errcode = ',ier
      CALL end_mpi()
    END IF

    !-- get the global number of processes --
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_SIZE: errcode = ',ier
      CALL end_mpi()
    END IF

    !-- record simumation start time --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    start_time = MPI_Wtime();   curnt_time = start_time

    !-- open paragems log file and write mpi data --
    CALL open_log();   CALL syncwrite_log_mpidata()

    RETURN

  END SUBROUTINE
! mpi_mod/start_mpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/end_mpi
!* SYNOPSIS
  SUBROUTINE end_mpi()
!* PURPOSE
!*   End MPI execution environment & stop ParaGEMS
!* SIDE EFFECTS
!*   - MPI execution environment ended and ParaGEMS stoped
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> End MPI execution environment & stop ParaGEMS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen) :: str_out          !> STOP output string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- end MPI execution environment and stop ParaGEMS --
    WRITE(str_out,'(A,A)') &
      'ParaGEMS: shutdown: the program terminated successfully'
    CALL syncwrite_log(str_out);   CALL syncwrite_log_time();   CALL close_log()
    CALL MPI_FINALIZE(ier)
    IF (ier /= MPI_SUCCESS) &
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_SIZE: errcode = ',ier
    STOP

    RETURN

  END SUBROUTINE
! mpi_mod/end_mpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f* mpi_mod/elm2proc
!* SYNOPSIS
  FUNCTION elm2proc(id)
!* PURPOSE
!*   Convert partitioning element to process number (rank)
!* INPUTS
!*   Name                    Description
!*   id                      global vertex index
!* OUTPUTS
!*   Name                    Description
!*   id                      global vertex index
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Convert partitioning element to process number (rank)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)   :: id         !>
    INTEGER               :: elm2proc   !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- identify process associated with given primary index --
    IF (extra_pelm==0 .OR. id<extra_pelm*num_pelm_pp(1)) THEN
      elm2proc = (id-1)/num_pelm_pp(1)
    ELSE
      elm2proc = (id-extra_pelm-1)/(num_pelm_pp(1)-1)
    END IF

    RETURN

  END FUNCTION
! mpi_mod/elm2proc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/open_log
!* SYNOPSIS
  SUBROUTINE open_log()
!* PURPOSE
!*   Open log file from root process and write ParaGEMS header
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Open log file from root process and write ParaGEMS header
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname          !> file name
    LOGICAL               :: err = .FALSE.  !> error variable
    integer               :: i=0            !> counter
    LOGICAL               :: fexists        !> logical for file existence

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process --
    IF (rank == root) THEN
      !-- Open log file and check for errors --
      fname = trim(log_prefix) // '.log'
      INQUIRE(FILE=fname, EXIST=fexists)
      IF (fexists) THEN
        i = 1
        DO WHILE(fexists)
          IF (i<10) THEN
            WRITE(fname,'(A,I1)') trim(log_prefix) // '.log.',i
          ELSEIF (i<100) THEN
            WRITE(fname,'(A,I2)') trim(log_prefix) // '.log.',i
          ELSEIF (i<100) THEN
            WRITE(fname,'(A,I3)') trim(log_prefix) // '.log.',i
          END IF
          INQUIRE(FILE=fname, EXIST=fexists)
          i = i+1
        END DO
      END IF
      OPEN(log_unit,FILE=fname,STATUS='NEW',ACTION='WRITE',IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
        err = .TRUE.
      END IF

      !-- write header --
      WRITE(log_unit,*) '======================================================'
      WRITE(log_unit,*) 'ParaGEMS: Parallel GEometric Mechanics of Solids'
      WRITE(log_unit,*) '------------------------------------------------------'
      WRITE(log_unit,*) 'Developed by: Pieter Boom'
      WRITE(log_unit,*) 'University of Manchester'
      WRITE(log_unit,*) '======================================================'
      CALL FLUSH(log_unit)
    END IF

    !-- check for errors
    CALL MPI_BCAST(err,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    IF (err) CALL end_mpi()

    RETURN

  END SUBROUTINE
! mpi_mod/open_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/close_log
!* SYNOPSIS
  SUBROUTINE close_log()
!* PURPOSE
!*   Close log file from root process
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Close log file from root process
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process: close log file --
    IF (rank == root) CLOSE(log_unit)

    RETURN

  END SUBROUTINE
! mpi_mod/close_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/rootwrite_log
!* SYNOPSIS
  SUBROUTINE rootwrite_log(msg)
!* PURPOSE
!*   Write unsynchronised message from root process to log file
!* INPUTS
!*   Name                    Description
!*   msg                     message string to be written to file
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Write unsynchronised message from root process to log file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)   :: msg   !> message string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- write message to log file --
    IF (rank==root) THEN
      WRITE(log_unit,*) msg;   CALL FLUSH(log_unit)
    END IF

    RETURN

  END SUBROUTINE
! mpi_mod/rootwrite_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/write_log
!* SYNOPSIS
  SUBROUTINE write_log(msg)
!* PURPOSE
!*   Write unsynchronised messages collated from processes to log file
!* INPUTS
!*   Name                    Description
!*   msg                     message string to be written to file
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Write unsynchronised messages collated from processes to log file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=slen), INTENT(IN)   :: msg          !> message string

    !-- local variables --
    INTEGER             :: i                          !> counter
    INTEGER             :: status(MPI_STATUS_SIZE)    !> mpi status variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IF (rank == root) THEN
      !-- write message to log file --
      WRITE(log_unit,'(A,I8,A,A)') 'rank: ',rank,' : ',msg
      DO i = 1,num_procs-1
        CALL MPI_RECV(msg,slen,MPI_CHARACTER,i,i,MPI_COMM_WORLD,status,ier)
        WRITE(log_unit,'(A,I8,A,A)') 'rank: ',i,' : ',msg
      END DO
      CALL FLUSH(log_unit)
    ELSE
      CALL MPI_SEND(msg,slen,MPI_CHARACTER,0,rank,MPI_COMM_WORLD,ier)
    END IF

    RETURN

  END SUBROUTINE
! mpi_mod/write_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/syncwrite_log
!* SYNOPSIS
  SUBROUTINE syncwrite_log(msg)
!* PURPOSE
!*   Write synchronised message from root process to log file
!* INPUTS
!*   Name                    Description
!*   msg                     message string to be written to file
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Write uynchronised message from root process to log file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)   :: msg          !> message string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- synchronise processes --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);

    !-- write message from root process --
    IF (rank == root) THEN
      WRITE(log_unit,*) msg;   CALL FLUSH(log_unit)
    END IF

    RETURN

  END SUBROUTINE
! mpi_mod/syncwrite_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/syncwrite_log_mpidata
!* SYNOPSIS
  SUBROUTINE syncwrite_log_mpidata()
!* PURPOSE
!*   Write a synchronised message with number of mpi ranks to log file
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Write a synchronised message with number of mpi ranks to log file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- synchronise processes --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);

    !-- write message from root process --
    IF (rank == root) THEN
      WRITE(log_unit,*) 'ParaGEMS: running: ', num_procs ,' processes'
      CALL FLUSH(log_unit)
    END IF

  END SUBROUTINE
! mpi_mod/syncwrite_log_mpidata
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/syncwrite_log_time
!* SYNOPSIS
  SUBROUTINE syncwrite_log_time()
!* PURPOSE
!*   Write a synchronised message to log file with timings
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Write a synchronised message to log file with timings
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)            :: log_msg       !> message string
    REAL(KIND=iwp)                 :: time          !> time variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- synchronise processes and get current time --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    time = MPI_WTIME()

    !-- write msg with difference in times --
    WRITE(log_msg,'(A,E10.3)') '   - TIME for process:     ', time-curnt_time
    CALL rootwrite_log(log_msg)
    WRITE(log_msg,'(A,E10.3)') '   - TIME from MPI start:  ', time-start_time
    CALL rootwrite_log(log_msg)

    !-- update current time --
    curnt_time = time

    RETURN

  END SUBROUTINE
! mpi_mod/syncwrite_log_time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/open_unsteady_log
!* SYNOPSIS
  SUBROUTINE open_unsteady_log()
!* PURPOSE
!*   Open log file for unsteady simulations from root process and write header
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Open log file for unsteady simulations from root process and write header
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname          !> file name
    LOGICAL               :: err = .FALSE.  !> error variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process --
    IF (rank == root) THEN
      !-- Open log file and check for errors --
      fname = trim(unstdy_prefix) // '.log'
      OPEN(ulog_unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
        err = .TRUE.
      END IF
      WRITE(ulog_unit,*) 'iter  max_flx indx   area   area^(3/2)   csum_area   csum_area^(3/2)'
      CALL FLUSH(ulog_unit)
    END IF

    !-- check for errors
    CALL MPI_BCAST(err,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    IF (err) CALL end_mpi()

    RETURN

  END SUBROUTINE
! mpi_mod/open_unsteady_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/close_unsteady_log
!* SYNOPSIS
  SUBROUTINE close_unsteady_log()
!* PURPOSE
!*   Close log file for unsteady simulations from root
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Close log file for unsteady simulations from root
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process: close log file --
    IF (rank == root) CLOSE(ulog_unit)

    RETURN

  END SUBROUTINE
! mpi_mod/close_unsteady_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! mpi_mod
!===============================================================================
!
