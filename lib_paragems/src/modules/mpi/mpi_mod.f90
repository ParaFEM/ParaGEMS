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
!*   start_mpi               start MPI execution environment & get details
!*   end_mpi                 end MPI execution environment & stop ParaGEMS
!*   chkerrMPI               check for errors, print msg and exit gracefully
!*   chkerrMPI_BCAST         complete check for errors and exit gracefully
!*   elm2proc                convert partitioning element to process num (rank)
!*   open_log                open log file from root process & write header
!*   close_log               close log file from root process
!*   rootwrite           write unsynchronised message to log from root
!*   syncwrite           write synchronised message to log from root
!*   syncwrite_time      write a synchronised message to log with timings
!*   open_unsteady_log       open log file for unsteady simulations from root
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
    CALL MPI_INIT(ier); CALL chkerr(ier,'start_mpi > MPI_INIT',6)

    !-- get local rank --
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
    CALL chkerrMPI(ier,'start_mpi > MPI_COMM_RANK',6)

    !-- get the global number of processes --
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ier)
    CALL chkerrMPI(ier,'start_mpi > MPI_COMM_SIZE',6)

    !-- record simumation start time --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    CALL chkerrMPI(ier,'start_mpi > MPI_BARRIER',6) !-- initial MPI check

    !-- open log file and get timing --
    CALL open_log()

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- end MPI execution environment and stop ParaGEMS --
    CALL syncwrite_time('ParaGEMS: shutdown - program terminated '//&
      'successfully',log_unit); CALL close_log()
    CALL MPI_FINALIZE(ier); CALL chkerr(ier,'MPI_FINALIZE',6); RETURN

  END SUBROUTINE
! mpi_mod/end_mpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/chkerrMPI
!* SYNOPSIS
  SUBROUTINE chkerrMPI(ierr,str,io_unit)
!* PURPOSE
!*   check for errors and exit gracefully from MPI environment
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> check for errors and exit gracefully from MPI environment
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)           :: ierr     !> error integer
    INTEGER, INTENT(IN)           :: io_unit  !> io unit
    CHARACTER(LEN=*), INTENT(IN)  :: str      !> error string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- check for errors and exit graacefully if any found --
    IF (ierr /= MPI_SUCCESS) THEN; WRITE(io_unit,'(A,I5)') 'Error :: '//&
      TRIM(str)//' : errcode = ',ierr; CALL end_mpi(); STOP; END IF; RETURN

  END SUBROUTINE
! chkerrMPI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/chkerrMPI
!* SYNOPSIS
  SUBROUTINE chkerrMPI_BCAST(err)
!* PURPOSE
!*   check for errors and exit gracefully from MPI environment
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/06/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/06/16
!> check for errors and exit gracefully from MPI environment
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    LOGICAL, INTENT(IN) :: err  !> error logical

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- check for errors and exit graacefully if any found --
    CALL MPI_BCAST(err,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    IF (err) THEN; CALL end_mpi(); STOP; END IF; RETURN

  END SUBROUTINE
! chkerrMPI_BCAST
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
!*   id                      global partitioning element index
!* OUTPUTS
!*   Name                    Description
!*   elm2proc                associated process
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
    INTEGER, INTENT(IN) :: id         !> global partitioning element index
    INTEGER             :: elm2proc   !> associated process

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- identify process associated with given primary index --
    IF (extra_pelm==0 .OR. id<extra_pelm*num_pelm_pp(1)) THEN
      elm2proc = (id-1)/num_pelm_pp(1)
    ELSE; elm2proc = (id-extra_pelm-1)/(num_pelm_pp(1)-1); END IF; RETURN

  END FUNCTION
! mpi_mod/elm2proc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f* mpi_mod/elm2proc_array
!* SYNOPSIS
  FUNCTION elm2proc_array(id)
!* PURPOSE
!*   Convert partitioning element to process number (rank)
!* INPUTS
!*   Name                    Description
!*   id                      global partitioning element index
!* OUTPUTS
!*   Name                    Description
!*   elm2proc                associated process
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
    INTEGER, INTENT(IN) :: id(:)  !> global partitioning element index
    INTEGER             :: elm2proc_array(size(id)) !> associated process

    !-- local variables --
    INTEGER :: i  !> counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- identify process associated with given primary index --
    DO i=1,size(id)
      IF (extra_pelm==0 .OR. id(i)<extra_pelm*num_pelm_pp(1)) THEN
        elm2proc_array(i) = (id(i)-1)/num_pelm_pp(1)
      ELSE; elm2proc_array(i) = (id(i)-extra_pelm-1)/(num_pelm_pp(1)-1); END IF
    END DO; RETURN

  END FUNCTION
! mpi_mod/elm2proc_array
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/open_log
!* SYNOPSIS
  SUBROUTINE open_log()
!* PURPOSE
!*   Open log file from root process and write ParaGEMS header
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
    CHARACTER(LEN=slen) :: fname          !> file name
    LOGICAL             :: err = .FALSE.  !> error variable
    integer             :: i              !> counter
    LOGICAL             :: fexists        !> logical for file existence

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process --
    IF (rank == root) THEN  !-- Open log file and check for errors --
      fname = TRIM(log_prefix) // '.log'; INQUIRE(FILE=fname, EXIST=fexists)
      IF (fexists) THEN; !-- don't overwrite previous logs
        i=0; DO WHILE(fexists); i = i+1
          IF (i<10) THEN;       WRITE(fname,'(A,I1)')TRIM(log_prefix)//'.log.',i
          ELSEIF (i<100) THEN;  WRITE(fname,'(A,I2)')TRIM(log_prefix)//'.log.',i
          ELSEIF (i<1000) THEN; WRITE(fname,'(A,I3)')TRIM(log_prefix)//'.log.',i
          ELSE; fname = TRIM(log_prefix)//'.log'; EXIT; END IF
          INQUIRE(FILE=fname, EXIST=fexists); END DO; END IF
      OPEN(log_unit,FILE=fname,STATUS='NEW',ACTION='WRITE',IOSTAT=ier)
      err = chkerr_log(ier,'OPEN - error opening '//TRIM(fname),6)
      IF (.NOT. err) THEN  !-- write header --
        WRITE(log_unit,'(A)') '================================================'
        WRITE(log_unit,'(A)') 'ParaGEMS: Parallel GEometric MechanicS'
        WRITE(log_unit,'(A)') '------------------------------------------------'
        WRITE(log_unit,'(A)') 'Developed by: Pieter Boom'
        WRITE(log_unit,'(A)') 'University of Manchester'
        WRITE(log_unit,'(A)') '================================================'
        WRITE(log_unit,'(A,I7,A)') 'ParaGEMS: running: ',num_procs,' processes'
        CALL FLUSH(log_unit); END IF; END IF
    CALL chkerrMPI_BCAST(err)  !-- complete check for errors

    !-- setup timing variables --
    start_time = MPI_Wtime(); curnt_time = start_time; RETURN

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
    IF (rank == root) CLOSE(log_unit); RETURN

  END SUBROUTINE
! mpi_mod/close_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/rootwrite
!* SYNOPSIS
  SUBROUTINE rootwrite(msg,io_unit)
!* PURPOSE
!*   Write unsynchronised message from root process to log file
!* INPUTS
!*   Name                    Description
!*   msg                     message string to be written to file
!*   io_unit                 io unit for file
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
    CHARACTER(LEN=*), INTENT(IN)  :: msg      !> message string
    INTEGER, INTENT(IN)           :: io_unit  !> io unit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- write message to log file --
    IF (rank==root) WRITE(io_unit,'(A)') TRIM(msg); CALL FLUSH(log_unit); RETURN

  END SUBROUTINE
! mpi_mod/rootwrite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/syncwrite
!* SYNOPSIS
  SUBROUTINE syncwrite(msg,io_unit)
!* PURPOSE
!*   Write synchronised message from root process to log file
!* INPUTS
!*   Name                    Description
!*   msg                     message string to be written to file
!*   io_unit                 io unit for file
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
    CHARACTER(LEN=*), INTENT(IN)  :: msg      !> message string
    INTEGER, INTENT(IN)           :: io_unit  !> io unit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- synchronise processes --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    !-- write message from root process --
    IF (rank == root) WRITE(io_unit,'(A)') TRIM(msg); CALL FLUSH(log_unit)
    RETURN

  END SUBROUTINE
! mpi_mod/syncwrite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/syncwrite_time
!* SYNOPSIS
  SUBROUTINE syncwrite_time(msg,io_unit)
!* PURPOSE
!*   Write a synchronised message to log file with timings
!* INPUTS
!*   Name                    Description
!*   msg                     message string to be written to file
!*   io_unit                 io unit for file
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

    !-- arguments --
    INTEGER, INTENT(IN)                   :: io_unit  !> io unit
    CHARACTER(LEN=*),OPTIONAL, INTENT(IN) :: msg      !> message string

    !-- local variables --
    CHARACTER(LEN=slen) :: str    !> message string
    REAL(KIND=PGMSiwp)      :: time   !> time variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- synchronise processes and get current time --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier); time = MPI_WTIME()

    !-- write msg with difference in times --
    IF (rank==root) THEN; IF (PRESENT(msg)) WRITE(io_unit,'(A)') TRIM(msg)
      WRITE(io_unit,'(A,E10.3,A,E10.3)') &
        '   - TIME for process:     ', time-curnt_time, NEW_LINE('A')//&
        '   - TIME from MPI start:  ', time-start_time; CALL FLUSH(log_unit)
    END IF

    !-- update current time --
    curnt_time = time; RETURN

  END SUBROUTINE
! mpi_mod/syncwrite_time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* mpi_mod/open_unsteady_log
!* SYNOPSIS
  SUBROUTINE open_unsteady_log()
!* PURPOSE
!*   Open log file for unsteady simulations from root process and write header
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
    CHARACTER(LEN=SLEN) :: fname          !> file name
    LOGICAL             :: err = .FALSE.  !> error variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process --
    IF (rank == root) THEN
      !-- Open log file and check for errors --
      fname = TRIM(unstdy_prefix) // '.log'
      OPEN(ulog_unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
      err = chkerr_log(ier,'OPEN - opening '//TRIM(fname),log_unit)
      IF (.NOT. err) WRITE(ulog_unit,'(A)') 'iter  max_flx indx   area   '//&
        'area^(3/2)   csum_area   csum_area^(3/2)'; CALL FLUSH(ulog_unit)
    END IF
    CALL chkerrMPI_BCAST(err)  !-- complete check for errors
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
    IF (rank == root) CLOSE(ulog_unit); RETURN

  END SUBROUTINE
! mpi_mod/close_unsteady_log
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! mpi_mod
!===============================================================================
!
