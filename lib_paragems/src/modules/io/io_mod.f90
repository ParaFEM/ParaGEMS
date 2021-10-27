!
!===============================================================================
!-- IO Module
!> Module contains subroutines for reading from and writing to disk
!===============================================================================
!/****/h* modules|io/io_mod
!* SYNOPSIS
MODULE io_mod
!* PURPOSE
!*   Module contains subroutines for reading from and writing to disk
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 ???
!*   dec_mod                 DEC functions
!*   ieee_arithmetic
!*   iso_fortran_env
!8   petsc.h
!* CONTAINS
!*   Subroutine              Purpose
!*   read_prml_elms          read local node indices of highest order primal elements
!*   read_glb_indx
!*   read_nodes_prml         read the nodal locations of the primal mesh
!*   root_open_file_read
!*   root_open_file_write
!*   write_solution_D0S
!*   write_solution_D0S2
!*   write_unsteady_D0S
!*   write_unsteady_D0S2
!*   read_crack_faces
!*   read_bndry_cond
!*   read_bndry_cond2
!*   init_write_MATLAB
!*   clean_PETSc_output
!*   write_centers_MATLAB2
!*   write_prml_volumes_MATLAB2
!*   write_solution_MATLAB2
!*   write_pressure_MATLAB2
!*   write_centers_MATLAB
!*   write_prml_volumes_MATLAB
!*   write_solution_MATLAB
!*   write_pressure_MATLAB
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Module contains subroutines for reading from and writing to disk
!===============================================================================
! TO DO:
! - MPI-IO?
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod
  USE dec_mod
  USE solver_mod
  USE, INTRINSIC :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
  USE, INTRINSIC :: iso_fortran_env, only: real32

  IMPLICIT NONE

  !-- include PETSc variables --
#if FORparaFEM
#else
#include <petsc/finclude/petsc.h>
#endif

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/read_prml_elms
!* SYNOPSIS
  SUBROUTINE read_prml_elms()
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: junk             !> junk IO variable
    INTEGER               :: k                !> simplicial order
    INTEGER               :: orient           !> orientation
    INTEGER               :: i,j,l            !> loop counters
    INTEGER               :: ptr              !> pointer/counter
    LOGICAL               :: err=.FALSE.      !> error variable
    CHARACTER(LEN=slen)   :: msg              !> error message
    INTEGER, ALLOCATABLE  :: int_buffer(:)    !> integer buffer for MPI comms
    INTEGER               :: buffer_size      !> size of MPI comm buffer
    INTEGER               :: iib, fib         !> buffer array indices
    INTEGER               :: num_read_iters   !> number of read iterations
    INTEGER               :: read_size        !> number of reads
    INTEGER               :: resid_read_size  !> reads in final loop
    INTEGER               :: stride           !> stride in buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- Open mesh file and check for errors --
    fname=trim(mesh_prefix)//".ele"; CALL root_open_file_read(fname,mesh_unit)

    IF (rank == root) THEN !-- read header and check for consistency --
      READ(mesh_unit,*) junk, k
      IF (glb_num_elm(k) /= junk) THEN; WRITE(msg,'(A,I5,A,I5)')'mesh file '//&
        TRIM(fname)//' - inconsistency: glb_num_elm(dim_cmplx+1)=',&
        glb_num_elm(k),'/=',junk; err=chkerr_log(-1,msg,log_unit)
      ELSEIF (dim_cmplx /= k - 1) THEN; WRITE(msg,'(A,I5,A,I5)')'mesh file '//&
        TRIM(fname)//' -  inconsistcy: dim_cmplx=',dim_cmplx,'/=',k-1
        err=chkerr_log(-1,msg,log_unit); END IF
    END IF; CALL chkerrMPI_BCAST(err) !-- complete error checking --

    !-- get some basic information for reuse --
    k=dim_cmplx+1; read_size=max_read_size; stride=k; ptr=0; fib=stride - 1
    buffer_size=read_size*stride; num_read_iters=glb_num_elm(k)/read_size
    resid_read_size=mod(glb_num_elm(k),read_size)
    ALLOCATE(int_buffer(read_size*(stride+1)))

    IF (rank == root) THEN !-- root process --
      DO i=0,num_read_iters !-- loop through read iterations --
        IF (i == num_read_iters) THEN !-- values for last read iteration --
          read_size=resid_read_size; buffer_size=read_size*stride; END IF
        iib=1; DO j=1, read_size  !-- read data for current iteration --
          READ(mesh_unit,*) junk, int_buffer(iib:iib+fib); iib=iib+stride
        END DO; int_buffer=int_buffer+indx_offset(1) !-- apply vertex offset --
        !-- broadcast vertex indices --
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        iib=1; DO j=1, read_size !-- count local elements --
          IF (minval(int_buffer(iib:iib+fib))<=num_pelm_pp(1)) ptr=ptr+1
          iib=iib+stride; END DO
      END DO
    ELSE !-- on neighbouring processes --
      DO i=0,num_read_iters !-- loop through read iterations --
        IF (i == num_read_iters) THEN !-- values for last read iteration --
          read_size=resid_read_size; buffer_size=read_size*stride; END IF
        !-- recieve broadcast vertex indices --
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        iib=1; DO j=1, read_size; DO l=0,dim_cmplx !-- count local elements --
            IF (elm2proc(int_buffer(iib+l))==rank) THEN
              ptr=ptr+1; EXIT; END IF !-- count and exit if local --
          END DO; iib=iib+stride; END DO
      END DO; END IF

    !-- allocate and assign primal element data to local structures --
    ALLOCATE(&
      lcl_complex(dim_cmplx+1)%glb_indx(ptr),&
      lcl_complex(dim_cmplx+1)%orientation(ptr),&
      lcl_complex(dim_cmplx+1)%node_indx(ptr,dim_cmplx+1))
    num_elm(dim_cmplx+1)=ptr

    !-- get some basic information for reuse --
    read_size=max_read_size; stride=stride+1; buffer_size=read_size*stride
    ptr=0; fib=fib+1

    IF (rank == root) THEN !-- on root process --
      REWIND(mesh_unit); READ(mesh_unit,*) junk !-- rewind and read header --
      DO i=0,num_read_iters !-- loop through read iterations --
        IF (i == num_read_iters) THEN !-- values for last read iteration --
          read_size=resid_read_size; buffer_size=read_size*stride; END IF
        iib=1; DO j=1, read_size  !-- read data for current iteration --
          READ(mesh_unit,*) int_buffer(iib:iib+fib)
          !-- add vertex offset
          int_buffer(iib+1:iib+fib)=int_buffer(iib+1:iib+fib)+indx_offset(1)
          int_buffer(iib)=i*max_read_size+j; iib=iib+stride; END DO
        !-- broadcast vertex indices --
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        iib=1; DO j=1,read_size !-- loop through data to find local values --
          IF (minval(int_buffer(iib+1:iib+fib))<=num_pelm_pp(1)) THEN; ptr=ptr+1
            !-- assign local values to local arrays and compute orientation --
            lcl_complex(dim_cmplx+1)%glb_indx(ptr)   =int_buffer(iib)
            CALL calc_orientation(int_buffer(iib+1:iib+fib),orient)
            lcl_complex(dim_cmplx+1)%node_indx(ptr,:)=int_buffer(iib+1:iib+fib)
            lcl_complex(dim_cmplx+1)%orientation(ptr)=orient
          END IF; iib=iib+stride; END DO
      END DO; CLOSE(mesh_unit)
    ELSE !-- on neighbouring processes --
      DO i=0,num_read_iters !-- loop through read iterations --
        IF (i == num_read_iters) THEN !-- values for last read iteration --
          read_size=resid_read_size; buffer_size=read_size*stride; END IF
        !-- recieve broadcast vertex indices --
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        iib=1; DO j=1,read_size; DO l=1,dim_cmplx+1 !-- find local values --
            IF (elm2proc(int_buffer(iib+l))==rank) THEN; ptr=ptr+1
              !-- assign local values to local arrays and compute orientation --
              lcl_complex(dim_cmplx+1)%glb_indx(ptr)   =int_buffer(iib)
              CALL calc_orientation(int_buffer(iib+1:iib+fib),orient)
              lcl_complex(dim_cmplx+1)%node_indx(ptr,:)=int_buffer(iib+1:iib+fib)
              lcl_complex(dim_cmplx+1)%orientation(ptr)=orient; EXIT
            END IF; END DO; iib=iib+stride; END DO
      END DO
    END IF

    !-- clean up and make sure all comms have completed --
    DEALLOCATE(int_buffer); CALL MPI_BARRIER(MPI_COMM_WORLD,ier); RETURN

  END SUBROUTINE
! io_mod|read_prml_elms
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/read_nodes_prml
!* SYNOPSIS
  SUBROUTINE read_nodes_prml()
!* PURPOSE
!*   Root process reads the location of primal nodes
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Root process reads the location of primal nodes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)         :: fname            !> file name
    INTEGER                     :: junk             !> junk IO variable
    INTEGER                     :: k                !> simplicial order
    INTEGER                     :: ip,ie            !> loop indices (procs,elems)
    INTEGER                     :: buffer_size      !> size of MPI comm buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE :: real_buffer(:)   !> buffer for MPI comms
    INTEGER                     :: iib,fib          !> initial/final index for MPI comm buffer
    INTEGER, ALLOCATABLE        :: status(:)        !> size of MPI comm buffer
    INTEGER                     :: ptr, ptrp        !>
    LOGICAL                     :: err=.FALSE.    !> error variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- open mesh file and check for errors --
    fname=trim(mesh_prefix) // ".node"; CALL root_open_file_read(fname,mesh_unit)

    !-- set buffer size and allocate buffer --
    buffer_size=maxval(num_pelm_pp)*(dim_embbd+1)
    ALLOCATE(real_buffer(buffer_size))

    !-- allocate storage for nodal locations --
    ALLOCATE(lcl_complex(1)%centers(num_elm(1),dim_embbd))

    IF (rank == root) THEN
      !-- on root process --

      !-- read header and check for consistency --
      READ(mesh_unit,*) junk, k
      IF (glb_num_elm(1) /= junk) THEN
        WRITE(log_unit,'(A,A,A,I5,A,I5)')'Error: mesh file',trim(fname),&
          ': inconsistency: glb_num_elm(1)=',glb_num_elm(1),'/=',junk
        err=.TRUE.
      ELSEIF (dim_embbd /= k) THEN
        WRITE(log_unit,'(A,A,A,I5,A,I5)')'Error: mesh file',trim(fname),&
          ': inconsistcy: dim_embbd=',dim_embbd,'/=',k
        err=.TRUE.
      END IF
    END IF

    !-- check for errors
    CALL MPI_BCAST(err,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    IF (err) CALL end_mpi()

    IF (rank == root) THEN
      !-- on root process --

      !-- advance to neighbour data --
      DO ie=1,num_pelm_pp(1)
        READ(mesh_unit,*)
      END DO

      !-- loop through neighbouring processes --
      DO ip=1,num_procs-1
        !-- loop through elements on neighbouring process --
        iib=1; fib=dim_embbd+1
        DO ie=1,num_pelm_pp(ip+1)
          !-- read nodal location of elements --
          READ(mesh_unit,*) real_buffer(iib:fib)

          !-- update indices --
          iib=iib+dim_embbd+1; fib=fib+dim_embbd+1
        END DO

        !-- send nodal locations (non-blocking) --
        buffer_size=num_pelm_pp(ip+1)*(dim_embbd+1)
        CALL MPI_SEND(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,ip,ip,&
          MPI_COMM_WORLD,ier)
      END DO

      !-- read local data --
      REWIND(mesh_unit)
      READ(mesh_unit,*) !-- header --

      !-- loop through local primary elements --
      ptr=1
      DO ie=1,num_pelm_pp(1)
        !-- read nodal location of elements --
        READ(mesh_unit,*) junk, lcl_complex(1)%centers(ptr,:)

        !-- find correct location in data structure --
        ptrp=ptr; DO
          IF (junk+indx_offset(1) == lcl_complex(1)%node_indx(ptrp,1)) EXIT
          ptrp=ptrp+1
        END DO

        !-- shuffle data and zero locations from neighbouring processes --
        IF (ptrp>ptr) THEN
          lcl_complex(1)%centers(ptrp,:)=lcl_complex(1)%centers(ptr,:)
          lcl_complex(1)%centers(ptr:ptrp-1,:)=0.d0
        END IF
        ptr=ptrp+1
      END DO
      !-- zero remaining locations from neighbouring processes --
      lcl_complex(1)%centers(ptr:num_elm(1),:)=0.d0

    !-- on neighbouring process --
    ELSE
      !-- allocate MPI comm status variable --
      ALLOCATE(status(MPI_STATUS_SIZE))

      !-- receive nodal indices (blocking) --
      buffer_size=num_pelm_pp(rank+1)*(dim_embbd+1)
      CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,0,rank,&
        MPI_COMM_WORLD,status,ier)

      !-- unpack nodal indices in to local element structure --
      iib=1; fib=dim_embbd+1; ptr=1
      DO ie=1,num_pelm_pp(rank+1)
        !-- read nodal indices of element and temporarily assign to data structure --
        junk=IDNINT(real_buffer(iib))
        lcl_complex(1)%centers(ptr,:)=real_buffer(iib+1:fib)

        !-- find correct location in data structure --
        ptrp=ptr; DO
          IF (junk+indx_offset(1) == lcl_complex(1)%node_indx(ptrp,1)) EXIT
          ptrp=ptrp+1
        END DO

        !-- shuffle data and zero locations from neighbouring processes --
        IF (ptrp>ptr) THEN
          lcl_complex(1)%centers(ptrp,:)=lcl_complex(1)%centers(ptr,:)
          lcl_complex(1)%centers(ptr:ptrp-1,:)=0.d0
        END IF

        !-- update indicies --
        iib=iib+dim_embbd+1; fib=fib+dim_embbd+1; ptr=ptrp+1
      END DO
      !-- zero remaining locations from neighbouring processes --
      lcl_complex(1)%centers(ptr:num_elm(1),:)=0.d0
    END IF

    !-- wait for MPI communication to finish, close file and clean up --
    IF (rank == root) THEN
      CLOSE(mesh_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF
      DEALLOCATE(real_buffer)
    ELSE
      DEALLOCATE(real_buffer,status)
    END IF

    !-- make sure all comms have completed --
    !-- all comms are blocking: no need for MPI_BARRIER --

    RETURN

  END SUBROUTINE
! io_mod|read_nodes_prml
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/read_glb_indx
!* SYNOPSIS
  SUBROUTINE read_glb_indx(k)
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER,INTENT(IN)    :: k                !> simplicial order

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: junk             !> junk IO variable
    INTEGER               :: orient           !> orientation
    INTEGER               :: i,j,l,n          !> loop counters
    INTEGER               :: ptr              !> pointer/counter
    LOGICAL               :: err=.FALSE.      !> error variable
    CHARACTER(LEN=slen)   :: msg              !> error message
    INTEGER, ALLOCATABLE  :: int_buffer(:)    !> integer buffer for MPI comms
    INTEGER               :: buffer_size      !> size of MPI comm buffer
    INTEGER               :: iib, fib         !> buffer array indices
    INTEGER               :: num_read_iters   !> number of read iterations
    INTEGER               :: read_size        !> number of reads
    INTEGER               :: resid_read_size  !> reads in final loop
    INTEGER               :: stride           !> stride in buffer
    INTEGER, ALLOCATABLE  :: prml_elm(:,:)    !> primal element list

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- build file name and open --
    SELECT CASE(k)
    CASE(1); fname=trim(mesh_prefix) // ".node"
    CASE(2); fname=trim(mesh_prefix) // ".edge"
    CASE(3); fname=trim(mesh_prefix) // ".face"
    CASE DEFAULT; WRITE(log_unit,'(A,I5)')'***Error: invalid geometric order'//&
      ' (0<k<dim_cmplx+1): ', k; CALL end_mpi(); END SELECT
    CALL root_open_file_read(fname,mesh_unit)

    !-- read header and broadcast global number of entities --
    IF (rank == root) READ(mesh_unit,*) glb_num_elm(k)
    CALL MPI_BCAST(glb_num_elm(k),1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    !-- get some basic information for reuse --
    read_size=max_read_size; stride=k+1; fib=k
    buffer_size=read_size*stride; num_read_iters=glb_num_elm(k)/read_size
    resid_read_size=mod(glb_num_elm(k),read_size)
    ALLOCATE(int_buffer(buffer_size))

    IF (rank == root) THEN !-- on root process --
      !-- first element to get index offset
      READ(mesh_unit,*) junk; indx_offset(k)=1-junk; BACKSPACE(mesh_unit)

      DO i=0,num_read_iters !-- loop through read iterations --
        IF (i == num_read_iters) THEN !-- values for last read iteration --
          read_size=resid_read_size; buffer_size=read_size*stride; END IF
        iib=1; DO j=1, read_size  !-- read data for current iteration --
          READ(mesh_unit,*) int_buffer(iib:iib+fib)
          !-- apply vertex offset, set global index and update buffer index --
          int_buffer(iib+1:iib+fib)=int_buffer(iib+1:iib+fib)+indx_offset(1)
          int_buffer(iib)=i*max_read_size+j; iib=iib+stride; END DO
        !-- broadcast vertex indices --
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        !--
        iib=1; DO j=1,read_size !-- sort each record --
          CALL int_insertion_sort(int_buffer(iib+1:iib+fib))
          IF (int_buffer(iib+1)>num_elm(1)) THEN
            iib=iib+stride; CYCLE; END IF !-- cycle if not local --
          DO l=1, num_elm(k) !-- search for matching element --
            IF (ALL(int_buffer(iib+1:iib+fib)==&
              lcl_complex(k)%node_indx(l,:))) THEN !-- assign to local array --
              lcl_complex(k)%glb_indx(l)=int_buffer(iib); EXIT; END IF
          END DO; iib=iib+stride
        END DO
      END DO; CLOSE(mesh_unit)
    ELSE !-- on neighbouring processes --
      DO i=0,num_read_iters !-- loop through read iterations --
        IF (i == num_read_iters) THEN !-- values for last read iteration --
          read_size=resid_read_size; buffer_size=read_size*stride; END IF
        !-- recieve broadcast vertex indices --
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        iib=1; DO j=1, read_size !-- sort each record --
          CALL int_insertion_sort(int_buffer(iib+1:iib+fib))
          IF ((elm2proc(int_buffer(iib+1))>rank) .OR. &
           (elm2proc(int_buffer(iib+k))<rank)) THEN;  !-- is local? --
           iib=iib+stride; CYCLE; END IF !-- cycle --
          DO l=1, num_elm(k) !-- search for matching element --
            IF (ALL(int_buffer(iib+1:iib+fib)==&
              lcl_complex(k)%node_indx(l,:))) THEN !-- assign to local array --
              lcl_complex(k)%glb_indx(l)=int_buffer(iib); EXIT;  END IF
          END DO; iib=iib+stride
        END DO
      END DO
    END IF
    DEALLOCATE(int_buffer)

    !-- broadcast index offset variable --
    CALL MPI_BCAST(indx_offset(k),1,MPI_INTEGER,0,MPI_COMM_WORLD,ier); RETURN

  END SUBROUTINE
! io_mod|read_glb_indx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/root_open_file_read
!* SYNOPSIS
  SUBROUTINE root_open_file_read(fname,unit)
!* PURPOSE
!*   Open given file on root process
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Open given file on root process
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: fname  !> file name
    INTEGER, INTENT(IN)           :: unit   !> file unit

    !-- local variables --
    LOGICAL :: err=.FALSE.  !> error variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process: --
    IF (rank == root) THEN !-- Open mesh file and check for errors --
      OPEN(unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
      err=chkerr_log(ier,'OPEN - error opening '//TRIM(fname),log_unit); END IF
    CALL chkerrMPI_BCAST(err) !-- complete error checking --
    RETURN

  END SUBROUTINE
! io_mod|root_open_file_read
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/root_open_file_write
!* SYNOPSIS
  SUBROUTINE root_open_file_write(fname,unit)
!* PURPOSE
!*   Open given file on root process
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Open given file on root process
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)    :: fname          !> file name
    INTEGER, INTENT(IN)             :: unit           !> file unit
    LOGICAL                         :: err=.FALSE.  !> error variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process: --
    IF (rank == root) THEN
      !-- Open mesh file and check for errors --
      OPEN(unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode=',ier
      err=.TRUE.
      END IF
    END IF

    !-- check for errors
    CALL MPI_BCAST(err,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    IF (err) CALL end_mpi()

    RETURN

  END SUBROUTINE
! io_mod|root_open_file_write
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_solution_D0S
!* SYNOPSIS
  SUBROUTINE write_solution_D0S(sclr_name,vctr_name)
!* PURPOSE
!*
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: sclr_name           !> scalar variable name (label)
    CHARACTER(LEN=*), INTENT(IN)  :: vctr_name           !> vector variable name (label)

    !-- local variables --
    CHARACTER(LEN=slen)           :: fname               !> file name
    INTEGER                       :: kmax                !>
    INTEGER                       :: i,j                 !>
    INTEGER                       :: iib,fib             !>
    INTEGER                       :: ptr                 !>
    INTEGER                       :: extra_comms         !>
    INTEGER,ALLOCATABLE           :: int_buffer(:)       !>
    INTEGER,ALLOCATABLE           :: junk_int_buffer(:)  !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: real_buffer(:)      !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: junk_real_buffer(:) !>
    INTEGER                       :: buffer_size         !>
    INTEGER                       :: junk                !>
    INTEGER, ALLOCATABLE          :: req(:)              !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE          :: status(:)           !> size of MPI comm buffer
    INTEGER, ALLOCATABLE          :: istatus(:,:)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- set commonly used variable --
    kmax=dim_cmplx+1

    !-- open solution file --
    fname=trim(sol_prefix) // '.vtk'; CALL root_open_file_write(fname,sol_unit)

    !-- gather total number of primal volumes accross all processes --
    CALL MPI_REDUCE(num_elm(kmax),extra_comms,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)

    !-- Open mesh file and check for errors --
    fname=trim(mesh_prefix) // ".ele"; CALL root_open_file_read(fname,mesh_unit)

    IF (rank==root) THEN
      !-- on root process --

      !-- compute extra comms that will be received --
      extra_comms=extra_comms - glb_num_elm(kmax)

      !-- write header --
      write(sol_unit,'(A)') '# vtk DataFile Version 3.0'
      write(sol_unit,'(A)') 'Unstructured Grid'
      write(sol_unit,'(A)') 'ASCII'
      write(sol_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'

      !-------------------------------------------------------------------------
      !-- write primal vertex locations --
      !-- write header --
      write(sol_unit,'(A,I10,A)') 'POINTS', glb_num_elm(1), ' double'
      !-- write local data --
      IF (dim_embbd==3) THEN
        DO i=1,num_pelm_pp(1)
          write(sol_unit,*) lcl_complex(1)%centers(i,:)
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=1,num_pelm_pp(1)
          write(sol_unit,*) lcl_complex(1)%centers(i,:), 0.d0
        END DO
      END IF
      !-- collect data from adjacent processes and write to file --
      !-- allocate MPI variables --
      ALLOCATE(real_buffer(3*num_pelm_pp(1)),status(MPI_STATUS_SIZE))
      IF (dim_embbd==3) THEN
        DO i=1,num_procs-1
          !-- recieve data from adjacent process --
          buffer_size=dim_embbd*num_pelm_pp(i+1)
          CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,i,0,&
            MPI_COMM_WORLD,status,ier)
          !-- unpack buffer and write to file --
          iib=1; fib=dim_embbd
          DO j=1,num_pelm_pp(i+1)
            write(sol_unit,*) real_buffer(iib:fib)
            iib=iib+dim_embbd; fib=fib+dim_embbd
          END DO
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=1,num_procs-1
          !-- recieve data from adjacent process --
          buffer_size=dim_embbd*num_pelm_pp(i+1)
          CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,i,0,&
            MPI_COMM_WORLD,status,ier)
          !-- unpack buffer and write to file --
          iib=1; fib=dim_embbd
          DO j=1,num_pelm_pp(i+1)
            write(sol_unit,*) real_buffer(iib:fib), 0.d0
            iib=iib+dim_embbd; fib=fib+dim_embbd
          END DO
        END DO
      END IF
      !-- clean up --
      DEALLOCATE(real_buffer)

      !-------------------------------------------------------------------------
      !-- write primal volume indices --
      !-- write header --
      write(sol_unit,'(A,I10,I10)') 'CELLS', glb_num_elm(kmax), (kmax+1)*glb_num_elm(kmax)

      !-- use input mesh file for node order --
      !-- read header --
      READ(mesh_unit,*) junk
      !-- process node indices --
      ALLOCATE(int_buffer(kmax))
      DO i=1,glb_num_elm(kmax)
        READ(mesh_unit,*) junk, int_buffer
        write(sol_unit,*) kmax, int_buffer+indx_offset(1) - 1
      END DO

      !-------------------------------------------------------------------------
      !-- write primal volumes types --
      !-- write header --
      write(sol_unit,'(A,I10)') 'CELL_TYPES', glb_num_elm(kmax)
      !-- write data --
      SELECT CASE(dim_cmplx)
      CASE(1)
        DO i=1,glb_num_elm(2); write(sol_unit,'(I1)') 3; END DO
      CASE(2)
        DO i=1,glb_num_elm(3); write(sol_unit,'(I1)') 5; END DO
      CASE(3)
        DO i=1,glb_num_elm(4); write(sol_unit,'(I2)') 10; END DO
      CASE DEFAULT
      END SELECT

      !-------------------------------------------------------------------------
      !-- write primal volumes scalar data --
      !-- write header --
      write(sol_unit,'(A,I10)') 'CELL_DATA', glb_num_elm(kmax)
      write(sol_unit,'(A,A,A)') 'SCALARS ',sclr_name,' double 1'
      write(sol_unit,'(A)') 'LOOKUP_TABLE default'
      !-- set local data pointer and allocate MPI variables --
      ptr=1
      ALLOCATE(real_buffer(1),junk_real_buffer(1))
      !-- wait for previous communication to end and deallocate old MPI variables --
      ! CALL MPI_WAITALL(extra_comms,req,istatus,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      ! DEALLOCATE(int_buffer,junk_int_buffer)
      !-- deallocate old and allocate new MPI variables --
      DEALLOCATE(int_buffer)
      ALLOCATE(req(extra_comms),istatus(MPI_STATUS_SIZE,extra_comms))
      !-- loop through all primal volumes and write data to file --
      DO i=1,glb_num_elm(kmax)
        IF (lcl_complex(kmax)%glb_indx(ptr)==i) THEN   !-- data is local --
          !-- write to file --
          write(sol_unit,*) lcl_complex(kmax)% dual_sol(ptr,1)
          !-- update local data pointer --
          ptr=min(ptr+1,num_elm(kmax))
        ELSE   !-- data is non-local --
          !-- recieve first communication to arrive for given volume --
          CALL MPI_RECV(real_buffer,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,i,&
            MPI_COMM_WORLD,status,ier)
          !-- write to file --
          write(sol_unit,*) real_buffer
        END IF
      END DO

      !-- recieve extras unneeded data to be discarded --
      DO i=1,extra_comms
        CALL MPI_IRECV(junk_real_buffer,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
          MPI_ANY_TAG,MPI_COMM_WORLD,req(i),ier)
      END DO

      !-------------------------------------------------------------------------
      !-- write cell center vector data --
      !-- write header --
      write(sol_unit,'(A,A,A)') 'VECTORS ',vctr_name,' double'
      !-- set local data pointer and allocate MPI variables --
      ptr=1
      DEALLOCATE(real_buffer,junk_real_buffer)
      ALLOCATE(real_buffer(dim_embbd),junk_real_buffer(dim_embbd))
      !-- wait for previous communication to end and deallocate old MPI variables --
      CALL MPI_WAITALL(extra_comms,req,istatus,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      !-- loop through all primal volumes and write data to file --
      DO i=1,glb_num_elm(kmax)
        IF (lcl_complex(kmax)%glb_indx(ptr)==i) THEN   !-- data is local --
          !-- write to file --
          write(sol_unit,*) lcl_complex(kmax)% whtny_sol(ptr,:)
          !-- update local data pointer --
          ptr=min(ptr+1,num_elm(kmax))
        ELSE   !-- data is non-local --
          !-- recieve first communication to arrive for given volume --
          CALL MPI_RECV(real_buffer,dim_embbd,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,i,&
            MPI_COMM_WORLD,status,ier)
          !-- write to file --
          write(sol_unit,*) real_buffer
        END IF
      END DO

      !-- recieve extras unneeded data to be discarded --
      DO i=1,extra_comms
        CALL MPI_IRECV(junk_real_buffer,dim_embbd,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
          MPI_ANY_TAG,MPI_COMM_WORLD,req(i),ier)
      END DO

      !-------------------------------------------------------------------------
      !-- wait for previous communication to end and deallocate variables --
      CALL MPI_WAITALL(extra_comms,req,istatus,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      DEALLOCATE(real_buffer,junk_real_buffer,req,status,istatus)

      !-- close files --
      CLOSE(sol_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF
      CLOSE(mesh_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF

    ELSE
      !-------------------------------------------------------------------------
      !-- find starting index of local primary vertices --
      DO i=1,num_elm(1)
        IF (lcl_complex(1)% glb_indx(i)==glb_offset+1) THEN
          ptr=i; EXIT
        END IF
      END DO
      buffer_size=dim_embbd*num_pelm_pp(rank+1); ALLOCATE(real_buffer(buffer_size))
      iib=1; fib=dim_embbd
      DO j=1,num_pelm_pp(rank+1)
        real_buffer(iib:fib)=lcl_complex(1)%centers(ptr+j-1,:)
        iib=iib+dim_embbd; fib=fib+dim_embbd
      END DO

      !-------------------------------------------------------------------------
      !-- send primal vertex locations --
      CALL MPI_ISEND(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,0,0,&
        MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)

      !-------------------------------------------------------------------------
      !-- send primal volume indices --
      ! DO i=1,num_elm(kmax)
      !   CALL MPI_ISEND(lcl_complex(kmax)%node_indx(i,1:kmax),kmax,MPI_INTEGER,0,&
      !     lcl_complex(kmax)%glb_indx(i),MPI_COMM_WORLD,junk,ier)
      !   CALL MPI_REQUEST_FREE(junk,ier)
      ! END DO

      !-- wait for previous communications to finish (free tags) --
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier); DEALLOCATE(real_buffer)

      !-------------------------------------------------------------------------
      !-- send primal volume values --
      DO i=1,num_elm(kmax)
        CALL MPI_ISEND(lcl_complex(kmax)%dual_sol(i,1),1,MPI_DOUBLE_PRECISION,0,&
          lcl_complex(kmax)%glb_indx(i),MPI_COMM_WORLD,junk,ier)
        CALL MPI_REQUEST_FREE(junk,ier)
      END DO

      !-- wait for previous communications to finish --
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      !-------------------------------------------------------------------------
      !-- send primal volume values --
      DO i=1,num_elm(kmax)
        CALL MPI_ISEND(lcl_complex(kmax)%whtny_sol(i,:),dim_embbd,MPI_DOUBLE_PRECISION,0,&
          lcl_complex(kmax)%glb_indx(i),MPI_COMM_WORLD,junk,ier)
        CALL MPI_REQUEST_FREE(junk,ier)
      END DO

      !-- wait for previous communications to finish --
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    END IF

    RETURN

  END SUBROUTINE
! io_mod|write_solution_D0S
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_solution_D0S2
!* SYNOPSIS
  SUBROUTINE write_solution_D0S2(sclr_name,vctr_name)
!* PURPOSE
!*
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: sclr_name           !> scalar variable name (label)
    CHARACTER(LEN=*), INTENT(IN)  :: vctr_name           !> vector variable name (label)

    !-- local variables --
    CHARACTER(LEN=slen)           :: fname               !> file name
    INTEGER                       :: kmax                !>
    INTEGER                       :: i,j,k               !>
    INTEGER                       :: iib,fib             !>
    INTEGER                       :: ptr                 !>
    INTEGER,ALLOCATABLE           :: int_buffer(:)       !>
    INTEGER,ALLOCATABLE           :: junk_int_buffer(:)  !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: real_buffer(:)      !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: junk_real_buffer(:) !>
    INTEGER                       :: buffer_size         !>
    INTEGER                       :: junk                !>
    INTEGER, ALLOCATABLE          :: req(:)              !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE          :: status(:)           !> size of MPI comm buffer
    INTEGER, ALLOCATABLE          :: istatus(:,:)        !> size of MPI comm buffer
    REAL(KIND=PGMSiwp)                :: tmp_time, tmp_time_int !> temporary timing variable
    INTEGER                       :: comm_size,max_comm_size,gindx,num_comm_iters, &
                                     strt_indx,resid_comm_size                 !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- set commonly used variable --
    kmax=dim_cmplx+1

    !-- open solution file --
    fname=trim(sol_prefix) // '.vtk'; CALL root_open_file_write(fname,sol_unit)

    !-- get some basic information for reuse --
    max_comm_size=10000

    IF (rank==root) THEN
      !-- on root process --

      !-------------------------------------------------------------------------
      !-- write header --
      write(sol_unit,'(A)') '# vtk DataFile Version 3.0'
      write(sol_unit,'(A)') 'Unstructured Grid'
      write(sol_unit,'(A)') 'ASCII'
      write(sol_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'

      !-------------------------------------------------------------------------
      !-- write primal vertex locations --
      !-- write header --
      write(sol_unit,'(A,I10,A)') 'POINTS', glb_num_elm(1), ' double'
      !-- write local data --
      IF (dim_embbd==3) THEN
        DO i=1,num_pelm_pp(1)
          write(sol_unit,*) lcl_complex(1)%centers(i,:)
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=1,num_pelm_pp(1)
          write(sol_unit,*) lcl_complex(1)%centers(i,:), 0.d0
        END DO
      END IF
      !-- collect data from adjacent processes and write to file --
      !-- allocate MPI variables --
      ALLOCATE(status(MPI_STATUS_SIZE))
      IF (dim_embbd==3) THEN
        DO k=1,num_procs-1

          num_comm_iters=num_pelm_pp(k+1)/max_comm_size
          resid_comm_size=mod(num_pelm_pp(k+1),max_comm_size)
          comm_size=max_comm_size
          buffer_size=comm_size*dim_embbd
          ALLOCATE(real_buffer(buffer_size))

          DO i=0,num_comm_iters
            IF (i == num_comm_iters) THEN
              comm_size=resid_comm_size
              buffer_size=comm_size*dim_embbd
            END IF

            !-- recieve data from adjacent process --
            CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,k,i,&
              MPI_COMM_WORLD,status,ier)

            !-- unpack buffer and write to file --
            iib=1; fib=dim_embbd
            DO j=1,comm_size
              write(sol_unit,*) real_buffer(iib:fib)
              iib=iib+dim_embbd; fib=fib+dim_embbd
            END DO
          END DO
          DEALLOCATE(real_buffer)
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO k=1,num_procs-1

          num_comm_iters=num_pelm_pp(k+1)/max_comm_size
          resid_comm_size=mod(num_pelm_pp(k+1),max_comm_size)
          comm_size=max_comm_size
          buffer_size=comm_size*dim_embbd
          ALLOCATE(real_buffer(buffer_size))

          DO i=0,num_comm_iters
            IF (i == num_comm_iters) THEN
              comm_size=resid_comm_size
              buffer_size=comm_size*dim_embbd
            END IF

            !-- recieve data from adjacent process --
            CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,k,i,&
              MPI_COMM_WORLD,status,ier)

            !-- unpack buffer and write to file --
            iib=1; fib=dim_embbd
            DO j=1,comm_size
              write(sol_unit,*) real_buffer(iib:fib), 0.d0
              iib=iib+dim_embbd; fib=fib+dim_embbd
            END DO
          END DO
          DEALLOCATE(real_buffer)
        END DO
      END IF
      !-- clean up --
      DEALLOCATE(status)

      !-------------------------------------------------------------------------
      !-- write primal volume indices --
      !-- write header --
      write(sol_unit,'(A,I10,I10)') 'CELLS', glb_num_elm(kmax), (kmax+1)*glb_num_elm(kmax)

      !--
      num_comm_iters=glb_num_elm(kmax)/max_comm_size
      resid_comm_size=mod(glb_num_elm(kmax),max_comm_size)
      comm_size=max_comm_size
      buffer_size=comm_size*kmax
      ALLOCATE(int_buffer(buffer_size),junk_int_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*kmax
        END IF

        int_buffer=0
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          iib=(gindx-1)*kmax+1;
          fib=iib+kmax-1
          IF (lcl_complex(kmax)%orientation(j) == 1) THEN
            int_buffer(iib:fib)=lcl_complex(kmax)%node_indx(j,1:kmax)-1
          ELSE
            int_buffer(iib)=lcl_complex(kmax)%node_indx(j,2)-1
            int_buffer(iib+1)=lcl_complex(kmax)%node_indx(j,1)-1
            int_buffer(iib+2:fib)=lcl_complex(kmax)%node_indx(j,3:kmax)-1
          END IF
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(int_buffer,junk_int_buffer,buffer_size,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)

        iib=1; fib=kmax-1
        DO j=1,comm_size
          WRITE(sol_unit,*) kmax, junk_int_buffer(iib:iib+fib)
          iib=iib+kmax
        END DO
      END DO

      !-- clean up --
      DEALLOCATE(int_buffer, junk_int_buffer)

      !-------------------------------------------------------------------------
      !-- write primal volumes types --
      !-- write header --
      write(sol_unit,'(A,I10)') 'CELL_TYPES', glb_num_elm(kmax)
      !-- write data --
      SELECT CASE(dim_cmplx)
      CASE(1)
        DO i=1,glb_num_elm(kmax); write(sol_unit,'(I1)') 3; END DO
      CASE(2)
        DO i=1,glb_num_elm(kmax); write(sol_unit,'(I1)') 5; END DO
      CASE(3)
        DO i=1,glb_num_elm(kmax); write(sol_unit,'(I2)') 10; END DO
      CASE DEFAULT
      END SELECT

      !-------------------------------------------------------------------------
      !-- write primal volumes scalar data --
      !-- write header --
      write(sol_unit,'(A,I10)') 'CELL_DATA', glb_num_elm(kmax)
      write(sol_unit,'(A,A,A)') 'SCALARS ',sclr_name,' double 1'
      write(sol_unit,'(A)') 'LOOKUP_TABLE default'

      !--
      comm_size=max_comm_size
      buffer_size=comm_size
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size
        END IF

        real_buffer=-large
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          real_buffer(gindx)=lcl_complex(kmax)% dual_sol(j,1)
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

        DO j=1,comm_size
          WRITE(sol_unit,*) junk_real_buffer(j)
        END DO
      END DO

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

      !-------------------------------------------------------------------------
      !-- write cell center vector data --
      !-- write header --
      write(sol_unit,'(A,A,A)') 'VECTORS ',vctr_name,' double'

      !--
      comm_size=max_comm_size
      buffer_size=comm_size*dim_embbd
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      IF (dim_embbd==3) THEN
        DO i=0,num_comm_iters
          IF (i == num_comm_iters) THEN
            comm_size=resid_comm_size
            buffer_size=comm_size*dim_embbd
          END IF

          real_buffer=-large
          DO j=strt_indx,num_elm(kmax)
            gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
            IF (gindx>comm_size) THEN
              strt_indx=j
              EXIT
            END IF
            iib=(gindx-1)*dim_embbd+1;
            fib=iib+dim_embbd-1
            real_buffer(iib:fib)=lcl_complex(kmax)%whtny_sol(j,1:dim_embbd)
            IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
          END DO
          CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

          iib=1; fib=dim_embbd-1
          DO j=1,comm_size
            WRITE(sol_unit,*) junk_real_buffer(iib:iib+fib)
            iib=iib+dim_embbd
          END DO
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=0,num_comm_iters
          IF (i == num_comm_iters) THEN
            comm_size=resid_comm_size
            buffer_size=comm_size*dim_embbd
          END IF

          real_buffer=-large
          DO j=strt_indx,num_elm(kmax)
            gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
            IF (gindx>comm_size) THEN
              strt_indx=j
              EXIT
            END IF
            iib=(gindx-1)*dim_embbd+1;
            fib=iib+dim_embbd-1
            real_buffer(iib:fib)=lcl_complex(kmax)%whtny_sol(j,1:dim_embbd)
            IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
          END DO
          CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

          iib=1; fib=dim_embbd-1
          DO j=1,comm_size
            WRITE(sol_unit,*) junk_real_buffer(iib:iib+fib), 0.d0
            iib=iib+dim_embbd
          END DO
        END DO
      END IF

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

      !-------------------------------------------------------------------------
      !-- close files --
      CLOSE(sol_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF

    ELSE

      !-- adjacent processes --
      !-------------------------------------------------------------------------
      !-- send header --

      !-------------------------------------------------------------------------
      !-- send primal vertex locations --

      !-- find starting index of local primary vertices --
      DO i=1,num_elm(1)
        IF (lcl_complex(1)% glb_indx(i)==glb_offset+1) THEN
          ptr=i-1; EXIT
        END IF
      END DO

      !-- set up and pack buffer --
      num_comm_iters=num_pelm_pp(rank+1)/max_comm_size
      resid_comm_size=mod(num_pelm_pp(rank+1),max_comm_size)
      comm_size=max_comm_size
      buffer_size=comm_size*dim_embbd
      ALLOCATE(real_buffer(buffer_size))

      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*dim_embbd
        END IF

        iib=1; fib=dim_embbd-1
        DO j=1,comm_size
          real_buffer(iib:iib+fib)=lcl_complex(1)%centers(ptr+j,:)
          iib=iib+dim_embbd
        END DO

        !-- send primal vertex locations --
        CALL MPI_SEND(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,0,i,MPI_COMM_WORLD,ier)
      END DO

      DEALLOCATE(real_buffer)

      !-------------------------------------------------------------------------
      !-- send primal volume indices --

      !--
      num_comm_iters=glb_num_elm(kmax)/max_comm_size
      resid_comm_size=mod(glb_num_elm(kmax),max_comm_size)
      comm_size=max_comm_size
      buffer_size=comm_size*kmax
      ALLOCATE(int_buffer(buffer_size),junk_int_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*kmax
        END IF

        int_buffer=0
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          iib=(gindx-1)*kmax+1;
          fib=iib+kmax-1
          IF (lcl_complex(kmax)%orientation(j) == 1) THEN
            int_buffer(iib:fib)=lcl_complex(kmax)%node_indx(j,1:kmax)-1
          ELSE
            int_buffer(iib)=lcl_complex(kmax)%node_indx(j,2)-1
            int_buffer(iib+1)=lcl_complex(kmax)%node_indx(j,1)-1
            int_buffer(iib+2:fib)=lcl_complex(kmax)%node_indx(j,3:kmax)-1
          END IF
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(int_buffer,junk_int_buffer,buffer_size,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)
      END DO

      !-- clean up --
      DEALLOCATE(int_buffer, junk_int_buffer)

      !-------------------------------------------------------------------------
      !-- send primal volumes types --

      !-------------------------------------------------------------------------
      !-- send primal volume scalar values --

      !--
      comm_size=max_comm_size
      buffer_size=comm_size
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size
        END IF

        real_buffer=-large
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          real_buffer(gindx)=lcl_complex(kmax)% dual_sol(j,1)
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

      END DO

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

      !-------------------------------------------------------------------------
      !-- send cell center vector data --

      !--
      comm_size=max_comm_size
      buffer_size=comm_size*dim_embbd
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*dim_embbd
        END IF

        real_buffer=-large
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          iib=(gindx-1)*dim_embbd+1;
          fib=iib+dim_embbd-1
          real_buffer(iib:fib)=lcl_complex(kmax)%whtny_sol(j,1:dim_embbd)
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

      END DO

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

    END IF

    !-- wait for previous communications to finish --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|write_solution_D0S2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_solution_D0S
!* SYNOPSIS
  SUBROUTINE write_unsteady_D0S(sclr_name,vctr_name,iter)
!* PURPOSE
!*
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: sclr_name           !> scalar variable name (label)
    CHARACTER(LEN=*), INTENT(IN)  :: vctr_name           !> vector variable name (label)
    INTEGER, INTENT(IN)           :: iter                !> integer iteration number

    !-- local variables --
    CHARACTER(LEN=slen)           :: fname               !> file name
    INTEGER                       :: prelen              !>
    INTEGER                       :: kmax                !>
    INTEGER                       :: i,j                 !>
    INTEGER                       :: iib,fib             !>
    INTEGER                       :: ptr                 !>
    INTEGER                       :: extra_comms         !>
    INTEGER,ALLOCATABLE           :: int_buffer(:)       !>
    INTEGER,ALLOCATABLE           :: junk_int_buffer(:)  !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: real_buffer(:)      !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: junk_real_buffer(:) !>
    INTEGER                       :: buffer_size         !>
    INTEGER                       :: junk                !>
    INTEGER, ALLOCATABLE          :: req(:)              !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE          :: status(:)           !> size of MPI comm buffer
    INTEGER, ALLOCATABLE          :: istatus(:,:)        !> size of MPI comm buffer
    REAL(KIND=PGMSiwp)        :: tmp_time, tmp_time_int !> temporary timing variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- set commonly used variable --
    kmax=dim_cmplx+1

    !-- open solution file --
    prelen=len_trim(unstdy_prefix)
    fname=trim(unstdy_prefix)//"_00000000.vtk"
    WRITE(fname(prelen+1:prelen+9),'(I9)') 100000000+iter
    fname(prelen+1:prelen+1)="_"
    CALL root_open_file_write(fname,unstdy_unit)

    !-- gather total number of primal volumes accross all processes --
    CALL MPI_REDUCE(num_elm(kmax),extra_comms,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)

    !-- Open mesh file and check for errors --
    fname=trim(mesh_prefix) // ".ele"; CALL root_open_file_read(fname,mesh_unit)

    IF (rank==root) THEN
      !-- on root process --

      !-- compute extra comms that will be received --
      extra_comms=extra_comms - glb_num_elm(kmax)

      !-- write header --
      write(unstdy_unit,'(A)') '# vtk DataFile Version 3.0'
      write(unstdy_unit,'(A)') 'Unstructured Grid'
      write(unstdy_unit,'(A)') 'ASCII'
      write(unstdy_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'

      !-------------------------------------------------------------------------
      !-- write primal vertex locations --
      !-- write header --
      write(unstdy_unit,'(A,I10,A)') 'POINTS', glb_num_elm(1), ' double'
      !-- write local data --
      IF (dim_embbd==3) THEN
        DO i=1,num_pelm_pp(1)
          write(unstdy_unit,*) lcl_complex(1)%centers(i,:)
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=1,num_pelm_pp(1)
          write(unstdy_unit,*) lcl_complex(1)%centers(i,:), 0.d0
        END DO
      END IF
      !-- collect data from adjacent processes and write to file --
      !-- allocate MPI variables --
      ALLOCATE(real_buffer(3*num_pelm_pp(1)),status(MPI_STATUS_SIZE))
      IF (dim_embbd==3) THEN
        DO i=1,num_procs-1
          !-- recieve data from adjacent process --
          buffer_size=dim_embbd*num_pelm_pp(i+1)
          CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,i,0,&
            MPI_COMM_WORLD,status,ier)
          !-- unpack buffer and write to file --
          iib=1; fib=dim_embbd
          DO j=1,num_pelm_pp(i+1)
            write(unstdy_unit,*) real_buffer(iib:fib)
            iib=iib+dim_embbd; fib=fib+dim_embbd
          END DO
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=1,num_procs-1
          !-- recieve data from adjacent process --
          buffer_size=dim_embbd*num_pelm_pp(i+1)
          CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,i,0,&
            MPI_COMM_WORLD,status,ier)
          !-- unpack buffer and write to file --
          iib=1; fib=dim_embbd
          DO j=1,num_pelm_pp(i+1)
            write(unstdy_unit,*) real_buffer(iib:fib), 0.d0
            iib=iib+dim_embbd; fib=fib+dim_embbd
          END DO
        END DO
      END IF
      !-- clean up --
      DEALLOCATE(real_buffer)

      !-------------------------------------------------------------------------
      !-- write primal volume indices --
      !-- write header --
      write(unstdy_unit,'(A,I10,I10)') 'CELLS', glb_num_elm(kmax), (kmax+1)*glb_num_elm(kmax)

      !-- use input mesh file for node order --
      !-- read header --
      READ(mesh_unit,*) junk
      !-- process node indices --
      ALLOCATE(int_buffer(kmax))
      DO i=1,glb_num_elm(kmax)
        READ(mesh_unit,*) junk, int_buffer
        write(unstdy_unit,*) kmax, int_buffer+indx_offset(1) - 1
      END DO

      ! !-- set local data pointer and allocate MPI variables --
      ! ptr=1
      ! ALLOCATE(&
      !   req(extra_comms),istatus(MPI_STATUS_SIZE,extra_comms),&
      !   int_buffer(kmax),junk_int_buffer(kmax))
      ! !-- wait for previous communication to end --
      ! CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      ! !-- loop through all primal volumes and write data to file --
      ! DO i=1,glb_num_elm(kmax)
      !   IF (lcl_complex(kmax)%glb_indx(ptr)==i) THEN !-- data is local --
      !     !-- write to file --
      !     write(unstdy_unit,*) kmax, lcl_complex(kmax)% node_indx(ptr,:)-1
      !     !-- update local data pointer --
      !     ptr=min(ptr+1,num_elm(kmax))
      !   ELSE !-- non-local data --
      !     !-- receive first communication to arrive for given volume --
      !     CALL MPI_RECV(int_buffer,kmax,MPI_INTEGER,MPI_ANY_SOURCE,i,&
      !       MPI_COMM_WORLD,status,ier)
      !     !-- write to file --
      !     write(unstdy_unit,*) kmax, int_buffer-1
      !   END IF
      ! END DO
      !
      ! !-- recieve extras unneeded data to be discarded --
      ! DO i=1,extra_comms
      !   CALL MPI_IRECV(junk_int_buffer,kmax,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,&
      !     MPI_COMM_WORLD,req(i),ier)
      ! END DO

      !-------------------------------------------------------------------------
      !-- write primal volumes types --
      !-- write header --
      write(unstdy_unit,'(A,I10)') 'CELL_TYPES', glb_num_elm(kmax)
      !-- write data --
      SELECT CASE(dim_cmplx)
      CASE(1)
        DO i=1,glb_num_elm(2); write(unstdy_unit,'(I1)') 3; END DO
      CASE(2)
        DO i=1,glb_num_elm(3); write(unstdy_unit,'(I1)') 5; END DO
      CASE(3)
        DO i=1,glb_num_elm(4); write(unstdy_unit,'(I2)') 10; END DO
      CASE DEFAULT
      END SELECT

      !-------------------------------------------------------------------------
      !-- write primal volumes scalar data --
      !-- write header --
      write(unstdy_unit,'(A,I10)') 'CELL_DATA', glb_num_elm(kmax)
      write(unstdy_unit,'(A,A,A)') 'SCALARS ',sclr_name,' double 1'
      write(unstdy_unit,'(A)') 'LOOKUP_TABLE default'
      !-- set local data pointer and allocate MPI variables --
      ptr=1
      ALLOCATE(real_buffer(1),junk_real_buffer(1))
      !-- wait for previous communication to end and deallocate old MPI variables --
      ! CALL MPI_WAITALL(extra_comms,req,istatus,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      ! DEALLOCATE(int_buffer,junk_int_buffer)
      !-- deallocate old and allocate new MPI variables --
      DEALLOCATE(int_buffer)
      ALLOCATE(req(extra_comms),istatus(MPI_STATUS_SIZE,extra_comms))
      !-- loop through all primal volumes and write data to file --
      DO i=1,glb_num_elm(kmax)
        IF (lcl_complex(kmax)%glb_indx(ptr)==i) THEN   !-- data is local --
          !-- write to file --
          write(unstdy_unit,*) lcl_complex(kmax)% dual_sol(ptr,1)
          !-- update local data pointer --
          ptr=min(ptr+1,num_elm(kmax))
        ELSE   !-- data is non-local --
          !-- recieve first communication to arrive for given volume --
          CALL MPI_RECV(real_buffer,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,i,&
            MPI_COMM_WORLD,status,ier)
          !-- write to file --
          write(unstdy_unit,*) real_buffer
        END IF
      END DO

      !-- recieve extras unneeded data to be discarded --
      DO i=1,extra_comms
        CALL MPI_IRECV(junk_real_buffer,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
          MPI_ANY_TAG,MPI_COMM_WORLD,req(i),ier)
      END DO

      !-------------------------------------------------------------------------
      !-- write cell center vector data --
      !-- write header --
      write(unstdy_unit,'(A,A,A)') 'VECTORS ',vctr_name,' double'
      !-- set local data pointer and allocate MPI variables --
      ptr=1
      !-- wait for previous communication to end and deallocate old MPI variables --
      CALL MPI_WAITALL(extra_comms,req,istatus,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      DEALLOCATE(real_buffer,junk_real_buffer)
      ALLOCATE(real_buffer(dim_embbd),junk_real_buffer(dim_embbd))
      !-- loop through all primal volumes and write data to file --
      DO i=1,glb_num_elm(kmax)
        IF (lcl_complex(kmax)%glb_indx(ptr)==i) THEN   !-- data is local --
          !-- write to file --
          write(unstdy_unit,*) lcl_complex(kmax)% whtny_sol(ptr,:)
          !-- update local data pointer --
          ptr=min(ptr+1,num_elm(kmax))
        ELSE   !-- data is non-local --
          !-- recieve first communication to arrive for given volume --
          CALL MPI_RECV(real_buffer,dim_embbd,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,i,&
            MPI_COMM_WORLD,status,ier)
          !-- write to file --
          write(unstdy_unit,*) real_buffer
        END IF
      END DO

      !-- recieve extras unneeded data to be discarded --
      DO i=1,extra_comms
        CALL MPI_IRECV(junk_real_buffer,dim_embbd,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
          MPI_ANY_TAG,MPI_COMM_WORLD,req(i),ier)
      END DO

      !-------------------------------------------------------------------------
      !-- wait for previous communication to end and deallocate variables --
      CALL MPI_WAITALL(extra_comms,req,istatus,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      DEALLOCATE(real_buffer,junk_real_buffer,req,status,istatus)

      !-- close files --
      CLOSE(unstdy_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF
      CLOSE(mesh_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF

    ELSE
      !-------------------------------------------------------------------------
      !-- find starting index of local primary vertices --
      DO i=1,num_elm(1)
        IF (lcl_complex(1)% glb_indx(i)==glb_offset+1) THEN
          ptr=i; EXIT
        END IF
      END DO
      buffer_size=dim_embbd*num_pelm_pp(rank+1); ALLOCATE(real_buffer(buffer_size))
      iib=1; fib=dim_embbd
      DO j=1,num_pelm_pp(rank+1)
        real_buffer(iib:fib)=lcl_complex(1)%centers(ptr+j-1,:)
        iib=iib+dim_embbd; fib=fib+dim_embbd
      END DO

      !-------------------------------------------------------------------------
      !-- send primal vertex locations --
      CALL MPI_ISEND(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,0,0,&
        MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)

      !-- wait for previous communications to finish --
      !CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      !-------------------------------------------------------------------------
      !-- send primal volume indices --
      ! DO i=1,num_elm(kmax)
      !   CALL MPI_ISEND(lcl_complex(kmax)%node_indx(i,1:kmax),kmax,MPI_INTEGER,0,&
      !     lcl_complex(kmax)%glb_indx(i),MPI_COMM_WORLD,junk,ier)
      !   CALL MPI_REQUEST_FREE(junk,ier)
      ! END DO

      !-- wait for previous communications to finish (free tags) --
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier); DEALLOCATE(real_buffer)

      !-------------------------------------------------------------------------
      !-- send primal volume values --
      DO i=1,num_elm(kmax)
        CALL MPI_ISEND(lcl_complex(kmax)%dual_sol(i,1),1,MPI_DOUBLE_PRECISION,0,&
          lcl_complex(kmax)%glb_indx(i),MPI_COMM_WORLD,junk,ier)
        CALL MPI_REQUEST_FREE(junk,ier)
      END DO

      !-- wait for previous communications to finish --
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      !-------------------------------------------------------------------------
      !-- send cell center vector data --
      DO i=1,num_elm(kmax)
        CALL MPI_ISEND(lcl_complex(kmax)%whtny_sol(i,1:dim_embbd),dim_embbd,MPI_DOUBLE_PRECISION,0,&
          lcl_complex(kmax)%glb_indx(i),MPI_COMM_WORLD,junk,ier)
        CALL MPI_REQUEST_FREE(junk,ier)
      END DO

      !-- wait for previous communications to finish --
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    END IF

    RETURN

  END SUBROUTINE
! io_mod|write_unsteady_D0S
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_solution_D0S2
!* SYNOPSIS
  SUBROUTINE write_unsteady_D0S2(sclr_name,vctr_name,iter)
!* PURPOSE
!*
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: sclr_name           !> scalar variable name (label)
    CHARACTER(LEN=*), INTENT(IN)  :: vctr_name           !> vector variable name (label)
    INTEGER, INTENT(IN)           :: iter                !> integer iteration number

    !-- local variables --
    CHARACTER(LEN=slen)           :: fname               !> file name
    INTEGER                       :: prelen              !>
    INTEGER                       :: kmax                !>
    INTEGER                       :: i,j,k               !>
    INTEGER                       :: iib,fib             !>
    INTEGER                       :: ptr                 !>
    INTEGER                       :: extra_comms         !>
    INTEGER,ALLOCATABLE           :: int_buffer(:)       !>
    INTEGER,ALLOCATABLE           :: junk_int_buffer(:)  !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: real_buffer(:)      !>
    REAL(KIND=PGMSiwp),ALLOCATABLE    :: junk_real_buffer(:) !>
    INTEGER                       :: buffer_size         !>
    INTEGER                       :: junk                !>
    INTEGER, ALLOCATABLE          :: req(:)              !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE          :: status(:)           !> size of MPI comm buffer
    INTEGER, ALLOCATABLE          :: istatus(:,:)        !> size of MPI comm buffer
    REAL(KIND=PGMSiwp)                :: tmp_time, tmp_time_int !> temporary timing variable
    INTEGER                       :: comm_size,max_comm_size,gindx,num_comm_iters, &
                                     strt_indx,resid_comm_size                 !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- set commonly used variable --
    kmax=dim_cmplx+1

    !-- open solution file --
    prelen=len_trim(unstdy_prefix)
    fname=trim(unstdy_prefix)//"_00000000.vtk"
    WRITE(fname(prelen+1:prelen+9),'(I9)') 100000000+iter
    fname(prelen+1:prelen+1)="_"
    CALL root_open_file_write(fname,unstdy_unit)

    !-- get some basic information for reuse --
    max_comm_size=10000

    IF (rank==root) THEN
      !-- on root process --

      !-------------------------------------------------------------------------
      !-- write header --
      write(unstdy_unit,'(A)') '# vtk DataFile Version 3.0'
      write(unstdy_unit,'(A)') 'Unstructured Grid'
      write(unstdy_unit,'(A)') 'ASCII'
      write(unstdy_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'

      !-------------------------------------------------------------------------
      !-- write primal vertex locations --
      !-- write header --
      write(unstdy_unit,'(A,I10,A)') 'POINTS', glb_num_elm(1), ' double'
      !-- write local data --
      IF (dim_embbd==3) THEN
        DO i=1,num_pelm_pp(1)
          write(unstdy_unit,*) lcl_complex(1)%centers(i,:)
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=1,num_pelm_pp(1)
          write(unstdy_unit,*) lcl_complex(1)%centers(i,:), 0.d0
        END DO
      END IF
      !-- collect data from adjacent processes and write to file --
      !-- allocate MPI variables --
      ALLOCATE(status(MPI_STATUS_SIZE))
      IF (dim_embbd==3) THEN
        DO k=1,num_procs-1

          num_comm_iters=num_pelm_pp(k+1)/max_comm_size
          resid_comm_size=mod(num_pelm_pp(k+1),max_comm_size)
          comm_size=max_comm_size
          buffer_size=comm_size*dim_embbd
          ALLOCATE(real_buffer(buffer_size))

          DO i=0,num_comm_iters
            IF (i == num_comm_iters) THEN
              comm_size=resid_comm_size
              buffer_size=comm_size*dim_embbd
            END IF

            !-- recieve data from adjacent process --
            CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,k,i,&
              MPI_COMM_WORLD,status,ier)

            !-- unpack buffer and write to file --
            iib=1; fib=dim_embbd
            DO j=1,comm_size
              write(unstdy_unit,*) real_buffer(iib:fib)
              iib=iib+dim_embbd; fib=fib+dim_embbd
            END DO
          END DO
          DEALLOCATE(real_buffer)
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO k=1,num_procs-1

          num_comm_iters=num_pelm_pp(k+1)/max_comm_size
          resid_comm_size=mod(num_pelm_pp(k+1),max_comm_size)
          comm_size=max_comm_size
          buffer_size=comm_size*dim_embbd
          ALLOCATE(real_buffer(buffer_size))

          DO i=0,num_comm_iters
            IF (i == num_comm_iters) THEN
              comm_size=resid_comm_size
              buffer_size=comm_size*dim_embbd
            END IF

            !-- recieve data from adjacent process --
            CALL MPI_RECV(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,k,i,&
              MPI_COMM_WORLD,status,ier)

            !-- unpack buffer and write to file --
            iib=1; fib=dim_embbd
            DO j=1,comm_size
              write(unstdy_unit,*) real_buffer(iib:fib), 0.d0
              iib=iib+dim_embbd; fib=fib+dim_embbd
            END DO
          END DO
          DEALLOCATE(real_buffer)
        END DO
      END IF
      !-- clean up --
      DEALLOCATE(status)

      !-------------------------------------------------------------------------
      !-- write primal volume indices --
      !-- write header --
      write(unstdy_unit,'(A,I10,I10)') 'CELLS', glb_num_elm(kmax), (kmax+1)*glb_num_elm(kmax)

      !--
      num_comm_iters=glb_num_elm(kmax)/max_comm_size
      resid_comm_size=mod(glb_num_elm(kmax),max_comm_size)
      comm_size=max_comm_size
      buffer_size=comm_size*kmax
      ALLOCATE(int_buffer(buffer_size),junk_int_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*kmax
        END IF

        int_buffer=0
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          iib=(gindx-1)*kmax+1;
          fib=iib+kmax-1
          IF (lcl_complex(kmax)%orientation(j) == 1) THEN
            int_buffer(iib:fib)=lcl_complex(kmax)%node_indx(j,1:kmax)-1
          ELSE
            int_buffer(iib)=lcl_complex(kmax)%node_indx(j,2)-1
            int_buffer(iib+1)=lcl_complex(kmax)%node_indx(j,1)-1
            int_buffer(iib+2:fib)=lcl_complex(kmax)%node_indx(j,3:kmax)-1
          END IF
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(int_buffer,junk_int_buffer,buffer_size,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)

        iib=1; fib=kmax-1
        DO j=1,comm_size
          WRITE(unstdy_unit,*) kmax, junk_int_buffer(iib:iib+fib)
          iib=iib+kmax
        END DO
      END DO

      !-- clean up --
      DEALLOCATE(int_buffer, junk_int_buffer)


      !-------------------------------------------------------------------------
      !-- write primal volumes types --
      !-- write header --
      write(unstdy_unit,'(A,I10)') 'CELL_TYPES', glb_num_elm(kmax)
      !-- write data --
      SELECT CASE(dim_cmplx)
      CASE(1)
        DO i=1,glb_num_elm(kmax); write(unstdy_unit,'(I1)') 3; END DO
      CASE(2)
        DO i=1,glb_num_elm(kmax); write(unstdy_unit,'(I1)') 5; END DO
      CASE(3)
        DO i=1,glb_num_elm(kmax); write(unstdy_unit,'(I2)') 10; END DO
      CASE DEFAULT
      END SELECT

      !-------------------------------------------------------------------------
      !-- write primal volumes scalar data --
      !-- write header --
      write(unstdy_unit,'(A,I10)') 'CELL_DATA', glb_num_elm(kmax)
      write(unstdy_unit,'(A,A,A)') 'SCALARS ',sclr_name,' double 1'
      write(unstdy_unit,'(A)') 'LOOKUP_TABLE default'

      !--
      comm_size=max_comm_size
      buffer_size=comm_size
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size
        END IF

        real_buffer=-large
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          real_buffer(gindx)=lcl_complex(kmax)% dual_sol(j,1)
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

        DO j=1,comm_size
          WRITE(unstdy_unit,*) junk_real_buffer(j)
        END DO
      END DO

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

      !-------------------------------------------------------------------------
      !-- write cell center vector data --
      !-- write header --
      write(unstdy_unit,'(A,A,A)') 'VECTORS ',vctr_name,' double'

      !--
      comm_size=max_comm_size
      buffer_size=comm_size*dim_embbd
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      IF (dim_embbd==3) THEN
        DO i=0,num_comm_iters
          IF (i == num_comm_iters) THEN
            comm_size=resid_comm_size
            buffer_size=comm_size*dim_embbd
          END IF

          real_buffer=-large
          DO j=strt_indx,num_elm(kmax)
            gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
            IF (gindx>comm_size) THEN
              strt_indx=j
              EXIT
            END IF
            iib=(gindx-1)*dim_embbd+1;
            fib=iib+dim_embbd-1
            real_buffer(iib:fib)=lcl_complex(kmax)%whtny_sol(j,1:dim_embbd)
            IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
          END DO
          CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

          iib=1; fib=dim_embbd-1
          DO j=1,comm_size
            WRITE(unstdy_unit,*) junk_real_buffer(iib:iib+fib)
            iib=iib+dim_embbd
          END DO
        END DO
      ELSEIF (dim_embbd==2) THEN
        DO i=0,num_comm_iters
          IF (i == num_comm_iters) THEN
            comm_size=resid_comm_size
            buffer_size=comm_size*dim_embbd
          END IF

          real_buffer=-large
          DO j=strt_indx,num_elm(kmax)
            gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
            IF (gindx>comm_size) THEN
              strt_indx=j
              EXIT
            END IF
            iib=(gindx-1)*dim_embbd+1;
            fib=iib+dim_embbd-1
            real_buffer(iib:fib)=lcl_complex(kmax)%whtny_sol(j,1:dim_embbd)
            IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
          END DO
          CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

          iib=1; fib=dim_embbd-1
          DO j=1,comm_size
            WRITE(unstdy_unit,*) junk_real_buffer(iib:iib+fib), 0.d0
            iib=iib+dim_embbd
          END DO
        END DO
      END IF

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

      !-------------------------------------------------------------------------
      !-- close files --
      CLOSE(unstdy_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF

    ELSE

      !-- adjacent processes --
      !-------------------------------------------------------------------------
      !-- send header --

      !-------------------------------------------------------------------------
      !-- send primal vertex locations --

      !-- find starting index of local primary vertices --
      DO i=1,num_elm(1)
        IF (lcl_complex(1)% glb_indx(i)==glb_offset+1) THEN
          ptr=i-1; EXIT
        END IF
      END DO

      !-- set up and pack buffer --
      num_comm_iters=num_pelm_pp(rank+1)/max_comm_size
      resid_comm_size=mod(num_pelm_pp(rank+1),max_comm_size)
      comm_size=max_comm_size
      buffer_size=comm_size*dim_embbd
      ALLOCATE(real_buffer(buffer_size))

      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*dim_embbd
        END IF

        iib=1; fib=dim_embbd-1
        DO j=1,comm_size
          real_buffer(iib:iib+fib)=lcl_complex(1)%centers(ptr+j,:)
          iib=iib+dim_embbd
        END DO

        !-- send primal vertex locations --
        CALL MPI_SEND(real_buffer,buffer_size,MPI_DOUBLE_PRECISION,0,i,MPI_COMM_WORLD,ier)
      END DO

      DEALLOCATE(real_buffer)

      !-------------------------------------------------------------------------
      !-- send primal volume indices --

      !--
      num_comm_iters=glb_num_elm(kmax)/max_comm_size
      resid_comm_size=mod(glb_num_elm(kmax),max_comm_size)
      comm_size=max_comm_size
      buffer_size=comm_size*kmax
      ALLOCATE(int_buffer(buffer_size),junk_int_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*kmax
        END IF

        int_buffer=0
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          iib=(gindx-1)*kmax+1;
          fib=iib+kmax-1
          IF (lcl_complex(kmax)%orientation(j) == 1) THEN
            int_buffer(iib:fib)=lcl_complex(kmax)%node_indx(j,1:kmax)-1
          ELSE
            int_buffer(iib)=lcl_complex(kmax)%node_indx(j,2)-1
            int_buffer(iib+1)=lcl_complex(kmax)%node_indx(j,1)-1
            int_buffer(iib+2:fib)=lcl_complex(kmax)%node_indx(j,3:kmax)-1
          END IF
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(int_buffer,junk_int_buffer,buffer_size,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)
      END DO

      !-- clean up --
      DEALLOCATE(int_buffer, junk_int_buffer)

      !-------------------------------------------------------------------------
      !-- send primal volumes types --

      !-------------------------------------------------------------------------
      !-- send primal volume scalar values --

      !--
      comm_size=max_comm_size
      buffer_size=comm_size
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size
        END IF

        real_buffer=-large
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          real_buffer(gindx)=lcl_complex(kmax)% dual_sol(j,1)
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

      END DO

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

      !-------------------------------------------------------------------------
      !-- send cell center vector data --

      !--
      comm_size=max_comm_size
      buffer_size=comm_size*dim_embbd
      ALLOCATE(real_buffer(buffer_size),junk_real_buffer(buffer_size))
      strt_indx=1
      DO i=0,num_comm_iters
        IF (i == num_comm_iters) THEN
          comm_size=resid_comm_size
          buffer_size=comm_size*dim_embbd
        END IF

        real_buffer=-large
        DO j=strt_indx,num_elm(kmax)
          gindx=lcl_complex(kmax)%glb_indx(j) - i*max_comm_size
          IF (gindx>comm_size) THEN
            strt_indx=j
            EXIT
          END IF
          iib=(gindx-1)*dim_embbd+1;
          fib=iib+dim_embbd-1
          real_buffer(iib:fib)=lcl_complex(kmax)%whtny_sol(j,1:dim_embbd)
          IF (j == num_elm(kmax)) strt_indx=num_elm(kmax)+1
        END DO
        CALL MPI_REDUCE(real_buffer,junk_real_buffer,buffer_size,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

      END DO

      !-- clean up --
      DEALLOCATE(real_buffer, junk_real_buffer)

    END IF

    !-- wait for previous communications to finish --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|write_unsteady_D0S2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/read_bndry_cond
!* SYNOPSIS
  SUBROUTINE read_bndry_cond()
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: junk             !> junk IO variable
    INTEGER, ALLOCATABLE  :: node_indices(:)  !> junk IO variable for node indicies

    INTEGER               :: ie,in,i,j        !> loop indices (procs,elems)

    INTEGER               :: buffer_size      !> size of MPI comm buffer
    INTEGER, ALLOCATABLE  :: int_buffer(:,:)  !> integer buffer for MPI comms
    INTEGER, ALLOCATABLE  :: req(:)           !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE  :: status(:)        !> size of MPI comm buffer

    INTEGER, ALLOCATABLE  :: prml_elm(:,:)    !> primal element list
    INTEGER               :: proc, proc_old   !> processor rank (current, past)
    INTEGER               :: ptr              !> pointer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- Open mesh file and check for errors --
    fname=trim(mesh_prefix) // ".face"; CALL root_open_file_read(fname,mesh_unit)

    !-- read header and broadcast global number of entities --
    IF (rank == root) READ(mesh_unit,*)

    !-- allocate boundary type --
    ALLOCATE(lcl_complex(dim_cmplx)%bc_type(num_elm(dim_cmplx)))
    lcl_complex(dim_cmplx)%bc_type=0

    !-- set buffer size and allocate temp. storage for primal elements --
    buffer_size=2
    ALLOCATE(prml_elm(num_elm(dim_cmplx),buffer_size))

    IF (rank == root) THEN
      !-- on root process --

      !-- allocate integer communication buffer --
      ALLOCATE(int_buffer(buffer_size,glb_num_elm(dim_cmplx)+1), &
        node_indices(dim_cmplx))

      !-- loop through primal elements --
      ptr=0;
      DO ie=1,glb_num_elm(dim_cmplx)
        !-- read nodal indices of primal element and sort --
        READ(mesh_unit,*) junk, node_indices,int_buffer(2,ie)
        IF (int_buffer(2,ie) == 0) CYCLE
        int_buffer(1,ie)=ie
        node_indices=node_indices+indx_offset(1)
        CALL int_insertion_sort(node_indices)

        !-- loop through dual elements --
        proc_old=-1
        DO in=1,dim_cmplx
          !-- get process of dual element --
          proc=elm2proc(node_indices(in))
          IF (proc==0 .and. proc/=proc_old) THEN
            !-- if local: assign to local array --
            ptr=ptr+1; prml_elm(ptr,:)=int_buffer(:,ie)
          ELSEIF (proc/=proc_old) THEN
            !-- otherwise: send data to neighbouring process --
            CALL MPI_ISEND(int_buffer(:,ie),buffer_size,MPI_INTEGER,proc,0,&
              MPI_COMM_WORLD,junk,ier)
            CALL MPI_REQUEST_FREE(junk,ier)
          END IF
          proc_old=proc
        END DO
      END DO

      !-- send communication termination message --
      int_buffer(:,ie)=-1
      DO proc=1,num_procs-1
        CALL MPI_ISEND(int_buffer(:,ie),buffer_size,MPI_INTEGER,proc,0,&
          MPI_COMM_WORLD,junk,ier)
        CALL MPI_REQUEST_FREE(junk,ier)
      END DO

    ELSE
      !-- on neighbouring process --

      !-- allocate integer communication buffer and status variable --
      ALLOCATE(int_buffer(buffer_size,1),status(MPI_STATUS_SIZE))

      !-- poll for incoming MPI communications --
      ptr=0
      DO
        !-- receive nodal indices (blocking) --
        CALL MPI_RECV(int_buffer(:,1),buffer_size,MPI_INTEGER,0,0,&
          MPI_COMM_WORLD,status,ier)
        IF (ALL(int_buffer==-1)) THEN
          !-- received termination message: stop polling for communication --
          EXIT
        ELSE
          !-- assign to local array --
          ptr=ptr+1; prml_elm(ptr,:)=int_buffer(:,1)
        END IF
      END DO
    END IF

    !-- assign primal element data to local structures --
    ptr_loop: DO i=1,ptr
      DO j=1,num_elm(dim_cmplx)
        IF (prml_elm(i,1)==lcl_complex(dim_cmplx)%glb_indx(j)) THEN
          lcl_complex(dim_cmplx)%bc_type(j)=prml_elm(i,2)
          CYCLE ptr_loop
        END IF
      END DO
      WRITE(log_unit,'(A,I5)')'Error : face not found=',prml_elm(i,1)
      CALL end_mpi()
    END DO ptr_loop

    !-- close file and clean up --
    IF (rank == root) THEN
      CLOSE(mesh_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF
      DEALLOCATE(prml_elm,int_buffer,node_indices)
    ELSE
      DEALLOCATE(prml_elm,int_buffer,status)
    END IF

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|read_bndry_cond
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/read_bndry_cond2
!* SYNOPSIS
  SUBROUTINE read_bndry_cond2()
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: junk             !> junk IO variable
    INTEGER, ALLOCATABLE  :: node_indices(:)  !> junk IO variable for node indicies

    INTEGER               :: ie,in,i,j        !> loop indices (procs,elems)

    INTEGER               :: buffer_size      !> size of MPI comm buffer
    INTEGER, ALLOCATABLE  :: int_buffer(:)    !> integer buffer for MPI comms
    INTEGER, ALLOCATABLE  :: req(:)           !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE  :: status(:)        !> size of MPI comm buffer

    INTEGER, ALLOCATABLE  :: prml_elm(:,:)    !> primal element list
    INTEGER               :: proc, proc_old   !> processor rank (current, past)
    INTEGER               :: ptr              !> pointer
    INTEGER               :: k, l, orient, iib, fib, num_read_iters, &
      resid_read_size, read_size, stride, strt_indx

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- Open mesh file and check for errors --
    IF (dim_cmplx == 3) then
      fname=trim(mesh_prefix) // ".face"; CALL root_open_file_read(fname,mesh_unit)
    ELSEIF (dim_cmplx == 2) then
      fname=trim(mesh_prefix) // ".edge"; CALL root_open_file_read(fname,mesh_unit)
    END IF
    !-- read header and broadcast global number of entities --
    IF (rank == root) READ(mesh_unit,*)

    !-- allocate boundary type --
    ALLOCATE(lcl_complex(dim_cmplx)%bc_type(num_elm(dim_cmplx)))
    lcl_complex(dim_cmplx)%bc_type=0

    !-- get some basic information for reuse --
    k=dim_cmplx
    read_size=max_read_size
    stride=2; buffer_size=read_size*stride
    num_read_iters=glb_num_elm(k)/read_size
    resid_read_size=mod(glb_num_elm(k),read_size)
    ptr=0; fib=1
    ALLOCATE(int_buffer(read_size*(stride+1)))


    IF (rank == root) THEN
      !-- on root process --

      !-- allocate integer communication buffer --
      ALLOCATE(node_indices(dim_cmplx))

      !--
      DO i=0,num_read_iters
        IF (i == num_read_iters) THEN
          read_size=resid_read_size
          buffer_size=read_size*stride
        END IF
        !--
        iib=1
        DO j=1, read_size
          READ(mesh_unit,*) junk, node_indices, int_buffer(iib+1)
          !--
          int_buffer(iib)=i*max_read_size+j
          iib=iib+stride
        END DO

        !--
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        !--
        iib=1
        DO j=1, read_size
          IF (int_buffer(iib+1)==0) THEN
            iib=iib+stride; CYCLE
          END IF
          DO l=1,num_elm(dim_cmplx)
            IF (int_buffer(iib)==lcl_complex(dim_cmplx)%glb_indx(l)) THEN
              lcl_complex(dim_cmplx)%bc_type(l)=int_buffer(iib+1)
              EXIT
            END IF
          END DO
          iib=iib+stride
        END DO
      END DO

      !-- close file --
      DEALLOCATE(node_indices)
      CLOSE(mesh_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF
    ELSE
      !-- on neighbouring processes --

      !--
      DO i=0,num_read_iters
        IF (i == num_read_iters) THEN
          read_size=resid_read_size
          buffer_size=read_size*stride
        END IF

        !--
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        !--
        iib=1
        DO j=1, read_size
          IF (int_buffer(iib+1)==0)  THEN
            iib=iib+stride; CYCLE
          END IF
          DO l=1,num_elm(dim_cmplx)
            IF (int_buffer(iib)==lcl_complex(dim_cmplx)%glb_indx(l)) THEN
              lcl_complex(dim_cmplx)%bc_type(l)=int_buffer(iib+1)
              EXIT
            END IF
          END DO
          iib=iib+stride
        END DO
      END DO
    END IF

    !-- clean up --
    DEALLOCATE(int_buffer)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|read_bndry_cond2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/read_bndry_cond3
!* SYNOPSIS
  SUBROUTINE read_bndry_cond3()
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: junk             !> junk IO variable
    INTEGER, ALLOCATABLE  :: node_indices(:)  !> junk IO variable for node indicies

    INTEGER               :: ie,in,i,j        !> loop indices (procs,elems)

    INTEGER               :: buffer_size, &
                             r_buffer_size    !> size of MPI comm buffer
    INTEGER, ALLOCATABLE  :: int_buffer(:)    !> integer buffer for MPI comms
    REAL(KIND=PGMSiwp), ALLOCATABLE  :: &
                             real_buffer(:)   !> integer buffer for MPI comms
    INTEGER, ALLOCATABLE  :: req(:)           !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE  :: status(:)        !> size of MPI comm buffer

    INTEGER, ALLOCATABLE  :: prml_elm(:,:)    !> primal element list
    INTEGER               :: proc, proc_old   !> processor rank (current, past)
    INTEGER               :: ptr              !> pointer
    INTEGER               :: k, l, orient, iib, fib, num_read_iters, &
      resid_read_size, read_size, stride, strt_indx

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- Open mesh file and check for errors --
    IF (dim_cmplx == 3) then
      fname=trim(mesh_prefix) // ".face.new"; CALL root_open_file_read(fname,mesh_unit)
    ELSEIF (dim_cmplx == 2) then
      fname=trim(mesh_prefix) // ".edge.new"; CALL root_open_file_read(fname,mesh_unit)
    END IF
    !-- read header and broadcast global number of entities --
    IF (rank == root) READ(mesh_unit,*)

    !-- allocate boundary type --
    ALLOCATE(lcl_complex(dim_cmplx)%bc_type(num_elm(dim_cmplx)), &
      lcl_complex(dim_cmplx)%bc_val(num_elm(dim_cmplx),1))
    lcl_complex(dim_cmplx)%bc_type=0; lcl_complex(dim_cmplx)%bc_val=0.d0

    !-- get some basic information for reuse --
    k=dim_cmplx
    read_size=max_read_size
    stride=2; buffer_size=read_size*stride
    r_buffer_size=read_size
    num_read_iters=glb_num_elm(k)/read_size
    resid_read_size=mod(glb_num_elm(k),read_size)
    ptr=0; fib=1
    ALLOCATE(int_buffer(buffer_size),real_buffer(r_buffer_size))


    IF (rank == root) THEN
      !-- on root process --

      !-- allocate integer communication buffer --
      ALLOCATE(node_indices(dim_cmplx))

      !--
      DO i=0,num_read_iters
        IF (i == num_read_iters) THEN
          read_size=resid_read_size
          buffer_size=read_size*stride
          r_buffer_size=read_size
        END IF
        !--
        iib=1
        DO j=1, read_size
          READ(mesh_unit,*) junk, node_indices, int_buffer(iib+1), real_buffer(j)
          !--
          int_buffer(iib)=i*max_read_size+j
          iib=iib+stride
        END DO

        !--
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST(real_buffer,r_buffer_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

        !--
        iib=1
        DO j=1, read_size
          IF (int_buffer(iib+1)==0) THEN
            iib=iib+stride; CYCLE
          END IF
          DO l=1,num_elm(dim_cmplx)
            IF (int_buffer(iib)==lcl_complex(dim_cmplx)%glb_indx(l)) THEN
              lcl_complex(dim_cmplx)%bc_type(l) =int_buffer(iib+1)
              lcl_complex(dim_cmplx)%bc_val(l,1)=real_buffer(j)
              EXIT
            END IF
          END DO
          iib=iib+stride
        END DO
      END DO

      !-- close file --
      DEALLOCATE(node_indices)
      CLOSE(mesh_unit,IOSTAT=ier)
      IF (ier /= 0) THEN
        WRITE(log_unit,'(A,A,A,I5)')'***Error: failed to close ',trim(fname),&
          ': errcode=',ier
        !err=.TRUE.
      END IF
    ELSE
      !-- on neighbouring processes --

      !--
      DO i=0,num_read_iters
        IF (i == num_read_iters) THEN
          read_size=resid_read_size
          buffer_size=read_size*stride
          r_buffer_size=read_size
        END IF

        !--
        CALL MPI_BCAST(int_buffer,buffer_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST(real_buffer,r_buffer_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

        !--
        iib=1
        DO j=1, read_size
          IF (int_buffer(iib+1)==0)  THEN
            iib=iib+stride; CYCLE
          END IF
          DO l=1,num_elm(dim_cmplx)
            IF (int_buffer(iib)==lcl_complex(dim_cmplx)%glb_indx(l)) THEN
              lcl_complex(dim_cmplx)%bc_type(l)=int_buffer(iib+1)
              lcl_complex(dim_cmplx)%bc_val(l,1)=real_buffer(j)
              EXIT
            END IF
          END DO
          iib=iib+stride
        END DO
      END DO
    END IF

    !-- clean up --
    DEALLOCATE(int_buffer,real_buffer)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|read_bndry_cond3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/init_write_MATLAB
!* SYNOPSIS
  SUBROUTINE init_write_MATLAB()
!* PURPOSE
!*   Setup PETSc for MATLAB output
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/23
!> Setup PETSc for MATLAB output
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- Set MATLAB as output format --
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_centers_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/clean_PETSc_output
!* SYNOPSIS
  SUBROUTINE clean_PETSc_output()
!* PURPOSE
!*   Clean up PETSc output format
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/23
!> Clean up PETSc output format
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- number of rows --
    CALL PetscViewerPopFormat(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|clean_PETSc_output
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_centers_MATLAB2
!* SYNOPSIS
  SUBROUTINE write_centers_MATLAB2()
!* PURPOSE
!*   Compute unsigned volume for given set of points
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Compute unsigned volume for given set of points
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER         :: i,id          !> loop and temporary indicies
    INTEGER         :: m             !> number of rows in global solution variables
    CHARACTER(LEN=slen)   :: fname   !> file name
    Vec             :: work          !> PETSc work array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- number of rows --
    m=glb_num_elm(dim_cmplx+1)+glb_num_elm(dim_cmplx)

    !-- write circumcenters to MATLAB file --
    CALL VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,m,work,petsc_ier)
    fname=trim(sol_prefix) // '_xyz.m'
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    DO id=1,3
      DO i=1,num_elm(dim_cmplx)
        CALL VecSetValues(work,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
          lcl_complex(dim_cmplx)%centers(i,id),INSERT_VALUES,petsc_ier)
      END DO
      DO i=1,num_elm(dim_cmplx+1)
        CALL VecSetValues(work,1,glb_num_elm(dim_cmplx)+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
          lcl_complex(dim_cmplx+1)%centers(i,id),INSERT_VALUES,petsc_ier)
      END DO
      CALL VecAssemblyBegin(work,petsc_ier); CALL VecAssemblyEnd(work,petsc_ier)
      SELECT CASE(id)
      CASE(1)
        CALL PetscObjectSetName(work,'x',petsc_ier)
      CASE(2)
        CALL PetscObjectSetName(work,'y',petsc_ier)
      CASE(3)
        CALL PetscObjectSetName(work,'z',petsc_ier)
      END SELECT
      CALL VecView(work,output_view,petsc_ier)
    END DO
    CALL PetscViewerDestroy(output_view,petsc_ier)
    CALL VecDestroy(work,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_centers_MATLAB2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_prml_volumes_MATLAB2
!* SYNOPSIS
  SUBROUTINE write_prml_volumes_MATLAB2()
!* PURPOSE
!*   Compute unsigned volume for given set of points
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Compute unsigned volume for given set of points
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER         :: i,id          !> loop and temporary indicies
    INTEGER         :: m             !> number of rows in global solution variables
    CHARACTER(LEN=slen)   :: fname   !> file name
    Vec             :: work          !> PETSc work array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- number of rows --
    m=glb_num_elm(dim_cmplx+1)+glb_num_elm(dim_cmplx)

    !-- write volumes to MATLAB file --
    CALL VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,m,work,petsc_ier)
    fname=trim(sol_prefix) // '_vol.m'
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    DO i=1,num_elm(dim_cmplx)
      CALL VecSetValues(work,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%prml_volume(i),INSERT_VALUES,petsc_ier)
    END DO
    DO i=1,num_elm(dim_cmplx+1)
      CALL VecSetValues(work,1,glb_num_elm(dim_cmplx)+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx+1)%prml_volume(i),INSERT_VALUES,petsc_ier)
    END DO
    CALL VecAssemblyBegin(work,petsc_ier); CALL VecAssemblyEnd(work,petsc_ier)
    CALL PetscObjectSetName(work,'vol',petsc_ier)
    CALL VecView(work,output_view,petsc_ier)
    CALL VecDestroy(work,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_prml_volumes_MATLAB2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_solution_MATLAB2
!* SYNOPSIS
  SUBROUTINE write_solution_MATLAB2(iter)
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN), OPTIONAL   :: iter   !> iteration for file name

    !-- local variables --
    INTEGER               :: prelen           !> string pointer
    CHARACTER(LEN=slen)   :: fname            !> file name

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IF (PRESENT(iter)) THEN
      prelen=len_trim(unstdy_prefix)
      fname=trim(unstdy_prefix)//"_00000000.m"
      WRITE(fname(prelen+1:prelen+9),'(I9)') 100000000+iter
      fname(prelen+1:prelen+1)="_"
    ELSE
      fname=trim(sol_prefix) // '.m'
    END IF

    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    CALL PetscObjectSetName(q,'flux',petsc_ier)
    CALL VecView(q,output_view,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_solution_MATLAB2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_pressure_MATLAB2
!* SYNOPSIS
  SUBROUTINE write_pressure_MATLAB2(iter)
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN), OPTIONAL   :: iter   !> iteration for file name

    !-- local variables --
    INTEGER               :: prelen           !> string pointer
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: i,j              !> loop index
    REAL(KIND=PGMSiwp),ALLOCATABLE :: sol_edge(:) !> solution variable for MATLAB output
    REAL(KIND=PGMSiwp),ALLOCATABLE :: cnt(:)      !> counter
    REAL(KIND=PGMSiwp)        :: dxyz             !> extrapolation variable
    REAL(KIND=PGMSiwp)        :: nan              !> nan variable
    INTEGER               :: indx             !> face index
    Vec                   :: work             !> PETSc work array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IF (PRESENT(iter)) THEN
      prelen=len_trim(unstdy_prefix)
      fname=trim(unstdy_prefix)//"_00000000_press.m"
      WRITE(fname(prelen+1:prelen+9),'(I9)') 100000000+iter
      fname(prelen+1:prelen+1)="_"
    ELSE
      fname=trim(sol_prefix) // '.m'
    END IF

    !-- extrapolate surface pressures from volumes --
    !-- interior faces --
    ALLOCATE(cnt(num_elm(dim_cmplx)),sol_edge(num_elm(dim_cmplx)))
    cnt=0.d0; sol_edge=0.d0
    DO i=1,num_elm(dim_cmplx+1)
      DO j=1,lcl_complex(dim_cmplx+1)%num_bndry(i)
        indx=lcl_complex(dim_cmplx+1)%bndry(i)%indx(j)
        dxyz=lcl_complex(dim_cmplx+1)%centers(i,3) - lcl_complex(dim_cmplx)%centers(indx,3)
        sol_edge(indx)=sol_edge(indx)+&
          lcl_complex(dim_cmplx)%prml_volume(indx)*( &
          lcl_complex(dim_cmplx+1)%dual_sol(i,1)+dxyz )
        cnt(indx)=cnt(indx)+1.d0
      END DO
    END DO

    !-- write extrapolated pressure solution to MATLAB file --
    CALL VecDuplicate(q,work,petsc_ier)
    DO i=1,num_elm(dim_cmplx)
      IF (cnt(i) < 1.d0) CYCLE
      CALL VecSetValues(work,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        sol_edge(i)/cnt(i),INSERT_VALUES,petsc_ier)
    END DO
    DO i=1,num_elm(dim_cmplx+1)
      CALL VecSetValues(work,1,glb_num_elm(dim_cmplx)+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx+1)%dual_sol(i,1),INSERT_VALUES,petsc_ier)
    END DO
    DEALLOCATE(cnt,sol_edge)
    CALL VecAssemblyBegin(work,petsc_ier); CALL VecAssemblyEnd(work,petsc_ier)
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    CALL PetscObjectSetName(work,'press',petsc_ier)
    CALL VecView(work,output_view,petsc_ier)
    CALL VecDestroy(work,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_pressure_MATLAB2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_centers_MATLAB
!* SYNOPSIS
  SUBROUTINE write_centers_MATLAB()
!* PURPOSE
!*   Compute unsigned volume for given set of points
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Compute unsigned volume for given set of points
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER         :: i,id          !> loop and temporary indicies
    CHARACTER(LEN=slen)   :: fname   !> file name
    Vec             :: work          !> PETSc work array
    Vec             :: wrk1, wrk2    !> sub arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- duplicate solution vector --
    CALL VecDuplicate(q,work,petsc_ier)

    !-- write circumcenters to MATLAB file --
    fname=trim(sol_prefix) // '_xyz.m'
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    DO id=1,dim_embbd
      CALL VecGetSubVector(work,isg(1),wrk1,petsc_ier)
      CALL VecGetSubVector(work,isg(2),wrk2,petsc_ier)
      DO i=1,num_elm(dim_cmplx)
        CALL VecSetValues(wrk1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
          lcl_complex(dim_cmplx)%centers(i,id),INSERT_VALUES,petsc_ier)
      END DO
      DO i=1,num_elm(dim_cmplx+1)
        CALL VecSetValues(wrk2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
          lcl_complex(dim_cmplx+1)%centers(i,id),INSERT_VALUES,petsc_ier)
      END DO
      CALL VecAssemblyBegin(wrk1,petsc_ier); CALL VecAssemblyEnd(wrk1,petsc_ier)
      CALL VecAssemblyBegin(wrk2,petsc_ier); CALL VecAssemblyEnd(wrk2,petsc_ier)
      CALL VecRestoreSubVector(work,isg(1),wrk1,petsc_ier)
      CALL VecRestoreSubVector(work,isg(2),wrk2,petsc_ier)

      SELECT CASE(id)
      CASE(1)
        CALL PetscObjectSetName(work,'x',petsc_ier)
      CASE(2)
        CALL PetscObjectSetName(work,'y',petsc_ier)
      CASE(3)
        CALL PetscObjectSetName(work,'z',petsc_ier)
      END SELECT
      CALL VecView(work,output_view,petsc_ier)
    END DO
    CALL PetscViewerDestroy(output_view,petsc_ier)
    CALL VecDestroy(work,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_centers_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_prml_volumes_MATLAB
!* SYNOPSIS
  SUBROUTINE write_prml_volumes_MATLAB()
!* PURPOSE
!*   Compute unsigned volume for given set of points
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Compute unsigned volume for given set of points
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER         :: i,id          !> loop and temporary indicies
    CHARACTER(LEN=slen)   :: fname   !> file name
    Vec             :: work          !> PETSc work array
    Vec             :: wrk1, wrk2    !> sub arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- duplicate solution vector --
    CALL VecDuplicate(q,work,petsc_ier)
    CALL VecGetSubVector(work,isg(1),wrk1,petsc_ier)
    CALL VecGetSubVector(work,isg(2),wrk2,petsc_ier)

    !-- write volumes to MATLAB file --
    fname=trim(sol_prefix) // '_vol.m'
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    DO i=1,num_elm(dim_cmplx)
      CALL VecSetValues(wrk1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%prml_volume(i),INSERT_VALUES,petsc_ier)
    END DO
    DO i=1,num_elm(dim_cmplx+1)
      CALL VecSetValues(wrk2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx+1)%prml_volume(i),INSERT_VALUES,petsc_ier)
    END DO
    CALL VecAssemblyBegin(wrk1,petsc_ier); CALL VecAssemblyEnd(wrk1,petsc_ier)
    CALL VecAssemblyBegin(wrk2,petsc_ier); CALL VecAssemblyEnd(wrk2,petsc_ier)
    CALL VecRestoreSubVector(work,isg(1),wrk1,petsc_ier)
    CALL VecRestoreSubVector(work,isg(2),wrk2,petsc_ier)
    CALL PetscObjectSetName(work,'vol',petsc_ier)
    CALL VecView(work,output_view,petsc_ier)
    CALL VecDestroy(work,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_prml_volumes_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_solution_MATLAB
!* SYNOPSIS
  SUBROUTINE write_solution_MATLAB(iter)
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN), OPTIONAL   :: iter   !> iteration for file name

    !-- local variables --
    INTEGER               :: prelen           !> string pointer
    CHARACTER(LEN=slen)   :: fname            !> file name

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IF (PRESENT(iter)) THEN
      prelen=len_trim(unstdy_prefix)
      fname=trim(unstdy_prefix)//"_00000000.m"
      WRITE(fname(prelen+1:prelen+9),'(I9)') 100000000+iter
      fname(prelen+1:prelen+1)="_"
    ELSE
      fname=trim(sol_prefix) // '.m'
    END IF

    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    CALL PetscObjectSetName(q,'flux',petsc_ier)
    CALL VecView(q,output_view,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_solution_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_pressure_MATLAB
!* SYNOPSIS
  SUBROUTINE write_pressure_MATLAB(iter)
!* PURPOSE
!*   Read local node indices of highest order primal elements
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Read local node indices of highest order primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN), OPTIONAL   :: iter   !> iteration for file name

    !-- local variables --
    INTEGER               :: prelen           !> string pointer
    CHARACTER(LEN=slen)   :: fname            !> file name
    INTEGER               :: i,j              !> loop index
    REAL(KIND=PGMSiwp),ALLOCATABLE :: sol_edge(:) !> solution variable for MATLAB output
    REAL(KIND=PGMSiwp),ALLOCATABLE :: cnt(:)      !> counter
    REAL(KIND=PGMSiwp)        :: dxyz             !> extrapolation variable
    REAL(KIND=PGMSiwp)        :: nan              !> nan variable
    INTEGER               :: indx             !> face index
    Vec                   :: work             !> PETSc work array
    Vec                   :: wrk1, wrk2       !> sub arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IF (PRESENT(iter)) THEN
      prelen=len_trim(unstdy_prefix)
      fname=trim(unstdy_prefix)//"_00000000_press.m"
      WRITE(fname(prelen+1:prelen+9),'(I9)') 100000000+iter
      fname(prelen+1:prelen+1)="_"
    ELSE
      fname=trim(sol_prefix) // '_press.m'
    END IF

    !-- extrapolate surface pressures from volumes --
    !-- interior faces --
    ALLOCATE(cnt(num_elm(dim_cmplx)),sol_edge(num_elm(dim_cmplx)))
    cnt=0.d0; sol_edge=0.d0
    DO i=1,num_elm(dim_cmplx+1)
      DO j=1,lcl_complex(dim_cmplx+1)%num_bndry(i)
        indx=lcl_complex(dim_cmplx+1)%bndry(i)%indx(j)
        dxyz=lcl_complex(dim_cmplx+1)%centers(i,dim_embbd) - lcl_complex(dim_cmplx)%centers(indx,dim_embbd)
        sol_edge(indx)=sol_edge(indx)+&
          lcl_complex(dim_cmplx)%prml_volume(indx)*( &
          lcl_complex(dim_cmplx+1)%dual_sol(i,1)+dxyz )
        cnt(indx)=cnt(indx)+1.d0
      END DO
    END DO

    !-- write extrapolated pressure solution to MATLAB file --
    CALL VecDuplicate(q,work,petsc_ier)
    CALL VecGetSubVector(work,isg(1),wrk1,petsc_ier)
    CALL VecGetSubVector(work,isg(2),wrk2,petsc_ier)
    DO i=1,num_elm(dim_cmplx)
      IF (cnt(i) < 1.d0) CYCLE
      CALL VecSetValues(wrk1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        sol_edge(i)/cnt(i),INSERT_VALUES,petsc_ier)
    END DO
    DO i=1,num_elm(dim_cmplx+1)
      CALL VecSetValues(wrk2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx+1)%dual_sol(i,1),INSERT_VALUES,petsc_ier)
    END DO
    DEALLOCATE(cnt,sol_edge)
    CALL VecAssemblyBegin(wrk1,petsc_ier); CALL VecAssemblyEnd(wrk1,petsc_ier)
    CALL VecAssemblyBegin(wrk2,petsc_ier); CALL VecAssemblyEnd(wrk2,petsc_ier)
    CALL VecRestoreSubVector(work,isg(1),wrk1,petsc_ier)
    CALL VecRestoreSubVector(work,isg(2),wrk2,petsc_ier)
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    CALL PetscObjectSetName(work,'press',petsc_ier)
    CALL VecView(work,output_view,petsc_ier)
    CALL VecDestroy(work,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_pressure_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_vec_MATLAB
!* SYNOPSIS
  SUBROUTINE write_vec_MATLAB(Vvar,vname,fname)
!* PURPOSE
!*   write PETSc vec to file in MATLAB format
!* INPUTS
!*   Name                    Description
!*   var                     PETSc vector to be written to file
!*   vname                   variable name for MATLAB output
!*   fname                   file name
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/03/15: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/03/15
!> write PETSc vec to file in MATLAB format
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: &
              fname, &       !> file name
              vname          !> variable name
    Vec, INTENT(IN) :: &
              Vvar           !> PETSc work array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- number of rows --
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    CALL PetscObjectSetName(Vvar,vname,petsc_ier)
    CALL VecView(Vvar,output_view,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_vec_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  io_mod/write_vec_MATLAB
!* SYNOPSIS
  SUBROUTINE write_mat_MATLAB(Mvar,vname,fname)
!* PURPOSE
!*   write PETSc vec to file in MATLAB format
!* INPUTS
!*   Name                    Description
!*   var                     PETSc vector to be written to file
!*   vname                   variable name for MATLAB output
!*   fname                   file name
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2021/03/15: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2021/03/15
!> write PETSc vec to file in MATLAB format
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    CHARACTER(LEN=*), INTENT(IN)  :: &
              fname, &       !> file name
              vname          !> variable name
    Mat, INTENT(IN) :: &
              Mvar           !> PETSc work array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- number of rows --
    CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,fname,output_view,petsc_ier)
    CALL PetscViewerPushFormat(output_view,PETSC_VIEWER_ASCII_MATLAB,petsc_ier)
    !CALL PetscObjectSetName(Mvar,vname,petsc_ier)
    CALL MatView(Mvar,output_view,petsc_ier)
    CALL PetscViewerDestroy(output_view,petsc_ier)

    RETURN

END SUBROUTINE
! io_mod|write_vec_MATLAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! io_mod
!===============================================================================
!
