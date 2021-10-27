!
!===============================================================================
!-- DEC Module
!> Module contains routines for partitioning parallel workload
!===============================================================================
!/****/h* modules|mpi/partition_mod
MODULE partition_mod
!* PURPOSE
!*   Module contains routines for partitioning parallel workload
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 general mpi routines (start, end, syncwrite,etc)
!*   io_mod                  file IO routines
!*   math_mod                basic math functions
!*   dec_mod                 DEC related functions
!* CONTAINS
!*   Subroutine              Purpose
!*   parallel_setup          sets up partitioning, read mesh and nodal
!*                           locations, and sets connectivity
!*   partition_dual          partitions dual volumes across available processes
!*   get_lcl_prml_elms       ??? get local node indices of highest order primal elements
!*   get_connectivity
!*   get_glb_indx_dual_vlm
!*   get_glb_indx
!*   exchange_prml_nodes
!*   exchange_glb_indx
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Module contains routines for partitioning parallel workload
!===============================================================================
! TO DO:
! - MPI+OpenMP?
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod
  USE io_mod
  USE math_mod
  USE dec_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/parallel_setup
!* SYNOPSIS
  SUBROUTINE parallel_setup()
!* PURPOSE
!*   Sets up partitioning, read mesh and nodal locations, and sets connectivity
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Sets up partitioning, read mesh and nodal locations, and sets connectivity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER :: k  !> simplicial order

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- partition dual elements (eq. primal nodes) across available processes
    !   and get local node indices of highest order primal elements --
    CALL setup_partitioning()

    !-- get local connectivity and setup parallel connectivity --
    DO k=dim_cmplx+1,2,-1; CALL calc_bndry_cobndry(k); END DO
    CALL setup_connectivity()

    !-- read primal node locations --
    CALL read_nodes_prml(); CALL exchange_prml_nodes(); RETURN

  END SUBROUTINE
! partition_mod|parallel_setup
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/setup_partitioning
!* SYNOPSIS
  SUBROUTINE setup_partitioning()
!* PURPOSE
!*   Sets up partitioning, read mesh and nodal locations, and sets connectivity
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Sets up partitioning, read mesh and nodal locations, and sets connectivity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER               :: i                      !> counter
    REAL(KIND=PGMSiwp)        :: tmp_time, tmp_time_int !> temporary timing variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- partition dual elements (eq. primal nodes) across available processes
    !   and get local node indices of highest order primal elements --
    CALL partition_dual(); CALL read_prml_elms(); RETURN

  END SUBROUTINE
! partition_mod|setup_partitioning
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/setup_connectivity
!* SYNOPSIS
  SUBROUTINE setup_connectivity()
!* PURPOSE
!*   Sets up partitioning, read mesh and nodal locations, and sets connectivity
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Sets up partitioning, read mesh and nodal locations, and sets connectivity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER :: k  !> simplicial order

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- get parallel connectivity --
    CALL get_glb_indx_dual_vlm(); CALL get_connectivity()
    DO k=2,dim_cmplx; CALL get_glb_indx(k); END DO; RETURN

  END SUBROUTINE
! partition_mod|setup_connectivity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/partition_dual
!* SYNOPSIS
  SUBROUTINE partition_dual()
!* PURPOSE
!*   Partitions dual volumes across available processes
!* ASSUMPTION
!*   Mesh is a simplicial complex in TetGen format
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Partitions dual volumes across available processes
!> Assumption: Mesh is a simplicial complex in TetGen format
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    CHARACTER(LEN=slen) :: fname          !> file name
    INTEGER             :: idx            !> index variable
    INTEGER             :: k              !> simplicial order
    INTEGER             :: int_buffer(6)  !> integer buffer for MPI comms
    INTEGER             :: junk           !> junk IO variable
    LOGICAL             :: err=.FALSE.    !> error variable
    CHARACTER(LEN=slen) :: msg            !> error message string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- get mesh data on root process --
    IF (rank == root) THEN !-- Open primary element file and check for errors --
      fname=trim(mesh_prefix) // ".ele"
      OPEN(mesh_unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
      err=chkerr_log(ier,'OPEN - error opening '//TRIM(fname),log_unit)
      IF (.NOT. err) THEN;
        !-- read header for number of primary elements and complex dimension --
        READ(mesh_unit,*) junk, k; dim_cmplx=k - 1; glb_num_elm(k)=junk
        !-- read first element to get index offset
        READ(mesh_unit,*) junk; indx_offset(k)=1-junk; CLOSE(mesh_unit)
        !-- check for minimum dimension of simplicial complex --
        IF (dim_cmplx<1) err=chkerr_log(-1,'element file '//TRIM(fname)//&
          ' - dimension of simplicial complex less than one',log_unit)
      END IF

      !-- Open primal node (eq. dual element) file and check for errors --
      fname=trim(mesh_prefix) // ".node"
      OPEN(mesh_unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
      err=chkerr_log(ier,'OPEN - error opening '//TRIM(fname),log_unit)
      IF (.NOT. err) THEN;
        !-- read header for number of dual elements and embedding dimension --
        READ(mesh_unit,*) glb_num_elm(1), dim_embbd
        !-- readfirst node to get index offset
        READ(mesh_unit,*) junk; indx_offset(1)=1-junk; CLOSE(mesh_unit)
        !-- Check parallel compatibility --
        IF (glb_num_elm(1) < num_procs)  THEN
          WRITE(msg,'(A,I5,A,I5)') 'number of available processes, ',&
            num_procs, ', greater than the number of dual elements '//&
            '(eq. primary nodes), ',glb_num_elm(1)
          err=chkerr_log(-1,msg,log_unit); END IF
      END IF

      !-- pack MPI buffer --
      int_buffer(1)=k; int_buffer(2)=dim_embbd
      int_buffer(3)=glb_num_elm(1); int_buffer(4)=glb_num_elm(k)
      int_buffer(5)=indx_offset(1); int_buffer(6)=indx_offset(k)
    END IF; CALL chkerrMPI_BCAST(err) !-- complete error checking --

    !-- broadcast data --
    CALL MPI_BCAST(int_buffer,6,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    IF (rank /= root) THEN !-- unpack data broadcast --
      k              =int_buffer(1); dim_embbd      =int_buffer(2)
      glb_num_elm(1) =int_buffer(3); glb_num_elm(k) =int_buffer(4)
      indx_offset(1) =int_buffer(5); indx_offset(k) =int_buffer(6)
      dim_cmplx=k - 1; END IF

    !-- calculate number of local dual elements per process, and the indices of
    !   the dual elements on each process --
    !-- partition dual elements (base number per node) --
    num_elm(1)=glb_num_elm(1) / num_procs

    !-- compute additional dual elements and partition also compute offset for
    !   ranges --
    extra_pelm=mod(glb_num_elm(1),num_procs); IF (rank < extra_pelm) THEN
      num_elm(1)=num_elm(1)+1; glb_offset=rank*num_elm(1)
    ELSE; glb_offset=extra_pelm+rank*num_elm(1); END IF

    !-- gather the number of elements per_node --
    IF (ALLOCATED(num_pelm_pp)) DEALLOCATE(num_pelm_pp)
    ALLOCATE(num_pelm_pp(num_procs)); CALL MPI_ALLGATHER(num_elm(1),1,&
      MPI_INTEGER,num_pelm_pp,1,MPI_INTEGER,MPI_COMM_WORLD,ier)

    !-- allocate local variables --
    ALLOCATE(lcl_complex(k)); ALLOCATE(lcl_complex(1)%glb_indx(num_elm(1)))
    lcl_complex(1)%glb_indx=(/ (idx, idx=1, num_elm(1)) /)+glb_offset
    RETURN

  END SUBROUTINE
! partition_mod|partition_dual
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/get_connectivity
!* SYNOPSIS
  SUBROUTINE get_connectivity()
!* PURPOSE
!*   Determine cross-process entity connections for all geometric orders
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
!> Determine cross-process entity connections for all geometric orders
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER               :: k,km             !> simplicial order
    INTEGER               :: i,j              !> loop indices
    INTEGER               :: junk             !> junk variable
    INTEGER               :: ptr,ptr2,ptr3    !> pointers
    INTEGER, ALLOCATABLE  :: req(:)           !> request var (non-blking comms)
    INTEGER, ALLOCATABLE  :: status(:,:)      !> status var (non-blking comms)
    INTEGER, ALLOCATABLE  :: sbuffer(:)       !> boundary index
    INTEGER, ALLOCATABLE  :: rbuffer(:)       !> boundary index

    CHARACTER(LEN=slen)   :: msg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- send expected number of boundaries from adjacent processes (for all
    !    geometric orders) then send boundary indices --
    !-- get number of adjacent processes --
    num_adj_proc=0; DO j=1,lcl_complex(2)%num_recv
      IF (j==1 .OR. lcl_complex(2)%recv_indx(j,1) /= &
          lcl_complex(2)%recv_indx(max(1,j-1),1)) num_adj_proc=num_adj_proc+1
    END DO; ALLOCATE(adj_proc(num_adj_proc))

    !-- allocate MPI arrays --
    ALLOCATE(num_send(dim_cmplx,num_adj_proc),num_recv(dim_cmplx,num_adj_proc),&
      req(dim_cmplx*num_adj_proc),status(MPI_STATUS_SIZE,dim_cmplx*num_adj_proc))
    num_send=0; num_recv=0

    !-- pack send buffer with number of required number of boundaries for each
    !   geometric order --
    ptr=0; DO j=1,lcl_complex(2)%num_recv
      IF (j==1 .OR. lcl_complex(2)%recv_indx(j,1) /= &
          lcl_complex(2)%recv_indx(max(1,j-1),1)) THEN; ptr=ptr+1
        adj_proc(ptr)=lcl_complex(2)%recv_indx(j,1); num_recv(1,ptr)=1
      ELSE; num_recv(1,ptr)=num_recv(1,ptr)+1; END IF
    END DO

    ptr2=lcl_complex(2)%num_recv
    DO k=3,dim_cmplx+1; km=k-1
      DO j=1,lcl_complex(k)%num_recv
        IF (j==1 .OR. lcl_complex(k)%recv_indx(j,1) /= &
            lcl_complex(k)%recv_indx(max(1,j-1),1)) THEN
          ptr=index_in_list(lcl_complex(k)%recv_indx(j,1),adj_proc,ptr)
          num_recv(km,ptr)=1
        ELSE; num_recv(km,ptr)=num_recv(km,ptr)+1; END IF
      END DO; ptr2=ptr2+lcl_complex(k)%num_recv*(k-1)
    END DO
    ALLOCATE(sbuffer(ptr2))

    !-- send/receive expected number of boundaries from adjacent processes
    !   then send boundary indices (for all geometric orders) --
    DO j=1,num_adj_proc; ptr=(j-1)*dim_cmplx; ptr2=ptr+dim_cmplx
      CALL MPI_ISEND(num_recv(:,j),dim_cmplx,MPI_INTEGER,adj_proc(j),0,&
          MPI_COMM_WORLD,junk,ier); CALL MPI_REQUEST_FREE(junk,ier)
      CALL MPI_IRECV(num_send(:,j),dim_cmplx,MPI_INTEGER,adj_proc(j),0,&
          MPI_COMM_WORLD,req(j),ier)
    END DO

    ptr=1
    DO k=2,dim_cmplx+1; km=k-1; ptr2=ptr
      DO j=1,lcl_complex(k)%num_recv; ptr3=ptr+km
        sbuffer(ptr:ptr3-1)=&
          lcl_complex(km)%node_indx(lcl_complex(k)%recv_indx(j,2),:)
        ptr=ptr3
      END DO

      DO j=1,num_adj_proc
        IF (num_recv(km,j) == 0)  CYCLE
        ptr3=ptr2+km*num_recv(km,j)
        CALL MPI_ISEND(sbuffer(ptr2:ptr3-1),km*num_recv(km,j),MPI_INTEGER,&
          adj_proc(j),k,MPI_COMM_WORLD,junk,ier); CALL MPI_REQUEST_FREE(junk,ier)
        ptr2=ptr3
      END DO
    END DO

    !-- receive receive boundary indices (for all geometric orders) --
    CALL MPI_WAITALL(num_adj_proc,req(1:num_adj_proc),status(:,1:num_adj_proc),&
      ier); req=0; status=0

    ptr=0
    DO k=2,dim_cmplx+1; km=k-1
      lcl_complex(k)%num_send=sum(num_send(km,:))
      ALLOCATE(lcl_complex(k)%send_indx(lcl_complex(k)%num_send,2))
      ptr=ptr+lcl_complex(k)%num_send*km
    END DO
    ALLOCATE(rbuffer(ptr))

    ptr=1 !-- recv buffer pointer --
    DO k=2,dim_cmplx+1; km=k-1; ptr2=1 !-- lcl_complex(k)%recv_indx pointer --
      DO j=1,num_adj_proc
        IF (num_send(km,j) == 0) THEN
          req((k-2)*num_adj_proc+j)=MPI_REQUEST_NULL; CYCLE; END IF
        ptr3=ptr2+num_send(km,j)
        lcl_complex(k)%send_indx(ptr2:ptr3-1,1)=adj_proc(j)
        ptr2=ptr3

        ptr3=ptr+km*num_send(km,j)
        CALL MPI_IRECV(rbuffer(ptr:ptr3-1),km*num_send(km,j),MPI_INTEGER,&
          adj_proc(j),k,MPI_COMM_WORLD,req((k-2)*num_adj_proc+j),ier)
        ptr=ptr3
      END DO
    END DO

    !-- get local indices for send_indx --
    CALL MPI_WAITALL(dim_cmplx*num_adj_proc,req,status,ier)

    ptr=1 !-- send buffer pointer
    DO k=2,dim_cmplx+1
      km=k-1
      DO j=1,lcl_complex(k)%num_send
        ptr2=ptr+km - 1
        DO i=1, num_elm(k-1)
          IF (ALL(rbuffer(ptr:ptr2) == lcl_complex(k-1)%node_indx(i,1:km))) THEN
            lcl_complex(k)%send_indx(j,2)=i
            EXIT
          END IF
        END DO
        ptr=ptr2+1
      END DO
    END DO

    !-- clean up --
    DEALLOCATE(req,status)

    DO k=2,dim_cmplx+1
      CALL MPI_REDUCE(lcl_complex(k)%num_send,i,1,MPI_INTEGER,MPI_SUM,0,&
        MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(lcl_complex(k)%num_recv,j,1,MPI_INTEGER,MPI_SUM,0,&
        MPI_COMM_WORLD,ier)
      write(msg,'(A,I1,A,I8,A,I8)') '   - k: ',k,' :number of sends: ',i,&
        ' :number of receives: ',j
      CALL rootwrite(msg,log_unit)
    END DO

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! partition_mod|get_connectivity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/get_connectivity2
!* SYNOPSIS
  SUBROUTINE get_connectivity2(k)
!* PURPOSE
!*   Determine cross-process entity connections for all geometric orders
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
!> Determine cross-process entity connections for all geometric orders
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER,INTENT(IN)    :: k                !> simplicial order

    !-- local variables --
    INTEGER               :: km               !> simplicial order
    INTEGER               :: i,j              !> loop indices
    INTEGER               :: junk             !> junk variable
    INTEGER               :: ptr,ptr2,ptr3    !> pointers
    INTEGER, ALLOCATABLE  :: req(:)           !> request var (non-blking comms)
    INTEGER, ALLOCATABLE  :: status(:,:)      !> status var (non-blking comms)
    INTEGER, ALLOCATABLE  :: sbuffer(:)       !> boundary index
    INTEGER, ALLOCATABLE  :: rbuffer(:)       !> boundary index

    CHARACTER(LEN=slen)   :: msg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- send expected number of boundaries from adjacent processes (for all
    !    geometric orders) then send boundary indices --
    !-- get number of adjacent processes --
    num_adj_proc=0; DO j=1,lcl_complex(k)%num_recv
      IF (j==1 .OR. lcl_complex(k)%recv_indx(j,1) /= &
          lcl_complex(k)%recv_indx(max(1,j-1),1)) num_adj_proc=num_adj_proc+1
    END DO; ALLOCATE(adj_proc(num_adj_proc))

    !-- allocate MPI arrays --
    ALLOCATE(num_send(dim_cmplx,num_adj_proc),num_recv(dim_cmplx,num_adj_proc))
    num_send=0; num_recv=0

    !-- pack send buffer with number of required number of boundaries for each
    !   geometric order --
    ptr=0; DO j=1,lcl_complex(k)%num_recv
      IF (j==1 .OR. lcl_complex(k)%recv_indx(j,1) /= &
          lcl_complex(k)%recv_indx(max(1,j-1),1)) THEN; ptr=ptr+1
        adj_proc(ptr)=lcl_complex(k)%recv_indx(j,1)
        num_recv(1,ptr)=1
      ELSE
        num_recv(1,ptr)=num_recv(1,ptr)+1; END IF
    END DO

    km=k-1
    DO j=1,lcl_complex(k)%num_recv
      IF (j==1 .OR. lcl_complex(k)%recv_indx(j,1) /= &
          lcl_complex(k)%recv_indx(max(1,j-1),1)) THEN
        ptr=index_in_list(lcl_complex(k)%recv_indx(j,1),adj_proc,ptr)
        num_recv(km,ptr)=1
      ELSE; num_recv(km,ptr)=num_recv(km,ptr)+1; END IF
    END DO

    CALL MPI_REDUCE(lcl_complex(k)%num_send,i,1,MPI_INTEGER,MPI_SUM,0,&
      MPI_COMM_WORLD,ier)
    CALL MPI_REDUCE(lcl_complex(k)%num_recv,j,1,MPI_INTEGER,MPI_SUM,0,&
      MPI_COMM_WORLD,ier)
    write(msg,'(A,I1,A,I8,A,I8)') '   - k: ',k,' :number of sends: ',i,&
      ' :number of receives: ',j
    CALL rootwrite(msg,log_unit)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! partition_mod|get_connectivity2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/get_glb_indx_dual_vlm
!* SYNOPSIS
  SUBROUTINE get_glb_indx_dual_vlm()
!* PURPOSE
!*   Set global index of nodes from predetermined node index
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Set global index of nodes from predetermined node index
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process: --
    IF (ALLOCATED(lcl_complex(1)%glb_indx)) DEALLOCATE(lcl_complex(1)%glb_indx)
    ALLOCATE(lcl_complex(1)%glb_indx(num_elm(1)))
    lcl_complex(1)%glb_indx=lcl_complex(1)%node_indx(:,1); RETURN

  END SUBROUTINE
! partition_mod|get_glb_indx_dual_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/get_glb_indx
!* SYNOPSIS
  SUBROUTINE get_glb_indx(k)
!* PURPOSE
!*   Get global entity indices from file
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Get global entity indices from file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER, INTENT(IN) :: k  !> simplicial order

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- on root process: --
    IF (.NOT.ALLOCATED(lcl_complex(k)%glb_indx)) &
        ALLOCATE(lcl_complex(k)%glb_indx(num_elm(k)))
    CALL read_glb_indx(k); CALL exchange_glb_indx(k); RETURN

  END SUBROUTINE
! partition_mod|get_glb_indx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/exchange_prml_nodes
!* SYNOPSIS
  SUBROUTINE exchange_prml_nodes()
!* PURPOSE
!*   Exchange nodal locations to adjacent process (ghost nodes)
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
!> Exchange nodal locations to adjacent process (ghost nodes)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER               :: ptr                !> pointer
    INTEGER               :: j                  !> loop index

    INTEGER               :: junk               !> junk comm variable
    INTEGER               :: buffer_size        !> size of comm buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE :: sbuffer(:)   !> send buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE :: rbuffer(:)   !> send buffer
    INTEGER, ALLOCATABLE  :: req(:)             !> request var (non-blking comm)
    INTEGER, ALLOCATABLE  :: status(:,:)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- pack send buffer --
    ALLOCATE(sbuffer(lcl_complex(2)%num_send*dim_embbd))
    DO j=1,lcl_complex(2)%num_send
      sbuffer((j-1)*dim_embbd+1:j*dim_embbd)=&
        lcl_complex(1)%centers(lcl_complex(2)%send_indx(j,2),1:dim_embbd)
    END DO

    !-- send data to adjacent processes --
    ptr=1
    DO j=1,num_adj_proc
      IF (num_send(1,j)==0) CYCLE
      buffer_size=num_send(1,j)*dim_embbd
      CALL MPI_ISEND(sbuffer(ptr:ptr+buffer_size-1),buffer_size,&
        MPI_DOUBLE_PRECISION,adj_proc(j),0,MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)
      ptr=ptr+buffer_size
    END DO

    !-- receive data from adjacent processes --
    ALLOCATE(&
      req(num_adj_proc),&
      status(MPI_STATUS_SIZE,num_adj_proc),&
      rbuffer(lcl_complex(2)%num_recv*dim_embbd))
    ptr=1
    DO j=1,num_adj_proc
      IF (num_recv(1,j)==0)  THEN
        req(j)=MPI_REQUEST_NULL
        CYCLE
      END IF
      buffer_size=num_recv(1,j)*dim_embbd
      CALL MPI_IRECV(rbuffer(ptr:ptr+buffer_size-1),buffer_size,&
        MPI_DOUBLE_PRECISION,adj_proc(j),0,MPI_COMM_WORLD,req(j),ier)
      ptr=ptr+buffer_size
    END DO

    !-- unpack receive buffer --
    CALL MPI_WAITALL(num_adj_proc,req,status,ier)
    DO j=1,lcl_complex(2)%num_recv
      lcl_complex(1)%centers(lcl_complex(2)%recv_indx(j,2),1:dim_embbd)=&
        rbuffer((j-1)*dim_embbd+1:j*dim_embbd)
    END DO

    !-- clean up --
    DEALLOCATE(req,status,sbuffer,rbuffer)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|exchange_prml_nodes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* partition_mod/exchange_glb_indx
!* SYNOPSIS
  SUBROUTINE exchange_glb_indx(k)
!* PURPOSE
!*   Exchange global indices to adjacent process (ghost nodes)
!* INPUTS
!*   Name                    Description
!*   k                       simplicial order
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
!> Exchange global indices to adjacent process (ghost nodes)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER,INTENT(IN)    :: k                  !> simplicial order

    !-- local variables --
    INTEGER               :: kp                 !> simplicial order
    INTEGER               :: ptr,ptrp           !> pointer
    INTEGER               :: i,j                !> loop index

    INTEGER               :: junk               !> junk comm variable
    INTEGER               :: buffer_size        !> size of comm buffer
    INTEGER, ALLOCATABLE  :: sbuffer(:)         !> send buffer
    INTEGER, ALLOCATABLE  :: rbuffer(:)         !> receive buffer
    INTEGER, ALLOCATABLE  :: req(:)             !> request var (non-blking comm)
    INTEGER, ALLOCATABLE  :: status(:,:)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- pack send buffer --
    kp=k+1
    ALLOCATE(sbuffer(lcl_complex(kp)%num_send))
    DO j=1,lcl_complex(kp)%num_send
      sbuffer(j)=lcl_complex(k)%glb_indx(lcl_complex(kp)%send_indx(j,2))
    END DO

    !-- send data to adjacent processes --
    ptr=1
    DO j=1,num_adj_proc
      IF (num_send(k,j)==0) CYCLE
      buffer_size=num_send(k,j)
      CALL MPI_ISEND(sbuffer(ptr:ptr+buffer_size-1),buffer_size,&
        MPI_INTEGER,adj_proc(j),0,MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)
      ptr=ptr+buffer_size
    END DO

    !-- receive data from adjacent processes --
    ALLOCATE(&
      req(num_adj_proc),&
      status(MPI_STATUS_SIZE,num_adj_proc),&
      rbuffer(lcl_complex(kp)%num_recv))
    ptr=1
    DO j=1,num_adj_proc
      IF (num_recv(k,j)==0) THEN
        req(j)=MPI_REQUEST_NULL
        CYCLE
      END IF
      buffer_size=num_recv(k,j)
      CALL MPI_IRECV(rbuffer(ptr:ptr+buffer_size-1),buffer_size,&
        MPI_INTEGER,adj_proc(j),0,MPI_COMM_WORLD,req(j),ier)
      ptr=ptr+buffer_size
    END DO

    !-- unpack receive buffer --
    CALL MPI_WAITALL(num_adj_proc,req,status,ier)
    DO j=1,lcl_complex(kp)%num_recv
      lcl_complex(k)%glb_indx(lcl_complex(kp)%recv_indx(j,2))=rbuffer(j)
    END DO

    !-- clean up --
    DEALLOCATE(req,status)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|exchange_glb_indx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!


END MODULE
! partition_mod
!===============================================================================
!
