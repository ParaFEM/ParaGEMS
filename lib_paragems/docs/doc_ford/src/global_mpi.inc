!
!===============================================================================
!/****v* modules|common/global_mpi
!* PURPOSE
!*   Define global MPI variables
!* CONTAINS
!*   Structure/Variable      Purpose
!*   general MPI info        - general MPI information
!*   MPI adjacent info       - MPI adjacent process (neighbour) information
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Define global MPI variables
!===============================================================================

  !> # general MPI information
  !> general MPI information
  !> ##
  INTEGER, PARAMETER    :: root = 0         !>- root process id
  INTEGER               :: rank             !>- local process id (indexed from 0)
  INTEGER               :: num_procs        !>- total # processes

  !> # MPI adjacent process (neighbour) information
  !> MPI adjacent process (neighbour) information
  !> ##
  INTEGER               :: num_adj_proc     !>- num adjacent processes
  INTEGER, ALLOCATABLE  :: adj_proc(:)      !>- ids of adjacent processes
  INTEGER, ALLOCATABLE  :: num_send(:,:)    !>- num sends for:
  INTEGER, ALLOCATABLE  :: num_recv(:,:)    !>- num receives for:
                                            !> 1 geometric order k;
                                            !> 2 adjacent process
! global_mpi
!===============================================================================
!
