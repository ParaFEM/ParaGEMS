!
!===============================================================================
!/****v* modules|common/global_io
!* PURPOSE
!*   Define global IO variables
!* CONTAINS
!*   Structure/Variable       Purpose
!*   MPI-IO                   - Use MPI-IO?
!*   file prefixes            - Default file prefexes
!*   IO units                 - Default IO units
!*   IO namelist              - IO namelist for input file
!*   index offsets            - Offsets to account indices not starting from 1
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Define global IO variables
!===============================================================================

  !> # MPI-IO
  !> Use MPI-IO?
  !> ##
  LOGICAL                 ::  use_mpiio = .FALSE. !>- use parallel IO

  !> # file prefexes
  !> Default file prefexes
  !> ##
  CHARACTER(LEN=slen)     ::  &                               !>- file prefixes:
                              input_file    ='input.param', & !>- input
                              mesh_prefix   ='mesh', &        !>- mesh
                              sol_prefix    ='solution', &    !>- solution
                              unstdy_prefix ='unsteady', &    !>- unsteady
                              log_prefix    ='paragems'       !>- log

  !> # IO units
  !> Default IO units
  !> ##
  INTEGER, PARAMETER      ::  &                   !>- IO units:
                              input_unit  = 10, & !>- input
                              mesh_unit   = 20, & !>- mesh files
                              node_unit   = 22, & !>- node files
                              edge_unit   = 24, & !>- edge files
                              face_unit   = 26, & !>- face files
                              volm_unit   = 28, & !>- volm files
                              sol_unit    = 30, & !>- solution
                              unstdy_unit = 40, & !>- unsteady
                              ulog_unit   = 45, & !>- unsteady log
                              log_unit    = 90    !>- logs

  !> # IO max read size
  !> max number of reads per iteration
  !> ##
  INTEGER, PARAMETER      ::  max_read_size=1000  !>- IO units:

  !> # IO namelist
  !> IO namelist for input file
  !> ##
  NAMELIST /io_param/ mesh_prefix, sol_prefix, unstdy_prefix, log_prefix

  !> # index offsets
  !> Offsets to account entity indices not starting from 1
  !> ##
  INTEGER                 ::  indx_offset(4)=0  !>- index offsets from 1 in mesh
                                                !> files
! global_io
!===============================================================================
!
