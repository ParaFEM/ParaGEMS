!
!===============================================================================
!/****v* modules|common/global_io
!* PURPOSE
!*   Define global IO variables
!* CONTAINS
!*   Structure/Variable      Purpose
!*   MPI-IO                  Use MPI-IO?
!*   file prefixes           Default file prefexes
!*   IO units                Default IO units
!*   IO namelist             IO namelist for input file
!*   index offsets           Offsets to account entity indices not starting from 1
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
  LOGICAL             :: use_mpiio = .FALSE.             !>- use parallel IO

  !> # file prefexes
  !> Default file prefexes
  !> ##
  CHARACTER(LEN=slen) :: input_file    ='input.param'    !>- file prefix: input
  CHARACTER(LEN=slen) :: mesh_prefix   ='mesh'           !>- file prefix: mesh
  CHARACTER(LEN=slen) :: sol_prefix    ='solution'       !>- file prefix: solution
  CHARACTER(LEN=slen) :: unstdy_prefix ='unsteady'       !>- file prefix: unsteady
  CHARACTER(LEN=slen) :: log_prefix    ='paragems'       !>- file prefix: log

  !> # IO units
  !> Default IO units
  !> ##
  INTEGER, PARAMETER :: input_unit  = 10                !>- IO unit: input
  INTEGER, PARAMETER :: mesh_unit   = 20                !>- IO unit: mesh files
  INTEGER, PARAMETER :: node_unit   = 22                !>- IO unit: node files
  INTEGER, PARAMETER :: edge_unit   = 24                !>- IO unit: edge files
  INTEGER, PARAMETER :: face_unit   = 26                !>- IO unit: face files
  INTEGER, PARAMETER :: volm_unit   = 28                !>- IO unit: volm files
  INTEGER, PARAMETER :: sol_unit    = 30                !>- IO unit: solution
  INTEGER, PARAMETER :: unstdy_unit = 40                !>- IO unit: unsteady
  INTEGER, PARAMETER :: ulog_unit   = 45                !>- IO unit: unsteady log
  INTEGER, PARAMETER :: log_unit    = 90                !>- IO unit: logs

  !> # IO namelist
  !> IO namelist for input file
  !> ##
  NAMELIST /io_param/ mesh_prefix, sol_prefix, unstdy_prefix, log_prefix

  !> # index offsets
  !> Offsets to account entity indices not starting from 1
  !> ##
  INTEGER               :: indx_offset(4)=0              !>- index offsets in files
! global_io
!===============================================================================
!
