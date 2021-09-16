!
!===============================================================================
!/****v* modules|common/global_mesh
!* PURPOSE
!*   Define global mesh variables
!* CONTAINS
!*   Structure/Variable      Purpose
!*   general mesh info       - dimensions of complex and embedding space
!*   global mesh info        - global mesh and partitioning variables
!*   local mesh info         - number of local mesh entities
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Define global mesh variables
!===============================================================================

  !> # general mesh information
  !> dimensions of complex and embedding space
  !> ##
  INTEGER               :: dim_cmplx          !>- dimension of complex
  INTEGER               :: dim_embbd          !>- dimension of embedding

  !> # global mesh information
  !> global mesh and partitioning variables
  !> ##
  INTEGER               :: glb_num_elm(4)=0   !>- global # of elements (p/e/f/v)
  INTEGER, ALLOCATABLE  :: num_pelm_pp(:)     !>- num partitioning elements per node
  INTEGER               :: glb_offset=0       !>- global offset of partitioning elements
  INTEGER               :: extra_pelm         !>- num extra partitioning elements
                                              !>- after initial basic partitioning

  !> # local mesh information
  !> number of local mesh entities
  !> ##
  INTEGER               :: num_elm(4)=0       !>- num local elements (p/e/f/v)

! global_mesh
!===============================================================================
!