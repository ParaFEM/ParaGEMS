!
!===============================================================================
!/****v* modules|common/global_element
!* PURPOSE
!*   Define global primary structures
!* CONTAINS
!*   Structure/Variable       Purpose
!*   index struct             - Structure for defining (co)boundary operators
!*   local complex struct     - Primary structure for simplicial complex and solution
!*        > global information
!*        > local information
!*        > DEC information
!*        > solution information
!*        > external and MPI comms information
!*        > surface information
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Define global primary structures
!===============================================================================

  !> # Index structure
  !> Structure for defining (co)boundary operator relations for an entity
  TYPE indx_struct
    INTEGER, ALLOCATABLE  ::  sgn(:), &     !>- sign
                              indx(:)       !>- index
  END TYPE

  !> # Local complex structure
  !> Primary structure for simplicial complex and solution
  TYPE lcl_complex_struct
    !> ## global information
    INTEGER, ALLOCATABLE  ::  glb_indx(:)   !>- global index of entity

    !> ## local information
    INTEGER, ALLOCATABLE  ::  node_indx(:,:), &     !>- global node indices
                              lcl_node_indx(:,:), & !>- local node indices
                              orientation(:)        !>- orientation

    !> ## boundary information
    INTEGER, ALLOCATABLE  ::  bc_type(:)    !>- boundary type of entity
    REAL(KIND=PGMSiwp), ALLOCATABLE :: &
                              bc_val(:,:)   !>- boundary value

    !> ## Geometric information
    REAL(KIND=PGMSiwp), ALLOCATABLE :: &
                              centers(:,:), &   !>- location of centers
                              prml_volume(:), & !>- primal volume
                              dual_volume(:), & !>- dual volume
                              prml_dir(:,:), &  !>- prml unit direction
                              dual_dir(:,:)    !>- dual unit direction

    !> ## DEC information
    INTEGER, ALLOCATABLE  ::  num_bndry(:), &   !>- num boundaries
                              num_cobndry(:)    !>- num co-boundaries
    TYPE(indx_struct), ALLOCATABLE  :: &
                              bndry(:), &       !>- boundary structure
                              cobndry(:)        !>- co-boundary structure
    REAL(KIND=PGMSiwp), ALLOCATABLE :: &
                              hdg_star(:), &    !>- hodge star
                              inv_hdg_star(:)   !>- hodge star inverse

    !> ## solution information
    REAL(KIND=PGMSiwp), ALLOCATABLE :: &
                              prml_sol(:,:), &  !>- primal solution
                              dual_sol(:,:), &  !>- dual solution
                              whtny_sol(:,:)    !>- Whitney interp. solution

    !> ## external (process non-local) and MPI comms information
    INTEGER               ::  num_send, &       !>- num boundaries to send
                              num_recv          !>- num boundaries to recv
    INTEGER, ALLOCATABLE  ::  send_indx(:,:), & !>- send info &
                              recv_indx(:,:)    !>- recieve info:
                                                !> 1 send/recv process rank
                                                !> 2 local boundary index

  END TYPE
  TYPE(lcl_complex_struct), ALLOCATABLE :: lcl_complex(:)

! global_element
!===============================================================================
!
