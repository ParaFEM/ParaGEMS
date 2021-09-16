!
!===============================================================================
!/****v* modules|common/global_element
!* PURPOSE
!*   Define global primary structures
!* CONTAINS
!*   Structure/Variable      Purpose
!*   index structure         - Structure for defining (co)boundary operators
!*   local complex structure - Primary structure for simplicial complex and solution
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
  !> Structure for defining (co)boundary operators
  TYPE indx_stuct
    INTEGER, ALLOCATABLE  :: sgn(:)             !>- sign
    INTEGER, ALLOCATABLE  :: indx(:)            !>- index
  END TYPE

  !> # Local complex structure
  !> Primary structure for simplicial complex and solution
  TYPE lcl_complex_struc
    !> ## global information
    INTEGER, ALLOCATABLE  :: glb_indx(:)        !>- global index

    !> ## local information
    INTEGER, ALLOCATABLE  :: node_indx(:,:)     !>- global node indices
    INTEGER, ALLOCATABLE  :: lcl_node_indx(:,:) !>- local node indices
    INTEGER, ALLOCATABLE  :: orientation(:)     !>- orientation
    INTEGER, ALLOCATABLE  :: deflt_order(:,:)   !>- default ordering

    !> ## boundary information
    INTEGER, ALLOCATABLE  :: bc_type(:)         !>- boundary type

    !> ## DEC information
    REAL(KIND=iwp), ALLOCATABLE :: centers(:,:)       !>- location of centers
    REAL(KIND=iwp), ALLOCATABLE :: prml_volume(:)     !>- primal volume
    REAL(KIND=iwp), ALLOCATABLE :: dual_volume(:)     !>- dual volume
    REAL(KIND=iwp), ALLOCATABLE :: hdg_star(:)        !>- hodge star
    REAL(KIND=iwp), ALLOCATABLE :: inv_hdg_star(:)    !>- hodge star inverse
    INTEGER, ALLOCATABLE        :: num_bndry(:)       !>- num boundaries
    INTEGER, ALLOCATABLE        :: num_cobndry(:)     !>- num co-boundaries
    TYPE(indx_stuct), ALLOCATABLE  :: &
                                   bndry(:)           !>- boundary indices
    TYPE(indx_stuct), ALLOCATABLE  :: &
                                   cobndry(:)         !>- co-boundary indices
    REAL(KIND=iwp), ALLOCATABLE :: prml_dir(:,:)      !>- prml unit direction
    REAL(KIND=iwp), ALLOCATABLE :: dual_dir(:,:)      !>- dual unit direction

    !> ## solution information
    REAL(KIND=iwp), ALLOCATABLE :: prml_sol(:,:)      !>- primal solution
    REAL(KIND=iwp), ALLOCATABLE :: dual_sol(:,:)      !>- dual solution
    REAL(KIND=iwp), ALLOCATABLE :: whtny_sol(:,:)     !>- Whitney interpolated solution

    !> ## external (process non-local) and MPI comms information
    INTEGER                     :: num_ext            !>- num external boundaries
    INTEGER, ALLOCATABLE        :: ext_indx(:,:)      !>- external info:
                                                      !> 1 sign;
                                                      !> 2 external boundary index;
                                                      !> 3 local entity index
    INTEGER                     :: num_send           !>- num boundaries to send
    INTEGER, ALLOCATABLE        :: send_indx(:,:)     !>- send info:
                                                      !> 1 send process rank;
                                                      !> 2 local boundary index
    INTEGER                     :: num_recv           !>- num boundaries to recv
    INTEGER, ALLOCATABLE        :: recv_indx(:,:)     !>- recieve info:
                                                      !> 1 recieve process rank;
                                                      !> 2 local boundary index

    !> ## surface information
    INTEGER                     :: num_surf=0         !>- num surface boundaries
    INTEGER, ALLOCATABLE        :: surf_indx(:,:)     !>- surface info:
                                                      !> 1 sign;
                                                      !> 2 surface boundary index;
                                                      !> 3 local entity index;
                                                      !> 4 boundary type


    !-- crack information --
    INTEGER                     :: num_crack=0        ! # cracked elements
    INTEGER, ALLOCATABLE        :: crack_indx(:)      ! index
  END TYPE
  TYPE(lcl_complex_struc), ALLOCATABLE :: lcl_complex(:)

! global_element
!===============================================================================
!
