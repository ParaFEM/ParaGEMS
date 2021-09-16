!
!===============================================================================
!/****v* modules|common/global_precision
!* PURPOSE
!*   Define real precision
!* CONTAINS
!*   Structure/Variable      Purpose
!*   real precision          - Define real precision
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Define real precision
!===============================================================================

  !> # real precision
  !> Define real precision
  !> ##
  INTEGER, PARAMETER          :: iwp = SELECTED_REAL_KIND(15,300)
                                                            !>- real precision

!END MODULE
! global_precision
!===============================================================================
!
