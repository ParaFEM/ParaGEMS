!
!===============================================================================
!/****v* modules|common/global_vars
!* PURPOSE
!*   Define general global variables and parameters
!* CONTAINS
!*   Structure/Variable      Purpose
!*   constants               - general physical constants
!*   reference values        - general reference values
!*   solution                - number of solution variables (prml,dual)
!*   timing                  - global timing variables
!*   error checking          - error variable
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Define general global variables and parameters
!===============================================================================

  !> # constants
  !> general physical constants
  !> ##
  REAL(KIND=PGMSiwp), PARAMETER :: &
                              pi = 3.14159265358979d0   !>- pi

  !> # reference values
  !> general reference values
  !> ##
  INTEGER, PARAMETER      ::  slen = 128    !>- default string length
  REAL(KIND=PGMSiwp), PARAMETER :: &
                              large = 1.d+15, & !>- large value (dp)
                              larges = 1.d+7, & !>- large value (sp)
                              largeh = 1.d+3, & !>- large value (hp)
                              small = 1.d-15, & !>- small value (dp)
                              smalls = 1.d-7, & !>- small value (sp)
                              smallh = 1.d-3    !>- small value (hp)

  !> # solution variables
  !> number of solution variables (prml,dual)
  !> ##
  INTEGER                 ::  num_vars(2,4) !>- num prml/dual solution variables
                                            !> for each geometric order

  !> # time dependent variables
  !> ##
  LOGICAL                 ::  time_dependent = .FALSE.  !>- is time dependent?

  !> # timing variables
  !> global timing variables
  !> ##
  REAL(KIND=PGMSiwp)          ::  start_time    !>- program start time
  REAL(KIND=PGMSiwp)          ::  curnt_time    !>- current program time

  !> # error checking variables
  !> error variable
  !> ##
  INTEGER                 ::  ier           !>- error flag

! global_vars
!===============================================================================
!
