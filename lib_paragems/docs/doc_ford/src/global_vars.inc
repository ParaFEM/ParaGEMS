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
  REAL(KIND=iwp), PARAMETER   :: &
                          pi = 3.14159265358979d0     !>- pi

  !> # reference values
  !> general reference values
  !> ##
  INTEGER, PARAMETER          :: slen = 128           !>- default string length
  REAL(KIND=iwp), PARAMETER   :: large = 1.d+15       !>- large value (dp)
  REAL(KIND=iwp), PARAMETER   :: larges = 1.d+8       !>- large value (sp)
  REAL(KIND=iwp), PARAMETER   :: largeh = 1.d+4       !>- large value (hp)
  REAL(KIND=iwp), PARAMETER   :: small = 1.d-15       !>- small value (dp)
  REAL(KIND=iwp), PARAMETER   :: smalls = 1.d-8       !>- small value (sp)
  REAL(KIND=iwp), PARAMETER   :: smallh = 1.d-4       !>- small value (hp)

  !> # solution variables
  !> number of solution variables (prml,dual)
  !> ##
  INTEGER                     :: num_vars(2,4)        !>- num prml/dual solution
                                                      !>- variables for each
                                                      !> geometric order

  !> # timing variables
  !> global timing variables
  !> ##
  REAL(KIND=iwp)              :: start_time           !>- program start time
  REAL(KIND=iwp)              :: curnt_time           !>- current time

  !> # error checking variables
  !> error variable
  !> ##
  INTEGER                     :: ier                  !>- error flag

! global_vars
!===============================================================================
!
