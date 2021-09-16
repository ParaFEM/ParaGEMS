!
!===============================================================================
!/****v* modules|time_marching/time_marching
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


  !> Global time marching variables
  !> ##
  REAL(KIND=iwp)            :: time = 0.d0          !>- time (simulation, not walltime)
  REAL(KIND=iwp)            :: dt = 0.d0            !>- unsteady step size
  INTEGER                   :: initial_step = 1, &
                               final_step = 1
  INTEGER                   :: step=0, &
                               stage=0

  !> # Index structure
  !> Structure for defining (co)boundary operator relations for an entity
  TYPE time_marching_struct
    LOGICAL                  :: singly = .FALSE.     !> single diagonal coefficient
    LOGICAL                  :: stifflyAcc = .FALSE. !> explicit first stage
    INTEGER                  :: implicit             !> -1 = Rosenbrock/W-method (not implemented)
                                                        !> 0 = explicit
                                                        !> 1 = diagonally-implicit
                                                        !> 2 = full implicit (not implemented)
    INTEGER                  :: p,q, &               !>- order and stage order
                                r,s,s_initial        !>- number of steps and stages
    INTEGER, ALLOCATABLE     :: q_step_indx(:)
    REAL(KIND=iwp), ALLOCATABLE :: &
                                A(:,:), &            !>- Coefficient matrix (internal function)
                                B(:,:), &            !>- Coefficient matrix (internal solution)
                                U(:,:), &            !>- Coefficient matrix (external function)
                                V(:,:), &            !>- Coefficient matrix (external solution)
                                c(:), dc(:), &       !>- abscissae and stage change in abscissae
                                LHS_tvec(:)          !>- LHS time term
          !>- abscissae and stage change in abscissae
  END TYPE
  TYPE(time_marching_struct) :: time_marching

  CHARACTER(len = slen)      :: time_marching_method = 'ie'


  !-- time marching namelist --
  NAMELIST /time_param/ time_marching_method, initial_step, final_step, dt, time

! global_element
!===============================================================================
!
