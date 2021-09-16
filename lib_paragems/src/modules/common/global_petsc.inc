!
!===============================================================================
!/****v* modules|common/global_petsc
!* PURPOSE
!*   Define global variables and parameters for PETSc
!* CONTAINS
!*   Structure/Variable       Purpose
!*   PETSc error variable     - PTEScerror variable
!*   PETSc solution variables - matrices and vectors for solving sys. of eqn's
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/07/08: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/07/08
!> Define global variables and parameters for PETSc
!===============================================================================

  !> # PETSc error variable
  !> ##
  PetscErrorCode          ::  petsc_ier     !>- PETSc error variable

  !> # PETSc solution variables
  !> ##
  Mat                     ::  A, Asub(4), & !>- LHS & sub-matrices
                              Sp            !>- Schur complement
  Vec                     ::  b, &          !>- Boundary conds. / source term
                              r, &          !>- negative spatial residual
                              q, dq, &      !>- solution & update
                              q_scl, scl    !>- scaled solution & scaling vec
  Vec, ALLOCATABLE        ::  q_step(:), &  !>- Unsteady previous solution
                              r_stage(:)    !>- Unsteady (neg) spatial residuals
  IS                      ::  isg(2)        !>- Index set of LHS

! global_vars
!===============================================================================
!
