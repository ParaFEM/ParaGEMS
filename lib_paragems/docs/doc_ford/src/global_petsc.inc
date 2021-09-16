!
!===============================================================================
!/****v* modules|common/global_petsc
!* PURPOSE
!*   Define global variables and parameters for PETSc
!* CONTAINS
!*   Structure/Variable      Purpose
!*   PETSc error variable
!*   PETSc solution variables
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
  PetscErrorCode  :: petsc_ier          !>- PETSc error variable

  !> # PETSc solution variables
  !> ##
  Mat             :: A, Asub(4), Sp     !>- Matries: LHS, sub-matrices, Schur complement
  Vec             :: b, sol             !>- Vectors: RHS and solution
  IS              :: isg(2)             !>- Index set: LHS w/ schur complement

! global_vars
!===============================================================================
!
