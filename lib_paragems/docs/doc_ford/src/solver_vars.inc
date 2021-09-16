!
!===============================================================================
!/****v* modules|solvers/solver_vars
!* PURPOSE
!*   Define global variables for solvers (linear and nonlinear)
!* CONTAINS
!*   Structure/Variable      Purpose
!*   PETSc error variable
!*   PETSc solution variables
!*   PETSc KSP variables
!*   PETSc viewer variables
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Define global variables for solvers (linear and nonlinear)
!===============================================================================

  !-----------------------------------------------------------------------------
  ! variables
  !-----------------------------------------------------------------------------
  !-- KSP variables --
  KSP             :: ksp_id                             ! linear solver id
  KSPType         :: ksp_type = KSPGMRES			          ! PETSC linear solver
  PetscReal       :: rel_tol = PETSC_DEFAULT_REAL, &	  ! Relative termination tolerance
                     abs_tol = PETSC_DEFAULT_REAL, &	  ! Absolute termination tolerance
                     div_tol = PETSC_DEFAULT_REAL		    ! divergence tolerance
  PetscInt        :: max_iter = PETSC_DEFAULT_REAL		  ! Maximum major iterations

  !-- use zero initial guess --
  LOGICAL        :: nz_init = .FALSE.                   ! use non-zero initial soluiton

  !-- Preconditioner variables --
  PC              :: pc_id                              ! preconditioner id
  PCType          :: pc_type = PCASM                    ! PETSC preconditioner (disabled)
  MatOrderingType :: pc_reorder = MATORDERINGRCM        ! PETSC reordering (disabled)

  !-- PETSc viewer variables --
  PetscViewer           :: output_view
  PetscViewerAndFormat  :: output_vfrmt

  !-- PETSc monitor output --
  LOGICAL         :: solver_monitor = .FALSE.           ! output convergence data to command line?

  !-- Solution output --
  LOGICAL         :: sol_output = .FALSE.               ! Output solution data to file?
  INTEGER         :: output_frqcy = 1                   ! Number of crack iterations between writing solution data to file

  !-----------------------------------------------------------------------------
  ! variables
  !-----------------------------------------------------------------------------
  NAMELIST /solver_param/ksp_type, rel_tol, abs_tol, div_tol, max_iter, &
    pc_type, pc_reorder, nz_init, solver_monitor, sol_output, output_frqcy

! solver_vars
!===============================================================================
!
