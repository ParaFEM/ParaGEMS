!
!===============================================================================
!-- Solver Module
!> Module contains general routines for solvers (linear and nonlinear)
!===============================================================================
!/****/h* modules|solvers/solver_mod
!* SYNOPSIS
MODULE solver_mod
!* PURPOSE
!*   Module contains general routines for solvers (linear and nonlinear)
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 ???
!*   ieee_arithmeti
!*   iso_fortran_env
!*   petsc.h
!*   solver_vars.inc
!* CONTAINS
!*   Subroutine              Purpose
!*
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!*   2019/11/21: KSP (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Module contains general routines for solvers (linear and nonlinear)
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod
  USE, INTRINSIC :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
  USE, INTRINSIC :: iso_fortran_env, only: real32

  IMPLICIT NONE

  !-- include solver variables --
#include <petsc/finclude/petsc.h>
#include "solver_vars.inc"

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/set_defaults_solver
 !* SYNOPSIS
  SUBROUTINE set_defaults_solver(solver)
!* PURPOSE
!*   Set default values for solver
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!*   2019/11/21: KSP (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Set default values for solver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  CHARACTER(LEN=*)  :: solver

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SELECT CASE(solver)
  CASE('KSP')
      pc_type = PCASM
  CASE DEFAULT
    WRITE(log_unit,'(A,A)')'***Error: invalid solver: name: ',solver
    CALL end_mpi()
  END SELECT

  RETURN

  END SUBROUTINE
! solver_mod|set_defaults_solver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/check_param_solver
 !* SYNOPSIS
  SUBROUTINE check_param_solver()
!* PURPOSE
!*   Check values set for solver
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Check values set for solver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  LOGICAL       :: error = .FALSE.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- Nothing checked (yet) --
  IF (error) CALL end_mpi()

  RETURN

  END SUBROUTINE
! solver_mod|check_param_solver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/start_petsc
 !* SYNOPSIS
  SUBROUTINE start_petsc()
!* PURPOSE
!*   Initialise PETSc library
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Initialise PETSc library
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> start_petsc()')

  !-- initialise PETSc and check for errors --
  CALL PetscInitialize(PETSC_NULL_CHARACTER,petsc_ier)
  IF (petsc_ier /= 0) THEN
    WRITE(log_unit,'(A,A,I5)')'***Error: failed to start PETSc: errcode = ',petsc_ier
    CALL end_mpi()
  END IF

  RETURN

  END SUBROUTINE
! solver_mod|start_petsc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/end_petsc
 !* SYNOPSIS
  SUBROUTINE end_petsc()
!* PURPOSE
!*   Finalise PETSc library
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Finalise PETSc library
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- finalise PETSc --
  CALL PetscFinalize(petsc_ier)

  RETURN

  END SUBROUTINE
! solver_mod|end_petsc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/solve_KSP
 !* SYNOPSIS
  SUBROUTINE solve_KSP()
!* PURPOSE
!*   Solve linear system using PETSc KSP
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Solve linear system using PETSc KSP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  CHARACTER(LEN=slen)   :: msg              !> message sting for log file
  INTEGER               :: iters            !> number of iterations
  REAL(KIND=iwp)        :: rnorm            !> final residual norm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> solve_KSP()')
  curnt_time = MPI_WTIME()

  !-- create KSP solver --
  CALL KSPCreate(MPI_COMM_WORLD,ksp_id,petsc_ier)
  CALL KSPSetType(ksp_id,ksp_type,petsc_ier)
  CALL KSPSetFromOptions(ksp_id,petsc_ier)
  CALL KSPSetTolerances(ksp_id,rel_tol,abs_tol,div_tol,max_iter,petsc_ier)
  IF (nz_init) CALL KSPSetInitialGuessNonzero(ksp_id,PETSC_TRUE,petsc_ier)

  !-- setup solver monitor --
  IF (solver_monitor) THEN
    CALL PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,output_vfrmt,petsc_ier)
    CALL KSPMonitorSet(ksp_id,KSPMonitorDefault,output_vfrmt,PetscViewerAndFormatDestroy,petsc_ier)
  END IF

  !-- set LHS --
  CALL KSPSetOperators(ksp_id,A,A,petsc_ier)

  !-- set preconditioner --
  CALL KSPGetPC(ksp_id,pc_id,petsc_ier)
  CALL PCFactorSetMatOrderingType(pc_id,pc_reorder,petsc_ier)
  CALL PCSetType(pc_id,pc_type,petsc_ier)

  !-- solve linear system and output details --
  CALL KSPSetUp(ksp_id,petsc_ier)
  CALL KSPSolve(ksp_id,b,sol,petsc_ier)
  CALL KSPView(ksp_id,PETSC_VIEWER_STDOUT_WORLD,petsc_ier)

  !-- write simulation details to log file --
  CALL KSPGetIterationNumber(ksp_id,iters,petsc_ier)
  WRITE(msg,*) '  - solution iterations: ', iters; CALL syncwrite_log(msg)
  CALL KSPGetResidualNorm(ksp_id,rnorm,petsc_ier)
  WRITE(msg,*) '  - final residual norm: ', rnorm; CALL syncwrite_log(msg)

  RETURN

  END SUBROUTINE
! solver_mod|end_petsc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/solve_KSP2
 !* SYNOPSIS
  SUBROUTINE solve_KSP2()
!* PURPOSE
!*   Solve linear system using PETSc KSP
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Solve linear system using PETSc KSP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  CHARACTER(LEN=slen)   :: msg              !> message sting for log file
  INTEGER               :: iters,m,n        !> number of iterations
  REAL(KIND=iwp)        :: rnorm            !> final residual norm

  CHARACTER(LEN=slen)   :: fname            !> file name
  KSPConvergedReason :: reason              !>
  KSP             :: subksp_id(2),sub2ksp_id(2)
  PC              :: subpc_id(2),sub2pc_id(2)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> solve_KSP()')
  curnt_time = MPI_WTIME()

  !-- create KSP solver --
  CALL KSPCreate(MPI_COMM_WORLD,ksp_id,petsc_ier)
  CALL KSPSetType(ksp_id,ksp_type,petsc_ier)
  CALL KSPSetFromOptions(ksp_id,petsc_ier)
  CALL KSPSetTolerances(ksp_id,rel_tol,abs_tol,div_tol,max_iter,petsc_ier)
  IF (nz_init) CALL KSPSetInitialGuessNonzero(ksp_id,PETSC_TRUE,petsc_ier)

  !-- setup solver monitor --
  IF (solver_monitor) THEN
    CALL PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,output_vfrmt,petsc_ier)
    CALL KSPMonitorSet(ksp_id,KSPMonitorDefault,output_vfrmt,PetscViewerAndFormatDestroy,petsc_ier)
  END IF

  !-- set LHS --
  CALL KSPSetOperators(ksp_id,A,A,petsc_ier)

  !-- set preconditioner --
  CALL KSPGetPC(ksp_id,pc_id,petsc_ier)
  CALL PCSetType(pc_id,PCFIELDSPLIT,petsc_ier)
  CALL PCFieldSplitSetIS(pc_id,"fcs",isg(1),petsc_ier)
  CALL PCFieldSplitSetIS(pc_id,"vls",isg(2),petsc_ier)
  CALL PCFieldSplitSetType(pc_id,PC_COMPOSITE_SCHUR,petsc_ier)
  ! - PC_COMPOSITE_ADDITIVE - Jacobi
  ! - PC_COMPOSITE_MULTIPLICATIVE (default) - Gauss-Seidel
  ! - PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE - Symmetric Gauss-Seidel
  ! - PC_COMPOSITE_SPECIAL - ???
  ! - PC_COMPOSITE_SCHUR - schur
  CALL PCFieldSplitSetSchurPre(pc_id,PC_FIELDSPLIT_SCHUR_PRE_SELFP,PETSC_NULL_Mat,petsc_ier)
  ! - custom - user
  ! - fuse a11 - a11 (a11=0)
  ! - full schur compliment - full
  ! - symbolic schur compliment - self
  ! - approximate schur compliment - selfp
  !!CALL MatSchurComplementSetAinvType(Mat S,MatSchurComplementAinvType ainvtype,petsc_ier)
  ! - for selfp
  ! -  MAT_SCHUR_COMPLEMENT_AINV_DIAG
  ! -  MAT_SCHUR_COMPLEMENT_AINV_LUMP
  ! -  MAT_SCHUR_COMPLEMENT_AINV_BLOCK_DIAG
  CALL PCFieldSplitSetSchurFactType(pc_id,PC_FIELDSPLIT_SCHUR_FACT_FULL,petsc_ier)
  ! - diag,lower,upper,full
  !!CALL PCFieldSplitSetSchurScale(pc_id,1.d0,petsc_ier)
  ! - default is -1.0 (for diag)

  ! CALL PCSetUseAmat(pc_id,PETSC_TRUE,petsc_ier)
! CALL PCFieldSplitSetDiagUseAmat(pc_id,PETSC_TRUE,petsc_ier)
! CALL PCFieldSplitSetOffDiagUseAmat(pc_id,PETSC_TRUE,petsc_ier)
  CALL PCSetUp(pc_id,petsc_ier)
  CALL PCFieldSplitGetSubKSP(pc_id,n,subksp_id,petsc_ier)

  CALL KSPSetType(subksp_id(2),KSPGMRES,petsc_ier);  CALL KSPGetPC(subksp_id(2),subpc_id(2),petsc_ier)
  CALL PCSetType(subpc_id(2),PCHYPRE,petsc_ier)
  !CALL PCHYPRESetType(subpc_id(2),"euclid",petsc_ier)
  CALL PCHYPRESetType(subpc_id(2),"parasails",petsc_ier)
  !CALL KSPSetTolerances(subksp_id(2),1.d-4,1.d-4,PETSC_DEFAULT_REAL,1000,petsc_ier)
  !CALL KSPMonitorSet(subksp_id(2),KSPMonitorDefault,output_vfrmt,PetscViewerAndFormatDestroy,petsc_ier)

  !-- solve linear system and output details --
  CALL KSPSetUp(ksp_id,petsc_ier)
  CALL KSPSolve(ksp_id,b,sol,petsc_ier)
  CALL KSPView(ksp_id,PETSC_VIEWER_STDOUT_WORLD,petsc_ier)

  !-- write simulation details to log file --
  CALL KSPGetIterationNumber(ksp_id,iters,petsc_ier)
  WRITE(msg,*) '  - solution iterations: ', iters; CALL syncwrite_log(msg)
  CALL KSPGetTotalIterations(subksp_id(2),iters,petsc_ier)
  WRITE(msg,*) '  - solution sub iters.: ', iters; CALL syncwrite_log(msg)
  CALL KSPGetResidualNorm(ksp_id,rnorm,petsc_ier)
  WRITE(msg,*) '  - final residual norm: ', rnorm; CALL syncwrite_log(msg)
  CALL KSPGetConvergedReason(ksp_id,reason,petsc_ier)
  WRITE(msg,*) '  - convergence reason : ',reason; CALL syncwrite_log(msg)

  CALL KSPDestroy(ksp_id,petsc_ier)

  RETURN

  END SUBROUTINE
! solver_mod|solve_KSP2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/extract_sol_KSP
 !* SYNOPSIS
  SUBROUTINE extract_sol_KSP()
!* PURPOSE
!*   extract local solution values from distributed PETSc vector
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> extract local solution values from distributed PETSc vector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER               :: i                !> loop index
  Vec                   :: lcl_sol          !> local solution vector
  VecScatter            :: scatter          !> scatter variable
  IS                    :: from, to         !> global and local vecor indices

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !-- get local solution values - primal faces --
  CALL VecCreateSeq(PETSC_COMM_SELF,num_elm(dim_cmplx),lcl_sol,petsc_ier)
  CALL ISCreateGeneral(MPI_COMM_WORLD,num_elm(dim_cmplx),&
    lcl_complex(dim_cmplx)%glb_indx-1,PETSC_COPY_VALUES,from,petsc_ier)
  CALL ISCreateGeneral(PETSC_COMM_SELF,num_elm(dim_cmplx),&
    (/(i,i=0,num_elm(dim_cmplx)-1)/),PETSC_COPY_VALUES,to,petsc_ier)
  CALL VecScatterCreate(sol,from,lcl_sol,to,scatter,petsc_ier)
  CALL VecScatterBegin(scatter,sol,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecScatterEnd(scatter,sol,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecGetValues(lcl_sol,num_elm(dim_cmplx),(/(i,i=0,num_elm(dim_cmplx)-1)/),&
    lcl_complex(dim_cmplx)%prml_sol(:,1),petsc_ier)
  CALL ISDestroy(from,petsc_ier);   CALL ISDestroy(to,petsc_ier)
  CALL VecScatterDestroy(scatter,petsc_ier)
  CALL VecDestroy(lcl_sol,petsc_ier)

  !-- get local solution values - dual vertices --
  CALL VecCreateSeq(PETSC_COMM_SELF,num_elm(dim_cmplx+1),lcl_sol,petsc_ier)
  CALL ISCreateGeneral(MPI_COMM_WORLD,num_elm(dim_cmplx+1),glb_num_elm(dim_cmplx)+&
    lcl_complex(dim_cmplx+1)%glb_indx-1,PETSC_COPY_VALUES,from,petsc_ier)
  CALL ISCreateGeneral(PETSC_COMM_SELF,num_elm(dim_cmplx+1),&
    (/(i,i=0,num_elm(dim_cmplx+1)-1)/),PETSC_COPY_VALUES,to,petsc_ier)
  CALL VecScatterCreate(sol,from,lcl_sol,to,scatter,petsc_ier)
  CALL VecScatterBegin(scatter,sol,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecScatterEnd(scatter,sol,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecGetValues(lcl_sol,num_elm(dim_cmplx+1),(/(i,i=0,num_elm(dim_cmplx+1)-1)/),&
    lcl_complex(dim_cmplx+1)%dual_sol(:,1),petsc_ier)
  CALL ISDestroy(from,petsc_ier);   CALL ISDestroy(to,petsc_ier)
  CALL VecScatterDestroy(scatter,petsc_ier)
  CALL VecDestroy(lcl_sol,petsc_ier)

  !-- log time to solve system --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! solver_mod|extract_sol_KSP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !/****/s* solver_mod/extract_sol_KSP2
 !* SYNOPSIS
  SUBROUTINE extract_sol_KSP2()
!* PURPOSE
!*   extract local solution values from distributed PETSc vector
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> extract local solution values from distributed PETSc vector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER               :: i                !> loop index
  Vec                   :: lcl_sol          !> local solution vector
  VecScatter            :: scatter          !> scatter variable
  IS                    :: from, to         !> global and local vecor indices
  Vec                   :: sol1, sol2       !>

  CHARACTER(LEN=slen)   :: fname          !> file name

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  CALL VecGetSubVector(sol,isg(1),sol1,petsc_ier);  CALL VecGetSubVector(sol,isg(2),sol2,petsc_ier)

  !-- get local solution values - primal faces --
  CALL VecCreateSeq(PETSC_COMM_SELF,num_elm(dim_cmplx),lcl_sol,petsc_ier)
  CALL ISCreateGeneral(MPI_COMM_WORLD,num_elm(dim_cmplx),&
    lcl_complex(dim_cmplx)%glb_indx-1,PETSC_COPY_VALUES,from,petsc_ier)
  CALL ISCreateGeneral(PETSC_COMM_SELF,num_elm(dim_cmplx),&
    (/(i,i=0,num_elm(dim_cmplx)-1)/),PETSC_COPY_VALUES,to,petsc_ier)
  CALL VecScatterCreate(sol1,from,lcl_sol,to,scatter,petsc_ier)
  CALL VecScatterBegin(scatter,sol1,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecScatterEnd(scatter,sol1,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecGetValues(lcl_sol,num_elm(dim_cmplx),(/(i,i=0,num_elm(dim_cmplx)-1)/),&
    lcl_complex(dim_cmplx)%prml_sol(:,1),petsc_ier)
  CALL ISDestroy(from,petsc_ier);   CALL ISDestroy(to,petsc_ier)
  CALL VecScatterDestroy(scatter,petsc_ier)
  CALL VecDestroy(lcl_sol,petsc_ier)

  !-- get local solution values - dual vertices --
  CALL VecCreateSeq(PETSC_COMM_SELF,num_elm(dim_cmplx+1),lcl_sol,petsc_ier)
  CALL ISCreateGeneral(MPI_COMM_WORLD,num_elm(dim_cmplx+1),&
    lcl_complex(dim_cmplx+1)%glb_indx-1,PETSC_COPY_VALUES,from,petsc_ier)
  CALL ISCreateGeneral(PETSC_COMM_SELF,num_elm(dim_cmplx+1),&
    (/(i,i=0,num_elm(dim_cmplx+1)-1)/),PETSC_COPY_VALUES,to,petsc_ier)
  CALL VecScatterCreate(sol2,from,lcl_sol,to,scatter,petsc_ier)
  CALL VecScatterBegin(scatter,sol2,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecScatterEnd(scatter,sol2,lcl_sol,INSERT_VALUES,SCATTER_FORWARD,petsc_ier)
  CALL VecGetValues(lcl_sol,num_elm(dim_cmplx+1),(/(i,i=0,num_elm(dim_cmplx+1)-1)/),&
    lcl_complex(dim_cmplx+1)%dual_sol(:,1),petsc_ier)
  CALL ISDestroy(from,petsc_ier);   CALL ISDestroy(to,petsc_ier)
  CALL VecScatterDestroy(scatter,petsc_ier)
  CALL VecDestroy(lcl_sol,petsc_ier)

  !-- log time to solve system --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! solver_mod|extract_sol_KSP2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!


END MODULE
! solver_mod
!===============================================================================
!
