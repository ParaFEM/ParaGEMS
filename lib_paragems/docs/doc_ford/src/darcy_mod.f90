!
!===============================================================================
!-- Darcy Module
!> Module contains routines specifically related to Darcy flow equations
!===============================================================================
!/****/h* modules|phy_darcy_flow/darcy_mod
MODULE darcy_mod
!* PURPOSE
!*   Module contains routines specifically related to Darcy flow equations
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 general mpi routines (start, end, syncwrite,etc)
!*   io_mod                  IO functions and routines
!*   solver_mod              Solver routines (PETSc)
!*   petsc.h                 PETSc variables and routines
!*   darcy_vars.inc          Darcy specific variables
!* CONTAINS
!*   Subroutine              Purpose
!*   read_input_darcy        read inputs, and check parameters (phys/io/solver)
!*   check_param_darcy       check parameters for Darcy flow simulations
!*   initialise_darcy
!*   finalise_darcy
!*   get_RHS_darcy
!*   get_LHS_darcy
!*   identify_crack
!*   identify_crack2
!*   identify_crack3
!*   exchange_bndry_cond
!*   initialise_darcy2
!*   finalise_darcy2
!*   get_RHS_darcy2
!*   get_LHS_darcy2
!*   identify_crack4
!*   identify_crack5
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Module contains routines specifically related to Darcy flow equations
!===============================================================================
! TO DO:
! - set SOLVER defaults
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod
  USE io_mod
  USE solver_mod

  IMPLICIT NONE

  !-- include Darcy flow specific variables --
#include <petsc/finclude/petsc.h>
#include "darcy_vars.inc"

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/read_input_darcy
!* SYNOPSIS
  SUBROUTINE read_input_darcy()
!* PURPOSE
!*   Read input file for user defined parameters, and perform checks to ensure
!*   reasonable values are set
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Read input filefor user defined parameters, and perform checks to ensure
!> reasonable values are set
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ??? read from root and distribute ???
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  CHARACTER(LEN=slen)   :: fname            !> file name
  CHARACTER(LEN=slen)   :: buffer           !> buffer for command line arguments
  CHARACTER(LEN=SLEN)   :: msg              !> error msg string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> read_input_darcy()')

  !-- check for user supplied input file name, otherwise use default --
  CALL GETARG(1,buffer)
  IF (LEN_TRIM(buffer) == 0) THEN
    fname = input_file
  ELSE
    READ(buffer,*) fname
  END IF

  !-- open input file --
  CALL rootwrite_log('   - input file name: '//fname)
  OPEN(input_unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(msg,'(A,A,A,I5)')'*** Error: failed to open ',trim(fname),&
      ': errcode = ',ier
    CALL rootwrite_log(msg);   CALL end_mpi()
  END IF

  !-- read input file --
  read(input_unit,nml=darcy_param)
  read(input_unit,nml=io_param)
  read(input_unit,nml=solver_param)

  !-- close input file --
  CLOSE(input_unit)

  !-- check parameters --
  CALL check_param_darcy()
  !CALL check_param_io()     !-- nothing checked (yet)
  !CALL check_param_solver() !-- nothing checked (yet)

  RETURN

  END SUBROUTINE
! darcy_mod/read_input_darcy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/check_param_darcy
!* SYNOPSIS
  SUBROUTINE check_param_darcy()
!* PURPOSE
!*   Check values set for solving Darcy flow simulation
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Check values set for solving Darcy flow simulation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  LOGICAL              :: err = .FALSE.    !> error flag
  CHARACTER(LEN=SLEN)  :: msg              !> error msg string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- physical parameters --
  IF (mu<0.d0) THEN
    WRITE(msg,'(A,E10.3)')'***Error: Invalid viscosity mu (must be greater ', &
      'than 0): value given ',mu
    CALL rootwrite_log(msg);   err = .TRUE.
  END IF
  IF (k<0.d0) THEN
    WRITE(log_unit,'(A,E10.3)')'***Error: Invalid permeability k (must be ', &
      'greater than 0): value given ',k
    CALL rootwrite_log(msg);   err = .TRUE.
  END IF
  IF (re<0.d0) THEN
    WRITE(log_unit,'(A,E10.3)')'***Error: Invalid Reynolds number re (must ', &
      'be greater than 0): value given ',re
    CALL rootwrite_log(msg);   err = .TRUE.
  END IF

  IF (err) CALL end_mpi()

  RETURN

  END SUBROUTINE
! darcy_mod/check_param_darcy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/initialise_darcy
!* SYNOPSIS
  SUBROUTINE initialise_darcy()
!* PURPOSE
!*   initialise matrices and vectors for Darcy simulations
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
!> initialise matrices and vectors for Darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER         :: m              !> number of rows/columns in global solution variables
  INTEGER         :: i              !> loop index
  INTEGER         :: iib,fib        !> loop index
  INTEGER         :: offset         !> loop index
  REAL(KIND=iwp), ALLOCATABLE :: o_nnz(:)  !>

  CHARACTER(LEN=slen)   :: fname    !> file name
  Vec             :: work           !> PETSc work array
  INTEGER         :: id             !> index

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> initialise_darcy()')
  curnt_time = MPI_WTIME()

  !-- allocate local solution variables --
  ALLOCATE(&
    lcl_complex(dim_cmplx)%prml_sol(num_elm(dim_cmplx),1),&
    lcl_complex(dim_cmplx+1)%whtny_sol(num_elm(dim_cmplx+1),dim_embbd),&
    lcl_complex(dim_cmplx+1)%dual_sol(num_elm(dim_cmplx+1),1))

  !-- allocate global solution variables --
  !-- number of rows/columns --
  m = glb_num_elm(dim_cmplx+1) + glb_num_elm(dim_cmplx)

  !-- vectors (RHS:b, solution:sol) --
  CALL VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,m,b,petsc_ier)
  CALL VecSetFromOptions(b,petsc_ier)
  CALL VecDuplicate(b,sol,petsc_ier)
  CALL VecDuplicate(b,work,petsc_ier)

  !-- get local range of vectors/matrix --
  CALL VecGetOwnershipRange(work,iib,fib,petsc_ier);   ALLOCATE(o_nnz(fib-iib))

  !-- system matrix (LHS:A) --
  CALL MatCreate(MPI_COMM_WORLD,A,petsc_ier)
  CALL MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m,petsc_ier)
  CALL MatSetFromOptions(A,petsc_ier)

  !-- preallocate system matrix --
  !-- calculate number of non-zeros per row for the equations for Darcy's law --
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%num_cobndry(i)>0) &
      CALL VecSetValues(work,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1.d0*lcl_complex(dim_cmplx)%num_cobndry(i),INSERT_VALUES,petsc_ier)
  END DO

  !-- calculate number of non-zeros per row for the continuity equations --
  offset = glb_num_elm(dim_cmplx)
  DO i = 1,num_elm(dim_cmplx+1)
    CALL VecSetValues(work,1,offset+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
      4.d0,INSERT_VALUES,petsc_ier)
  END DO
  !-- assemble vector --
  CALL VecAssemblyBegin(work,petsc_ier);   CALL VecAssemblyEnd(work,petsc_ier)
  CALL VecGetValues(work,fib-iib,(/(i,i=iib,fib-1)/),o_nnz,petsc_ier)

  !-- preallocate system matrix --
  IF (num_procs==1) THEN
    CALL MatSeqAIJSetPreallocation(A,1,int(o_nnz+1.d0),petsc_ier)
  ELSEIF (iib<glb_num_elm(dim_cmplx) .AND. fib>=glb_num_elm(dim_cmplx)) THEN
    CALL MatMPIAIJSetPreallocation(A,1,int(min(o_nnz+1.d0,real(fib-iib))),1,int(o_nnz+1.d0),petsc_ier)
  ELSE
    CALL MatMPIAIJSetPreallocation(A,1,PETSC_NULL_INTEGER,1,int(o_nnz),petsc_ier)
  END IF

  !-- clean up --
  CALL VecDestroy(work,petsc_ier)
  DEALLOCATE(o_nnz)

  !-- log time to initialise the darcy flow problem --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/initialise_darcy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/finalise_darcy
!* SYNOPSIS
  SUBROUTINE finalise_darcy()
!* PURPOSE
!*   clean up matrices and vectors for Darcy simulations
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
!> clean up matrices and vectors for Darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !--deallocate vectors, matrices, and solver --
  !-- vectors --
  CALL VecDestroy(b,petsc_ier)
  CALL VecDestroy(sol,petsc_ier)

  !-- matrices --
  CALL MatDestroy(A,petsc_ier)

  !-- linear solver --
  CALL KSPDestroy(ksp_id,petsc_ier)

  RETURN

  END SUBROUTINE
! darcy_mod/finalise_darcy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/get_RHS_darcy
!* SYNOPSIS
  SUBROUTINE get_RHS_darcy()
!* PURPOSE
!*   setup right hand side vector and boundary conditions for darcy simulations
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> setup right hand side vector and boundary conditions for darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER :: i,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> get_RHS_darcy()');  curnt_time = MPI_WTIME()

  CALL syncwrite_log('> get_RHS_darcy(): read_bndry_cond()');
  CALL read_bndry_cond2()
  !CALL exchange_bndry_cond(dim_cmplx)
  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  !--  set known values for boundary pressures and fluxes, as well as initial condition --
  num_flx = 0
  ALLOCATE(flx_indx(num_elm(dim_cmplx)))
  DO i = 1,num_elm(dim_cmplx)
    type_bc = lcl_complex(dim_cmplx)%bc_type(i)
    IF (type_bc == -1) THEN
      !-- impermeable boundary --
      CALL VecSetValues(sol,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      num_flx = num_flx + 1
      flx_indx(num_flx) = lcl_complex(dim_cmplx)%glb_indx(i) - 1
    ELSEIF (type_bc == -2) THEN
      !-- flux boundary --
      CALL VecSetValues(sol,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-sum(vel*&
        lcl_complex(dim_cmplx)%dual_dir(i,1:3))*lcl_complex(dim_cmplx)%prml_volume(i),&
        INSERT_VALUES,petsc_ier)
      num_flx = num_flx + 1
      flx_indx(num_flx) = lcl_complex(dim_cmplx)%glb_indx(i) - 1
    ELSEIF (type_bc>=0 .AND. nz_init) THEN
      CALL VecSetValues(sol,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-sum(vel*&
        lcl_complex(dim_cmplx)%dual_dir(i,1:3))*lcl_complex(dim_cmplx)%prml_volume(i),&
        INSERT_VALUES,petsc_ier)
    ! ELSE
    !   WRITE(log_unit,'(A,I5)')'Error : undefined boundary type = ',type_bc
    !   CALL end_mpi()
    END IF
  END DO

  IF (nz_init)  THEN
    DO i = 1,num_elm(dim_cmplx+1)
      CALL VecSetValues(sol,1,glb_num_elm(dim_cmplx)+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        -sum(vel*lcl_complex(dim_cmplx+1)%centers(i,1:3)) + p_ref,INSERT_VALUES,petsc_ier)
      END DO
  END IF
  !-- assemble vector of initial condition --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)

  !-- set known values for boundary pressures and fluxes --
  CALL syncwrite_log('> get_RHS_darcy(): assign bndry_cond()');
  l_use_pref = .TRUE.
  DO i = 1,lcl_complex(dim_cmplx+1)%num_surf
    indx = lcl_complex(dim_cmplx+1)%surf_indx(i,2)
    type_bc = lcl_complex(dim_cmplx+1)%surf_indx(i,4)
    IF (type_bc>0) THEN
      !-- pressure boundary condition --
      CALL VecSetValues(b,1,lcl_complex(dim_cmplx)%glb_indx(indx)-1,&
        -lcl_complex(dim_cmplx+1)%surf_indx(i,1)*bc_press(type_bc),&
        INSERT_VALUES,petsc_ier)
      l_use_pref = .FALSE.
    END IF
  END DO

  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  !-- set reference pressure if no boundary pressure specified --
  CALL syncwrite_log('> get_RHS_darcy(): pressure bndry_cond() and assembly');
  CALL MPI_REDUCE(l_use_pref,use_pref,1,MPI_LOGICAL,MPI_LAND,0,MPI_COMM_WORLD,petsc_ier)
  IF (rank==root .and. .not. nz_init .and. use_pref) &
    CALL VecSetValues(sol,1,glb_num_elm(dim_cmplx),p_ref,INSERT_VALUES,petsc_ier)

  !-- assemble vector of known values --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
  CALL VecAssemblyBegin(b,petsc_ier);   CALL VecAssemblyEnd(b,petsc_ier)
  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  CALL syncwrite_log('> get_RHS_darcy(): update LHS');
  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  IF (rank==root) THEN
    IF (use_pref) THEN
      CALL MatZeroRowsColumns(A,num_flx+1,(/ flx_indx(1:num_flx),&
      glb_num_elm(dim_cmplx) /),1.d0,sol,b,petsc_ier)
    ELSE
      CALL MatZeroRowsColumns(A,num_flx,(/ flx_indx(1:num_flx)/),1.d0,sol,b,petsc_ier)
    END IF
  ELSE
    CALL MatZeroRowsColumns(A,num_flx,(/ flx_indx(1:num_flx)/),1.d0,sol,b,petsc_ier)
  END IF
  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  !-- clean up --
  DEALLOCATE(flx_indx)

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/get_RHS_darcy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/get_LHS_darcy
!* SYNOPSIS
  SUBROUTINE get_LHS_darcy()
!* PURPOSE
!*   setup left hand side matrix for darcy simulations
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> setup left hand side matrix for darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER :: i,j,offset     !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> get_LHS_darcy()')
  curnt_time = MPI_WTIME()

  !-- build LHS (A matrix) --
  !-- equations for Darcy's law --
  offset = glb_num_elm(dim_cmplx)
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%num_cobndry(i)>0) &
      CALL MatSetValues(A,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i)+1,(/ lcl_complex(dim_cmplx)%glb_indx(i),&
        offset+lcl_complex(dim_cmplx+1)%glb_indx(lcl_complex(dim_cmplx)%cobndry(i)%indx) /)-1,&
        (/ -mu/k*lcl_complex(dim_cmplx)%hdg_star(i),real(lcl_complex(dim_cmplx)%cobndry(i)%sgn,iwp) /),&
        INSERT_VALUES,petsc_ier)
    ! CALL MatSetValues(A,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
    !   lcl_complex(dim_cmplx)%num_cobndry(i)+1,(/ lcl_complex(dim_cmplx)%glb_indx(i),&
    !   offset+lcl_complex(dim_cmplx+1)%glb_indx(lcl_complex(dim_cmplx)%cobndry(i)%indx) /)-1,&
    !   (/ -mu/k,real(lcl_complex(dim_cmplx)%cobndry(i)%sgn,iwp) /)*lcl_complex(dim_cmplx)%inv_hdg_star(i),&
    !   INSERT_VALUES,petsc_ier)
  END DO

  !-- continuity equations: dual edges corresponding to internal primal faces --
  DO i = 1,num_elm(dim_cmplx+1)
    CALL MatSetValues(A,1,offset+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
      lcl_complex(dim_cmplx+1)%num_bndry(i)+1,(/ offset+lcl_complex(dim_cmplx+1)%glb_indx(i),&
      lcl_complex(dim_cmplx)%glb_indx(lcl_complex(dim_cmplx+1)%bndry(i)%indx) /)-1,&
      (/ 0.d0, real(lcl_complex(dim_cmplx+1)%bndry(i)%sgn,iwp) /),INSERT_VALUES,petsc_ier)
    ! CALL MatSetValues(A,1,offset+lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
    !   lcl_complex(dim_cmplx+1)%num_bndry(i)+1,(/ offset+lcl_complex(dim_cmplx+1)%glb_indx(i),&
    !   lcl_complex(dim_cmplx)%glb_indx(lcl_complex(dim_cmplx+1)%bndry(i)%indx) /)-1,&
    !   (/ 0.d0, real(lcl_complex(dim_cmplx+1)%bndry(i)%sgn,iwp)*&
    !   lcl_complex(dim_cmplx)%inv_hdg_star(lcl_complex(dim_cmplx+1)%bndry(i)%indx) /),&
    !   INSERT_VALUES,petsc_ier)
  END DO

  !-- continuity equations: dual edges corresponding to surface primal faces --
  DO i = 1,lcl_complex(dim_cmplx+1)%num_surf
    CALL MatSetValues(A,&
      1,offset+lcl_complex(dim_cmplx+1)%glb_indx(lcl_complex(dim_cmplx+1)%surf_indx(i,3))-1,&
      1,lcl_complex(dim_cmplx)%glb_indx(lcl_complex(dim_cmplx+1)%surf_indx(i,2))-1,&
      -real(lcl_complex(dim_cmplx+1)%surf_indx(i,1),iwp),INSERT_VALUES,petsc_ier)
    ! CALL MatSetValues(A,&
    !   1,offset+lcl_complex(dim_cmplx+1)%glb_indx(lcl_complex(dim_cmplx+1)%surf_indx(i,3))-1,&
    !   1,lcl_complex(dim_cmplx)%glb_indx(lcl_complex(dim_cmplx+1)%surf_indx(i,2))-1,&
    !   -real(lcl_complex(dim_cmplx+1)%surf_indx(i,1),iwp)*&
    !   lcl_complex(dim_cmplx)%inv_hdg_star(lcl_complex(dim_cmplx+1)%surf_indx(i,2)),&
    !   INSERT_VALUES,petsc_ier)
  END DO

  !-- assemble system matrix --
  CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- log time to build LHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/get_LHS_darcy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack
!* SYNOPSIS
  SUBROUTINE identify_crack(exit_cond)
!* PURPOSE
!*   identify faces to crack based on threshold value
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> identify faces to crack based on threshold value
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  LOGICAL, INTENT(OUT)  :: exit_cond          !>

  !-- local variables --
  INTEGER :: i                                !> loop and temporary indicies
  INTEGER :: type_bc                          !> type of boundary condition
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices
  LOGICAL  :: lcl_exit_cond = .FALSE.         !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  !-- set new cracks --
  num_flx = 0
  ALLOCATE(flx_indx(num_elm(dim_cmplx)))
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0 .AND. &
      ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/&
      lcl_complex(dim_cmplx)%prml_volume(i)) > crck_thrshld) THEN

      !-- impermeable boundary --
      CALL VecSetValues(sol,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      num_flx = num_flx + 1
      flx_indx(num_flx) = lcl_complex(dim_cmplx)%glb_indx(i) - 1
      lcl_complex(dim_cmplx)%bc_type(i) = -1
    END IF
  END DO
  !-- assemble solution vector  --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)

  CALL syncwrite_log('> identify_max_flux(): update LHS');
  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  CALL MatZeroRowsColumns(A,num_flx,(/ flx_indx(1:num_flx)/),1.d0,sol,b,petsc_ier)
  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  !-- clean up --
  DEALLOCATE(flx_indx)

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  IF (num_flx==0) lcl_exit_cond = .TRUE.
  CALL MPI_ALLREDUCE(lcl_exit_cond,exit_cond,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ier)

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack2
!* SYNOPSIS
  SUBROUTINE identify_crack2(exit_cond)
!* PURPOSE
!*   identify faces to crack based on threshold value
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> identify faces to crack based on threshold value
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  LOGICAL, INTENT(OUT)  :: exit_cond          !>

  !-- local variables --
  INTEGER :: i,j,indx                         !> loop and temporary indicies
  INTEGER :: type_bc                          !> type of boundary condition
  INTEGER :: num_flx,lcl_num_flx,glb_num_flx  !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices
  REAL(KIND=iwp), ALLOCATABLE :: flxs(:)      !> flx boundary condition indices
  REAL(KIND=iwp) :: flx, tmp_crck_thrshld     !> flx boundary condition indices

  LOGICAL  :: lcl_exit_cond = .FALSE.         !>
  CHARACTER(LEN=slen) :: msg                  !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  !-- set new cracks --
  lcl_num_flx = 0
  ALLOCATE(flx_indx(num_elm(dim_cmplx)),flxs(num_elm(dim_cmplx)))
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0) THEN
      flx = ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/lcl_complex(dim_cmplx)%prml_volume(i))
      IF (flx > crck_thrshld) THEN
        lcl_num_flx = lcl_num_flx + 1
        flx_indx(lcl_num_flx) = i
        flxs(lcl_num_flx) = flx
      END IF
    END IF
  END DO

  CALL MPI_ALLREDUCE(lcl_num_flx,glb_num_flx,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  write(msg,*) glb_num_flx,crck_thrshld
  CALL syncwrite_log(msg)

  tmp_crck_thrshld = 1.05*crck_thrshld
  DO WHILE (glb_num_flx > crcks_pstep)
    num_flx = lcl_num_flx
    lcl_num_flx = 0
    DO j = 1,num_flx
      i = flx_indx(j)
      IF (flxs(j) > tmp_crck_thrshld) THEN
        lcl_num_flx = lcl_num_flx + 1
        flx_indx(lcl_num_flx) = flx_indx(j)
        flxs(lcl_num_flx) = flxs(j)
      END IF
    END DO

    CALL MPI_ALLREDUCE(lcl_num_flx,glb_num_flx,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
    write(msg,*) glb_num_flx,tmp_crck_thrshld
    CALL syncwrite_log(msg)
    IF (glb_num_flx == 0) THEN
      tmp_crck_thrshld = 0.99*tmp_crck_thrshld
      glb_num_flx = crcks_pstep + 1
      lcl_num_flx = num_flx
    ELSE
      tmp_crck_thrshld = 1.05*tmp_crck_thrshld
    END IF
  END DO

  DO j = 1,lcl_num_flx
    !-- impermeable boundary --
    i = flx_indx(j)
    CALL VecSetValues(sol,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
      INSERT_VALUES,petsc_ier)
    flx_indx(j) = lcl_complex(dim_cmplx)%glb_indx(i) - 1
    lcl_complex(dim_cmplx)%bc_type(i) = -1
  END DO

  !-- assemble solution vector  --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)

  CALL syncwrite_log('> identify_max_flux(): update LHS');
  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  CALL MatZeroRowsColumns(A,lcl_num_flx,(/ flx_indx(1:lcl_num_flx)/),1.d0,sol,b,petsc_ier)
  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  !-- clean up --
  DEALLOCATE(flx_indx,flxs)

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  IF (lcl_num_flx==0) lcl_exit_cond = .TRUE.
  CALL MPI_ALLREDUCE(lcl_exit_cond,exit_cond,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ier)

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack3
!* SYNOPSIS
  SUBROUTINE identify_crack3(iter,exit_cond)
!* PURPOSE
!*   identify faces to crack based on max value
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> identify faces to crack based on max value
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  INTEGER, INTENT(IN)   :: iter              !>
  LOGICAL, INTENT(OUT)  :: exit_cond         !>

  !-- local variables --
  INTEGER :: i,junk,max_indx(2)              !> loop and temporary indicies
  REAL(KIND=iwp) :: glb_max_flx(3)           !>
  REAL(KIND=iwp) :: max_flx(3)               !>
  REAL(KIND=iwp) :: flx                      !>
  REAL(KIND=iwp) :: area                     !>
  CHARACTER(LEN=slen) :: msg                 !>
  INTEGER        :: status(MPI_STATUS_SIZE)  !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  !-- find new crack --
  max_flx = (/ 0.d0, dble(rank), 0.d0 /)
  max_indx = -1
  area = 0.d0
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0) THEN
      area = lcl_complex(dim_cmplx)%prml_volume(i)
      flx = ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/area)
      IF (flx > max_flx(1)) THEN
        max_indx = (/ i, lcl_complex(dim_cmplx)%glb_indx(i)-1 /)
        max_flx(1) = flx
        max_flx(2) = dble(max_indx(2))
        max_flx(3) = area
      END IF
    END IF
  END DO

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,max_flx(1:2),2,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,ier)

  !-- impermeable boundary --
  IF (max_flx(1)>crck_thrshld) THEN
    IF (max_indx(2)==int(max_flx(2))) THEN
      CALL VecSetValues(sol,1,max_indx(2),0.d0,INSERT_VALUES,petsc_ier)
      CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
      CALL MatZeroRowsColumns(A,1,(/ max_indx(2) /),1.d0,sol,b,petsc_ier)
      lcl_complex(dim_cmplx)%bc_type(max_indx(1)) = -1
    ELSE
      CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
      CALL MatZeroRowsColumns(A,0,(/0/),1.d0,sol,b,petsc_ier)
      max_flx=0.d0
    END IF
    CALL MPI_REDUCE(max_flx,glb_max_flx,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)
    IF (rank==root) THEN
      area = glb_max_flx(3)
      crck_area = crck_area + area
      WRITE(ulog_unit,*) iter, glb_max_flx(1), int(glb_max_flx(2)), area, area**(3.d0/2.d0), crck_area, crck_area**(3.d0/2.d0)
      CALL FLUSH(ulog_unit)
    END IF

    exit_cond = .FALSE.
  ELSE
    exit_cond = .TRUE.
  END IF

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/exchange_bndry_cond
!* SYNOPSIS
  SUBROUTINE exchange_bndry_cond(k)
!* PURPOSE
!*   exchange boundary conditions
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> exchange boundary conditions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER,INTENT(IN)     :: k                 !> simplicial order

    !-- local variables --
    INTEGER               :: kp                 !> simplicial order
    INTEGER               :: ptr,ptrp           !> pointer
    INTEGER               :: i,j                !> loop index

    INTEGER               :: junk               !> junk comm variable
    INTEGER               :: buffer_size        !> size of comm buffer
    INTEGER, ALLOCATABLE  :: sbuffer(:)         !> send buffer
    INTEGER, ALLOCATABLE  :: rbuffer(:)         !> receive buffer
    INTEGER, ALLOCATABLE  :: req(:)             !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE  :: status(:,:)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- pack send buffer --
    kp = k+1
    ALLOCATE(sbuffer(lcl_complex(kp)%num_send))
    DO j = 1,lcl_complex(kp)%num_send
      sbuffer(j) = lcl_complex(k)%bc_type(lcl_complex(kp)%send_indx(j,2))
    END DO

    !-- send data to adjacent processes --
    ptr = 1
    DO j = 1,num_adj_proc
      IF (num_send(k,j)==0) CYCLE
      buffer_size = num_send(k,j)
      CALL MPI_ISEND(sbuffer(ptr:ptr+buffer_size-1),buffer_size,&
        MPI_INTEGER,adj_proc(j),0,MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)
      ptr = ptr + buffer_size
    END DO

    !-- receive data from adjacent processes --
    ALLOCATE(&
      req(num_adj_proc),&
      status(MPI_STATUS_SIZE,num_adj_proc),&
      rbuffer(lcl_complex(kp)%num_recv))
    ptr = 1
    DO j = 1,num_adj_proc
      IF (num_recv(k,j)==0) THEN
        req(j) = MPI_REQUEST_NULL
        CYCLE
      END IF
      buffer_size = num_recv(k,j)
      CALL MPI_IRECV(rbuffer(ptr:ptr+buffer_size-1),buffer_size,&
        MPI_INTEGER,adj_proc(j),0,MPI_COMM_WORLD,req(j),ier)
      ptr = ptr + buffer_size
    END DO

    !-- unpack receive buffer --
    CALL MPI_WAITALL(num_adj_proc,req,status,ier)
    DO j = 1,lcl_complex(kp)%num_recv
      lcl_complex(k)%bc_type(lcl_complex(kp)%recv_indx(j,2)) = rbuffer(j)
    END DO

    !-- clean up --
    DEALLOCATE(req,status)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|exchange_bndry_cond
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/initialise_darcy2
!* SYNOPSIS
  SUBROUTINE initialise_darcy2()
!* PURPOSE
!*   initialise matrices and vectors for Darcy simulations
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> initialise matrices and vectors for Darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER         :: m        !> number of rows/columns in global solution variables
  INTEGER         :: offset   !> loop index

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> initialise_darcy()')
  curnt_time = MPI_WTIME()

  !-- allocate local solution variables --
  ALLOCATE(&
    lcl_complex(dim_cmplx)%prml_sol(num_elm(dim_cmplx),1),&
    lcl_complex(dim_cmplx+1)%whtny_sol(num_elm(dim_cmplx+1),dim_embbd),&
    lcl_complex(dim_cmplx+1)%dual_sol(num_elm(dim_cmplx+1),1))

  !-- log time to initialise the darcy flow problem --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/initialise_darcy2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/finalise_darcy2
!* SYNOPSIS
  SUBROUTINE finalise_darcy2()
!* PURPOSE
!*   clean up matrices and vectors for Darcy simulations
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
!> clean up matrices and vectors for Darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !--deallocate vectors, matrices, and solver --
  !-- vectors --
  CALL VecDestroy(b,petsc_ier)
  CALL VecDestroy(sol,petsc_ier)

  !-- matrices --
  CALL MatDestroy(Asub(1),petsc_ier)
  CALL MatDestroy(Asub(2),petsc_ier)
  CALL MatDestroy(Asub(3),petsc_ier)
  CALL MatDestroy(Asub(4),petsc_ier)
  CALL MatDestroy(A,petsc_ier)
  CALL MatDestroy(Sp,petsc_ier)

  !-- index set --
  CALL ISDestroy(isg(1),petsc_ier)
  CALL ISDestroy(isg(2),petsc_ier)

  !-- linear solver --
  CALL KSPDestroy(ksp_id,petsc_ier)

  RETURN

  END SUBROUTINE
! darcy_mod/finalise_darcy2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/get_RHS_darcy2
!* SYNOPSIS
  SUBROUTINE get_RHS_darcy2()
!* PURPOSE
!*   setup right hand side vector and boundary conditions for darcy simulations
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> setup right hand side vector and boundary conditions for darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER :: i,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=iwp)  :: flx                      !>
  Vec             :: b1,b2,sol1,sol2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> get_RHS_darcy()');  curnt_time = MPI_WTIME()

  CALL syncwrite_log('> get_RHS_darcy(): read_bndry_cond()');
  ! CALL read_bndry_cond()
  ! CALL exchange_bndry_cond(dim_cmplx)
  CALL read_bndry_cond2()
  CALL syncwrite_log_time();    curnt_time = MPI_WTIME()

  !-- vectors (RHS:b, solution:sol) --
  CALL MatGetOwnershipRange(A,m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(dim_cmplx+1) + glb_num_elm(dim_cmplx),sol,petsc_ier)
  CALL VecSetFromOptions(sol,petsc_ier);            CALL VecDuplicate(sol,b,petsc_ier)
  CALL VecGetSubVector(b,isg(1),b1,petsc_ier);      CALL VecGetSubVector(b,isg(2),b2,petsc_ier)
  CALL VecGetSubVector(sol,isg(1),sol1,petsc_ier);  CALL VecGetSubVector(sol,isg(2),sol2,petsc_ier)

  CALL MatGetOwnershipRange(Asub(1),m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(dim_cmplx),wrk,petsc_ier)
  CALL VecSetFromOptions(wrk,petsc_ier)

  !--  set known values for initial and boundary conditions (fluxes) --
  num_flx = 0
  ALLOCATE(flx_indx(num_elm(dim_cmplx)))
  DO i = 1,num_elm(dim_cmplx)
    type_bc = lcl_complex(dim_cmplx)%bc_type(i)
    IF (type_bc == -1) THEN
      !-- impermeable boundary --
      CALL VecSetValues(sol1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-1.d0,INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),&
        lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),&
        INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc == -3) THEN
      !-- impermeable boundary --
      CALL VecSetValues(sol1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-1.d0,INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),&
        lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),&
        INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc == -2) THEN
      !-- flux boundary --
      flx = -sum(vel*lcl_complex(dim_cmplx)%dual_dir(i,1:3))*lcl_complex(dim_cmplx)%prml_volume(i)
      CALL VecSetValues(sol1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(wrk,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(b1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc>=0 .AND. nz_init) THEN
      CALL VecSetValues(sol1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-sum(vel*&
        lcl_complex(dim_cmplx)%dual_dir(i,1:3))*lcl_complex(dim_cmplx)%prml_volume(i),&
        INSERT_VALUES,petsc_ier)
    ! ELSE
    !   WRITE(log_unit,'(A,A,A,I5)')'Error : undefined boundary type = ',type_bc
    !   CALL end_mpi()
    END IF
  END DO

  !--  set known values for initial conditions (pressures) --
  IF (nz_init)  THEN
    DO i = 1,num_elm(dim_cmplx+1)
      CALL VecSetValues(sol2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        -sum(vel*lcl_complex(dim_cmplx+1)%centers(i,1:3)) + p_ref,INSERT_VALUES,petsc_ier)
      END DO
  END IF

  !-- set known values for boundary conditions (pressures) --
  CALL syncwrite_log('> get_RHS_darcy(): assign bndry_cond()');
  DO i = 1,lcl_complex(dim_cmplx+1)%num_surf
    indx = lcl_complex(dim_cmplx+1)%surf_indx(i,2)
    type_bc = lcl_complex(dim_cmplx+1)%surf_indx(i,4)
    IF (type_bc>0) THEN
      !-- pressure boundary condition --
      CALL VecSetValues(b1,1,lcl_complex(dim_cmplx)%glb_indx(indx)-1,&
        -lcl_complex(dim_cmplx+1)%surf_indx(i,1)*bc_press(type_bc),&
        INSERT_VALUES,petsc_ier)
    END IF
  END DO

  !-- assemble sub-vectors --
  CALL VecAssemblyBegin(sol1,petsc_ier); CALL VecAssemblyEnd(sol1,petsc_ier)
  CALL VecAssemblyBegin(sol2,petsc_ier); CALL VecAssemblyEnd(sol2,petsc_ier)
  CALL VecAssemblyBegin(b1,petsc_ier); CALL VecAssemblyEnd(b1,petsc_ier)
  CALL VecAssemblyBegin(b2,petsc_ier); CALL VecAssemblyEnd(b2,petsc_ier)
  CALL VecAssemblyBegin(wrk,petsc_ier); CALL VecAssemblyEnd(wrk,petsc_ier)

  CALL syncwrite_log('> get_RHS_darcy(): update LHS');
  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  !CALL MatMultAdd(Asub(1),wrk,b1,b1,petsc_ier)
  CALL MatMultAdd(Asub(3),wrk,b2,b2,petsc_ier)

  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_INITIAL_MATRIX,Asub(3),petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- rebuild matrices and vectors --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  CALL VecRestoreSubVector(b,isg(1),b1,petsc_ier);      CALL VecRestoreSubVector(b,isg(2),b2,petsc_ier)
  CALL VecRestoreSubVector(sol,isg(1),sol1,petsc_ier);  CALL VecRestoreSubVector(sol,isg(2),sol2,petsc_ier)
  CALL VecDestroy(wrk,petsc_ier)

  !-- clean up --
  DEALLOCATE(flx_indx)

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/get_RHS_darcy2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/get_LHS_darcy2
!* SYNOPSIS
  SUBROUTINE get_LHS_darcy2()
!* PURPOSE
!*   setup left hand side matrix for darcy simulations
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> setup left hand side matrix for darcy simulations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- local variables --
  INTEGER         :: m          !> number of rows/columns in global solution variables
  INTEGER         :: i,j,offset !>
  REAL(KIND=iwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> get_LHS_darcy()')
  curnt_time = MPI_WTIME()

  !-- predefine some useful variables --
  offset = glb_num_elm(dim_cmplx)
  m = glb_num_elm(dim_cmplx+1) + glb_num_elm(dim_cmplx)

  !-- system matrix - block 00 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(1),petsc_ier)
  CALL MatSetSizes(Asub(1),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(dim_cmplx),glb_num_elm(dim_cmplx),petsc_ier)
  CALL MatSetType(Asub(1), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(1),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetFromOptions(Asub(1),petsc_ier)


  ! DO i = 1,num_elm(dim_cmplx+1)
  !   lcl_complex(dim_cmplx+1)%prml_sol(i,1) = -mu/k*(rand()-0.5d0)
  ! END DO
  ! CALL exchange_prml_sol(dim_cmplx+1)
  ! DO i = 1,num_elm(dim_cmplx)
  !   IF (lcl_complex(dim_cmplx)%num_cobndry(i) == 1) THEN
  !     lcl_complex(dim_cmplx)%hdg_star(i) = lcl_complex(dim_cmplx)%hdg_star(i) * &
  !       lcl_complex(dim_cmplx+1)%prml_sol(lcl_complex(dim_cmplx)%cobndry(i)%indx(1),1)
  !   ELSE
  !     indx1 = lcl_complex(dim_cmplx)%cobndry(i)%indx(1)
  !     pt1 = lcl_complex(dim_cmplx+1)%centers(indx1,:)
  !     s1 = lcl_complex(dim_cmplx+1)%prml_sol(indx1,1)
  !
  !     pt2 = lcl_complex(dim_cmplx)%centers(i,:)
  !
  !     indx2 = lcl_complex(dim_cmplx)%cobndry(i)%indx(2)
  !     pt3 = lcl_complex(dim_cmplx+1)%centers(indx2,:)
  !     s2 = lcl_complex(dim_cmplx+1)%prml_sol(,1)
  !
  !     dx1 = pt2 - pt1
  !     dx2 = pt3 - pt2
  !     dx1 = sqrt(dot_product(dx1,dx1))
  !     dx1 = sqrt(dot_product(dx2,dx2))
  !
  !     lcl_complex(dim_cmplx)%hdg_star(i) = (dx1*s1 + dx2*s2) / (dx1+dx2)
  !
  !   END IF
  !
  ! END DO

  !-- equations for Darcy's law --
  CALL RANDOM_SEED()
  DO i = 1,num_elm(dim_cmplx)
    CALL RANDOM_NUMBER(rnd)
    IF (lcl_complex(dim_cmplx)%num_cobndry(i)>0) &
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-(1.d0+rnd_max*2.0d0*(rnd-0.5d0))*&
        mu/k*sign(1.d0,lcl_complex(dim_cmplx)%hdg_star(i))*&
        max(abs(lcl_complex(dim_cmplx)%hdg_star(i)),smallh),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(Asub(1),MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- system matrix - block 10 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(3),petsc_ier)
  CALL MatSetSizes(Asub(3),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(dim_cmplx+1),glb_num_elm(dim_cmplx),petsc_ier)
  CALL MatSetType(Asub(3), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(3),4,PETSC_NULL_INTEGER,4,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetOption(Asub(3),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(3),petsc_ier)
  !-- continuity equations: dual edges corresponding to internal primal faces --
  DO i = 1,num_elm(dim_cmplx+1)
    CALL MatSetValues(Asub(3),1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
      lcl_complex(dim_cmplx+1)%num_bndry(i),&
      lcl_complex(dim_cmplx)%glb_indx(lcl_complex(dim_cmplx+1)%bndry(i)%indx)-1,&
      real(lcl_complex(dim_cmplx+1)%bndry(i)%sgn,iwp),INSERT_VALUES,petsc_ier)
  END DO

  !-- continuity equations: dual edges corresponding to surface primal faces --
  DO i = 1,lcl_complex(dim_cmplx+1)%num_surf
    CALL MatSetValues(Asub(3),&
      1,lcl_complex(dim_cmplx+1)%glb_indx(lcl_complex(dim_cmplx+1)%surf_indx(i,3))-1,&
      1,lcl_complex(dim_cmplx)%glb_indx(lcl_complex(dim_cmplx+1)%surf_indx(i,2))-1,&
      -real(lcl_complex(dim_cmplx+1)%surf_indx(i,1),iwp),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- system matrix - block 01 --
  CALL MatTranspose(Asub(3),MAT_INITIAL_MATRIX,Asub(2),petsc_ier)
  CALL MatSetOption(Asub(2),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(2),petsc_ier)

  !-- system matrix - block 11 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(4),petsc_ier)
  CALL MatSetSizes(Asub(4),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(dim_cmplx+1),glb_num_elm(dim_cmplx+1),petsc_ier)
  CALL MatSetType(Asub(4), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(4),0,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetFromOptions(Asub(4),petsc_ier)
  CALL MatAssemblyBegin(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(Asub(4),MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- system matrix (LHS:A) --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  CALL MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)
  CALL MatNestGetISs(A,isg,PETSC_NULL_IS,petsc_ier)

  !-- system Schur Complement --
  ! CALL MatCreate(MPI_COMM_WORLD,Sp,petsc_ier)
  ! CALL MatSetSizes(Sp,PETSC_DECIDE,PETSC_DECIDE,&
  !   glb_num_elm(dim_cmplx)+glb_num_elm(dim_cmplx+1),&
  !   glb_num_elm(dim_cmplx)+glb_num_elm(dim_cmplx+1),petsc_ier)
  ! CALL MatSetType(Sp, MATMPIAIJ,petsc_ier)
  ! CALL MatMPIAIJSetPreallocation(Sp,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,petsc_ier)
  ! CALL MatSetFromOptions(Sp,petsc_ier)

  !-- log time to build LHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/get_LHS_darcy2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack4
!* SYNOPSIS
  SUBROUTINE identify_crack4(iter,exit_cond)
!* PURPOSE
!*   identify faces to crack based on threshold value
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> identify faces to crack based on threshold value
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  INTEGER, INTENT(IN)   :: iter               !>
  LOGICAL, INTENT(OUT)  :: exit_cond          !>

  !-- local variables --
  INTEGER :: i,j                              !> loop and temporary indicies
  INTEGER :: type_bc                          !> type of boundary condition
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices
  LOGICAL  :: lcl_exit_cond = .FALSE.         !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  !-- set new cracks --
  num_flx = 0
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0 .AND. &
      ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/&
      lcl_complex(dim_cmplx)%prml_volume(i)) > crck_thrshld) THEN

      !-- impermeable boundary --
      CALL VecSetValues(sol,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (small, j=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    END IF
  END DO
  !-- assemble solution vector  --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)

  !-- assemble solution vector  --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  !-- clean up --
  DEALLOCATE(flx_indx)

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  IF (num_flx==0) lcl_exit_cond = .TRUE.
  CALL MPI_ALLREDUCE(lcl_exit_cond,exit_cond,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ier)

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack5
!* SYNOPSIS
  SUBROUTINE identify_crack5(iter,exit_cond)
!* PURPOSE
!*    identify faces to crack based on max value
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!> identify faces to crack based on max value
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  INTEGER, INTENT(IN)   :: iter      !>
  LOGICAL, INTENT(OUT)  :: exit_cond      !>

  !-- local variables --
  INTEGER :: i,junk,max_indx(2)              !> loop and temporary indicies
  REAL(KIND=iwp) :: glb_max_flx(3)      !>
  REAL(KIND=iwp) :: max_flx(3)      !>
  REAL(KIND=iwp) :: flx      !>
  REAL(KIND=iwp) ::  area      !>
  CHARACTER(LEN=slen) :: msg      !>
  INTEGER        :: status(MPI_STATUS_SIZE)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  !-- find new crack --
  max_flx = (/ 0.d0, dble(rank), 0.d0 /)
  max_indx = -1
  area = 0.d0
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0) THEN
      area = lcl_complex(dim_cmplx)%prml_volume(i)
      flx = ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/area)
      IF (flx > max_flx(1)) THEN
        max_indx = (/ i, lcl_complex(dim_cmplx)%glb_indx(i)-1 /)
        max_flx(1) = flx
        max_flx(2) = dble(max_indx(2))
        max_flx(3) = area
      END IF
    END IF
  END DO

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,max_flx(1:2),2,MPI_2DOUBLE_PRECISION,&
    MPI_MAXLOC,MPI_COMM_WORLD,ier)

  !-- impermeable boundary --
  IF (max_flx(1)>crck_thrshld) THEN
    IF (max_indx(2)==int(max_flx(2))) THEN
      CALL VecSetValues(sol,1,max_indx(2),0.d0,INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,max_indx(2),1,max_indx(2),1.d0,INSERT_VALUES,&
        petsc_ier)
      CALL MatSetValues(Asub(2),1,max_indx(2),&
        lcl_complex(dim_cmplx)%num_cobndry(max_indx(1)),&
        lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(max_indx(1))%indx)-1,&
        (/ (small, i=1,lcl_complex(dim_cmplx)%num_cobndry(max_indx(1))) /),&
        INSERT_VALUES,petsc_ier)
      lcl_complex(dim_cmplx)%bc_type(max_indx(1)) = -1
    ELSE
      max_flx=0.d0
    END IF
    CALL MPI_REDUCE(max_flx,glb_max_flx,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
      MPI_COMM_WORLD,ier)
    IF (rank==root) THEN
      area = glb_max_flx(3)
      crck_area = crck_area + area
      WRITE(ulog_unit,*) iter, glb_max_flx(1), int(glb_max_flx(2)), area, &
        area**(3.d0/2.d0), crck_area, crck_area**(3.d0/2.d0)
      CALL FLUSH(ulog_unit)
    END IF

    !-- assemble solution vector  --
    CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
    CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
    CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
    ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

    exit_cond = .FALSE.
  ELSE
    exit_cond = .TRUE.
  END IF

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack6
!* SYNOPSIS
  SUBROUTINE identify_crack6(iter,exit_cond)
!* PURPOSE
!*   Introduce random crack to face not on the boundary or already cracked
!* INPUTS
!*   Name                    Description
!*   iter
!*   exit_cond
!* OUTPUTS
!*   Name                    Description
!*
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/12/09: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/12/09
!>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  INTEGER, INTENT(IN)   :: iter      !>
  LOGICAL, INTENT(OUT)  :: exit_cond      !>

  !-- local variables --
  INTEGER :: i,j,junk,indx             !> loop and temporary indicies
  REAL(KIND=iwp) ::  area      !>
  CHARACTER(LEN=slen) :: msg      !>
  INTEGER        :: status(MPI_STATUS_SIZE)        !> size of MPI comm buffer
  REAL(KIND=iwp)  :: rnd
  LOGICAL :: search

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  search = .TRUE.
  area = 0.d0
  DO WHILE(search)
    !-- find random face --
    IF (rank==root) THEN
      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(rnd)
      indx = nint(rnd*(glb_num_elm(dim_cmplx)-1))+1
    END IF

    !-- check if it is an interface --
    CALL MPI_BCAST(indx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    !-- if local and interface, then end search --
    DO i = 1,num_elm(dim_cmplx)
      IF (lcl_complex(dim_cmplx)%glb_indx(i) == indx .AND. &
        lcl_complex(dim_cmplx)%bc_type(i) == 0) THEN
        area = lcl_complex(dim_cmplx)%prml_volume(i)

        CALL VecSetValues(sol,1,indx-1,0.d0,INSERT_VALUES,petsc_ier)
        CALL MatSetValues(Asub(1),1,indx-1,1,indx-1,1.d0,INSERT_VALUES,&
          petsc_ier)
        CALL MatSetValues(Asub(2),1,indx-1,&
          lcl_complex(dim_cmplx)%num_cobndry(i),&
          lcl_complex(dim_cmplx+1)%glb_indx(&
          lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
          (/ (0.d0, j=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),&
          INSERT_VALUES,petsc_ier)

        lcl_complex(dim_cmplx)%bc_type(i) = -1
        search = .FALSE.
      END IF
    END DO

    !-- end search? --
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,search,1,MPI_LOGICAL,MPI_LAND,&
      MPI_COMM_WORLD,ier)
  END DO

  IF (rank==root) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,area,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
      MPI_COMM_WORLD,ier)
    crck_area = crck_area + area
    WRITE(ulog_unit,*) iter, 0.d0, indx, area, area**(3.d0/2.d0), &
      crck_area, crck_area**(3.d0/2.d0)
    CALL FLUSH(ulog_unit)
  ELSE
    CALL MPI_REDUCE(area,area,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
      MPI_COMM_WORLD,ier)
  END IF

  !-- assemble solution vector  --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  exit_cond = .FALSE.

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack6
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* darcy_mod/identify_crack7
!* SYNOPSIS
  SUBROUTINE identify_crack7(iter,exit_cond)
!* PURPOSE
!*
!* INPUTS
!*   Name                    Description
!*
!* OUTPUTS
!*   Name                    Description
!*
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
!>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments --
  INTEGER, INTENT(IN)   :: iter      !>
  LOGICAL, INTENT(OUT)  :: exit_cond      !>

  !-- local variables --
  INTEGER :: i,junk,max_indx(2)              !> loop and temporary indicies
  REAL(KIND=iwp) :: glb_max_flx(3)      !>
  REAL(KIND=iwp) :: max_flx(3)      !>
  REAL(KIND=iwp) :: flx      !>
  REAL(KIND=iwp) ::  area      !>
  CHARACTER(LEN=slen) :: msg      !>
  INTEGER        :: status(MPI_STATUS_SIZE)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL syncwrite_log('> update_LRHS_darcy_crack()');  curnt_time = MPI_WTIME()

  !-- find new crack --
  max_flx = (/ large, dble(rank), 0.d0 /)
  max_indx = -1
  area = 0.d0
  write(*,*) "----------------------------"
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0) THEN
      area = lcl_complex(dim_cmplx)%prml_volume(i)
      flx = ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/area)
      IF (flx < max_flx(1)) THEN
        max_indx = (/ i, lcl_complex(dim_cmplx)%glb_indx(i)-1 /)
        max_flx(1) = flx
        max_flx(2) = dble(max_indx(2))
        max_flx(3) = area
      END IF
    END IF
  END DO

  write(*,*) max_flx, max_indx

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,max_flx(1:2),2,MPI_2DOUBLE_PRECISION,&
    MPI_MINLOC,MPI_COMM_WORLD,ier)
  write(*,*) max_flx, max_indx

  !-- impermeable boundary --
  IF (max_indx(2)==int(max_flx(2))) THEN
    CALL VecSetValues(sol,1,max_indx(2),0.d0,INSERT_VALUES,petsc_ier)
    CALL MatSetValues(Asub(1),1,max_indx(2),1,max_indx(2),1.d0,INSERT_VALUES,&
      petsc_ier)
    CALL MatSetValues(Asub(2),1,max_indx(2),&
      lcl_complex(dim_cmplx)%num_cobndry(max_indx(1)),&
      lcl_complex(dim_cmplx+1)%glb_indx(&
      lcl_complex(dim_cmplx)%cobndry(max_indx(1))%indx)-1,&
      (/ (small, i=1,lcl_complex(dim_cmplx)%num_cobndry(max_indx(1))) /),&
      INSERT_VALUES,petsc_ier)
    lcl_complex(dim_cmplx)%bc_type(max_indx(1)) = -1
  ELSE
    max_flx=large
  END IF

  CALL MPI_REDUCE(max_flx,glb_max_flx,3,MPI_DOUBLE_PRECISION,MPI_MIN,0,&
    MPI_COMM_WORLD,ier)
  write(*,*) max_flx, max_indx

  IF (rank==root) THEN
    area = glb_max_flx(3)
    crck_area = crck_area + area
    WRITE(ulog_unit,*) iter, glb_max_flx(1), int(glb_max_flx(2)), area, &
    area**(3.d0/2.d0), crck_area, crck_area**(3.d0/2.d0)
    CALL FLUSH(ulog_unit)
  END IF

  !-- assemble solution vector  --
  CALL VecAssemblyBegin(sol,petsc_ier); CALL VecAssemblyEnd(sol,petsc_ier)
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  exit_cond = .FALSE.

  IF (glb_max_flx(1)==large) exit_cond = .TRUE.

  !-- log time to construct RHS --
  CALL syncwrite_log_time()

  RETURN

  END SUBROUTINE
! darcy_mod/identify_crack7
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! darcy_mod
!===============================================================================
!
