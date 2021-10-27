!
!===============================================================================
!-- Diffusion Module
!> Module contains routines specifically related to diffusion equations
!===============================================================================
!/****/h* modules|phy_diffusion_flow/diffusion_mod
MODULE diffusion_mod
!* PURPOSE
!*   Module contains routines specifically related to diffusion equations
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 general mpi routines (start, end, syncwrite,etc)
!*   io_mod                  IO functions and routines
!*   solver_mod              Solver routines (PETSc)
!*   petsc.h                 PETSc variables and routines
!*   diffusion_vars.inc      diffusion specific variables
!* CONTAINS
!*   Subroutine              Purpose
!*   read_input_diffusion    read inputs, and check parameters (phys/io/solver)
!*   check_param_diffusion   check parameters for diffusion simulations
!*   initialise_diffusion
!*   finalise_diffusion
!*   get_RHS_diffusion
!*   get_LHS_diffusion
!*   exchange_bndry_cond
!*   identify_threshold_crack
!*   identify_random_crack
!*   identify_minFD_crack
!*   identify_maxFD_crack
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Module contains routines specifically related to diffusion equations
!===============================================================================
! TO DO:
! - set SOLVER defaults
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod
  USE io_mod
  USE solver_mod
  USE time_marching_mod

  IMPLICIT NONE

  !-- include diffusion specific variables --
#include <petsc/finclude/petsc.h>
#include "diffusion_vars.inc"

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/read_input_diffusion
!* SYNOPSIS
  SUBROUTINE read_input_diffusion()
!* PURPOSE
!*   Read input file for user defined parameters, and perform checks to ensure
!*   reasonable values are set
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

  IMPLICIT NONE

  !-- local variables --
  INTEGER               :: i,k            !> counters
  CHARACTER(LEN=SLEN)   :: fname          !> file name
  CHARACTER(LEN=SLEN)   :: buffer         !> buffer for command line arguments
  LOGICAL               :: err = .FALSE.  !> error logical
  INTEGER               :: mystat(5)      !> namelist io stat
  REAL(KIND=PGMSiwp)        :: realbuf(42)    !> real MPI buffer
  INTEGER               :: intbuf(8)      !> integer MPI buffer
  LOGICAL               :: logbuf(3)      !> logical MPI buffer
  CHARACTER(LEN=4*SLEN) :: charbuf        !> character MPI buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- check for user supplied input file name, otherwise use default --
  IF (rank == root) THEN  !-- get input filename --
    CALL GETARG(1,buffer); IF (LEN_TRIM(buffer) == 0) THEN; fname = input_file
    ELSE; READ(buffer,*) fname; END IF
    CALL rootwrite('   - input file name: '//fname,log_unit)

    !-- open input file --
    OPEN(input_unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
    err = chkerr_log(ier,'OPEN - error opening '//TRIM(fname),log_unit)

    IF (.NOT. err) THEN; !-- read input file --
      REWIND(input_unit); READ(input_unit,nml=diff_param,iostat=mystat(1))
      REWIND(input_unit); READ(input_unit,nml=crack_param,iostat=mystat(2))
      REWIND(input_unit); READ(input_unit,nml=io_param,iostat=mystat(3))
      REWIND(input_unit); READ(input_unit,nml=time_param,iostat=mystat(4))
      REWIND(input_unit); READ(input_unit,nml=solver_param,iostat=mystat(5))

      !-- check namelists in input file are correct --
      IF (ANY(mystat>0)) THEN; err = chkerr_log(MAXVAL(mystat),&
        'READ - namelist error in '//TRIM(fname),log_unit)
      ELSEIF (ANY(mystat<-1)) THEN; err = chkerr_log(MINVAL(mystat),&
        'READ - namelist error in '//TRIM(fname),log_unit); END IF

      !-- close input file --
      CLOSE(input_unit)
    END IF; CALL chkerrMPI_BCAST(err)  !-- complete check for errors --

    !-- pack and broadcast input data --
    charbuf = trim(mesh_prefix)//','//trim(sol_prefix)//','//trim( &
      unstdy_prefix)//','//trim(log_prefix)//','//trim(time_marching_method)//&
      ','//trim(ksp_type)//','//trim(pc_type)//','//trim(pc_reorder)
    realbuf = (/ Dk, ref_val, blk_dir, bc_vals, crck_thrshld, rnd_max, dt, &
      time, rel_tol, abs_tol, div_tol /)
    intbuf = (/ crck_type, max_crcks, crcks_pstep, initial_step, final_step, &
      max_iter, output_frqcy, LEN_TRIM(charbuf) /)
    logbuf = (/ nz_init, solver_monitor, sol_output /)
    CALL MPI_BCAST(intbuf,8,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(realbuf,42,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(logbuf,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(charbuf,intbuf(8),MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  ELSE !-- non-root processes
    CALL chkerrMPI_BCAST(err)  !-- complete check for errors --

    !-- recieve broadcast input data and unpack --
    CALL MPI_BCAST(intbuf,8,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(realbuf,42,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(logbuf,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(charbuf,intbuf(8),MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
    Dk              = realbuf(1);   ref_val         = realbuf(2)
    blk_dir         = realbuf(3:5); bc_vals         = realbuf(6:35)
    crck_thrshld    = realbuf(36);  rnd_max         = realbuf(37)
    dt              = realbuf(38);  time            = realbuf(39)
    rel_tol         = realbuf(40);  abs_tol         = realbuf(41)
    div_tol         = realbuf(42)
    crck_type       = intbuf(1);    max_crcks       = intbuf(2)
    crcks_pstep     = intbuf(3);    initial_step    = intbuf(4)
    final_step      = intbuf(5);    max_iter        = intbuf(6)
    output_frqcy    = intbuf(7)
    nz_init         = logbuf(1);    solver_monitor  = logbuf(2);
    sol_output      = logbuf(3);
    i=1;   k = INDEX(charbuf(i:),","); mesh_prefix          = charbuf(i:i+k-2)
    i=i+k; k = INDEX(charbuf(i:),","); sol_prefix           = charbuf(i:i+k-2)
    i=i+k; k = INDEX(charbuf(i:),","); unstdy_prefix        = charbuf(i:i+k-2)
    i=i+k; k = INDEX(charbuf(i:),","); log_prefix           = charbuf(i:i+k-2)
    i=i+k; k = INDEX(charbuf(i:),","); time_marching_method = charbuf(i:i+k-2)
    i=i+k; k = INDEX(charbuf(i:),","); ksp_type             = charbuf(i:i+k-2)
    i=i+k; k = INDEX(charbuf(i:),","); pc_type              = charbuf(i:i+k-2)
    i=i+k; k = LEN_TRIM(charbuf);      pc_reorder           = charbuf(i:k)
  END IF; RETURN

  END SUBROUTINE
! diffusion_mod/read_input_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/initialise_diffusion
!* SYNOPSIS
  SUBROUTINE initialise_diffusion()
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

  !-- local variables --
  INTEGER         :: m        !> number of rows/columns in global solution variables
  INTEGER         :: offset   !> loop index

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- allocate local solution variables --
  ALLOCATE(&
    lcl_complex(dim_cmplx)%prml_sol(num_elm(dim_cmplx),1),&
    lcl_complex(dim_cmplx+1)%dual_sol(num_elm(dim_cmplx+1),1),&
    lcl_complex(dim_cmplx+1)%whtny_sol(num_elm(dim_cmplx+1),dim_embbd)); RETURN

  END SUBROUTINE
! diffusion_mod/initialise_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/initialise_diffusion_d
!* SYNOPSIS
  SUBROUTINE initialise_diffusion_d()
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

  !-- local variables --
  INTEGER         :: m        !> number of rows/columns in global solution variables
  INTEGER         :: offset   !> loop index

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- allocate local solution variables --
  ALLOCATE(&
    lcl_complex(2)%dual_sol(num_elm(2),1),&
    lcl_complex(1)%prml_sol(num_elm(1),1)); RETURN

  END SUBROUTINE
! diffusion_mod/initialise_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/finalise_diffusion
!* SYNOPSIS
  SUBROUTINE finalise_diffusion()
!* PURPOSE
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

  INTEGER :: i

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !--deallocate vectors, matrices, and solver --
  !-- vectors --
  CALL VecDestroy(b,petsc_ier)
  CALL VecDestroy(q,petsc_ier); CALL VecDestroy(dq,petsc_ier)

  IF (time_dependent) THEN
    DO i = 1,time_marching%r; CALL VecDestroy(q_step(i),petsc_ier); END DO
    DO i = 1,time_marching%s; CALL VecDestroy(r_stage(i),petsc_ier); END DO
    CALL VecDestroy(q_scl,petsc_ier); CALL VecDestroy(scl,petsc_ier)
  END IF

  !-- matrices --
  CALL MatDestroy(Asub(1),petsc_ier); CALL MatDestroy(Asub(2),petsc_ier)
  CALL MatDestroy(Asub(3),petsc_ier); CALL MatDestroy(Asub(4),petsc_ier)
  CALL MatDestroy(A,petsc_ier);       !CALL MatDestroy(Sp,petsc_ier)

  !-- index set --
  CALL ISDestroy(isg(1),petsc_ier); CALL ISDestroy(isg(2),petsc_ier)

  !-- linear solver --
  CALL KSPDestroy(ksp_id,petsc_ier); RETURN

  END SUBROUTINE
! diffusion_mod/finalise_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_RHS_diffusion
!* SYNOPSIS
  SUBROUTINE setup_RHS_diffusion()
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

  !-- local variables --
  INTEGER :: i,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- vectors --
  CALL MatGetOwnershipRange(A,m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(dim_cmplx+1) + glb_num_elm(dim_cmplx),q,petsc_ier)
  CALL VecSetFromOptions(q,petsc_ier); CALL VecDuplicate(q,dq,petsc_ier)
  CALL VecDuplicate(q,r,petsc_ier);    CALL VecDuplicate(q,b,petsc_ier)
  IF (time_dependent) THEN
    IF (ALLOCATED(q_step)) DEALLOCATE(q_step)
    IF (ALLOCATED(r_stage)) DEALLOCATE(r_stage)
    ALLOCATE(q_step(time_marching%r),r_stage(time_marching%s))
    DO i=1,time_marching%r; CALL VecDuplicate(q,q_step(i),petsc_ier); END DO
    DO i=1,time_marching%s; CALL VecDuplicate(q,r_stage(i),petsc_ier); END DO
  END IF

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_RHS_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_RHS_diffusion
!* SYNOPSIS
  SUBROUTINE setup_IC_BC_diffusion()
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

  !-- local variables --
  INTEGER :: i,j,k,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL read_bndry_cond2()

  !-- vectors --
  CALL VecGetSubVector(b,isg(1),b1,petsc_ier)
  CALL VecGetSubVector(b,isg(2),b2,petsc_ier)
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier)
  CALL VecGetSubVector(q,isg(2),q2,petsc_ier)

  CALL MatGetOwnershipRange(Asub(1),m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(dim_cmplx),wrk,petsc_ier)
  CALL VecSetFromOptions(wrk,petsc_ier)

  !--  set known values for initial and boundary conditions (fluxes) --
  num_flx = 0
  ALLOCATE(flx_indx(num_elm(dim_cmplx)))
  DO i = 1,num_elm(dim_cmplx)
    type_bc = lcl_complex(dim_cmplx)%bc_type(i)
    IF (type_bc == -2) THEN
      !-- flux boundary --
      flx = -sum(blk_dir*lcl_complex(dim_cmplx)%dual_dir(i,1:3))*lcl_complex(dim_cmplx)%prml_volume(i)
      lcl_complex(dim_cmplx)%prml_sol(i,1) = flx
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(wrk,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(b1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      lcl_complex(dim_cmplx)%hdg_star(i) = 1.d0
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc == -1) THEN
      !-- impermeable boundary --
      lcl_complex(dim_cmplx)%prml_sol(i,1) = 0.d0
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      lcl_complex(dim_cmplx)%hdg_star(i) = 1.d0
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc>0) THEN
      !-- pressure boundary condition --
      DO j=1,lcl_complex(dim_cmplx)%num_cobndry(i)
        indx = lcl_complex(dim_cmplx)%cobndry(i)%indx(j)
        DO k = 1,lcl_complex(dim_cmplx+1)%num_bndry(indx)
          IF (lcl_complex(dim_cmplx+1)%bndry(indx)%indx(k)==i) THEN

            CALL VecSetValues(b1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
              lcl_complex(dim_cmplx+1)%bndry(indx)%sgn(k)*bc_vals(type_bc),&
              INSERT_VALUES,petsc_ier)

            CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
              lcl_complex(dim_cmplx)%prml_sol(i,1),INSERT_VALUES,petsc_ier)
          END IF
        END DO
      END DO
    ELSE
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%prml_sol(i,1),INSERT_VALUES,petsc_ier)
    END IF
  END DO

  !--  set known values for initial conditions (pressures) --
  DO i = 1,num_elm(dim_cmplx+1)
    CALL VecSetValues(q2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
      lcl_complex(dim_cmplx+1)%dual_sol(i,1),INSERT_VALUES,petsc_ier)
  END DO

  !-- assemble sub-vectors --
  CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
  CALL VecAssemblyBegin(q2,petsc_ier); CALL VecAssemblyEnd(q2,petsc_ier)
  CALL VecAssemblyBegin(b1,petsc_ier); CALL VecAssemblyEnd(b1,petsc_ier)
  CALL VecAssemblyBegin(b2,petsc_ier); CALL VecAssemblyEnd(b2,petsc_ier)
  CALL VecAssemblyBegin(wrk,petsc_ier); CALL VecAssemblyEnd(wrk,petsc_ier)

  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatMultAdd(Asub(3),wrk,b2,b2,petsc_ier)

  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_INITIAL_MATRIX,Asub(3),petsc_ier)
  CALL MatScale(Asub(3),-1.d0,petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- rebuild matrices and vectors --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  CALL VecRestoreSubVector(b,isg(1),b1,petsc_ier);   CALL VecRestoreSubVector(b,isg(2),b2,petsc_ier)
  CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier);   CALL VecRestoreSubVector(q,isg(2),q2,petsc_ier)
  CALL VecDestroy(wrk,petsc_ier)

  !-- set previous solution vectors (TODO: proper multistep startup) --
  IF (time_dependent) THEN
    CALL VecDuplicate(q,scl,petsc_ier);  CALL VecDuplicate(q,q_scl,petsc_ier)
    CALL VecGetSubVector(scl,isg(1),q1,petsc_ier)
    CALL VecGetSubVector(scl,isg(2),q2,petsc_ier)
    CALL VecSet(q1,0.d0,petsc_ier)
    DO i = 1,num_elm(dim_cmplx+1)
      CALL VecSetValues(q2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx+1)%inv_hdg_star(i),INSERT_VALUES,petsc_ier)
    END DO
    CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
    CALL VecAssemblyBegin(q2,petsc_ier); CALL VecAssemblyEnd(q2,petsc_ier)
    CALL VecRestoreSubVector(scl,isg(1),q1,petsc_ier)
    CALL VecRestoreSubVector(scl,isg(2),q2,petsc_ier)

    CALL VecPointwiseMult(q_scl,q,scl,petsc_ier)

    DO i=1,time_marching%r
      CALL VecCopy(q_scl,q_step(i),petsc_ier)
    END DO
  END IF

  !-- clean up --
  DEALLOCATE(flx_indx)

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_RHS_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_IC_BC_diffusion2
!* SYNOPSIS
  SUBROUTINE setup_IC_BC_diffusion2()
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

  !-- local variables --
  INTEGER :: i,j,k,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL read_bndry_cond3()

  !-- vectors --
  CALL VecGetSubVector(b,isg(1),b1,petsc_ier);   CALL VecGetSubVector(b,isg(2),b2,petsc_ier)
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier);   CALL VecGetSubVector(q,isg(2),q2,petsc_ier)

  CALL MatGetOwnershipRange(Asub(1),m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(dim_cmplx),wrk,petsc_ier)
  CALL VecSetFromOptions(wrk,petsc_ier)

  !--  set known values for initial and boundary conditions (fluxes) --
  num_flx = 0
  ALLOCATE(flx_indx(num_elm(dim_cmplx)))
  DO i = 1,num_elm(dim_cmplx)
    type_bc = lcl_complex(dim_cmplx)%bc_type(i)
    IF (type_bc == -10) THEN
      !--
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%prml_sol(i,1),INSERT_VALUES,petsc_ier)
      !--
      lcl_complex(dim_cmplx)%hdg_star(i) = lcl_complex(dim_cmplx)%hdg_star(i)/&
        lcl_complex(dim_cmplx)%bc_val(i,1)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,lcl_complex(dim_cmplx)%hdg_star(i),&
        INSERT_VALUES,petsc_ier)
      !--
      lcl_complex(dim_cmplx)%bc_type(i) = 0
    ELSEIF (type_bc == -2) THEN
      !-- flux boundary --
      flx = -sum(blk_dir*lcl_complex(dim_cmplx)%dual_dir(i,1:3))*lcl_complex(dim_cmplx)%prml_volume(i)
      lcl_complex(dim_cmplx)%prml_sol(i,1) = flx
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(wrk,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(b1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      lcl_complex(dim_cmplx)%hdg_star(i) = 1.d0
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc == -1) THEN
      !-- impermeable boundary --
      lcl_complex(dim_cmplx)%prml_sol(i,1) = 0.d0
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      lcl_complex(dim_cmplx)%hdg_star(i) = 1.d0
      CALL MatSetValues(Asub(2),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%num_cobndry(i),lcl_complex(dim_cmplx+1)%glb_indx(&
        lcl_complex(dim_cmplx)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(dim_cmplx)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc>0) THEN
      !-- pressure boundary condition --
      DO j=1,lcl_complex(dim_cmplx)%num_cobndry(i)
        indx = lcl_complex(dim_cmplx)%cobndry(i)%indx(j)
        DO k = 1,lcl_complex(dim_cmplx+1)%num_bndry(indx)
          IF (lcl_complex(dim_cmplx+1)%bndry(indx)%indx(k)==i) THEN

            CALL VecSetValues(b1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
              lcl_complex(dim_cmplx+1)%bndry(indx)%sgn(k)*bc_vals(type_bc),&
              INSERT_VALUES,petsc_ier)

            CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
              lcl_complex(dim_cmplx)%prml_sol(i,1),INSERT_VALUES,petsc_ier)
          END IF
        END DO
      END DO
    ELSE
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx)%prml_sol(i,1),INSERT_VALUES,petsc_ier)
    END IF
  END DO

  !--  set known values for initial conditions (pressures) --
  DO i = 1,num_elm(dim_cmplx+1)
    CALL VecSetValues(q2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
      lcl_complex(dim_cmplx+1)%dual_sol(i,1),INSERT_VALUES,petsc_ier)
  END DO

  !-- assemble sub-vectors --
  CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
  CALL VecAssemblyBegin(q2,petsc_ier); CALL VecAssemblyEnd(q2,petsc_ier)
  CALL VecAssemblyBegin(b1,petsc_ier); CALL VecAssemblyEnd(b1,petsc_ier)
  CALL VecAssemblyBegin(b2,petsc_ier); CALL VecAssemblyEnd(b2,petsc_ier)
  CALL VecAssemblyBegin(wrk,petsc_ier); CALL VecAssemblyEnd(wrk,petsc_ier)

  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatMultAdd(Asub(3),wrk,b2,b2,petsc_ier)

  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_INITIAL_MATRIX,Asub(3),petsc_ier)
  CALL MatScale(Asub(3),-1.d0,petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- rebuild matrices and vectors --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  CALL VecRestoreSubVector(b,isg(1),b1,petsc_ier);   CALL VecRestoreSubVector(b,isg(2),b2,petsc_ier)
  CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier);   CALL VecRestoreSubVector(q,isg(2),q2,petsc_ier)
  CALL VecDestroy(wrk,petsc_ier)

  !-- set previous solution vectors (TODO: proper multistep startup) --
  IF (time_dependent) THEN
    CALL VecDuplicate(q,scl,petsc_ier);  CALL VecDuplicate(q,q_scl,petsc_ier)
    CALL VecGetSubVector(scl,isg(1),q1,petsc_ier)
    CALL VecGetSubVector(scl,isg(2),q2,petsc_ier)
    CALL VecSet(q1,0.d0,petsc_ier)
    DO i = 1,num_elm(dim_cmplx+1)
      CALL VecSetValues(q2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        lcl_complex(dim_cmplx+1)%inv_hdg_star(i),INSERT_VALUES,petsc_ier)
    END DO
    CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
    CALL VecAssemblyBegin(q2,petsc_ier); CALL VecAssemblyEnd(q2,petsc_ier)
    CALL VecRestoreSubVector(scl,isg(1),q1,petsc_ier)
    CALL VecRestoreSubVector(scl,isg(2),q2,petsc_ier)

    CALL VecPointwiseMult(q_scl,q,scl,petsc_ier)

    DO i=1,time_marching%r
      CALL VecCopy(q_scl,q_step(i),petsc_ier)
    END DO
  END IF

  !-- clean up --
  DEALLOCATE(flx_indx)

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_RHS_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_LHS_diffusion2
!* SYNOPSIS
  SUBROUTINE setup_LHS_diffusion()
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

  !-- local variables --
  INTEGER         :: m          !> number of rows/columns in global solution variables
  INTEGER         :: i,j,offset !>
  REAL(KIND=PGMSiwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- predefine some useful variables --
  offset = glb_num_elm(dim_cmplx)
  m = glb_num_elm(dim_cmplx+1) + glb_num_elm(dim_cmplx)

  !-- system matrix - block 00 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(1),petsc_ier)
  CALL MatSetSizes(Asub(1),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(dim_cmplx),glb_num_elm(dim_cmplx),petsc_ier)
  CALL MatSetType(Asub(1), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(1),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetFromOptions(Asub(1),petsc_ier)

  !-- equations for diffusion's law --
  DO i = 1,num_elm(dim_cmplx)
    CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
      1,lcl_complex(dim_cmplx)%glb_indx(i)-1,lcl_complex(dim_cmplx)%hdg_star(i),&
      INSERT_VALUES,petsc_ier)
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
      -real(lcl_complex(dim_cmplx+1)%bndry(i)%sgn,PGMSiwp),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- system matrix - block 01 --
  CALL MatTranspose(Asub(3),MAT_INITIAL_MATRIX,Asub(2),petsc_ier)
  CALL MatScale(Asub(2),-1.d0,petsc_ier)
  CALL MatSetOption(Asub(2),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(2),petsc_ier)

  !-- system matrix - block 11 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(4),petsc_ier)
  CALL MatSetSizes(Asub(4),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(dim_cmplx+1),glb_num_elm(dim_cmplx+1),petsc_ier)
  CALL MatSetType(Asub(4), MATMPIAIJ,petsc_ier)
  IF (time_dependent) THEN
    CALL MatMPIAIJSetPreallocation(Asub(4),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)

    !-- set memory place holder --
    DO i = 1,num_elm(dim_cmplx+1)
      CALL MatSetValues(Asub(4),1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,small,INSERT_VALUES,petsc_ier)
    END DO
  ELSE
    CALL MatMPIAIJSetPreallocation(Asub(4),0,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  END IF
  CALL MatSetFromOptions(Asub(4),petsc_ier)
  CALL MatAssemblyBegin(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(Asub(4),MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- system matrix (LHS:A) --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  !CALL MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)
  CALL MatNestGetISs(A,isg,PETSC_NULL_IS,petsc_ier)

  !-- system Schur Complement --
  ! CALL MatCreate(MPI_COMM_WORLD,Sp,petsc_ier)
  ! CALL MatSetSizes(Sp,PETSC_DECIDE,PETSC_DECIDE,&
  !   glb_num_elm(dim_cmplx)+glb_num_elm(dim_cmplx+1),&
  !   glb_num_elm(dim_cmplx)+glb_num_elm(dim_cmplx+1),petsc_ier)
  ! CALL MatSetType(Sp, MATMPIAIJ,petsc_ier)
  ! CALL MatMPIAIJSetPreallocation(Sp,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,petsc_ier)
  ! CALL MatSetFromOptions(Sp,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_LHS_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/update_LHS_unsteady
!* SYNOPSIS
  SUBROUTINE update_LHS_unsteady()
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

  !-- local variables --
  INTEGER :: i,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- A11 --
  DO i = 1,num_elm(dim_cmplx+1)
    CALL MatSetValues(Asub(4),1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,&
      1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1, lcl_complex(dim_cmplx+1)%inv_hdg_star(i)* &
      time_marching%LHS_tvec(stage),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/update_LHS_unsteady
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/get_RHS_diffusion
!* SYNOPSIS
  SUBROUTINE get_RHS_diffusion()
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

  !-- local variables --
  INTEGER :: i,j,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: r1,r2      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- vectors --
  CALL VecGetSubVector(r,isg(1),r1,petsc_ier);   CALL VecGetSubVector(r,isg(2),r2,petsc_ier)

  if (rank>0) CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  !-- equations for diffusion's law --
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%num_cobndry(i)==0) CYCLE
    !-- contribution from primal faces --
    IF (lcl_complex(dim_cmplx)%bc_type(i)>=0) THEN
      flx = lcl_complex(dim_cmplx)%hdg_star(i)*lcl_complex(dim_cmplx)%prml_sol(i,1)
      DO j = 1,lcl_complex(dim_cmplx)%num_cobndry(i)
        flx = flx + lcl_complex(dim_cmplx)%cobndry(i)%sgn(j)* &
          lcl_complex(dim_cmplx+1)%dual_sol(lcl_complex(dim_cmplx)%cobndry(i)%indx(j),1)
      END DO
    ELSEIF (lcl_complex(dim_cmplx)%bc_type(i)==-1) THEN
      flx = 0.d0
    END IF

    !-- add contribution to (negative) residual --
    CALL VecSetValues(r1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-flx,&
      INSERT_VALUES,petsc_ier)
  END DO
  if (rank==0) CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  !-- continuity equations --
  DO i = 1,num_elm(dim_cmplx+1)
    !-- contribution from primal faces --
    flx = 0
    DO j = 1,lcl_complex(dim_cmplx+1)%num_bndry(i)
      flx = flx + lcl_complex(dim_cmplx+1)%bndry(i)%sgn(j)* &
        lcl_complex(dim_cmplx)%prml_sol(lcl_complex(dim_cmplx+1)%bndry(i)%indx(j),1)
    END DO

    !-- add contribution to (negative) residual --
    CALL VecSetValues(r2,1,lcl_complex(dim_cmplx+1)%glb_indx(i)-1,flx,&
      INSERT_VALUES,petsc_ier)
  END DO

  !-- assemble sub-vectors and reassemble full vector --
  CALL VecAssemblyBegin(r1,petsc_ier); CALL VecAssemblyEnd(r1,petsc_ier)
  CALL VecAssemblyBegin(r2,petsc_ier); CALL VecAssemblyEnd(r2,petsc_ier)
  CALL VecRestoreSubVector(r,isg(1),r1,petsc_ier);   CALL VecRestoreSubVector(r,isg(2),r2,petsc_ier)

  !-- add boundary conditions to (negative) residual --
  CALL VecAXPY(r,1.d0,b,petsc_ier)
  CALL VecAssemblyBegin(r,petsc_ier); CALL VecAssemblyEnd(r,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/get_RHS_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/get_RHS_DIMRK
!* SYNOPSIS
  SUBROUTINE get_RHS_DIMRK()
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

  !-- local variables --
  INTEGER :: i          !> loop and temporary indicies

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- store current stage (negative) spatial residual --
  !r_stage(time_marching%r_stage_indx(stage)) = r_spatial
  CALL VecCopy(r,r_stage(stage),petsc_ier);

  !-- add current (scaled) solution to full residual --
  !r = r + time_marching%LHS_tvec(stage)*q
  CALL VecPointwiseMult(q_scl,q,scl,petsc_ier)
  CALL VecAXPY(r,-time_marching%LHS_tvec(stage),q_scl,petsc_ier)

  !-- add past (scaled) step solutions to full residual --
  DO i=1,time_marching%r;   IF (time_marching%U(stage,i)==0) CYCLE
    !r = r + (-)time_marching%U(stage,i)*q_step(time_marching%q_step_indx(i))
    CALL VecAXPY(r,time_marching%U(stage,i),q_step(time_marching%q_step_indx(i)),petsc_ier)
  END DO

  !-- add (scaled) stage residuals to full residual --
  DO i=1,stage-1;   IF (time_marching%A(stage,i)==0) CYCLE
    !r = r + time_marching%A(stage,i)*r_stage(time_marching%r_stage_indx(i))
    CALL VecAXPY(r,time_marching%A(stage,i),r_stage(i),petsc_ier)
  END DO

  RETURN

  END SUBROUTINE
! diffusion_mod/get_RHS_DIMRK
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/get_RHS_unsteady
!* SYNOPSIS
  SUBROUTINE update_time_step()
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

  !-- local variables --
  INTEGER :: i          !> loop and temporary indicies

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IF (time_marching%stifflyAcc) THEN
    !-- update previous solution indices and vector --
    DO i=1,time_marching%r
      time_marching%q_step_indx(i) = time_marching%q_step_indx(i) - 1
      IF (time_marching%q_step_indx(i) < 1) time_marching%q_step_indx(i) = time_marching%r
    END DO

    CALL VecPointwiseMult(q_step(time_marching%q_step_indx(1)),q,scl,petsc_ier)
  ELSE
    !-- zero solution and add past (scaled) step solutions to full residual --
    CALL VecScale(q_scl,0.d0,petsc_ier)
    DO i=1,time_marching%r;   IF (time_marching%V(1,i)==0) CYCLE
      !q = q + time_marching%V(1,i)*q_step(time_marching%q_step_indx(i))
      CALL VecAXPY(q_scl,time_marching%V(1,i),q_step(time_marching%q_step_indx(i)),petsc_ier)
    END DO

    !-- add (scaled) stage residuals to full residual --
    DO i=1,time_marching%s;   IF (time_marching%B(1,i)==0) CYCLE
      !q = q + (-)time_marching%B(1,i)*r_stage(i)
      CALL VecAXPY(q_scl,time_marching%B(1,i),r_stage(i),petsc_ier)
    END DO

    !-- update previous solution indices and vector --
    DO i=1,time_marching%r
      time_marching%q_step_indx(i) = time_marching%q_step_indx(i) - 1
      IF (time_marching%q_step_indx(i) < 1) time_marching%q_step_indx(i) = time_marching%r
    END DO

    CALL VecCopy(q_scl,q_step(time_marching%q_step_indx(1)),petsc_ier)
    CALL VecPointwiseDivide(q,q_scl,scl,petsc_ier)
  END IF

  RETURN

  END SUBROUTINE
! diffusion_mod/update_time_step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/get_LHS_diffusion2
!* SYNOPSIS
  SUBROUTINE get_LHS_diffusion()
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

  !-- local variables --
  INTEGER         :: m          !> number of rows/columns in global solution variables
  INTEGER         :: i,j,offset !>
  REAL(KIND=PGMSiwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- predefine some useful variables --
  offset = glb_num_elm(dim_cmplx)
  m = glb_num_elm(dim_cmplx+1) + glb_num_elm(dim_cmplx)

  !-- system matrix - block 00 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(1),petsc_ier)
  CALL MatSetSizes(Asub(1),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(dim_cmplx),glb_num_elm(dim_cmplx),petsc_ier)
  CALL MatSetType(Asub(1), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(1),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetFromOptions(Asub(1),petsc_ier)

  !-- equations for diffusion's law --
  CALL RANDOM_SEED()
  DO i = 1,num_elm(dim_cmplx)
    CALL RANDOM_NUMBER(rnd)
    IF (lcl_complex(dim_cmplx)%num_cobndry(i)>0) &
      CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
        1,lcl_complex(dim_cmplx)%glb_indx(i)-1,lcl_complex(dim_cmplx)%hdg_star(i),INSERT_VALUES,petsc_ier)
      ! CALL MatSetValues(Asub(1),1,lcl_complex(dim_cmplx)%glb_indx(i)-1,&
      !   1,lcl_complex(dim_cmplx)%glb_indx(i)-1,-(1.d0+rnd_max*2.0d0*(rnd-0.5d0))*&
      !   1.d0/Dk*sign(1.d0,lcl_complex(dim_cmplx)%hdg_star(i))*&
      !   max(abs(lcl_complex(dim_cmplx)%hdg_star(i)),smallh),INSERT_VALUES,petsc_ier)
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
      real(lcl_complex(dim_cmplx+1)%bndry(i)%sgn,PGMSiwp),INSERT_VALUES,petsc_ier)
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

  RETURN

  END SUBROUTINE
! diffusion_mod/get_LHS_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/identify_threshold_crack
!* SYNOPSIS
  SUBROUTINE identify_threshold_crack(iter,exit_cond)
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
  INTEGER, INTENT(IN)   :: iter               !>
  LOGICAL, INTENT(OUT)  :: exit_cond          !>

  !-- local variables --
  INTEGER :: i,j                              !> loop and temporary indicies
  INTEGER :: type_bc                          !> type of boundary condition
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices
  LOGICAL  :: lcl_exit_cond = .FALSE.         !>
  Vec             :: q1,q2      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- get face sub solution vector --
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier)

  !-- set new cracks --
  num_flx = 0
  DO i = 1,num_elm(dim_cmplx)
    IF (lcl_complex(dim_cmplx)%bc_type(i) == 0 .AND. &
      ABS(lcl_complex(dim_cmplx)%prml_sol(i,1)/&
      lcl_complex(dim_cmplx)%prml_volume(i)) > crck_thrshld) THEN

      !-- impermeable boundary --
      lcl_complex(dim_cmplx)%prml_sol(i,1) = 0.d0
      CALL VecSetValues(q1,1,lcl_complex(dim_cmplx)%glb_indx(i)-1,0.d0,&
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
  CALL VecAssemblyBegin(q,petsc_ier); CALL VecAssemblyEnd(q,petsc_ier)

  !-- assemble solution vector  --
  CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
  CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier)
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
  ! CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  !-- clean up --
  DEALLOCATE(flx_indx)

  IF (num_flx==0) lcl_exit_cond = .TRUE.
  CALL MPI_ALLREDUCE(lcl_exit_cond,exit_cond,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/identify_threshold_crack
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/identify_maxFD_crack
!* SYNOPSIS
  SUBROUTINE identify_maxFD_crack(iter,exit_cond)
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
  REAL(KIND=PGMSiwp) :: glb_max_flx(3)      !>
  REAL(KIND=PGMSiwp) :: max_flx(3)      !>
  REAL(KIND=PGMSiwp) :: flx      !>
  REAL(KIND=PGMSiwp) ::  area      !>
  CHARACTER(LEN=slen) :: msg      !>
  INTEGER        :: status(MPI_STATUS_SIZE)        !> size of MPI comm buffer
  Vec             :: q1,q2      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- get face sub solution vector --
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier)

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
      lcl_complex(dim_cmplx)%prml_sol(max_indx(1),1) = 0.d0
      CALL VecSetValues(q1,1,max_indx(2),0.d0,INSERT_VALUES,petsc_ier)
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
    CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
    CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier)
    CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
    CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
    ! CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
    ! CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
    ! CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
    ! CALL MatScale(Asub(3),-1.d0,petsc_ier)
    ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

    exit_cond = .FALSE.
  ELSE
    exit_cond = .TRUE.
  END IF

  RETURN

  END SUBROUTINE
! diffusion_mod/identify_maxFD_crack
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/identify_random_crack
!* SYNOPSIS
  SUBROUTINE identify_random_crack(iter,exit_cond)
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
  REAL(KIND=PGMSiwp) ::  area      !>
  CHARACTER(LEN=slen) :: msg      !>
  INTEGER        :: status(MPI_STATUS_SIZE)        !> size of MPI comm buffer
  REAL(KIND=PGMSiwp)  :: rnd
  LOGICAL :: search
  Vec             :: q1,q2      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- get face sub solution vector --
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier)

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

        CALL VecSetValues(q1,1,indx-1,0.d0,INSERT_VALUES,petsc_ier)
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
  CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
  CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier)
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
  ! CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  exit_cond = .FALSE.

  RETURN

  END SUBROUTINE
! diffusion_mod/identify_random_crack
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/identify_minFD_crack
!* SYNOPSIS
  SUBROUTINE identify_minFD_crack(iter,exit_cond)
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
  REAL(KIND=PGMSiwp) :: glb_max_flx(3)      !>
  REAL(KIND=PGMSiwp) :: max_flx(3)      !>
  REAL(KIND=PGMSiwp) :: flx      !>
  REAL(KIND=PGMSiwp) ::  area      !>
  CHARACTER(LEN=slen) :: msg      !>
  INTEGER        :: status(MPI_STATUS_SIZE)        !> size of MPI comm buffer
  Vec             :: q1,q2      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- get face sub solution vector --
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier)

  !-- find new crack --
  max_flx = (/ large, dble(rank), 0.d0 /)
  max_indx = -1
  area = 0.d0
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

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,max_flx(1:2),2,MPI_2DOUBLE_PRECISION,&
    MPI_MINLOC,MPI_COMM_WORLD,ier)

  !-- impermeable boundary --
  IF (max_indx(2)==int(max_flx(2))) THEN
    CALL VecSetValues(q1,1,max_indx(2),0.d0,INSERT_VALUES,petsc_ier)
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

  IF (rank==root) THEN
    area = glb_max_flx(3)
    crck_area = crck_area + area
    WRITE(ulog_unit,*) iter, glb_max_flx(1), int(glb_max_flx(2)), area, &
    area**(3.d0/2.d0), crck_area, crck_area**(3.d0/2.d0)
    CALL FLUSH(ulog_unit)
  END IF

  !-- assemble solution vector  --
  CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
  CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier)
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatTranspose(Asub(2),MAT_REUSE_MATRIX,Asub(3),petsc_ier)
  ! CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  ! CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  exit_cond = .FALSE.

  IF (glb_max_flx(1)==large) exit_cond = .TRUE.

  RETURN

  END SUBROUTINE
! diffusion_mod/identify_minFD_crack
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* time_marching_mod/set_IC_diffusion
!* SYNOPSIS
  SUBROUTINE set_IC_diffusion()
!* PURPOSE
!*   Deallocate global structures/variables
!* SIDE EFFECTS
!*   Global structure (lcl_complex) deallocated
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Deallocate global structures/variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER         :: i
    REAL(KIND=PGMSiwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- deallocate primary elements --
    DO i = 1,num_elm(dim_cmplx)
      lcl_complex(dim_cmplx)%prml_sol(i,1) = &
        sum(blk_dir(1:dim_embbd)*&
        lcl_complex(dim_cmplx)%dual_dir(i,1:dim_embbd))*&
        lcl_complex(dim_cmplx)%prml_volume(i)
    END DO
    DO i = 1,num_elm(dim_cmplx+1)
      lcl_complex(dim_cmplx+1)%dual_sol(i,1) = &
        -sum(blk_dir(1:dim_embbd)*&
        lcl_complex(dim_cmplx+1)%centers(i,1:dim_embbd)) + ref_val
    END DO
    ! lcl_complex(dim_cmplx)%prml_sol = 0.d0
    ! lcl_complex(dim_cmplx+1)%dual_sol = 0.d0

    IF (rnd_max>0.d0) THEN
      CALL RANDOM_SEED()
      DO i = 1,num_elm(dim_cmplx)
        CALL RANDOM_NUMBER(rnd)
        lcl_complex(dim_cmplx)%hdg_star(i) = (1.d0+rnd_max*2.0d0*(rnd-0.5d0))*&
          1.d0/Dk*sign(1.d0,lcl_complex(dim_cmplx)%hdg_star(i))*&
          max(abs(lcl_complex(dim_cmplx)%hdg_star(i)),smalls)
      END DO
    ELSE
      DO i = 1,num_elm(dim_cmplx)
        lcl_complex(dim_cmplx)%hdg_star(i) = &
          1.d0/Dk*sign(1.d0,lcl_complex(dim_cmplx)%hdg_star(i))*&
          max(abs(lcl_complex(dim_cmplx)%hdg_star(i)),smalls)
          !max(lcl_complex(dim_cmplx)%hdg_star(i),smallh)
      END DO
    END IF

    RETURN

  END SUBROUTINE
! set_IC_diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* time_marching_mod/set_IC_diffusion_d
!* SYNOPSIS
  SUBROUTINE set_IC_diffusion_d()
!* PURPOSE
!*   Deallocate global structures/variables
!* SIDE EFFECTS
!*   Global structure (lcl_complex) deallocated
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Deallocate global structures/variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER         :: i
    REAL(KIND=PGMSiwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- deallocate primary elements --
    lcl_complex(1)%prml_sol = 0.d0
    lcl_complex(2)%dual_sol = 0.d0

    IF (rnd_max>0.d0) THEN
      CALL RANDOM_SEED()
      DO i = 1,num_elm(2)
        CALL RANDOM_NUMBER(rnd)
        lcl_complex(2)%inv_hdg_star(i) = (1.d0+rnd_max*2.0d0*(rnd-0.5d0))*&
          1.d0/Dk*sign(1.d0,lcl_complex(2)%inv_hdg_star(i))*&
          max(abs(lcl_complex(2)%inv_hdg_star(i)),smalls)
      END DO
    ELSE
      DO i = 1,num_elm(2)
        lcl_complex(2)%inv_hdg_star(i) = 1.d0/Dk*&
          sign(1.d0,lcl_complex(2)%inv_hdg_star(i))*&
          max(abs(lcl_complex(2)%inv_hdg_star(i)),smalls)
      END DO
    END IF

    RETURN

  END SUBROUTINE
! set_IC_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_RHS_diffusion_d
!* SYNOPSIS
  SUBROUTINE setup_RHS_diffusion_d()
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

  !-- local variables --
  INTEGER :: i,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- vectors --
  CALL MatGetOwnershipRange(A,m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(1) + glb_num_elm(2),q,petsc_ier)
  CALL VecSetFromOptions(q,petsc_ier);  CALL VecDuplicate(q,dq,petsc_ier)
  CALL VecDuplicate(q,r,petsc_ier);     CALL VecDuplicate(q,b,petsc_ier)
  IF (time_dependent) THEN
    IF (ALLOCATED(q_step)) DEALLOCATE(q_step)
    IF (ALLOCATED(r_stage)) DEALLOCATE(r_stage)

    ALLOCATE(q_step(time_marching%r),r_stage(time_marching%s))

    DO i=1,time_marching%r;  CALL VecDuplicate(q,q_step(i),petsc_ier);  END DO
    DO i=1,time_marching%s;  CALL VecDuplicate(q,r_stage(i),petsc_ier);  END DO
  END IF

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_RHS_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_RHS_diffusion_d
!* SYNOPSIS
  SUBROUTINE setup_IC_BC_diffusion_d()
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

  !-- local variables --
  INTEGER :: i,j,k,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !CALL read_bndry_cond_d()

  !-- vectors --
  CALL VecGetSubVector(b,isg(1),b1,petsc_ier);   CALL VecGetSubVector(b,isg(2),b2,petsc_ier)
  CALL VecGetSubVector(q,isg(1),q1,petsc_ier);   CALL VecGetSubVector(q,isg(2),q2,petsc_ier)

  CALL MatGetOwnershipRange(Asub(1),m,n,petsc_ier)
  CALL VecCreateMPI(MPI_COMM_WORLD,n-m,glb_num_elm(2),wrk,petsc_ier)
  CALL VecSetFromOptions(wrk,petsc_ier)

  !--  set known values for initial and boundary conditions (fluxes) --
  num_flx = 0
  ALLOCATE(flx_indx(num_elm(2)))
  DO i = 1,num_elm(2)
    type_bc = lcl_complex(2)%bc_type(i)
    IF (type_bc == -2) THEN
      !-- flux boundary --
      flx = -sum(blk_dir*lcl_complex(2)%prml_dir(i,1:3))*lcl_complex(2)%dual_volume(i)
      lcl_complex(2)%dual_sol(i,1) = flx
      CALL VecSetValues(q1,1,lcl_complex(2)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(wrk,1,lcl_complex(2)%glb_indx(i)-1,-flx,&
        INSERT_VALUES,petsc_ier)
      CALL VecSetValues(b1,1,lcl_complex(2)%glb_indx(i)-1,flx,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(2)%glb_indx(i)-1,&
        1,lcl_complex(2)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      lcl_complex(2)%inv_hdg_star(i) = 1.d0
      CALL MatSetValues(Asub(2),1,lcl_complex(2)%glb_indx(i)-1,&
        lcl_complex(2)%num_cobndry(i),lcl_complex(1)%glb_indx(&
        lcl_complex(2)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(2)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc == -1) THEN
      !-- impermeable boundary --
      lcl_complex(2)%dual_sol(i,1) = 0.d0
      CALL VecSetValues(q1,1,lcl_complex(2)%glb_indx(i)-1,0.d0,&
        INSERT_VALUES,petsc_ier)
      CALL MatSetValues(Asub(1),1,lcl_complex(2)%glb_indx(i)-1,&
        1,lcl_complex(2)%glb_indx(i)-1,1.d0,INSERT_VALUES,petsc_ier)
      lcl_complex(2)%inv_hdg_star(i) = 1.d0
      CALL MatSetValues(Asub(2),1,lcl_complex(2)%glb_indx(i)-1,&
        lcl_complex(2)%num_cobndry(i),lcl_complex(1)%glb_indx(&
        lcl_complex(2)%cobndry(i)%indx)-1,&
        (/ (0.d0, m=1,lcl_complex(2)%num_cobndry(i)) /),INSERT_VALUES,petsc_ier)
    ELSEIF (type_bc>0) THEN
      !-- pressure boundary condition --
      DO j=1,lcl_complex(2)%num_cobndry(i)
        indx = lcl_complex(2)%cobndry(i)%indx(j)
        DO k = 1,lcl_complex(1)%num_bndry(indx)
          IF (lcl_complex(1)%bndry(indx)%indx(k)==i) THEN

            CALL VecSetValues(b1,1,lcl_complex(2)%glb_indx(i)-1,&
              lcl_complex(1)%bndry(indx)%sgn(k)*bc_vals(type_bc),&
              INSERT_VALUES,petsc_ier)
          END IF
        END DO
      END DO
    ELSE
      CALL VecSetValues(q1,1,lcl_complex(2)%glb_indx(i)-1,&
        lcl_complex(2)%dual_sol(i,1),INSERT_VALUES,petsc_ier)
    END IF
  END DO

  !--  set known values for initial conditions (pressures) --
  DO i = 1,num_elm(1)
    CALL VecSetValues(q2,1,lcl_complex(1)%glb_indx(i)-1,&
      lcl_complex(1)%prml_sol(i,1),INSERT_VALUES,petsc_ier)
  END DO

  !-- assemble sub-vectors --
  CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
  CALL VecAssemblyBegin(q2,petsc_ier); CALL VecAssemblyEnd(q2,petsc_ier)
  CALL VecAssemblyBegin(b1,petsc_ier); CALL VecAssemblyEnd(b1,petsc_ier)
  CALL VecAssemblyBegin(b2,petsc_ier); CALL VecAssemblyEnd(b2,petsc_ier)
  CALL VecAssemblyBegin(wrk,petsc_ier); CALL VecAssemblyEnd(wrk,petsc_ier)

  !-- set RHS (b vector) from known values and remove relevant equations from LHS (A matrix) --
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatMultAdd(Asub(3),wrk,b2,b2,petsc_ier)

  CALL MatAssemblyBegin(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(2),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatTranspose(Asub(2),MAT_INITIAL_MATRIX,Asub(3),petsc_ier)
  CALL MatScale(Asub(3),-1.d0,petsc_ier)
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- rebuild matrices and vectors --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  CALL VecRestoreSubVector(b,isg(1),b1,petsc_ier);   CALL VecRestoreSubVector(b,isg(2),b2,petsc_ier)
  CALL VecRestoreSubVector(q,isg(1),q1,petsc_ier);   CALL VecRestoreSubVector(q,isg(2),q2,petsc_ier)
  CALL VecDestroy(wrk,petsc_ier)

  !-- set previous solution vectors (TODO: proper multistep startup) --
  IF (time_dependent) THEN
    CALL VecDuplicate(q,scl,petsc_ier);  CALL VecDuplicate(q,q_scl,petsc_ier)
    CALL VecGetSubVector(scl,isg(1),q1,petsc_ier)
    CALL VecGetSubVector(scl,isg(2),q2,petsc_ier)
    CALL VecSet(q1,0.d0,petsc_ier)
    DO i = 1,num_elm(1)
      CALL VecSetValues(q2,1,lcl_complex(1)%glb_indx(i)-1,&
        lcl_complex(1)%hdg_star(i),INSERT_VALUES,petsc_ier)
    END DO
    CALL VecAssemblyBegin(q1,petsc_ier); CALL VecAssemblyEnd(q1,petsc_ier)
    CALL VecAssemblyBegin(q2,petsc_ier); CALL VecAssemblyEnd(q2,petsc_ier)
    CALL VecRestoreSubVector(scl,isg(1),q1,petsc_ier)
    CALL VecRestoreSubVector(scl,isg(2),q2,petsc_ier)

    CALL VecPointwiseMult(q_scl,q,scl,petsc_ier)

    DO i=1,time_marching%r
      CALL VecCopy(q_scl,q_step(i),petsc_ier)
    END DO
  END IF

  !-- clean up --
  DEALLOCATE(flx_indx)

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_RHS_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/setup_LHS_diffusion_d
!* SYNOPSIS
  SUBROUTINE setup_LHS_diffusion_d()
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

  !-- local variables --
  INTEGER         :: m          !> number of rows/columns in global solution variables
  INTEGER         :: i,j,offset !>
  REAL(KIND=PGMSiwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- predefine some useful variables --
  offset = glb_num_elm(2)
  m = glb_num_elm(1) + glb_num_elm(2)

  !-- system matrix - block 00 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(1),petsc_ier)
  CALL MatSetSizes(Asub(1),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(2),glb_num_elm(2),petsc_ier)
  CALL MatSetType(Asub(1), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(1),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetFromOptions(Asub(1),petsc_ier)

  !-- equations for diffusion's law --
  DO i = 1,num_elm(2)
    CALL MatSetValues(Asub(1),1,lcl_complex(2)%glb_indx(i)-1,&
      1,lcl_complex(2)%glb_indx(i)-1,lcl_complex(2)%inv_hdg_star(i),&
      INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(Asub(1),MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- system matrix - block 10 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(3),petsc_ier)
  CALL MatSetSizes(Asub(3),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(1),glb_num_elm(2),petsc_ier)
  CALL MatSetType(Asub(3), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(3),4,PETSC_NULL_INTEGER,4,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetOption(Asub(3),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(3),petsc_ier)
  !-- continuity equations: prml edges corresponding to internal primal faces --
  DO i = 1,num_elm(1)
    CALL MatSetValues(Asub(3),1,lcl_complex(1)%glb_indx(i)-1,&
      lcl_complex(1)%num_bndry(i),&
      lcl_complex(2)%glb_indx(lcl_complex(1)%bndry(i)%indx)-1,&
      -real(lcl_complex(1)%bndry(i)%sgn,PGMSiwp),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- system matrix - block 01 --
  CALL MatTranspose(Asub(3),MAT_INITIAL_MATRIX,Asub(2),petsc_ier)
  CALL MatScale(Asub(2),-1.d0,petsc_ier)
  CALL MatSetOption(Asub(2),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(2),petsc_ier)

  !-- system matrix - block 11 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(4),petsc_ier)
  CALL MatSetSizes(Asub(4),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(1),glb_num_elm(1),petsc_ier)
  CALL MatSetType(Asub(4), MATMPIAIJ,petsc_ier)
  IF (time_dependent) THEN
    CALL MatMPIAIJSetPreallocation(Asub(4),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)

    !-- set memory place holder --
    DO i = 1,num_elm(1)
      CALL MatSetValues(Asub(4),1,lcl_complex(1)%glb_indx(i)-1,&
        1,lcl_complex(1)%glb_indx(i)-1,small,INSERT_VALUES,petsc_ier)
    END DO
  ELSE
    CALL MatMPIAIJSetPreallocation(Asub(4),0,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  END IF
  CALL MatSetFromOptions(Asub(4),petsc_ier)
  CALL MatAssemblyBegin(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(Asub(4),MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- system matrix (LHS:A) --
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)
  !CALL MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)
  CALL MatNestGetISs(A,isg,PETSC_NULL_IS,petsc_ier)

  !-- system Schur Complement --
  ! CALL MatCreate(MPI_COMM_WORLD,Sp,petsc_ier)
  ! CALL MatSetSizes(Sp,PETSC_DECIDE,PETSC_DECIDE,&
  !   glb_num_elm(2)+glb_num_elm(1),&
  !   glb_num_elm(2)+glb_num_elm(1),petsc_ier)
  ! CALL MatSetType(Sp, MATMPIAIJ,petsc_ier)
  ! CALL MatMPIAIJSetPreallocation(Sp,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,petsc_ier)
  ! CALL MatSetFromOptions(Sp,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/setup_LHS_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/update_LHS_unsteady_d
!* SYNOPSIS
  SUBROUTINE update_LHS_unsteady_d()
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

  !-- local variables --
  INTEGER :: i,indx                           !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: b1,b2,q1,q2,wrk      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- A11 --
  DO i = 1,num_elm(1)
    CALL MatSetValues(Asub(4),1,lcl_complex(1)%glb_indx(i)-1,&
      1,lcl_complex(1)%glb_indx(i)-1, lcl_complex(1)%hdg_star(i)* &
      time_marching%LHS_tvec(stage),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier);  CALL MatAssemblyEnd(Asub(4),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatCreateNest(MPI_COMM_WORLD,2,PETSC_NULL_IS,2,PETSC_NULL_IS,Asub,A,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/update_LHS_unsteady_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/get_RHS_diffusion_d
!* SYNOPSIS
  SUBROUTINE get_RHS_diffusion_d()
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

  !-- local variables --
  INTEGER :: i,j,indx                         !> loop and temporary indicies
  INTEGER :: offset                           !> vector offset for volumes (after faces)

  INTEGER :: type_bc                          !> type of boundary condition
  LOGICAL :: use_pref,l_use_pref              !> use reference pressure (if no pressure bc used)
  INTEGER :: num_flx                          !> number of flux boundary faces
  INTEGER, ALLOCATABLE :: flx_indx(:)         !> flx boundary condition indices

  CHARACTER(LEN=slen)   :: fname              !> file name

  INTEGER         :: m,n                      !>
  REAL(KIND=PGMSiwp)  :: flx                      !>
  Vec             :: r1,r2      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- vectors --
  CALL VecGetSubVector(r,isg(1),r1,petsc_ier);   CALL VecGetSubVector(r,isg(2),r2,petsc_ier)

  if (rank>0) CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  !-- equations for diffusion's law --
  DO i = 1,num_elm(2)
    IF (lcl_complex(2)%num_cobndry(i)==0) CYCLE
    !-- contribution from primal faces --

    IF (lcl_complex(2)%bc_type(i)>=0) THEN
      flx = lcl_complex(2)%inv_hdg_star(i)*lcl_complex(2)%dual_sol(i,1)
      DO j = 1,lcl_complex(2)%num_cobndry(i)
        flx = flx + lcl_complex(2)%cobndry(i)%sgn(j)* &
          lcl_complex(1)%prml_sol(lcl_complex(2)%cobndry(i)%indx(j),1)
      END DO
    ELSEIF (lcl_complex(2)%bc_type(i)==-1) THEN
      flx = 0.d0
    END IF

    !-- add contribution to (negative) residual --
    CALL VecSetValues(r1,1,lcl_complex(2)%glb_indx(i)-1,-flx,&
      INSERT_VALUES,petsc_ier)
  END DO
  if (rank==0) CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  !-- continuity equations --
  DO i = 1,num_elm(1)
    !-- contribution from primal faces --
    flx = 0
    DO j = 1,lcl_complex(1)%num_bndry(i)
      flx = flx + lcl_complex(1)%bndry(i)%sgn(j)* &
        lcl_complex(2)%dual_sol(lcl_complex(1)%bndry(i)%indx(j),1)
    END DO

    !-- add contribution to (negative) residual --
    CALL VecSetValues(r2,1,lcl_complex(1)%glb_indx(i)-1,flx,&
      INSERT_VALUES,petsc_ier)
  END DO

  !-- assemble sub-vectors and reassemble full vector --
  CALL VecAssemblyBegin(r1,petsc_ier); CALL VecAssemblyEnd(r1,petsc_ier)
  CALL VecAssemblyBegin(r2,petsc_ier); CALL VecAssemblyEnd(r2,petsc_ier)
  CALL VecRestoreSubVector(r,isg(1),r1,petsc_ier);   CALL VecRestoreSubVector(r,isg(2),r2,petsc_ier)

  !-- add boundary conditions to (negative) residual --
  CALL VecAXPY(r,1.d0,b,petsc_ier)
  CALL VecAssemblyBegin(r,petsc_ier); CALL VecAssemblyEnd(r,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/get_RHS_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* diffusion_mod/get_LHS_diffusion_d
!* SYNOPSIS
  SUBROUTINE get_LHS_diffusion_d()
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

  !-- local variables --
  INTEGER         :: m          !> number of rows/columns in global solution variables
  INTEGER         :: i,j,offset !>
  REAL(KIND=PGMSiwp)  :: rnd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !-- predefine some useful variables --
  offset = glb_num_elm(2)
  m = glb_num_elm(1) + glb_num_elm(2)

  !-- system matrix - block 00 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(1),petsc_ier)
  CALL MatSetSizes(Asub(1),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(2),glb_num_elm(2),petsc_ier)
  CALL MatSetType(Asub(1), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(1),1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetFromOptions(Asub(1),petsc_ier)

  !-- equations for diffusion's law --
  CALL RANDOM_SEED()
  DO i = 1,num_elm(2)
    CALL RANDOM_NUMBER(rnd)
    IF (lcl_complex(2)%num_cobndry(i)>0) &
      CALL MatSetValues(Asub(1),1,lcl_complex(2)%glb_indx(i)-1,&
        1,lcl_complex(2)%glb_indx(i)-1,lcl_complex(2)%inv_hdg_star(i),INSERT_VALUES,petsc_ier)
      ! CALL MatSetValues(Asub(1),1,lcl_complex(2)%glb_indx(i)-1,&
      !   1,lcl_complex(2)%glb_indx(i)-1,-(1.d0+rnd_max*2.0d0*(rnd-0.5d0))*&
      !   1.d0/Dk*sign(1.d0,lcl_complex(2)%inv_hdg_star(i))*&
      !   max(abs(lcl_complex(2)%inv_hdg_star(i)),smallh),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(1),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatSetOption(Asub(1),MAT_SYMMETRIC,PETSC_TRUE,petsc_ier)

  !-- system matrix - block 10 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(3),petsc_ier)
  CALL MatSetSizes(Asub(3),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(1),glb_num_elm(2),petsc_ier)
  CALL MatSetType(Asub(3), MATMPIAIJ,petsc_ier)
  CALL MatMPIAIJSetPreallocation(Asub(3),4,PETSC_NULL_INTEGER,4,PETSC_NULL_INTEGER,petsc_ier)
  CALL MatSetOption(Asub(3),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(3),petsc_ier)
  !-- continuity equations: prml edges corresponding to internal primal faces --
  DO i = 1,num_elm(1)
    CALL MatSetValues(Asub(3),1,lcl_complex(1)%glb_indx(i)-1,&
      lcl_complex(1)%num_bndry(i),&
      lcl_complex(2)%glb_indx(lcl_complex(1)%bndry(i)%indx)-1,&
      real(lcl_complex(1)%bndry(i)%sgn,PGMSiwp),INSERT_VALUES,petsc_ier)
  END DO
  CALL MatAssemblyBegin(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)
  CALL MatAssemblyEnd(Asub(3),MAT_FINAL_ASSEMBLY,petsc_ier)

  !-- system matrix - block 01 --
  CALL MatTranspose(Asub(3),MAT_INITIAL_MATRIX,Asub(2),petsc_ier)
  CALL MatSetOption(Asub(2),MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,petsc_ier)
  CALL MatSetFromOptions(Asub(2),petsc_ier)

  !-- system matrix - block 11 --
  CALL MatCreate(MPI_COMM_WORLD,Asub(4),petsc_ier)
  CALL MatSetSizes(Asub(4),PETSC_DECIDE,PETSC_DECIDE,glb_num_elm(1),glb_num_elm(1),petsc_ier)
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
  !   glb_num_elm(2)+glb_num_elm(1),&
  !   glb_num_elm(2)+glb_num_elm(1),petsc_ier)
  ! CALL MatSetType(Sp, MATMPIAIJ,petsc_ier)
  ! CALL MatMPIAIJSetPreallocation(Sp,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,petsc_ier)
  ! CALL MatSetFromOptions(Sp,petsc_ier)

  RETURN

  END SUBROUTINE
! diffusion_mod/get_LHS_diffusion_d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! diffusion_mod
!===============================================================================
!
