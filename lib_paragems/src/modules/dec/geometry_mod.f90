!
!===============================================================================
!-- DEC Module
!> Module containing routines for performing DEC operations
!===============================================================================
!/****/h* modules|dec/geometry_mod
!* SYNOPSIS
MODULE geometry_mod
!* PURPOSE
!*   Module containing routines for performing DEC operations
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 ???
!*   math_mod                basic math functions
!* CONTAINS
!*   Subroutine              Purpose
!*   initialise_geo          initialise geometric quantities for the given mesh
!*   get_lcl_node_indx       create map between global and local node indices
!*   calc_circumcenters      compute circumcenter of given elements
!*   calc_prml_sgnd_vlm      compute signed volume of primal elements
!*   calc_prml_unsgnd_vlm    compute unsigned volume of primal elements
!*   calc_dual_vlm           compute volume for dual elements of all geometric order
!*   calc_dual_vlm_i         recursively add to dual volume calculation
!*   add_points_i            recursively adds points for dual volume calculation
!*   calc_unsgnd_vlm         compute unsigned volume for given set of points
!*   exchange_dual_vlm       exchange external dual volumes between adjacent processes
!*   exchange_dual_dir       exchange external dual edge directions between adjacent processes
!*   calc_prml_unsgnd_vlm    compute unsigned volume of primal elements
!*   calc_prml_dir           compute the unit direction of primal edges
!*   calc_dual_dir           compute the unit direction of dual edges
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2020/09/16
!> Module containing routines for performing DEC operations
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod
  USE math_mod
  USE dec_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/initialise_geo
!* SYNOPSIS
  SUBROUTINE initialise_geo(ki,km)
!* PURPOSE
!*   Initialise geometric quantities for the given mesh
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
!> Initialise geometric quantities for the given mesh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER                 :: ki,km               !> loop counters

    !-- local variables --
    INTEGER                 :: i,k               !> loop counters
    REAL(KIND=PGMSiwp)          :: minV,maxV,absminV !> min/max volumes
    CHARACTER(LEN=slen)     :: str               !> string for writing to log file

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- log element information --
    DO k = ki,km
      write(str,'(A,I1,A,I10)') '   - number of C',k-1,': ',glb_num_elm(k)
      CALL rootwrite(str,log_unit); END DO

    !-- get local node indices --
    CALL get_lcl_node_indx()

    !-- compute circumcenters --
    DO k = max(ki,2),km; CALL calc_circumcenters(k); END DO

    !-- compute signed primal volumes --
    IF (ki==1) THEN
      IF (.NOT. ALLOCATED(lcl_complex(1)%prml_volume)) &
        ALLOCATE(lcl_complex(1)%prml_volume(num_elm(1)))
        lcl_complex(1)%prml_volume = 1; END IF
    DO k=max(ki,2),km
      IF (k-1==dim_embbd) THEN; CALL calc_prml_sgnd_vlm(k)
      ELSE; CALL calc_prml_unsgnd_vlm(k); END IF; END DO

    !-- compute signed dual volumes --
    CALL calc_dual_vlm(ki)
    DO k =ki,min(km,dim_cmplx); CALL exchange_dual_vlm(k); END DO

    !-- assess maximumm and minimum volumes --
    DO k =ki,km
      CALL MPI_REDUCE(minval(lcl_complex(k)%prml_volume),minV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(minval(abs(lcl_complex(k)%prml_volume)),absminV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(maxval(lcl_complex(k)%prml_volume),maxV,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)
      IF (rank == root) &
        write(str,'(A,I1,A,ES10.3,A,ES10.3,A,ES10.3)') &
          '   - The min/abs(min)/max C',k-1,&
          ' volume is: ',minV,'/',absminV,'/',maxV
      CALL syncwrite(str,log_unit)

      CALL MPI_REDUCE(minval(lcl_complex(k)%dual_volume),minV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(minval(abs(lcl_complex(k)%dual_volume)),absminV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(maxval(lcl_complex(k)%dual_volume),maxV,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)
      IF (rank == root) &
        write(str,'(A,I1,A,ES10.3,A,ES10.3,A,ES10.3)') &
          '   - The min/abs(min)/max D',dim_cmplx+1-k,&
          ' volume is: ',minV,'/',absminV,'/',maxV
      CALL syncwrite(str,log_unit)
    END DO

    !-- compute primal and dual edge directions --
    if (ki<3 .and. km>1) CALL calc_prml_dir()
    if (ki<dim_cmplx+1 .and. km>dim_cmplx-1) CALL calc_dual_dir()

    RETURN

  END SUBROUTINE
! geometry_mod/initialise_geo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/get_lcl_node_indx
!* SYNOPSIS
  SUBROUTINE get_lcl_node_indx()
!* PURPOSE
!*   Create map between global and local node indices
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
!> Create map between global and local node indices
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER                 :: i,k,m,n        !> loop counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate array for map and assign default value --
    DO k = 2, dim_cmplx+1
      IF (.NOT. ALLOCATED(lcl_complex(k)%lcl_node_indx)) &
        ALLOCATE(lcl_complex(k)%lcl_node_indx(num_elm(k),k))
      lcl_complex(k)%lcl_node_indx = -1
    END DO

    !-- build map --
    DO k = 2, dim_cmplx+1     !-- simplicial order --
      DO i=1,num_elm(1)       !-- local node index --
        DO n = 1,k            !-- element node index --
          DO m=1,num_elm(k)   !-- element index --
            IF (lcl_complex(k)%node_indx(m,n) == lcl_complex(1)%glb_indx(i)) &
              lcl_complex(k)%lcl_node_indx(m,n) = i
          END DO
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE
! geometry_mod/get_lcl_node_indx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_circumcenters
!* SYNOPSIS
  SUBROUTINE calc_circumcenters(k)
!* PURPOSE
!*   Compute circumcenter of given elements
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
!> Compute circumcenter of given elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER                     :: k                !> simplicial order

    !-- local variables --
    INTEGER                     :: kp               !> simplicial order plus one
    INTEGER                     :: i,j              !> loop indices
    REAL(KIND=PGMSiwp), ALLOCATABLE :: A(:,:),b(:)      !> solution variables
    REAL(KIND=PGMSiwp), ALLOCATABLE :: pts(:,:)         !> bounding points
    REAL(KIND=PGMSiwp), ALLOCATABLE :: ipiv(:),work(:)  !> work arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate centers and temp. variables --
    kp=k+1
    IF (.NOT. ALLOCATED(lcl_complex(k)%centers)) &
      ALLOCATE(lcl_complex(k)%centers(num_elm(k),dim_embbd))
    ALLOCATE(pts(k,dim_embbd),A(kp,kp),b(kp),ipiv(kp),work(kp))

    !-- compute circumcenters for all elements --
    !-- loop through all elements of this order --
    DO i = 1, num_elm(k)
      !-- get bounding points of the element --
      DO j = 1,k
        pts(j,:) = lcl_complex(1)%centers(lcl_complex(k)%lcl_node_indx(i,j),:)
      END DO

      !-- setup system to compute circumcenter --
      A(1:k,1:k) = 2*MATMUL(pts,TRANSPOSE(pts))
      A(1:k,kp) = 1;  A(kp,1:k) = 1;  A(kp,kp) = 0
      b(1:k) = sum(pts*pts,DIM=2)
      b(kp) = 1

      !-- solve system to compute circumcenter --
      CALL DSYSV('L',kp,1,A,kp,ipiv,b,kp,work,kp,ier)

      !-- assign circumcenter --
      lcl_complex(k)%centers(i,:) = MATMUL(b(1:k),pts)
    END DO

    !-- clean up --
    DEALLOCATE(pts,A,b,ipiv,work)

    RETURN

  END SUBROUTINE
! geometry_mod/calc_circumcenters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_circumcenters_WBC
!* SYNOPSIS
  SUBROUTINE calc_circumcenters_WBC(k)
!* PURPOSE
!*   Compute circumcenter of given elements and store associated barycentric
!*   coordinates
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
!> Compute circumcenter of given elements and store associated barycentric
!>   coordinates
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER                     :: k                !> simplicial order

    !-- local variables --
    INTEGER                     :: kp               !> simplicial order plus one
    INTEGER                     :: i,j              !> loop indices
    REAL(KIND=PGMSiwp), ALLOCATABLE :: A(:,:),b(:)      !> solution variables
    REAL(KIND=PGMSiwp), ALLOCATABLE :: pts(:,:)         !> bounding points
    REAL(KIND=PGMSiwp), ALLOCATABLE :: ipiv(:),work(:)  !> work arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate centers and temp. variables --
    kp=k+1
    IF (.NOT. ALLOCATED(lcl_complex(k)%centers)) &
      ALLOCATE(lcl_complex(k)%centers(num_elm(k),dim_embbd), &
               lcl_complex(k)%b_coord(num_elm(k),k))
    ALLOCATE(pts(k,dim_embbd),A(kp,kp),b(kp),ipiv(kp),work(kp))

    !-- compute circumcenters for all elements --
    !-- loop through all elements of this order --
    DO i = 1, num_elm(k)
      !-- get bounding points of the element --
      DO j = 1,k
        pts(j,:) = lcl_complex(1)%centers(lcl_complex(k)%lcl_node_indx(i,j),:)
      END DO

      !-- setup system to compute circumcenter --
      A(1:k,1:k) = 2*MATMUL(pts,TRANSPOSE(pts))
      A(1:k,kp) = 1;  A(kp,1:k) = 1;  A(kp,kp) = 0
      b(1:k) = sum(pts*pts,DIM=2)
      b(kp) = 1

      !-- solve system to compute circumcenter --
      CALL DSYSV('L',kp,1,A,kp,ipiv,b,kp,work,kp,ier)

      !-- assign circumcenter --
      lcl_complex(k)%centers(i,:) = MATMUL(b(1:k),pts)
      lcl_complex(k)%b_coord(i,:) = b(1:k)
    END DO

    !-- clean up --
    DEALLOCATE(pts,A,b,ipiv,work)

    RETURN

  END SUBROUTINE
! geometry_mod/calc_circumcenters_WBC
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_prml_sgnd_vlm
!* SYNOPSIS
  SUBROUTINE calc_prml_sgnd_vlm(k)
!* PURPOSE
!*   Compute signed volume of primal elements
!* INPUTS
!*   Name                    Description
!*   k
!* OUTPUTS
!*   Name                    Description
!*   k
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
!> Compute signed volume of primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER                     :: k            !> simplicial order

    !-- local variables --
    INTEGER                     :: i,j          !> loop counters
    INTEGER                     :: fac          !> factorial
    REAL(KIND=PGMSiwp), ALLOCATABLE :: A(:,:)       !> solution variable
    REAL(KIND=PGMSiwp), ALLOCATABLE :: pts(:)       !> bounding points

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate volumes and temp. variables --
    IF (.NOT. ALLOCATED(lcl_complex(k)%prml_volume)) &
      ALLOCATE(lcl_complex(k)%prml_volume(num_elm(k)))
    ALLOCATE(pts(dim_embbd),A(dim_embbd,dim_embbd))

    !-- compute factorial of embedding dimension --
    fac = 1; DO i=2,dim_embbd; fac = fac*i; END DO

    !-- compute signed volumes --
    !-- loop through all elements --
    DO i = 1, num_elm(k)
      !-- get a bounding point of the element --
      pts(:) = lcl_complex(1)%centers(lcl_complex(k)%lcl_node_indx(i,1),:)

      !-- build system matrix from bounding points of the element --
      DO j = 1,k-1
        A(j,:) = lcl_complex(1)%centers(lcl_complex(k)%lcl_node_indx(i,j+1),:) - pts
      END DO

      !-- compute and assign signed volume --
      lcl_complex(k)%prml_volume(i) = lcl_complex(k)%orientation(i)*RCSV_determinant(A,dim_embbd)/fac
    END DO

    !-- clean up --
    DEALLOCATE(pts,A); RETURN

  END SUBROUTINE
! geometry_mod/calc_prml_sgnd_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_prml_unsgnd_vlm
!* SYNOPSIS
  SUBROUTINE calc_prml_unsgnd_vlm(k)
!* PURPOSE
!*   Compute unsigned volume of primal elements
!* INPUTS
!*   Name                    Description
!*   k
!* OUTPUTS
!*   Name                    Description
!*   k
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
!> Compute unsigned volume of primal elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER                     :: k            !> loop counters

    !-- local variables --
    INTEGER                     :: i,j          !> loop counters
    INTEGER                     :: fac          !> factorial
    REAL(KIND=PGMSiwp), ALLOCATABLE :: A(:,:)       !> solution variable
    REAL(KIND=PGMSiwp), ALLOCATABLE :: pts(:)       !> bounding points

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate volumes and temp. variables --
    IF (.NOT. ALLOCATED(lcl_complex(k)%prml_volume)) &
      ALLOCATE(lcl_complex(k)%prml_volume(num_elm(k)))
    ALLOCATE(pts(dim_embbd),A(k-1,dim_embbd))

    !-- compute necessary factorial --
    fac = 1; DO i=2,k-1; fac = fac*i; END DO

    !-- compute unsigned volumes --
    !-- loop for all elements --
    DO i = 1, num_elm(k)
      !-- get a bounding point of the element --
      pts(:) = lcl_complex(1)%centers(lcl_complex(k)%lcl_node_indx(i,1),:)

      !-- build system matrix from remaining bounding points of the element --
      DO j = 1,k-1
        A(j,:) = lcl_complex(1)%centers(lcl_complex(k)%lcl_node_indx(i,j+1),:) - pts
      END DO

      !-- compute and assign unsigned volume --
      lcl_complex(k)%prml_volume(i) = SQRT(ABS(RCSV_determinant(MATMUL(A,TRANSPOSE(A)),k-1)))/fac
    END DO

    !-- clean up --
    DEALLOCATE(pts,A)

    RETURN

  END SUBROUTINE
! geometry_mod/calc_prml_unsgnd_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_dual_vlm
!* SYNOPSIS
  SUBROUTINE calc_dual_vlm(ki)
!* PURPOSE
!*   Compute volume for dual elements of all geometric order
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
!> Compute volume for dual elements of all geometric order
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER                 :: ki               !> loop counters

    !-- local variables --
    INTEGER                 :: k                !> simplicial order
    INTEGER                 :: i,indx           !> loop counters
    REAL(KIND=PGMSiwp), ALLOCATABLE :: pts(:,:)     !> bounding points
    REAL(KIND=PGMSiwp), ALLOCATABLE :: sgn(:)       !> sign of volume elements

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate volumes and temp. variables --
    DO k = ki,dim_cmplx+1
      IF (.NOT. ALLOCATED(lcl_complex(k)%dual_volume)) &
        ALLOCATE(lcl_complex(k)%dual_volume(num_elm(k)))
      lcl_complex(k)%dual_volume = 0.d0
    END DO
    ALLOCATE(pts(dim_cmplx+1,dim_embbd),sgn(dim_cmplx))
    sgn = 1.d0

    !-- start recursive computation of dual volumes on the interior --
    k = dim_cmplx+1
    DO i = 1, num_elm(k)
      CALL calc_dual_vlm_i(pts,sgn,i,i,k,ki)
    END DO

    !-- clean up --
    DEALLOCATE(pts,sgn)

    RETURN

  END SUBROUTINE
! geometry_mod/calc_dual_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_dual_vlm_i
!* SYNOPSIS
  RECURSIVE SUBROUTINE calc_dual_vlm_i(pts,sgn,indx,p_indx,k,ki)
!* PURPOSE
!*   Recursively add to dual volume calculation
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
!> Recursively add to dual volume calculation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)           :: k,ki         !> simplicial order
    INTEGER, INTENT(IN)           :: indx,p_indx  !> element and parent index
    REAL(KIND=PGMSiwp), INTENT(INOUT) :: pts(:,:)     !> bounding points
    REAL(KIND=PGMSiwp), INTENT(INOUT) :: sgn(:)       !> sign of volume elements

    !-- local variables --
    INTEGER                       :: i            !> loop counters
    INTEGER                       :: t_indx       !> tmp index
    INTEGER                       :: n            !> number of points

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- add to array of bounding points --
    pts(k,:) = lcl_complex(k)%centers(indx,:)

    !-- get sign --
    IF (k>ki .AND. k<=dim_cmplx) THEN
      DO i = 1,k+1
        t_indx = lcl_complex(k+1)% node_indx(p_indx,i)
        IF (.NOT. ANY(lcl_complex(k)% node_indx(indx,:)==t_indx)) THEN
          pts(k-1,:) = lcl_complex(1)%centers(lcl_complex(k+1)%lcl_node_indx(p_indx,i),:)
          sgn(k-1) = sign(1.d0,dot_product(pts(k+1,:)-pts(k,:),pts(k-1,:)-pts(k,:)))
          EXIT
        END IF
      END DO
    END IF

    !-- compute number of points to be used in volume calculation --
    n = dim_cmplx + 1 - k

    !-- add to k^th order dual volume
    lcl_complex(k)%dual_volume(indx) = lcl_complex(k)%dual_volume(indx) + &
      product(sgn)*calc_unsgnd_vlm(pts(k:,:),n)

    !-- recursively add to dual volumes from lower order elements --
    IF (k>ki) THEN
      DO i = 1, lcl_complex(k)% num_bndry(indx)
        CALL calc_dual_vlm_i(pts,sgn,lcl_complex(k)% bndry(indx)% indx(i),indx,k-1,ki)
      END DO
    END IF

    IF (k>ki .AND. k<=dim_cmplx) sgn(k-1) = 1.d0

    RETURN

  END SUBROUTINE
! geometry_mod/calc_dual_vlm_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f* geometry_mod/calc_unsgnd_vlm
!* SYNOPSIS
  FUNCTION calc_unsgnd_vlm(pts,n)
!* PURPOSE
!*   Compute unsigned volume for given set of points
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
!> Compute unsigned volume for given set of points
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN)  :: pts(:,:)         !> bounding points
    INTEGER, INTENT(IN)         :: n                !> # points to use

    !-- local variables --
    INTEGER                     :: i                !> loop counters
    INTEGER                     :: fac              !> factorial
    REAL(KIND=PGMSiwp), ALLOCATABLE :: A(:,:)           !> solution variable
    REAL(KIND=PGMSiwp)              :: calc_unsgnd_vlm  !> function value

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IF (n==0) THEN
      !-- assign default value when no points are used --
      calc_unsgnd_vlm = 1.d0
    ELSE
      !-- allocate temp. variables --
      ALLOCATE(A(n,dim_embbd))

      !-- compute necessary factorial and build system matrix --
      fac = 1
      DO i=1,n
        fac = fac*i
        A(i,:) = pts(i+1,:) - pts(1,:)
      END DO

      !-- compute and assign unsigned volume --
      calc_unsgnd_vlm = SQRT(ABS(RCSV_determinant(MATMUL(A,TRANSPOSE(A)),n)))/fac

      !-- clean up --
      DEALLOCATE(A)
    END IF

    RETURN

  END FUNCTION
! geometry_mod/calc_unsgnd_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/exchange_dual_vlm
!* SYNOPSIS
  SUBROUTINE exchange_dual_vlm(k)
!* PURPOSE
!*   Exchange external dual volumes between adjacent processes
!* INPUTS
!*   Name                    Description
!*   k
!* OUTPUTS
!*   Name                    Description
!*   k
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/22: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/22
!> Exchange external dual volumes between adjacent processes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER,INTENT(IN)     :: k                !> simplicial order

    !-- local variables --
    INTEGER               :: kp                 !> simplicial order plus one
    INTEGER               :: cnt                !> counter
    INTEGER               :: i                  !> loop index

    INTEGER               :: junk               !> junk comm variable
    INTEGER               :: buffer_size        !> size of comm buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE  :: sbuffer(:)  !> send buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE  :: rbuffer(:)  !> receive buffer
    INTEGER, ALLOCATABLE  :: req(:)             !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE  :: status(:,:)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- allocate and pack send buffer --
    kp = k+1
    ALLOCATE(sbuffer(lcl_complex(kp)%num_send))
    DO i = 1,lcl_complex(kp)%num_send
      sbuffer(i) = lcl_complex(k)%dual_volume(lcl_complex(kp)%send_indx(i,2))
    END DO

    !-- send data to adjacent processes --
    cnt = 1
    DO i = 1,num_adj_proc
      IF (num_send(k,i)==0) CYCLE
      buffer_size = num_send(k,i)
      CALL MPI_ISEND(sbuffer(cnt:cnt+buffer_size-1),buffer_size,&
        MPI_DOUBLE_PRECISION,adj_proc(i),0,MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)
      cnt = cnt + buffer_size
    END DO

    !-- allocate recieve buffer and receive data from adjacent processes --
    ALLOCATE(rbuffer(lcl_complex(kp)%num_recv),&
      req(num_adj_proc),status(MPI_STATUS_SIZE,num_adj_proc))
    cnt = 1
    DO i = 1,num_adj_proc
      IF (num_recv(k,i)==0) THEN
        req(i) = MPI_REQUEST_NULL
        CYCLE
      END IF
      buffer_size = num_recv(k,i)
      CALL MPI_IRECV(rbuffer(cnt:cnt+buffer_size-1),buffer_size,&
        MPI_DOUBLE_PRECISION,adj_proc(i),0,MPI_COMM_WORLD,req(i),ier)
      cnt = cnt + buffer_size
    END DO

    !-- wait for communication to finish, then unpack receive buffer --
    CALL MPI_WAITALL(num_adj_proc,req,status,ier)
    DO i = 1,lcl_complex(kp)%num_recv
      lcl_complex(k)%dual_volume(lcl_complex(kp)%recv_indx(i,2)) = rbuffer(i)
    END DO

    !-- clean up --
    DEALLOCATE(req,status)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|exchange_dual_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/exchange_dual_dir
!* SYNOPSIS
  SUBROUTINE exchange_dual_dir(k)
!* PURPOSE
!*   Exchange external dual edge directions between adjacent processes
!* INPUTS
!*   Name                    Description
!*   k
!* OUTPUTS
!*   Name                    Description
!*   k
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/22: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/22
!> Exchange external dual edge directions between adjacent processes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER,INTENT(IN)    :: k                  !> simplicial order

    !-- local variables --
    INTEGER               :: kp                 !> simplicial order
    INTEGER               :: cnt                !> counter
    INTEGER               :: i                  !> loop index

    INTEGER               :: junk               !> junk comm variable
    INTEGER               :: buffer_size        !> size of comm buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE  :: sbuffer(:)  !> send buffer
    REAL(KIND=PGMSiwp), ALLOCATABLE  :: rbuffer(:)  !> receive buffer
    INTEGER, ALLOCATABLE  :: req(:)             !> request variable for non-blocking comms
    INTEGER, ALLOCATABLE  :: status(:,:)        !> size of MPI comm buffer

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- allocate and pack send buffer --
    kp = k+1
    ALLOCATE(sbuffer(dim_embbd*lcl_complex(kp)%num_send))
    DO i = 1,lcl_complex(kp)%num_send
      sbuffer((i-1)*dim_embbd+1:i*dim_embbd) = &
        lcl_complex(k)%dual_dir(lcl_complex(kp)%send_indx(i,2),:)
    END DO

    !-- send data to adjacent processes --
    cnt = 1
    DO i = 1,num_adj_proc
      IF (num_send(k,i)==0) CYCLE
      buffer_size = dim_embbd*num_send(k,i)
      CALL MPI_ISEND(sbuffer(cnt:cnt+buffer_size-1),buffer_size,&
        MPI_DOUBLE_PRECISION,adj_proc(i),0,MPI_COMM_WORLD,junk,ier)
      CALL MPI_REQUEST_FREE(junk,ier)
      cnt = cnt + buffer_size
    END DO

    !-- allocate recieve buffer and receive data from adjacent processes --
    ALLOCATE(&
      req(num_adj_proc),&
      status(MPI_STATUS_SIZE,num_adj_proc),&
      rbuffer(dim_embbd*lcl_complex(kp)%num_recv))
    cnt = 1
    DO i = 1,num_adj_proc
      IF (num_recv(k,i)==0) THEN
        req(i) = MPI_REQUEST_NULL
        CYCLE
      END IF
      buffer_size = dim_embbd*num_recv(k,i)
      CALL MPI_IRECV(rbuffer(cnt:cnt+buffer_size-1),buffer_size,&
        MPI_DOUBLE_PRECISION,adj_proc(i),0,MPI_COMM_WORLD,req(i),ier)
      cnt = cnt + buffer_size
    END DO

    !-- wait for communication to finish, then unpack receive buffer --
    CALL MPI_WAITALL(num_adj_proc,req,status,ier)
    DO i = 1,lcl_complex(kp)%num_recv
      lcl_complex(k)%dual_dir(lcl_complex(kp)%recv_indx(i,2),:) = &
        rbuffer((i-1)*dim_embbd+1:i*dim_embbd)
    END DO

    !-- clean up --
    DEALLOCATE(req,status)

    !-- make sure all comms have completed --
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE
! io_mod|exchange_dual_dir
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_prml_dir
!* SYNOPSIS
  SUBROUTINE calc_prml_dir()
!* PURPOSE
!*   Compute the unit direction of primal edges
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
!> Compute the unit direction of primal edges
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER           :: i        !> loop counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate direction variables --
    IF (.NOT. ALLOCATED(lcl_complex(2)%prml_dir)) &
      ALLOCATE(lcl_complex(2)%prml_dir(num_elm(2),dim_embbd))

    !-- compute edge direction --
    DO i = 1,num_elm(2)
      !-- add contribution from end-points (primal vertices) --
      lcl_complex(2)%prml_dir(i,:) = &
        lcl_complex(1)%centers(lcl_complex(2)%lcl_node_indx(i,1),:) - &
        lcl_complex(1)%centers(lcl_complex(2)%lcl_node_indx(i,2),:)

      !-- normalise --
      lcl_complex(2)%prml_dir(i,:) = lcl_complex(2)%prml_dir(i,:) / &
        max(sqrt(dot_product(lcl_complex(2)%prml_dir(i,:),lcl_complex(2)%prml_dir(i,:))),small)
    END DO

    RETURN

  END SUBROUTINE
! geometry_mod/calc_prml_dir
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_dual_dir
!* SYNOPSIS
  SUBROUTINE calc_dual_dir()
!* PURPOSE
!*   Compute the unit direction of dual edges
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
!> Compute the unit direction of dual edges
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER           :: i,j        !> loop counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate direction variables and zero --
    IF (.NOT. ALLOCATED(lcl_complex(dim_cmplx)%dual_dir)) &
      ALLOCATE(lcl_complex(dim_cmplx)%dual_dir(num_elm(dim_cmplx),dim_embbd))
    lcl_complex(dim_cmplx)%dual_dir = 0.d0

    !-- compute edge direction for each local dual edge --
    DO i = 1,num_elm(dim_cmplx)
      !-- add contribution from internal end-points (dual vertices) --
      IF (lcl_complex(dim_cmplx)%num_cobndry(i) == 1) THEN
        lcl_complex(dim_cmplx)%dual_dir(i,:) = lcl_complex(dim_cmplx)%cobndry(i)%sgn(1)*&
          (lcl_complex(dim_cmplx+1)%centers(lcl_complex(dim_cmplx)%cobndry(i)%indx(1),:) - &
          lcl_complex(dim_cmplx)%centers(i,:))
      ELSEIF (lcl_complex(dim_cmplx)%num_cobndry(i) == 2) THEN
        lcl_complex(dim_cmplx)%dual_dir(i,:) = lcl_complex(dim_cmplx)%cobndry(i)%sgn(1)*( &
          lcl_complex(dim_cmplx+1)%centers(lcl_complex(dim_cmplx)%cobndry(i)%indx(1),:) - &
          lcl_complex(dim_cmplx+1)%centers(lcl_complex(dim_cmplx)%cobndry(i)%indx(2),:) )
      END IF
    END DO

    !-- exchange external dual edge directions to adjacent processes --
    CALL exchange_dual_dir(dim_cmplx)

    !-- normalise direction vectors --
    DO i = 1,num_elm(dim_cmplx)
      lcl_complex(dim_cmplx)%dual_dir(i,:) = lcl_complex(dim_cmplx)%dual_dir(i,:) / &
        max(sqrt(dot_product(lcl_complex(dim_cmplx)%dual_dir(i,:),lcl_complex(dim_cmplx)%dual_dir(i,:))),small)
    END DO

    RETURN

  END SUBROUTINE
! geometry_mod/calc_dual_dir
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_barycentric_grad
!* SYNOPSIS
  SUBROUTINE calc_barycentric_grad(pts,n)
!* PURPOSE
!*   Compute barycentric gradients
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
!> Compute barycentric gradients
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(INOUT) :: pts(:,:) !> vertices of simplex
    INTEGER, INTENT(IN)           :: n        !> number of vertices

    !-- local variables --
    INTEGER                       :: i,j      !> loop counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute n-1 edge gradients --
    DO i =1, n-1;   pts(i,:) = pts(i,:) - pts(n,:);   END DO
    SELECT CASE (n)
    CASE(2)
      !-- do nothing: already finished --
    CASE(3)
      pts(1:n-1,:) = MATMUL(matinv2(MATMUL(pts(1:n-1,:),TRANSPOSE(pts(1:n-1,:)))),pts(1:n-1,:))
    CASE(4)
      pts(1:n-1,:) = MATMUL(matinv3(MATMUL(pts(1:n-1,:),TRANSPOSE(pts(1:n-1,:)))),pts(1:n-1,:))
    CASE DEFAULT
      !-- not setup for for other dimension simplices --
    END SELECT
    !-- compute last gradient using the fact that the sum of gradients must be zero --
    pts(n,:) = -SUM(pts(1:n-1,:),DIM=1)

    RETURN

  END SUBROUTINE
! geometry_mod/calc_barycentric_grad
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/calc_whitney_C2_BC
!* SYNOPSIS
  SUBROUTINE calc_whitney_C2_BC()
!* PURPOSE
!*   Compute the Whitney interpolation of primal faces to primal volume barycentric centers
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
!> Compute the Whitney interpolation of primal faces to primal volume barycentric centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER           :: i,j,k      !> loop indices
    INTEGER           :: vrts(dim_cmplx+1)    !> volume vertex indices
    INTEGER           :: f_vrts(dim_cmplx)  !> face vertex indices
    INTEGER           :: indx       !> face index
    LOGICAL           :: not_found  !> is boundary face index found
    REAL(KIND=PGMSiwp)    :: grads(dim_cmplx+1,dim_embbd) !> barycentric gradients

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- zero solution array --
    lcl_complex(dim_cmplx+1)% whtny_sol = 0.d0

    !-- loop through volumes computing whitney interpolation --
    DO i=1,num_elm(dim_cmplx+1)
      !-- get local vertex indices of volume --
      vrts = lcl_complex(dim_cmplx+1)% lcl_node_indx(i,:)
      !-- compute barycentric gradients of volume --
      grads = lcl_complex(1)%centers(vrts,:)
      CALL calc_barycentric_grad(grads,dim_cmplx+1)

      !-- compute contribution from each face of volume --
      DO j = 1,dim_cmplx+1
        !-- get face vertices --
        f_vrts = (/vrts(1:dim_cmplx+1-j), vrts(dim_cmplx+1-j+2:dim_cmplx+1)/)
        !-- find face local index of face --
        not_found = .TRUE.
        !-- check local faces --
        DO k = 1,lcl_complex(dim_cmplx+1)% num_bndry(i)
          indx = lcl_complex(dim_cmplx+1)% bndry(i)% indx(k)
          IF (ALL(f_vrts == lcl_complex(dim_cmplx)% lcl_node_indx(indx,:))) THEN
            !-- match found; exit search --
            not_found = .FALSE.;   EXIT
          END IF
        END DO

        !-- compute contribution from face at barycenter using Whitney interpolation --
        !-- build index array for computation
        f_vrts = (/(k,k=1,dim_cmplx+1-j), (k,k=dim_cmplx+1-j+2,dim_cmplx+1)/)
        !-- add contribution --
        IF (dim_cmplx==3) THEN
          lcl_complex(dim_cmplx+1)% whtny_sol(i,:) = lcl_complex(dim_cmplx+1)% whtny_sol(i,:) + &
              0.5d0*lcl_complex(dim_cmplx)% prml_sol(indx,1)*(&
              cross_product(grads(f_vrts(1),:),grads(f_vrts(3),:)) + &
              cross_product(grads(f_vrts(2),:),grads(f_vrts(1),:)) + &
              cross_product(grads(f_vrts(3),:),grads(f_vrts(2),:)) )
        ELSEIF (dim_cmplx==2) THEN
          lcl_complex(dim_cmplx+1)% whtny_sol(i,1) = lcl_complex(dim_cmplx+1)% whtny_sol(i,1) + &
              lcl_complex(dim_cmplx)% prml_sol(indx,1)*(grads(f_vrts(2),2) - grads(f_vrts(1),2))/3.d0
          lcl_complex(dim_cmplx+1)% whtny_sol(i,2) = lcl_complex(dim_cmplx+1)% whtny_sol(i,2) + &
              lcl_complex(dim_cmplx)% prml_sol(indx,1)*(grads(f_vrts(1),1) - grads(f_vrts(2),1))/3.d0
        END IF
      END DO
      !write(*,*) '  -',i,lcl_complex(dim_cmplx+1)% whtny_sol(i,:)
    END DO

    RETURN

  END SUBROUTINE
! geometry_mod/calc_whitney_C2_BC
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* geometry_mod/elm2smplx
!* SYNOPSIS
  SUBROUTINE elm2smplx(n_smplx,smplx,pts,element,nod)
!* PURPOSE
!*   Compute the Whitney interpolation of primal faces to primal volume barycentric centers
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
!> Compute the Whitney interpolation of primal faces to primal volume barycentric centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(OUT)              :: n_smplx, nod !>
    INTEGER, ALLOCATABLE, INTENT(OUT) :: smplx(:,:) !>
    REAL(KIND=PGMSiwp), INTENT(IN)    :: pts(:,:) !> vertices of element
    CHARACTER(LEN=*), INTENT(IN)      :: element  !> element type

    !-- local variables --
    INTEGER           :: ori
    LOGICAL           :: sqr(6), to_split     !> loop indices

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- zero solution array --
    IF (ALLOCATED(smplx)) DEALLOCATE(smplx)

    !--
    SELECT CASE (element)
    CASE('line')
      IF (nod /= 2) CALL end_mpi()
      n_smplx=1;  ALLOCATE(smplx(1,2));  smplx(1,:) = (/ 1, 2 /)
    CASE('triangle')
      IF (nod /= 3) CALL end_mpi()
      n_smplx=1;  ALLOCATE(smplx(1,3));  smplx(1,:) = (/ 1, 2, 3 /)
    CASE('tetrahedron')
      IF (nod /= 4) CALL end_mpi()
      n_smplx=1;  ALLOCATE(smplx(1,4));  smplx(1,:) = (/ 1, 2, 3, 4 /)
    CASE('quadrilateral')
      IF (nod /= 4) CALL end_mpi()
      n_smplx=2; ALLOCATE(smplx(2,3))
      IF ( dist(pts(1,:),pts(3,:)) < dist(pts(2,:),pts(4,:)) ) THEN
        smplx(1,:) = (/ 1, 2, 3 /); smplx(2,:) = (/ 1, 3, 4 /)
      ELSE;  smplx(1,:) = (/ 1, 2, 4 /); smplx(2,:) = (/ 2, 3, 4 /);  END IF
    CASE('hexahedron')
      IF (nod /= 8) CALL end_mpi()

      sqr = .FALSE.;  ori = 0
      CALL get_diag((/1,2,3,4/),pts,sqr,ori,1)      ![1 2 3 4] zeta
      CALL get_diag((/1,5,6,2/),pts,sqr,ori,2)      ![1 2 6 5] eta
      CALL get_diag((/1,4,8,5/),pts,sqr,ori,3)      ![4 1 5 8] xi
      CALL get_diag((/7,6,5,8/),pts,sqr,ori,4)      ![5 6 7 8] zeta
      CALL get_diag((/7,8,4,3/),pts,sqr,ori,5)      ![3 4 8 7] eta
      CALL get_diag((/7,3,2,6/),pts,sqr,ori,6)      ![2 3 7 6] xi

      to_split = .TRUE.
      DO WHILE(to_split)
        to_split = .FALSE.
        SELECT CASE(ori)
        CASE(0);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 3, 7  ,  1, 2, 6, 7  ,  1, 5, 6, 7  , &
          1, 3, 4, 7  ,  1, 4, 8, 7  ,  1, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(1);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 4, 7  ,  1, 4, 8, 7  ,  2, 3, 4, 7  , &
          1, 2, 6, 7  ,  1, 5, 6, 7  ,  1, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(2);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 3, 7  ,  1, 2, 5, 7  ,  2, 5, 6, 7  , &
          1, 3, 4, 7  ,  1, 4, 8, 7  ,  1, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(3);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 4, 7  ,  1, 4, 8, 7  ,  2, 3, 4, 7  , &
          1, 2, 5, 7  ,  2, 5, 6, 7  ,  1, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(4);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 3, 7  ,  1, 2, 6, 7  ,  1, 5, 6, 7  , &
          1, 3, 4, 7  ,  1, 4, 5, 7  ,  4, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(5);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 4, 5, 7  ,  1, 5, 6, 7  ,  4, 5, 8, 7  , &
          1, 2, 4, 7  ,  1, 2, 6, 7  ,  2, 3, 4, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(6);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 3, 7  ,  1, 2, 5, 7  ,  2, 5, 6, 7  , &
          1, 3, 4, 7  ,  1, 4, 5, 7  ,  4, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(7);  n_smplx=5;  ALLOCATE(smplx(5,4));  smplx = RESHAPE( (/ &
          1, 2, 4, 5  ,  2, 5, 6, 7  ,  2, 3, 4, 7  , &
          4, 5, 7, 8  ,  2, 4, 5, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(8);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 3, 4, 7  ,  1, 4, 8, 7  ,  1, 2, 3, 7  , &
          1, 2, 6, 7  ,  1, 6, 8, 7  ,  1, 5, 6, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(9);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 4, 7  ,  1, 4, 8, 7  ,  2, 3, 4, 7  , &
          1, 2, 6, 7  ,  1, 5, 6, 8  ,  1, 6, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(10);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori + 1;
          ELSEIF (sqr(2)) THEN;  ori = ori - 2;
          ELSEIF (sqr(4)) THEN;  ori = ori - 8;
          ELSEIF (sqr(5)) THEN;  ori = ori + 16;
          ELSE;  ori = ori + 100;  END IF
        CASE(11);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 1, 4, 8  ,  2, 4, 7, 8  ,  2, 3, 4, 7  , &
          2, 1, 5, 8  ,  2, 5, 6, 8  ,  2, 6, 7, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(12);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori + 1;
          ELSEIF (sqr(3)) THEN;  ori = ori - 4;
          ELSEIF (sqr(4)) THEN;  ori = ori - 8;
          ELSEIF (sqr(6)) THEN;  ori = ori + 32;
          ELSE;  ori = ori + 100; END IF
        CASE(13);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
            4, 5, 8, 6  ,  4, 1, 5, 6  ,  4, 7, 8, 6  , &
            4, 1, 2, 6  ,  4, 2, 7, 6  ,  4, 2, 3, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(14);  to_split = .TRUE.
          IF (sqr(1)) THEN;                   ori = ori + 1;
          ELSEIF (sqr(4)) THEN;               ori = ori - 8;
          ELSEIF (sqr(2) .AND. sqr(3)) THEN;  ori = ori - 6;
          ELSEIF (sqr(2) .AND. sqr(6)) THEN;  ori = ori + 30;
          ELSEIF (sqr(3) .AND. sqr(5)) THEN;  ori = ori + 12;
          ELSEIF (sqr(5) .AND. sqr(6)) THEN;  ori = ori + 48;
          ELSE;  ori = ori + 100;  END IF
        CASE(15);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 5, 8, 6  ,  4, 2, 5, 6  ,  4, 1, 2, 5  , &
          4, 7, 8, 6  ,  4, 2, 7, 6  ,  4, 2, 3, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(16);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 3, 7  ,  1, 3, 8, 7  ,  1, 3, 4, 8  , &
          1, 5, 6, 7  ,  1, 2, 6, 7  ,  1, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(17);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori - 1;
          ELSEIF (sqr(2)) THEN;  ori = ori + 2;
          ELSEIF (sqr(4)) THEN;  ori = ori + 8;
          ELSEIF (sqr(5)) THEN;  ori = ori - 16;
          ELSE;  ori = ori + 100;  END IF
        CASE(18);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 5, 7  ,  1, 2, 3, 7  ,  2, 5, 6, 7  , &
          1, 5, 8, 7  ,  1, 3, 8, 7  ,  1, 3, 4, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(19);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 3, 7, 8  ,  2, 1, 4, 8  ,  2, 3, 4, 8  , &
          2, 1, 5, 8  ,  2, 5, 7, 8  ,  2, 5, 6, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(20);  to_split = .TRUE.
          IF (sqr(2)) THEN;      ori = ori + 2;
          ELSEIF (sqr(3)) THEN;  ori = ori - 4;
          ELSEIF (sqr(5)) THEN;  ori = ori - 16;
          ELSEIF (sqr(6)) THEN;  ori = ori + 32;
          ELSE;  ori = ori + 100;  END IF
        CASE(21);  to_split = .TRUE.
          IF (sqr(2)) THEN;                   ori = ori + 2;
          ELSEIF (sqr(5)) THEN;               ori = ori - 16;
          ELSEIF (sqr(1) .AND. sqr(3)) THEN;  ori = ori - 5;
          ELSEIF (sqr(1) .AND. sqr(6)) THEN;  ori = ori + 31;
          ELSEIF (sqr(3) .AND. sqr(4)) THEN;  ori = ori + 4;
          ELSEIF (sqr(4) .AND. sqr(6)) THEN;  ori = ori + 40;
          ELSE;  ori = ori + 100; END IF
        CASE(22);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 1, 2, 5  ,  3, 2, 7, 5  ,  2, 6, 7, 5  , &
          3, 1, 4, 5  ,  3, 4, 8, 5  ,  3, 7, 8, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(23);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 2, 7, 5  ,  3, 7, 8, 5  ,  2, 6, 7, 5  , &
          3, 2, 4, 5  ,  3, 4, 8, 5  ,  1, 2, 4, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(24);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 2, 3, 7  ,  1, 3, 8, 7  ,  1, 3, 4, 8  , &
          1, 2, 6, 7  ,  1, 6, 8, 7  ,  1, 5, 6, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(25);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 3, 7, 8  ,  2, 3, 4, 8  ,  2, 1, 4, 8  , &
          2, 1, 6, 8  ,  2, 6, 7, 8  ,  1, 5, 6, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(26);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 3, 7, 8  ,  2, 1, 3, 8  ,  1, 3, 4, 8  , &
          2, 5, 6, 8  ,  2, 6, 7, 8  ,  2, 1, 5, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(27);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 3, 7, 8  ,  2, 3, 4, 8  ,  2, 1, 4, 8  , &
          2, 5, 6, 8  ,  2, 6, 7, 8  ,  2, 1, 5, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(28);  to_split = .TRUE.
          IF (sqr(3)) THEN;                   ori = ori - 4;
          ELSEIF (sqr(6)) THEN;               ori = ori + 32;
          ELSEIF (sqr(1) .AND. sqr(2)) THEN;  ori = ori + 3;
          ELSEIF (sqr(1) .AND. sqr(4)) THEN;  ori = ori - 15;
          ELSEIF (sqr(2) .AND. sqr(4)) THEN;  ori = ori - 6;
          ELSEIF (sqr(4) .AND. sqr(5)) THEN;  ori = ori - 24;
          ELSE;  ori = ori + 100;  END IF
        CASE(29);  to_split = .TRUE.
          IF (sqr(2)) THEN;      ori = ori + 2;
          ELSEIF (sqr(3)) THEN;  ori = ori - 4;
          ELSEIF (sqr(5)) THEN;  ori = ori - 16;
          ELSEIF (sqr(6)) THEN;  ori = ori + 32;
          ELSE;  ori = ori + 100;  END IF
        CASE(30);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori + 1;
          ELSEIF (sqr(3)) THEN;  ori = ori - 4;
          ELSEIF (sqr(4)) THEN;  ori = ori - 8;
          ELSEIF (sqr(6)) THEN;  ori = ori + 32;
          ELSE;  ori = ori + 100; END IF
        CASE(31);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 3, 7, 8  ,  2, 3, 4, 8  ,  2, 6, 7, 8  , &
          2, 5, 6, 8  ,  2, 4, 5, 8  ,  2, 1, 4, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(32);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 5, 6, 7  ,  1, 3, 6, 7  ,  1, 2, 3, 6  , &
          1, 3, 4, 7  ,  1, 4, 8, 7  ,  1, 5, 8, 7 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(33);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori - 1;
          ELSEIF (sqr(3)) THEN;  ori = ori + 4;
          ELSEIF (sqr(4)) THEN;  ori = ori + 8;
          ELSEIF (sqr(6)) THEN;  ori = ori - 32;
          ELSE;  ori = ori + 100; END IF
        CASE(34);  to_split = .TRUE.
          IF (sqr(2)) THEN;      ori = ori - 2;
          ELSEIF (sqr(3)) THEN;  ori = ori + 4;
          ELSEIF (sqr(5)) THEN;  ori = ori + 16;
          ELSEIF (sqr(6)) THEN;  ori = ori - 32;
          ELSE;  ori = ori + 100;  END IF
        CASE(35);  to_split = .TRUE.
          IF (sqr(3)) THEN;                   ori = ori + 4;
          ELSEIF (sqr(6)) THEN;               ori = ori - 32;
          ELSEIF (sqr(1) .AND. sqr(2)) THEN;  ori = ori - 3;
          ELSEIF (sqr(1) .AND. sqr(5)) THEN;  ori = ori + 15;
          ELSEIF (sqr(2) .AND. sqr(4)) THEN;  ori = ori + 6;
          ELSEIF (sqr(4) .AND. sqr(5)) THEN;  ori = ori + 24;
          ELSE;  ori = ori + 100;  END IF
        CASE(36);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 6, 7, 5  ,  3, 4, 7, 5  ,  4, 7, 8, 5  , &
          3, 1, 4, 5  ,  3, 1, 6, 5  ,  3, 1, 2, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(37);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 3, 7, 6  ,  4, 5, 7, 6  ,  4, 5, 7, 8  , &
          4, 1, 2, 6  ,  4, 2, 3, 6  ,  4, 1, 5, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(38);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 6, 7, 5  ,  3, 4, 7, 5  ,  4, 7, 8, 5  , &
          3, 1, 2, 5  ,  3, 2, 6, 5  ,  3, 1, 4, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(39);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 3, 7, 6  ,  4, 5, 7, 6  ,  4, 5, 7, 8  , &
          4, 2, 3, 6  ,  4, 2, 5, 6  ,  4, 1, 2, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(40);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 1, 8, 6  ,  4, 7, 8, 6  ,  1, 5, 8, 6  , &
          4, 1, 3, 6  ,  4, 3, 7, 6  ,  1, 2, 3, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(41);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 1, 8, 6  ,  4, 7, 8, 6  ,  1, 5, 8, 6  , &
          4, 3, 7, 6  ,  4, 2, 3, 6  ,  4, 1, 2, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(42);  to_split = .TRUE.
          IF (sqr(2)) THEN;                   ori = ori - 2;
          ELSEIF (sqr(5)) THEN;               ori = ori + 16;
          ELSEIF (sqr(1) .AND. sqr(3)) THEN;  ori = ori + 5;
          ELSEIF (sqr(1) .AND. sqr(6)) THEN;  ori = ori - 31;
          ELSEIF (sqr(3) .AND. sqr(4)) THEN;  ori = ori - 4;
          ELSEIF (sqr(4) .AND. sqr(6)) THEN;  ori = ori - 40;
          ELSE;  ori = ori + 100;  END IF
        CASE(43);  to_split = .TRUE.
          IF (sqr(2)) THEN;      ori = ori - 2;
          ELSEIF (sqr(3)) THEN;  ori = ori + 4;
          ELSEIF (sqr(5)) THEN;  ori = ori + 16;
          ELSEIF (sqr(6)) THEN;  ori = ori - 32;
          ELSE;  ori = ori + 100;  END IF
        CASE(44);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 1, 5, 6  ,  4, 5, 8, 6  ,  4, 7, 8, 6  , &
          4, 1, 3, 6  ,  4, 3, 7, 6  ,  1, 2, 3, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(45);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 3, 7, 6  ,  4, 5, 8, 6  ,  4, 7, 8, 6  , &
          4, 1, 2, 6  ,  4, 2, 3, 6  ,  4, 1, 5, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(46);  to_split = .TRUE.
          IF (sqr(2)) THEN;      ori = ori + 2;
          ELSEIF (sqr(3)) THEN;  ori = ori - 4;
          ELSEIF (sqr(5)) THEN;  ori = ori - 16;
          ELSEIF (sqr(6)) THEN;  ori = ori + 32;
          ELSE;  ori = ori + 100;  END IF
        CASE(47);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 3, 7, 6  ,  4, 5, 8, 6  ,  4, 7, 8, 6  , &
          4, 2, 3, 6  ,  4, 2, 5, 6  ,  4, 1, 2, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(48);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          1, 5, 6, 7  ,  1, 3, 6, 7  ,  1, 2, 3, 6  , &
          1, 5, 8, 7  ,  1, 3, 8, 7  ,  1, 3, 4, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(49);  to_split = .TRUE.
          IF (sqr(1)) THEN;                   ori = ori - 1;
          ELSEIF (sqr(4)) THEN;               ori = ori + 8;
          ELSEIF (sqr(2) .AND. sqr(3)) THEN;  ori = ori + 6;
          ELSEIF (sqr(2) .AND. sqr(6)) THEN;  ori = ori - 30;
          ELSEIF (sqr(3) .AND. sqr(5)) THEN;  ori = ori - 12;
          ELSEIF (sqr(5) .AND. sqr(6)) THEN;  ori = ori - 48;
          ELSE;  ori = ori + 100;  END IF
        CASE(50);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 6, 7, 5  ,  3, 2, 6, 5  ,  3, 1, 2, 5  , &
          3, 7, 8, 5  ,  3, 1, 8, 5  ,  3, 1, 4, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(51);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori - 1;
          ELSEIF (sqr(3)) THEN;  ori = ori + 4;
          ELSEIF (sqr(4)) THEN;  ori = ori + 8;
          ELSEIF (sqr(6)) THEN;  ori = ori - 16;
          ELSE;  ori = ori + 100;  END IF
        CASE(52);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 1, 6, 5  ,  3, 6, 7, 5  ,  3, 1, 2, 6  , &
          3, 1, 4, 5  ,  3, 4, 8, 5  ,  3, 7, 8, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(53);  to_split = .TRUE.
          IF (sqr(1)) THEN;      ori = ori - 1;
          ELSEIF (sqr(2)) THEN;  ori = ori + 2;
          ELSEIF (sqr(4)) THEN;  ori = ori + 8;
          ELSEIF (sqr(5)) THEN;  ori = ori - 32;
          ELSE;  ori = ori + 100;  END IF
        CASE(54);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 1, 2, 5  ,  3, 2, 6, 5  ,  3, 6, 7, 5  , &
          3, 1, 4, 5  ,  3, 4, 8, 5  ,  3, 7, 8, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(55);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 2, 6, 5  ,  3, 2, 4, 5  ,  1, 2, 4, 5  , &
          3, 4, 8, 5  ,  3, 7, 8, 5  ,  3, 6, 7, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(56);  n_smplx=5;  ALLOCATE(smplx(5,4));  smplx = RESHAPE( (/ &
          1, 5, 6, 8  ,  1, 2, 3, 6  ,  1, 3, 4, 8  , &
          3, 6, 7, 8  ,  1, 3, 6, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(57);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 3, 8, 6  ,  4, 2, 3, 6  ,  3, 7, 8, 6  , &
          4, 1, 2, 6  ,  4, 1, 8, 6  ,  1, 5, 8, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(58);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 2, 6, 5  ,  3, 6, 8, 5  ,  3, 6, 7, 8  , &
          3, 1, 8, 5  ,  3, 1, 2, 5  ,  3, 1, 4, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(59);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          2, 1, 5, 8  ,  2, 5, 6, 8  ,  2, 1, 4, 8  , &
          2, 3, 4, 8  ,  2, 3, 6, 8  ,  3, 6, 7, 8 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(60);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 6, 8, 5  ,  3, 4, 8, 5  ,  3, 6, 7, 8  , &
          3, 1, 4, 5  ,  3, 1, 6, 5  ,  3, 1, 2, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(61);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          4, 3, 8, 6  ,  4, 2, 3, 6  ,  3, 7, 8, 6  , &
          4, 1, 2, 6  ,  4, 1, 5, 6  ,  4, 5, 8, 6 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(62);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 6, 8, 5  ,  3, 4, 8, 5  ,  3, 6, 7, 8  , &
          3, 1, 4, 5  ,  3, 1, 2, 5  ,  3, 2, 6, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE(63);  n_smplx=6;  ALLOCATE(smplx(6,4));  smplx = RESHAPE( (/ &
          3, 6, 8, 5  ,  3, 4, 8, 5  ,  3, 6, 7, 8  , &
          3, 2, 6, 5  ,  3, 2, 4, 5  ,  1, 2, 4, 5 /) , SHAPE(smplx), ORDER=(/2,1/) )
        CASE DEFAULT
          CALL end_mpi();
        END SELECT
      END DO
    CASE DEFAULT
      CALL end_mpi();
    END SELECT

    RETURN

  END SUBROUTINE
! geometry_mod/elm2smplx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

SUBROUTINE get_diag(vrts,pts,sqr,ori,ptr)
  LOGICAL, INTENT(INOUT) :: sqr(:)
  INTEGER, INTENT(INOUT) :: ori
  INTEGER, INTENT(IN)    :: vrts(4), ptr
  REAL(KIND=PGMSiwp), INTENT(IN) :: pts(:,:)
  REAL(KIND=PGMSiwp) :: d1, d2

  d1 = dist(pts(vrts(1),:),pts(vrts(3),:))
  d2 = dist(pts(vrts(2),:),pts(vrts(4),:))

  IF ( ABS(d1 - d2) < small ) THEN;  sqr(ptr) = .TRUE.
  ELSE;  sqr(ptr) = .FALSE.
    IF ( d1 > d2 ) ori = ori + 2**(ptr-1)
  END IF

  RETURN
END SUBROUTINE

FUNCTION dist(pts1,pts2)
  REAL(KIND=PGMSiwp) :: pts1(:), pts2(:), dist
  dist = sqrt(sum((pts2-pts1)**2))
END FUNCTION

END MODULE
! geometry_mod
!===============================================================================
!
