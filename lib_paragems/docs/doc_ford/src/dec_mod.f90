!
!===============================================================================
!-- DEC Module
!> Module containing routines for performing DEC operations
!===============================================================================
!/****/h* modules|dec/dec_mod
!* SYNOPSIS
MODULE dec_mod
!* PURPOSE
!*   Module containing routines for performing DEC operations
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 ???
!*   math_mod                basic math functions
!* CONTAINS
!*   Subroutine              Purpose
!*   calc_orientation        sort nodal indices and compute ±1 orientation
!*   calc_bndry_cobndry      recursively compute element (co-)boundaries
!*   build_bndry_work_array  build boundary data working array
!*   count_bndry_cobndry     count [co-]boundaries: internal, external, surface
!*   set_bndry_cobndry       set [co-]boundaries: internal, external, surface
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
!*   calc_hodge_star         compute hodge star and it's inverse from primal and dual volumes
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

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_orientation
!* SYNOPSIS
  SUBROUTINE calc_orientation(elm,orientation)
!* PURPOSE
!*   Sort nodal indices and compute ± orientation of an element
!* INPUTS
!*   Name                    Description
!*   elem
!*   orientation
!* OUTPUTS
!*   Name                    Description
!*   elem
!*   orientation
!* SIDE EFFECTS
!*   - the nodal indices of elm are sorted numerically in ascending order
!*   - orientation contains the ±1 orientation of the element
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Procedure adapted from:
!*   https://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Sort nodal indices and compute ± orientation of an element
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(INOUT) :: elm(:)           !> nodal indices of elements
    INTEGER, INTENT(INOUT) :: orientation      !> orientations of the elements

    !-- local variables --
    INTEGER                :: i,j              !> loop counters
    INTEGER                :: length           !> size of elm array
    INTEGER                :: work             !> work variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- check size of argument --
    length  = SIZE(elm,DIM=1)

    !-------------------------------------------------------------------------
    !-- sort indices using insertion sort and count the number of swaps --
    !-------------------------------------------------------------------------
    orientation = 0
    !-- loop through all indices, starting from the second --
    DO i = 2, length
      !-- get next value to insert --
      work = elm(i)
      j = i - 1
      !-- find location to insert given value --
      DO WHILE (j >= 1)
        !-- ??? found location ??? --
        IF (elm(j) <= work) EXIT
        !-- shuffle (swap) values --
        elm(j + 1) = elm(j)
        j = j - 1
        !-- count swaps --
        orientation = orientation + 1
      END DO
      !-- insert value --
      elm(j + 1) = work
    END DO

    !-- compute orientation from number of swaps --
    orientation = 1 - 2*MOD(orientation,2)

    RETURN

  END SUBROUTINE
! dec_mod/calc_orientation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE calc_bndry_cobndry()
!* PURPOSE
!*   Recursively compute element (co-)boundaries from highest to lowest geometric order
!* ASSUMPTION
!*   Mesh is a simplicial complex in TetGen format
!* SIDE EFFECTS
!*   - the nodal indices of elm are sorted numerically in ascending order
!*   - orientation contains the ±1 orientation of the element
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Recursively compute element (co-)boundaries from highest to lowest geometric order
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- local variables --
    INTEGER               :: k              !> simplicial order

    INTEGER, ALLOCATABLE  :: bndry(:,:)     !> temp. boundary work array
    INTEGER, ALLOCATABLE  :: bndry_cnt(:)   !> counter local boundaries
    INTEGER, ALLOCATABLE  :: cobndry_cnt(:) !> counter co-boundaries

    INTEGER               :: cnt            !> counter boundaries
    INTEGER               :: ext_cnt(2)     !> counter external boundaries
    INTEGER               :: surf_cnt       !> counter surface boundaries

    INTEGER, ALLOCATABLE  :: work(:,:)      !> work array for merge sort
    LOGICAL, ALLOCATABLE  :: lwork(:)       !> work array to determine locality

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- calculate signed adjacency [co-boundary operator] --
    !-- loop for each geometric order in complex in descending order --
    DO k=dim_cmplx+1,2,-1
      CALL syncwrite_log('===============================================================')
      !-- initial guess for the number of boundaries --
      num_elm(k-1) = k*num_elm(k)

      !-- build working array with element/boundary data --
      CALL syncwrite_log('> parallel_setup() - calc_bndry_cobndry '//&
        '- build_bndry_work_array')
      ALLOCATE(bndry(num_elm(k-1),k+2))
      CALL build_bndry_work_array(bndry,k);   CALL syncwrite_log_time()

      !-- sort working array --
      CALL syncwrite_log('> parallel_setup() - calc_bndry_cobndry '//&
        '- int_merge_sort_rows')
      ALLOCATE(work(num_elm(k-1),k+2))
      CALL int_merge_sort_rows(bndry,num_elm(k-1),1,k-1,work)
      DEALLOCATE(work);   CALL syncwrite_log_time()

      !-- count number of unique (local/external/surface) [co-]boundaries --
      CALL syncwrite_log('> parallel_setup() - calc_bndry_cobndry '//&
        '- count_bndry_cobndry')
      ALLOCATE(bndry_cnt(num_elm(k)),cobndry_cnt(num_elm(k-1)),&
        lwork(num_elm(k)*k))
      CALL count_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,&
        surf_cnt,lwork);   CALL syncwrite_log_time()

      !-- allocate [co-]boundary structures and variables --
      CALL syncwrite_log('> parallel_setup() - calc_bndry_cobndry '//&
        '- allocate_bndry_cobndry')
      CALL allocate_bndry_cobndry(k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,surf_cnt)
      CALL syncwrite_log_time()

      !-- setup boundaries external to the current process (not local),
      !   setup parallel mapping, setup node indices for (k-1)^th order element
      !   structure, and setup co-boundaries for (k-1)^th order element
      !   structure --
      CALL syncwrite_log('> parallel_setup() - calc_bndry_cobndry '//&
        '- set_bndry_cobndry')
      CALL set_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,&
        surf_cnt,lwork);   CALL syncwrite_log_time()

      !-- clean up temporary variables --
      DEALLOCATE(bndry,bndry_cnt,cobndry_cnt,lwork)
    END DO
    CALL syncwrite_log('===============================================================')

    RETURN

  END SUBROUTINE
! dec_mod/calc_bndry_cobndry
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/build_bndry_work_array
!* SYNOPSIS
  SUBROUTINE build_bndry_work_array(bndry,k)
!* PURPOSE
!*   Build boundary data working array
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
!> Build boundary data working array
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(INOUT)  :: bndry(:,:)     !> boundary work array
    INTEGER, INTENT(IN)     :: k              !> simplicial order

    !-- local variables --
    INTEGER                 :: i,j            !> loop index
    INTEGER                 :: iib,fib        !> array initial/final index

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    iib = 1;   fib = num_elm(k)
    DO i=1,k
      !-- get boundary node indices --
      bndry(iib:fib,1:i-1) = lcl_complex(k)%node_indx(:,1:i-1)
      bndry(iib:fib,i:k-1) = lcl_complex(k)%node_indx(:,i+1:k)

      !-- get simplex indices --
      bndry(iib:fib,k)   = (/ (j, j=1,num_elm(k)) /)
      bndry(iib:fib,k+1) = i

      !-- get orientation --
      IF (k == dim_cmplx+1) THEN
        bndry(iib:fib,k+2)  = (-1)**i*lcl_complex(k)%orientation
      ELSE
        bndry(iib:fib,k+2)  = (-1)**i
      END IF

      !-- update array indices --
      iib = iib + num_elm(k)
      fib = fib + num_elm(k)
    END DO

    RETURN

  END SUBROUTINE
! dec_mod/build_bndry_work_array
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/count_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE count_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,&
    surf_cnt,local)
!* PURPOSE
!*   Count [co-]boundaries: internal, external, surface
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
!> Count [co-]boundaries: internal, external, surface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)     :: k                !> simplicial order
    INTEGER, INTENT(IN)     :: bndry(:,:)       !> boundary work array
    INTEGER, INTENT(INOUT)  :: bndry_cnt(:)     !> counter local boundaries
    INTEGER, INTENT(INOUT)  :: cobndry_cnt(:)   !> counter co-boundaries
    INTEGER, INTENT(INOUT)  :: cnt              !> counter (k-1) element
    INTEGER, INTENT(INOUT)  :: ext_cnt(2)       !> counter external boundaries
    INTEGER, INTENT(INOUT)  :: surf_cnt         !> counter surface boundaries
    LOGICAL, INTENT(INOUT)  :: local(:)         !> locality work array

    !-- local variables --
    INTEGER                 :: i                !> counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- reset counters --
    cnt = 0; bndry_cnt = 0; cobndry_cnt = 0; ext_cnt = 0; surf_cnt = 0


    !-- loop for all boundaries --
    DO i = 1, num_elm(k)*k
      IF (i == 1 .OR. ANY(bndry(i,1:k-1) /= bndry(max(1,i-1),1:k-1))) THEN
        !-- boundary different from previous: new boundary --
        cnt = cnt + 1
        IF (any_element_in_list(bndry(i,1:k-1),lcl_complex(1)%glb_indx)) THEN
          !-- boundary is local --
          local(i) = .TRUE.
          IF (ANY(bndry(i,1:k-1) /= bndry(min(i+1,num_elm(k)*k),1:k-1)) .OR. &
            i == num_elm(k)*k) THEN
            !-- boundary is unique: count surface boundaries --
            surf_cnt = surf_cnt + 1
          ELSE
            !-- boundary is NOT unique: count local boundaries --
            bndry_cnt(bndry(i,k)) = bndry_cnt(bndry(i,k)) + 1
          END IF
          !-- count number of co-boundaries for each local boundary --
          cobndry_cnt(cnt) = cobndry_cnt(cnt) + 1
        ELSE
          local(i) = .FALSE.
          !-- boundary is non-local: count external boundaries --
          ext_cnt = ext_cnt + 1
        END IF
      ELSE
        !-- boundary same as previous --
        IF (any_element_in_list(bndry(i,1:k-1),lcl_complex(1)%glb_indx)) THEN
          !-- boundary is local --
          local(i) = .TRUE.
          !-- count local boundaries --
          bndry_cnt(bndry(i,k)) = bndry_cnt(bndry(i,k)) + 1
          !-- count number of co-boundaries for each local boundary --
          cobndry_cnt(cnt) = cobndry_cnt(cnt) + 1
        ELSE
          local(i) = .FALSE.
          !-- boundary is non-local: count external boundaries --
          ext_cnt(1) = ext_cnt(1) + 1
        END IF
      END IF
    END DO

    RETURN

  END SUBROUTINE
! dec_mod/count_bndry_cobndry
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/allocate_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE allocate_bndry_cobndry(k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,&
    surf_cnt)
!* PURPOSE
!*   Allocate [co-]boundaries structures/variables
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
!> Allocate [co-]boundaries structures/variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)     :: k                !> simplicial order
    INTEGER, INTENT(IN)     :: bndry_cnt(:)     !> temp # local boundaries
    INTEGER, INTENT(IN)     :: cobndry_cnt(:)   !> temp # co-boundaries
    INTEGER, INTENT(IN)     :: cnt              !> counter (k-1) element
    INTEGER, INTENT(IN)     :: ext_cnt(2)       !> counter external boundaries
    INTEGER, INTENT(IN)     :: surf_cnt         !> counter surface boundaries

    !-- local variables --
    INTEGER                 :: i                !> counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate boundary structures and variables: k^th order --
    !-- # external boundaries and # expect MPI receives --
    lcl_complex(k)%num_ext = ext_cnt(1)
    lcl_complex(k)%num_recv = ext_cnt(2)
    ALLOCATE(&
      lcl_complex(k)%ext_indx(ext_cnt(1),3),&
      lcl_complex(k)%recv_indx(ext_cnt(2),2))
    !-- # surface boundaries --
    lcl_complex(k)%num_surf = surf_cnt
    ALLOCATE(lcl_complex(k)%surf_indx(surf_cnt,4))
    !-- # internal boundaries --
    ALLOCATE(&
      lcl_complex(k)%num_bndry(num_elm(k)),&
      lcl_complex(k)%bndry(num_elm(k)))
    lcl_complex(k)%num_bndry = bndry_cnt
    DO i=1,num_elm(k)
      ALLOCATE(&
        lcl_complex(k)%bndry(i)%sgn(lcl_complex(k)%num_bndry(i)),&
        lcl_complex(k)%bndry(i)%indx(lcl_complex(k)%num_bndry(i)))
    END DO

    !-- allocate [co-]boundary structures and variables: (k-1))^th order --
    !-- # (k-1))^th order elements --
    num_elm(k-1) = cnt
    ALLOCATE(lcl_complex(k-1)%node_indx(num_elm(k-1),k-1))
    !-- # (k-1))^th order coboundaries --
    ALLOCATE(&
      lcl_complex(k-1)%num_cobndry(num_elm(k-1)),&
      lcl_complex(k-1)%cobndry(num_elm(k-1)))
    lcl_complex(k-1)%num_cobndry = cobndry_cnt(1:num_elm(k-1))
    DO i=1,num_elm(k-1)
      ALLOCATE(&
        lcl_complex(k-1)%cobndry(i)%sgn(cobndry_cnt(i)),&
        lcl_complex(k-1)%cobndry(i)%indx(cobndry_cnt(i)))
    END DO

    RETURN

  END SUBROUTINE
! dec_mod/allocate_bndry_cobndry
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/set_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE set_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,surf_cnt,local)
!* PURPOSE
!*   Set [co-]boundaries: internal, external, surface
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
!> Set [co-]boundaries: internal, external, surface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)     :: k                !> simplicial order
    INTEGER, INTENT(INOUT)  :: bndry(:,:)       !> boundary work array
    INTEGER, INTENT(INOUT)  :: bndry_cnt(:)     !> temp # local boundaries
    INTEGER, INTENT(INOUT)  :: cobndry_cnt(:)   !> temp # co-boundaries
    INTEGER, INTENT(INOUT)  :: cnt              !> counter (k-1) element
    INTEGER, INTENT(INOUT)  :: ext_cnt(2)       !> counter external boundaries
    INTEGER, INTENT(INOUT)  :: surf_cnt         !> counter surface boundaries
    LOGICAL, INTENT(INOUT)  :: local(:)         !> locality work array

    !-- local variables --
    INTEGER                 :: i                !> counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- reset counters --
    cnt = 0;  bndry_cnt = 0;  cobndry_cnt = 0;  ext_cnt = 0;  surf_cnt = 0

    !-- loop for all boundaries --
    DO i = 1, num_elm(k)*k
      IF (i == 1 .OR. ANY(bndry(i,1:k-1) /= bndry(max(1,i-1),1:k-1))) THEN
        !-- boundary different from previous: new boundary --
        cnt = cnt + 1
        IF (local(i)) THEN
          !-- boundary is local --
          IF (ANY(bndry(i,1:k-1) /= bndry(min(i+1,num_elm(k)*k),1:k-1)) .OR. &
            i == num_elm(k)*k) THEN
            !-- boundary is unique: set signed adjacency for surface boundaries --
            surf_cnt = surf_cnt + 1
            lcl_complex(k)%surf_indx(surf_cnt,:) = (/ -bndry(i,k+2), cnt, bndry(i,k) /)
          ELSE
            !-- boundary is NOT unique: set signed adjacency for local boundaries --
            bndry_cnt(bndry(i,k)) = bndry_cnt(bndry(i,k)) + 1
            lcl_complex(k)%bndry(bndry(i,k))%sgn(bndry_cnt(bndry(i,k))) = bndry(i,k+2)
            lcl_complex(k)%bndry(bndry(i,k))%indx(bndry_cnt(bndry(i,k))) = cnt
          END IF
          !-- set signed adjacency for co-boundaries for each local boundary --
          cobndry_cnt(cnt) = 1
          lcl_complex(k-1)%cobndry(cnt)%sgn(1) = bndry(i,k+2)
          lcl_complex(k-1)%cobndry(cnt)%indx(1) = bndry(i,k)
        ELSE
          !-- boundary is non-local: set signed adjacency for external boundaries --
          ext_cnt = ext_cnt + 1
          lcl_complex(k)%ext_indx(ext_cnt(1),:) = (/ bndry(i,k+2), cnt, bndry(i,k) /)
          lcl_complex(k)%recv_indx(ext_cnt(2),:) = (/ elm2proc(bndry(i,1)), cnt /)                          !????? find actual global index
        END IF
        !-- set node indices for (k-1)^th order elements --
        lcl_complex(k-1)%node_indx(cnt,:) = bndry(i,1:k-1)
      ELSE
        !-- boundary same as previous --
        IF (local(i)) THEN
          !-- boundary is local: set signed adjacency for local boundaries --
          bndry_cnt(bndry(i,k)) = bndry_cnt(bndry(i,k)) + 1
          lcl_complex(k)%bndry(bndry(i,k))%sgn(bndry_cnt(bndry(i,k))) = bndry(i,k+2)
          lcl_complex(k)%bndry(bndry(i,k))%indx(bndry_cnt(bndry(i,k))) = cnt
          !-- set signed adjacency for co-boundaries for each local boundary --
          cobndry_cnt(cnt) = cobndry_cnt(cnt) + 1
          lcl_complex(k-1)%cobndry(cnt)%sgn(cobndry_cnt(cnt)) = bndry(i,k+2)
          lcl_complex(k-1)%cobndry(cnt)%indx(cobndry_cnt(cnt)) = bndry(i,k)
        ELSE
          !-- boundary is non-local: set signed adjacency for external boundaries --
          ext_cnt(1) = ext_cnt(1) + 1
          lcl_complex(k)%ext_indx(ext_cnt(1),:) = (/ bndry(i,k+2), cnt, bndry(i,k) /)
        END IF
      END IF
    END DO

    RETURN

  END SUBROUTINE
! dec_mod/set_bndry_cobndry
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/initialise_geo
!* SYNOPSIS
  SUBROUTINE initialise_geo()
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

    !-- local variables --
    INTEGER                 :: i,k               !> loop counters
    REAL(KIND=iwp)          :: minV,maxV,absminV !> min/max volumes
    CHARACTER(LEN=slen)     :: str               !> string for writing to log file

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL syncwrite_log('> initialise_geo()')
    curnt_time = MPI_WTIME()

    !-- log element information --
    DO k =1,dim_cmplx+1
      write(str,'(A,I1,A,I10)') '   - number of C',k-1,': ',glb_num_elm(k)
      CALL rootwrite_log(str)
    END DO

    !-- get local node indices --
    CALL get_lcl_node_indx()

    !-- compute circumcenters --
    DO k =2,dim_cmplx+1
      CALL calc_circumcenters(k)
    END DO

    !-- compute signed primal volumes --
    IF (.NOT. ALLOCATED(lcl_complex(1)%prml_volume)) &
      ALLOCATE(lcl_complex(1)%prml_volume(num_elm(1)))
    lcl_complex(1)%prml_volume = 1
    DO k =2,dim_cmplx+1
      IF (k-1==dim_embbd) THEN
        CALL calc_prml_sgnd_vlm(k)
      ELSE
        CALL calc_prml_unsgnd_vlm(k)
      END IF
    END DO

    !-- compute signed dual volumes --
    CALL calc_dual_vlm()
    DO k =1,dim_cmplx
      CALL exchange_dual_vlm(k)
    END DO

    !-- assess maximumm and minimum volumes --
    DO k =1,dim_cmplx+1
      CALL MPI_REDUCE(minval(lcl_complex(k)%prml_volume),minV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(minval(abs(lcl_complex(k)%prml_volume)),absminV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(maxval(lcl_complex(k)%prml_volume),maxV,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)
      IF (rank == root) &
        write(str,'(A,I1,A,ES10.3,A,ES10.3,A,ES10.3)') '   - The min/abs(min)/max C',k-1,&
          ' volume is: ',minV,'/',absminV,'/',maxV
      CALL syncwrite_log(str)

      CALL MPI_REDUCE(minval(lcl_complex(k)%dual_volume),minV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(minval(abs(lcl_complex(k)%dual_volume)),absminV,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ier)
      CALL MPI_REDUCE(maxval(lcl_complex(k)%dual_volume),maxV,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)
      IF (rank == root) &
        write(str,'(A,I1,A,ES10.3,A,ES10.3,A,ES10.3)') '   - The min/abs(min)/max D',dim_cmplx+1-k,&
          ' volume is: ',minV,'/',absminV,'/',maxV
      CALL syncwrite_log(str)
    END DO

    !-- compute Hodge stars --
    DO k = 1,dim_cmplx+1
      CALL calc_hodge_star(k)
    END DO

    !-- compute primal and dual edge directions --
    CALL calc_prml_dir()
    CALL calc_dual_dir()

    !-- log time to initialise the geometric information --
    CALL syncwrite_log_time()

    RETURN

  END SUBROUTINE
! dec_mod/initialise_geo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/get_lcl_node_indx
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
! dec_mod/get_lcl_node_indx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_circumcenters
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
    REAL(KIND=iwp), ALLOCATABLE :: A(:,:),b(:)      !> solution variables
    REAL(KIND=iwp), ALLOCATABLE :: pts(:,:)         !> bounding points
    REAL(KIND=iwp), ALLOCATABLE :: ipiv(:),work(:)  !> work arrays

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
! dec_mod/calc_circumcenters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_prml_sgnd_vlm
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
    REAL(KIND=iwp), ALLOCATABLE :: A(:,:)       !> solution variable
    REAL(KIND=iwp), ALLOCATABLE :: pts(:)       !> bounding points

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate volumes and temp. variables --
    IF (.NOT. ALLOCATED(lcl_complex(k)%prml_volume)) &
      ALLOCATE(lcl_complex(k)%prml_volume(num_elm(k)))
    ALLOCATE(pts(dim_embbd),A(dim_embbd,dim_embbd))

    !-- compute factorial of embedding dimension --
    fac = 1
    DO i=2,dim_embbd
      fac = fac*i
    END DO

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
      lcl_complex(k)%prml_volume(i) = lcl_complex(k)%orientation(i)*determinant(A,dim_embbd)/fac
    END DO

    !-- clean up --
    DEALLOCATE(pts,A)

    RETURN

  END SUBROUTINE
! dec_mod/calc_prml_sgnd_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_prml_unsgnd_vlm
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
    REAL(KIND=iwp), ALLOCATABLE :: A(:,:)       !> solution variable
    REAL(KIND=iwp), ALLOCATABLE :: pts(:)       !> bounding points

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate volumes and temp. variables --
    IF (.NOT. ALLOCATED(lcl_complex(k)%prml_volume)) &
      ALLOCATE(lcl_complex(k)%prml_volume(num_elm(k)))
    ALLOCATE(pts(dim_embbd),A(k-1,dim_embbd))

    !-- compute necessary factorial --
    fac = 1
    DO i=2,k-1
      fac = fac*i
    END DO

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
      lcl_complex(k)%prml_volume(i) = SQRT(ABS(determinant(MATMUL(A,TRANSPOSE(A)),k-1)))/fac
    END DO

    !-- clean up --
    DEALLOCATE(pts,A)

    RETURN

  END SUBROUTINE
! dec_mod/calc_prml_unsgnd_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_dual_vlm
!* SYNOPSIS
  SUBROUTINE calc_dual_vlm()
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

    !-- local variables --
    INTEGER                 :: k                !> simplicial order
    INTEGER                 :: i,indx           !> loop counters
    REAL(KIND=iwp), ALLOCATABLE :: pts(:,:)     !> bounding points
    REAL(KIND=iwp), ALLOCATABLE :: sgn(:)       !> sign of volume elements

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate volumes and temp. variables --
    DO k = 1,dim_cmplx+1
      IF (.NOT. ALLOCATED(lcl_complex(k)%dual_volume)) &
        ALLOCATE(lcl_complex(k)%dual_volume(num_elm(k)))
      lcl_complex(k)%dual_volume = 0.d0
    END DO
    ALLOCATE(pts(dim_cmplx+1,dim_embbd),sgn(dim_cmplx))
    sgn = 1.d0

    !-- start recursive computation of dual volumes on the interior --
    k = dim_cmplx+1
    DO i = 1, num_elm(k)
      CALL calc_dual_vlm_i(pts,sgn,i,i,k)
    END DO

    !-- start recurve addition of contributions from the surfaces --
    k = dim_cmplx+1
    DO i = 1, lcl_complex(k)%num_surf
      indx = lcl_complex(k)% surf_indx(i,2)
      CALL add_points_i(pts,indx,k-1)
      CALL calc_dual_vlm_i(pts,sgn,indx,lcl_complex(k)% surf_indx(i,3),k-1)
    END DO

    !-- clean up --
    DEALLOCATE(pts,sgn)

    RETURN

  END SUBROUTINE
! dec_mod/calc_dual_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_dual_vlm_i
!* SYNOPSIS
  RECURSIVE SUBROUTINE calc_dual_vlm_i(pts,sgn,indx,p_indx,k)
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
    INTEGER, INTENT(IN)           :: k            !> simplicial order
    INTEGER, INTENT(IN)           :: indx,p_indx  !> element and parent index
    REAL(KIND=iwp), INTENT(INOUT) :: pts(:,:)     !> bounding points
    REAL(KIND=iwp), INTENT(INOUT) :: sgn(:)       !> sign of volume elements

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
    IF (k>1 .AND. k<=dim_cmplx) THEN
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
    IF (k>1) THEN
      DO i = 1, lcl_complex(k)% num_bndry(indx)
        CALL calc_dual_vlm_i(pts,sgn,lcl_complex(k)% bndry(indx)% indx(i),indx,k-1)
      END DO
    END IF

    IF (k>1 .AND. k<=dim_cmplx) sgn(k-1) = 1.d0

    RETURN

  END SUBROUTINE
! dec_mod/calc_dual_vlm_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/add_points_i
!* SYNOPSIS
  RECURSIVE SUBROUTINE add_points_i(pts,indx,k)
!* PURPOSE
!*   Recursively adds points for dual volume calculation
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
!> Recursively adds points for dual volume calculation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)           :: k          !> simplicial order
    INTEGER, INTENT(IN)           :: indx       !> element index
    REAL(KIND=iwp), INTENT(INOUT) :: pts(:,:)   !> bounding points

    !-- local variables --
    INTEGER                       :: i          !> loop counters
    INTEGER                       :: n          !> number of points

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- add to array of points --
    pts(k,:) = lcl_complex(k)%centers(indx,:)

    !-- recursively add points from lower order elements --
    IF (k<dim_cmplx+1) THEN
      DO i = 1, lcl_complex(k)% num_cobndry(indx)
        CALL add_points_i(pts,lcl_complex(k)%cobndry(indx)%indx(i),k+1)
      END DO
    END IF

    RETURN

  END SUBROUTINE
! dec_mod/add_points_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f* dec_mod/calc_unsgnd_vlm
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
    REAL(KIND=iwp), INTENT(IN)  :: pts(:,:)         !> bounding points
    INTEGER, INTENT(IN)         :: n                !> # points to use

    !-- local variables --
    INTEGER                     :: i                !> loop counters
    INTEGER                     :: fac              !> factorial
    REAL(KIND=iwp), ALLOCATABLE :: A(:,:)           !> solution variable
    REAL(KIND=iwp)              :: calc_unsgnd_vlm  !> function value

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
      calc_unsgnd_vlm = SQRT(ABS(determinant(MATMUL(A,TRANSPOSE(A)),n)))/fac

      !-- clean up --
      DEALLOCATE(A)
    END IF

    RETURN

  END FUNCTION
! dec_mod/calc_unsgnd_vlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/exchange_dual_vlm
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
    REAL(KIND=iwp), ALLOCATABLE  :: sbuffer(:)  !> send buffer
    REAL(KIND=iwp), ALLOCATABLE  :: rbuffer(:)  !> receive buffer
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
!/****/s* dec_mod/exchange_dual_dir
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
    REAL(KIND=iwp), ALLOCATABLE  :: sbuffer(:)  !> send buffer
    REAL(KIND=iwp), ALLOCATABLE  :: rbuffer(:)  !> receive buffer
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
!/****/s* dec_mod/calc_hodge_star
!* SYNOPSIS
  SUBROUTINE calc_hodge_star(k)
!* PURPOSE
!*   Compute hodge star and it's inverse from primal and dual volumes
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
!> Compute hodge star and it's inverse from primal and dual volumes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)           :: k        !> geometric order

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate and compute hodge stars and it's invers --
    IF (.NOT. ALLOCATED(lcl_complex(k)%hdg_star)) ALLOCATE(lcl_complex(k)%hdg_star(num_elm(k)))
    IF (.NOT. ALLOCATED(lcl_complex(k)%inv_hdg_star)) ALLOCATE(lcl_complex(k)%inv_hdg_star(num_elm(k)))
    lcl_complex(k)%hdg_star = lcl_complex(k)%dual_volume / lcl_complex(k)%prml_volume
    lcl_complex(k)%inv_hdg_star = lcl_complex(k)%prml_volume / lcl_complex(k)%dual_volume

    RETURN

  END SUBROUTINE
! dec_mod/calc_hodge_star
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_prml_dir
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
! dec_mod/calc_prml_dir
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_dual_dir
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
      DO j = 1,lcl_complex(dim_cmplx)%num_cobndry(i)
        lcl_complex(dim_cmplx)%dual_dir(i,:) = lcl_complex(dim_cmplx)%dual_dir(i,:) + &
          lcl_complex(dim_cmplx)%cobndry(i)%sgn(j)*&
          lcl_complex(dim_cmplx+1)%centers(lcl_complex(dim_cmplx)%cobndry(i)%indx(j),:)
      END DO
    END DO

    !-- add contribution from surface end-points --
    DO i = 1,lcl_complex(dim_cmplx+1)%num_surf
      lcl_complex(dim_cmplx)%dual_dir(lcl_complex(dim_cmplx+1)%surf_indx(i,2),:) = &
        lcl_complex(dim_cmplx)%dual_dir(lcl_complex(dim_cmplx+1)%surf_indx(i,2),:) + &
        lcl_complex(dim_cmplx+1)%surf_indx(i,1)*&
        lcl_complex(dim_cmplx)%centers(lcl_complex(dim_cmplx+1)%surf_indx(i,2),:)
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
! dec_mod/calc_dual_dir
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_whitney_C2_BC
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
    INTEGER           :: vrts(4)    !> volume vertex indices
    INTEGER           :: f_vrts(3)  !> face vertex indices
    INTEGER           :: indx       !> face index
    LOGICAL           :: not_found  !> is boundary face index found
    REAL(KIND=iwp)    :: grads(4,3) !> barycentric gradients

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
      CALL calc_barycentric_grad(grads,4)

      !-- compute contribution from each face of volume --
      DO j = 1,4
        !-- get face vertices --
        f_vrts = (/vrts(1:4-j), vrts(4-j+2:4)/)
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
        IF (not_found) THEN
          !-- check non-local internal faces --
          DO k = 1,lcl_complex(dim_cmplx+1)% num_ext
            indx = lcl_complex(dim_cmplx+1)% ext_indx(k,2)
            IF (ALL(f_vrts == lcl_complex(dim_cmplx)% lcl_node_indx(indx,:))) THEN
              !-- match found; exit search --
              not_found = .FALSE.;   EXIT
            END IF
          END DO
          IF (not_found) THEN
            !-- check surface faces --
            DO k = 1,lcl_complex(dim_cmplx+1)% num_surf
              indx = lcl_complex(dim_cmplx+1)% surf_indx(k,2)
              IF (ALL(f_vrts == lcl_complex(dim_cmplx)% lcl_node_indx(indx,:))) THEN
                !-- match found; exit search --
                not_found = .FALSE.;   EXIT
              END IF
            END DO
          END IF
        END IF
        !-- compute contribution from face at barycenter using Whitney interpolation --
        !-- build index array for computation
        f_vrts = (/(k,k=1,4-j), (k,k=4-j+2,4)/)
        ! write(*,*) '  -',i,j,indx,lcl_complex(dim_cmplx+1)% whtny_sol(i,:),&
        !   0.5d0*lcl_complex(dim_cmplx)% prml_sol(indx,1)*(&
        !   cross_product(grads(f_vrts(1),:),grads(f_vrts(3),:)) + &
        !   cross_product(grads(f_vrts(2),:),grads(f_vrts(1),:)) + &
        !   cross_product(grads(f_vrts(3),:),grads(f_vrts(2),:)) )
        !-- add contribution --
        lcl_complex(dim_cmplx+1)% whtny_sol(i,:) = lcl_complex(dim_cmplx+1)% whtny_sol(i,:) + &
            0.5d0*lcl_complex(dim_cmplx)% prml_sol(indx,1)*(&
            cross_product(grads(f_vrts(1),:),grads(f_vrts(3),:)) + &
            cross_product(grads(f_vrts(2),:),grads(f_vrts(1),:)) + &
            cross_product(grads(f_vrts(3),:),grads(f_vrts(2),:)) )
      END DO
      !write(*,*) '  -',i,lcl_complex(dim_cmplx+1)% whtny_sol(i,:)
    END DO

    RETURN

  END SUBROUTINE
! dec_mod/calc_whitney_C2_BC
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_barycentric_grad
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
    REAL(KIND=iwp), INTENT(INOUT) :: pts(:,:) !> vertices of simplex
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
      !-- not setup for higher dimension simplices --
    END SELECT
    !-- compute last gradient using the fact that the sum of gradients must be zero --
    pts(n,:) = -SUM(pts(1:n-1,:),DIM=1)

    RETURN

  END SUBROUTINE
! dec_mod/calc_barycentric_grad
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! dec_mod
!===============================================================================
!
