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
!*   count_bndry_cobndry     count [co-]boundaries
!*   set_bndry_cobndry       set [co-]boundaries
!*   calc_hodge_star         compute hodge star and it's inverse from primal and dual volumes
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
  SUBROUTINE calc_orientation(elm,orient)
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
!*  -the nodal indices of elm are sorted numerically in ascending order
!*  -orientation contains the ±1 orientation of the element
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
    INTEGER, INTENT(INOUT)  :: elm(:)   !> nodal indices of elements
    INTEGER, INTENT(INOUT)  :: orient   !> orientations of the elements

    !-- local variables --
    INTEGER                :: i,j              !> loop counters
    INTEGER                :: length           !> size of elm array
    INTEGER                :: work             !> work variable

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-------------------------------------------------------------------------
    !-- sort indices using insertion sort and count the number of swaps --
    !-------------------------------------------------------------------------
    orient=0; length=SIZE(elm,DIM=1)
    DO i=2, length !-- loop through all indices, starting from the second --
      !-- get next value to insert and set loop index --
      work=elm(i); j=i-1
      DO WHILE (j >= 1) !-- find location to insert given value --
        IF (elm(j) <= work) EXIT !-- ??? found location ??? --
        !-- shuffle (swap) values, update loop index and count swaps --
        elm(j+1)=elm(j); j=j-1; orient=orient+1; END DO
      elm(j+1)=work; END DO !-- insert value --
    !-- compute orientation from number of swaps --
    orient=1-2*MOD(orient,2); RETURN

  END SUBROUTINE
! dec_mod/calc_orientation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/calc_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE calc_bndry_cobndry(k)
!* PURPOSE
!*   Recursively compute element (co-)boundaries from highest to lowest geometric order
!* ASSUMPTION
!*   Mesh is a simplicial complex in TetGen format
!* SIDE EFFECTS
!*  -the nodal indices of elm are sorted numerically in ascending order
!*  -orientation contains the ±1 orientation of the element
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

    !-- arguments --
    INTEGER, INTENT(IN) :: k  !> simplicial order

    !-- local variables --
    INTEGER, ALLOCATABLE  :: bndry(:,:)     !> temp. boundary work array
    INTEGER, ALLOCATABLE  :: bndry_cnt(:)   !> counter local boundaries
    INTEGER, ALLOCATABLE  :: cobndry_cnt(:) !> counter co-boundaries
    INTEGER               :: cnt            !> counter boundaries
    INTEGER               :: ext_cnt(2)     !> counter external boundaries
    INTEGER, ALLOCATABLE  :: work(:,:)      !> work array for merge sort
    INTEGER, ALLOCATABLE  :: lwork(:)       !> work array to determine locality

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- calculate signed adjacency [co-boundary operator] --
    !-- initial guess for the number of boundaries --
    num_elm(k-1)=k*num_elm(k)

    !-- build working array with element/boundary data --
    ALLOCATE(bndry(num_elm(k-1),k+2)); CALL build_bndry_work_array(bndry,k)

    !-- sort working array --
    ALLOCATE(work(num_elm(k-1),k+2))
    CALL int_merge_sort_rows(bndry,num_elm(k-1),1,k-1,work); DEALLOCATE(work)

    !-- count number of unique (local/external) [co-]boundaries --
    ALLOCATE(bndry_cnt(num_elm(k)),cobndry_cnt(num_elm(k-1)),&
      lwork(num_elm(k)*k))
    CALL count_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,lwork)

    !-- allocate [co-]boundary structures and variables --
    CALL allocate_bndry_cobndry(k,cnt,bndry_cnt,cobndry_cnt,ext_cnt)

    !-- setup boundaries external to the current process (not local),
    !   setup parallel mapping, setup node indices for (k-1)^th order element
    !   structure, and setup co-boundaries for (k-1)^th order element
    !   structure --
    CALL set_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,lwork)

    !-- clean up temporary variables --
    DEALLOCATE(bndry,bndry_cnt,cobndry_cnt,lwork); RETURN

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
    iib=1; fib=num_elm(k); DO i=1,k
      !-- get boundary node indices --
      bndry(iib:fib,1:i-1)=lcl_complex(k)%node_indx(:,1:i-1)
      bndry(iib:fib,i:k-1)=lcl_complex(k)%node_indx(:,i+1:k)

      !-- get simplex indices --
      bndry(iib:fib,k)  =(/ (j, j=1,num_elm(k)) /); bndry(iib:fib,k+1)=i

      !-- get orientation --
      IF (k == dim_cmplx+1) THEN
        bndry(iib:fib,k+2) =(-1)**i*lcl_complex(k)%orientation
      ELSE; bndry(iib:fib,k+2) =(-1)**i; END IF

      !-- update array indices --
      iib=iib+num_elm(k); fib=fib+num_elm(k)
    END DO; RETURN

  END SUBROUTINE
! dec_mod/build_bndry_work_array
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/count_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE count_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,local)
!* PURPOSE
!*   Count [co-]boundaries: internal, external
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
!> Count [co-]boundaries: internal, external
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)     :: k                !> simplicial order
    INTEGER, INTENT(IN)     :: bndry(:,:)       !> boundary work array
    INTEGER, INTENT(INOUT)  :: bndry_cnt(:)     !> counter local boundaries
    INTEGER, INTENT(INOUT)  :: cobndry_cnt(:)   !> counter co-boundaries
    INTEGER, INTENT(INOUT)  :: cnt              !> counter (k-1) element
    INTEGER, INTENT(INOUT)  :: ext_cnt(2)       !> counter external boundaries
    INTEGER, INTENT(INOUT)  :: local(:)         !> locality work array

    !-- local variables --
    INTEGER :: i,j,np       !> counter
    INTEGER :: procs(k-1)   !> proc ids
    INTEGER :: procs2(30)   !> proc ids

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- reset counters --
    cnt=0; bndry_cnt=0; cobndry_cnt=0; ext_cnt=0

    !-- loop for all boundaries --
    DO i=1, num_elm(k)*k
      bndry_cnt(bndry(i,k))=bndry_cnt(bndry(i,k))+1

      IF (i == 1 .OR. ANY(bndry(i,1:k-1) /= bndry(max(1,i-1),1:k-1))) THEN
        !-- boundary different from previous: new boundary --
        cnt=cnt+1
        ! procs = elm2proc_array(bndry(i,1:k-1))
        ! IF (procs(1)==rank) THEN
        !   !-- boundary is local & sender: count co-boundaries --
        !   cobndry_cnt(cnt)=cobndry_cnt(cnt)+1; np=1
        !   procs2(1) = elm2proc(lcl_complex(k)%node_indx(bndry(i,k),bndry(i,k+1)))
        !   IF (ALL(procs/=procs2(1))) THEN; local(i)=-2; ext_cnt(2)=ext_cnt(2)+1
        !   ELSE; local(i)=-1; END IF
        ! ELSEIF (ANY(procs(2:k-1)==rank)) THEN
        IF (ANY(elm2proc_array(bndry(i,1:k-1))==rank)) THEN
          !-- boundary is local: count co-boundaries --
          local(i)=1; cobndry_cnt(cnt)=cobndry_cnt(cnt)+1
        ELSE !-- boundary is non-local (recv'r): count external boundaries --
          local(i)=0; ext_cnt(1)=ext_cnt(1)+1
        END IF
      ELSE
        !-- boundary same as previous --
        ! IF (local(i-1)<0) THEN !-- boundary is local: count co-boundaries --
        !   cobndry_cnt(cnt)=cobndry_cnt(cnt)+1; np=np+1
        !   procs2(np) = elm2proc(lcl_complex(k)%node_indx(bndry(i,k),bndry(i,k+1)))
        !   IF (ALL(procs/=procs2(np)) .AND. ALL(procs2(1:np-1)/=procs2(np))) THEN
        !     local(i)=-2; ext_cnt(2)=ext_cnt(2)+1
        !   ELSE; local(i)=-1; END IF
        ! ELSEIF (local(i-1)==1) THEN
        IF (local(i-1)/=0) THEN
          local(i)=1; cobndry_cnt(cnt)=cobndry_cnt(cnt)+1
        ELSE; local(i)=0; END IF
      END IF
    END DO; RETURN

  END SUBROUTINE
! dec_mod/count_bndry_cobndry
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/allocate_bndry_cobndry
!* SYNOPSIS
  SUBROUTINE allocate_bndry_cobndry(k,cnt,bndry_cnt,cobndry_cnt,ext_cnt)
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

    !-- local variables --
    INTEGER                 :: i                !> counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- allocate boundary structures and variables: k^th order --
    !-- # external boundaries and # expect MPI receives --
    lcl_complex(k)%num_recv=ext_cnt(1)
    ! lcl_complex(k)%num_send=ext_cnt(2)
    ALLOCATE(lcl_complex(k)%recv_indx(ext_cnt(1),2))
    ! ALLOCATE(lcl_complex(k)%send_indx(ext_cnt(2),2))

    !-- # internal boundaries --
    ALLOCATE(&
      lcl_complex(k)%num_bndry(num_elm(k)),&
      lcl_complex(k)%bndry(num_elm(k)))
    lcl_complex(k)%num_bndry=bndry_cnt
    DO i=1,num_elm(k); ALLOCATE(&
        lcl_complex(k)%bndry(i)%sgn(lcl_complex(k)%num_bndry(i)),&
        lcl_complex(k)%bndry(i)%indx(lcl_complex(k)%num_bndry(i)))
    END DO

    !-- allocate [co-]boundary structures and variables: (k-1))^th order --
    !-- # (k-1))^th order elements --
    num_elm(k-1)=cnt; ALLOCATE(lcl_complex(k-1)%node_indx(num_elm(k-1),k-1))
    !-- # (k-1))^th order coboundaries --
    ALLOCATE(&
      lcl_complex(k-1)%num_cobndry(num_elm(k-1)),&
      lcl_complex(k-1)%cobndry(num_elm(k-1)))
    lcl_complex(k-1)%num_cobndry=cobndry_cnt(1:num_elm(k-1))
    DO i=1,num_elm(k-1); ALLOCATE(&
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
  SUBROUTINE set_bndry_cobndry(bndry,k,cnt,bndry_cnt,cobndry_cnt,ext_cnt,local)
!* PURPOSE
!*   Set [co-]boundaries: internal, external
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
!> Set [co-]boundaries: internal, external
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)     :: k                !> simplicial order
    INTEGER, INTENT(INOUT)  :: bndry(:,:)       !> boundary work array
    INTEGER, INTENT(INOUT)  :: bndry_cnt(:)     !> temp # local boundaries
    INTEGER, INTENT(INOUT)  :: cobndry_cnt(:)   !> temp # co-boundaries
    INTEGER, INTENT(INOUT)  :: cnt              !> counter (k-1) element
    INTEGER, INTENT(INOUT)  :: ext_cnt(2)       !> counter external boundaries
    INTEGER, INTENT(INOUT)  :: local(:)         !> locality work array

    !-- local variables --
    INTEGER                 :: i                !> counter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- reset counters --
    cnt=0; bndry_cnt=0; cobndry_cnt=0; ext_cnt=0

    !-- loop for all boundaries --
    DO i=1, num_elm(k)*k
      IF (i == 1 .OR. ANY(bndry(i,1:k-1) /= bndry(max(1,i-1),1:k-1))) THEN
        !-- boundary different from previous: new boundary --
        cnt=cnt+1

        IF (local(i)/=0) THEN
          !-- boundary is local: set signed adjacency for co-boundaries --
          cobndry_cnt(cnt)=1
          lcl_complex(k-1)%cobndry(cnt)%sgn(1)=bndry(i,k+2)
          lcl_complex(k-1)%cobndry(cnt)%indx(1)=bndry(i,k)
        ELSE
          !-- boundary is non-local: set signed adjacency for external boundaries --
          ext_cnt(1)=ext_cnt(1)+1
          lcl_complex(k)%recv_indx(ext_cnt(1),:)=(/ elm2proc(bndry(i,1)), cnt /)                          !????? find actual global index
        END IF
        !-- set node indices for (k-1)^th order elements --
        lcl_complex(k-1)%node_indx(cnt,:)=bndry(i,1:k-1)
      ELSEIF (local(i)/=0) THEN !-- boundary same as previous and local --
        !-- set signed adjacency for co-boundaries for each local boundary --
        cobndry_cnt(cnt)=cobndry_cnt(cnt)+1
        lcl_complex(k-1)%cobndry(cnt)%sgn(cobndry_cnt(cnt))=bndry(i,k+2)
        lcl_complex(k-1)%cobndry(cnt)%indx(cobndry_cnt(cnt))=bndry(i,k)
      END IF

      ! IF (local(i)==-2) THEN
      !   ext_cnt(2)=ext_cnt(2)+1
      !   lcl_complex(k)%send_indx(ext_cnt(2),:)=&
      !     (/ elm2proc(lcl_complex(k)%node_indx(bndry(i,k),bndry(i,k+1))), cnt /)
      ! END IF

      bndry_cnt(bndry(i,k))=bndry_cnt(bndry(i,k))+1
      lcl_complex(k)%bndry(bndry(i,k))%sgn(bndry_cnt(bndry(i,k)))=bndry(i,k+2)
      lcl_complex(k)%bndry(bndry(i,k))%indx(bndry_cnt(bndry(i,k)))=cnt
    END DO

    RETURN

  END SUBROUTINE
! dec_mod/set_bndry_cobndry
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* dec_mod/initialise_hodge_star
!* SYNOPSIS
  SUBROUTINE initialise_hodge_star
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
    INTEGER :: k        !> geometric order

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute Hodge stars --
    DO k=1,dim_cmplx+1; CALL calc_hodge_star(k); END DO; RETURN

  END SUBROUTINE
! dec_mod/initialise_hodge_star
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
    !-- allocate and compute hodge stars and it's inverse --
    IF (.NOT. ALLOCATED(lcl_complex(k)%hdg_star))&
      ALLOCATE(lcl_complex(k)%hdg_star(num_elm(k)))
    IF (.NOT. ALLOCATED(lcl_complex(k)%inv_hdg_star))&
      ALLOCATE(lcl_complex(k)%inv_hdg_star(num_elm(k)))
    lcl_complex(k)%hdg_star=&
      lcl_complex(k)%dual_volume/lcl_complex(k)%prml_volume
    lcl_complex(k)%inv_hdg_star=&
      lcl_complex(k)%prml_volume/lcl_complex(k)%dual_volume
    RETURN

  END SUBROUTINE
! dec_mod/calc_hodge_star
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! dec_mod
!===============================================================================
!
