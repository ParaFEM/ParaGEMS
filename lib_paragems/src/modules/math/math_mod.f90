!
!===============================================================================
!-- Math Module
!> Module for common math functions and routines
!===============================================================================
!/****/h* modules|math/math_mod
!* SYNOPSIS
MODULE math_mod
!* PURPOSE
!*   Module for common math functions and routines
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!*   mpi_mod                 ???
!* CONTAINS
!*   Subroutine              Purpose
!*   int_merge_sort_rows()   merge sort rows of integer array in ascending order
!*   int_merge_rows()        merge rows of integer array in ascending order
!*   int_insertion_sort()    insertion sort for integer array
!*   any_element_in_list()   test if any element of one array is in another
!*   num_element_in_list()   determine how many elements of one array are in another
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/23: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Module for common math functions and routines
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  math_mod/int_merge_sort_rows
!* SYNOPSIS
  RECURSIVE SUBROUTINE int_merge_sort_rows(A,length,irow,frow,work)
!* PURPOSE
!*   Merge sort rows of integer array in ascending order
!* ASSUMPTION
!*   Rows are already sorted in ascending order
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
!* NOTES
!*   ADAPTED FROM: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Merge sort rows of integer array in ascending order
!> Assumption: rows are already sorted in ascending order
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(INOUT)  :: A(:,:)      !>
    INTEGER                 :: length      !>
    INTEGER, INTENT(IN)     :: irow, frow  !>
    INTEGER, INTENT(INOUT)  :: work(:,:)   !>

    !-- local variables --
    INTEGER                 :: half        !>
    INTEGER                 :: ir,i        !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- sort directly if there are two elements --
    !-- NOTE: if there is only one row, then it is already sorted! --
    IF (length==2) THEN
      !-- if out of order, then swap
      DO ir = irow, frow
        IF (A(1,ir) > A(2,ir)) THEN
          work(1,:) = A(1,:);  A(1,:) = A(2,:);  A(2,:) = work(1,:);  EXIT
        ELSEIF (A(1,ir) < A(2,ir)) THEN;  EXIT
        END IF
      END DO

    !-- for more than two rows, split the array and call merge sort on each half --
    ELSEIF (length>2) THEN
      !-- compute array midpoint --
      half = (length+1)/2
      !-- merge sort half arrays --
      CALL int_merge_sort_rows(A(: half,:),  half,       irow,frow,work)
      CALL int_merge_sort_rows(A(half+1 :,:),length-half,irow,frow,work)
      !-- merge sorted half arrays --
      DO ir = irow, frow
        IF (A(half,ir) > A(half+1,ir)) THEN
          work(1 : half,:) = A(1 : half,:)
          CALL int_merge_rows(work(1 : half,:),A(half+1 : length,:),A,irow,frow)
          EXIT
        ELSEIF (A(half,ir) < A(half+1,ir)) THEN;  EXIT
        END IF
      END DO
    END IF
  END SUBROUTINE
!  math_mod/int_merge_sort_rows
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  math_mod/int_merge_rows
!* SYNOPSIS
  SUBROUTINE int_merge_rows(A,B,C,irow,frow)
!* PURPOSE
!*   Merge rows of an array in ascending order
!* ASSUMPTION
!*   Rows are already sorted in ascending order
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
!* NOTES
!*   ADAPTED FROM: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/23
!> Merge rows of an array in ascending order
!> Assumption: rows are already sorted in ascending order
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)     :: A(:,:),B(:,:)   !>
    INTEGER, INTENT(INOUT)  :: C(:,:)          !>
    INTEGER, INTENT(IN)     :: irow, frow      !>

    !-- local variables --
    LOGICAL                 :: is_le           !>
    INTEGER                 :: i,j,k           !>
    INTEGER                 :: ir              !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- dimension check --
    IF (SIZE(A,DIM=1) + SIZE(B,DIM=1) > SIZE(C,DIM=1)) THEN
      WRITE(*,'(A)')'Error in merge: Arrays sizes are not compatible'
      CALL end_mpi()
    END IF

    !-- loop for all elements of merged array --
    i = 1; j = 1
    DO k = 1, SIZE(C,DIM=1)
      IF (i <= SIZE(A,DIM=1) .and. j <= SIZE(B,DIM=1)) THEN
        is_le = .TRUE.

        DO ir = irow, frow
          IF (A(i,ir) < B(j,ir)) THEN;  EXIT
          ELSEIF (A(i,ir) > B(j,ir)) THEN;  is_le = .FALSE.;  EXIT;  END IF
        END DO

        IF (is_le) THEN;  C(k,:) = A(i,:);  i = i + 1
        ELSE;  C(k,:) = B(j,:);  j = j + 1;  END IF
      ELSEIF (i <= SIZE(A,DIM=1)) THEN;  C(k,:) = A(i,:);  i = i + 1
      ELSEIF (j <= SIZE(B,DIM=1)) THEN;  C(k,:) = B(j,:);  j = j + 1;  END IF
    END DO

  END SUBROUTINE
!  math_mod/int_merge_rows
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s*  math_mod/int_insertion_sort
!* SYNOPSIS
  SUBROUTINE int_insertion_sort(A)
!* PURPOSE
!*   Insertion sort of integer array
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Insertion sort of integer array
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(INOUT)  :: A(:)       !>

    !-- local variables --
    INTEGER                 :: i,j        !>
    INTEGER                 :: work       !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DO i = 2, SIZE(A)
      work = A(i);  j = i - 1
      DO WHILE (j >= 1)
          IF (A(j) <= work) EXIT
          A(j+1) = A(j);  j = j - 1
      END DO
      A(j+1) = work
    END DO

    RETURN

  END SUBROUTINE
!  math_mod/int_insertion_sort
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/any_element_in_list
!* SYNOPSIS
  FUNCTION any_element_in_list(elm,node_list)
!* PURPOSE
!*   Test if any element of one array is in another
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Test if any element of one array is in another
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)   :: elm(:)                !>
    INTEGER, INTENT(IN)   :: node_list(:)          !>
    LOGICAL               :: any_element_in_list   !>

    !-- local variables --
    INTEGER               :: in, jn, end_array     !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- assume it is true --
    any_element_in_list = .TRUE.

    !-- search second array for match --
    end_array = SIZE(node_list,DIM=1)
    DO in = 1,SIZE(elm,DIM=1)
      DO jn = 1,CEILING(end_array/1000.d0)-1
        IF (ANY(elm(in)==node_list(jn*1000-999:jn*1000))) RETURN
      END DO
      IF (ANY(elm(in)==node_list(jn*1000-999:end_array))) RETURN
    END DO

    !-- it is false --
    any_element_in_list = .FALSE.

    RETURN

  END FUNCTION
!  math_mod/any_element_in_list
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/num_element_in_list
!* SYNOPSIS
  FUNCTION num_element_in_list(elm,node_list)
!* PURPOSE
!*   Count number of elements of one array is in another
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Test if any element of one array is in another
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)   :: elm(:)                  !>
    INTEGER, INTENT(IN)   :: node_list(:)            !>
    INTEGER               :: num_element_in_list     !>

    !-- local variables --
    INTEGER               :: in                      !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- search second array for matches --
    num_element_in_list = 0
    DO in = 1,SIZE(elm,DIM=1)
      IF (ANY(elm(in)==node_list)) num_element_in_list = num_element_in_list + 1
    END DO

    RETURN

  END FUNCTION
!  math_mod/num_element_in_list
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/index_in_list
!* SYNOPSIS
  FUNCTION index_in_list(elm,list,init_guess)
!* PURPOSE
!*   Index of element in array
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Test if any element of one array is in another
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    INTEGER, INTENT(IN)   :: elm, list(:)        !>
    INTEGER, INTENT(IN)   :: init_guess          !>
    INTEGER               :: index_in_list       !>

    !-- local variables --
    INTEGER               :: in                  !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- search second array for matches --
    DO in = init_guess,SIZE(list,DIM=1)
      IF (elm==list(in)) THEN;  index_in_list=in;  RETURN;  END IF
    END DO
    DO in = 1,init_guess-1
      IF (elm==list(in)) THEN;  index_in_list=in;  RETURN;  END IF
    END DO
    index_in_list=-1

    RETURN

  END FUNCTION
!  math_mod/index_in_list
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/RCSV_determinant
!* SYNOPSIS
  RECURSIVE FUNCTION RCSV_determinant(A,n) RESULT(accumulation)
!* PURPOSE
!*   Compute determinant
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from:
!*   http://fortranwiki.org/fortran/show/Matrix+inversion
!*   https://rosettacode.org/wiki/Determinant_and_permanent#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Compute determinant
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN) :: A(n,n)    !>
    INTEGER, INTENT(IN)        :: n         !>

    !-- local variables --
    INTEGER               :: i, sgn         !>
    REAL(KIND=PGMSiwp)        :: B(n-1,n-1)     !>
    REAL(KIND=PGMSiwp)        :: accumulation   !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- search second array for matches --
    SELECT CASE (n)
    CASE (1)
      accumulation = A(1,1)
    CASE (2)
      accumulation = (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    CASE (3)
      accumulation = (A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1)))
    CASE DEFAULT
      accumulation = 0;   sgn = 1
      DO i=1, n
        B(:, :(i-1)) = A(2:, :i-1);   B(:, i:) = A(2:, i+1:)
        accumulation = accumulation + sgn * A(1, i) * RCSV_determinant(B, n-1)
        sgn = -sgn
      END DO
    END SELECT

    RETURN

  END FUNCTION
!  math_mod/RCSV_determinant
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/matinv2
!* SYNOPSIS
  FUNCTION matinv2(A) RESULT(B)
!* PURPOSE
!*   Performs a direct calculation of the inverse of a 2x2 matrix.
!* INPUTS
!*   Name                    Description
!* OUTPUTS
!*   Name                    Description
!* SIDE EFFECTS
!*   -
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: http://fortranwiki.org/fortran/show/Matrix+inversion
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Performs a direct calculation of the inverse of a 2x2 matrix.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN) :: A(2,2)   !> Matrix

    !-- local variables --
    REAL(KIND=PGMSiwp)             :: B(2,2)   !> Inverse matrix
    REAL(KIND=PGMSiwp)             :: detinv   !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- Calculate the inverse determinant of the matrix --
    detinv = 1.d0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    !-- Calculate the inverse of the matrix --
    B(1,1) = +detinv * A(2,2);  B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2);  B(2,2) = +detinv * A(1,1)

    RETURN

  END FUNCTION
!  math_mod/matinv2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/matinv3
!* SYNOPSIS
  FUNCTION matinv3(A) RESULT(B)
!* PURPOSE
!*   Performs a direct calculation of the inverse of a 3×3 matrix.
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: http://fortranwiki.org/fortran/show/Matrix+inversion
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Performs a direct calculation of the inverse of a 3×3 matrix.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN) :: A(3,3)   !> Matrix

    !-- local variables --
    REAL(KIND=PGMSiwp)             :: B(3,3)   !> Inverse matrix
    REAL(KIND=PGMSiwp)             :: detinv   !>

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- Calculate the inverse determinant of the matrix --
    detinv = 1.d0/(A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))&
      - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))&
      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1)))

    !-- Calculate the inverse of the matrix --
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    RETURN

  END FUNCTION
!  math_mod/matinv3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/cross_product
!* SYNOPSIS
  FUNCTION cross_product(a, b)
!* PURPOSE
!*   Performs a direct calculation of vector cross product
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: https://rosettacode.org/wiki/Vector_products#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Performs a direct calculation of vector cross product
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN)  :: a(3), b(3)       !> input vectors

    !-- local variables --
    REAL(KIND=PGMSiwp)              :: cross_product(3) !> resulting vector

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute the cross product --
    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - b(1)*a(2)

    RETURN

  END FUNCTION
!  math_mod/cross_product
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/kron_eye
!* SYNOPSIS
  FUNCTION kron_eye(A, n)
!* PURPOSE
!*   Performs a direct calculation of the kronnecker product
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: https://rosettacode.org/wiki/Vector_products#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Performs a direct calculation of the kronnecker product
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN)  :: A(:,:)             !> input matrices
    INTEGER, INTENT(IN)             :: n                  !> input eye size

    !-- local variables --
    REAL(KIND=PGMSiwp)              :: kron_eye(n*size(A,1),n*size(A,2))  !> resulting matrix
    INTEGER                         :: i,j,k,ii,jj        !> counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute the kronnecker product --
    kron_eye = 0.d0
    DO j = 1,size(A,dim=2)
      jj = n*(j-1)
      DO i = 1,size(A,dim=1)
        IF (A(i,j)==0) CYCLE
        ii = n*(i-1)
        DO k = 1,n
          kron_eye(ii+k,jj+k) = A(i,j)
        END DO
      END DO
    END DO

    RETURN

  END FUNCTION
!  math_mod/kron_eye
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/eye_kron
!* SYNOPSIS
  FUNCTION eye_kron(A, n)
!* PURPOSE
!*   Performs a direct calculation of the kronnecker product
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: https://rosettacode.org/wiki/Vector_products#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> Performs a direct calculation of the kronnecker product
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN)  :: A(:,:)             !> input matrices
    INTEGER, INTENT(IN)             :: n                  !> input eye size

    !-- local variables --
    REAL(KIND=PGMSiwp)              :: eye_kron(n*size(A,1),n*size(A,2))  !> resulting matrix
    INTEGER                         :: i,j,k,m,ii,jj        !> counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute the kronnecker product --
    eye_kron = 0.d0
    DO j = 1,n
      jj = size(A,dim=2)*(j-1)
      ii = size(A,dim=1)*(j-1)
      DO m = 1,size(A,dim=2)
        DO k = 1,size(A,dim=1)
          eye_kron(ii+k,jj+m) = A(k,m)
        END DO
      END DO
    END DO

    RETURN

  END FUNCTION
!  math_mod/eye_kron
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/kron_product
!* SYNOPSIS
  FUNCTION norm(v)
!* PURPOSE
!*   calculates norm of vector
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: https://rosettacode.org/wiki/Vector_products#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> calculates norm of vector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN)  :: v(:)             !> input vector

    !-- local variables --
    REAL(KIND=PGMSiwp)              :: norm  !> resulting norm
    INTEGER                         :: i        !> counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute the kronnecker product --
    norm = 0.d0
    DO i = 1,size(v)
      norm = norm + v(i)**2
    END DO
    norm = sqrt(norm)

    RETURN

  END FUNCTION
!  math_mod/norm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/f*  math_mod/kron_product
!* SYNOPSIS
  FUNCTION diagM(v)
!* PURPOSE
!*   creates diagonal matrix
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
!*   2019/08/29: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!* NOTES
!*   Adapted from: https://rosettacode.org/wiki/Vector_products#Fortran
!******/
!> author: Pieter Boom
!> date: 2019/08/29
!> creates diagonal matrix
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=PGMSiwp), INTENT(IN)  :: v(:)             !> input vector

    !-- local variables --
    REAL(KIND=PGMSiwp)              :: diagM(size(v),size(v))  !> resulting matrix
    INTEGER                         :: i        !> counters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- compute the kronnecker product --
    diagM = 0.d0
    DO i = 1,size(v)
      diagM(i,i) = v(i)
    END DO

    RETURN

  END FUNCTION
!  math_mod/diagM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! math_mod
!===============================================================================
!
