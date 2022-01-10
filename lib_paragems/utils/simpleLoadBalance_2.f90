!
!===============================================================================
! /****p*/ /utils/simpleLoadBalance
PROGRAM simpleLoadBalance
!
! PURPOSE:  miniapp to to reorder nodes file for 'better' partitioning/load balancing
!
! PRE:      - Mesh files exists (TetGen tetrahedral mesh)
!
! POST:     - Mesh files reordered for 'better' partitioning/load balancing
!
! UPDATES:  created (PDB) :: 2020/06/02
!
!===============================================================================

  !-- implicit none --
  IMPLICIT NONE

  !-- real precision --
  INTEGER, PARAMETER          :: iwp = SELECTED_REAL_KIND(15,300)

  !-- local variables --
  CHARACTER(LEN=128)          :: mesh_prefix      ! file name prefix
  CHARACTER(LEN=128)          :: fname            ! file name
  INTEGER                     :: unit=10          ! IO unit
  INTEGER                     :: junk             ! junk IO variable
  INTEGER                     :: offset           ! junk IO variable
  INTEGER                     :: ie,je,ke,kn      ! loop indices
  INTEGER, ALLOCATABLE        :: elems(:,:)       ! node indices of elements
  REAL(KIND=iwp), ALLOCATABLE :: nodes(:,:)       ! xyz coordinates of nodes
  REAL(KIND=iwp), ALLOCATABLE :: centers(:,:)     ! xyz coordinates of element baryceneters
  INTEGER                     :: num_elm, num_node! number of elements and nodes
  INTEGER                     :: ier              ! error status variable
  REAL(KIND=iwp), ALLOCATABLE :: work(:,:)        ! work array for merge sort
  INTEGER, ALLOCATABLE        :: work_int(:,:)        ! work array for merge sort
  INTEGER, ALLOCATABLE        :: indx(:)          ! indx

!===============================================================================
! MAIN EXECUTION
!===============================================================================
  !-- write miniapp details --
  write(*,*) 'ParaGEMS utility - reorder_ele: reorder .ele file for better partition/load balancing'

  !-- get file prefix --
  write(*,*) 'What is the mesh prefix for the .node and .ele files?'
  read(*,*) mesh_prefix

  !-- open node file and check for errors --
  fname = trim(mesh_prefix) // ".node"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- read header --
  write(*,*) "reading file: ",fname
  READ(unit,*) num_node, kn
  !-- get offset --
  READ(unit,*) offset
  BACKSPACE(unit)
  !-- read data --
  ALLOCATE(nodes(num_node,kn+1))
  DO ie=1,num_node
    READ(unit,*) nodes(ie,:)
  END DO
  !-- close file --
  CLOSE(unit)

  !-- sort nodes --
  write(*,*) "sorting nodes"
  ALLOCATE(work(num_node,kn+1))
  CALL merge_sort_rows(nodes,num_node,2,kn+1,work)
  DEALLOCATE(work)

  !-- rename files --
  fname = trim(mesh_prefix)//".node"
  write(*,*) "renaming file: ",fname
  CALL RENAME(trim(mesh_prefix) // ".node",trim(mesh_prefix) // ".node.old",ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error renaming ',trim(mesh_prefix) // ".node",': errcode = ',ier
    STOP
  END IF

  !-- write new node file --
  fname = trim(mesh_prefix)//".node"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- write header --
  write(*,*) "writing file: ",fname
  WRITE(unit,*) num_node, kn, 0, 0
  !-- write data --
  DO ie=1,num_node
    !write(*,*) ie, nodes(ie,2:kn+1)
    WRITE(unit,*) ie, nodes(ie,2:kn+1)
  END DO

  !-- get inverse node sort index --
  write(*,*) "inverse sort"
  ALLOCATE(indx(num_node))
  DO ie=1,num_node
    indx(NINT(nodes(ie,1))) = ie
  END DO

  !-- open element file and check for errors --
  fname = trim(mesh_prefix) // ".ele"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- read header --
  write(*,*) "reading file: ",fname
  READ(unit,*) num_elm, ke
  !-- read data --
  ALLOCATE(elems(num_elm,ke))
  DO ie=1,num_elm
    READ(unit,*) junk, elems(ie,:)
  END DO
  IF (offset/=1) elems = elems - (offset-1)
  !-- close file --
  CLOSE(unit)

  !-- update element indicies --
  write(*,*) "index update"
  DO ie=1,num_elm
    DO je=1,ke
      elems(ie,je) = indx(elems(ie,je))
    END DO
  END DO

  !-- sort elements --
  write(*,*) "sorting elements"
  ALLOCATE(work_int(num_elm,ke))
  CALL merge_sort_rows_int(elems,num_elm,1,ke,work_int)
  DEALLOCATE(work_int)


  !-- rename files --
  write(*,*) "renaming file: ",fname
  CALL RENAME(trim(mesh_prefix) // ".ele",trim(mesh_prefix) // ".ele.old",ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error renaming ',trim(mesh_prefix) // ".ele",': errcode = ',ier
    STOP
  END IF

  !-- write new ele file --
  fname = trim(mesh_prefix) // ".ele"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- write header --
  write(*,*) "writing file: ",fname
  WRITE(unit,*) num_elm, ke, 0
  !-- write data --
  DO ie=1,num_elm
    WRITE(unit,*) ie, elems(ie,:)
  END DO

  !-- close file  and clean up --
  write(*,*) "close and clean"
  CLOSE(unit)
  DEALLOCATE(elems)

  !-- open element file and check for errors --
  fname = trim(mesh_prefix) // ".face"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- read header --
  write(*,*) "reading file: ",fname
  READ(unit,*) num_elm, ke
  ke=3
  !-- read data --
  ALLOCATE(elems(num_elm,ke+1))
  DO ie=1,num_elm
    READ(unit,*) junk, elems(ie,:)
  END DO
  IF (offset/=1) elems(:,1:ke) = elems(:,1:ke) - (offset-1)
  !-- close file --
  CLOSE(unit)

  !-- update element indicies --
  write(*,*) "index update"
  DO ie=1,num_elm
    DO je=1,ke
      elems(ie,je) = indx(elems(ie,je))
    END DO
  END DO

  !-- sort faces --
  write(*,*) "sorting faces"
  ALLOCATE(work_int(num_elm,ke+1))
  CALL merge_sort_rows_int(elems,num_elm,1,ke,work_int)
  DEALLOCATE(work_int)

  !-- rename files --
  write(*,*) "renaming file: ",fname
  CALL RENAME(trim(mesh_prefix) // ".face",trim(mesh_prefix) // ".face.old",ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error renaming ',trim(mesh_prefix) // ".face",': errcode = ',ier
    STOP
  END IF

  !-- write new ele file --
  fname = trim(mesh_prefix) // ".face"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- write header --
  write(*,*) "writing file: ",fname
  WRITE(unit,*) num_elm, 1
  !-- write data --
  DO ie=1,num_elm
    WRITE(unit,*) ie, elems(ie,:)
  END DO

  !-- close file  and clean up --
  write(*,*) "close and clean"
  CLOSE(unit)
  DEALLOCATE(elems)

  !-- open element file and check for errors --
  fname = trim(mesh_prefix) // ".edge"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- read header --
  write(*,*) "reading file: ",fname
  READ(unit,*) num_elm, ke
  ke = 2
  !-- read data --
  ALLOCATE(elems(num_elm,ke+1))
  DO ie=1,num_elm
    READ(unit,*) junk, elems(ie,:)
  END DO
  IF (offset/=1) elems(:,1:ke) = elems(:,1:ke) - (offset-1)
  !-- close file --
  CLOSE(unit)

  !-- update element indicies --
  write(*,*) "index update"
  DO ie=1,num_elm
    DO je=1,ke
      elems(ie,je) = indx(elems(ie,je))
    END DO
  END DO

  !-- sort edges --
  write(*,*) "sorting edges"
  ALLOCATE(work_int(num_elm,ke+1))
  CALL merge_sort_rows_int(elems,num_elm,1,ke,work_int)
  DEALLOCATE(work_int)

  !-- rename files --
  write(*,*) "renaming file: ",fname
  CALL RENAME(trim(mesh_prefix) // ".edge",trim(mesh_prefix) // ".edge.old",ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error renaming ',trim(mesh_prefix) // ".edge",': errcode = ',ier
    STOP
  END IF

  !-- write new ele file --
  fname = trim(mesh_prefix) // ".edge"
  write(*,*) "opening file: ",fname
  OPEN(unit,FILE=fname,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ier)
  IF (ier /= 0) THEN
    WRITE(*,'(A,A,A,I5)')'Error opening ',trim(fname),': errcode = ',ier
    STOP
  END IF

  !-- write header --
  write(*,*) "writing file: ",fname
  WRITE(unit,*) num_elm, 1
  !-- write data --
  DO ie=1,num_elm
    WRITE(unit,*) ie, elems(ie,:)
  END DO

  !-- close file  and clean up --
  write(*,*) "close and clean"
  CLOSE(unit)
  DEALLOCATE(elems)

  !-- clean up --
  write(*,*) "clean"
  DEALLOCATE(nodes,indx)

!===============================================================================
CONTAINS
!===============================================================================
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! /****s*/ /src/modules/math_mod|int_merge_sort_rows
  RECURSIVE SUBROUTINE merge_sort_rows(A,length,irow,frow,work)
  !
  ! PURPOSE:  merge sort rows of array in ascending order
  !
  ! PRE:      rows are already sorted in ascending order
  !
  ! UPDATES:  created (PDB) :: 2019/08/23
  !
  ! ADAPTED FROM: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=iwp), INTENT(INOUT)  :: A(:,:)
    INTEGER, INTENT(IN)     :: length
    INTEGER, INTENT(IN)     :: irow, frow
    REAL(KIND=iwp), INTENT(INOUT)  :: work(:,:)

    !-- local variables --
    INTEGER                 :: half
    INTEGER                 :: ir,i

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! MAIN EXECUTION
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !-- sort directly if there are two elements --
    !-- NOTE: if there is only one row, then it is already sorted! --
    IF (length==2) THEN
      !-- if out of order, then swap
      DO ir = irow, frow
        IF (A(1,ir) > A(2,ir)) THEN
          work(1,:) = A(1,:)
          A(1,:) = A(2,:)
          A(2,:) = work(1,:)
          EXIT
        ELSEIF (A(1,ir) < A(2,ir)) THEN
          EXIT
        END IF
      END DO

    !-- for more than two rows, split the array and call merge sort on each half --
    ELSEIF (length>2) THEN
      !-- compute array midpoint --
      half = (length+1)/2

      !-- merge sort half arrays --
      CALL merge_sort_rows(A(: half,:),  half,       irow,frow,work)
      CALL merge_sort_rows(A(half+1 :,:),length-half,irow,frow,work)

      !-- merge sorted half arrays --
      DO ir = irow, frow
        IF (A(half,ir) > A(half+1,ir)) THEN
          work(1 : half,:) = A(1 : half,:)
          CALL merge_rows(work(1 : half,:),A(half+1 : length,:),A,irow,frow)
          EXIT
        ELSEIF (A(half,ir) < A(half+1,ir)) THEN
          EXIT
        END IF
      END DO
    END IF
  END SUBROUTINE
  ! math_mod|int_merge_sort_rows
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! /****s*/ /src/modules/math_mod|int_merge_rows
  SUBROUTINE merge_rows(A,B,C,irow,frow)
  !
  ! PURPOSE:  merge rows of an array in ascending order
  !
  ! PRE:      rows are already sorted in ascending order
  !
  ! UPDATES:  created (PDB) :: 2019/08/23
  !
  ! ADAPTED FROM: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE

    !-- arguments --
    REAL(KIND=iwp), INTENT(IN)     :: A(:,:),B(:,:)
    REAL(KIND=iwp), INTENT(INOUT)  :: C(:,:)
    INTEGER, INTENT(IN)     :: irow, frow

    !-- local variables --
    LOGICAL                 :: is_le
    INTEGER                 :: i,j,k
    INTEGER                 :: ir

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! MAIN EXECUTION
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- dimension check --
    IF (SIZE(A,DIM=1) + SIZE(B,DIM=1) > SIZE(C,DIM=1)) THEN
      STOP
    END IF

    !-- loop for all elements of merged array --
    i = 1; j = 1
    DO k = 1, SIZE(C,DIM=1)
      IF (i <= SIZE(A,DIM=1) .and. j <= SIZE(B,DIM=1)) THEN
        is_le = .TRUE.

        DO ir = irow, frow
          IF (A(i,ir) < B(j,ir)) THEN
            EXIT
          ELSEIF (A(i,ir) > B(j,ir)) THEN
            is_le = .FALSE.
            EXIT
          END IF
        END DO

        IF (is_le) THEN
          C(k,:) = A(i,:)
          i = i + 1
        ELSE
          C(k,:) = B(j,:)
          j = j + 1
        END IF
      ELSEIF (i <= SIZE(A,DIM=1)) THEN
        C(k,:) = A(i,:)
        i = i + 1
      ELSEIF (j <= SIZE(B,DIM=1)) THEN
        C(k,:) = B(j,:)
        j = j + 1
      END IF
    END DO

  END SUBROUTINE
  ! math_mod|int_merge_rows
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
    !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! /****s*/ /src/modules/math_mod|int_merge_sort_rows
    RECURSIVE SUBROUTINE merge_sort_rows_int(A,length,irow,frow,work)
    !
    ! PURPOSE:  merge sort rows of array in ascending order
    !
    ! PRE:      rows are already sorted in ascending order
    !
    ! UPDATES:  created (PDB) :: 2019/08/23
    !
    ! ADAPTED FROM: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
    !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE

      !-- arguments --
      INTEGER, INTENT(INOUT)  :: A(:,:)
      INTEGER, INTENT(IN)     :: length
      INTEGER, INTENT(IN)     :: irow, frow
      INTEGER, INTENT(INOUT)  :: work(:,:)

      !-- local variables --
      INTEGER                 :: half
      INTEGER                 :: ir,i

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! MAIN EXECUTION
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !-- sort directly if there are two elements --
      !-- NOTE: if there is only one row, then it is already sorted! --
      IF (length==2) THEN
        IF (minval(A(1,:)) > minval(A(2,:))) THEN
          work(1,:) = A(1,:)
          A(1,:) = A(2,:)
          A(2,:) = work(1,:)
        ELSEIF (minval(A(1,:)) == minval(A(2,:))) THEN
          !-- if out of order, then swap
          DO ir = irow, frow
            IF (A(1,ir) > A(2,ir)) THEN
              work(1,:) = A(1,:)
              A(1,:) = A(2,:)
              A(2,:) = work(1,:)
              EXIT
            ELSEIF (A(1,ir) < A(2,ir)) THEN
              EXIT
            END IF
          END DO
        END IF

      !-- for more than two rows, split the array and call merge sort on each half --
      ELSEIF (length>2) THEN
        !-- compute array midpoint --
        half = (length+1)/2

        !-- merge sort half arrays --
        CALL merge_sort_rows_int(A(: half,:),  half,       irow,frow,work)
        CALL merge_sort_rows_int(A(half+1 :,:),length-half,irow,frow,work)

        !-- merge sorted half arrays --
        IF (minval(A(half,:)) > minval(A(half+1,:))) THEN
          work(1 : half,:) = A(1 : half,:)
          CALL merge_rows_int(work(1 : half,:),A(half+1 : length,:),A,irow,frow)
        ELSEIF (minval(A(half,:)) == minval(A(half+1,:))) THEN
          DO ir = irow, frow
            IF (A(half,ir) > A(half+1,ir)) THEN
              work(1 : half,:) = A(1 : half,:)
              CALL merge_rows_int(work(1 : half,:),A(half+1 : length,:),A,irow,frow)
              EXIT
            ELSEIF (A(half,ir) < A(half+1,ir)) THEN
              EXIT
            END IF
          END DO
        END IF
      END IF
    END SUBROUTINE
    ! math_mod|int_merge_sort_rows
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! /****s*/ /src/modules/math_mod|int_merge_rows
    SUBROUTINE merge_rows_int(A,B,C,irow,frow)
    !
    ! PURPOSE:  merge rows of an array in ascending order
    !
    ! PRE:      rows are already sorted in ascending order
    !
    ! UPDATES:  created (PDB) :: 2019/08/23
    !
    ! ADAPTED FROM: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
    !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE

      !-- arguments --
      INTEGER, INTENT(IN)     :: A(:,:),B(:,:)
      INTEGER, INTENT(INOUT)  :: C(:,:)
      INTEGER, INTENT(IN)     :: irow, frow

      !-- local variables --
      LOGICAL                 :: is_le
      INTEGER                 :: i,j,k
      INTEGER                 :: ir

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! MAIN EXECUTION
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !-- dimension check --
      IF (SIZE(A,DIM=1) + SIZE(B,DIM=1) > SIZE(C,DIM=1)) THEN
        STOP
      END IF

      !-- loop for all elements of merged array --
      i = 1; j = 1
      DO k = 1, SIZE(C,DIM=1)
        IF (i <= SIZE(A,DIM=1) .and. j <= SIZE(B,DIM=1)) THEN
          is_le = .TRUE.

          IF (minval(A(i,:)) > minval(B(j,:))) THEN
            is_le = .FALSE.
          ELSEIF (minval(A(i,:)) == minval(B(j,:))) THEN
            DO ir = irow, frow
              IF (A(i,ir) < B(j,ir)) THEN
                EXIT
              ELSEIF (A(i,ir) > B(j,ir)) THEN
                is_le = .FALSE.
                EXIT
              END IF
            END DO
          END IF

          IF (is_le) THEN
            C(k,:) = A(i,:)
            i = i + 1
          ELSE
            C(k,:) = B(j,:)
            j = j + 1
          END IF
        ELSEIF (i <= SIZE(A,DIM=1)) THEN
          C(k,:) = A(i,:)
          i = i + 1
        ELSEIF (j <= SIZE(B,DIM=1)) THEN
          C(k,:) = B(j,:)
          j = j + 1
        END IF
      END DO

    END SUBROUTINE
    ! math_mod|int_merge_rows
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !
END PROGRAM
! simpleLoadBalance
!===============================================================================
!
