!
!===============================================================================
!-- time_marching Module
!> Module time_marching structures and code
!===============================================================================
!/****/h* modules|time_marching/time_marching_mod
!* SYNOPSIS
MODULE time_marching_mod
!* PURPOSE
!*   Module for global variables, included library, and time_marching code
!* INCLUDES
!*   Name                  Purpose
!* CONTAINS
!*   Subroutine              Purpose
!*
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/20: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/20
!> Module time_marching structures and code
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod
  USE mpi_mod

  IMPLICIT NONE

#include "time_marching.inc"

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* time_marching_mod/initialise_time_marching
!* SYNOPSIS
  SUBROUTINE initialise_time_marching()
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !-- deallocate primary elements --

    IF (ALLOCATED(time_marching%A)) DEALLOCATE(time_marching%A)
    IF (ALLOCATED(time_marching%B)) DEALLOCATE(time_marching%B)
    IF (ALLOCATED(time_marching%U)) DEALLOCATE(time_marching%U)
    IF (ALLOCATED(time_marching%V)) DEALLOCATE(time_marching%V)
    IF (ALLOCATED(time_marching%c)) DEALLOCATE(time_marching%c)
    IF (ALLOCATED(time_marching%dc)) DEALLOCATE(time_marching%dc)
    IF (ALLOCATED(time_marching%LHS_tvec)) DEALLOCATE(time_marching%LHS_tvec)
    IF (ALLOCATED(time_marching%q_step_indx)) DEALLOCATE(time_marching%q_step_indx)

    SELECTCASE(time_marching_method)

    CASE('explicit_euler','ee')
      time_marching%implicit = 0
      time_marching%p = 1;    time_marching%q = 1
      time_marching%r = 1;    time_marching%s = 1;    time_marching%s_initial = 1

      ALLOCATE(time_marching%A(1,1),time_marching%U(1,1),time_marching%B(1,1),&
        time_marching%V(1,1),time_marching%c(1),time_marching%q_step_indx(1))

      time_marching%A = transpose(reshape((/ 0.d0 /),(/ time_marching%s,time_marching%s /)))
      time_marching%B = transpose(reshape((/ 1.d0 /),(/ time_marching%s,time_marching%r /)))
      time_marching%U = transpose(reshape((/ 1.d0 /),(/ time_marching%r,time_marching%s /)))
      time_marching%V = transpose(reshape((/ 1.d0 /),(/ time_marching%r,time_marching%r /)))
      time_marching%c = (/ 0.d0 /)

    CASE('implicit_euler','ie')
      time_marching%implicit = 1
      time_marching%stifflyAcc = .TRUE.
      time_marching%p = 1;    time_marching%q = 1
      time_marching%r = 1;    time_marching%s = 1;    time_marching%s_initial = 1

      ALLOCATE(time_marching%A(1,1),time_marching%U(1,1),time_marching%B(1,1),&
        time_marching%V(1,1),time_marching%c(1),time_marching%q_step_indx(1))

      time_marching%A = transpose(reshape((/ 1.d0 /),(/ time_marching%s,time_marching%s /)))
      time_marching%B = transpose(reshape((/ 1.d0 /),(/ time_marching%s,time_marching%r /)))
      time_marching%U = transpose(reshape((/ 1.d0 /),(/ time_marching%r,time_marching%s /)))
      time_marching%V = transpose(reshape((/ 1.d0 /),(/ time_marching%r,time_marching%r /)))
      time_marching%c = (/ 1.d0 /)

    CASE('bdf2')
      time_marching%implicit = 1
      time_marching%stifflyAcc = .TRUE.
      time_marching%p = 2;  time_marching%q = 2
      time_marching%r = 2;  time_marching%s = 1;  time_marching%s_initial = 1

      !call AllocateUnsteady(unsteady)

      ALLOCATE(time_marching%A(1,1),time_marching%U(1,2),time_marching%B(2,1),&
        time_marching%V(2,2),time_marching%c(1),time_marching%q_step_indx(2))

      time_marching% A(1,1) = 2.d0/3.d0
      time_marching% B = transpose(reshape((/ &
           2.d0/3.d0, &
           0.d0 &
           /),(/ time_marching% s, time_marching% r /)))
      time_marching% U = transpose(reshape((/ &
           4.d0/3.d0, -1.d0/3.d0 &
           /),(/ time_marching% r, time_marching% s /)))
      time_marching% V = transpose(reshape((/ &
           4.d0/3.d0, -1.d0/3.d0, &
           1.d0, 0.d0 &
           /),(/ time_marching% r, time_marching% r /)))
           time_marching%c = (/ 1.d0 /)

       CASE('sdirk3')
          time_marching%implicit = 1
          time_marching%stifflyAcc = .TRUE.
          time_marching%p = 3;  time_marching%q = 1
          time_marching%r = 1;  time_marching%s = 4;  time_marching%s_initial = 1

          !call AllocateUnsteady(unsteady)

          ALLOCATE(time_marching%A(4,4),time_marching%U(4,1),time_marching%B(1,4),&
            time_marching%V(1,1),time_marching%c(4),time_marching%q_step_indx(1))

          time_marching%A = transpose(reshape( (/ &
               5.7281606245015908d-1,0.0000000000000000d0,0.0000000000000000d0,0.0000000000000000d0, &
               -4.2974075486898089d-2,5.7281606245015908d-1,0.0000000000000000d0,0.0000000000000000d0, &
               -3.8405397546217372d0,3.6079433958911875d0,5.7281606245015908d-1,0.0000000000000000d0, &
               -9.9166439219916533d0,1.1013593846890206d1,-6.6976598734871517d-1,5.7281606245015908d-1 &
               /),(/ time_marching%s, time_marching%s /)))
          time_marching%B = transpose(reshape( (/ &
               -9.9166439219916533d0,1.1013593846890206d1,-6.6976598734871517d-1,5.7281606245015908d-1 &
               /),(/ time_marching%s, time_marching%r /)))
          time_marching%U = transpose(reshape( (/ &
               1.0000000000000000d0, &
               1.0000000000000000d0, &
               1.0000000000000000d0, &
               1.0000000000000000d0 &
               /),(/ time_marching%r, time_marching%s /)))
          time_marching%V = transpose(reshape( (/ &
               1.0000000000000000d0 &
               /),(/ time_marching%r, time_marching%r /)))
           time_marching% c = SUM(time_marching% A,DIM=2)

       CASE('sdirk4')
          time_marching%implicit = 1
          time_marching%stifflyAcc = .FALSE.
          time_marching%p = 4;  time_marching%q = 1
          time_marching%r = 1;  time_marching%s = 3;  time_marching%s_initial = 1

          !call AllocateUnsteady(unsteady)

          ALLOCATE(time_marching%A(3,3),time_marching%U(3,1),time_marching%B(1,3),&
            time_marching%V(1,1),time_marching%c(3),time_marching%q_step_indx(1))

          time_marching%A = transpose(reshape( (/ &
               1.0685790213016289d0,0.0000000000000000d0,0.0000000000000000d0, &
               -5.6857902130162898d-1,1.0685790213016289d0,0.0000000000000000d0, &
               2.1371580426032577d0,-3.2743160852065154d0,1.0685790213016289d0 &
               /),(/ time_marching%s, time_marching%s /)))
          time_marching%B = transpose(reshape( (/ &
               1.2888640051572038d-1,7.4222719896855924d-1,1.2888640051572038d-1 &
               /),(/ time_marching%s, time_marching%r /)))
          time_marching%U = transpose(reshape( (/ &
               1.0000000000000000d0, &
               1.0000000000000000d0, &
               1.0000000000000000d0 &
               /),(/ time_marching%r, time_marching%s /)))
          time_marching%V = transpose(reshape( (/ &
               1.0000000000000000d0 &
               /),(/ time_marching%r, time_marching%r /)))
           time_marching% c = SUM(time_marching% A,DIM=2)

       CASE('trap')
         time_marching%implicit = 1
         time_marching%stifflyAcc = .TRUE.
         time_marching%p = 2;  time_marching%q = 2
         time_marching%r = 1;  time_marching%s = 2;  time_marching%s_initial = 2

         !call AllocateUnsteady(unsteady)

         ALLOCATE(time_marching%A(2,2),time_marching%U(2,1),time_marching%B(1,2),&
           time_marching%V(1,1),time_marching%c(2),time_marching%q_step_indx(1))

         time_marching%A = transpose(reshape( (/ &
              0.0000000000000000d0,0.0000000000000000d0, &
              5.0000000000000000d-1,5.0000000000000000d-1 &
              /),(/ time_marching%s, time_marching%s /)))
         time_marching%B = transpose(reshape( (/ &
              5.0000000000000000d-1,5.0000000000000000d-1 &
              /),(/ time_marching%s, time_marching%r /)))
         time_marching%U = transpose(reshape( (/ &
              1.0000000000000000d0, &
              1.0000000000000000d0 &
              /),(/ time_marching%r, time_marching%s /)))
         time_marching%V = transpose(reshape( (/ &
              1.0000000000000000d0 &
              /),(/ time_marching%r, time_marching%r /)))
         time_marching% c = SUM(time_marching% A,DIM=2)

       CASE('esdirk4')
         time_marching%implicit = 1
         time_marching%stifflyAcc = .TRUE.
         time_marching%p = 4;  time_marching%q = 2
         time_marching%r = 1;  time_marching%s = 6;  time_marching%s_initial = 2

         !call AllocateUnsteady(unsteady)

         ALLOCATE(time_marching%A(6,6),time_marching%U(6,1),time_marching%B(1,6),&
           time_marching%V(1,1),time_marching%c(6),time_marching%q_step_indx(1))

         time_marching% A = transpose(reshape( (/ &
              0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
              1.d0/4.d0, 1.d0/4.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
              8611.d0/62500.d0, -1743.d0/31250.d0, 1.d0/4.d0, 0.d0, 0.d0, 0.d0, &
              5012029.d0/34652500.d0, -654441.d0/2922500.d0, 174375.d0/388108.d0, &
                1.d0/4.d0, 0.d0, 0.d0, &
              15267082809.d0/155376265600.d0, -71443401.d0/120774400.d0, &
                730878875.d0/902184768.d0, 2285395.d0/8070912.d0, 1.d0/4.d0, 0.d0, &
              82889.d0/524892.d0, 0.d0, 15625.d0/83664.d0, 69875.d0/102672.d0, &
                -2260.d0/8211.d0, 1.d0/4.d0 &
              /),(/ time_marching% s, time_marching% s /)))
         time_marching% B = transpose(reshape((/ &
              82889.d0/524892.d0, 0.d0, 15625.d0/83664.d0, 69875.d0/102672.d0, &
                -2260.d0/8211.d0, 1.d0/4.d0 &
              /),(/ time_marching% s, time_marching% r /)))
         time_marching% U = transpose(reshape((/ &
              1.d0, &
              1.d0, &
              1.d0, &
              1.d0, &
              1.d0, &
              1.d0 &
              /),(/ time_marching% r, time_marching% s /)))
         time_marching% V(1,1) = 1.d0
         time_marching% c = SUM(time_marching% A,DIM=2)

    CASE DEFAULT
      CALL end_mpi()
    END SELECT

    DO i=1,time_marching%r;   time_marching%q_step_indx(i) = i;    END DO

    IF (time_marching%implicit == 0) THEN;  time_marching%A(i,:) = dt*time_marching%A(i,:)
    ELSEIF (time_marching%implicit == 1) THEN
      ALLOCATE(time_marching%LHS_tvec(time_marching%s))
      DO i = time_marching%s_initial,time_marching%s
        time_marching%LHS_tvec(i) = 1.d0/(time_marching%A(i,i) * dt)
        time_marching%U(i,:) = time_marching%U(i,:) * time_marching%LHS_tvec(i)
        time_marching%A(i,:) = time_marching%A(i,:) / time_marching%A(i,i)
      END DO
      ALLOCATE(time_marching%dc(time_marching%s+1))
      time_marching%dc(1) = time_marching%c(1)
      DO i = 2,time_marching%s
        time_marching%dc(i) = time_marching%c(i) - time_marching%c(i-1)
      END DO
      time_marching%dc(i) = 1.d0 - time_marching%c(i-1)
    END IF
    time_marching%B = dt*time_marching%B

    RETURN

  END SUBROUTINE
! initialise_time_marching
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! time_marching_mod
!===============================================================================
!
