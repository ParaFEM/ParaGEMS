!
!===============================================================================
!-- Elasticity Module
!> Module contains routines specifically related to elasticity equations
!===============================================================================
!/****/h* modules|phy_elasticity/elasticity_mod
MODULE elasticity_mod
!* PURPOSE
!*   Module contains routines specifically related to elasticity equations
!* INCLUDES
!*   Name                    Purpose
!*   common_mod              variable definitions
!* CONTAINS
!*   Subroutine              Purpose
!*   deemat_linelasticity    create material matrix for linear elasticity
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2022/04/30: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2022/04/30
!> Module contains routines specifically related to elasticity equations
!===============================================================================

  !-- use statements and implicit none --
  USE common_mod

  IMPLICIT NONE

!===============================================================================
CONTAINS
!===============================================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!/****/s* elasticity_mod/deemat_linelasticity
!* SYNOPSIS
  SUBROUTINE deemat_linelasticity(dee,e,v)
!* PURPOSE
!*   Compute the material matrix (D-matrix) for a linear elastic material.
!* AUTHOR
!*   Pieter Boom (adapted from parafem deemat subroutine)
!* MODIFICATION HISTORY
!*   2020/09/16: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom (adapted from parafem deemat subroutine)
!> date: 2020/09/16
!> Compute the material matrix (D-matrix) for a linear elastic material.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !-- arguments and local variables --
  REAL(KIND=PGMSiwp),INTENT(IN)    :: e,v
  REAL(KIND=PGMSiwp),INTENT(INOUT) :: dee(:,:)
  REAL(KIND=PGMSiwp)               :: v2,vv
  REAL(KIND=PGMSiwp),PARAMETER     :: zero=0.d0,one=1.d0,two=2.d0
  INTEGER                 :: i

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAIN EXECUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    dee = zero; vv = v*e/(one-2*v)/(one+v); v2 = e/(two+two*v)
    DO i=2,4; dee(i,i)=v2; END DO
    DO i=6,8; dee(i,i)=v2; END DO
    DO i=1,9,4
      dee(i,(/ 1,5,9 /)) = vv
      dee(i,i) = dee(i,i) + 2*v2
    END DO
    dee(2,4) = v2; dee(3,7) = v2; dee(4,2) = v2
    dee(6,8) = v2; dee(7,3) = v2; dee(8,6) = v2

  END SUBROUTINE
! elasticity_mod/deemat_linelasticity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

END MODULE
! elasticity_mod
!===============================================================================
!
