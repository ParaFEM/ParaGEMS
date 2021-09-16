!
!===============================================================================
!/****v* modules|phy_darcy_flow/darcy_vars
!* PURPOSE
!*   Define global variables for Darcy flow equations
!* CONTAINS
!*   Structure/Variable      Purpose
!*   Darcy parameters
!*   Darcy namelist
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Define global variables for Darcy flow equations
!===============================================================================

  !-----------------------------------------------------------------------------
  ! variables
  !-----------------------------------------------------------------------------
  !-- Darcy parameters --
  REAL(KIND=iwp) :: mu = 1.d0                         ! viscosity
  REAL(KIND=iwp) :: k  = 1.d0                         ! permeability
  REAL(KIND=iwp) :: re = 1.d0                         ! Reynold number
  REAL(KIND=iwp) :: vel(3) = (/ 1.d0, 0.d0, 0.d0 /)   ! velocity
  REAL(KIND=iwp) :: bc_press(30) = 0.d0               ! boundary press (max 30)
  REAL(KIND=iwp) :: p_ref = 1.d0                      ! reference pressure

  !-- cracking variables --
  INTEGER        :: crck_type = 0                     ! cracking process
  INTEGER        :: max_crcks = 1                     ! max cracks
  INTEGER        :: crcks_pstep = 1                   ! max faces cracked / step
  REAL(KIND=iwp) :: crck_thrshld = 0.0d0              ! min vel for cracking
  REAL(KIND=iwp) :: crck_area = 0.0d0                 ! cumulative crack area
  REAL(KIND=iwp) :: rnd_max = 0.0d0                   ! face variation

  !-- Darcy namelist --
  NAMELIST /darcy_param/ mu, k, re, vel, p_ref, bc_press, &
    crck_type, max_crcks, crcks_pstep, crck_thrshld, rnd_max

! darcy_vars
!===============================================================================
!
