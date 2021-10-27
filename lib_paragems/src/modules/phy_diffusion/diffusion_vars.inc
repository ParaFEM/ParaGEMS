!
!===============================================================================
!/****v* modules|phy_diffusion_flow/diffusion_vars
!* PURPOSE
!*   Define global variables for diffusion equations
!* CONTAINS
!*   Structure/Variable      Purpose
!*   diffusion parameters
!*   diffusion namelist
!* AUTHOR
!*   Pieter Boom
!* MODIFICATION HISTORY
!*   2019/08/21: Created (PB)
!* COPYRIGHT
!*   (c) University of Manchester
!******/
!> author: Pieter Boom
!> date: 2019/08/21
!> Define global variables for diffusion equations
!===============================================================================

  !-----------------------------------------------------------------------------
  ! variables
  !-----------------------------------------------------------------------------
  !-- diffusion parameters --
  REAL(KIND=PGMSiwp) :: Dk = 1.d0                           ! diffusion coefficient
  REAL(KIND=PGMSiwp) :: blk_dir(3) = (/ 1.d0, 0.d0, 0.d0 /) ! bulk direction
  REAL(KIND=PGMSiwp) :: ref_val = 1.d0                      ! reference value
  REAL(KIND=PGMSiwp) :: bc_vals(30) = 0.d0                  ! boundary value (max 30)

  !-- Darcy's Law (fluid) parameters --
  !-- p - pressure; Re*k/mu - Reynolds * permeability / viscosity; q - flow rate

  !-- Ohm's Law (elec) parameters --
  ! q - Charge; R * k - resistance*Coulomb constant; I - Current

  !-- Fourier's Law (heat) parameters --
  ! T - temperature; k - conductivity; q heat flux

  !-- Fick's Law (mass) parameters --
  ! phi - concentrations; D - diffusion coefficient; J - flux

  !-- cracking variables --
  INTEGER        :: crck_type = 0                       ! cracking process
  INTEGER        :: max_crcks = 1                       ! max cracks
  INTEGER        :: crcks_pstep = 1                     ! max faces cracked / step
  REAL(KIND=PGMSiwp) :: crck_thrshld = 0.0d0                ! min vel for cracking
  REAL(KIND=PGMSiwp) :: crck_area = 0.0d0                   ! cumulative crack area
  REAL(KIND=PGMSiwp) :: rnd_max = 0.0d0                     ! face variation

  !-- diffusion namelist --
  NAMELIST /diff_param/ Dk, blk_dir, ref_val, bc_vals
  NAMELIST /crack_param/ crck_type, max_crcks, crcks_pstep, crck_thrshld, rnd_max

! diffusion_vars
!===============================================================================
!
