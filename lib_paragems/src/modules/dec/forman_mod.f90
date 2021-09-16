!
!===============================================================================
!-- forman Module
!> Module containing routines for performing forman operations
!===============================================================================
!/****/h* modules|forman/forman_mod
!* SYNOPSIS
MODULE forman_mod
!* PURPOSE
!*   Module containing routines for performing forman operations
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
!> Module containing routines for performing forman operations
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

END MODULE
! forman_mod
!===============================================================================
!
