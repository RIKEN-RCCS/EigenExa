!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_house_mod
!
! Purpose
! =======
!
! eigen_house_mod manages the modules for householder transformation.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_house_mod

      implicit none

      ! Since not sure of allocatables in module,
      ! thus, replaced into pointer attribute.
      real(8), pointer       ::  u0_z(:), v0_z(:)
      real(8), pointer       ::  u1_z(:), v1_z(:)
      integer                ::  offset1, offset2
      integer                ::  offset3, offset4

      integer, parameter     ::  LOOP_INFO_X_START = 1
      integer, parameter     ::  LOOP_INFO_X_END   = 2
      integer, parameter     ::  LOOP_INFO_X_OINDX = 3
      integer, parameter     ::  LOOP_INFO_X_ROOT  = 4
      integer, parameter     ::  LOOP_INFO_X_POS   = 5

      integer, parameter     ::  LOOP_INFO_Y_START = 6
      integer, parameter     ::  LOOP_INFO_Y_END   = 7
      integer, parameter     ::  LOOP_INFO_Y_OINDX = 8
      integer, parameter     ::  LOOP_INFO_Y_ROOT  = 9
      integer, parameter     ::  LOOP_INFO_Y_POS   = 10

      integer, parameter     ::  LOOP_INFO_V_WIDTH = 11
      integer, parameter     ::  LOOP_INFO_V_CURNT = 12
      integer, parameter     ::  LOOP_INFO_I_BASE  = 13

      integer, parameter     ::  LOOP_INFO_L_RANK  = 15
      integer, parameter     ::  LOOP_INFO_L_SIZE  = 16

      end module eigen_house_mod
