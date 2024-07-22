!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_dc_mod
!
! Purpose
! =======
!
! eigen_dc_mod manages the modules related for the DC methods.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_dc_mod

      implicit none

      real(8) ::  flops, dgemm_time, dgemm_dummy(1:2)
      real(8) ::  p_time0 = 0d0, p_time2 = 0d0, p_time3 = 0d0
      real(8) ::  p_times = 0d0, p_timez = 0d0, p_timer = 0d0

      end module eigen_dc_mod
