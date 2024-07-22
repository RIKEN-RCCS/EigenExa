!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trbak_mod
!
! Purpose
! =======
!
! eigen_trbak_mod manages the modules of the back-transformation.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_trbak_mod

!$    use omp_lib

      implicit none
      public

      real(8) :: trbk_time_bcast
      real(8) :: trbk_time_reduc
      real(8) :: trbk_time_fr
      real(8) :: trbk_time_trbk1
      real(8) :: trbk_time_trbk1_
      real(8) :: trbk_time_trbk1x
      real(8) :: trbk_time_trbk1x_
      real(8) :: trbk_time_trbk2

      real(8) :: trbk_time_reduc_overhead_x

!$    integer(kind=omp_lock_kind) :: trbk_lock
!$    integer(8)                  :: trbk_mask(1:2)

      integer :: do_overlap_bcast_level
      integer :: trbk_time_counter
      integer :: trbk_time_interval
      integer :: trbk_time_next
      integer :: trbk_switched

!     This must be attributed with pointer for Fortran 200x later
      real(8), pointer :: trbk_buf(:)

      end module eigen_trbak_mod
