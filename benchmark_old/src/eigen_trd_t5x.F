!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t5x_mod
!
! Purpose
! =======
!
! eigen_trd_t5_mod manages the modules for the special update to the
! next target vector for creating a Householder reflector.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_trd_t5x_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_local_2update0

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_local_2update0
!
! Purpose
! =======
!
! The special updater to the next target vector for creating a
! Householder reflector.
!
!
! Arguments
! =========
!
! uxx    (input/output) real(8) array, dimension(*)
!        uxx returns the next target reflector vector updated.
!
! ux     (input) real(8) array, dimension(*)
!        ux contains the reflector vector.
!
! vx     (input) real(8) array, dimension(*)
!        vx contains the counterpart vector to ux.
!
! m      (input) integer
!        Block factor
!
! i      (input) integer
!        The global loop counter
!
! u_n    (input) real(8)
!        u_n contains replicated elements of ux(i,i:i+m-1)
!
! v_n    (input) real(8)
!        v_n contains replicated elements of vx(i,i:i+m-1)
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_local_2update0(
     &     uxx, ux, vx,
     &     m, i, u_n, v_n)

      real(8), intent(inout) ::  uxx(*)
      real(8), intent(in)    ::  ux(*)
      real(8), intent(in)    ::  vx(*)
      integer, intent(in)    ::  m
      integer, intent(in)    ::  i
      real(8), intent(in)    ::  u_n
      real(8), intent(in)    ::  v_n

      integer                ::  j, k, L
      integer                ::  j_1, j_2, j_3

      include 'CSTAB.h'


      if (m < 1) return
     
      L = i - 1

      j_2 = eigen_loop_start(1, 'X')
      j_3 = eigen_loop_end  (L, 'X')

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
      do j_1=j_2,j_3
        uxx(j_1) = uxx(j_1)
     &           - (ux(j_1)*v_n + vx(j_1)*u_n)
      end do                    ! j_1

      return

      end subroutine eigen_trd_local_2update0

      end module eigen_trd_t5x_mod
