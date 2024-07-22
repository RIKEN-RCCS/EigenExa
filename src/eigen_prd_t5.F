!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_prd_t5_mod
!
! Purpose
! =======
!
! eigen_prd_t5_mod manages the modules for the local update.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_prd_t5_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_prd_local_2update

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_local_2update
!
! Purpose
! =======
!
! The local updater.
!
!
! Arguments
! =========
!
! w       (input/output) real(8) array, dimension(nm,*)
!         a contains the elements of the target panel.
!         On exit, it returns the panel updated.
!
! nm      (input) integer 
!         The leading dimension of the array w.
!
! ux      (input) real(8) array, dimension(nv,*)
!         ux contanis the reflector vectors.
!
! uy      (input) real(8) array, dimension(nv,*)
!         uy contanis the reflector vectors.
!
! vx      (input) real(8) array, dimension(nv,*)
!         vx contanis the counterpart vectors.
!
! vy      (input) real(8) array, dimension(nv,*)
!         vy contanis the counterpart vectors.
!
! nv      (input) integer 
!         The leading dimension of the working arrays, u_x, u_y, ...
!
! i_base  (input) integer 
!         The offset pointer
!
! i       (input) integer 
!         The global loop counter
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_prd_local_2update(
     &     w, nm,
     &     ux, uy, vx, vy, nv,
     &     i_base, i)

      real(8), intent(inout) ::  w(1:nm,*)
      integer, intent(in)    ::  nm
      real(8), intent(in)    ::  ux(1:nv,*)
      real(8), intent(in)    ::  uy(1:nv,*)
      real(8), intent(in)    ::  vx(1:nv,*)
      real(8), intent(in)    ::  vy(1:nv,*)
      integer, intent(in)    ::  nv
      integer, intent(in)    ::  i_base
      integer, intent(in)    ::  i

      include 'CSTAB.h'

      integer                :: iloop_sta, iloop_end
      integer                :: jloop_sta, jloop_end
      integer                :: i_1, i_4
      integer                :: j, j_1, j_2, j_3
      integer                :: l, l_1
      integer                :: jj_1
      integer                :: k_1, kk_1, kk_2, kk_3
      integer                :: LX

      real(8)                :: u_x0, v_x0
      real(8)                :: u_x1, v_x1
      real(8)                :: uy0_0, vy0_0
      real(8)                :: uy1_0, vy1_0
      real(8)                :: uy0_1, vy0_1
      real(8)                :: uy1_1, vy1_1
      real(8)                :: uy0_2, vy0_2
      real(8)                :: uy1_2, vy1_2
      real(8)                :: uy0_3, vy0_3
      real(8)                :: uy1_3, vy1_3
      real(8)                :: uy0_4, vy0_4
      real(8)                :: uy1_4, vy1_4
      real(8)                :: uy0_5, vy0_5
      real(8)                :: uy1_5, vy1_5
      real(8)                :: w0
      real(8)                :: w1
      real(8)                :: w2
      real(8)                :: w3
      real(8)                :: w4
      real(8)                :: w5


      k_1 = i - i_base
      if (k_1 <= 1) return

      LX = 64

      iloop_sta = eigen_loop_start(i_base+1, y_nnod,y_inod)
      iloop_end = eigen_loop_end  (i-MBAND,  y_nnod,y_inod)
      i_4 = mod(iloop_end-iloop_sta+1,6)+iloop_sta
      if (iloop_sta > iloop_end) return

      L = i - MBAND
      jloop_sta = eigen_loop_start(1, x_nnod,x_inod)
      jloop_end = eigen_loop_end  (L, x_nnod,x_inod)

      do jj_1=jloop_sta,jloop_end,LX
        j_2 = jj_1; j_3 = min(jj_1+LX-1,jloop_end)

        do i_1=iloop_sta,i_4-1  ! 0
          j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
          l_1 = j-i_base
          uy0_0 = uy(i_1+0,k_1-0)
          vy0_0 = vy(i_1+0,k_1-0)
          uy1_0 = uy(i_1+0,k_1-1)
          vy1_0 = vy(i_1+0,k_1-1)

!DIR$ IVDEP
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
          do j_1=j_2,j_3
            u_x0 = ux(j_1,k_1-0)
            v_x0 = vx(j_1,k_1-0)
            u_x1 = ux(j_1,k_1-1)
            v_x1 = vx(j_1,k_1-1)

            w0 = w(j_1,l_1+0*y_nnod)
            w0 = w0
     &           - (u_x0 * vy0_0)
     &           - (v_x0 * uy0_0)
            w0 = w0
     &           - (u_x1 * vy1_0)
     &           - (v_x1 * uy1_0)
            w(j_1,l_1+0*y_nnod) = w0
          end do                ! j_1

        end do                  ! l_1

      end do                  ! l_1

      do jj_1=jloop_sta,jloop_end,LX
        j_2 = jj_1; j_3 = min(jj_1+LX-1,jloop_end)

        do i_1=i_4,iloop_end,6  ! 5
          j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
          l_1 = j-i_base
          uy0_0 = uy(i_1+0,k_1-0)
          vy0_0 = vy(i_1+0,k_1-0)
          uy1_0 = uy(i_1+0,k_1-1)
          vy1_0 = vy(i_1+0,k_1-1)
          uy0_1 = uy(i_1+1,k_1-0)
          vy0_1 = vy(i_1+1,k_1-0)
          uy1_1 = uy(i_1+1,k_1-1)
          vy1_1 = vy(i_1+1,k_1-1)
          uy0_2 = uy(i_1+2,k_1-0)
          vy0_2 = vy(i_1+2,k_1-0)
          uy1_2 = uy(i_1+2,k_1-1)
          vy1_2 = vy(i_1+2,k_1-1)
          uy0_3 = uy(i_1+3,k_1-0)
          vy0_3 = vy(i_1+3,k_1-0)
          uy1_3 = uy(i_1+3,k_1-1)
          vy1_3 = vy(i_1+3,k_1-1)
          uy0_4 = uy(i_1+4,k_1-0)
          vy0_4 = vy(i_1+4,k_1-0)
          uy1_4 = uy(i_1+4,k_1-1)
          vy1_4 = vy(i_1+4,k_1-1)
          uy0_5 = uy(i_1+5,k_1-0)
          vy0_5 = vy(i_1+5,k_1-0)
          uy1_5 = uy(i_1+5,k_1-1)
          vy1_5 = vy(i_1+5,k_1-1)

          kk_2=j_2; kk_3=j_3

!DIR$ IVDEP
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
          do j_1=kk_2,kk_3
            u_x0 = ux(j_1,k_1-0)
            v_x0 = vx(j_1,k_1-0)
            u_x1 = ux(j_1,k_1-1)
            v_x1 = vx(j_1,k_1-1)

            w0 = w(j_1,l_1+0*y_nnod)
            w1 = w(j_1,l_1+1*y_nnod)
            w0 = w0
     &           - (u_x0 * vy0_0)
     &           - (v_x0 * uy0_0)
            w1 = w1
     &           - (u_x0 * vy0_1)
     &           - (v_x0 * uy0_1)
            w0 = w0
     &           - (u_x1 * vy1_0)
     &           - (v_x1 * uy1_0)
            w1 = w1
     &           - (u_x1 * vy1_1)
     &           - (v_x1 * uy1_1)
            w(j_1,l_1+0*y_nnod) = w0
            w(j_1,l_1+1*y_nnod) = w1

            w2 = w(j_1,l_1+2*y_nnod)
            w3 = w(j_1,l_1+3*y_nnod)
            w2 = w2
     &           - (u_x0 * vy0_2)
     &           - (v_x0 * uy0_2)
            w3 = w3
     &           - (u_x0 * vy0_3)
     &           - (v_x0 * uy0_3)
            w2 = w2
     &           - (u_x1 * vy1_2)
     &           - (v_x1 * uy1_2)
            w3 = w3
     &           - (u_x1 * vy1_3)
     &           - (v_x1 * uy1_3)
            w(j_1,l_1+2*y_nnod) = w2
            w(j_1,l_1+3*y_nnod) = w3


            w4 = w(j_1,l_1+4*y_nnod)
            w5 = w(j_1,l_1+5*y_nnod)
            w4 = w4
     &           - (u_x0 * vy0_4)
     &           - (v_x0 * uy0_4)
            w5 = w5
     &           - (u_x0 * vy0_5)
     &           - (v_x0 * uy0_5)
            w4 = w4
     &           - (u_x1 * vy1_4)
     &           - (v_x1 * uy1_4)
            w5 = w5
     &           - (u_x1 * vy1_5)
     &           - (v_x1 * uy1_5)
            w(j_1,l_1+4*y_nnod) = w4
            w(j_1,l_1+5*y_nnod) = w5
          end do              ! j_1

        end do                  ! l_1

      end do                    ! jj_1

      return

      end subroutine eigen_prd_local_2update

      end module eigen_prd_t5_mod
