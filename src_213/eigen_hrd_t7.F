      subroutine eigen_hrd_load(a, w, nm,
     &              d_t,
     &              u_x, u_y, v_x, v_y, nv,
     &              m0, i_base, i_block)
!$    use omp_lib
      use eigen_devel_mod
      use eigen_libs_mod
      implicit NONE

      integer, intent(in)    ::  nm, nv, m0, i_base, i_block
      complex(8), intent(inout) ::  a(1:nm,*)
      complex(8), intent(out)   ::  w(1:nm,*)
      complex(8), intent(out)   ::  u_x(1:nv,*), u_y(1:nv,*)
      complex(8), intent(out)   ::  v_x(1:nv,*), v_y(1:nv,*)
      real(8), intent(out)   ::  d_t(*)

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j, L
      real(8) :: s


      i_2 = eigen_loop_start(i_base+1,  'Y')
      i_3 = eigen_loop_end  (i_base+m0, 'Y')
      j_2 = eigen_loop_start(1,         'X')
      j_3 = eigen_loop_end  (i_base+m0, 'X')

!$    if ( omp_get_thread_num() == 0 ) then

      do i_1=i_2,i_3
         j = eigen_translate_l2g(i_1, 'Y')
         do j_1=j_2,j_3
            w(j_1,j-i_base) = a(j_1,i_1)
         end do! j_1
      end do! i_1

      if ( x_nnod /= y_nnod ) then
         if ( diag_0 > 0 ) then
            call eigen_diag_loop_info( 1, i_base+m0-1,
     &                 k_2, k_3, i_2, i_4, j_2, j_4 )
            if (k_2 <= k_3) then
               do  k_1 = k_2, k_3
                  i_1 = i_2 + k_1 * i_4
                  j_1 = j_2 + k_1 * j_4
                  d_t(i_1)   = DBLE(a(j_1,i_1))
                  a(j_1,i_1) = ZEROZ
               end do! k_1
            end if
         end if
      else
         if ( x_inod == y_inod ) then
            i_2 = eigen_loop_start(1,           'Y')
            i_3 = eigen_loop_end  (i_base+m0-1, 'Y')
            if (i_2 <= i_3) then
               do i_1=i_2,i_3
                  d_t(i_1)   = DBLE(a(i_1,i_1))
                  a(i_1,i_1) = ZEROZ
               end do! i_1
            end if
         end if
      end if

!$    end if

!$    if ( omp_get_num_threads() == 1 .OR.
!$   $     omp_get_thread_num() == 1 ) then

      i = i_base+m0
      L = i - 2
      y_pos  = eigen_translate_g2l(i, 'Y')
      x_pos  = eigen_translate_g2l(i, 'X')

      k_2 = m0
      k_3 = MAX(1, 3*(2-i_block))
      do k_1=1,k_2
         do j_1=1,x_pos
            u_x(j_1,k_1) = ZEROZ
            v_x(j_1,k_1) = ZEROZ
         end do! j_1
         do j_1=1,y_pos
            u_y(j_1,k_1) = ZEROZ
            v_y(j_1,k_1) = ZEROZ
         end do! j_1
      end do! k_1

!$    end if

      return
      end subroutine eigen_hrd_load


      subroutine eigen_hrd_store(a, w, nm,
     &              d_t,
     &              m0, i_base)
!$    use OMP_LIB
      use eigen_devel_mod
      use eigen_libs_mod
      implicit NONE

      integer, intent(in)    ::  nm, m0, i_base
      complex(8), intent(out)   ::  a(1:nm,*)
      complex(8), intent(in)    ::  w(1:nm,*)
      real(8), intent(in)    ::  d_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j, L


      i_2 = eigen_loop_start(i_base+1,  'Y')
      i_3 = eigen_loop_end  (i_base+m0, 'Y')
      do i_1=i_2,i_3
         j   = eigen_translate_l2g(i_1, 'Y')
         j_3 = eigen_loop_end  (j, 'X')
         do j_1=1,j_3
            a(j_1,i_1) = w(j_1,j-i_base)
         end do! j_1
      end do! i_1

      if ( x_nnod /= y_nnod ) then
         if ( diag_0 > 0 ) then
            call eigen_diag_loop_info( 1, i_base,
     &                 k_2, k_3, i_2, i_4, j_2, j_4 )
            if (k_2 <= k_3) then
               do  k_1 = k_2, k_3
                  i_1 = i_2 + k_1 * i_4
                  j_1 = j_2 + k_1 * j_4
                  a(j_1, i_1) = DCMPLX(d_t(i_1),ZERO)
               end do! k_1
            end if
         end if
      else
         if ( x_inod == y_inod ) then
            i_2 = eigen_loop_start(1,      'Y')
            i_3 = eigen_loop_end  (i_base, 'Y')
            if ( i_2 <= i_3 ) then
               do i_1=i_2,i_3
                  a(i_1,i_1) = DCMPLX(d_t(i_1),ZERO)
               end do! i_1
            end if
         end if
      end if

      return
      end subroutine eigen_hrd_store

