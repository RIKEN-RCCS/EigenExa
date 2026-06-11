      subroutine  eigen_hrd_init(a, nm, n,
     &              d_out, e_out,
     &              u_t, v_t, nv)
!$    use omp_lib
      use eigen_devel_mod
      use eigen_libs_mod
      use eigen_house_h_mod
      use CSTAB_mod
      implicit NONE

      integer, intent(in)    ::  nm, n, nv
      complex(8), intent(inout) ::  a(1:nm, *)
      real(8), intent(out)   ::  d_out(*)
      real(8), intent(out)   ::  e_out(*)
      complex(8), intent(in)    ::  u_t(*), v_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j

      integer                ::  local_size, local_rank

      include 'CSTAB.h'


      d_out(1:n) = ZERO
      e_out(1:n) = ZERO

      if ( diag_0 > 0 ) then
         call eigen_diag_loop_info( 1, n,
     &              k_2, k_3, i_2, i_4, j_2, j_4 )
         if ( k_2 <= k_3 ) then
            do  k_1 = k_2, k_3
               i_1 = i_2 + k_1 * i_4
               j_1 = j_2 + k_1 * j_4
               j   = (i_1-1)*y_nnod+y_inod
               d_out(j) = DBLE(a(j_1, i_1))
            end do ! k_1
         end if
      end if

      i_2 = eigen_loop_start(1, 'Y')
      i_3 = eigen_loop_end  (n, 'Y')
      do  i_1 = i_2, i_3
         j   = eigen_translate_l2g(i_1, 'Y')
         j_2 = eigen_loop_start   (j+1, 'X')
         if ( j <= n ) then
            a(j_2:nm, i_1) = ZEROZ
         else
            a(1:nm, i_1) = ZEROZ
         end if
      end do ! i_1

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( local_size > 1 ) then

         allocate(u0_z_h(nv*local_size+n_columns),
     &            v0_z_h(nv*local_size+n_columns))
         u0_z_h(:) = ZEROZ
         v0_z_h(:) = ZEROZ
         call CSTAB_adjust_base(u0_z_h, u_t, offset1)
         call CSTAB_adjust_base(v0_z_h, v_t, offset2)

      end if


      return
      end subroutine  eigen_hrd_init


      subroutine  eigen_hrd_final(a, nm, n, d_out, e_out, u_t)
!$    use omp_lib
      use eigen_devel_mod
      use eigen_libs_mod
      use eigen_house_h_mod
      use comm_mod
      use, intrinsic :: iso_c_binding
      implicit NONE

      integer, intent(in)    ::  nm, n
      complex(8), intent(inout) ::  a(1:nm, *)
      real(8), intent(out)   ::  d_out(n)
      real(8), intent(out)   ::  e_out(n)
      complex(8), intent(out), target ::  u_t(n)

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j, k, L
      complex(8)                ::  t
      complex(8)  :: a_n
      real(8)  :: g_n
      real(8), pointer ::  u_r(:)

      integer                ::  local_size, local_rank

      include 'CSTAB.h'


      if ( n >= 2 ) then
         i = 2; L = i-1
         i_1 = eigen_owner_index(i, 'Y')
         j_1 = eigen_owner_index(L, 'X')
         if ( i_1 > 0 .and. j_1 > 0 ) then
            a_n = a(j_1, i_1)
            g_n   = -SIGN(ABS(a_n), DBLE(a_n))
            e_out(i)    = g_n
            a(j_1, i_1) = a(j_1, i_1) - g_n
         end if
         y_owner_nod = eigen_owner_node (i, 'Y')
         x_owner_nod = eigen_owner_node (L, 'X')
         if ( y_owner_nod == y_inod ) then
         call bcast_dbl( e_out(i:), 1, x_owner_nod, 1, x_COMM_WORLD )
         end if
         call bcast_dbl( e_out(i:), 1, y_owner_nod, 1, y_COMM_WORLD )
      end if

      if ( diag_0 > 0 ) then
        call eigen_diag_loop_info( 1, n,
     &             k_2, k_3, i_2, i_4, j_2, j_4 )
         if ( k_2 <= k_3 ) then
            do  k_1 = k_2, k_3
               i_1 = i_2 + k_1 * i_4
               j_1 = j_2 + k_1 * j_4
               j   = (i_1-1)*y_nnod+y_inod
               t           = DCMPLX(d_out(j), ZERO)
               d_out(j)    = DBLE(a(j_1, i_1))
               a(j_1, i_1) = t
            end do ! k_1
         end if
      end if

      u_t(1) = ZEROZ
      call c_f_pointer(c_loc(u_t),u_r,shape=[n])
      call reduce_dbl(d_out, u_r, n, 1, x_COMM_WORLD)
      call reduce_dbl(d_out, u_r, n, 1, y_COMM_WORLD)
      u_r=>null()

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( local_size > 1 ) then

         deallocate ( u0_z_h )
         deallocate ( v0_z_h )

      end if

      return
      end subroutine eigen_hrd_final

