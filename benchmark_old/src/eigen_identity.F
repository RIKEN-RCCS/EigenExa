      subroutine eigen_identity( n, x, ldx )

      use eigen_libs_mod
      use eigen_devel_mod
      implicit none

      integer, intent(in)    :: n
      real(8), intent(inout) :: x(ldx,*)
      integer, intent(in)    :: ldx

      integer :: i, i_1, i_2, i_3, i_4
      integer :: j, j_1, j_2, j_3, j_4
      integer :: k, k_1, k_2, k_3, k_4

      call eigen_loop_info( 1, n, i_2, i_3, 'Y' )
      call eigen_loop_info( 1, n, j_2, j_3, 'X' )

#if 0
!$OMP PARALLEL DO PRIVATE(i,i_1,j_1)
      do i_1 = i_2, i_3
        x(1:ldx, i_1) = ZERO
        i   = eigen_translate_l2g(i_1, 'Y')
        j_1 = eigen_owner_index(i, 'X')
        if ( j_1 >= 1 ) then
          x(j_1, i_1) = ONE
        end if
      end do
!$OMP END PARALLEL DO
#else
      if ( x_nnod /= y_nnod ) then
!$OMP PARALLEL DO
        do i_1 = i_2, i_3
          x(1:ldx, i_1) = ZERO
        end do
!$OMP END PARALLEL DO
        if ( diag_0 > 0 ) then
          call eigen_diag_loop_info( 1, n,
     &             k_2, k_3, i_2, i_4, j_2, j_4 )
          if ( k_2 <= k_3 ) then
!$OMP PARALLEL DO PRIVATE(i_1,j_1)
            do k_1 = k_2, k_3
              i_1 = i_2 + k_1 * i_4
              j_1 = j_2 + k_1 * j_4
              x(j_1, i_1) = ONE
            end do
!$OMP END PARALLEL DO
          end if
        end if
      else
!$OMP DO
        do i_1 = i_2, i_3
          x(1:ldx, i_1) = ZERO
          if ( x_inod == y_inod ) then
            x(i_1, i_1) = ONE
          end if
        end do
!$OMP END DO
      end if
#endif

      return
      end subroutine eigen_identity

