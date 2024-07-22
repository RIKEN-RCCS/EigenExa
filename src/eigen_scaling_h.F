      subroutine eigen_scaling_h( n, a, lda, x )
!$    use OMP_LIB
      use MPI
      use eigen_libs_mod
      use eigen_devel_mod
      use, intrinsic :: ieee_arithmetic
      implicit NONE

      integer, intent(IN)    :: lda, n
!hermite
!      real(8), intent(INOUT) :: a(lda, *)
      complex(8), intent(INOUT) :: a(lda, *)
      real(8), intent(OUT)   :: x

      real(8)                :: SAFMIN, EPS, SMLNUM, BIGNUM
      real(8)                :: RMIN, RMAX, ANRM, SIGMA
      complex(8)             :: t
      real(8)                :: DNAN, tt(2), ss(2)
      real(8), external      :: DLAMCH

      integer                :: ERR_MPI
      integer                :: i, i_1, i_2, i_3
      integer                :: j, j_1, j_2, j_3


         SAFMIN = DLAMCH( 'Safe minimum' )
         EPS    = DLAMCH( 'Precision' )
         SMLNUM = SAFMIN / EPS
         BIGNUM = ONE / SMLNUM
         RMIN = SQRT( SMLNUM )
         RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*-
         i_2 = eigen_loop_start   ( 1, 'Y' )
         i_3 = eigen_loop_end     ( n, 'Y' )

         ANRM = ZERO
         DNAN = ZERO
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1,j_2,j_3)
!$OMP+         PRIVATE(t)
!$OMP+         REDUCTION(max: ANRM, DNAN)
         do i_1=i_2,i_3

            i   = eigen_translate_l2g( i_1, 'Y' )
            j   = i

            j_2 = eigen_loop_start   ( 1, 'X' )
            j_3 = eigen_loop_end     ( j, 'X' )

            do j_1=j_2,j_3
               t = a(j_1, i_1)
               if ( ieee_is_finite( DBLE(t) ) .and.
     &              ieee_is_finite( DIMAG(t) ) ) then
                  ANRM = MAX( ANRM, ABS( t ) )
               else
                  DNAN = ONE
               end if
            end do

         end do
!$OMP END PARALLEL DO

         tt(1) = ANRM
         tt(2) = DNAN
         call MPI_Allreduce( tt, ss,
     &                       2, MPI_DOUBLE_PRECISION,
     &                       MPI_MAX, TRD_COMM_WORLD, ERR_MPI )
         ANRM = ss(1)
         DNAN = ss(2)
         if ( DNAN /= ZERO ) then
            x = get_constant_nan()
            return
         end if
*-
!
!   SIGMA is calcurated in order to avoid Under/over-flow
!
         SIGMA = ONE
         IF( ANRM /= ZERO .AND. ANRM < RMIN ) THEN
            SIGMA = RMIN / ANRM
         ELSE IF( ANRM > RMAX ) THEN
            SIGMA = RMAX / ANRM
         END IF
*-
         x = SIGMA
         if ( x == ONE ) return
*-
         i_2 = eigen_loop_start   ( 1, 'Y' )
         i_3 = eigen_loop_end     ( n, 'Y' )

!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1)
         do i_1=i_2,i_3

            i = eigen_translate_l2g( i_1, 'Y' )
            j = eigen_loop_end     ( i,   'X' )

            call zdscal(j, x, a(1, i_1), 1)

         end do
!$OMP END PARALLEL DO

      return
      end subroutine  eigen_scaling_h

