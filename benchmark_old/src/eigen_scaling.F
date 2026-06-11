!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_scaling_mod
!
! Purpose
! =======
!
! eigen_scaling_mod manages the modules for scaling of the input
! matrix.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_scaling_mod

      use eigen_libs_mod
      use eigen_devel_mod
      use mpi
!$    use omp_lib
      use, intrinsic :: ieee_arithmetic

      implicit none
      private


      public :: eigen_scaling


      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_scaling
!
! Purpose
! =======
!
! eigen_scaling scales the input matrix in order to avoid overflow and
! underflow numerically.
!
!
! Arguments
! =========
!
! n       (input) integer
!         The dimension of the matrix. N >= 0.
!
! a       (input/output) real(8) array, dimension(lda,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! lda     (input) integer
!         The leading dimension of the array a.
!
! x       (output) real(8)
!         The scaling factor
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_scaling(n, a, lda, x)

      integer, intent(in)    :: n
      real(8), intent(inout) :: a(lda, *)
      integer, intent(in)    :: lda
      real(8), intent(inout) :: x

      real(8)                :: SAFMIN, EPS, SMLNUM, BIGNUM
      real(8)                :: RMIN, RMAX, ANRM, SIGMA
      real(8)                :: DNAN, t, tt(2), ss(2)
      real(8), external      :: DLAMCH

      integer                :: ERR_MPI
      integer                :: i, i_1, iloop_sta, iloop_end
      integer                :: j, j_1, jloop_sta, jloop_end


      SAFMIN = DLAMCH('Safe minimum')
      EPS    = DLAMCH('Precision')
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = sqrt(SMLNUM)
      RMAX   = min(sqrt(BIGNUM), ONE / sqrt(sqrt(SAFMIN)))

      iloop_sta = eigen_loop_start(1, 'Y')
      iloop_end = eigen_loop_end  (n, 'Y')

      ANRM = ZERO
      DNAN = ZERO
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1,jloop_sta,jloop_end)
!$OMP+         PRIVATE(t)
!$OMP+         REDUCTION(max: ANRM, DNAN)
      do i_1=iloop_sta,iloop_end

        i   = eigen_translate_l2g(i_1, 'Y')
        j   = i

        jloop_sta = eigen_loop_start(1, 'X')
        jloop_end = eigen_loop_end  (j, 'X')

        do j_1=jloop_sta,jloop_end
          t = a(j_1, i_1)
          if ( ieee_is_finite( t ) ) then
            ANRM = max(ANRM, abs( t ))
          else
            DNAN = ONE
          end if
        end do

      end do
!$OMP END PARALLEL DO

      tt(1) = ANRM
      tt(2) = DNAN
      call MPI_Allreduce(tt, ss,
     &     2, MPI_DOUBLE_PRECISION,
     &     MPI_MAX, TRD_COMM_WORLD, ERR_MPI)
      ANRM = ss(1)
      DNAN = ss(2)
      if ( DNAN /= ZERO ) then
         x = get_constant_nan()
         return
      end if
!-
!     
!     SIGMA is calcurated in order to avoid Under/over-flow
!     
      SIGMA = ONE
      if (ANRM /= ZERO .and. ANRM < RMIN) then
         SIGMA = RMIN / ANRM
      else if (ANRM > RMAX) then
         SIGMA = RMAX / ANRM
      end if

      x = SIGMA
      if (x == ONE) return

      iloop_sta = eigen_loop_start(1, 'Y')
      iloop_end = eigen_loop_end  (n, 'Y')

!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1)
      do i_1=iloop_sta,iloop_end

        i = eigen_translate_l2g(i_1, 'Y')
        j = eigen_loop_end  (i, 'X')

        call dscal(j, x, a(1, i_1), 1)

      end do
!$OMP END PARALLEL DO

      return

      end subroutine eigen_scaling

      end module eigen_scaling_mod
