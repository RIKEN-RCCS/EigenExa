      subroutine zero_padding( n, nbb, a, lda )
!$    use OMP_LIB
      use eigen_libs_mod
      implicit NONE

      integer, intent(IN)    :: n, nbb, lda
      real(8), intent(INOUT) :: a(*)
      real(8), parameter     :: zero = 0.0D0

      integer                :: nprocs, npcol, nprow
      integer                :: iam, mycol, myrow
      integer                :: i_1, i_2, i_3
      integer                :: j_1, j_2, j_3


      call eigen_get_procs( nprocs, nprow, npcol )
      call eigen_get_id   ( iam,    myrow, mycol )

      ! zero padding
      i_2 = eigen_loop_start( 1, 'Y' )
      i_3 = eigen_loop_end  ( n, 'Y' )
      j_2 = eigen_loop_end  ( n, 'X' ) + 1
      j_3 = nbb
!$OMP PARALLEL DO PRIVATE(j_1)
      do i_1=i_2,i_3
!OCL XFILL
         do j_1=j_2,j_3
            a(j_1+lda*(i_1-1)) = zero
         end do
      end do
!$OMP END PARALLEL DO
      i_2 = eigen_loop_end  ( n, 'Y' ) + 1
      i_3 = nbb
      j_2 = 1
      j_3 = nbb
!$OMP PARALLEL DO PRIVATE(j_1)
      do i_1=i_2,i_3
!OCL XFILL
         do j_1=j_2,j_3
            a(j_1+lda*(i_1-1)) = zero
         end do
      end do
!$OMP END PARALLEL DO


      return
      end subroutine zero_padding

      subroutine diag_mult( n, a, lda, w, b, ldb, mode )
!$    use OMP_LIB
      use eigen_libs_mod
      implicit NONE

      integer, intent(IN)    :: n, lda, ldb
      real(8), intent(IN)    :: a(*), w(*)
      real(8), intent(INOUT) :: b(*)
      real(8), intent(IN), optional :: mode

      integer                :: i, m
      integer                :: nprocs, npcol, nprow
      integer                :: iam, mycol, myrow
      integer                :: i_1, i_2, i_3
      integer                :: j_1, j_2, j_3
      real(8)                :: r


      if ( PRESENT(mode) ) then
         m = INT(mode * 2)
      else
         m = 2
      end if
      if ( m == 0 ) return

      call eigen_get_procs( nprocs, nprow, npcol )
      call eigen_get_id   ( iam,    myrow, mycol )

      i_2 = eigen_loop_start( 1, 'Y' )
      i_3 = eigen_loop_end  ( n, 'Y' )
      j_2 = eigen_loop_start( 1, 'X' )
      j_3 = eigen_loop_end  ( n, 'X' )
!$OMP PARALLEL DO PRIVATE(i,j_1,r)
      do i_1=i_2,i_3
         i = eigen_translate_l2g( i_1, 'Y' )

         if ( m == -2 ) then
            r = 1d0/w(i)
         else if ( m == -1 ) then
            r = 1d0/SQRT(w(i))
         else if ( m == 1 ) then
            r = SQRT(w(i))
         else if ( m == 2 ) then
            r = w(i)
         else
            r = ONE
         end if

         do j_1=j_2,j_3
            b(j_1+ldb*(i_1-1)) = a(j_1+lda*(i_1-1)) * r
         end do
      end do
!$OMP END PARALLEL DO

      return
      end subroutine diag_mult

      subroutine mat_lcopy( n, a, lda, b, ldb )
!$    use OMP_LIB
      use eigen_libs_mod
      implicit NONE

      integer, intent(IN)    :: n, lda, ldb
      real(8), intent(IN)    :: a(*)
      real(8), intent(INOUT) :: b(*)

      integer                :: i
      integer                :: nprocs, npcol, nprow
      integer                :: iam, mycol, myrow
      integer                :: i_1, i_2, i_3
      integer                :: j_1, j_2, j_3


      call eigen_get_procs( nprocs, nprow, npcol )
      call eigen_get_id   ( iam,    myrow, mycol )

      i_2 = eigen_loop_start( 1, 'Y' )
      i_3 = eigen_loop_end  ( n, 'Y' )
      j_2 = eigen_loop_start( 1, 'X' )
      j_3 = eigen_loop_end  ( n, 'X' )
!$OMP PARALLEL DO PRIVATE(i,j_1)
      do i_1=i_2,i_3
         do j_1=j_2,j_3
            b(j_1+ldb*(i_1-1)) = a(j_1+lda*(i_1-1))
         end do
      end do
!$OMP END PARALLEL DO

      return
      end subroutine mat_lcopy

      subroutine trpos_utol( n, a, lda, b, ldb )
!$    use OMP_LIB
      use eigen_libs_mod
      implicit NONE

      integer, intent(IN)    :: n, lda, ldb
      real(8), intent(IN)    :: a(*)
      real(8), intent(INOUT) :: b(*)

      integer                :: i
      integer                :: nprocs, npcol, nprow
      integer                :: iam, mycol, myrow
      integer                :: i_1, i_2, i_3
      integer                :: j_1, j_2, j_3


      call eigen_get_procs( nprocs, nprow, npcol )
      call eigen_get_id   ( iam,    myrow, mycol )

      i_2 = eigen_loop_start( 1, 'Y' )
      i_3 = eigen_loop_end  ( n, 'Y' )
!$OMP PARALLEL DO PRIVATE(i,j_1,j_2,j_3)
      do i_1=i_2,i_3
         i = eigen_translate_l2g( i_1, 'Y' )
         j_2 = eigen_loop_start ( i+1, 'X' )
         j_3 = eigen_loop_end   ( n,   'X' )
         do j_1=j_2,j_3
            b(j_1+(i_1-1)*ldb) = a(j_1+(i_1-1)*lda)
         end do
      end do
!$OMP END PARALLEL DO

      return
      end subroutine trpos_utol

