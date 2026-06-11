!--------+---------+---------+---------+---------+---------+---------+---------|
      subroutine KMATH_EIGEN_GEV( n, a, lda, b, ldb, w, z, ldz )
!
!
!   integer(In)    n        :: Matrix dimension
!   real(8)(InOut) a(lda,*) :: An input matrix to be solved for GEV of LHS.
!                              A is real symmetric matrix,
!                              and an upper triangular format distributed on
!                              a 2D cyclic-cyclic fashion.
!                              On exit, data is destroyed.
!   integer(In)    lda      :: Leading Dimension of the array A
!   real(8)(InOut) b(ldb,*) :: An input matrix to be solved for GEV of RHS.
!                              B is real positive definite symmetric matrix,
!                              and an upper triangular format distributed on
!                              a 2D cyclic-cyclic fashion.
!                              The transformation factor to SEP,
!                                 eig.vec(B) * diag(eig.val(B))^{-1/2}.
!                              is returned on exit.
!   integer(In)    ldb      :: Leading Dimension of the array B
!   real(8)(Out)   w(*)     :: An array where eigenvalues computed are stored
!   real(8)(Out)   z(ldz,*) :: An output matrix where eigenvectors are stored
!   integer(In)    ldz      :: Leading Dimension of the array Z
!
!
      use MPI
!$    use OMP_LIB
      use eigen_libs_mod
      use eigen_devel_mod
      implicit none

      integer, intent(InOut)    ::  n, lda, ldb, ldz
      real(8), intent(InOut) ::  a(*), b(*)
      real(8), intent(InOut) ::  w(*), z(*)

      integer :: nprocs1, nprow1,nprow, npcol1,npcol

      call eigen_get_procs( nprocs1, nprow1, npcol1 )
      nprow = nprow1; npcol = npcol1

      if ( nprow == npcol ) then

         call KMATH_EIGEN_GEV_1( n, a, lda, b, ldb, w, z, ldz )

      else

         call KMATH_EIGEN_GEV_2( n, a, lda, b, ldb, w, z, ldz )

      end if

      return
      contains
!--------+---------+---------+---------+---------+---------+---------+---------|

      include 'KMATH_EIGEN_GEV_1.F'

      include 'KMATH_EIGEN_GEV_2.F'

      include 'KMATH_EIGEN_GEV_misc.F'


!--------+---------+---------+---------+---------+---------+---------+---------|
      end subroutine
!--------+---------+---------+---------+---------+---------+---------+---------|

