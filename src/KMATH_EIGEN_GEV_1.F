      subroutine KMATH_EIGEN_GEV_1( n, a, lda, b, ldb, w, z, ldz )

      use MPI
!$    use OMP_LIB
      use eigen_libs_mod
      use eigen_blacs_mod
      implicit none

      integer, intent(InOut) ::  n, lda, ldb, ldz
      real(8), intent(InOut) ::  a(*), b(*)
      real(8), intent(InOut) ::  w(*), z(*)

      real(8),parameter :: one = 1.0D0
      real(8),parameter :: zero = 0.0D0

      integer :: i,j,k
      integer :: n2, nx, nm, larray
      real(8), allocatable :: c(:)

      integer,parameter :: nbb  = 1  ! NB for eigen_s
      integer           :: nbb2      ! NB for block distributed matrix
      integer,parameter :: nbb3 = 64 ! NB for reditributed matrix

      integer :: ictxt
      integer :: iam, nprocs, ierr
      integer :: nprow, npcol, myrow, mycol
      integer :: j_1,j_2,j_3,i_1,i_2,i_3

      integer, parameter :: DESC_DIM = 9
      integer, dimension(DESC_DIM) :: a0desc, z0desc, b1desc, z1desc
      integer, dimension(DESC_DIM) :: FBdesc

      real(8) :: t1,t2,tx(20)


      tx = ZERO


      ictxt = eigen_get_blacs_context()

      call eigen_get_procs( nprocs, nprow, npcol )
      call eigen_get_id   ( iam, myrow, mycol )

      nbb2   = eigen_translate_g2l( n, npcol, mycol )
      n2     = nbb2 * npcol

      call eigen_get_matdims( n, nm, nx )


      call descinit( FBdesc, n2,n2, nbb2,nbb2, 0,0, ictxt, nm,  ierr )
      call descinit( a0desc, n, n,  nbb, nbb,  0,0, ictxt, lda, ierr )
      call descinit( z0desc, n, n,  nbb, nbb,  0,0, ictxt, ldz, ierr )
      call descinit( b1desc, n2,n2, nbb2,nbb2, 0,0, ictxt, ldb, ierr )
      call descinit( z1desc, n2,n2, nbb2,nbb2, 0,0, ictxt, ldz, ierr )


      t1 = MPI_Wtime()
      call pdtran( n, n, one, a, 1, 1, a0desc, zero, z, 1, 1, z0desc )
      call trpos_utol( n, z, ldz, a, lda )
      t2 = MPI_Wtime()
      tx(1) = t2-t1


      t1 = MPI_Wtime()
!      call eigen_sx( n, n, b(1), ldb, w(1), z(1), ldz, mode='X' )
      call eigen_s( n, n, b(1), ldb, w(1), z(1), ldz, mode='X' )
      t2 = MPI_Wtime()
      tx(2) = t2-t1


      if ( w(1) <= zero ) then
         if ( iam == 1 ) then
            print*,"Matrix B is not positive definite!"
            return
         end if
      end if


      larray = nm * nx
      allocate(c(larray), stat=ierr)
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error [GEV].", 1 )
      end if


      ! z = z * eigenvalue(w)^(-1/2)
      if ( nm /= lda ) then
!         call diag_mult( n, z, ldz, w, b, ldb, mode=-one/2 )
         call diag_mult( n, z, ldz, w, b, ldb, -one/2 )
         call mat_lcopy( n, a, lda, c, nm )
         call mat_lcopy( n, c, nm, a, nm )
      else
!         call diag_mult( n, z, ldz, w, b, ldb, mode=-one/2 )
         call diag_mult( n, z, ldz, w, b, ldb, -one/2 )
      end if
      call zero_padding( n, nbb2, b, ldb )
      call zero_padding( n, nbb2, a, nm )


      !calculate a' = b^T * a * b
      t1 = MPI_Wtime()
      call pdgemm('N','N', n2, n2, n2,
     &      one,
     &         a, 1, 1, FBdesc,
     &         b, 1, 1, b1desc,
     &      zero,
     &         c, 1, 1, FBdesc)
      t2 = MPI_Wtime()
      tx(3) = t2-t1
      t1 = MPI_Wtime()
      call pdgemm('T','N',n2,n2,n2,
     &      one,
     &         b, 1, 1, b1desc,
     &         c, 1, 1, FBdesc,
     &      zero,
     &         z, 1, 1, FBdesc)
      t2 = MPI_Wtime()
      tx(4) = t2-t1


      deallocate(c)


      ! eigen_s routine for a'
      t1 = MPI_Wtime()
!      call eigen_sx( n, n, z(1), nm, w(1), a(1), nm, mode='X' )
      call eigen_s( n, n, z(1), nm, w(1), a(1), nm, mode='X' )
      t2 = MPI_Wtime()
      tx(5) = t2-t1


      call zero_padding( n, nbb2, a, nm )


      ! Back-transform eigenvector
      t1 = MPI_Wtime()
      call pdgemm('N','N',n2,n2,n2,
     &      one,
     &         b, 1, 1, b1desc,
     &         a, 1, 1, FBdesc,
     &      zero,
     &         z, 1, 1, z1desc)
      t2 = MPI_Wtime()
      tx(6) = t2-t1


      if ( iam == 1 ) then
         print*,1,'PDTRAN    ',tx(1)
         print*,2,'Eigen_sx1 ',tx(2)
         print*,3,'PDGEMM1a  ',tx(3)
         print*,4,'PDGEMM1b  ',tx(4)
         print*,5,'Eigen_sx2 ',tx(5)
         print*,6,'PDGEMM2   ',tx(6)
      end if


      return
      end subroutine

