      subroutine KMATH_EIGEN_GEV_2( n, a, lda, b, ldb, w, z, ldz )

      use MPI
!$    use OMP_LIB
      use eigen_libs_mod
      use eigen_blacs_mod
      implicit none

      integer, intent(InOut)    ::  n, lda, ldb, ldz
      real(8), intent(InOut) ::  a(*), b(*)
      real(8), intent(InOut) ::  w(*), z(*)

      real(8),parameter :: one = 1.0D0
      real(8),parameter :: zero = 0.0D0

      integer :: i,j
      integer :: nx,nm,NB,nmz,nmw,larray
      integer :: larray1
      real(8),allocatable :: c(:)

      logical :: sflag
      integer,parameter :: nbb = 1  ! eigen_s's block size
      integer,parameter :: nbb3 = eigen_NB ! block size for reditributed matrix
      integer :: ictxt
      integer :: iam,nprocs
      integer :: nprow,npcol, myrow, mycol
      integer :: j_1,j_2,j_3,i_1,i_2,i_3
      integer :: ierr

      integer,parameter :: DESC_DIM = 9
      integer, dimension(DESC_DIM) :: FBdesc, adesc, bdesc, zdesc


      real(8) :: t1,t2,tx(20)


      tx = ZERO


      ictxt = eigen_get_blacs_context()

      call eigen_get_procs( nprocs, nprow, npcol )
      call eigen_get_id( iam, myrow, mycol )

      NB  = nbb3

      nmz = ((n-1)/nprow+1)
      nmz = ((nmz-1)/NB+1)*NB
      nx  = nmz
      nmz = ((n-1)/NB+1)
      nmz = ((nmz-1)/nprow+1)*NB
      nmz = MAX(nmz, nx)

      nmw = ((n-1)/npcol+1)
      nmw = ((nmw-1)/NB+1)*NB
      nx  = nmz
      nmw = ((n-1)/NB+1)
      nmw = ((nmw-1)/nprow+1)*NB
      nmw = MAX(nmw, nx)

      nm  = nmz
      larray = nmz * nmw

      call descinit( FBdesc, n, n,  nbb3,nbb3, 0,0, ictxt, nm,  ierr )
      call descinit( adesc,  n, n,  nbb, nbb,  0,0, ictxt, lda, ierr )
      call descinit( bdesc,  n, n,  nbb, nbb,  0,0, ictxt, ldb, ierr )
      call descinit( zdesc,  n, n,  nbb, nbb,  0,0, ictxt, ldz, ierr )


      t1 = MPI_Wtime()
      call pdtran( n, n, one, a, 1, 1, adesc, zero, z, 1, 1, zdesc )
      call mat_lcopy( n, z, ldz, a, lda )
      t2 = MPI_Wtime()
      tx(1) = t2-t1


      t1 = MPI_Wtime()
!      call eigen_sx( n, n, b(1), ldb, w(1), z(1), ldz )
      call eigen_s( n, n, b(1), ldb, w(1), z(1), ldz )
      t2 = MPI_Wtime()
      tx(2) = t2-t1


      if ( w(1) <= zero ) then
         if ( iam == 1 ) then
            print*,"Matrix B is not positive definite!"
            return
         end if
      end if


!      call diag_mult( n, z, ldz, w, b, ldb, mode=-one/2 )
      call diag_mult( n, z, ldz, w, b, ldb, -one/2 )


      ! data redistribution
      t1 = MPI_Wtime()
      call pdgemr2d( n, n, b, 1, 1, bdesc,
     &                        z, 1, 1, FBdesc, ictxt)
      t2 = MPI_Wtime()
      tx(3) = t2-t1


      allocate(c(larray))

      t1 = MPI_Wtime()
      call pdgemr2d( n, n, a, 1, 1, adesc,
     &                        c, 1, 1, FBdesc, ictxt)
      t2 = MPI_Wtime()
      tx(4) = t2-t1


      !calculate  = z^T * c * z
      t1 = MPI_Wtime()
      call pdgemm('N','N',n,n,n,
     &      one,
     &         c, 1, 1, FBdesc,
     &         z, 1, 1, FBdesc,
     &      zero,
     &         a, 1, 1, FBdesc)
      t2 = MPI_Wtime()
      tx(5) = t2-t1


      t1 = MPI_Wtime()
      call pdgemm('T','N',n,n,n,
     &      one,
     &         z, 1, 1, FBdesc,
     &         a, 1, 1, FBdesc,
     &      zero,
     &         c, 1, 1, FBdesc)
      t2 = MPI_Wtime()
      tx(6) = t2-t1


      ! back to original form
      t1 = MPI_Wtime()
      call pdgemr2d( n, n, c, 1, 1, FBdesc,
     &                        a, 1, 1, zdesc, ictxt)
      t2 = MPI_Wtime()
      tx(7) = t2-t1


      deallocate(c)

      t1 = MPI_Wtime()
!      call eigen_sx( n, n, a(1), ldz, w(1), z(1), ldz )
      call eigen_s( n, n, a(1), ldz, w(1), z(1), ldz )
      t2 = MPI_Wtime()
      tx(8) = t2-t1


      t1 = MPI_Wtime()
      call pdgemr2d( n, n, z, 1, 1, zdesc,
     &                        a, 1, 1, FBdesc, ictxt)
      t2 = MPI_Wtime()
      tx(9) = t2-t1


      t1 = MPI_Wtime()
      call pdgemr2d( n, n, b, 1, 1, bdesc,
     &                        z, 1, 1, FBdesc, ictxt)
      t2 = MPI_Wtime()
      tx(10) = t2-t1


      allocate(c(larray))

      t1 = MPI_Wtime()
      call pdgemm('N','N',n,n,n,
     &      one,
     &         z, 1, 1, FBdesc,
     &         a, 1, 1, FBdesc,
     &      zero,
     &         c, 1, 1, FBdesc)
      t2 = MPI_Wtime()
      tx(11) = t2-t1


      t1 = MPI_Wtime()
!      call pdgemr2d( n, n, c, 1, 1, FBdesc,
!     &                        z, 1, 1, zdesc, ictxt)

!      call dc_redist1C( n, NB, c, a, nm, z, nbb*nbb )
!      call dc_redist2C( n, NB, a, nm, z, ldz, c, nbb*nbb )

      call dc_redist1( n, n, NB, c, a, nm, z, nbb*nbb )
      call dc_redist2( n, n, NB, a, nm, z, ldz, c, nbb*nbb )


      t2 = MPI_Wtime()
      tx(12) = t2-t1


      deallocate(c)


      if ( iam == 1 ) then
         print*,1,'PDTRAN    ',tx(1)
         print*,2,'Eigen_sx1 ',tx(2)
         print*,3,'PDGEMR2D1 ',tx(3)
         print*,4,'PDGEMR2D2 ',tx(4)
         print*,5,'PDGEMM1a  ',tx(5)
         print*,6,'PDGEMM1b  ',tx(6)
         print*,7,'PDGEMR2D3 ',tx(7)
         print*,8,'Eigen_sx2 ',tx(8)
         print*,9,'PDGEMR2D4 ',tx(9)
         print*,10,'PDGEMR2D5 ',tx(10)
         print*,11,'PDGEMM2   ',tx(11)
         print*,12,'PDGEMR2D6 ',tx(12)
      end if


      return
      end subroutine

