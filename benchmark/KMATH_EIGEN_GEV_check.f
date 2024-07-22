      subroutine KMATH_EIGEN_GEV_check(N,A,lda,B,ldb,w,z,ldz,
     &     ifl,nfl)
      use MPI
!     $    use OMP_LIB
      use eigen_libs_mod
      use eigen_blacs_mod
      implicit none

      integer :: N
      real(8),intent(in) :: w(*)
      real(8),intent(inout) :: z(*)
      real(8),intent(in) :: A(*),B(*)
      integer,intent(in) :: lda,ldb,ldz,nfl,ifl

      real(8),parameter :: one = 1.0D0
      real(8),parameter :: zero = 0.0D0

      integer :: NPROW,NPCOL,larray
      integer :: ictxt
      integer :: iam1,nprocs1
      integer :: nprow_i,npcol_i,myrow1,mycol1
      integer :: numr_row,numr_col
      integer :: ierr,inod
      integer :: i,j
      real(8) :: s_time,e_time
      integer :: nm, nx, NB

      integer :: numroc,iceil,indxl2g
      external :: numroc,iceil,indxl2g
      real(8) :: pdlange
      external :: pdlange

      real(8) :: norm2
      integer :: gcol
      real(8) :: DUMWORK
      real(8), allocatable :: Y1(:), Y2(:)
      real(8), allocatable :: AB(:), BB(:), ZB(:)

      integer, parameter :: DESC_DIM = 9
      integer, dimension(DESC_DIM) :: a0desc, b0desc, z0desc
      integer, dimension(DESC_DIM) :: z1desc

      integer :: i_1,i_2,i_3
      integer :: j_1,j_2,j_3
      integer :: iam, myrow, mycol      

!     for pdsygvx
      call eigen_get_id   ( iam, myrow, mycol )            

      call eigen_get_procs( nprocs1, NPROW, NPCOL )            

      nprow_i=NPROW
      npcol_i=NPCOL
     
      call mpi_barrier(MPI_COMM_WORLD, ierr )
      call eigen_get_id( inod, myrow1, mycol1 )      
      call mpi_barrier(MPI_COMM_WORLD, ierr )
      
      inod=inod-1; myrow1=myrow1-1; mycol1=mycol1-1
      ictxt = eigen_get_blacs_context()

      call eigen_get_matdims( n, nm, nx )

      larray = nm * nx

*     
      call descinit( a0desc, n, n, 1, 1, 0, 0, ictxt, lda, ierr )
      call descinit( b0desc, n, n, 1, 1, 0, 0, ictxt, ldb, ierr )
      call descinit( z0desc, n, n, 1, 1, 0, 0, ictxt, ldz, ierr )

      NB = 64
      nm = (((n-1)/NB)/NPROW+1)*NB

      call descinit( z1desc, n, n, NB, NB, 0, 0, ictxt, nm, ierr )

      allocate(AB(larray), ZB(larray))


      call pdgemr2d( n, n, a, 1, 1, a0desc,
     &     ab, 1, 1, z1desc, ictxt)
      call pdgemr2d( n, n, z, 1, 1, z0desc,
     &     zb, 1, 1, z1desc, ictxt)

      allocate(Y1(larray))

!     |AX-BXW|
      call pdgemm('N','N',n,n,n,
     &     one,
     &     ab, 1, 1, z1desc,
     &     zb, 1, 1, z1desc,
     &     zero,
     &     Y1, 1, 1, z1desc)

      
      deallocate( AB )


      i_2 = eigen_loop_start( 1, 'Y' )
      i_3 = eigen_loop_end  ( n, 'Y' )
!$OMP PARALLEL DO PRIVATE(i,j_1,j_2,j_3)
      do i_1=i_2,i_3
         i = eigen_translate_l2g( i_1, 'Y' )
         j_2 = eigen_loop_start( 1, 'X' )
         j_3 = eigen_loop_end  ( n, 'X' )
         do j_1=j_2,j_3
            z(j_1+(i_1-1)*ldz) = z(j_1+(i_1-1)*ldz) * w(i)
         end do
      end do
!$OMP END PARALLEL DO

      allocate(Y2(larray))

      call pdgemr2d( n, n, z, 1, 1, z0desc,
     &     Y2, 1, 1, z1desc, ictxt)

      allocate(BB(larray))

      call pdgemr2d( n, n, b, 1, 1, b0desc,
     &     bb, 1, 1, z1desc, ictxt)

      call pdgemm('N','N',n,n,n,
     &     one,
     &     bb, 1, 1, z1desc,
     &     Y2, 1, 1, z1desc,
     &     -one,
     &     Y1, 1, 1, z1desc)

      norm2 = pdlange( 'F', n, n, Y1, 1, 1, z1desc, DUMWORK )

      if ( inod == 0 ) then
         print*,'|AX-BXW|_F=',norm2
      end if

!     |XBX-I|

      call pdgemm('T','N',n,n,n,
     &     one,
     &     zb, 1, 1, z1desc,
     &     bb, 1, 1, z1desc,
     &     zero,
     &     Y1, 1, 1, z1desc)

      call pdgemm('N','N',n,n,n,
     &     one,
     &     Y1, 1, 1, z1desc,
     &     zb, 1, 1, z1desc,
     &     zero,
     &     Y2, 1, 1, z1desc)

      call pdgemr2d( n, n, Y2, 1, 1, z1desc,
     &     z, 1, 1, z0desc, ictxt)

      i_2 = eigen_loop_start( 1, 'Y' )
      i_3 = eigen_loop_end  ( n, 'Y' )
!     $OMP PARALLEL DO PRIVATE(i,j_1,j_2,j_3)
      do i_1=i_2,i_3
         i = eigen_translate_l2g( i_1, 'Y' )
         j_1 = eigen_owner_index( i, 'X' )
         if (j_1 > 0) then
            z(j_1+(i_1-1)*ldz) = z(j_1+(i_1-1)*ldz) - one
         end if
      end do
!     $OMP END PARALLEL DO

      call pdgemr2d( n, n, z, 1, 1, z0desc,
     &     Y2, 1, 1, z1desc, ictxt)

      
      norm2 = pdlange( 'F', n, n, Y2, 1, 1, z1desc, DUMWORK )

      if ( inod == 0 ) then
         print*,'|XBX-I|_F=',norm2
      end if
      call mpi_barrier(MPI_COMM_WORLD, ierr )
      deallocate( BB, ZB )
      deallocate( Y1, Y2 )

      return
      end subroutine

