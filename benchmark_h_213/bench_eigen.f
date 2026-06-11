      use MPI
      use eigen_libs_mod
      implicit none


      real(8), allocatable :: a(:,:),z(:,:)
      real(8),    allocatable :: w(:)
!
      real(8), allocatable :: a_(:,:),z_(:,:)
      real(8),    allocatable :: w_(:)
!
      integer :: x_nnod, y_nnod, x_inod, y_inod
      integer :: i_inod, i_nnod, inod, nnod

      integer :: n, nm, ny 
      integer :: ierr, istat
!
      integer :: i,j

      real(8) :: w_err_max, z_err_max
      real(8) :: d1, d2, d3, d4


      call mpi_init_thread( MPI_THREAD_MULTIPLE, i, ierr )
      call mpi_comm_rank( mpi_comm_world, i_inod, ierr )
      call mpi_comm_size( mpi_comm_world, i_nnod, ierr )
*
      if ( i_inod == 0 ) call system( "date" )
*
      call eigen_init( )

      call eigen_get_procs( nnod, x_nnod, y_nnod )
      call eigen_get_id   ( inod, x_inod, y_inod )
*
      if ( i_inod == 0 ) then
         open(10,file="IN")
      endif

      do
*
      if ( i_inod == 0 ) then
         read(10,*) n
      endif
      call MPI_Bcast( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
      if ( n < 1 ) exit

      call eigen_get_matdims( n, nm, ny )

      if ( i_inod == 0 ) then
         print*," n=",n," nm=",nm," ny=",ny
      end if
      
      istat = 0
      allocate(
     &     a(nm, ny), z(nm, ny), w(n),
     &     stat=istat)
      if ( istat /= 0 ) then
         print*,"Memory exhausted"
         call flush(6)
         goto 999
      end if
*     
      call mat_set( n, a, nm )
*-
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      d1 = MPI_WTIME( )
*
      call eigen_s( n, n, a, nm, w, z, nm,
     &     m_forward=48, m_backward=128 )
*
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      d2 = MPI_WTIME( )
*
      if ( i_inod == 0 ) then
         print*,"Matrix dimension = ",n
         print*,"Elapsed time for eigen_S = ",d2-d1," [sec]"
      end if
*-
      ! reproducible test
      istat = 0
      allocate(
     &     a_(nm, ny), z_(nm, ny), w_(n),
     &     stat=istat)
      if ( istat /= 0 ) then
         print*,"Memory exhausted"
         call flush(6)
         goto 999
      end if
*
      call mat_set( n, a_, nm )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      d3 = MPI_WTIME( )
      call eigen_s( n, n, a_, nm, w_, z_, nm,
     &     m_forward=48, m_backward=128 )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      d4 = MPI_WTIME( )
      if ( i_inod == 0 ) then
         print*,"                         = ",d4-d3," [sec]"
      end if
*-
      if ( i_inod == 0 ) then
         w_err_max = 0d0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:w_err_max)
         do i=1,n
            w_err_max = max( w_err_max, ABS(w(i)-w_(i)) )
         end do
!$OMP END PARALLEL DO
         print*," Repro test : max(w-w_)=",w_err_max
      end if
      z_err_max = 0d0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(max:z_err_max)
      do i=1,n
      if ( MOD(i-1,y_nnod)+1 == y_inod ) then
      do j=1,n
      if ( MOD(j-1,x_nnod)+1 == x_inod ) then
         z_err_max = max(z_err_max,
     &     ABS( z ( (j-1)/x_nnod+1, (i-1)/y_nnod+1 ) -
     &          z_( (j-1)/x_nnod+1, (i-1)/y_nnod+1 ) ) )
      end if
      end do
      end if
      end do
!$OMP END PARALLEL DO
      call MPI_Allreduce(MPI_IN_PLACE,z_err_max,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,MPI_COMM_WORLD, ierr)
      if ( i_inod == 0 ) then
         print*," Repro test : max(z-z_)=",z_err_max
      end if
      deallocate( a_, z_, w_ )
*-
      call mat_set( n, a, nm )
      call test_s(n, a, nm, w, z, nm)

      deallocate( a )
      deallocate( z )
      deallocate( w )

      enddo

      if ( i_inod == 0 ) then
         close(10)
      endif

  999 continue

      call eigen_free( )
*
      call MPI_Finalize( ierr )
      end
*
      subroutine test_s(n, a, lda, w, z, ldz)
      use MPI
      use eigen_libs_mod
      implicit none
* 
      integer :: n,lda,ldz, ny,lda_,ldz_
      real(8) :: a(lda,*),z(ldz,*)
      real(8) :: w(*)
*
      integer :: i,j
      integer :: i_1,i_2,i_3,i_4,j_1,j_2,j_3,j_4, k_1
      integer :: i_inod,nnod,x_nnod,y_nnod,inod,x_inod,y_inod
      integer :: w_comm, x_comm, y_comm
      integer :: NB, ierr

      integer :: lwork, lrwork, liwork, info
      real(8),allocatable :: a_(:,:),z_(:,:),work(:)
      real(8),allocatable :: w_(:), rwork(:)
      integer,allocatable :: iwork(:)
      real(8),allocatable :: a__(:,:),z__(:,:)
      real(8),allocatable :: w__(:)
      real(8) :: w_err_max,z_err_max, rdummy
      real(8) :: z1, z2, zdummy, z_err, signz
      real(8) :: rad1, rad2, sign1,sign2, d1,d2,d3,d4

      integer :: ictxt
      integer, parameter :: DESC_DIM = 9
      integer, dimension(DESC_DIM) :: desca, descz, desca_, descz_
      integer :: NP,NQ

      call mpi_comm_rank( MPI_COMM_WORLD, i_inod, ierr )
      call eigen_get_comm ( w_comm, x_comm, y_comm )
      call eigen_get_procs( nnod, x_nnod, y_nnod )
      call eigen_get_id   ( inod, x_inod, y_inod )


      NB = eigen_NB
      lda_ = max((((n-1)/NB)/x_nnod+1)*NB, lda)
      ny   =     (((n-1)/NB)/y_nnod+1)*NB
      ldz_ = lda_

      ictxt = eigen_get_blacs_context()
      call descinit( desca, n, n, 1, 1, 0, 0, ictxt, lda, ierr )
      call descinit( descz, n, n, 1, 1, 0, 0, ictxt, ldz, ierr )
      call descinit( desca_, n, n, NB, NB, 0, 0, ictxt, lda_, ierr )
      call descinit( descz_, n, n, NB, NB, 0, 0, ictxt, ldz_, ierr )

      ierr=0
      allocate( a_(lda_,ny), z_(ldz_,ny), w_(n), stat=ierr )
      if ( ierr /= 0 ) then
         call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
      endif
      ierr=0
      allocate( a__(lda_,ny), z__(ldz_,ny), w__(n), stat=ierr )
      if ( ierr /= 0 ) then
         call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
      endif

      a_ (:,:) = 0d0
      z_ (:,:) = 0d0
      a__(:,:) = 0d0
      z__(:,:) = 0d0
      w_ (:)   = 0d0
      w__(:)   = 0d0

! A:cyclic -> A:block(NB)
      call pdgemr2d( n,n, a,1,1,desca, a_,1,1,desca_, ictxt )
      a__(:,:) = a_(:,:)

      NP = (((n-1)/NB)/x_nnod+1)*NB
      NQ = (((n-1)/NB)/y_nnod+1)*NB
      lrwork = 3*n+max(NB*(NP+1), 3*NB)
      lwork  = max(1+6*n+2*NP*NQ, lrwork) + 2*n + NP*NQ
      liwork = 7*n+8*max(x_nnod,y_nnod)+2

      allocate( work(lwork), iwork(liwork), stat=ierr )
      if ( ierr /= 0 ) then
         call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
      endif

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      d1 = MPI_Wtime()

      call pdsyevd( 'V', 'U',
     &        n, a_,  1, 1, desca_, w_, z_, 1, 1, descz_,
     &        work, lwork, iwork, liwork, info)

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      d2 = MPI_Wtime()

      if ( inod == 1 ) then
         print*,"Elapsed time for PDSYEVD = ",d2-d1," [sec]"
      end if

      z__(:,:) = z_(:,:)
      z_ (:,:) = 0d0
      w__(:)   = w_(:)
      w_ (:)   = 0d0

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      d3 = MPI_Wtime()
      call pdsyevd( 'V', 'U',
     &        n, a__, 1, 1, desca_, w_, z_, 1, 1, descz_,
     &        work, lwork, iwork, liwork, info)
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      d4 = MPI_Wtime()
      if ( inod == 1 ) then
         print*,"                         = ",d4-d3," [sec]"
      end if

      deallocate(work, iwork)

      if ( inod == 1 ) then
         w_err_max = 0d0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:w_err_max)
         do i=1,n
            w_err_max = max( w_err_max, ABS(w_(i)-w__(i)) )
         end do
!$OMP END PARALLEL DO
         print*," Repro test : max(w-w_)=",w_err_max
      end if
      z_err_max = 0d0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(max:z_err_max)
      do i=1,ny
      do j=1,ldz_
         z_err_max = max(z_err_max, ABS(z_(j,i)-z__(j,i)))
      end do
      end do
!$OMP END PARALLEL DO
      call MPI_Allreduce(MPI_IN_PLACE,z_err_max,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,MPI_COMM_WORLD, ierr)
      if ( i_inod == 0 ) then
         print*," Repro test : max(z-z_)=",z_err_max
      end if
      deallocate( a__, z__, w__ )

      if ( inod == 1 ) then
         w_err_max = 0d0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:w_err_max)
         do i=1,n
            w_err_max = max(w_err_max,dabs(w(i) - w_(i)))
         end do
!$OMP END PARALLEL DO
         print*,"MAX Difference of Eigen Value = ",w_err_max
         if ( w_err_max > 1d-10 ) then
         print*,"*************** FAILED **********************"
         call MPI_abort(MPI_COMM_WORLD, 1, ierr)
         endif
      end if

! A:cyclic <- Z:block(NB)
      call pdgemr2d( n,n, z_,1,1,descz_, a,1,1,desca, ictxt )

      deallocate( a_, z_, w_ )

      j_2 = eigen_loop_start( 1, x_nnod, x_inod )
      j_3 = eigen_loop_end  ( n, x_nnod, x_inod )
      i_2 = eigen_loop_start( 1, y_nnod, y_inod )
      i_3 = eigen_loop_end  ( n, y_nnod, y_inod )

      lwork = max(1,i_3-i_2+1)
      allocate( work(lwork), stat=ierr )
      work(:) = 0d0

      if ( x_inod == 1 ) then
        if ( i_2 <= i_3 .and. j_2 <= j_3 ) then
          do i_1=i_2,i_3
! signz : a unitary rotation factor to a computed eigenvector
            z_err_max = 0d0; k_1 = -1
            do j_1=j_2,j_3
              z1 = z(j_1,i_1)
              z2 = a(j_2,i_1)
              rad1 = ABS( z1 )
              rad2 = ABS( z2 )
              if ( rad1 * rad2 > z_err_max ) then
                k_1 = j_1
                z_err_max = rad1 * rad2
              end if
            end do
            j_1 = k_1
            z1 = z(j_1,i_1)
            z2 = a(j_1,i_1)
            rad1 = ABS( z1 )
            rad2 = ABS( z2 )
            if ( rad1*rad2 == 0d0 ) then
               signz = DCMPLX(1d0,0d0)
            else
               z1 = z1 / rad1
               z2 = z2 / rad2
               signz = z2 / z1
            endif
            work(i_1) = signz
          end do
        end if
      end if

      if ( i_2 <= i_3 ) then
        call MPI_Bcast(work,i_3-i_2+1,MPI_DOUBLE_PRECISION,
     &       0, x_comm, ierr )
      end if

      z_err_max = 0d0
      if ( i_2 <= i_3 .and. j_2 <= j_3 ) then
!$OMP PARALLEL DO PRIVATE(i_1,j_1,rdummy) REDUCTION(max:z_err_max)
        do i_1=i_2,i_3
          do j_1=j_2,j_3
            rdummy = ABS(work(i_1) * z(j_1,i_1) - a(j_1,i_1))
            z_err_max = max(z_err_max, rdummy)
          end do
        end do
!$OMP END PARALLEL DO
      end if

      deallocate(work)

      call MPI_Allreduce( MPI_IN_PLACE, z_err_max,
     &     1, MPI_DOUBLE_PRECISION,
     &     MPI_MAX, MPI_COMM_WORLD, ierr )

      if ( inod == 1 ) then
         print*,"MAX Difference of Eigen Vector= ",z_err_max
         if ( z_err_max > 1d-10 ) then
         print*,"*************** FAILED **********************"
         print*,"  or would be near-by clustered eigenvalues  "
         print*,"*************** +SKIP+ **********************"
         endif
      end if

      end subroutine
