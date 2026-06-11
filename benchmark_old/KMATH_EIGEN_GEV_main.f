      use MPI
      use eigen_libs_mod
      implicit NONE

      logical,parameter :: check_flag = .true.
!      logical,parameter :: check_flag = .false.
      real(8),parameter :: one = 1.0D0
      real(8),allocatable :: a(:),b(:),z(:),w(:)
      integer :: ierr,inod,nnod
      integer :: n
      integer :: NPROW,NPCOL,nx,nm,nmw, ifl,nfl
      integer :: larray_a,larray_b,larray_z
      integer :: numroc,iceil
      external :: numroc,iceil
      real(8) :: s_time,e_time
      integer :: lda,ldb,ldz
!
!      call MPI_Init_thread( MPI_THREAD_MULTIPLE, nx, ierr )
      call MPI_Init_thread( MPI_THREAD_SERIALIZED, nx, ierr )
      call MPI_Comm_rank( MPI_COMM_WORLD, inod, ierr )
      call MPI_Comm_size( MPI_COMM_WORLD, nnod, ierr )

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
*
      call eigen_init( )
      call eigen_get_procs( ierr, NPROW, NPCOL )
*
      if ( inod == 0 ) then
         open(10,FILE='IN_GEV')
      end if
*
      do

         if ( inod == 0 ) then
            read(10,*) n
         end if
         call MPI_Bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if ( n < 1 ) exit
         if(inod==0)then
            write(*,*)"Matrix dimension = ",n
         end if

         call eigen_get_matdims( n, nm, nmw )

         lda = nm
         ldb = nm
         ldz = nm

         larray_a = lda*nmw
         larray_b = ldb*nmw
         larray_z = ldz*nmw

         allocate(a(larray_a), b(larray_b))
         allocate(z(larray_z), w(n))

         s_time = MPI_Wtime( )
         call mat_set( n, a, lda, 2 )
         call mat_set( n, b, ldb, 10)
         e_time = MPI_Wtime( )
         if(inod==0)then
            write(*,*)"Matset time:",e_time-s_time
         end if
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierr )
         s_time = MPI_Wtime( )
*
         call KMATH_EIGEN_GEV( n, a, lda, b, ldb, w, z, ldz )
*
         call MPI_Barrier( MPI_COMM_WORLD, ierr )
         e_time = MPI_Wtime( )
*-
         if(inod==0)then
            write(*,*)"Solver time:",e_time-s_time
         end if
*
         if(check_flag) then
         s_time = MPI_Wtime( )
            call mat_set( n, a, lda, 2 )
            call mat_set( n, b, ldb, 10)
         e_time = MPI_Wtime( )
         if(inod==0)then
            write(*,*)"Matset time:",e_time-s_time
         end if
            ifl=0
            nfl=0
         s_time = MPI_Wtime( )
            call KMATH_EIGEN_GEV_check(n,a,lda,b,ldb,w,z,ldz,ifl,nfl)
         e_time = MPI_Wtime( )
         if(inod==0)then
            write(*,*)"Checking time:",e_time-s_time
         end if
         end if

         deallocate(a,b,z,w)

      end do

      if ( inod == 0 ) then
         close(10)
      end if

      call MPI_Finalize( ierr )

      end

