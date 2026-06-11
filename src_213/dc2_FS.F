!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: dc2_mod
!
! Purpose
! =======
!
! dc2_mod manages the modules for divide and conquer algrithm for
! solving a symmetric tridiagonal matrix.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module dc2_FS_mod

      use eigen_devel_mod
      use eigen_dc_mod
      use FS_EDC_mod
      use eigen_libs_mod
!$    use omp_lib

      implicit none
      private

      public :: dc2_FS

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_dc2
!
! Purpose
! =======
!
! eigen_dc2 invokes the main body of the divide and conquer solver,
! eigen_dc2_body, to solve the eigenpairs of the symmetric tridiagonal
! matrix.
!
!
! Arguments
! =========
!
! n       (input) integer
!         The dimension of the symmetric tridiagonal matrix. N >= 0.
!
! nvec    (input) integer
!         The number of eigenmodes to be computed. N >= NVEC >= 0.
!
! d       (input/output) real(8) array, dimension(n)
!         On entry, d contains the diagonal elements of the symmetric
!         tridiagonal matrix.
!         On exit, d contains eigenvalues of the input matrix.
!
! e       (input/output) real(8) array, dimension(n-1)
!         On entry, e contains the off-diagonal elements of the
!         symmetric tridiagonal matrix.
!         On exit, values has been destroyed.
!
! z       (output) real(8) array, dimension(ldz,(n-1)/y_nnod+1)
!         z returns the eigenvectors of the input matrix.
!
! ldz     (input) integer
!         The leading dimension of the array z. ldz >= ceil(N/x_nnod).
!
! info    (output) integer
!         = 0: successful exit
!         < 0: error status as same as scalapack
!         > 0: error status as same as scalapack
!
! ret     (output) real(8)
!         The number of floating point operations.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine dc2_FS(n, nvec, d, e, z, ldz, info, ret)
      implicit none
      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      real(8), intent(inout) :: d(1:n)
      real(8), intent(inout) :: e(1:n-1)
      real(8), intent(out)   :: z(1:ldz,*)
      integer, intent(in)    :: ldz
      integer(8), intent(out)   :: info
      real(8), intent(out)   :: ret

!     Parameters BLACS array descritor(the position of entry tags), etc
      integer, parameter     :: BLOCK_CYCLIC_2D = 1
      integer, parameter     :: DLEN_  = 9
      integer, parameter     :: DTYPE_ = 1
      integer, parameter     :: CTXT_  = 2
      integer, parameter     :: M_     = 3
      integer, parameter     :: N_     = 4
      integer, parameter     :: MB_    = 5
      integer, parameter     :: NB_    = 6
      integer, parameter     :: RSRC_  = 7
      integer, parameter     :: CSRC_  = 8
      integer, parameter     :: LLD_   = 9

      logical, parameter     :: USE_MY_REDIST = .true.

      integer                :: descz(DLEN_)
      integer                :: descw(DLEN_)

      integer                :: i, j, k, nx, nxv, istat, nb, lddz, lddw
      integer                :: np, nq, npcol, nprow
      integer                :: nprocs, iam, mycol, myrow
      integer                :: ictxt, ierr
      integer(8)             :: lwork, liwork
      integer(8)             :: lwork_, liwork_
      integer                :: eigen_comm, eigen_x_comm, eigen_y_comm            


      real(8), allocatable   :: work(:)
      integer, allocatable   :: iwork(:)

      integer                :: local_size, local_rank
      integer                :: i_2, i_3

#if defined(__INTEL_COMPILER) && USE_MKL
      integer                :: mkl_mode
      integer, external      :: mkl_get_dynamic
#endif

      integer, external      :: NUMROC

      type(FS_prof)          :: prof
!

      flops = 0d0
      dgemm_time = 0d0
      p_time0 = 0d0
      p_timer = 0d0
      p_time2 = 0d0
      p_time3 = 0d0
      p_times = 0d0
      p_timez = 0d0

      call eigen_timer_reset(1, 0, 0, 0)

      nprocs = TRD_nnod
      iam    = TRD_inod-1

!      nprow = x_nnod
!      npcol = y_nnod
!      myrow = x_inod-1
!      mycol = y_inod-1

!      call FS_WorkSize(n, lwork, liwork)
      call eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm)
      call FS_WorkSize(n, lwork_, liwork_)
      if(.not.FS_COMM_MEMBER)then
         lwork_  = 0
         liwork_ = 0
      endif
      call MPI_Allreduce(lwork_,lwork,1,MPI_INTEGER8,
     &     MPI_MAX,eigen_comm,ierr)
      call MPI_Allreduce(liwork_,liwork,1,MPI_INTEGER8,
     &     MPI_MAX,eigen_comm,ierr)

      allocate(work(lwork), iwork(liwork), stat=istat)
      if (istat /= 0) then
        call eigen_abort( "Memory exhausted [DC-FS]", 1 )
      end if

#if defined(__INTEL_COMPILER) && USE_MKL
      mkl_mode = mkl_get_dynamic()
      call MKL_SET_DYNAMIC(0)
#endif
#if TIMER_PRINT
      call FS_prof_init(prof)
#endif
      call FS_EDC(n, d(1), e(1), z(1,1), ldz, 
     &     work(1), lwork, iwork(1), liwork, info, prof)
#if TIMER_PRINT
      call FS_prof_finalize(prof)
#endif
#if TIMER_PRINT>1
      p_time0 = prof%region_time(21)
      p_timer = prof%region_time(70)
      p_time2 = prof%region_time(50)
      p_time3 = prof%region_time(60)
      dgemm_time = prof%region_time(67)
      p_timez = prof%region_time(40)
#endif
#if defined(__INTEL_COMPILER) && USE_MKL
      call MKL_SET_DYNAMIC(mkl_mode)
#endif

#if TIMER_PRINT>1
      if (iam == 0) print*,"FS_EDC     ", prof%region_time(10)
#endif

!     freeing working arrays
      deallocate(work)
      deallocate(iwork)

      comm_time_dc =
     &     eigen_timer_print('EigenExa (D-and-C for tri-diag)')

#if TIMER_PRINT>1
      if (iam==0) print*,"FS_dividing", p_time0
      if (iam==0) print*,"FS_PDLASRT ", p_timer
      if (iam==0) print*,"FS_PDLAED2 ", p_time2
      if (iam==0) print*,"FS_PDLAED3 ", p_time3
      if (iam==0) print*,"FS_PDLAEDZ ", p_timez
      if (iam==0) print*,"DGEMM      ", dgemm_time
#endif

      call MPI_Allreduce(flops,ret,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,eigen_comm,ierr)
      ! dgemm_time ! flops/dgemm_time

      return

      end subroutine

      end module dc2_FS_mod
