!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: dcx_mod
!
! Purpose
! =======
!
! dcx_mod manages the modules for divide and conquer algrithm for
! solving a symmetric penta-diagonal matrix.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module dcx_mod

      use my_pdsxedc_mod
      use eigen_libs_mod
      use comm_mod
      use eigen_blacs_mod, only : eigen_get_blacs_context
      use eigen_devel_mod
      use eigen_dc_mod
      use mpi
!$    use omp_lib

      implicit none
      private

      public :: eigen_dcx

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_dcx
!
! Purpose
! =======
!
! eigen_dcx invokes the main body of the divide and conquer solver,
! eigen_dcx_body, to solve the eigenpairs of the symmetric
! penta-diagonal matrix.
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
! e       (input/output) real(8) array, dimension(2*lde)
!         On entry, e contains the off-diagonal elements of the
!         symmetric tridiagonal matrix.
!         On exit, values has been destroyed.
!
! lde     (input) integer
!         The leading dimension of the array z. lde >= N.
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

      subroutine eigen_dcx(n, nvec, d, e, lde, z, ldz, info, ret)

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      real(8), intent(inout) :: d(1:n)
      real(8), intent(inout) :: e(1:2*lde)
      integer, intent(in)    :: lde
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

      integer :: info_i4

      real(8), allocatable   :: work(:)
      integer, allocatable   :: iwork(:)

      integer                :: local_size, local_rank
      integer                :: i_2, i_3

      external               :: descinit
      external               :: dc_redist1, dc_redist2
#if defined(__INTEL_COMPILER) && USE_MKL
      integer                :: mkl_mode
      integer, external      :: mkl_get_dynamic
#endif
      integer, external      :: NUMROC

#if TIMER_PRINT
      real(8)                :: d1,d2
#endif


      flops = 0d0
      dgemm_time = 0d0
      
      call eigen_timer_reset(1, 1, 1, 1)

#if TIMER_PRINT>1
      d1 = eigen_get_wtime()
#endif

      nprocs = TRD_nnod
      iam    = TRD_inod - 1

      nprow = x_nnod
      npcol = y_nnod

      ictxt = eigen_get_blacs_context()

      myrow = x_inod - 1
      mycol = y_inod - 1


!     BLACS array registration
      nb = eigen_nb
      nb = min(nb, N)

      np = NUMROC(n, nb, myrow, 0, nprow)
      nq = NUMROC(n, nb, mycol, 0, npcol)

      lddz = (n-1)/nprow+1
      lddz = ((lddz-1)/nb+1)*nb
      lddw = (n-1)/npcol+1
      lddw = ((lddw-1)/nb+1)*nb
      nx   = (n-1)/npcol+1
      nxv  = (nvec-1)/npcol+1

      call descinit(descz, n, n, nb, nb, 0, 0, ictxt, lddz, INFO_i4)

!     preparing working arrays
      lwork  = max(1+6*N+2*INT(np,8)*INT(nq+max(nq,nb),8),
     &         INT(lddz,8)*INT(lddw,8),
     &         INT(ldz,8)*INT(nx,8))
      liwork = 2+7*n+8*npcol
      allocate(work(lwork), iwork(liwork), stat=istat)
      if (istat /= 0) then
        call eigen_abort( "Memory exhausted [DC-x]", 1 )
      end if

!     Somehow, Z must be nullified (Originally next loop is not required.)
!$OMP PARALLEL DO
      do i = 1, lddz*lddw
         z(i,1) = ZERO
      end do
!$OMP END PARALLEL DO
#if TIMER_PRINT>1
      d2 = eigen_get_wtime()
      if (iam == 0) print*,"before PDSTEDC",d2-d1
#endif

#if defined(__INTEL_COMPILER) && USE_MKL
      mkl_mode = mkl_get_dynamic()
      call mkl_set_dynamic(0)
#endif

#if TIMER_PRINT>1
      d1 = eigen_get_wtime()
#endif
      call my_pdsxedc('I',2, n, d(1), e(1), lde, z(1,1), 1, 1, descz,
     &     WORK(1), lwork, IWORK(1), liwork, INFO)
#if TIMER_PRINT>1
      d2 = eigen_get_wtime()
#endif

#if defined(__INTEL_COMPILER) && USE_MKL
      call mkl_set_dynamic(mkl_mode)
#endif

#if TIMER_PRINT>1
      if (iam == 0) print*,"PDSTEDC",d2-d1
#endif

      if (nb == 1 .or. TRD_nnod == 1) then
        if (lddz /= ldz) then
          k = min(lddz,ldz)
!$OMP PARALLEL
!$OMP+         PRIVATE(i,j)
!$OMP DO
          do i = 1, lddz*lddw
            work(i) = z(i,1)
          end do
!$OMP ENDDO
!$OMP BARRIER
!$OMP DO
          do i = 1,nxv
            j = (i-1)*lddz
            z(1:k,i) = work(1+j:k+j)
          end do
!$OMP ENDDO
          if (lddz < ldz) then
!$OMP BARRIER
!$OMP DO
            do i = 1,nxv
              z(lddz+1:ldz,i) = ZERO
            end do
!$OMP ENDDO
          end if
!$OMP END PARALLEL
        end if
        call bcast_dbl(d(1), n, 1, 1, TRD_COMM_WORLD)

      else

        call barrier(TRD_COMM_WORLD)

        if (USE_MY_REDIST) then
#if TIMER_PRINT>1
          d1 = eigen_get_wtime()
#endif
          call dc_redist1(n, nvec, nb,
     &         z, work, lddz, iwork, liwork/2)
#if TIMER_PRINT>1
          d2 = eigen_get_wtime()
          if (iam == 0) print*,"MY-REDIST1",d2-d1
          d1 = eigen_get_wtime()
#endif
          call dc_redist2(n, nvec, nb,
     &         work, lddz, z, ldz, iwork, liwork/2)
#if TIMER_PRINT>1
          d2 = eigen_get_wtime()
          if (iam == 0) print*,"MY-REDIST2",d2-d1
#endif
        else
#if TIMER_PRINT>1
          d1 = eigen_get_wtime()
#endif
          call descinit(descw, n, n, 1, 1, 0, 0, ictxt, ldz, INFO_I4)
          call pdgemr2d(n, n, z, 1, 1, descz, work, 1, 1, descw,
     &         ictxt)
#if TIMER_PRINT>1
          d2 = eigen_get_wtime()
          if (iam == 0) print*,"PDGEMR2D",d2-d1
#endif
        end if                  ! if (USE_MY_REDIST) .. else ..

#if TIMER_PRINT>1
        d1 = eigen_get_wtime()
#endif
!$OMP PARALLEL
!$OMP+         PRIVATE(i,i_2,i_3,local_size,local_rank)
!$      local_size = omp_get_num_threads()
!$      local_rank = omp_get_thread_num()
!$      if (local_rank == 0) then
          call bcast_dbl(d(1), n, 1, 1, TRD_COMM_WORLD)
!$      end if
        if (.not.USE_MY_REDIST)then
          i_2=1
          i_3=nxv
!$        if (local_size > 1) then
!$        if (local_rank > 0) then
!$          i_2 = nxv * (local_rank-1)
!$          i_2 = i_2 / (local_size-1) + 1
!$          i_3 = nxv * (local_rank  )
!$          i_3 = i_3 / (local_size-1)
!$        else
!$          i_3 = 0
!$        end if
!$        end if
          do i=i_2, i_3
!OCL XFILL
          do j=1,ldz
            z(j,i)=work(j+(i-1)*ldz)
          end do
          end do
        end if
!$OMP END PARALLEL
#if TIMER_PRINT>1
        d2 = eigen_get_wtime()
        if (IAM==0) print*,"RERE1",d2-d1
#endif
      end if                    ! if (NB == 1) .. else ..

!     freeing working arrays
      deallocate(work)
      deallocate(iwork)

      comm_time_dc =
     &     eigen_timer_print('EigenExa (D-and-C for penta-diag)')


#if TIMER_PRINT>1
      if (iam == 0) then
        print*,"PDGEMM", dgemm_time, 
     &       1D-9*flops/dgemm_time,"GFLOPS",flops
      end if
#endif

      ret = flops

      return

      end subroutine  eigen_dcx

      end module dcx_mod
