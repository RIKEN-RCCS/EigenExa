
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_sx
!
! Purpose
! =======
!
! implementation of interface eigen_sx (see eigen_lib.F module)
!
!
! Arguments
! =========
!
! n, nvec, lda, ldz    (input) integer
!  lda and ldz are the leading dimensions corresponding to a and z,
!  respectively. These values must be equal or larger than one obtained 
!  by calling eigen_get_matdims.
!
! m_forward, m_backward  (input) integer, optional
!  default values are m_forward = 48 and m_backward = 128.
!     
! mode  (input) character(*), optional
!  = 'A' : eigenvalues and corresponding eigenvectors (default)
!  = 'N' : eigenvalues only
!  = 'X' : mode 'A' + improve accuracy of eigenvalues
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_sx(n, nvec, a, lda, w, z, ldz,
     &     m_forward, m_backward, mode)

      !
      ! Trick to avoid a conflict between func-name and interface
      !
      use eigen_libs_mod, eigen_sx_interface => eigen_sx
      use eigen_devel_mod
      use eigen_scaling_mod
      use eigen_prd_mod
      use bisect2_mod
      use dcx_mod
      use trbakwy4_mod
      use mpi, only : MPI_COMM_NULL
!$    use omp_lib
      use, intrinsic :: ieee_arithmetic
      use, intrinsic :: iso_c_binding

      implicit none

      integer,   intent(in)           :: n
      integer,   intent(in)           :: nvec
      real(8),   intent(inout)        :: a(lda,*)
      integer,   intent(in)           :: lda
      real(8),   intent(out)          :: w(1:n)
      real(8),   intent(out)          :: z(ldz,*)
      integer,   intent(in)           :: ldz
      integer,   intent(in), optional :: m_forward
      integer,   intent(in), optional :: m_backward
      character(*), intent(in), optional :: mode

      real(8), allocatable   :: d(:), e(:), e2(:)

      integer                :: m_f, m_b
      character              :: mode_

      real(8)                :: SIGMA
      integer                :: nm, ny, nb, nme
      integer                :: my_rank, world_size, ierr
      integer(8)             :: INFO
      integer                :: i
      logical                :: flag

      real(8)                :: hs0, hs1, s0, s1
      real(8)                :: ret, ret_1, ret_2, ret_3
!-
!-------------------------------------------------------------------
!-
      ERROR_INFO = 0
!-
!-------------------------------------------------------------------
!-
      call eigen_get_initialized(flag)
      if (.not. flag) then
        return
      end if
!-
!-------------------------------------------------------------------
!-
      if (TRD_COMM_WORLD == MPI_COMM_NULL) then
        return
      end if
!-
!-------------------------------------------------------------------
!-
      if (n <= 0) then
        print*,"Warining: Negative dimesion is invalid!"
        return
      end if
!-
!-------------------------------------------------------------------
!-
      if (present(mode)) then
        mode_ = mode(1:1)
      else
        mode_ = 'A'
      end if
      if ( nvec == 0 ) then
        mode_ = 'N'
      end if

      if (present(m_forward)) then
        m_f = m_forward
      else
        m_f = eigen_NB_f
      end if

      if (present(m_backward)) then
        m_b = m_backward
      else
        m_b = eigen_NB_b
      end if
!-
!-------------------------------------------------------------------
!-
      call eigen_get_matdims0(n, nm, ny, m_f, m_b)
      if (nm <= 0 .or. ny <= 0) then
        print*,"Warining: Problem size is too large for"
     &       //  " 32bit fortarn integer binary."
        return
      end if
!-
!-------------------------------------------------------------------
!-
      hs0 = eigen_get_wtime()
      ret_1 = ZERO; ret_2 = ZERO; ret_3 = ZERO
!-
!-------------------------------------------------------------------
!-
      nme = ((n-1)/4+1)*4
      allocate( d(1:n), e(1:nme*2), e2(1:nme*2), stat=ierr )
      if ( ierr /= 0 ) then
        call eigen_abort( "Memory allocation error [eigen_sx].", 1 )
      end if

      world_size = TRD_nnod
      my_rank    = TRD_inod-1
!-
!-------------------------------------------------------------------
!-
      SIGMA = ONE
      call eigen_scaling(n, a, lda, SIGMA)
      if ( ieee_is_nan( SIGMA ) ) then
         w(1:n) = SIGMA
         return
      end if
!-
!-------------------------------------------------------------------
!-
      if ( mode_ /= 'R' ) then
#if TIMER_PRINT
      s0 = eigen_get_wtime()
#endif

      call eigen_prd(n, a, lda, d, e, nme, m_f)
      ret_1 = dble(n)**3*4/3

#if TIMER_PRINT
      s1 = eigen_get_wtime()
      if (my_rank == 0) then
        print 10000,"TRD-BLK ",
     &  n, s1-s0, ret_1, 1D-9*ret_1/(s1-s0), "GFLOPS"
      end if
      call flush(6)
#endif
      else
      if (my_rank == 0) then
         print*,"Reading D.data"
         open(20,FILE='D.data',FORM="unformatted",STATUS="old")
         do i=1,n; read(20) D(i); enddo
         close(20)
         print*,"Reading E.data"
         open(21,FILE='E.data',FORM="unformatted",STATUS="old")
         do i=1,n; read(21) E(i); enddo
         close(21)
         print*,"Reading F.data"
         open(22,FILE='F.data',FORM="unformatted",STATUS="old")
         do i=1,n; read(22) E(i+nme); enddo
         close(22)
      endif
      call MPI_Bcast(D,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i)
      call MPI_Bcast(E,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i)
      call MPI_Bcast(E(1+nme),n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i)
      end if
!-
!-------------------------------------------------------------------
!-
#if TIMER_PRINT
      s0 = eigen_get_wtime()
#endif

      if (mode_ == 'A' .or. mode_ == 'X' .or.
     &    mode_ == 'T' .or. mode_ == 'R') then
        e2(0*nme+1:0*nme+N-1) = e(0*nme+2:0*nme+N)
        e2(0*nme+N) = 0
        e2(1*nme+1:1*nme+N-2) = e(1*nme+3:1*nme+N)
        e2(1*nme+N-1) = 0
        e2(1*nme+N) = 0
        w(1:n)=d(1:n)
        call eigen_dcx(n, ABS(nvec), w, e2, nme, z, ldz, INFO, ret_2)
        ERROR_INFO = INFO
        if (mode_ == 'X') then
          call eigen_bisect2(d, e, e(1+nme:), w, n, 1)
        end if
      end if
      if (mode_ == 'S' .or. mode_ == 'C') then
        call eigen_identity(n, z, ldz)
        call eigen_bisect2(d, e, e(1+nme:), w, n, 0)
      end if
      if (mode_ == 'N') then
        call eigen_bisect2(d, e, e(1+nme:), w, n, 0)
      end if

      deallocate( e2 )
      deallocate( d )

#if TIMER_PRINT
      s1 = eigen_get_wtime()
      if (my_rank == 0) then
        print 10000,"D&C     ",
     &  n, s1-s0, ret_2, 1D-9*ret_2/(s1-s0), "GFLOPS"
      end if
      call flush(6)
#endif
      if (mode_ == 'N') then
        goto 99999
      end if
!-
!-------------------------------------------------------------------
!-
      if ( mode_ /= 'T' .and. mode_ /= 'C' .and. mode_ /= 'R' ) then
#if TIMER_PRINT
      s0 = eigen_get_wtime()
#endif

      nb = MBAND
      call eigen_common_trbakwy(n, ABS(nvec), a, lda, z, ldz,
     &     e(1+nme:), m_b, nb)
      ret_3 = 2*dble(ABS(nvec))*dble(n)**2

#if TIMER_PRINT
      s1 = eigen_get_wtime()
      if (my_rank == 0) then
        print 10000,"TRDBAK  ",
     &  n, s1-s0, ret_3, 1D-9*ret_3/(s1-s0), "GFLOPS"
      end if
      call flush(6)
#endif
      end if
!-
!-------------------------------------------------------------------
!-
      if ( SIGMA /= ONE .and. SIGMA /= ZERO ) then
        SIGMA = ONE / SIGMA
        call dscal(n, SIGMA, w, 1)
      end if
!-
!-------------------------------------------------------------------
!-
#if TIMER_PRINT
      hs1 = eigen_get_wtime()
      ret = ret_1 + ret_2 + ret_3
      if (my_rank == 0) then
        print*,"Total(eigen_sx)",
     &  n, hs1-hs0, ret, 1D-9*ret/(hs1-hs0), "GFLOPS"
      end if
      call flush(6)
#endif
!-
!-------------------------------------------------------------------
!-
99999 continue
!-
!-------------------------------------------------------------------
!-
      hs1 = eigen_get_wtime()
      ret = ret_1 + ret_2 + ret_3
      if (ret_2 == 0) ret = -ret
      a(1, 1) = ret
      a(2, 1) = hs1-hs0
#if TIMER_PRINT
      a(3, 1) = comm_time_reduction
     &     + comm_time_dc
     &     + comm_time_backtrafo
#else
      a(3, 1) = -1d0
#endif
!-
!-------------------------------------------------------------------
!-
      deallocate( e )
!-
!-------------------------------------------------------------------
!-
10000 format (X, A8, I8, 3E25.16e2, X, A)
!-
!-------------------------------------------------------------------
!-
      end subroutine eigen_sx
