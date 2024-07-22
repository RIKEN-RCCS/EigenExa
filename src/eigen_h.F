
!--------*---------*---------*---------*---------*---------*---------*-*
!
! Subroutine: eigen_h
!
! Purpose
! =======
!
! implementation of interface eigen_h (see eigen_lib.F module)
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

         subroutine eigen_h( n, nvec, a, lda, w, z, ldz,
     $                       m_forward, m_backward, mode )
*-
         use eigen_libs_mod, eigen_h_interface => eigen_h
         use eigen_devel_mod
         use bisect_mod
         use dc2_mod
         use dc2_FS_mod
         use MPI, only : MPI_COMM_NULL
!$       use OMP_LIB
         use, intrinsic :: ieee_arithmetic
         use, intrinsic :: iso_c_binding

         implicit NONE
*-
         integer, intent(in)    :: n, nvec, lda, ldz
         complex(8), intent(inout) :: a(lda,*)
         complex(8), intent(out), target :: z(ldz,*)
         real(8), intent(out)      :: w(*)
         integer, intent(in), optional :: m_forward, m_backward
         character(*), intent(in), optional :: mode
*-
         real(8), allocatable     :: d(:), e(:), e2(:)
         complex(8), allocatable  :: beta(:)
*-
         integer                :: m_f, m_b
         character              :: mode_

         real(8)                :: SIGMA
         integer                :: nm, ny, nb
         integer                :: my_rank, world_size, ierr
         integer(8)             :: INFO
         logical                :: flag
*-
         real(8), pointer       :: z_r(:,:)
         logical                :: use_DC2FS

         real(8)                :: hs0, hs1, s0, s1
         real(8)                :: ret, ret_1, ret_2, ret_3

*-
*-------------------------------------------------------------------
*-
      call eigen_get_initialized(flag)
      if (.not. flag) then
        return
      end if
         if ( TRD_COMM_WORLD == MPI_COMM_NULL ) then
            return
         end if
*-
*-------------------------------------------------------------------
*-
      if (n <= 0) then
        print*,"Warining: Negative dimesion is invalid!"
        return
      end if
*-
*-------------------------------------------------------------------
*-
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
      m_f = min(m_f, n)
      m_f = max(m_f, 1)

      if (present(m_backward)) then
        m_b = m_backward
      else
        m_b = eigen_NB_b
      end if
      m_b = min(m_b, n)
      m_b = max(m_b, 1)
*-
*-------------------------------------------------------------------
*-
         call eigen_get_matdims0( n, nm, ny, m_f, m_b )
         if (nm <= 0 .or. ny <= 0) then
        print*,"Warining: Problem size is too large for"
     &       // " 32bit fortarn integer binary."
        return
      end if
*-
*-------------------------------------------------------------------
*-
         hs0 = eigen_get_wtime( )
         ret_1 = 0D0; ret_2 = 0D0; ret_3 = 0D0
*-
*-------------------------------------------------------------------
*-
         allocate( d(1:n), e(1:n), stat=ierr )
      if ( ierr /= 0 ) then
        call eigen_abort( "Memory allocation error [eigen_h].", 1 )
      end if
#if DEBUG
      if ( TRD_inod == 1 ) print*,"global data allocation done"
#endif

         world_size = TRD_nnod
         my_rank    = TRD_inod-1
*-
*-------------------------------------------------------------------
*-
         SIGMA = ONE
         call eigen_scaling_h( n, a, lda, SIGMA )
         if ( ieee_is_nan( SIGMA ) ) then
            w(1:n) = SIGMA
            return
         end if
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         s0 = eigen_get_wtime( )
#endif

         if ( mode_ == 'N' .or. mode_ == 'S' ) then
            call eigen_hrd( n, a, lda, d, e, m_f )
         else
            call eigen_hrd( n, a, lda, w, e, m_f )
         end if
         ret_1 = dble(n)**3*4/3
        
#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if(my_rank==0)then
            print*,"HRD-BLK", n, s1-s0, 1D-9*ret_1/(s1-s0), "GFLOPS"
         endif
         call flush(6)
#endif
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         s0 = eigen_get_wtime( )
#endif

         use_DC2FS = (n > TRD_nnod) .and. (TRD_nnod >= 4)
         call C_F_Pointer( C_LOC(z), z_r, shape=[ldz,ny] )

         if ( mode_ == 'X' ) then
           allocate( e2(1:n), stat=ierr )
           if ( ierr /= 0 ) then
             call eigen_abort( "Memory allocation error [eigen_s].", 1 )
             return
           end if
           d (1:n) = w(1:n)
           e2(1:n) = e(1:n)
           if ( use_DC2FS ) then
           call dc2_FS(n, ABS(nvec), w, e2(2:), z_r, ldz, INFO,ret_2)
           else
           call eigen_dc2(n, ABS(nvec), w, e2(2:), z_r, ldz, INFO,ret_2)
           end if
           deallocate( e2 )
           call eigen_bisect(d, e, w, n, 1)
           call convert_DtoZ( ny, z_r, ldz, z, ldz, d, n )
         end if
         if ( mode_ == 'A' ) then
           d(1:n) = e(1:n)
           if ( use_DC2FS ) then
           call dc2_FS( n, ABS(nvec), w, d(2:), z_r, ldz, INFO,ret_2 )
           else
           call eigen_dc2(n, ABS(nvec), w, d(2:), z_r, ldz, INFO,ret_2)
           end if
           call convert_DtoZ( ny, z_r, ldz, z, ldz, d, n )
         end if
         if ( mode_ == 'S' ) then
           call eigen_identity(n, z_r, ldz)
           call eigen_bisect( d, e, w, n, 0 )
           call convert_DtoZ( ny, z_r, ldz, z, ldz, d, n )
         end if
         if ( mode_ == 'N' ) then
           call eigen_bisect( d, e, w, n, 0 )
         end if
         deallocate( d )

         z_r=>null()

#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if ( my_rank == 0 ) then
            print*,"D&C", n, s1-s0, 1D-9*ret_2/(s1-s0), 'GFLOPS'
         end if
         call flush(6)
#endif
         if ( mode_ == 'N' ) then
            goto 99999
         end if

*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         s0 = eigen_get_wtime( )
#endif
         allocate( beta(1:n), stat=ierr )
         if ( ierr /= 0 ) then
             call eigen_abort( "Memory allocation error [eigen_s].", 1 )
             return
         end if
         beta(1:n) = DCMPLX(e(1:n))
         deallocate( e )

         nb = 1
         call eigen_hrbakwyx( n, ABS(nvec), a, lda, z, ldz, beta,
     &        m_b, nb )
         ret_3 = 2*dble(ABS(nvec))*dble(n)**2
            
#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if ( my_rank == 0 ) then
            print*,"HRDBAK", n, s1-s0, 1D-9*ret_3/(s1-s0), "GFLOPS"
         end if
         call flush(6)
#endif
*-
*-------------------------------------------------------------------
*-
         IF( SIGMA /= ONE .and. SIGMA /= ZERO ) THEN
            SIGMA = ONE / SIGMA
            call dscal( n, SIGMA, w, 1 )
         END IF
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         hs1 = eigen_get_wtime( )
         ret = ret_1 + ret_2 + ret_3
         if ( my_rank == 0 ) then
            print*,"Total", hs1-hs0, 1D-9*ret/(hs1-hs0), "GFLOPS"
         end if
         call flush(6)
#endif
*-
*-------------------------------------------------------------------
*-
99999    continue
*-
*-------------------------------------------------------------------
*-
         hs1 = eigen_get_wtime( )
         ret = ret_1 + ret_2 + ret_3
         a(1, 1) = DCMPLX(ret, ZERO)
         a(2, 1) = DCMPLX(hs1-hs0, ZERO)
*-
*-------------------------------------------------------------------
*-
         deallocate( beta )
*-
*-------------------------------------------------------------------
*-
         contains

         subroutine convert_DtoZ( n, a, lda, z, ldz, buf, lbuf )
         integer,    intent(in)  :: n, lda, ldz, lbuf
         real(8),    intent(in)  :: a(1:lda,*)
         complex(8), intent(out) :: z(1:ldz,*)
         real(8),    intent(out) :: buf(1:lbuf)

         integer :: nm, i, j, k, l

         nm = min(lda,ldz,lbuf)
!$OMP PARALLEL PRIVATE(i,j)
         do i=n,1,-1
!$OMP DO
           do j=1,nm
              buf(j) = a(j,i)
           end do
!$OMP END DO
!$OMP DO
           do j=1,nm
               z(j,i) = DCMPLX(buf(j), ZERO)
           end do
!$OMP END DO
         end do
!$OMP END PARALLEL

         return
         end subroutine  convert_DtoZ

         end subroutine  eigen_h

