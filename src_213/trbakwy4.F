!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: trbakwy4_mod
!
! Purpose
! =======
!
! trbakwy4_mod manages the modules for back-transformation routines.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module trbakwy4_mod

      use trbakwy4_body_mod
      use comm_mod
      use eigen_libs_mod
      use eigen_devel_mod
      use eigen_trbak_mod
      use CSTAB_mod
      use mpi
!$    use omp_lib

      implicit none
      private
      

      public  :: eigen_common_trbakwy
      private :: eigen_trbakwy_body

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_common_trbakwy
!
! Purpose
! =======
!
! eigen_common_trbakwy is the main body of the Householder back
! transformation routine by taking use of the WY representation
! algorithm.
!
!
! Arguments
! =========
!
! n       (input) integer
!         The dimension of the input matrix. N >= 0.
!
! nvec    (input) integer
!         The number of vectors to be transformed. N >= NVEC >= 0.
!         If nvec is zero, the function returns immediately.
!
! a       (input) real(8) array, dimension(*)
!         a contains the data of the householder reflectors.
!
! nma0    (input) integer
!         The leading dimension of the array a.
!
! z       (output) real(8) array, dimension(*)
!         z contains the elements of the eigenvectors.
!
! nmz0    (input) integer
!         The leading dimension of the array z.
!
! beta    (in/output) real(8) array, dimension(*)
!         beta contains the data of the householder transformation.
!
! m0      (input) integer
!         Block factor
!
! iblk    (input) integer
!         Auxiliary data for back-transformation (tri-/penta-/...).
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_common_trbakwy(n, nvec, a, nma0, z, nmz0,
     &     beta, m0, iblk)

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      real(8), intent(in)    :: a(*)
      integer, intent(in)    :: nma0
      real(8), intent(out)   :: z(*)
      integer, intent(in)    :: nmz0
      real(8), intent(inout) :: beta(*)
      integer, intent(in)    :: m0
      integer, intent(in)    :: iblk

      real(8),  allocatable  :: v1(:), ss1(:)
      real(8),  allocatable  :: v2(:), ss2(:)
      real(8),  allocatable  :: v3(:), ss3(:)
      real(8),  allocatable  :: tt(:)

      include 'CSTAB.h'

      integer                :: nma, nmz, m
      integer                :: na, nm, kx
      integer                :: i_v1, i_v2, i_v3
      integer                :: i_s1, i_s2, i_s3
      integer                :: i_t
      integer                :: i_stat, ierr


      if ( nvec == 0 ) return

      nma = nma0
      nmz = nmz0
      m   = min(nsm, m0)
      if (m < 1) m = 1

      call eigen_timer_reset(2, 4, 0, 1)

      na   = (n-1)/y_nnod+1
      na   = na + mod(na - 1,2)
      call CSTAB_get_optdim(nma,9,16*4,16*6,nm)

#if BOOST_BY_CACHE_ALIGNMENT
      allocate(
     &     v1(1:max(nm*m,n)+n_columns),
     &     ss1(1:na*m+ns0+n_columns),
     &     v2(1:max(nm*m,n)+n_columns),
     &     ss2(1:na*m+ns0+n_columns),
     &     v3(1:max(nm*m,n)+n_columns),
     &     ss3(1:na*m+ns0+n_columns),
     &     tt(1:na*m+ns0+n_columns),
     &     stat=i_stat)
#else
      allocate(
     &     v1(1:max(nm*m,n)),
     &     ss1(1:na*m+ns0),
     &     v2(1:max(nm*m,n)),
     &     ss2(1:na*m+ns0),
     &     v3(1:max(nm*m,n)),
     &     ss3(1:na*m+ns0),
     &     tt(1:na*m+ns0),
     &     stat=i_stat)
#endif
      if (i_stat /= 0) then
        call eigen_abort( "Memory allocation error [TBK].", 1 )
      end if

      v1 = ZERO
      ss1 = ZERO
      v2 = ZERO
      ss2 = ZERO
      v3 = ZERO
      ss3 = ZERO
      tt = ZERO

#if BOOST_BY_CACHE_ALIGNMENT
      call CSTAB_adjust_base(v1,  z, i_v1)
      call CSTAB_adjust_base(ss1, z, i_s1)
      call CSTAB_adjust_base(v2,  z, i_v2)
      call CSTAB_adjust_base(ss2, z, i_s2)
      call CSTAB_adjust_base(v3,  z, i_v3)
      call CSTAB_adjust_base(ss3, z, i_s3)
      call CSTAB_adjust_base(tt,  z, i_t)
      kx = (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &     +(L1_LSIZE)
     &     +(L2_LSIZE/8)
      i_v1 = i_v1 + kx*5
      i_s1 = i_s1 + kx*1
      i_v2 = i_v2 + kx*5
      i_s2 = i_s2 + kx*1
      i_v3 = i_v3 + kx*5
      i_s3 = i_s3 + kx*1
      i_t = i_t + kx*1
      call CSTAB_round_offset(i_v1)
      call CSTAB_round_offset(i_s1)
      call CSTAB_round_offset(i_v2)
      call CSTAB_round_offset(i_s2)
      call CSTAB_round_offset(i_v3)
      call CSTAB_round_offset(i_s3)
      call CSTAB_round_offset(i_t)
#else
      i_v1=0
      i_s1=0
      i_v2=0
      i_s2=0
      i_v3=0
      i_s3=0
      i_t=0
#endif

      allocate(TRBK_buf(1:m*nm), stat=i_stat)
      if (i_stat /= 0) then
        call eigen_abort( "Memory allocation error [TBK].", 1 )
      end if
     
      TRBK_TIME_REDUC_OVERHEAD_x = Reduce_cont_Overhead_x

      call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP PARALLEL
      call eigen_trbakwy_body( n, nvec,
     &     a, nma,
     &     z, nmz,
     &     beta,
     &     v1(1+i_v1:), v2(1+i_v2:), v3(1+i_v3:), nm, m,
     &     ss1(1+i_s1:), ss2(1+i_s2:), ss3(1+i_s3:),
     &     tt(1+i_t:), iblk, na )
!$OMP END PARALLEL
      call MPI_Barrier(TRD_COMM_WORLD,ierr)

      deallocate(TRBK_buf)

      deallocate(v1)
      deallocate(ss1)
      deallocate(v2)
      deallocate(ss2)
      deallocate(v3)
      deallocate(ss3)
      deallocate(tt)

      comm_time_backtrafo =
     &     eigen_timer_print('EigenExa (Back-transformation)')

      return

      end subroutine eigen_common_trbakwy

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine eigen_trbakwy_body(n, nvec, a, nma, z, nmz, beta,
     &     v1, v2, v3, nm, m, ss1, ss2, ss3, tt, iblk, nss)

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      real(8), intent(in)    :: a(1:nma,*)
      integer, intent(in)    :: nma
      real(8), intent(inout) :: z(1:nmz,*)
      integer, intent(in)    :: nmz
      real(8), intent(inout) :: beta(1:n)
      real(8), intent(out)   :: v1(1:nm,*)
      real(8), intent(out)   :: v2(1:nm,*)
      real(8), intent(out)   :: v3(1:nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: m
      real(8), intent(out)   :: ss1(*)
      real(8), intent(out)   :: ss2(*)
      real(8), intent(out)   :: ss3(*)
      real(8), intent(out)   :: tt(*)
      integer, intent(in)    :: iblk
      integer, intent(in)    :: nss

      real(8), allocatable   :: wk(:)
      integer                :: nodes(0:nsm-1)
      integer                :: nx, ierr
      integer(8)             :: lwk
      integer                :: i, L, mode
      integer                :: x_root
      integer                :: iloop_sta, iloop_end
      integer                :: jloop_sta, jloop_end
      integer                :: i_1, i_4, j_1
      real(8)                :: s0, s1, s2, s3
      real(8)                :: d0, d1, d2, ds, de, dcom, dx, dy, dz


#if TIMER_PRINT
      dx = ZERO; dy = ZERO; dz = ZERO
#endif

!$OMP MASTER
      call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
#if TIMER_PRINT
      d1 = eigen_get_wtime()
#endif
      dcom = ZERO
!-
#if AT_BCAST_OVERLAP
!$    if(.TRUE.) DO_OVERLAP_BCAST_LEVEL  = 2
!$    TRBK_TIME_COUNTER  = 0
!$    TRBK_TIME_INTERVAL = 1
!$    TRBK_TIME_NEXT     = 2
!$    TRBK_SWITCHED      = 0
!$    TRBK_TIME_TRBK1    = ZERO
!$    TRBK_TIME_TRBK1_   = ZERO
!$    TRBK_TIME_TRBK1x   = ZERO
!$    TRBK_TIME_TRBK1x_  = ZERO
!$    TRBK_TIME_TRBK2    = ZERO
!$    TRBK_TIME_BCAST    = ZERO
!$    TRBK_TIME_REDUC    = ONE
!$    TRBK_TIME_FR       = ZERO
!===  
!$    call sync_other_than_master_init(TRBK_LOCK, TRBK_MASK)
!===  
#else
!$    if(.true.) DO_OVERLAP_BCAST_LEVEL  = 0
#endif
!-
!$OMP END MASTER

      nx = min(mod(n-(1+iblk)+1,m)+(1+iblk)-1,n)
!$OMP MASTER
      lwk=int((m-1)/y_nnod+1,8)*int((n-1)/x_nnod+1, 8)
      lwk=max(lwk, int(n,8))
      allocate(wk(lwk), stat=ierr)
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error [TBK].", 1 )
      end if
!$OMP END MASTER

!$OMP DO
      do i=1,nx
        l = i-iblk
        if ( i >= 1+iblk .and.
     &       mod(i-1, y_nnod)+1 == y_inod .and.
     &       mod(l-1, x_nnod)+1 == x_inod ) then
          i_1 = (i-y_inod)/y_nnod+1
          j_1 = (l-x_inod)/x_nnod+1
          beta(i) = a(j_1,i_1)*beta(i)
        else
          beta(i) = ZERO
        end if
      end do
!$OMP ENDDO
!$OMP MASTER
      call reduce_dbl(beta, wk, nx, 1, TRD_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
      do i=1+iblk,nx
        if (beta(i) == ZERO) then
          s0 = ONE
        else
          s0 = ONE/beta(i)
        end if
        beta(i) = s0
      end do                    ! i
!$OMP ENDDO
!
!----------------------------------------------------------------

      iloop_sta = eigen_loop_start(1,    'Y')
      iloop_end = eigen_loop_end  (nvec, 'Y')

      nx = min(mod(n-(1+iblk)+1,m)+(1+iblk)-1,n)

      do i=(1+iblk),nx

!$OMP BARRIER

        if (beta(i) == ZERO) cycle

        jloop_sta = eigen_loop_start(1,      'X')
        jloop_end = eigen_loop_end  (i-iblk, 'X')
        i_4=mod(iloop_end-iloop_sta+1,4)+iloop_sta

!$OMP MASTER
        ds = eigen_get_wtime()
        nodes(0) = eigen_owner_node(i, 'Y')
        if (nodes(0) == y_inod) then
          i_1 = eigen_translate_g2l(i, 'Y')
          do j_1=jloop_sta,jloop_end
            v1(j_1,1) = a(j_1,i_1)
          end do                ! j_1
        end if
        call bcast_dbl(v1(jloop_sta,1), jloop_end-jloop_sta+1,
     &       nodes(0), 1, y_COMM_WORLD)

        de = eigen_get_wtime()
        dcom = dcom + (de - ds)
!$OMP END MASTER

!$OMP BARRIER


#define	USE_BLAS	1

!$OMP MASTER
        call barrier(y_COMM_WORLD)
        call barrier(x_COMM_WORLD)
!$OMP END MASTER
        if (jloop_end >= jloop_sta) then
          i_4=mod(iloop_end-iloop_sta+1,8)+iloop_sta
          if (i_4 > iloop_sta) then
!$OMP DO SCHEDULE(STATIC)
            do i_1=iloop_sta,i_4-1
#if USE_BLAS
              ss1(i_1+0)=ddot(jloop_end-jloop_sta+1, z(jloop_sta,i_1),1,
     &             v1(jloop_sta,1),1)
#else
              s0 = ZERO
              do j_1=jloop_sta,jloop_end
                s0=s0+v1(j_1,1)*z(j_1,i_1+0)
              end do            ! j_1
              ss1(i_1+0)=s0
#endif
            end do
!$OMP END DO
          end if
!$OMP DO SCHEDULE(STATIC)
#if USE_BLAS
          do i_1=i_4,iloop_end,8
            call dgemv('T', jloop_end-jloop_sta+1, 8,
     &           ONE, z(jloop_sta,i_1), nmz, v1(jloop_sta,1), 1,
     &           ZERO, ss1(i_1), 1)
          end do                ! i_1
#else
          do i_1=i_4,iloop_end,4
            s0 = ZERO
            s1 = ZERO
            s2 = ZERO
            s3 = ZERO
            do j_1=jloop_sta,jloop_end
              s0=s0+v1(j_1,1)*z(j_1,i_1+0)
              s1=s1+v1(j_1,1)*z(j_1,i_1+1)
              s2=s2+v1(j_1,1)*z(j_1,i_1+2)
              s3=s3+v1(j_1,1)*z(j_1,i_1+3)
            end do              ! j_1
            ss1(i_1+0)=s0
            ss1(i_1+1)=s1
            ss1(i_1+2)=s2
            ss1(i_1+3)=s3
          end do                ! i_1
#endif
!$OMP END DO
        else
!$OMP DO
          do i_1=iloop_sta,iloop_end
            ss1(i_1) = ZERO
          end do
!$OMP END DO
        end if

!$OMP BARRIER

!$OMP MASTER
        ds = eigen_get_wtime()
        call reduce_dbl(ss1(iloop_sta),tt, iloop_end-iloop_sta+1, 2, 
     &       x_COMM_WORLD)
        de = eigen_get_wtime()
        dcom = dcom + (de - ds)

#if USE_BLAS
        call dscal(iloop_end-iloop_sta+1, beta(i), ss1(iloop_sta), 1)
#else
        s0 = beta(i)
        do i_1=iloop_sta,iloop_end
          ss1(i_1) = ss1(i_1) * s0
        end do                  ! i_1
#endif
!$OMP END MASTER

!$OMP BARRIER

        jloop_sta = eigen_loop_start(1,      'X')
        jloop_end = eigen_loop_end  (i-iblk, 'X')
        i_4=mod(iloop_end-iloop_sta+1,8)+iloop_sta

        if (jloop_end >= jloop_sta) then
          if (i_4 > iloop_sta) then
!$OMP DO SCHEDULE(STATIC)
            do i_1=iloop_sta, i_4-1
#if USE_BLAS
              call daxpy(jloop_end-jloop_sta+1,
     &             ss1(i_1), v1(jloop_sta,1), 1, z(jloop_sta,i_1), 1)
#else
              s0 = ss1(i_1+0)
              do j_1=jloop_sta,jloop_end
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
              end do            ! j_1
#endif
            end do
!$OMP END DO
          end if
!$OMP DO SCHEDULE(STATIC)
#if USE_BLAS
          do i_1=i_4,iloop_end,8
            call dger(jloop_end-jloop_sta+1, 8,
     &           ONE, v1(jloop_sta,1), 1, 
     &           ss1(i_1), 1, z(jloop_sta,i_1), nmz)
          end do                ! i_1
#else
          do i_1=i_4,iloop_end,4
            s0 = ss1(i_1+0)
            s1 = ss1(i_1+1)
            s2 = ss1(i_1+2)
            s3 = ss1(i_1+3)
            do j_1=jloop_sta,jloop_end
              z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
              z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
              z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v1(j_1,1)
              z(j_1,i_1+3) = z(j_1,i_1+3) + s3 * v1(j_1,1)
            end do              ! j_1
          end do                ! i_1
#endif
!$OMP ENDDO
        end if

!$OMP BARRIER

      end do

!----------------------------------------------------------------

      d2 = eigen_get_wtime()
      d0 = d2 -d1

#define	USE_BCASTW	1

      mode = 0; i = nx+1

!$OMP MASTER
      if (nx+1 <= n) then
        ds = eigen_get_wtime()
        call trbakwy_datacast(iloop_end, m, i, a(1,1), nma,
     &       v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     &       wk(1),
#endif
     &       iblk)
        de = eigen_get_wtime()
        dcom = dcom + (de - ds)
      end if
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
      if (nx+1+m <= n) then
        ds = eigen_get_wtime()
        call trbakwy_datacast(iloop_end, m, i+m, a(1,1), nma,
     &       v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     &       wk(1),
#endif
     &       iblk)
        de = eigen_get_wtime()
        dcom = dcom + (de - ds)
      end if
!$OMP END MASTER

      do i=nx+1, n, m

        if (mode == 0) then
          call eigen_trbakwy_block_body(iloop_end, z(1,1), nmz,
#if 0
     &         d(1), e(1),
#endif
     &         v1(1,1), nm, m,
     &         i, ss1(1), tt, nss, iblk,
     &         dcom, dx, dy, dz)
        end if
        if (mode == 1) then
          call eigen_trbakwy_block_body(iloop_end, z(1,1), nmz,
#if 0
     &         d(1), e(1),
#endif
     &         v2(1,1), nm, m,
     &         i, ss2(1), tt, nss, iblk,
     &         dcom, dx, dy, dz)
        end if
        if (mode == 2) then
          call eigen_trbakwy_block_body(iloop_end, z(1,1), nmz,
#if 0
     &         d(1), e(1),
#endif
     &         v3(1,1), nm, m,
     &         i, ss3(1), tt, nss, iblk,
     &         dcom, dx, dy, dz)
        end if

!$OMP MASTER
        if (i+2*m <= n) then
          ds = eigen_get_wtime()
          if (mode == 0) then
            call trbakwy_datacast(iloop_end, m, i+2*m, a(1,1), nma,
     &           v3(1,1), nm, ss3(1),
#if !USE_BCASTW
     &           wk(1),
#endif
     &           iblk)
          end if
          if (mode == 1) then
            call trbakwy_datacast(iloop_end, m, i+2*m, a(1,1), nma,
     &           v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     &           wk(1),
#endif
     &           iblk)
          end if
          if (mode == 2) then
            call trbakwy_datacast(iloop_end, m, i+2*m, a(1,1), nma,
     &           v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     &           wk(1),
#endif
     &           iblk)
          end if
          de = eigen_get_wtime()
          dcom = dcom + (de - ds)
        end if
!$OMP END MASTER

        mode = mod(mode + 1, 3)

      end do

!----------------------------------------------------------------

!$    if (omp_get_num_threads() > 1) then
!$OMP BARRIER
!$       if (omp_get_thread_num() == 1) then
!$          ss1(1) = dx; ss1(2) = dy
!$       end if
!$OMP BARRIER
!$       if (omp_get_thread_num() == 0) then
!$          dx = ss1(1); dy = ss1(2)
!$       end if
!$    end if
!$OMP BARRIER
!     
#if AT_BCAST_OVERLAP
!===  
!$OMP MASTER
!$    call sync_other_than_master_finalize(TRBK_LOCK)
!$OMP END MASTER
!===  
#endif
!$OMP MASTER
      call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
      deallocate(wk)

#if TIMER_PRINT>1
      d2 = eigen_get_wtime()

      if (1 == TRD_inod) then
        print*,"TRBAK=",(d2-d1)
        print*,"COMM=",dcom
        print*,"   ",d2-d1,(2d0*nvec*n*n)/(d2-d1)*1d-9,"GFLOPS"
        if (dx > 0)
     &       print*,"   ",dx,(1d0*nvec*n*n)/(dx)*1d-9,"GFLOPS"
        if (dy > 0)
     &       print*,"   ",dy,(1d0*nvec*n*n)/(dy)*1d-9,"GFLOPS"
        print*,"   ",d0,dz
      end if
#endif
!$OMP END MASTER

      return

      end subroutine eigen_trbakwy_body

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine trbakwy_datacast(iloop_end, m, i, a, nma, v, nm,
     &     ss,
#if !USE_BCASTW
     &     wk,
#endif
     &     iblk)

      integer, intent(in)   :: iloop_end
      integer, intent(in)   :: m
      integer, intent(in)   :: i
      real(8), intent(in)   :: a(nma,*)
      integer, intent(in)   :: nma
      real(8), intent(out)  :: v(nm,*)
      integer, intent(in)   :: nm
      real(8), intent(out)  :: ss(*)
#if !USE_BCASTW
! size of wk should be larger than ((m-1)/y_nnod+1)*((n-1)/x_nnod+1)
      real(8), intent(out)  :: wk(*)
#endif
      integer, intent(in)   :: iblk

      integer               :: j, iy, k0
      integer               :: i_1, j_1, j_4, j_5
      integer               :: jloop_sta, jloop_end
      integer               :: nodes(0:nsm-1)

!$    real(8)               :: ds, de


!$    ds = omp_get_wtime()

      ss(1:iloop_end*m+ns0) = ZERO

#if USE_BCASTW
      jloop_end = min(eigen_loop_end(i+m-1-iblk, 'X'), nm)

      if (mod(m, y_nnod) == 0) then

        do j=0,m-1
          if (y_inod == eigen_owner_node(i+j, 'Y')) then
            iy  = eigen_owner_node(i+j, 'Y')
            i_1 = eigen_translate_g2l(i+j, 'Y', iy)
            k0  = (j/y_nnod)*jloop_end
            v(1+k0:jloop_end+k0, 1) = a(1:jloop_end, i_1)
          end if
        end do

        k0  = (m/y_nnod)*jloop_end
        call allgather_dbl(v, TRBK_buf, k0, 1, y_COMM_WORLD)

        do j=0,m-1
          iy = eigen_owner_node(i+j, 'Y')
          k0 = ((iy-1)*(m/y_nnod)+(j/y_nnod))*jloop_end
          v(1:jloop_end, j+1) = TRBK_buf(1+k0:jloop_end+k0)
          jloop_sta = eigen_loop_start(i+j, 'X')
          if (jloop_sta <= jloop_end)
     &         v(jloop_sta:jloop_end, j+1) = ZERO
        end do

      else

        do j=0,m-1
          if (y_inod == eigen_owner_node(i+j, 'Y')) then
            iy  = eigen_owner_node(i+j, 'Y')
            i_1 = eigen_translate_g2l(i+j, 'Y', iy)
            v(1:jloop_end, j+1) = a(1:jloop_end, i_1)
            jloop_sta = eigen_loop_start(i+j, 'X')
            if (jloop_sta <= jloop_end)
     &           v(jloop_sta:jloop_end, j+1) = ZERO
          end if
        end do
        do j=0,m-1,y_nnod
          iy  = eigen_owner_node(i+j, 'Y')
          k0  = min(m-j, y_nnod)
          call bcastw_dbl(v(1,j+1), jloop_end, iy,
     &         nm, k0, TRBK_buf, 2, y_COMM_WORLD)
        end do

      end if
#else

      if (m > y_nnod .and. y_nnod > 1) then

        do j=0,m-1
          nodes(j) = eigen_owner_node(i+j, 'Y')
        end do

        do iy=1,y_nnod

          jloop_sta = eigen_loop_start(1,         'X')
          jloop_end = eigen_loop_end  (i+m-1-iblk,'X')

          k0=0
          do j=0,m-1
            if (nodes(j) == iy) then
              i_1 = eigen_translate_g2l(i+j, 'Y', iy)
              if (y_inod == iy) then
                do j_1=jloop_sta,jloop_end
                  wk(k0+j_1) = a(j_1, i_1)
                end do          ! k
              end if
              k0=k0+(jloop_end-jloop_sta+1)
            end if
          end do

          call bcast_dbl(wk, k0, iy, 2, y_COMM_WORLD)

          k0=0
          do j=0,m-1
            if (nodes(j) == iy) then
              do j_1=jloop_sta,jloop_end
                v(j_1, j+1) = wk(k0+j_1)
              end do            ! k
              k0=k0+(jloop_end-jloop_sta+1)
            end if
          end do

          do j=0,m-1
            if (nodes(j) == iy) then
              j_4 = eigen_loop_start(i+j,       'X')
              j_5 = eigen_loop_end  (i+m-1-iblk,'X')
              do j_1=j_4,j_5
                v(j_1, j+1) = ZERO
              end do
            end if
          end do

        end do

      else

        do j=0,m-1
          nodes(j) = eigen_owner_node(i+j, 'Y')
          if (nodes(j) == y_inod) then
            i_1 = eigen_translate_g2l(i+j, 'Y')
            jloop_sta = eigen_loop_start(1,         'X')
            jloop_end = eigen_loop_end  (i+m-1-iblk,'X')
            do j_1=jloop_sta,jloop_end
              v(j_1, j+1) = a(j_1, i_1)
            end do              ! k
            jloop_sta = eigen_loop_start(i+j,       'X')
            jloop_end = eigen_loop_end  (i+m-1-iblk,'X')
!     CDIR NOVECTOR
            do j_1=jloop_sta,jloop_end
              v(j_1, j+1) = ZERO
            end do
          end if
        end do
        jloop_sta = eigen_loop_start(1,          'X')
        jloop_end = eigen_loop_end  (i+m-1-iblk, 'X')
        do j=0,m-1
          call bcast_dbl(v(1,j+1), jloop_end-jloop_sta+1,
     &         nodes(j), 2, y_COMM_WORLD)
        end do

      end if

#endif

!$    de = omp_get_wtime()
!$    TRBK_TIME_BCAST = de - ds

      return

      end subroutine

      end module trbakwy4_mod
