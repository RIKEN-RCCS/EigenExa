!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_prd_mod
!
! Purpose
! =======
!
! eigen_prd_mod manages the modules of the penta-diagonalization of
! eigen_sx.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_prd_mod

      use eigen_prd_t8_mod
      use eigen_prd_t7_mod
      use eigen_prd_t6_3_mod
      use eigen_prd_t5_mod
      use eigen_prd_t4x_mod
      use eigen_prd_t2_mod
      use eigen_t1_mod
      use comm_mod
      use eigen_devel_mod
      use CSTAB_mod
      use mpi
!$    use omp_lib

      implicit none
      private


      public  :: eigen_prd
      private :: eigen_prd_stub
      private :: eigen_prd_body

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd
!
! Purpose
! =======
!
! eigen_prd invokes the stub function of the penta-diagonalization,
! eigen_prd_stub, to transform a dense matrix to a penta-diagonal form.
!
!
! Arguments
! =========
!
! n      (input) integer
!        The dimension of the input matrix. N >= 0.
!
! a      (input/output) real(8) array, dimension(*)
!        On entry, a contains the input real symmetric matrix to be
!        transformed.
!        On exit, a contains the data of the hoseholder reflectors.
!
! nma0   (input) integer
!        The leading dimension of the input array a.
!        nma0 >= ceil(N/x_nnod)
!
! d_out  (output) real(8) array, dimension(*)
!        d_out contains diagonal elements of the output penta-diagonal
!        matrix.
!
! e_out  (output) real(8) array, dimension(*)
!        e_out contains off-diagonal elements of the output
!        penta-diagonal matrix.
!
! nme0   (input) integer
!        The leading dimension of the output array e_out. nme0 >= N.
!
! m_orig (input) integer
!        Block factor.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_prd(n, a, nma0,
     &     d_out, e_out, nme0, m_orig)

      integer, intent(in)    :: n
      real(8), intent(inout) :: a(*)
      integer, intent(in)    :: nma0
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(*)
      integer, intent(in)    :: nme0
      integer, intent(in)    :: m_orig

      include 'CSTAB.h'

      integer                :: nm, ne, m, ierr, lda


      nm = nma0
      ne = nme0
      m  = m_orig

      call MPI_Barrier(TRD_COMM_WORLD, ierr)
      call eigen_timer_reset(2, 6, 2, 2)

      lda = nma0

      call eigen_prd_stub(a, lda, d_out, e_out, ne, n, m)

      call flush(6)

      call MPI_Barrier(TRD_COMM_WORLD, ierr)
      comm_time_reduction =
     &     eigen_timer_print('EigenExa(Penta-diagonalization)')

      return

      end subroutine eigen_prd

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_stub
!
! Purpose
! =======
!
! eigen_prd_stub invokes the main body of the penta-diagonalization,
! eigen_prd_body, to transform a dense matrix to a penta-diagonal form.
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(nm,*)
!        On entry, a contains the input real symmetric matrix to be
!        transformed.
!        On exit, a contains the data of the hoseholder reflectors.
!
! nm     (input) integer
!        The leading dimension of the input array a.
!        nm >= ceil(N/x_nnod)
!
! d_out  (output) real(8) array, dimension(n)
!        d_out contains diagonal elements of the output penta-diagonal
!        matrix.
!
! e_out  (output) real(8) array, dimension(ne,*)
!        e_out contains off-diagonal elements of the output
!        penta-diagonal matrix.
!
! ne     (input) integer
!        The leading dimension of the output array e_out. ne >= N.
!
! n      (input) integer
!        The dimension of the input matrix. N >= 0.
!
! m_orig (input) integer
!        Block factor.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_prd_stub(a, nm, d_out, e_out, ne, n, m_orig)

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(out)   :: d_out(1:n)
      real(8), intent(out)   :: e_out(1:ne,*)
      integer, intent(in)    :: ne
      integer, intent(in)    :: n
      integer, intent(in)    :: m_orig

      include 'CSTAB.h'

      integer, parameter     :: nm_max_L1 = 16*4
      integer, parameter     :: nm_max_L2 = 16*6

      real(8), allocatable   :: u_t(:), v_t(:), d_t(:)
      real(8), allocatable   :: w(:)
      real(8), allocatable   :: ux_vx(:)
      real(8), allocatable   :: uy_vy(:)

      integer                :: nv, nx, ny, nz
      integer                :: ierr, kx
      integer                :: offset1, offset2, offset3
      integer                :: offset4, offset5, offset6
      integer                :: offset7

      integer(8) :: l_array_1, l_array_2


      nx = (n-1)/x_nnod+1 +2
      nv = nm

      call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv)

      nx = (n-1)/x_nnod+1
      kx = y_nnod / n_common
      nz = ((nx-1)/kx+1)
      ny = nz*MAX(x_nnod,y_nnod)
!
      l_array_1 = MAX(MAX(m_orig,2)*ny, nx+4*m_orig+6)
      l_array_2 = MAX(MAX(m_orig,2)*nz, nx+4*m_orig+6)

      allocate(
     &     w(1:nm*m_orig+n_columns),
     &     u_t(1:l_array_1+n_columns),
     &     v_t(1:l_array_2+n_columns),
     &     d_t(1:nv+n_columns),
     &     ux_vx(1:nv*2*m_orig+2*n_columns),
     &     uy_vy(1:nv*2*m_orig+2*n_columns),
     &     stat=ierr)
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error.", 1 )
      end if

      w = ZERO
      v_t = ZERO
      u_t = ZERO
      d_t = ZERO
      ux_vx = ZERO
      uy_vy = ZERO

      kx = nv*m_orig+n_columns
      call CSTAB_adjust_base(ux_vx(1:), a(1:,1),offset1)
      call CSTAB_adjust_base(ux_vx(1+kx:), a(1:,1),offset3)
      call CSTAB_adjust_base(uy_vy(1:), a(1:,1),offset2)
      call CSTAB_adjust_base(uy_vy(1+kx:), a(1:,1),offset4)
      call CSTAB_adjust_base(u_t(1:),a(1:,1),offset5)
      call CSTAB_adjust_base(v_t(1:),a(1:,1),offset6)
      call CSTAB_adjust_base(w(1:),a(1:,1),offset7)
      kx = (L1_LSIZE/8) + (L2_LSIZE/8)
      offset1 = offset1+kx*1
      offset2 = offset2+kx*2
      offset3 = offset3+kx*3
      offset4 = offset4+kx*4
      offset5 = offset5+kx*5
      offset6 = offset6+kx*6
      offset7 = offset7+kx*7
      call CSTAB_round_offset(offset1)
      call CSTAB_round_offset(offset2)
      call CSTAB_round_offset(offset3)
      call CSTAB_round_offset(offset4)
      call CSTAB_round_offset(offset5)
      call CSTAB_round_offset(offset6)
      call CSTAB_round_offset(offset7)

      kx = nv*m_orig+n_columns

!$OMP PARALLEL
      call eigen_prd_body(
     &     a, nm, d_out, e_out, ne,
     &     n, m_orig,
     &     w(1+offset7),
     &     ux_vx(1   +offset1), ! u_x(1:nv,m)
     &     uy_vy(1   +offset2), ! u_y(1:nv,m)
     &     ux_vx(1+kx+offset3), ! v_x(1:nv,m)
     &     uy_vy(1+kx+offset4), ! v_y(1:nv,m)
     &     u_t(1+offset5),
     &     v_t(1+offset6),
     &     d_t(1),
     &     nv)
!$OMP END PARALLEL

      deallocate(w)
      deallocate(v_t)
      deallocate(u_t)
      deallocate(d_t)
      deallocate(ux_vx)
      deallocate(uy_vy)

      return

      end subroutine eigen_prd_stub

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_body
!
! Purpose
! =======
!
! eigen_prd_body is the main body of the penta-diagonalization to
! transform a dense matrix to a penta-diagonal form.
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(nm,*)
!        On entry, a contains the input real symmetric matrix to be
!        transformed.
!        On exit, a contains the data of the hoseholder reflectors.
!
! nm     (input) integer
!        The leading dimension of the input array a.
!        nm >= ceil(N/x_nnod)
!
! d_out  (output) real(8) array, dimension(n)
!        d_out contains diagonal elements of the output penta-diagonal
!        matrix.
!
! e_out  (output) real(8) array, dimension(ne,*)
!        e_out contains off-diagonal elements of the output
!        penta-diagonal matrix.
!
! ne     (input) integer
!        The leading dimension of the output array e_out. ne >= N.
!
! n      (input) integer
!        The dimension of the input matrix. N >= 0.
!
! m_orig (input) integer
!        Block factor.
!
! w       (input/output) real(8) array, dimension(nm,*)
!         the panel buffer.
!
! u_x     (input/output) real(8) array, dimension(nv,*)
!         the working buffer.
!
! u_y     (input/output) real(8) array, dimension(nv,*)
!         the working buffer.
!
! v_x     (input/output) real(8) array, dimension(nv,*)
!         the working buffer.
!
! v_y     (input/output) real(8) array, dimension(nv,*)
!         the working buffer.
!
! u_t     (input/output) real(8) array, dimension(2*nv)
!         the working buffer.
!
! v_t     (input/output) real(8) array, dimension(2*nv)
!         the working buffer.
!
! d_t     (input/output) real(8) array, dimension(nv)
!         the working buffer.
!
! nv      (input) integer
!         The leading dimension of the working buffers.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_prd_body(a, nm, d_out, e_out, ne,
     &     n, m_orig,
     &     w, u_x, u_y, v_x, v_y, u_t, v_t, d_t, nv)

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(out)   :: d_out(1:n)
      real(8), intent(out)   :: e_out(1:ne,*)
      integer, intent(in)    :: ne
      integer, intent(in)    :: n
      integer, intent(in)    :: m_orig
      real(8)                :: w(1:nm,*)
      real(8)                :: u_x(1:nv,*)
      real(8)                :: u_y(1:nv,*)
      real(8)                :: v_x(1:nv,*)
      real(8)                :: v_y(1:nv,*)
      real(8)                :: u_t(1:2*nv)
      real(8)                :: v_t(1:2*nv)
      real(8)                :: d_t(1:nv)
      integer, intent(in)    :: nv


      real(8), save          :: c(MBAND,MBAND) ! c is used in t4x and t6
!     _3
      integer                :: i, m, m0, mm
      integer                :: k_1, k_2, k_3, k_4
      integer                :: i_block, i_base
      real(8)                :: d1, d2, dd(100)

!     
!     initialization
!     
#if TIMER_PRINT>1
      dd = ZERO
#endif
      m  = min(m_orig, n)

!$OMP MASTER

#if TIMER_PRINT
      if (TRD_inod == 1) then
        print*,"NUM.OF.PROCESS =",TRD_nnod,"(",x_nnod,y_nnod,")"
!$      print*,"NUM.OF.THREADS =",omp_get_num_threads()
      end if
      call flush(6)
#endif

      call eigen_prd_init(a(1,1), nm, n,
     &     d_out(1), e_out(1,1), ne,
     &     u_t(1), v_t(1), nv)

!$OMP END MASTER

!     
!     main iteration to reduce the input matrix into a penta-diagonal
!     
      mm = ((n-(MBAND+mod(n,MBAND)))-1)/m+1+1

!$OMP BARRIER

      do i_block=mm,2,-1


!$OMP BARRIER

        i_base = (i_block-2)*m+(MBAND+MOD(n,MBAND))
        m0     = min(m,n-i_base)
        if (m0 < 1) cycle

!$OMP MASTER

        call eigen_prd_panel_load(a, w, nm,
     &       d_t,
     &       u_x, u_y, v_x, v_y, nv,
     &       m0, i_base, i_block)

!$OMP END MASTER

        k_2 = m0
        k_3 = max(1,3*(2-i_block))

!$OMP BARRIER

        do k_1=k_2,k_3,-MBAND; k_4=k_1-MBAND+1

!$OMP BARRIER

          i = i_base+k_1
!     
!     u=...
!     
#if TIMER_PRINT>1
!$OMP MASTER
          d1=eigen_get_wtime()
!$OMP END MASTER
#endif

          call eigen_prd_compute_u(
     &         w(1:,k_4), nm,
     &         u_x(1:,k_4), u_y(1:,k_4), nv,
     &         u_t, v_t, i,
     &         c, e_out, ne)

#if TIMER_PRINT>1
!$OMP MASTER
          d2=eigen_get_wtime()
          dd(1)=dd(1)+(d2-d1)
!$OMP END MASTER
#endif
!$OMP BARRIER
#if TIMER_PRINT>1
!$OMP MASTER
          d1=eigen_get_wtime()
!$OMP END MASTER
#endif
!     
!     v:=Au
!     
          call eigen_prd_au(
     &         a, nm,
     &         u_x(1:,k_4), u_y(1:,k_4), v_x(1:,k_4), nv,
     &         u_t, v_t, d_t,
     &         i, i_base, m0)

#if TIMER_PRINT>1
!$OMP MASTER
          d2=eigen_get_wtime()
          dd(2)=dd(2)+(d2-d1)
!$OMP END MASTER
#endif
!$OMP BARRIER
#if TIMER_PRINT>1
!$OMP MASTER
          d1=eigen_get_wtime()
!$OMP END MASTER
#endif
!     
!     v=v-(UV+VU)u
!     v':= v-((u,v)/2|u|^2)u
!     
          call eigen_prd_compute_v(
     &         u_x, v_x, v_y, nv,
     &         u_t, v_t, c,
     &         i, i_base, m0)

#if TIMER_PRINT>1
!$OMP MASTER
          d2=eigen_get_wtime()
          dd(5)=dd(5)+(d2-d1)
!$OMP END MASTER
#endif
!$OMP BARRIER
#if TIMER_PRINT>1
!$OMP MASTER
          d1=eigen_get_wtime()
!$OMP END MASTER
#endif
!     
!     w':= w-uv-vu
!     
!$OMP MASTER
          call eigen_prd_local_2update(
     &         w, nm,
     &         u_x, u_y, v_x, v_y, nv,
     &         i_base, i)
!$OMP END MASTER

#if TIMER_PRINT>1
!$OMP MASTER
          d2=eigen_get_wtime()
          dd(6)=dd(6)+(d2-d1)
!$OMP END MASTER
#endif

        end do                  ! k_1

!$OMP BARRIER

!$OMP MASTER
        call eigen_prd_panel_store(a, w, nm,
     &       d_t,
     &       m0, i_base)
!$OMP END MASTER

!$OMP BARRIER

#if TIMER_PRINT>1
!$OMP MASTER
        d1=eigen_get_wtime()
!$OMP END MASTER
#endif
!     
!     A:=A-v^Tu-uT^v
!     
!     if (i_block >= 1) then
        call eigen_common_2update(
     &       a, nm,
     &       u_x, u_y, v_x, v_y, nv,
     &       m0, i_base)
!     end if

!$OMP BARRIER

#if TIMER_PRINT>1
!$OMP MASTER
        d2=eigen_get_wtime()
        dd(3)=dd(3)+(d2-d1)
!$OMP END MASTER
#endif

!$OMP BARRIER

      end do                    ! i_block

!$OMP BARRIER

!$OMP MASTER

      call eigen_prd_final(a, nm, n,
     &     d_out, e_out, ne, u_t)

#if TIMER_PRINT>1
      if (TRD_inod==1) then
        print*,"calc (u,beta)    ",dd(1)
        print*,"mat-vec (Au)     ",dd(2),dble(n)**3*2d-9/3/dd(2)
        print*,"2update (A-uv-vu)",dd(3),dble(n)**3*2d-9/3/dd(3)
!       print*,"calc v           ",dd(4)
        print*,"v=v-(UV+VU)u     ",dd(5)
        print*,"UV post reduction",dd(6)
      end if
      call flush(6)
#endif

!$OMP END MASTER

      return

      end subroutine            ! trd_body

      end module eigen_prd_mod
