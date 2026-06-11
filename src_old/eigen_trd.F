!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_mod
!
! Purpose
! =======
!
! eigen_trd_mod manages the modules of the tri-diagonalization of
! eigen_s.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_trd_mod

      use eigen_t1_mod
      use eigen_trd_t2_mod
      use eigen_trd_t4_mod
      use eigen_trd_t5_mod
      use eigen_trd_t5x_mod
      use eigen_trd_t6_3_mod
      use eigen_trd_t7_mod
      use eigen_trd_t8_mod

      use eigen_house_mod

      use comm_mod
      use eigen_devel_mod
      use eigen_libs_mod

      use CSTAB_mod
      use mpi
!$    use omp_lib

      implicit none
      private

      public  :: eigen_trd
      private :: eigen_trd_stub
      private :: eigen_trd_body

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd
!
! Purpose
! =======
!
! eigen_trd invokes the stub function of the tri-diagonalization,
! eigen_trd_stub, to transform a dense matrix to a tri-diagonal form.
!
!
! Arguments
! =========
!
! n       (input) integer
!         The dimension of the input matrix. N >= 0.
!
! a       (input/output) real(8) array, dimension(*)
!         On entry, a contains the input real symmetric matrix to be
!         transformed.
!         On exit, a contains the data of the hoseholder reflectors.
!
! nma0    (input) integer
!         The leading dimension of the input array a.
!         nma0 >= ceil(N/x_nnod)
!
! d_out   (output) real(8) array, dimension(n)
!         d_out contains diagonal elements of the output tri-diagonal
!         matrix.
!
! e_out   (output) real(8) array, dimension(n)
!         e_out contains off-diagonal elements of the output
!         tri-diagonal matrix.
!
! m0      (input) integer
!         Block factor
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd(n, a, nma0, d_out, e_out, m0)

      integer, intent(in)    :: n
      real(8), intent(inout) :: a(*)
      integer, intent(in)    :: nma0
      real(8), intent(inout) :: d_out(1:n)
      real(8), intent(inout) :: e_out(1:n)
      integer, intent(in)    :: m0

      include 'CSTAB.h'

      integer                :: nm, m, ierr, lda


      nm = nma0
      m  = m0

      call MPI_Barrier(TRD_COMM_WORLD, ierr)
      call eigen_timer_reset(BCAST=2, REDUCE=3, REDIST=2, GATHER=2)

      lda = nma0

      call eigen_trd_stub(a, d_out, e_out, n, lda, m)

      call MPI_Barrier(TRD_COMM_WORLD,ierr)

      comm_time_reduction =
     &     eigen_timer_print('EigenExa (Tri-diagonalization)')

      return

      end subroutine eigen_trd

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_stub
!
! Purpose
! =======
!
! eigen_trd_stub invokes the main body of the tri-diagonalization,
! eigen_trd_body, to transform a dense matrix to a tri-diagonal form.
!
!
! Arguments
! =========
!
! a       (input/output) real(8) array, dimension(*)
!         On entry, a contains the input real symmetric matrix to be
!         transformed.
!         On exit, a contains the data of the hoseholder reflectors.
!
! d_out   (output) real(8) array, dimension(n)
!         d_out contains diagonal elements of the output tri-diagonal
!         matrix.
!
! e_out   (output) real(8) array, dimension(n)
!         e_out contains off-diagonal elements of the output
!         tri-diagonal matrix.
!
! n       (input) integer
!         The dimension of the input matrix. N >= 0.
!
! nm      (input) integer
!         The leading dimension of the input array a.
!         nm >= ceil(N/x_nnod)
!     
! m       (input) integer
!         Block factor
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_stub(a, d_out, e_out, n, nm, m)

      real(8), intent(inout) :: a(nm,*)
      real(8), intent(inout) :: d_out(1:n)
      real(8), intent(inout) :: e_out(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: nm
      integer, intent(in)    :: m

      include 'CSTAB.h'

      integer, parameter     :: nm_max_L1 = 16*4
      integer, parameter     :: nm_max_L2 = 16*4*2

      real(8), allocatable   :: u_t(:), v_t(:)
      real(8), allocatable   :: ux_vx(:)
      real(8), allocatable   :: uy_vy(:)
      real(8), allocatable   :: uvw_nn(:)

      integer                :: nv, nx, ny, nz
      integer                :: ierr, kx
      integer                :: offset1, offset2, offset3
      integer                :: offset4, offset5, offset6
      integer                :: offset7

      integer(8)             :: l_array_1, l_array_2


      nx = (n-1)/x_nnod+1 +2
      nv = nm

      call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv)

      if (x_nnod /= y_nnod) then
        nx = (n-1)/x_nnod+1
        kx = y_nnod / n_common
        nz = ((nx-1)/kx+1)
        ny = nz*MAX(x_nnod,y_nnod)
      else
        nx = (n-1)/x_nnod+1
        nz = nx
        ny = nz*x_nnod
      end if
!
      l_array_1 = MAX(MAX(m,2)*ny, nx+4*m+8)
      l_array_2 = MAX(MAX(m,2)*nz, nx+4*m+8)

#if BOOST_BY_CACHE_ALIGNMENT
      allocate(
     &     u_t(1:l_array_1+n_columns),
     &     v_t(1:l_array_2+n_columns),
     &     ux_vx(1:nv*2*m+2*n_columns),
     &     uy_vy(1:nv*2*m+2*n_columns),
     &     uvw_nn(m*3),
     &     stat=ierr)
#else
      allocate(
     &     u_t(1:l_array_1),
     &     v_t(1:l_array_2),
     &     ux_vx(1:nv*2*m),
     &     uy_vy(1:nv*2*m),
     &     uvw_nn(m*3),
     &     stat=ierr)
#endif
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error [TRD].", 1 )
      end if

      v_t = ZERO
      u_t = ZERO
      ux_vx = ZERO
      uy_vy = ZERO

#if BOOST_BY_CACHE_ALIGNMENT
      kx = nv*m+n_columns
      call CSTAB_adjust_base(ux_vx(1:), a(1:,1),offset1)
      call CSTAB_adjust_base(ux_vx(1+kx:), a(1:,1),offset3)
      call CSTAB_adjust_base(uy_vy(1:), a(1:,1),offset2)
      call CSTAB_adjust_base(uy_vy(1+kx:), a(1:,1),offset4)
      call CSTAB_adjust_base(u_t(1:),a(1:,1),offset5)
      call CSTAB_adjust_base(v_t(1:),a(1:,1),offset6)
      kx = (L1_LSIZE/8) + (L2_LSIZE/8)
      offset1 = offset1+kx*1
      offset2 = offset2+kx*3
      offset3 = offset3+kx*5
      offset4 = offset4+kx*2
      offset5 = offset5+kx*4
      offset6 = offset6+kx*6
      call CSTAB_round_offset(offset1)
      call CSTAB_round_offset(offset2)
      call CSTAB_round_offset(offset3)
      call CSTAB_round_offset(offset4)
      call CSTAB_round_offset(offset5)
      call CSTAB_round_offset(offset6)
      kx = nv*m+n_columns
#else
      kx = nv*m
      offset1=0
      offset2=0
      offset3=0
      offset4=0
      offset5=0
      offset6=0
#endif

      call MPI_Barrier(TRD_COMM_WORLD, ierr)

!$OMP PARALLEL
      call eigen_trd_body(a, nm, d_out, e_out, n, nv, m,
     &     ux_vx(1   +offset1), ! u_x(1:nv, 1:m)
     &     uy_vy(1   +offset2), ! u_y(1:nv, 1:m)
     &     ux_vx(1+kx+offset3), ! v_x(1:nv, 1:m)
     &     uy_vy(1+kx+offset4), ! v_y(1:nv, 1:m)
     &     u_t(1+offset5),      ! u_t(1:2*nv)
     &     v_t(1+offset6),      ! v_t(1:2*nv)
     &     uvw_nn(1:),          ! u_nn(1:m)
     &     uvw_nn(1+1*m:),      ! v_nn(1:m)
     &     uvw_nn(1+2*m:))      ! w_nn(1:m)
!$OMP END PARALLEL

      call MPI_Barrier(TRD_COMM_WORLD, ierr)

      deallocate(v_t)
      deallocate(u_t)
      deallocate(ux_vx)
      deallocate(uy_vy)
      deallocate(uvw_nn)

      return

      end subroutine eigen_trd_stub

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_body
!
! Purpose
! =======
!
! eigen_trd_body is the main body of the tri-diagonalization to
! transform a dense matrix to a tri-diagonal form.
!
!
! Arguments
! =========
!
! a       (input/output) real(8) array, dimension(*)
!         On entry, a contains the input real symmetric matrix to be
!         transformed.
!         On exit, a contains the data of the hoseholder reflectors.
!
! nm      (input) integer
!         The leading dimension of the input array a.
!         nm >= ceil(N/x_nnod)
!     
! d_out   (output) real(8) array, dimension(n)
!         d_out contains diagonal elements of the output tri-diagonal
!         matrix.
!
! e_out   (output) real(8) array, dimension(n)
!         e_out contains off-diagonal elements of the output
!         tri-diagonal matrix.
!
! n       (input) integer
!         The dimension of the input matrix. N >= 0.
!
! nv      (input) integer
!         The leading dimension of the working buffers.
!
! m_orig  (input) integer
!         Block factor.
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
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_body(a, nm, d_out, e_out, n, nv, m_orig,
     &     u_x, u_y, v_x, v_y, u_t, v_t, u_nn, v_nn, w_nn)
      implicit none
      include 'CSTAB.h'

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(inout) :: d_out(1:n)
      real(8), intent(inout) :: e_out(1:n)

      integer, intent(in)    :: n
      integer, intent(in)    :: nv
      integer, intent(in)    :: m_orig
      real(8)                :: u_x(1:nv,*)
      real(8)                :: u_y(1:nv,*)
      real(8)                :: v_x(1:nv,*)
      real(8)                :: v_y(1:nv,*)
      real(8)                :: u_t(1:2*nv)
      real(8)                :: v_t(1:2*nv)
      real(8)                :: u_nn(*)
      real(8)                :: v_nn(*)
      real(8)                :: w_nn(*)

      integer, parameter     :: MBAND_ = 1

      real(8), save          :: c(MBAND_,MBAND_)
      real(8), save          :: u_n, v_n
      integer                :: i, k_1, k_2, k_3, k_4, iref, iref_prev
      integer                :: m0, mm, m
      integer                :: i_block, i_base
      integer                :: local_rank, local_size
#if 1
      integer                :: loop_info(16)
#endif
      real(8)                :: d1, d2, dd(100)

#if TIMER_PRINT>1
      dd = ZERO
#endif
      m = min(m_orig, n)
     
!
!     initialization
!     

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

#if 1
      loop_info(LOOP_INFO_L_RANK) = local_rank
      loop_info(LOOP_INFO_L_SIZE) = local_size
#endif

!$OMP MASTER

#if TIMER_PRINT>1
      if (TRD_inod == 1) then
!$      if (local_rank == 0) then
          print*,"NUM.OF.PROCESS =",TRD_nnod,"(",x_nnod,y_nnod,")"
!$        print*,"NUM.OF.THREADS =",local_size
!$      end if
      end if
      call flush(6)
#endif

      call eigen_trd_init(
     &        a, nm, n,
     &        u_t, v_t, nv )

!$OMP END MASTER

      mm = ((n-mod(n,MBAND_))-1)/m+1

      do i_block = mm, max(1,3*(2-m)),-1

!$OMP BARRIER

        i_base = (i_block-1)*m
        m0     = min(m,n-i_base)

        call eigen_trd_panel_load(
     &           a, nm, u_x, nv,
     &           d_out, v_x,
     &           m0, i_base, i_block )

!$OMP MASTER
#if 0
         k_2 = eigen_translate_g2l (i_base+m0, 'X')
         k_3 = eigen_translate_g2l (i_base+m0, 'Y')
         flag_oncache = (k_2 * k_3 / 2 * 2 < 5 * (L2_SIZE/8))
#else
         flag_oncache = .TRUE.
#endif
!$OMP END MASTER

!$OMP BARRIER

        k_2 = m0
        k_3 = max(1,3*(2-i_block))

        do k_1=k_2,k_3,-1

          iref      = k_1-MBAND_+1
          iref_prev = iref+1

          i = i_base+k_1
#if DEBUG
          if ( TRD_inod == 1 ) print*,"TRD{",i,"}"
#endif


          loop_info(LOOP_INFO_X_START) =
     &         eigen_loop_start    (1,   'X')
          loop_info(LOOP_INFO_X_END  ) =
     &         eigen_loop_end      (i-1, 'X')
          loop_info(LOOP_INFO_X_OINDX) =
     &         eigen_owner_index   (i-1, 'X')
          loop_info(LOOP_INFO_X_ROOT ) =
     &         eigen_owner_node    (i-1, 'X')
          loop_info(LOOP_INFO_X_POS  ) =
     &         eigen_translate_g2l (i-1, 'X')

          loop_info(LOOP_INFO_Y_START) =
     &         eigen_loop_start    (1,   'Y')
          loop_info(LOOP_INFO_Y_END  ) =
     &         eigen_loop_end      (i-1, 'Y')
          loop_info(LOOP_INFO_Y_OINDX) =
     &         eigen_owner_index   (i-1, 'Y')
          loop_info(LOOP_INFO_Y_ROOT ) =
     &         eigen_owner_node    (i-1, 'Y')
          loop_info(LOOP_INFO_Y_POS  ) =
     &         eigen_translate_g2l (i-1, 'Y')

          loop_info(LOOP_INFO_V_WIDTH ) = m0
          loop_info(LOOP_INFO_V_CURNT ) = k_1
          loop_info(LOOP_INFO_I_BASE  ) = i_base

!$        if (local_rank == 0) then
          if (k_1 < k_2) then

#if TIMER_PRINT>1
            d1=eigen_get_wtime()
#endif
!     
!     w':= w-uv-vu
!     
            call eigen_trd_local_2update0(
     &           u_x(1:,iref),
     &           u_x(1:,iref_prev), v_x(1:,iref_prev),
     &           i-i_base, i+1, u_n, v_n)

#if TIMER_PRINT>1
            d2=eigen_get_wtime()
            dd(6)=dd(6)+(d2-d1)
#endif

            u_nn(1:k_2) = w_nn(1:k_2)

          end if

!     
!     u = Householder reflactor(...)
!
     
#if TIMER_PRINT>1
          d1=eigen_get_wtime()
#endif

          call eigen_trd_compute_u(
     &         d_out(i),
     &         u_x(1:,iref), u_y(1:,iref), nv,
     &         u_t, v_t, i,
     &         (k_2/=k_1), (k_1==k_3),
     &         v_x(1:,iref_prev), v_y(1:,iref_prev))

#if DEBUG
!$OMP MASTER
          if ( TRD_inod == 1 ) print*,"TRD t4 passed{",i,"}"
!$OMP END MASTER
#endif

#if TIMER_PRINT>1
          d2=eigen_get_wtime()
          dd(1)=dd(1)+(d2-d1)
#endif

!$          else
!$OMP BARRIER
!$          end if

#if TIMER_PRINT>1
!$OMP MASTER
          d1=eigen_get_wtime()
!$OMP END MASTER
#endif
!     
!     v:=Au
!     

          call eigen_trd_au(
     &         a, nm,
     &         u_x(1:,iref), u_y(1:,iref), v_x(1:,iref),
     &         u_x, v_x, nv,
     &         u_t, v_t, d_out,
     &         loop_info,
     &         d_out(i), e_out(i), c(1,1),
     &         u_n, v_n, w_nn, v_nn)


!$        if (local_size == 1 .or.
!$   &        local_rank > 0) then
            if (k_1 < k_2) then

              call eigen_trd_local_2update(
     &             u_x, nv,
     &             u_x(1:,iref_prev), v_x(1:,iref_prev),
     &             i-i_base-1, i+1, u_nn, v_nn,
     &             local_size, local_rank )

            end if
!$        end if

!$OMP BARRIER

#if TIMER_PRINT>1
!$OMP MASTER
          d2=eigen_get_wtime()
          dd(2)=dd(2)+(d2-d1)
!$OMP END MASTER
#endif

#if DEBUG
!$OMP MASTER
          if ( TRD_inod == 1 ) print*,"TRD au passed{",i,"}"
!$OMP END MASTER
#endif

!     
!     v=v-(UV+VU)u
!     v':= v-((u,v)/2|u|^2)u
!     

#if TIMER_PRINT>1
!$OMP MASTER
          d1=eigen_get_wtime()
!$OMP END MASTER
#endif

          call eigen_trd_compute_v(
     &         u_x(1:,iref), v_x(1:,iref), v_y(1:,iref),
     &         u_x, v_x, nv,
     &         u_t, v_t, c(1,1), u_n, v_n,
     &         loop_info, (k_3/=k_1), v_nn )

#if TIMER_PRINT>1
!$OMP MASTER
          d2=eigen_get_wtime()
          dd(4)=dd(4)+(d2-d1)
!$OMP END MASTER
#endif

        end do                  ! k_1

!$OMP BARRIER

        if (i_base == 0) then
          k_1=k_3; iref=k_1-MBAND_+1
          i = k_1

#if TIMER_PRINT>1
          d1=eigen_get_wtime()
#endif
!     
!     w':= w-uv-vu
!     
!$OMP MASTER
          call reduce_dbl(v_nn, u_t, i-1, 1, x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER

#if DEBUG
          if ( TRD_inod == 1 ) print*,"TRD 2update{",i,"}"
#endif

          call eigen_trd_local_2update(
     &         u_x, nv,
     &         u_x(1:,iref), v_x(1:,iref),
     &         i-i_base-1, i, w_nn, v_nn,
     &         local_size, local_rank)

#if TIMER_PRINT>1
          d2=eigen_get_wtime()
          dd(5)=dd(5)+(d2-d1)
#endif

        end if

!$OMP BARRIER

        call eigen_trd_panel_restore(
     &       a, nm, u_x, nv,
     &       d_out,
     &       m0, i_base)

!$OMP BARRIER

#if TIMER_PRINT>1
!$OMP MASTER
        d1=eigen_get_wtime()
!$OMP END MASTER
#endif
!     
!     A:=A-v^Tu-uT^v
!     

        if (i_block > 1) then
          call eigen_common_2update(
     &         a, nm,
     &         u_x,u_y,v_x,v_y, nv,
     &         m0, i_base)
        end if

!$OMP BARRIER

#if TIMER_PRINT>1
!$OMP MASTER
        d2=eigen_get_wtime()
        dd(3)=dd(3)+(d2-d1)
!$OMP END MASTER
!$OMP BARRIER
#endif

      end do                    ! i_block

!$OMP BARRIER

!$OMP MASTER
      call eigen_trd_final(
     &        a, nm, n,
     &        d_out, e_out, u_t)
!$OMP END MASTER


#if TIMER_PRINT>1
      if (TRD_inod==1) then

!$OMP BARRIER
!$      if (local_size == 1 .or.
!$   &      local_rank == 1) then
          u_t(1) = dd(5)
!$      end if
!$OMP BARRIER
!$OMP MASTER
        dd(5) = u_t(1)
!$OMP END MASTER
!$OMP BARRIER

!$OMP MASTER
        print*,"calc (u,beta)    ", dd(1)
        print*,"mat-vec (Au)     ", dd(2),dble(n)**3*2d-9/3/dd(2)
        print*,"2update (A-uv-vu)", dd(3),dble(n)**3*2d-9/3/dd(3)
        print*,"v=v-(UV+VU)u     ", dd(4)
        print*,"UV post reduction", dd(5)
!$OMP END MASTER

      end if
      call flush(6)
#endif

      return

      end subroutine eigen_trd_body

      end module eigen_trd_mod
