!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t2_mod
!
! Purpose
! =======
!
! eigen_trd_t2_mod manages the modules for the PDSYMV routine called
! from eigen_trd.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_trd_t2_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use eigen_house_mod
      use mpi
!$    use omp_lib

      implicit none
      private

      public  :: eigen_trd_au
      private :: eigen_trd_au_sub_UV_ux
      private :: eigen_trd_au_sub_UV_uy
      private :: eigen_trd_au_body0
      private :: eigen_trd_au_body1
      private :: eigen_trd_au_body2
      private :: eigen_trd_au_body4
      private :: eigen_trd_au_body4_subXX
      private :: eigen_trd_au_body4_subX1
      private :: eigen_trd_au_body4_sub11

#if __FUJITSU
!     for FX10 and K
#   define	DO_UNROLL	5
#   define	K_LOOP		0
#   define	V_LOOP		0
#   define	DO_ITR		1
#   define	VLEN		(16*11)
#endif

#if __IBM_REGISTER_VARS
!     for BlueGene/Q
#   define	DO_UNROLL	4
#   define	K_LOOP		1
#   define	V_LOOP		1
#   define	DO_ITR		1
#   define	VLEN		(16*15)
#endif

#if SX
!     for SX-ACE
#   define	DO_UNROLL	7
#   define	K_LOOP		0
#   define	V_LOOP		0
#   define	DO_ITR		1
#   define	VLEN		(16*16*16)
#endif

#if __INTEL_COMPILER
!     for x86_64
#   define	DO_UNROLL	4
#   define	K_LOOP		0
#   define	V_LOOP		1
#   define	DO_ITR		4
#   define	VLEN		(16*10)
#endif

#if !defined(DO_UNROLL) || DO_UNROLL <= 1
!     for other platform
#   define	DO_UNROLL	6
#   define	K_LOOP		1
#   define	V_LOOP		1
#   define	DO_ITR		1
#   define	VLEN		(16*15)
#endif

#define	SPLIT_LOOP	0

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au
!
! Purpose
! =======
!
! paralleli SYMV, v := A u
!
!
! Arguments
! =========
!
! a       (input/output) real(8) array, dimension(nm,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! nm      (input) integer 
!         The leading dimension of the array a.
!
! u_x     (input/output) real(8) array, dimension(nv)
!         On entry, u_x contains the specfic row-vector from a.
!         On exit, u_x returns the reflector vector.
!
! u_y     (input/output) real(8) array, dimension(nv)
!         On entry, u_y contains the specfic row-vector from a.
!         On exit, u_y returns the reflector vector.
!
! v_x     (input/output) real(8) array, dimension(nv)
!         On entry, v_x is a zero vector.
!         On exit, v_x returns A * u_x.
!
! ux      (input) real(8) array, dimension(nv,*)
!         ux contains U, in which previously calculated u's are stored.
!
! vx      (input) real(8) array, dimension(nv,*)
!         vx contains V, in which previously calculated v's are stored.
!
! nv      (input) integer 
!         The leading dimension of the working arrays, u_x, u_y, ...
!
! u_t     (input/output) real(8) array, dimension(*)
!         working buffer, the size must be larger than nv+4*m+8
!
! v_t     (input/output) real(8) array, dimension(*)
!         working buffer, the size must be larger than nv+4*m+8
!
! d_t     (input) real(8) array, dimension(*)
!         d_t contains the diagonal elements of A.
!
! loop_info (input) integer array, dimension(*)
!         Information regarding loop and so on.
!
! d_out   (input/output) real(8)
!         diagonal part of the tridiagonal output matrix.
!
! e_out   (input/output) real(8)
!         off-diagonal part of the tridiagonal output matrix.
!
! beta    (input/output) real(8)
!         beta contains the sum of square of the reflector vector u_x.
!
! u_n     (output) real(8)
!         n-th element of u_x
!
! v_n     (output) real(8)
!         n-th element of v_x
!
! u_nn    (output) real(8) array, dimension(*)
!         u_nn contains a part of u_x(:) which must be commonly shared.
!
! v_nn    (output) real(8) array, dimension(*)
!         v_nn contains a part of v_x(:) which must be commonly shared.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_trd_au(
     &     a, nm,
     &     u_x, u_y, v_x,
     &     ux, vx, nv,
     &     u_t, v_t, d_t,
     &     loop_info,
     &     d_out, e_out, beta, u_n, v_n, u_nn, v_nn)

      real(8), intent(inout) ::  a(1:nm,*)
      integer, intent(in)    ::  nm
      real(8), intent(inout) ::  u_x(1:nv)
      real(8), intent(inout) ::  u_y(1:nv)
      real(8), intent(inout) ::  v_x(1:nv)
      real(8), intent(in)    ::  ux(1:nv,*)
      real(8), intent(in)    ::  vx(1:nv,*)
      integer, intent(in)    ::  nv
      real(8), intent(inout) ::  u_t(*)
      real(8), intent(inout) ::  v_t(*)
      real(8), intent(in)    ::  d_t(*)
      integer, intent(in)    ::  loop_info(*)
      real(8), intent(inout) ::  d_out
      real(8), intent(inout) ::  e_out
      real(8), intent(inout) ::  beta
      real(8), intent(out)   ::  u_n
      real(8), intent(out)   ::  v_n
      real(8), intent(out)   ::  u_nn(*)
      real(8), intent(out)   ::  v_nn(*)

      integer                ::  i, i_base, m
      integer                ::  x_pos, y_pos, x_root, y_root
      integer                ::  num_v_done, num_v_ptr
      integer                ::  n1, n2
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  L, j
      integer                ::  ierr

      real(8)                :: anorm2, a_n, g_n
      real(8)                :: prod_uv, prod_ua, d_n
      real(8)                :: w0, u0, v0, ux0, vx0

      real(8)                :: d1, d2

      integer                ::  local_size, local_rank, local_root

      real(8)                :: dcx(1:200)
      real(8)                :: dcy(1:200)


!---------------------------------------------------------------
!---------------------------------------------------------------

      local_rank = loop_info(LOOP_INFO_L_RANK)
      local_size = loop_info(LOOP_INFO_L_SIZE)

      num_v_ptr  = loop_info(LOOP_INFO_V_CURNT)
      i_base     = loop_info(LOOP_INFO_I_BASE)
      m          = loop_info(LOOP_INFO_V_WIDTH)
      num_v_done = m - num_v_ptr

      i = i_base + num_v_ptr
      L = i - 1

      i_2 = loop_info(LOOP_INFO_Y_START)
      i_3 = loop_info(LOOP_INFO_Y_END)
      i_4 = loop_info(LOOP_INFO_Y_OINDX)

      j_2 = loop_info(LOOP_INFO_X_START)
      j_3 = loop_info(LOOP_INFO_X_END)
      j_4 = loop_info(LOOP_INFO_X_OINDX)

      y_root = loop_info(LOOP_INFO_Y_ROOT)
      y_pos  = loop_info(LOOP_INFO_Y_POS )

      x_root = loop_info(LOOP_INFO_X_ROOT)
      x_pos  = loop_info(LOOP_INFO_X_POS )


      if (num_v_done == 0) then
        call eigen_trd_au_body0(
     &       a, nm,
     &       i_2, i_3,
     &       local_rank, local_size)
!$OMP BARRIER
      end if

#if DEBUG
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-body0 passed"
#endif

!---------------------------------------------------------------
!---------------------------------------------------------------

      if (local_size > 1) then

        if (flag_overlap) then
          local_root = 1
        else
          local_root = 0
        end if

        if (local_rank >= local_root) then

          if (.not. flag_oncache) then
          if (local_rank == local_root) then
          if (num_v_ptr == 1) then
            timer_t2 = eigen_get_wtime()
          end if
          end if
          end if

          n1 = offset1+nv*(local_rank-local_root)
          n2 = offset2+nv*(local_rank-local_root)

          u0_z(1+n1:x_pos+n1) = ZERO
          v0_z(1+n2:y_pos+n2) = ZERO

          call eigen_trd_au_body1(
     &         a, nm,
     &         u_x, u_y, u0_z(1+n1:), v0_z(1+n2:), nv,
     &         i_2, i_3, x_pos, y_pos
     &         ,local_rank-local_root, local_size-local_root )

          if (.not. flag_oncache) then
          if (local_rank == local_root) then
          if (num_v_ptr == 1) then
            timer_t2 = eigen_get_wtime() - timer_t2
          end if
          end if
          end if

        end if

!$OMP BARRIER

        call eigen_trd_au_body2(
     &       v_x, u0_z(1+offset1:),
     &       v_t, v0_z(1+offset2:),
     &       nv, x_pos, y_pos
     &       ,local_size-local_root )

      else

        v_x(1:x_pos) = ZERO
        v_t(1:y_pos) = ZERO

        call eigen_trd_au_body1(
     &       a, nm,
     &       u_x, u_y, v_x, v_t, nv,
     &       i_2, i_3, x_pos, y_pos
     &       ,local_rank, local_size )

#if DEGBU
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-body1 passed"
#endif

      end if

!---------------------------------------------------------------
!---------------------------------------------------------------

!$OMP BARRIER

      if (local_size == 1 .or. local_rank >= 1) then

        if (m > num_v_ptr) then
        if (num_v_done > 0) then
          if ( y_inod == y_root ) then
            call eigen_trd_au_sub_UV_ux(
     &           loop_info, num_v_done,
     &           u_x,
     &           ux(1:, num_v_ptr+1), vx(1:, num_v_ptr+1), nv,
     &           u_t(1), local_size, local_rank )
          else
            u_t(1:4*num_v_done) = ZERO
          end if
        end if
        end if

#if DEBUG
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-body UV passed"
#endif

      end if

      if (local_rank == 0) then

        anorm2  = ZERO          ! norm(u0, 2)^2
        a_n     = ZERO          ! u0_{L} = a_{L,i}
        prod_uv = ZERO          ! u^tAu
        d_n     = ZERO          ! a_{L,L}
        prod_ua = ZERO          ! u0^t A_{:,L}


                                ! u = u0 - s * e_L, s = - sign(u0_{L})*
                                !                         norm(u0,2)
                                ! -> 
                                ! u^t A u = (u0-se)^tA(u0-se)
                                !         = u0^t A u0 - 2sa_{:,L}^tu0 
                                !                     + s^2a_{L,L}
                                !         = u0^t A u0 - 2*g_n*prod_ua
                                !                     + g_n^2*d_n
                                ! prod_uv = prod_uv + g_n * (g_n * d_n -
                                !                       2 * prod_ua)


        prod_uv = ddot(y_pos, v_t, 1, u_y, 1)

        if (i_4 > 0) then
!DIR$ IVDEP
!OCL NOEVAL
          do  j_1=1,x_pos
            w0  = u_x(j_1)
            anorm2  = anorm2  + w0**2
            prod_uv = prod_uv + w0 * v_x(j_1)
            prod_ua = prod_ua + w0 * a(j_1, i_4)
          end do
          if (j_4 > 0) then
            d_n = d_t(i_4)
            prod_ua = prod_ua + d_n * u_x(j_4)
          end if
        else

          prod_uv = prod_uv + ddot(x_pos, v_x, 1, u_x, 1)

        end if

        if (j_4 > 0) then
          a_n = u_x(j_4)
        end if


                                !  :: v_x += diag(A) * u_y   if diag(<L)
        call  eigen_trd_au_body4(
     &       u_y, v_x, d_t,
     &       L-1, nv, prod_uv)

#if DEBUG
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-body4 passed"
#endif

        u_nn(1:m) = ZERO
        if (num_v_ptr-1 > 0) then
          k_2 = eigen_loop_start (i_base+1, 'X')
          k_3 = eigen_loop_end   (i-1,      'X')
          do k_1=k_2,k_3
            j = eigen_translate_l2g(k_1, 'X')
            u_nn(j-i_base) = u_x(k_1)
          end do
        end if
      end if

!---------------------------------------------------------------
!---------------------------------------------------------------

!$OMP BARRIER

      if (local_rank == 0) then
        if (x_nnod > 1) then

          j_1 = y_pos
          call pack1_dbl(anorm2,  v_t, j_1)
          call pack1_dbl(a_n,     v_t, j_1)
          call pack1_dbl(prod_uv, v_t, j_1)
          call pack1_dbl(d_n,     v_t, j_1)
          call pack1_dbl(prod_ua, v_t, j_1)
          call pack1_dbl(d_out,   v_t, j_1)
          if (num_v_ptr-1 > 0) then
            call pack_dbl(u_nn, num_v_ptr-1, v_t, j_1)
          end if
          if (num_v_done > 0 .and. num_v_ptr > 0) then
            call pack_dbl(v_nn, num_v_ptr, v_t, j_1)
          end if
          if (.not. flag_oncache) then
            if (num_v_ptr == 1) then
              call pack1_dbl(timer_t1, v_t, j_1)
              call pack1_dbl(timer_t2, v_t, j_1)
            end if
          end if
          if (num_v_done > 0) then
            call pack_dbl(u_t, 4*num_v_done, v_t, j_1)
          end if
 
! nv + 6 + (m1-1) + m1 + 2+ 4*(m-m1) = nv + 4*m - 2*m1 + 7 <= nv+4*m+7
          call reduce_dbl(v_t, u_t, j_1, 1, x_COMM_WORLD)

#if DEBUG
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-reduce(X) passed"
#endif

          j_1 = y_pos
          call unpack1_dbl(anorm2,  v_t, j_1)
          call unpack1_dbl(a_n,     v_t, j_1)
          call unpack1_dbl(prod_uv, v_t, j_1)
          call unpack1_dbl(d_n,     v_t, j_1)
          call unpack1_dbl(prod_ua, v_t, j_1)
          call unpack1_dbl(d_out,   v_t, j_1)
          if ( num_v_ptr-1 > 0 ) then
            call unpack_dbl(u_nn, num_v_ptr-1, v_t, j_1)
          end if
          if (num_v_done > 0 .and. num_v_ptr > 0) then
            call unpack_dbl(v_nn, num_v_ptr, v_t, j_1)
          end if
          if (.not. flag_oncache) then
            if (num_v_ptr == 1) then
              call unpack1_dbl(timer_t1, v_t, j_1)
              call unpack1_dbl(timer_t2, v_t, j_1)
            end if
          end if
          if (num_v_done > 0) then
            call unpack_dbl(u_t, 4*num_v_done, v_t, j_1)
          end if

        end if

!---------------------------------------------------------------
!---------------------------------------------------------------

        if (i_4 > 0) then
          if (anorm2 /= ZERO) then
            anorm2  =  sqrt(anorm2)
            g_n     = -sign(anorm2, a_n) ! s
            u_n     =  a_n - g_n
          else
            g_n     =  ZERO
            u_n     =  ZERO
          end if

          prod_uv = prod_uv + d_n * u_y(i_4)**2
          u_y(i_4) =  u_n

          ! v = Au = A(u0-se) = v0 - g_n * a_{:,L}
          call daxpy(x_pos, -g_n, a(1, i_4), 1, v_x(1), 1)

          if (j_4 > 0) then
            v_x(j_4) = v_x(j_4) + d_t(i_4) * u_y(i_4)
          end if

          v_n  = v_t(i_4)
        else
          v_n     = ZERO
          anorm2  = ZERO
        end if

                                !  :: v_x += v_t   if diag
        call  eigen_trd_au_body3(
     &       v_x, v_t,
     &       L-0, nv
     &       )

#if DEBUG
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-body3 passed"
#endif

      end if

!---------------------------------------------------------------
!---------------------------------------------------------------

!$OMP BARRIER

      if (local_rank == 0) then

        if (y_nnod > 1) then

          j_1 = 4*num_v_done
          call pack1_dbl(anorm2,         u_t, j_1)
          call pack1_dbl(prod_uv,        u_t, j_1)
          call pack1_dbl(prod_ua,        u_t, j_1)
          call pack1_dbl(d_n,            u_t, j_1)
          call pack1_dbl(v_n,            u_t, j_1)
          call pack_dbl (v_x,     x_pos, u_t, j_1)
          if (.not. flag_oncache) then
            if (num_v_ptr == 1) then
              call pack1_dbl(timer_t1,   u_t, j_1)
              call pack1_dbl(timer_t2,   u_t, j_1)
            end if
          end if

! 1 + 4*m + 5 + nv + 2 = nv + 4*m + 8
          call reduce_dbl(u_t, v_t, j_1, 2, y_COMM_WORLD)

#if DEBUG
      if ( local_rank ==0 .and. TRD_inod == 1 )
     &  print*,"TRD t2-reduce(Y) passed"
#endif

          j_1 = 4*num_v_done
          call unpack1_dbl(anorm2,         u_t, j_1)
          call unpack1_dbl(prod_uv,        u_t, j_1)
          call unpack1_dbl(prod_ua,        u_t, j_1)
          call unpack1_dbl(d_n,            u_t, j_1)
          call unpack1_dbl(v_n,            u_t, j_1)
          call unpack_dbl (v_x,     x_pos, u_t, j_1)
          if (.not. flag_oncache) then
            if (num_v_ptr == 1) then
              call unpack1_dbl(timer_t1,   u_t, j_1)
              call unpack1_dbl(timer_t2,   u_t, j_1)
            end if
          end if

        end if

!---------------------------------------------------------------
!---------------------------------------------------------------

        if (anorm2 /= ZERO) then
          g_n      = -sign(anorm2, a_n) ! s
          u_n      =  a_n - g_n
!          beta     =  anorm2 * abs(u_n) ! norm(u)^2 / 2
          beta     = -u_n * g_n
        else
          g_n      =  ZERO
          u_n      =  ZERO
          beta     =  ONE
        end if
        e_out      =  g_n

        if (num_v_ptr-1 > 0) then
          u_nn(num_v_ptr-1)  =  u_n
        end if

        if (j_4 > 0) then
          u_x(j_4) =  u_n
        end if

        v_n        = v_n + d_n * u_n

        if (num_v_done > 0) then
          !
          ! here u_t is updated as u_t(:,1) += u_n * u(:,4)
          !                     as u_t(:,3) += u_n * u(:,2)
          !
          call daxpy(num_v_done, u_n,
     &         u_t(1+3*num_v_done), 1, u_t(1+0*num_v_done), 1)
          call daxpy(num_v_done, u_n,
     &         u_t(1+2*num_v_done), 1, u_t(1+1*num_v_done), 1)
          !
          v_n      = v_n
     &       - ddot( 2*num_v_done,
     &         u_t(1+0*num_v_done), 1, u_t(1+2*num_v_done), 1)
        end if

        if ( j_4 > 0 ) then
          prod_uv  = prod_uv + g_n * (g_n * d_n - 2 * prod_ua)
          u_t(1+2*num_v_done)  = prod_uv
        endif

!---------------------------------------------------------------
!---------------------------------------------------------------

        if (.not. flag_oncache) then
        if (num_v_ptr == 1) then
          w0 = ONE
          v0 = ONE
          if (x_pos * y_pos * 8 / 2 < 5 * 1024 * 1024) then
            ! if on-cache
            if ( local_size > 1 ) then
              ! if multi-threaded
              if ( flag_overlap ) then
                ! if currently overlapped
                v0 = DBLE(local_size-1) / local_size
              else
                ! if currently non-overlapped
                w0 = DBLE(local_size) / (local_size-1)
              end if
            end if
          end if
          ! if overlapped
          d1 = max(timer_t2*w0, timer_t1)
          ! if not-overlapped
          d2 = timer_t2*v0 + timer_t1 - Bcast_cont_Overhead_x
          flag_overlap = (d1 < d2)
        end if
        end if

      end if

!---------------------------------------------------------------
!---------------------------------------------------------------

      return

      end subroutine  eigen_trd_au

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_sub_UV_ux
!
! Purpose
! =======
!
! compute  (U^t,V^t) u_x
!
!
! Arguments
! =========
!
! loop_info (input) integer array, dimension(*)
!         Information regarding loop and so on.
!
! m       (input) integer 
!         The number of vectors stored in U and V.
!
! u_x     (input) real(8) array, dimension(*)
!         contains the reflector vector
!
! ux      (input) real(8) array, dimension(nv,*)
!         U(:,1:m)
!
! vx      (input) real(8) array, dimension(nv,*)
!         V(:,1:m)
!
! nv      (input) integer
!         The leading dimension of the buffers.
!
! u_t     (output) real(8) array, dimension(m,*)
!         u_t contains the (U^t,V^t) u_x
!
! local_rank(input) integer 
!         Thread ID
!
! local_size(input) integer 
!         The number of threads in the main team
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_sub_UV_ux(
     &     loop_info, m,
     &     u_x, ux, vx, nv, u_t,
     &     local_size, local_rank)

      integer, intent(in)   :: loop_info(*)
      integer, intent(in)   :: m
      real(8), intent(in)   :: u_x(*)
      real(8), intent(in)   :: ux(1:nv, *)
      real(8), intent(in)   :: vx(1:nv, *)
      integer, intent(in)   :: nv
      real(8), intent(out)  :: u_t(1:m, *)
      integer, intent(in)   :: local_size
      integer, intent(in)   :: local_rank
      
      integer               :: i_1, i_2, i_3, i_4
      integer               :: j_1, j_2, j_3, j_4
      real(8)                :: u0, v0, w0, ux0, vx0


      if (m < 1) return

      j_2 = loop_info(LOOP_INFO_X_START)
      j_3 = loop_info(LOOP_INFO_X_END)
      j_4 = loop_info(LOOP_INFO_X_OINDX)
      if (j_4 > 0) j_3 = j_3-1

      i_2 = 1+max(local_rank-1,0)
      i_3 = m
      i_4 = max(local_size-1,1)

      do i_1=i_2,i_3,i_4
        u0 = ZERO
        v0 = ZERO
        do j_1=j_2,j_3
          w0 = u_x(j_1)
          u0 = u0 + ux(j_1, i_1) * w0
          v0 = v0 + vx(j_1, i_1) * w0
        end do
        ux0 = ZERO
        vx0 = ZERO
        if (j_4 > 0) then
          ux0 = ux(j_4, i_1)
          vx0 = vx(j_4, i_1)
        end if
        u_t(i_1, 1) = v0
        u_t(i_1, 2) = u0
        u_t(i_1, 3) = ux0
        u_t(i_1, 4) = vx0
      end do

      return

      end subroutine eigen_trd_au_sub_UV_ux

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_sub_UV_uy
!
! Purpose
! =======
!
! compute  (U^t,V^t) u_y
!
!
! Arguments
! =========
!
! loop_info (input) integer array, dimension(*)
!         Information regarding loop and so on.
!
! m       (input) integer 
!         The number of vectors stored in U and V.
!
! u_y     (input) real(8) array, dimension(*)
!         contains the reflector vector
!
! uy      (input) real(8) array, dimension(nv,*)
!         U(:,1:m)
!
! vy      (input) real(8) array, dimension(nv,*)
!         V(:,1:m)
!
! nv      (input) integer
!         The leading dimension of the buffers.
!
! u_t     (output) real(8) array, dimension(m,*)
!         u_t contains the (U^t,V^t) u_y
!
! local_rank(input) integer 
!         Thread ID
!
! local_size(input) integer 
!         The number of threads in the main team
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_sub_UV_uy(
     &     loop_info, m,
     &     u_y, uy, vy, nv, u_t,
     &     local_size, local_rank)

      integer, intent(in)   :: loop_info(*)
      integer, intent(in)   :: m
      real(8), intent(in)   :: u_y(*)
      real(8), intent(in)   :: uy(1:nv, *)
      real(8), intent(in)   :: vy(1:nv, *)
      integer, intent(in)   :: nv
      real(8), intent(out)  :: u_t(1:m, *)
      integer, intent(in)   :: local_size
      integer, intent(in)   :: local_rank
      
      integer               :: i_1, i_2, i_3, i_4
      integer               :: j_1, j_2, j_3, j_4
      real(8)               :: u0, v0, w0, ux0, vx0


      if (m < 1) return

      i_2 = loop_info(LOOP_INFO_Y_START)
      i_3 = loop_info(LOOP_INFO_Y_END)
      i_4 = loop_info(LOOP_INFO_Y_OINDX)

      j_2 = 1+max(local_rank-1,0)
      j_3 = m
      j_4 = max(local_size-1,1)

      do j_1=j_2,j_3,j_4
        u0 = ZERO
        v0 = ZERO
        do i_1=i_2,i_3
          w0 = u_y(i_1)
          u0 = u0 + uy(i_1, j_1) * w0
          v0 = v0 + vy(i_1, j_1) * w0
        end do
        ux0 = ZERO
        vx0 = ZERO
        if (i_4 > 0) then
          ux0 = uy(i_4, j_1)
          vx0 = vy(i_4, j_1)
        end if
        u_t(j_1, 1) = v0
        u_t(j_1, 2) = u0
        u_t(j_1, 3) = ux0
        u_t(j_1, 4) = vx0
      end do

      return

      end subroutine eigen_trd_au_sub_UV_uy

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body0
!
! Purpose
! =======
!
! Fill zeros for lower triangular part.
!
!
! Arguments
! =========
!
! a       (input/output) real(8) array, dimension(nm,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! nm      (input) integer 
!         The leading dimension of the array a.
!
! n1      (input) integer 
!         initial of loop interval
!
! n2      (input) integer 
!         termination of loop interval
!
! local_rank(input) integer 
!         Thread ID
!
! local_size(input) integer 
!         The number of threads in the main team
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body0(
     &     a, nm,
     &     n1, n2
     &     ,local_rank, local_size)

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: n1
      integer, intent(in)    :: n2
      integer                :: local_size
      integer                :: local_rank

      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4, j


      i_2 = n1
      i_3 = n2

      do i_1=i_2+local_rank,i_3,local_size
        j   = eigen_translate_l2g(i_1, 'Y')
        j_3 = eigen_loop_end  (j, 'X')
        j   = j+y_nnod*(DO_UNROLL*DO_ITR-1)
        j_4 = eigen_loop_end  (j, 'X')
        do j_1=j_3+1,min(j_4,nm)
          a(j_1,i_1) = ZERO
        end do                  ! j_1
      end do                    ! i_1

      return

      end subroutine eigen_trd_au_body0

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body1
!
! Purpose
! =======
!
! < purpose of this subroutine ... >
!
!
! Arguments
! =========
!
! a       (input/output) real(8) array, dimension(nm,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! nm      (input) integer 
!         The leading dimension of the array a.
!
! u_x     (input/output) real(8) array, dimension(nv)
!         On entry, u_x contains the row-vector from a.
!         On exit, u_x returns the reflector vector.
!
! u_y     (input/output) real(8) array, dimension(nv)
!         On entry, u_y contains the row-vector from a.
!         On exit, u_y returns the reflector vector.
!
! u_t     (input/output) real(8) array, dimension(nv)
!         returns upper(A) * u
!
! v_t     (input/output) real(8) array, dimension(nv)
!         returns lower(A) * u
!
! nv      (input) integer 
!         leading dimension of u_[xyt] and v_t
!
! n1      (input) integer 
!         initial of loop interval
!
! n2      (input) integer 
!         termination of loop interval
!
! x_pos   (input) integer 
!         local index corresponding to the source of the reflector.
!
! y_pos   (input) integer 
!         local index corresponding to the source of the reflector.
!
! local_rank(input) integer 
!         Thread ID
!
! local_size(input) integer 
!         The number of threads in the main team
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body1(
     &     a, nm,
     &     u_x, u_y, u_t, v_t, nv, 
     &     n1, n2, x_pos, y_pos,
     &     local_rank, local_size)

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(inout) :: u_x(1:nv)
      real(8), intent(inout) :: u_y(1:nv)
      real(8), intent(inout) :: u_t(1:nv)
      real(8), intent(inout) :: v_t(1:nv)
      integer, intent(in)    :: nv
      integer, intent(in)    :: n1
      integer, intent(in)    :: n2
      integer, intent(in)    :: x_pos
      integer, intent(in)    :: y_pos
      integer, intent(in)    :: local_size
      integer, intent(in)    :: local_rank

      integer                :: i_0
      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: k_1, k_2, k_3, k_4
      integer                :: l_1, l_2, l_3, l_4
      integer                :: i, j, k

      real(8)                :: v0, u0
      real(8)                :: a0_0
      real(8)                :: a0_1
      real(8)                :: a0_2
      real(8)                :: a0_3
      real(8)                :: a0_4
      real(8)                :: a0_5
      real(8)                :: a0_6
      real(8)                :: a0_7
      real(8)                :: w0_(0:7)
#define	w0_(x)		w_0_ x
      real(8)                :: w_0_0
      real(8)                :: w_0_1
      real(8)                :: w_0_2
      real(8)                :: w_0_3
      real(8)                :: w_0_4
      real(8)                :: w_0_5
      real(8)                :: w_0_6
      real(8)                :: w_0_7
      real(8)                :: v_t0, u_y0
      real(8)                :: v_t1, u_y1
      real(8)                :: v_t2, u_y2
      real(8)                :: v_t3, u_y3
      real(8)                :: v_t4, u_y4
      real(8)                :: v_t5, u_y5
      real(8)                :: v_t6, u_y6
      real(8)                :: v_t7, u_y7

      integer                :: LX, LY
      integer                :: ii_1, ii_2, ii_3, ii_4, ii_5
      integer                :: jj_1, jj_2, jj_3, jj_4, jj_5
      integer                :: kk_1, kk_2, kk_3, kk_4, kk_5


      i_2 = n1
      i_3 = n2
!     
!     v:= Au
!     
      LX = max(VLEN*4,96)

      l_2 = i_2; l_3 = i_3
      l_1 = l_2; l_4 = l_3

      k_2 = 1
      k   = eigen_translate_l2g(l_4, 'Y')
      k_3 = eigen_loop_end(k-1, 'X')

#if K_LOOP
      do k_1=k_2,k_3,LX; k_4 = min(k_3,k_1+LX-1)
#else
        k_1=k_2; k_4 = k_3
#endif

        j    = eigen_translate_l2g(k_1, 'X')
        ii_2 = eigen_loop_start(j, 'Y')
        ii_2 = max(l_1,ii_2)
        ii_3 = l_4
#if K_LOOP
        if (ii_2 > ii_3) cycle
#else
        if (ii_2 > ii_3) goto 99999
#endif
        ii_4 = mod(ii_3-ii_2+1,DO_UNROLL*DO_ITR)+ii_2


        if (ii_2 < ii_4) then
!-------------------------------------------------------------
          do i_1=ii_2+local_rank,ii_4-1,local_size
            j    = eigen_translate_l2g(i_1, 'Y')
            j    = j+(1-1)*y_nnod
            jj_2 = k_1
            jj_3 = eigen_loop_end  (j-1, 'X')
            jj_3 = min(k_4,jj_3)
            if (jj_2 > jj_3) cycle
!=============================================================
#if V_LOOP
            do kk_1=jj_2,jj_3,VLEN
              kk_4=min(kk_1+VLEN-1,jj_3)
#else
              kk_1=jj_2; kk_4=jj_3
#endif

              w0_(0) = v_t(i_1+0)

#if SPLIT_LOOP
!GCC$ UNROLL 8
#if __INTEL_COMPILER
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
#endif
#if __IBM_REGISTER_VARS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
#endif
#if __FUJITSU
#if SPLIT_LOOP
!OCL LOOP_NOFUSION
#endif
!OCL UNROLL(8)
!OCL SIMD
#endif
              do j_1=kk_1,kk_4
                u0 = u_x(j_1+0)
                a0_0 = a(j_1+0,i_1+0)
                w0_(0) = w0_(0)
     &               + (a0_0*u0)
              end do            ! j_1
              v_t(i_1+0) = w0_(0)
#endif

              u_y0 = u_y(i_1+0)
!GCC$ UNROLL 8
#if __INTEL_COMPILER
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
#endif
#if __IBM_REGISTER_VARS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
#endif
#if __FUJITSU
#if SPLIT_LOOP
!OCL LOOP_NOFUSION
#endif
!OCL UNROLL(8)
!OCL SIMD
#endif
              do j_1=kk_1,kk_4

                a0_0 = a(j_1+0,i_1+0)

#if !SPLIT_LOOP
                u0 = u_x(j_1+0)
                w0_(0) = w0_(0)
     &               + (a0_0*u0)
#endif

                u_t(j_1+0) = u_t(j_1+0)
     &               + (a0_0*u_y0)

              end do            ! j_1
#if !SPLIT_LOOP
              v_t(i_1+0) = w0_(0)
#endif

#if V_LOOP
            end do              ! kk_1
#endif
!=============================================================
          end do                ! i_1
!-------------------------------------------------------------
        end if


        if (ii_4 <= ii_3) then
!-------------------------------------------------------------
          do i_0 = ii_4+local_rank*DO_UNROLL*DO_ITR,
     &         ii_3,
     &         local_size*DO_UNROLL*DO_ITR

            j    = eigen_translate_l2g(i_0, 'Y')
            j    = j+(DO_UNROLL*DO_ITR-1)*y_nnod
            jj_2 = k_1
            jj_3 = eigen_loop_end(j-1, 'X')
            jj_3 = min(k_4,jj_3)
            if (jj_2 > jj_3) cycle
!=============================================================
#if V_LOOP
            do kk_1=jj_2,jj_3,VLEN
              kk_4=min(kk_1+VLEN-1,jj_3)
#else
              kk_1=jj_2; kk_4=jj_3
#endif

#if DO_ITR>1
              do i_1 = i_0,
     &             i_0+DO_UNROLL*DO_ITR-1,
     &             DO_UNROLL
#else
                i_1 = i_0
#endif

                w0_(0) = v_t(i_1+0)
                w0_(1) = v_t(i_1+1)

#if SPLIT_LOOP
!GCC$ UNROLL 8
#if __INTEL_COMPILER
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
!DIR$ UNROLL(8)
#endif
#if __IBM_REGISTER_VARS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
#endif
#if __FUJITSU
#if SPLIT_LOOP
!OCL LOOP_NOFUSION
#endif
!OCL UNROLL(8)
!OCL SIMD
#endif
                do j_1=kk_1,kk_4
                  u0 = u_x(j_1+0)
                  a0_0 = a(j_1+0,i_1+0)
                  a0_1 = a(j_1+0,i_1+1)
                  w0_(0) = w0_(0)
     &                 + (a0_0*u0)
                  w0_(1) = w0_(1)
     &                 + (a0_1*u0)
                end do
                v_t(i_1+0) = w0_(0)
                v_t(i_1+1) = w0_(1)
#endif

#               if DO_UNROLL>=3
                w0_(2) = v_t(i_1+2)
#               endif
#               if DO_UNROLL>=4
                w0_(3) = v_t(i_1+3)
#               endif
#               if DO_UNROLL>=5
                w0_(4) = v_t(i_1+4)
#               endif
#               if DO_UNROLL>=6
                w0_(5) = v_t(i_1+5)
#               endif
#               if DO_UNROLL>=7
                w0_(6) = v_t(i_1+6)
#               endif
#               if DO_UNROLL>=8
                w0_(7) = v_t(i_1+7)
#               endif

                u_y0 = u_y(i_1+0)
                u_y1 = u_y(i_1+1)
#               if DO_UNROLL>=3
                u_y2 = u_y(i_1+2)
#               endif
#               if DO_UNROLL>=4
                u_y3 = u_y(i_1+3)
#               endif
#               if DO_UNROLL>=5
                u_y4 = u_y(i_1+4)
#               endif
#               if DO_UNROLL>=6
                u_y5 = u_y(i_1+5)
#               endif
#               if DO_UNROLL>=7
                u_y6 = u_y(i_1+6)
#               endif
#               if DO_UNROLL>=8
                u_y7 = u_y(i_1+7)
#               endif

!GCC$ UNROLL 4
#if __INTEL_COMPILER
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
!DIR$ UNROLL(4)
#endif
#if __IBM_REGISTER_VARS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
#endif
#if __FUJITSU
#if SPLIT_LOOP
!OCL LOOP_NOFUSION
#endif
!OCL UNROLL(4)
!OCL SIMD
#endif
                do j_1=kk_1,kk_4

                  a0_0 = a(j_1+0,i_1+0)
                  a0_1 = a(j_1+0,i_1+1)

                  u0 = u_x(j_1+0)
                  v0 = ZERO

#if !SPLIT_LOOP
                  w0_(0) = w0_(0)
     &                 + (a0_0*u0)
                  w0_(1) = w0_(1)
     &                 + (a0_1*u0)
#endif

                  v0 = v0
     &                 + (a0_0*u_y0)
     &                 + (a0_1*u_y1)

#                 if DO_UNROLL>=3
                  a0_2 = a(j_1+0,i_1+2)
#                 endif
#                 if DO_UNROLL>=4
                  a0_3 = a(j_1+0,i_1+3)
#                 endif
#                 if DO_UNROLL>=3
                  w0_(2) = w0_(2)
     &                 + (a0_2*u0)
#                 endif
#                 if DO_UNROLL>=4
                  w0_(3) = w0_(3)
     &                 + (a0_3*u0)
#                 endif
#                 if DO_UNROLL>=3
                  v0 = v0
     &                 + (a0_2*u_y2)
#                 endif
#                 if DO_UNROLL>=4
     &                 + (a0_3*u_y3)
#                 endif

#                 if DO_UNROLL>=5
                  a0_4 = a(j_1+0,i_1+4)
#                 endif
#                 if DO_UNROLL>=6
                  a0_5 = a(j_1+0,i_1+5)
#                 endif
#                 if DO_UNROLL>=5
                  w0_(4) = w0_(4)
     &                 + (a0_4*u0)
#                 endif
#                 if DO_UNROLL>=6
                  w0_(5) = w0_(5)
     &                 + (a0_5*u0)
#                 endif
#                 if DO_UNROLL>=5
                  v0 = v0
     &                 + (a0_4*u_y4)
#                 endif
#                 if DO_UNROLL>=6
     &                 + (a0_5*u_y5)
#                 endif

#                 if DO_UNROLL>=7
                  a0_6 = a(j_1+0,i_1+6)
#                 endif
#                 if DO_UNROLL>=8
                  a0_7 = a(j_1+0,i_1+7)
#                 endif
#                 if DO_UNROLL>=7
                  w0_(6) = w0_(6)
     &                 + (a0_6*u0)
#                 endif
#                 if DO_UNROLL>=8
                  w0_(7) = w0_(7)
     &                 + (a0_7*u0)
#                 endif
#                 if DO_UNROLL>=7
                  v0 = v0
     &                 + (a0_6*u_y6)
#                 endif
#                 if DO_UNROLL>=8
     &                 + (a0_7*u_y7)
#                 endif

                  u_t(j_1+0) = v0 + u_t(j_1+0)

                end do          ! j_1
#if !SPLIT_LOOP
                v_t(i_1+0) = w0_(0)
                v_t(i_1+1) = w0_(1)
#endif
#               if DO_UNROLL>=3
                v_t(i_1+2) = w0_(2)
#               endif
#               if DO_UNROLL>=4
                v_t(i_1+3) = w0_(3)
#               endif
#               if DO_UNROLL>=5
                v_t(i_1+4) = w0_(4)
#               endif
#               if DO_UNROLL>=6
                v_t(i_1+5) = w0_(5)
#               endif
#               if DO_UNROLL>=7
                v_t(i_1+6) = w0_(6)
#               endif
#               if DO_UNROLL>=8
                v_t(i_1+7) = w0_(7)
#               endif

#if DO_ITR>1
              end do            ! i_1
#endif

#if V_LOOP
            end do              ! kk_1
#endif
!=============================================================
          end do                ! i_0
!-------------------------------------------------------------
        end if
99999   continue
#if K_LOOP
      end do                    ! k_1
#endif

      return

      end subroutine eigen_trd_au_body1

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body2
!
! Purpose
! =======
!
! Sum up the results by each thread.
! v_x <= (u_t in this routine) upper(A)*u
! v_t <= (v_t in this routine) lower(A)*u
!
!
! Arguments
! =========
!
! u_t     (output) real(8) array, dimension(nv)
!         returns upper(A) * u
!
! u_z     (input) real(8) array, dimension(nv,*)
!         contains results of upper(A) * u done by each thread
!
! v_t     (output) real(8) array, dimension(nv)
!         returns lower(A) * u
!
! v_z     (input) real(8) array, dimension(nv,*)
!         contains results of lower(A) * u done by each thread
!
! nv      (input) integer 
!         The leading dimension of the buffers.
!
! x_pos   (input) integer 
!         local index corresponding to the source of the reflector.
!
! y_pos   (input) integer 
!         local index corresponding to the source of the reflector.
!
! local_rank(input) integer 
!         Thread ID
!
! local_size(input) integer 
!         The number of threads in the main team
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body2(
     &     u_t, u_z, v_t, v_z, nv,
     &     x_pos, y_pos,
     &     local_size)

      real(8), intent(out)   :: u_t(1:nv)
      real(8), intent(in)    :: u_z(1:nv,*)
      real(8), intent(out)   :: v_t(1:nv)
      real(8), intent(in)    :: v_z(1:nv,*)
      integer, intent(in)    :: nv
      integer, intent(in)    :: x_pos
      integer, intent(in)    :: y_pos
      integer, intent(in)    :: local_size

      integer, parameter     :: LX = 32

      integer                :: j_1, j_3, j_4, j
      integer                :: jj_1, jj_2, jj_3, jj_4

      integer                :: ll_size, ll_rank


!$OMP DO
      do jj_1=1,x_pos,LX
        j_3=jj_1; j_4=min(jj_1+LX-1,x_pos)

        if (mod(local_size,2) == 1) then
          do j_1=j_3,j_4
            u_t(j_1) = u_z(j_1,1)
          end do
        else
          do j_1=j_3,j_4
            u_t(j_1) = u_z(j_1,1)+u_z(j_1,2)
          end do
        end if
        if (local_size > 2) then
          do j=3-mod(local_size,2),local_size,2
            do j_1=j_3,j_4
              u_t(j_1) = u_t(j_1)+u_z(j_1,j+0)+u_z(j_1,j+1)
            end do
          end do
        end if

      end do
!$OMP END DO
!$OMP DO
      do jj_1=1,y_pos,LX
        j_3=jj_1; j_4=min(jj_1+LX-1,y_pos)

        if (mod(local_size,2) == 1) then
          do j_1=j_3,j_4
            v_t(j_1) = v_z(j_1,1)
          end do
        else
          do j_1=j_3,j_4
            v_t(j_1) = v_z(j_1,1)+v_z(j_1,2)
          end do
        end if
        if (local_size > 2) then
          do j=3-mod(local_size,2),local_size,2
            do j_1=j_3,j_4
              v_t(j_1) = v_t(j_1)+v_z(j_1,j+0)+v_z(j_1,j+1)
            end do
          end do
        end if

      end do
!$OMP END DO

      return

      end subroutine eigen_trd_au_body2

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body3
!
! Purpose
! =======
!
! Sumup the rest part of SYMV,
!   v_x:=Upper(A)*u+diag(A)*u
!       =v_x(:) + diag(A)*u(:)
!
!
! Arguments
! =========
!
! v_x     (input/output) real(8) array, dimension(nv)
!         v_x contains Upper(A)*u(:)
!
! v_t     (input) real(8) array, dimension(*)
!         v_t contains Lower(A)*u(:)
!
! n       (input) integer
!         loop length
!
! nv      (input) integer 
!         The leading dimension of the buffers.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body3(
     &     v_x, v_t,
     &     n, nv)

      real(8), intent(inout) :: v_x(1:nv)
      real(8), intent(in)    :: v_t(*)
      integer, intent(in)    :: n
      integer, intent(in)    :: nv

      integer                :: i0, j0, k0
      integer                :: nx, nm1, nm2


      if (x_nnod /= y_nnod) then
        call eigen_diag_loop_info( 1, n,
     &             k0, nx, i0, nm2, j0, nm1 )
        nx = nx + 1
        call daxpy(nx, ONE, v_t(i0), nm2, v_x(j0), nm1)
      else
        if (x_inod==y_inod)then
          nx = eigen_loop_end  (n, 'X')
          call daxpy(nx, ONE, v_t, 1, v_x, 1)
        end if
      end if

      return

      end subroutine  eigen_trd_au_body3

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body4
!
! Purpose
! =======
!
! Sumup the rest part of SYMV,
!   v:=Upper(A)*u+Lower(A)*u+diag(A)*u
!     =v_x(:,1:2) + v_t(:,1:2)
!   prod_uv:=(u, v)
!
!
! Arguments
! =========
!
! u_y     (input) real(8) array, dimension(nv)
!         u_y contains u(:)
!
! v_x     (input/output) real(8) array, dimension(nv)
!         v_x contains Upper(A)*u(:)
!
! d_t     (input) real(8) array, dimension(nv)
!         d_t contains the diagonal elements of A.
!
! n       (input) integer
!         loop length
!
! nv      (input) integer 
!         The leading dimension of the buffers.
!
! prod_uv (input/output) real(8)
!         prod_uv returns (u, v).
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body4(
     &     u_y, v_x, d_t,
     &     n, nv, prod_uv)

      real(8), intent(in)    :: u_y(1:nv)
      real(8), intent(inout) :: v_x(1:nv)
      real(8), intent(in)    :: d_t(1:nv)
      integer, intent(in)    :: n
      integer, intent(in)    :: nv
      real(8), intent(inout) :: prod_uv

      integer                :: i0, j0, k0
      integer                :: nx, nm1, nm2


      if (x_nnod /= y_nnod) then
        call eigen_diag_loop_info( 1, n,
     &             k0, nx, i0, nm2, j0, nm1 )
        nx = nx + 1
        if (nx > 0) then
          if (nm2 == 1) then
          if (nm1 == 1) then
            call eigen_trd_au_body4_sub11(
     &               v_x(j0), d_t(i0), u_y(i0),
     &               nx, nm1, nm2, prod_uv)
          else
            call eigen_trd_au_body4_subX1(
     &               v_x(j0), d_t(i0), u_y(i0),
     &               nx, nm1, nm2, prod_uv)
          end if
          else
            call eigen_trd_au_body4_subXX(
     &             v_x(j0), d_t(i0), u_y(i0),
     &             nx, nm1, nm2, prod_uv)
          end if
        end if
      else
        if (x_inod == y_inod) then
          nx = eigen_loop_end  (n, 'X')
          call eigen_trd_au_body4_sub11(
     &               v_x, d_t, u_y, nx, 1, 1, prod_uv)
        end if
      end if

      return

      end subroutine eigen_trd_au_body4

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body4_subXX
!
! Purpose
! =======
!
! subroutine for body4
!
!
! Arguments
! =========
!
! v_x     (input/output) real(8) array, dimension(nm1,*)
!         v_x contains Upper(A)*u(:)
!
! d_t     (input) real(8) array, dimension(nm2,*)
!         d_t contains the diagonal elements of A.
!
! u_y     (input) real(8) array, dimension(nm2,*)
!         u_y contains u(:)
!
! n       (input) integer
!         loop length
!
! nm1     (input) integer
! nm2     (input) integer
!
! prod_uv (input/output) real(8)
!         prod_uv returns (u, v).
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body4_subXX(v_x,d_t,u_y,n,nm1,nm2,
     &     prod_uv)

      real(8), intent(inout) :: v_x(nm1,*)
      real(8), intent(in)    :: d_t(nm2,*)
      real(8), intent(in)    :: u_y(nm2,*)
      integer, intent(in)    :: n
      integer, intent(in)    :: nm1
      integer, intent(in)    :: nm2
      real(8), intent(inout) :: prod_uv

      real(8)                :: s, t
      integer                :: i


!DIR$ VECTOR ALWAYS
      do i=1,n
        t = u_y(1,i)
        s = d_t(1,i) * t
        v_x(1,i) = v_x(1,i) + s
        prod_uv  = prod_uv  + s * t
      end do                    ! i


      return

      end subroutine eigen_trd_au_body4_subXX

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body4_subX1
!
! Purpose
! =======
!
! subroutine for body4
!
!
! Arguments
! =========
!
! v_x     (input/output) real(8) array, dimension(nm1,*)
!         v_x contains Upper(A)*u(:)
!
! d_t     (input) real(8) array, dimension(*)
!         d_t contains the diagonal elements of A.
!
! u_y     (input) real(8) array, dimension(*)
!         u_y contains u(:)
!
! n       (input) integer
!         loop length
!
! nm1     (input) integer
! nm2     (input) integer
!
! prod_uv (input/output) real(8)
!         prod_uv returns (u, v).
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body4_subX1(v_x,d_t,u_y,n,nm1,nm2,
     &     prod_uv)

      real(8), intent(inout) :: v_x(nm1,*)
      real(8), intent(in)    :: d_t(*)
      real(8), intent(in)    :: u_y(*)
      integer, intent(in)    :: n
      integer, intent(in)    :: nm1
      integer, intent(in)    :: nm2
      real(8), intent(inout) :: prod_uv

      real(8)                :: s, t
      integer                :: i


!DIR$ VECTOR ALWAYS
      do i=1,n
        t = u_y(i)
        s = d_t(i) * t
        v_x(1,i) = v_x(1,i) + s
        prod_uv  = prod_uv  + s * t
      end do                    ! i

      return

      end subroutine eigen_trd_au_body4_subX1

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_au_body4_sub11
!
! Purpose
! =======
!
! subroutine for body4
!
!
! Arguments
! =========
!
! v_x     (input/output) real(8) array, dimension(n)
!         v_x contains Upper(A)*u(:)
!
! d_t     (input) real(8) array, dimension(n)
!         d_t contains the diagonal elements of A.
!
! u_y     (input) real(8) array, dimension(n)
!         u_y contains u(:)
!
! n       (input) integer
!         loop length
!
! nm1     (input) integer
! nm2     (input) integer
!
! prod_uv (input/output) real(8)
!         prod_uv returns (u, v).
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_au_body4_sub11(v_x, d_t,u_y,n,nm1,nm2,
     &     prod_uv)

      real(8), intent(inout) :: v_x(1:n)
      real(8), intent(in)    :: d_t(1:n)
      real(8), intent(in)    :: u_y(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: nm1
      integer, intent(in)    :: nm2
      real(8), intent(inout) :: prod_uv

      real(8)                :: s, t
      integer                :: i


!DIR$ VECTOR ALWAYS
      do i=1,n
        t = u_y(i)
        s = d_t(i) * t
        v_x(i)  = v_x(i) + s
        prod_uv = prod_uv + s * t
      end do                    ! i

      return

      end subroutine  eigen_trd_au_body4_sub11


      subroutine xx_check( x, y, n, comm, name )
      real(8), intent(in)  :: x
      real(8), intent(out) :: y(n)
      integer, intent(in)  :: n, comm
      character(*), intent(in) :: name
      integer :: i, nz
      real(8) :: dummy(1)

      dummy(1) = x
      call allgather_dbl( dummy, y, 1, 0, comm )

      nz = 0
      do i=1,n
        if ( y(i) /= 0 ) nz = nz + 1
      enddo
      if ( nz<= 1 ) then
         print*,"Bcast-replacement by allreduce: ",name,"@",
     &              x_inod, y_inod, "::", nz
      endif

      end subroutine

      subroutine yy_check( x, m, y, n, comm, name )
      real(8), intent(in)  :: x(m)
      real(8), intent(out) :: y(n)
      integer, intent(in)  :: m, n, comm
      character(*), intent(in) :: name
      integer :: ii, i, nz, nzz

      nzz = 0
      do ii=1,m
        call allgather_dbl( x(ii:), y, 1, 0, comm )
        nz = 0
        do i=1,n
          if ( y(i) /= 0 ) nz = nz + 1
        enddo
        if ( nz > 1 ) nzz = nzz + 1
      enddo

      if ( nzz == 0 ) then
         print*,"Bcast-replacement by allreduce: ",name,"@",
     &              x_inod, y_inod
      endif

      end subroutine

      subroutine co_check( x, y, n, comm, name )
      real(8), intent(in)  :: x(n)
      real(8), intent(out) :: y(n)
      integer, intent(in)  :: n, comm
      character(*), intent(in) :: name
      integer :: i

      if ( n >= 1 ) then
        y(1:n) = x(1:n)
        call bcast_dbl( y, n, 1, 0, comm )
        do i=1,n
          if ( y(i) /= x(i) ) then
             print*,"Not identical: ",name,"(",i,")@",
     &              x_inod,y_inod,y(i),x(i)
             exit
          endif
        enddo
      endif

      end subroutine

      end module eigen_trd_t2_mod
