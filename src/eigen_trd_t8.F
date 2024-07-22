!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t8_mod
!
! Purpose
! =======
!
! eigen_trd_t8_mod manages the modules for initializer and finalizer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_trd_t8_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use eigen_house_mod
      use CSTAB_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_init
      public :: eigen_trd_final

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_init
!
! Purpose
! =======
!
! initializer
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(nm,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! nm     (input) integer
!        The leading dimension of the array a.
!
! n      (input) integer
!        The dimension of matrix A
!
! u_t    (input) real(8) array, dimension(*)
!        working buffer
!
! v_t    (input) real(8) array, dimension(*)
!        working buffer
!
! nv     (input) integer
!        The leading dimension of the working arrays.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_init(a, nm, n,
     &     u_t, v_t, nv)

      real(8), intent(inout) ::  a(1:nm, *)
      integer, intent(in)    ::  nm
      integer, intent(in)    ::  n
      real(8), intent(in)    ::  u_t(*)
      real(8), intent(in)    ::  v_t(*)
      integer, intent(in)    ::  nv

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j

      integer                ::  local_size, local_rank, ierr

      include 'CSTAB.h'


      i_2 = eigen_loop_start(1, 'Y')
      i_3 = eigen_loop_end  (n, 'Y')
      do  i_1 = i_2, i_3
        j   = eigen_translate_l2g(i_1, 'Y')
        j_2 = eigen_loop_start   (j+1, 'X')
        if (j <= n) then
          a(j_2:nm, i_1) = ZERO
        else
          a(1:nm, i_1) = ZERO
        end if
      end do                    ! i_1

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if (local_size > 1) then

#if BOOST_BY_CACHE_ALIGNMENT
        allocate(u0_z(nv*local_size+n_columns),
     &           v0_z(nv*local_size+n_columns), stat=ierr)
#else
        allocate(u0_z(nv*local_size),
     &           v0_z(nv*local_size), stat=ierr)
#endif
        if (ierr /= 0) then
          call eigen_abort( "Memory allocation error [TRD/t8].", 1 )
          return
        end if
#if BOOST_BY_CACHE_ALIGNMENT
        call CSTAB_adjust_base(u0_z, u_t, offset1)
        call CSTAB_adjust_base(v0_z, v_t, offset2)
#else
        offset1 = 0
        offset2 = 0
#endif

      end if

      flag_overlap = .true.
      timer_t1 = ZERO
      timer_t2 = ONE

      return

      end subroutine eigen_trd_init

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_final
!
! Purpose
! =======
!
! finalizer
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(nm,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! nm     (input) integer
!        The leading dimension of the array a.
!
! n      (input) integer
!        The dimension of matrix A
!
! d_out  (output) real(8) array, dimension(*)
!        d_out contains the diagonal elements of the tridiagonal matrix.
!
! e_out  (output) real(8) array, dimension(*)
!        e_out contains the off-diagonal elements of the tridiagonal
!        matrix.
!
! u_t    (output) real(8) array, dimension(*)
!        working buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_final(a, nm, n, d_out, e_out, u_t)

      real(8), intent(inout) ::  a(1:nm, *)
      integer, intent(in)    ::  nm
      integer, intent(in)    ::  n
      real(8), intent(out)   ::  d_out(1:2)
      real(8), intent(out)   ::  e_out(1:2)
      real(8), intent(inout) ::  u_t(1:6)
      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j, k, L
      real(8)                ::  t

      integer                ::  local_size, local_rank

      include 'CSTAB.h'


      if (n >= 2) then

        i = 2; L = i-1
        i_1 = eigen_owner_index(i, 'Y')
        j_1 = eigen_owner_index(L, 'X')
        if (i_1 > 0 .and. j_1 > 0) then
          t           =   a(j_1, i_1)
          u_t(1)      = - t
          a(j_1, i_1) = 2*t
        else
          u_t(1)      = ZERO
        end if

        do i=1,2
          i_1 = eigen_owner_index(i, 'Y')
          j_1 = eigen_owner_index(i, 'X')
          if (i_1 > 0 .and. j_1 > 0) then
            u_t(1+i) = a(j_1, i_1)
          else
            u_t(1+i) = ZERO
          end if
        end do

        call reduce_dbl(u_t(1), u_t(4), 3, 2, x_COMM_WORLD)
        call reduce_dbl(u_t(1), u_t(4), 3, 2, y_COMM_WORLD)

        e_out(1) = ZERO
        e_out(2) = u_t(1)
        d_out(1) = u_t(2)
        d_out(2) = u_t(3)

      end if
      
      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if (local_size > 1) then

        deallocate(u0_z)
        deallocate(v0_z)

      end if

      return

      end subroutine eigen_trd_final

      end module eigen_trd_t8_mod
