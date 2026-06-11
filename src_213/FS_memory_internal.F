!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function : FS_memory_internal
!
! Purpose
! =======
!
! Returns size of the required memory internally while proceeding the
! FS library.
! 
! Note
! =======
!
! This routine is modified from EigenExa-2.4b eigen_memory_internal
!
! =========

      integer(8) function FS_memory_internal(n, lda, ldz,
     &     m1_opt, m0_opt)
     &     result( byte )
      use FS_libs_mod
      use comm_mod
      use eigen_blacs_mod
      use eigen_devel_mod
      use CSTAB_mod
      use eigen_libs_mod
      use mpi
!$    use omp_lib
      implicit none

      integer, intent(in)           :: n, lda, ldz
      integer, intent(in), optional :: m1_opt, m0_opt

      integer, parameter :: nm_max_L1 = 16*4
      integer, parameter :: nm_max_L2 = 16*6
      integer, parameter :: REAL_SIZE = 8_8 ! must be in integer(8)
      integer, parameter :: INT_SIZE  = 4_8 ! must be in integer(8)

      integer    :: m1, m0
      integer    :: MYROW, MYCOL, NPROW, NPCOL, NP, NQ, NB
      integer    :: local_rank, local_size
      integer    :: m , nv, nm, ierr

      integer    :: lddz, lddw
      integer(8) :: LWORK, LIWORK
      integer    :: na, nx, ny, nz, kx
      integer    :: l_array_1, l_array_2

      integer(8) :: nax(2)
      integer(8) :: byte_tridi
      integer(8) :: byte_dcx
      integer(8) :: byte_trbak

      integer, external  :: NUMROC

      include 'CSTAB.h'

      nm = 0

      if (n <= 0) then
        byte = -1
        return
      end if

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if (present(m1_opt)) then
        m1 = m1_opt
      else
        m1 = 48
      end if

      if (present(m0_opt)) then
        m0 = m0_opt
      else
        m0 = 128
      end if

!     ===== for TRIDI ====

      nx = (n-1)/x_nnod+1 +2
      call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv)

      nx = (n-1)/x_nnod+1
      kx = y_nnod / n_common
      nz = ((nx-1)/kx+1)
      ny = nz*MAX(x_nnod,y_nnod)
!
      m = 0
      l_array_1 = MAX(MAX(m,2)*ny, nx+4*m+6)
      l_array_2 = MAX(MAX(m,2)*nz, nx+4*m+6)

      LWORK =
     &       (l_array_1+n_columns)          ! u_t
     &     + (l_array_2+n_columns)          ! v_t
     &     + MAX(3*m1, (nm*m1+n_columns))   ! w
     &     + (nv+n_columns)                 ! d_t
     &     + 2*(nv*2*m1+2*n_columns)        ! u_x, v_x, u_y, v_y
     &     + 4*(nv*local_size+n_columns)    ! u0_z, v0_z

      byte_tridi = LWORK * REAL_SIZE

!     ===== for FS ====

      call FS_WorkSize(n, LWORK, LIWORK)
      byte_dcx = LWORK * REAL_SIZE + LIWORK * INT_SIZE

!     ===== for TRBAK ====

      m    = min(nsm, m0)
      if (m < 1) m = 1

      na   = (n-1)/y_nnod+1
      na   = na  + mod(na  -1,2)
      call CSTAB_get_optdim(lda, 9, 16*4, 16*6, nm)

      LWORK = n
     &     + 3*(max(nm*m,n)+n_columns)
     &     + 4*(na*m+ns0+n_columns)
     &     + m*nm
      LIWORK = max(((m-1)/y_nnod+1)*((n-1)/x_nnod+1),n)

      byte_trbak = LWORK * REAL_SIZE + LIWORK * INT_SIZE


      nax(1) = max(byte_tridi, byte_dcx, byte_trbak)
      call MPI_Allreduce(nax(1), nax(2), 1, MPI_INTEGER8,
     &     MPI_MAX, TRD_COMM_WORLD, ierr)
      byte = nax(2)

      return

      end function FS_memory_internal

