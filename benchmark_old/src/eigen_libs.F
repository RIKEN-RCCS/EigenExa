!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_libs_mod
!
! Purpose
! =======
!
! eigen_libs_mod manages the public modules for the EigenExa library
! specificall, it combines eigen_libs0 and FS_libs, and organizes
! significant APIs
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_libs_mod

      use eigen_libs0_mod
      use FS_libs_mod
      use eigen_blacs_mod

      public


      ! interfaces for the API of external functions
      ! which enable optional arguments.

      interface

       ! eigen_sx is defined in outside of the module chain
       subroutine eigen_sx(n, nvec, a, lda, w, z, ldz,
     &        m_forward, m_backward, mode)
       integer,   intent(in)           :: n
       integer,   intent(in)           :: nvec
       real(8),   intent(inout)        :: a(1:lda, *)
       integer,   intent(in)           :: lda
       real(8),   intent(out)          :: w(1:n)
       real(8),   intent(out)          :: z(1:ldz, *)
       integer,   intent(in)           :: ldz
       integer,   intent(in), optional :: m_forward
       integer,   intent(in), optional :: m_backward
       character(*), intent(in), optional :: mode
       end subroutine eigen_sx

       ! eigen_h is defined in outside of the module chain
       subroutine eigen_h(n, nvec, a, lda, w, z, ldz,
     &        m_forward, m_backward, mode)
       integer,   intent(in)           :: n
       integer,   intent(in)           :: nvec
       complex(8),intent(inout)        :: a(1:lda, *)
       integer,   intent(in)           :: lda
       real(8),   intent(out)          :: w(1:n)
       complex(8),intent(out)          :: z(1:ldz, *)
       integer,   intent(in)           :: ldz
       integer,   intent(in), optional :: m_forward
       integer,   intent(in), optional :: m_backward
       character(*), intent(in), optional :: mode
       end subroutine eigen_h

      end interface

      public  :: eigen_init
      public  :: eigen_get_matdims
      public  :: eigen_free
      
      public  :: eigen_sx
      public  :: eigen_s

      contains


      subroutine eigen_init(comm,order)

      use mpi,only : MPI_COMM_WORLD
      use eigen_libs0_mod,only : eigen_init0
      use FS_libs_mod,only : FS_init
      implicit none

      integer, intent(in), optional :: comm
      character(*), intent(in), optional :: order

      character :: order0
      integer     :: comm0

      if (present(comm)) then
         comm0=comm
      else
         comm0=MPI_COMM_WORLD
      endif
      if (present(order)) then
         order0 = order(1:1)
      else
         order0 = 'C'
      endif
      if (order0 == 'R' .or. order0 == 'r') then
         order0 = 'R'
      else
         order0 = 'C'
      end if
#if CODE_AKASHI
      call eigen_init0(comm=comm0, order=order0)
#else
      call FS_init(comm0,order0)
#endif

      end subroutine eigen_init

      subroutine eigen_get_matdims(n, nx, ny,
     &           m_forward, m_backward,  mode)

      use eigen_libs0_mod,only : eigen_get_matdims0
      use FS_libs_mod,only : FS_get_matdims
      implicit none

      integer,   intent(in)           :: n
      integer,   intent(out)          :: nx, ny
      integer,   intent(in), optional :: m_forward
      integer,   intent(in), optional :: m_backward
      character, intent(in), optional :: mode

      integer    nx0, ny0
      integer    m_f, m_b
      character              :: mode_

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
      if (present(mode)) then
         mode_ = mode
      else
         mode_ = 'O'
      end if

#if CODE_AKASHI
      call eigen_get_matdims0(n, nx, ny, m_f, m_b, mode_)
#else
      call eigen_get_matdims0(n, nx, ny, m_f, m_b, mode_)
      call FS_get_matdims(n, nx0, ny0)
      nx=max(nx,nx0)
      ny=max(ny,ny0)
#endif

      end subroutine eigen_get_matdims

      subroutine eigen_s(n, nvec, a, lda, w, z, ldz,
     & m_forward, m_backward, mode)

      use eigen_libs0_mod,only : eigen_s0
      use FS_libs_mod,only : eigen_FS
      implicit none

      integer,   intent(in)           :: n
      integer,   intent(in), optional :: nvec
      real(8),   intent(inout)        :: a(lda,*)
      integer,   intent(in)           :: lda
      real(8),   intent(out)          :: w(*)
      real(8),   intent(out)          :: z(ldz,*)
      integer,   intent(in)           :: ldz
      integer,   intent(in), optional :: m_forward
      integer,   intent(in), optional :: m_backward
      character, intent(in), optional :: mode

      integer   :: nvec0
      integer   :: m_forward0
      integer   :: m_backward0
      character :: mode0

      if (present(nvec)) then
         nvec0=nvec
      else
         nvec0=n
      endif
      if (present(m_forward)) then
         m_forward0 = m_forward
      else
         m_forward0 = 48
      endif
      if (present(m_backward)) then
         m_backward0 = m_backward
      else
         m_backward0 = 128
      endif
      if (present(mode)) then
         mode0 = mode
      else
         mode0 = 'A'
      endif

#if CODE_AKASHI
      call eigen_s0(n, nvec0, a, lda, w, z, ldz,
     &   m_forward0, m_backward0, mode0)
#else
      call eigen_FS(n, nvec0, a, lda, w, z, ldz,
     &   m_forward0, m_backward0, mode0)
#endif

      end subroutine eigen_s

      subroutine eigen_free()

        use eigen_libs0_mod,only : eigen_free0
        use FS_libs_mod,only : FS_free
       implicit none

#if CODE_AKASHI
        call eigen_free0()
#else
        call FS_free()
#endif

       end subroutine eigen_free

      end module eigen_libs_mod
