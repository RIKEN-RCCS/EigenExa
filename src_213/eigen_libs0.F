!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_libs0_mod
!
! Purpose
! =======
!
! eigen_libs_mod manages the public modules for the EigenExa library
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_libs0_mod

      use comm_mod
      use eigen_blacs_mod
      use eigen_devel_mod
      use CSTAB_mod
      use mpi
!$    use omp_lib

      use, intrinsic :: iso_c_binding
      use, intrinsic :: ieee_arithmetic

      implicit none   

      ! default attribute is private
      private

      character(32)   ::  CodeNAME = 'EigenExa'
      type :: version_t
        integer       :: Major_Version ! Major version
        integer       :: Minor_Version ! Minor version
        integer       :: Patch_Level   ! Patchlevel 0=none,1=a,2=b,...
        character(32) :: date          ! Release date
        character(32) :: vcode         ! Version code name
      end type  version_t

      type(version_t), public, parameter :: Eigen_Version 
     &     = version_t (
     &     2, 13, 0,                   ! Major, Minor, Patchlevel
     &     'July 23, 2024',            ! Release date
#if CODE_AKASHI
     &     'tamakazura/akashi'         ! Version code
#else
     &     'tamakazura'                ! Version code
#endif
     &     )

      integer, public, parameter ::  eigen_NB    = 64
      integer, public, parameter ::  eigen_NB_f  = 48
      integer, public, parameter ::  eigen_NB_b  = 128

      ! interfaces for the API of external functions
      ! which enable optional arguments.

      interface

       ! eigen_s0 is defined outside of the module chain
       subroutine eigen_s0(n, nvec, a, lda, w, z, ldz,
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
       end subroutine  eigen_s0

      end interface

      ! interfaces for generalized API,
      ! which is called by the same funcname,
      ! but cannot relaized by optional arguments

      interface eigen_loop_info
        module procedure  eigen_loop_info_, eigen_loop_info_XY
      end interface eigen_loop_info

      interface eigen_loop_start
        module procedure eigen_loop_start_, eigen_loop_start_XY 
      end interface eigen_loop_start

      interface eigen_loop_end
        module procedure eigen_loop_end_, eigen_loop_end_XY
      end interface eigen_loop_end
     
      interface eigen_translate_l2g
        module procedure eigen_translate_l2g_, eigen_translate_l2g_XY
      end interface eigen_translate_l2g

      interface eigen_translate_g2l
        module procedure eigen_translate_g2l_, eigen_translate_g2l_XY
      end interface eigen_translate_g2l

      interface eigen_owner_node
        module procedure eigen_owner_node_, eigen_owner_node_XY
      end interface eigen_owner_node

      interface eigen_owner_index
        module procedure eigen_owner_index_, eigen_owner_index_XY
      end interface eigen_owner_index

      ! significant drivers
      public  :: eigen_s0
      public  :: eigen_init0
      public  :: eigen_free0

      ! information query
      public  :: eigen_get_version
      public  :: eigen_show_version
      public  :: eigen_initialized
      public  :: eigen_get_matdims0
      public  :: eigen_memory_internal
      public  :: eigen_get_procs
      public  :: eigen_get_id
      public  :: eigen_get_comm
      public  :: eigen_get_errinfo

      ! utility functions regarding counters/indices conversion
      public  :: eigen_loop_info
      public  :: eigen_loop_start
      public  :: eigen_loop_end
      public  :: eigen_translate_l2g
      public  :: eigen_translate_g2l
      public  :: eigen_owner_node
      public  :: eigen_owner_index
      public  :: eigen_convert_ID_xy2w
      public  :: eigen_convert_ID_w2xy
      public  :: eigen_vector_zeropad_x
      public  :: eigen_vector_zeropad_y

      ! special constants
      public  :: get_constant_eps
      public  :: get_constant_nan
      public  :: get_constant_pai
      public  :: get_constant_2pai
      public  :: get_constant_pai_2

      public  :: eigen_diag_loop_info

      ! extra routines to check NaN
      public  :: eigen_check_mat_NaN
      public  :: eigen_check_vec_NaN

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_version_
!
! Purpose
! =======
!
! Get version information, body part
!
!
! Arguments
! =========
!
! version    (output) integer
!            current version number
!
! date       (output) character(*)
!            date string
!
! vcode      (output) character(*)
!            version code string
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_version(version, date, vcode)

      integer,      intent(out) :: version
      character(*), intent(out), optional :: date
      character(*), intent(out), optional :: vcode

      version = Eigen_Version%Major_Version * 10000
     &     + Eigen_Version%Minor_Version * 100
     &     + Eigen_Version%Patch_Level

      if ( present(date) ) date = Eigen_Version%date
      if ( present(vcode) ) vcode = Eigen_Version%vcode

      return

      end subroutine eigen_get_version

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_show_version
!
! Purpose
! =======
!
! Show version information
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_show_version()

      character(256) :: fmt
      character(256) :: version
      character(1  ) :: patchlevel
      integer        :: i


      i = min(26, Eigen_Version%Patch_Level) + 1
      patchlevel = (" abcdefghijklmnopqrstuvwxyz*" (i:i))

      fmt = '(I1,A,I1,A)'
      i = int(log10(dble(Eigen_Version%Major_Version)))+1
      fmt(3:3) = ("123456789*" (i:i))
      i = int(log10(dble(Eigen_Version%Minor_Version)))+1
      fmt(8:8) = ("123456789*" (i:i))

      write(version, fmt)
     &     Eigen_Version%Major_Version,
     &     '.',Eigen_Version%Minor_Version, trim(patchlevel)

      if (TRD_inod == 1) then
        print*, "## EigenExa version (", trim(version),
     &       ") / (", trim(Eigen_Version%date),
     &       ") / (", trim(Eigen_Version%vcode), ")"
      end if

      return

      end subroutine eigen_show_version

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_initialized
!
! Purpose
! =======
!
! Get the flag that means initialized.
!
!
! Arguments
! =========
!
! flag   (output) logical
!        flag that means whether or not initialized.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_initialized(flag)

      logical, intent(out) :: flag


      call eigen_get_initialized(flag)

      return

      end subroutine eigen_initialized

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_init
!
! Purpose
! =======
!
! Initialize the utility of the EigenExa library.
!
!
! Arguments
! =========
!
! comm   (input) integer, optional
!        communicator
!
! order  (input) character(*), optional
!        Grid major
!        'R': row-major
!        'C': column-major
!
! scalapack_context  (input) integer, optional
!         This option will be supported in future version.
!
! gridmap (input) integer, optional
!         This option will be supported in a future version.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_init0(comm, order,
     &     scalapack_context, gridmap)

      integer,      intent(in), optional :: comm
      character(*), intent(in), optional :: order
      integer,      intent(in), optional :: scalapack_context
      integer,      intent(in), optional :: gridmap(:,:)

      logical          :: flag
      integer          :: n1, n2, n3, i, j, k, ierr
      integer          :: topo_type, cart_dim, dims(2), coords(2)
      logical          :: periods(2)
      integer          :: old_grp, new_grp

      integer          :: comm0
      integer          :: group0, group1, kk0(1), kk1(1)

      integer          :: local_size
!$    integer          :: th0(2), th1(2)

      real(8)          :: s1, s2, d1, d2, buff(4)
      real(8), allocatable :: buff1(:), buff2(:)

      character(1)     :: GRID_major = 'C'

!--------------------

      call eigen_timer_reset(0, 0, 0, 0)

!--------------------

      call eigen_initialized(flag)

      if (flag) then
        if (TRD_inod == 1) then
          print*,"*************"
          print*,"** CAUTION **"
          print*,"*************"
          print*,"You are going to initialize EigenExa,"
          print*,"while EigenExa was not freed at last call."
          print*,"EigenExa restarts again by itself."
        end if
        call eigen_free0()
      end if

!---- Comm setup ---

      call eigen_init_comm_setup()

!---- OpenMP thread check ---

      call eigen_init_omp_setup()

!---- Cartesian check ---

      call eigen_init_cartesian_check()

!---- BLACS setup ---

      call eigen_init_blacs_setup()
#if DEBUG
      if ( TRD_inod == 1 ) print*,"BLACS done"
#endif

!---- Preliminary Sampling Collective comms. ---

      call eigen_init_collective_comms()
#if DEBUG
      if ( TRD_inod == 1 ) print*,"Collective sampling done"
#endif

!---- Final setup ---

      call eigen_set_initialized()
#if DEBUG
      if ( TRD_inod == 1 ) print*,"initialization done"
#endif

!--------------------

      return

      contains

!--------------------

      subroutine eigen_init_comm_setup()
#if defined(__INTEL_COMPILER)
      integer :: info
#endif

      if (present(comm)) then
        comm0 = comm
      else
        comm0 = MPI_COMM_WORLD
      end if
      call MPI_Comm_Test_Inter(comm0, flag, ierr)
      if (flag .or. ierr /= MPI_SUCCESS) then
        print*,"*************"
        print*,"** CAUTION **"
        print*,"*************"
        print*,"You are going to initialize EigenExa with"
        print*,"an invalid communicator."
        print*,"EigenExa terminates this run."
        call flush(6)
        call sleep(1)
        call flush(6)
        call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
      end if
      if (comm0 == MPI_COMM_NULL) then
        TRD_COMM_WORLD = MPI_COMM_NULL
        x_COMM_WORLD   = MPI_COMM_NULL
        y_COMM_WORLD   = MPI_COMM_NULL
        z_COMM_WORLD   = MPI_COMM_NULL
        w_COMM_WORLD   = MPI_COMM_NULL
        TRD_nnod = 0; TRD_inod = 0
        x_nnod   = 0; x_inod   = 0
        y_nnod   = 0; y_inod   = 0
        z_inod   = 0; z_nnod   = 0
        w_inod   = 0; w_nnod   = 0
      else
        call MPI_Comm_dup(comm0, TRD_COMM_WORLD, ierr)
#if defined(__INTEL_COMPILER)
        call MPI_Info_create(info, ierr)
        call MPI_Info_set(info, "I_MPI_CBWR", "2", ierr)
        call MPI_Comm_set_info(TRD_COMM_WORLD, info, ierr)
        call MPI_Info_free(info, ierr)
#endif
      end if

      return

      end subroutine eigen_init_comm_setup

!--------------------

      subroutine eigen_init_omp_setup()

!$    if (TRD_COMM_WORLD /= MPI_COMM_NULL) then
!$      call MPI_Query_thread(MPI_THREAD_MODE, ierr)
!$      local_size = 1
!$OMP PARALLEL
!$OMP MASTER
!$      local_size = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$      th0(1) = local_size
!$      th0(2) = -th0(1)
!$      call MPI_Allreduce(th0, th1, 2, MPI_INTEGER, MPI_MAX,
!$   &        TRD_COMM_WORLD, ierr)
!$      j = th1(1) + th1(2)
!$      if (j /= 0) then
!$        call MPI_Barrier(TRD_COMM_WORLD, ierr)
!$        call flush(6)
!$        if (TRD_inod == 1) then
!$          print*,"*************"
!$          print*,"** CAUTION **"
!$          print*,"*************"
!$          print*,"EigenExa supports only homogeneous thread"
!$   &           // " setting !"
!$          print*,"EigenExa terminates this run."
!$          print*,"*************"
!$        end if
!$        do j=1,2
!$          call MPI_Barrier(TRD_COMM_WORLD, ierr)
!$          call flush(6)
!$        end do
!$        call MPI_Barrier(TRD_COMM_WORLD, ierr)
!$        call flush(6)
!$        call sleep(1)
!$        call flush(6)
!$        call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
!$      end if
!$    end if

      return

      end subroutine eigen_init_omp_setup

!--------------------

      subroutine eigen_init_cartesian_check()
      integer :: ierr

      if (TRD_COMM_WORLD /= MPI_COMM_NULL) then

!---- Cartesian check ---

        call MPI_Comm_size(TRD_COMM_WORLD, TRD_nnod, ierr)
        call MPI_Comm_rank(TRD_COMM_WORLD, TRD_inod, ierr)
        TRD_inod = TRD_inod+1

        call MPI_Topo_test(TRD_COMM_WORLD, topo_type, ierr)
        if (topo_type == MPI_CART) then
          call MPI_Cartdim_get(TRD_COMM_WORLD, cart_dim, ierr)
        else
          cart_dim = 1
        end if

!---- 2D Cartesian setup ---

        if (cart_dim == 2) then

          call MPI_Cart_get(TRD_COMM_WORLD, cart_dim,
     &         dims, periods, coords, ierr)
          x_nnod = dims(1)
          y_nnod = dims(2)
          x_inod = coords(1) +1 ! fixed
          y_inod = coords(2) +1 ! fixed

          GRID_major = 'R'
          if (present(order)) then
            if (x_inod == 1 .and. y_inod == 1) then
              if (order(1:1) == 'C' .or. order(1:1) == 'c') then
                print*,"*************"
                print*,"** CAUTION **"
                print*,"*************"
                print*,"The MPI_CART you specified is based on"
     &               // " R-major,"
                print*,"but you also specified C-major option."
                print*,"EigenExa solve this conflict by taking"
     &               // " R-major,"
              end if
            end if
          end if

        else                    ! if (cart_dim == 2

!---- Setup 2D process map ---

          x_nnod = int(sqrt(dble(TRD_nnod)))
          i = 1                 ! minimum factor, x_nnod must be
!     multiple of k
          if (mod(TRD_nnod, i) == 0) then
            k = i
          else
            k = 1
          end if
          do
            if (x_nnod <= k) exit
            if (mod(x_nnod, k) == 0 .and.
     &          mod(TRD_nnod, x_nnod) == 0) exit
            x_nnod = x_nnod-1
          end do                !!
          y_nnod = TRD_nnod/x_nnod

          if (present(order)) then
            GRID_major = order(1:1)
          else
            GRID_major = 'C'
          end if

          if (GRID_major == 'R' .or. GRID_major == 'r') then
            GRID_major = 'R'
          else
            GRID_major = 'C'
          end if

          call eigen_set_grid_major(GRID_major)

          if (GRID_major == 'R') then
!     row-major
            x_inod =    (TRD_inod-1)/y_nnod +1
            y_inod = mod(TRD_inod-1, y_nnod)+1
          else
!     column-major
!     ** EigenK adopts column-major in default
!     The process ordering on the Earth Simulator was done in the way
!     such as incrmental from internal processors to external nodes.
!     We want to make the processes, which possess a distributed vector,
!     close physically on the network connection.
!     These are the historical reason to adopt column-major.
            x_inod = mod(TRD_inod-1, x_nnod)+1
            y_inod =    (TRD_inod-1)/x_nnod +1
          end if

        end if                  ! if (cart_dim == 2) then

!---- x-Comm & y-Comm setup ---

#if DEBUG
        if ( TRD_inod==1 ) print*,"Comm split[X]"
#endif
        call MPI_Comm_split(TRD_COMM_WORLD, y_inod, x_inod,
     &       x_COMM_WORLD, ierr)
#if DEBUG
        if ( TRD_inod==1 ) print*,"Comm split[Y]"
#endif
        call MPI_Comm_split(TRD_COMM_WORLD, x_inod, y_inod,
     &       y_COMM_WORLD, ierr)

!---- Parameter setup for datacast_dbl (row-column vector redistribution) ---

      if (x_nnod /= y_nnod) then

        n1 = max(x_nnod, y_nnod)
        n2 = min(x_nnod, y_nnod)
        do
          if (n1 == n2) then
            n_common = n1
            exit
          end if
          n3 = n1-n2
          n1 = max(n2, n3)
          n2 = min(n2, n3)
        end do                  !!

        if (associated(p0_)) then
          deallocate(p0_)
        end if
        allocate(p0_(1:max(x_nnod, y_nnod)), stat=ierr)
        if (ierr /= 0) then
          call eigen_abort( "Memory allocation error [init].", 1 )
        end if

        if (associated(q0_)) then
          deallocate(q0_)
        end if
        allocate(q0_(1:max(x_nnod, y_nnod)), stat=ierr)
        if (ierr /= 0) then
          call eigen_abort( "Memory allocation error [init].", 1 )
        end if

        p0_(:) = -1
        q0_(:) = -1
        do i = 1, x_nnod
          if (mod(i-1, n_common) == mod(y_inod-1, n_common)) then
            n1 = y_inod-i
            if (n1 >= 0) then
              do j = 1, x_nnod
                k = +n1+(j-1)*y_nnod
                if (mod(k, x_nnod) == 0) then
                  p0_(i) = k/x_nnod
                  q0_(i) = (j-1)
                  exit
                end if
              end do            ! j
            else
              do j = 1, y_nnod
                k = -n1+(j-1)*x_nnod
                if (mod(k, y_nnod) == 0) then
                  q0_(i) = k/y_nnod
                  p0_(i) = (j-1)
                  exit
                end if
              end do            ! j
            end if
          end if
        end do                  ! i
        p0_(:) = p0_(:)+1
        q0_(:) = q0_(:)+1


        diag_0 = 0
        diag_1 = 0
        do i = 1, y_nnod/n_common
          j = (i-1)*y_nnod+y_inod
          k = mod(j-1, x_nnod)+1
          if (k == x_inod) then
            diag_0 = i
            diag_1 = (j-1)/x_nnod+1
            exit
          end if
        end do                  ! i_1

#if DEBUG
        if ( TRD_inod==1 ) print*,"p0_, q0_, diag_ ", n_common
#endif

      else

        n_common = x_nnod

        if (associated(p0_)) then
          deallocate(p0_)
        end if
        if (associated(q0_)) then
          deallocate(q0_)
        end if

        if (y_inod == x_inod) then
          diag_0 = 1
          diag_1 = 1
        else
          diag_0 = 0
          diag_1 = 0
        end if
          
      end if

!---- extra-Comm setup for datacast_dbl ---

        if (n_common > 1) then
        if (x_nnod /= y_nnod) then
          w_inod =     (x_inod-1)/ n_common +1
          z_inod = mod((x_inod-1), n_common)+1
          w_nnod = x_nnod/n_common
          z_nnod = n_common
#if DEBUG
        if ( TRD_inod==1 ) print*,"Split[Z]"
#endif
          call MPI_Comm_split(x_COMM_WORLD, w_inod-1, z_inod-1,
     &         z_COMM_WORLD, ierr)
#if DEBUG
        if ( TRD_inod==1 ) print*,"Split[W]"
#endif
          call MPI_Comm_split(x_COMM_WORLD, z_inod-1, w_inod-1,
     &         w_COMM_WORLD, ierr)
        else
          w_COMM_WORLD = MPI_COMM_SELF
          z_COMM_WORLD = x_COMM_WORLD
          w_inod = 1; w_nnod = 1
          z_inod = x_inod; z_nnod = x_nnod
        end if
        else
          z_COMM_WORLD = MPI_COMM_SELF
          w_COMM_WORLD = x_COMM_WORLD
          z_inod = 1; z_nnod = 1
          w_inod = x_inod; w_nnod = x_nnod
        end if

!---- Exceptional setup ---

      else                      !  if (TRD_COMM_WORLD /= MPI_COMM_NULL

        x_COMM_WORLD = MPI_COMM_NULL
        y_COMM_WORLD = MPI_COMM_NULL
        z_COMM_WORLD = MPI_COMM_NULL
        w_COMM_WORLD = MPI_COMM_NULL
        x_nnod = 0; x_inod = 0
        y_nnod = 0; y_inod = 0
        z_inod = 0; z_nnod = 0
        w_inod = 0; w_nnod = 0

        GRID_major = 'C'
        call eigen_set_grid_major(GRID_major)

      end if                    !  if (TRD_COMM_WORLD /= MPI_COMM_NULL

#if DEBUG
        if ( TRD_inod==1 ) print*,"Cartesian check is OK"
#endif

      return

      end subroutine eigen_init_cartesian_check

!--------------------

      subroutine eigen_init_blacs_setup()

#if TIMER_PRINT>1
      if (TRD_inod <= 1) then
        print*,"GRID major ",GRID_major," is specified."
      end if
#endif

      if (TRD_COMM_WORLD == MPI_COMM_NULL) then
        comm0  = MPI_COMM_SELF
        x_nnod = 1; x_inod = 1
        y_nnod = 1; y_inod = 1
      else
        comm0  = TRD_COMM_WORLD
      end if

      call Eigen_BLACS_Init(comm0, x_nnod, y_nnod, GRID_Major)

      if (TRD_COMM_WORLD == MPI_COMM_NULL) then
        x_nnod = 0; x_inod = 0
        y_nnod = 0; y_inod = 0
      end if

      return

      end subroutine eigen_init_blacs_setup

!--------------------

      subroutine eigen_init_collective_comms()
      integer :: ierr

      call repro_check()

      if (TRD_COMM_WORLD /= MPI_COMM_NULL) then

!---- Measure the overhead of collective communications ---

        call MPI_Barrier(TRD_COMM_WORLD,ierr)
        call MPI_Barrier(TRD_COMM_WORLD,ierr)

!$OMP PARALLEL
!$OMP+          PRIVATE(i,j,d1,d2,buff1,buff2)

        allocate(buff1(2048), buff2(2048), stat=ierr)
        if (ierr /= 0) then
          call eigen_abort( "Memory allocation error [init].", 1 )
        end if

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do  i = 1, 10
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        call MPI_Barrier(TRD_COMM_WORLD,ierr)
        call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        call MPI_Barrier(TRD_COMM_WORLD,ierr)
        call MPI_Barrier(TRD_COMM_WORLD,ierr)

        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Barrier_Overhead_x = buff(3) / TRD_nnod
        Barrier_Overhead_y = buff(4) / TRD_nnod
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),2048,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),2048,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Reduce_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
        Reduce_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call reduce_dbl(buff2(1),buff1(1),1024,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call reduce_dbl(buff2(1),buff1(1),1024,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Reduce_cont_Overhead_x = buff(3) / TRD_nnod
     &       - Barrier_Overhead_x - Reduce_Overhead_x
        Reduce_cont_Overhead_y = buff(4) / TRD_nnod
     &       - Barrier_Overhead_y - Reduce_Overhead_y
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),2048,1,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),2048,1,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, TRD_COMM_WORLD, ierr)
        Bcast_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
        Bcast_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call bcast_dbl(buff2(1),1024,1,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do  i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call bcast_dbl(buff2(1),1024,1,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Bcast_cont_Overhead_x = buff(3) / TRD_nnod
     &       - Barrier_Overhead_x - Bcast_Overhead_x
        Bcast_cont_Overhead_y = buff(4) / TRD_nnod
     &       - Barrier_Overhead_y - Bcast_Overhead_y
!$OMP END MASTER

        deallocate(buff1, buff2)

!$OMP END PARALLEL

!--------------------

        if (Reduce_cont_Overhead_x < ZERO) then
          Reduce_cont_Overhead_x = ZERO
        end if
        if (Reduce_cont_Overhead_y < ZERO) then
          Reduce_cont_Overhead_y = ZERO
        end if
        if (Bcast_cont_Overhead_x  < ZERO) then
          Bcast_cont_Overhead_x  = ZERO
        end if
        if (Bcast_cont_Overhead_y  < ZERO) then
          Bcast_cont_Overhead_y  = ZERO
        end if

#if TIMER_PRINT>0
        if (TRD_inod == 1) then
          print*,";----- MPI in-(x/y)-group overhead -----"
          print*,"Overhead Barrier(x) :: ", Barrier_Overhead_x
          print*,"Overhead Barrier(y) :: ", Barrier_Overhead_y
          print*,"Overhead Bcast(x)   :: ", Bcast_Overhead_x
          print*,"Overhead Bcast(y)   :: ", Bcast_Overhead_y
          print*,"Overhead Reduce(x)  :: ", Reduce_Overhead_x
          print*,"Overhead Reduce(y)  :: ", Reduce_Overhead_y
          print*,"Overhead Bcast_cont(x)   :: ",
     &         Bcast_cont_Overhead_x
          print*,"Overhead Bcast_cont(y)   :: ",
     &         Bcast_cont_Overhead_y
          print*,"Overhead Reduce_cont(x)  :: ",
     &         Reduce_cont_Overhead_x
          print*,"Overhead Reduce_cont(y)  :: ",
     &         Reduce_cont_Overhead_y
        end if
#endif

!--------------------

      end if

      return

      end subroutine  eigen_init_collective_comms

      subroutine repro_check()
      character(128) :: val
      integer        :: tmp(2), res(2)

      repro_reduce = .TRUE.

#if defined(__INTEL_COMPILER)
      call getenv("I_MPI_ADJUST_ALLREDUCE", val)
      if ( val(1:1) == '4' .or. val(1:1) == '6' ) then
        repro_reduce = .FALSE.
      end if

      call getenv("I_MPI_ADJUST_REDUCE", val)
      if ( val(1:1) == '3' .or. val(1:1) == '4'
     &      .or. val(1:1) == '6' ) then
        repro_reduce = .FALSE.
      end if

      if ( repro_reduce ) then
        tmp(1) =  1
        tmp(2) = -1
      else
        tmp(1) = 0
        tmp(2) = 0
      end if
      call MPI_Allreduce(tmp, res, 2, MPI_INTEGER, MPI_MAX,
     &     TRD_COMM_WORLD, ierr)
      if ( res(1)*res(2) == -1 ) then
         repro_reduce = .TRUE.
      else
         repro_reduce = .FALSE.
      end if
#endif

      return

      end subroutine  repro_check

      end subroutine  eigen_init0

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_free
!
! Purpose
! =======
!
! Free the utility of the EigenExa library
!
!
! Arguments
! =========
!
! flag   (input) integer, optional
!        TRUE:  prit out the stats of timer table.
!        FALSE: 
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_free0(flag)

      integer, intent(in), optional ::  flag

      integer                ::  ierr
      logical                ::  local_flag
      real(8)                ::  tmp


      call eigen_get_initialized(local_flag)

      if (.not. local_flag) then
        return
      end if

      if (TRD_COMM_WORLD /= MPI_COMM_NULL) then

        if (n_common > 1) then
        if (x_nnod /= y_nnod) then
          call MPI_Comm_free(z_COMM_WORLD, ierr)
          call MPI_Comm_free(w_COMM_WORLD, ierr)
        end if
        end if

        call MPI_Comm_free(x_COMM_WORLD, ierr)
        call MPI_Comm_free(y_COMM_WORLD, ierr)

        call MPI_Comm_free(TRD_COMM_WORLD, ierr)

        if (x_nnod /= y_nnod) then
          deallocate(p0_)
          deallocate(q0_)
        end if

        if (present(flag)) then
          if (flag == 1) then
            tmp = eigen_timer_print('EigenExa(finalized)')
          end if
        end if
        call eigen_timer_reset(0, 0, 0, 0)

      end if

!     call Eigen_BLACS_Exit()

      TRD_COMM_WORLD = MPI_COMM_WORLD

      call eigen_unset_initialized()

      return

      end subroutine eigen_free0

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_matdims0
!
! Purpose
! =======
!
! Returns the recommended array size for the distributed matrix,
! such as (1:nx, 1:ny).
!
!
! Arguments
! =========
!
! n      (input) integer
!        Global matrix dimension
!
! nx     (output) integer
!        Array size for the first index
!
! ny     (output) integer
!        Array size for the second index
!
! m_forward  (input) integer
! m_backward (input) integer
!        Block size to be specified in calling eigenvalue solver
!
! mode   (input) character
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_matdims0(n, nx, ny,
     &             m_forward, m_backward, mode)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_get_matdims0

      integer, intent(in)    :: n
      integer, intent(out)   :: nx, ny
      integer, intent(in), optional :: m_forward
      integer, intent(in), optional :: m_backward
      character, intent(in), optional :: mode

      integer :: NPROW, NPCOL, NB
      integer :: n1, nm, nmz, nmw, nn, larray

      integer :: m_f, m_b
      integer :: lddw, lddz, nxx

      integer(8) :: LX1, LX2
      character              :: mode_

      if (n <= 0) then
        nx = -1; ny = -1
        return
      end if

      if (present(mode)) then
         mode_ = mode
      else
         mode_ = 'O'
      end if
      If ( mode_ == 'M' ) then
         nx = (n-1)/x_nnod+1
         ny = (n-1)/y_nnod+1
         return
      endif
      if ( mode_ == 'L' ) then
         nx =(n-1)/x_nnod+1
         nx =((nx-1)/32+1)*32
         ny =(n-1)/y_nnod+1
         return
      endif
!
!     for mode == 'O' : optimized for performance
!
      NPROW = x_nnod
      NPCOL = y_nnod

      n1 = ((n-1)/NPROW+1)
      call CSTAB_get_optdim(n1, 6, 16*4, 16*4*2, nm)

      if (present(m_forward)) then
        m_f  = m_forward
      else
        m_f  = eigen_NB_f
      end if
      if (present(m_backward)) then
        m_b  = m_backward
      else
        m_b  = eigen_NB_b
      end if
#if 0
      NB = max( m_f, m_b, eigen_NB )
#else
      NB = max( m_b, eigen_NB )
#endif

      nmz = ((n-1)/NPROW+1)
      nmz = ((nmz-1)/NB+1)*NB+1
      nn  = nmz
      nmz = (n-1)/NB+1
      nmz = ((nmz-1)/NPROW+1)*NB
!     Fix on version 2.2b
!     to avoid unexpected SIGSEGV,
!     use the maximum of nn and nmz.
      nmz = max(nn, nmz)

      nmw = ((n-1)/NPCOL+1)
      nmw = ((nmw-1)/NB+1)*NB+1
      nn  = nmw
      nmw = (n-1)/NB+1
      nmw = ((nmw-1)/NPCOL+1)*NB
!     Fix on version 2.2b
!     to avoid unexpected SIGSEGV,
!     use the maximum of nn and nmw.
      nmw = max(nn, nmw)


      larray = max(nmz, nm)*nmw

      nx = nm
      ny = (larray-1)/nm+1

#if 1
!     
!     check whether intermediate datasize is valid with in integer(4).
!     
      NB = eigen_NB
      NB = min(NB, N)

      lddz = (n-1)/min(NPCOL,NPROW)+1
      lddz = ((lddz-1)/NB+1)*NB
      nxx  = (N-1)/min(NPCOL,NPROW)+1

      LX1 = INT(lddz,8); LX1 = LX1 * lddz
      LX2 = INT(nxx,8);  LX2 = LX2 * nxx

      if (LX1 >= (2_8)**31 .OR.
     &    LX2 >= (2_8)**31) then
        print*, "Warning :: oversized problem !!"
        larray = -1
        nx     = -1
        ny     = -1
      end if
#endif

      return

!      end subroutine  eigen_get_matdims
      end subroutine  eigen_get_matdims0

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function : eigen_memory_internal
!
! Purpose
! =======
!
! Returns size of the required memory internally while proceeding the
! EigenExa library.
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*
!
!     ******* ! important ! API has been changed ! *********************
!     Since the size of intermediate data may overfolow integer(4),
!     version 2.3b has changed the datatype of the function  result
!     from integer to integer(8).
!     ******************************************************************
!
      integer(8) function eigen_memory_internal(n, lda, ldz,
     &     m1_opt, m0_opt)
     &     result( byte )
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_memory_internal

      integer, intent(in)           :: n, lda, ldz
      integer, intent(in), optional :: m1_opt, m0_opt

      integer, parameter :: nm_max_L1 = 16*4
      integer, parameter :: nm_max_L2 = 16*6
      integer(8), parameter :: REAL_SIZE = 8_8 ! must be in integer(8)
      integer(8), parameter :: INT_SIZE  = 4_8 ! must be in integer(8)

      integer    :: m1, m0
      integer    :: MYROW, MYCOL, NPROW, NPCOL, NP, NQ, NB
      integer    :: local_rank, local_size
      integer    :: m , nv, nm, ierr

      integer(8) :: lddz, lddw, lwk
      integer(8) :: LWORK, LIWORK
      integer    :: na, nx, nxv, ny, nz, kx
      integer(8) :: l_array_1, l_array_2

      integer(8) :: nax(2)
      integer(8) :: padding
      integer(8) :: byte_tridi
      integer(8) :: byte_dcx
      integer(8) :: byte_trbak

      integer, external  :: NUMROC

      include 'CSTAB.h'


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
        m1 = eigen_NB_f
      end if

      if (present(m0_opt)) then
        m0 = m0_opt
      else
        m0 = eigen_NB_b
      end if

#if BOOST_BY_CACHE_ALIGNMENT
      padding = n_columns
#else
      padding = 0
#endif

!     ===== for TRIDI ====

      m    = min(m1, n)
      if (m < 1) m = 1

      nx = (n-1)/x_nnod+1 +2
      nv = lda
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

      LWORK =
     &     + (lda*m           +  padding)  ! w for prd
     &     + (l_array_1       +  padding)  ! u_t
     &     + (l_array_2       +  padding)  ! v_t
     &     + (nv              +2*padding)  ! d_t for prd
     &     + 2*(nv*2*m1       +2*padding)  ! ux_vx & uy_vy
     &     + 3*(m)                         ! uvw_nn for trd
     &     + 4*(nv*local_size +  padding)  ! u0_z, v0_z
      LIWORK = 0

      byte_tridi = LWORK * REAL_SIZE + LIWORK * INT_SIZE

!     ===== for DCX ====

      NB = eigen_NB
      NB = min(NB, N)

      NPROW = x_nnod
      NPCOL = y_nnod
      MYROW = x_inod-1
      MYCOL = y_inod-1
      NP = NUMROC(n, NB, MYROW, 0, NPROW)
      NQ = NUMROC(n, NB, MYCOL, 0, NPCOL)

      lddz = (n-1)/NPROW+1
      lddz = ((lddz-1)/NB+1)*NB
      lddw = (n-1)/NPCOL+1
      lddw = ((lddw-1)/NB+1)*NB

      nx     = (N-1)/NPCOL+1
      nxv    = (N-1)/NPCOL+1
      LWORK  = max(1+6*N+2*NP*(NQ+max(NQ,NB)), lddz*lddw, ldz*nx)
     &       + 2*n*local_size
      LIWORK = 2+7*n+8*NPCOL

      byte_dcx = LWORK * REAL_SIZE + LIWORK * INT_SIZE

!     ===== for TRBAK ====

      m    = min(nsm, m0, n)
      if (m < 1) m = 1

      na   = (n-1)/y_nnod+1
      na   = na  + mod(na  -1,2)
      call CSTAB_get_optdim(lda, 9, 16*4, 16*6, nm)

      nx = min(mod(n-1,m)+1,n)
      lwk=int((m-1)/y_nnod+1,8)*int((n-1)/x_nnod+1, 8)
      lwk=max(lwk, int(n,8))

      LWORK =
     &     + 3*(max(nm*m,n)+padding)     ! v1, v2, v3
     &     + 3*(na*m+ns0   +padding)     ! ss1, ss2, ss3
     &     + (na*m+ns0     +padding)     ! tt
     &     + m*nm                        ! trbk_buf
     &     + lwk                         ! wk
      LIWORK = 0

      byte_trbak = LWORK * REAL_SIZE + LIWORK * INT_SIZE


      nax(1) = max(byte_tridi, byte_dcx, byte_trbak)
      call MPI_Allreduce(nax(1), nax(2), 1, MPI_INTEGER8,
     &     MPI_MAX, TRD_COMM_WORLD, ierr)

      byte = nax(2) + 5*n*REAL_SIZE      ! d, e, e2

      return

      end function eigen_memory_internal

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_procs
!
! Purpose
! =======
!
! Returns the number of processes corresponding to three comms.
!
!
! Arguments
! =========
!
! procs  (output) integer
!        Num of Global communicator
!
! x_procs (output) integer
!        Num of x- communicator
!
! y_procs (output) integer
!        Num of y- communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_procs(procs, x_procs, y_procs)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_get_procs

      integer, intent(out)   ::  procs
      integer, intent(out)   ::  x_procs
      integer, intent(out)   ::  y_procs


      procs   = TRD_nnod
      x_procs = x_nnod
      y_procs = y_nnod

      return

      end subroutine  eigen_get_procs

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_id
!
! Purpose
! =======
!
! Returns process ID corresponding to three comms.
!
!
! Arguments
! =========
!
! id     (output) integer
!        Global ID
!
! x_id   (output) integer
!        x-ID
!
! y_id   (output) integer
!        y-ID
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_get_id(id, x_id, y_id)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_get_id

      integer, intent(out)   ::  id
      integer, intent(out)   ::  x_id
      integer, intent(out)   ::  y_id


      id   = TRD_inod
      x_id = x_inod
      y_id = y_inod

      return

      end subroutine  eigen_get_id

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_comm
!
! Purpose
! =======
!
! Returns communicators.
!
!
! Arguments
! =========
!
! comm   (output) integer
!        global communicator
!
! x_comm (output) integer
!        x-dir communicator, column-group
!
! y_comm (output) integer
!        y-dir communicator, row-group
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_comm(comm, x_comm, y_comm)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_get_comm

      integer, intent(out)   ::  comm
      integer, intent(out)   ::  x_comm
      integer, intent(out)   ::  y_comm


      comm   = TRD_COMM_WORLD
      x_comm = x_COMM_WORLD
      y_comm = y_COMM_WORLD

      return

      end subroutine eigen_get_comm

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_errinfo
!
! Purpose
! =======
!
! Returns error information.
!
!
! Arguments
! =========
!
! info   (output) integer
!        local error status
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_errinfo(info)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_get_errinfo

      integer, intent(out)   ::  info

      info   = ERROR_INFO

      return

      end subroutine eigen_get_errinfo

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_loop_info
!
! Purpose
! =======
!
! Return the local loop info, such as the initial and termination values
! of the specific global loop interval
!
!
! Arguments
! =========
!
! lstart (output) integer
!        local loop initial value
!
! lend (output) integer
!        local loop termination value
!
! istart (input) integer
!        global loop initial value
!
! iend (input) integer
!        global loop termination value
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'X' or 'Y' (not case sensitive)
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_loop_info_(istart, iend, lstart, lend, nnod, inod
     & )
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_loop_info_

      integer, intent(in)    ::  istart
      integer, intent(in)    ::  iend
      integer, intent(out)   ::  lstart
      integer, intent(out)   ::  lend
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      lstart = eigen_loop_start( istart, nnod, inod )
      lend   = eigen_loop_end  ( iend  , nnod, inod )

      return

      end subroutine eigen_loop_info_

      subroutine eigen_loop_info_XY(istart, iend, lstart, lend, pdir,
     & inod)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_loop_info_XY

      integer, intent(in)           ::  istart
      integer, intent(in)           ::  iend
      integer, intent(out)          ::  lstart
      integer, intent(out)          ::  lend
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      lstart = eigen_loop_start( istart, pdir, inod )
      lend   = eigen_loop_end  ( iend  , pdir, inod )

      return

      end subroutine eigen_loop_info_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_loop_start
!
! Purpose
! =======
!
! Return the local loop initial value of the specific global loop
! interval
!
!
! Arguments
! =========
!
! istart (input) integer
!        global loop initial value
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'W', 'X' or 'Y' (not case sensitive)
!        'W' : TRD_nnod, TRD_inod
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_loop_start_(istart, nnod, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_loop_start_

      integer, intent(in)    ::  istart
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod


      ret = (istart + nnod - 1 - inod) / nnod + 1

      return

      end function eigen_loop_start_

      integer function eigen_loop_start_XY(istart, pdir, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_loop_start_XY

      integer, intent(in)           ::  istart
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      select case (pdir(1:1))
      case( 'W', 'w', 'T', 't' )
        if ( present(inod) ) then
          ret = eigen_loop_start_(istart, TRD_nnod, inod)
        else
          ret = eigen_loop_start_(istart, TRD_nnod, TRD_inod)
        end if
      case( 'X', 'x', 'R', 'r' )
        if ( present(inod) ) then
          ret = eigen_loop_start_(istart, x_nnod, inod)
        else
          ret = eigen_loop_start_(istart, x_nnod, x_inod)
        end if
      case( 'Y', 'y', 'C', 'c' )
        if ( present(inod) ) then
          ret = eigen_loop_start_(istart, y_nnod, inod)
        else
          ret = eigen_loop_start_(istart, y_nnod, y_inod)
        end if
      case default
        ret = 0
      end select

      return

      end function eigen_loop_start_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_loop_end
!
! Purpose
! =======
!
! Return the local loop temination value of the specific global loop
! interval
!
!
! Arguments
! =========
!
! iend   (input) integer
!        global loop termination value
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'W', 'X' or 'Y' (not case sensitive)
!        'W' : TRD_nnod, TRD_inod
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_loop_end_(iend, nnod, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_loop_end_

      integer, intent(in)    ::  iend
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = (iend + nnod - 0 - inod) / nnod + 0

      return

      end function eigen_loop_end_

      integer function eigen_loop_end_XY(iend, pdir, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_loop_end_XY

      integer, intent(in)           ::  iend
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      select case (pdir(1:1))
      case( 'W', 'w', 'T', 't' )
        if ( present(inod) ) then
          ret = eigen_loop_end_(iend, TRD_nnod, inod)
        else
          ret = eigen_loop_end_(iend, TRD_nnod, TRD_inod)
        end if
      case( 'X', 'x', 'R', 'r' )
        if ( present(inod) ) then
          ret = eigen_loop_end_(iend, x_nnod, inod)
        else
          ret = eigen_loop_end_(iend, x_nnod, x_inod)
        end if
      case( 'Y', 'y', 'C', 'c' )
        if ( present(inod) ) then
          ret = eigen_loop_end_(iend, y_nnod, inod)
        else
          ret = eigen_loop_end_(iend, y_nnod, y_inod)
        end if
      case default
        ret = -1
      end select

      return

      end function eigen_loop_end_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_translate_l2g
!
! Purpose
! =======
!
! Translate the local index to the global index
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'W', 'X' or 'Y' (not case sensitive)
!        'W' : TRD_nnod, TRD_inod
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_translate_l2g_(ictr, nnod, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_translate_l2g_

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = (ictr-1)*nnod+inod

      return

      end function eigen_translate_l2g_

      integer function eigen_translate_l2g_XY(ictr, pdir, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_translate_l2g_XY

      integer, intent(in)           ::  ictr
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      select case (pdir(1:1))
      case( 'W', 'w', 'T', 't' )
        if ( present(inod) ) then
          ret = eigen_translate_l2g_(ictr, TRD_nnod, inod)
        else
          ret = eigen_translate_l2g_(ictr, TRD_nnod, TRD_inod)
        end if
      case( 'X', 'x', 'R', 'r' )
        if ( present(inod) ) then
          ret = eigen_translate_l2g_(ictr, x_nnod, inod)
        else
          ret = eigen_translate_l2g_(ictr, x_nnod, x_inod)
        end if
      case( 'Y', 'y', 'C', 'c' )
        if ( present(inod) ) then
          ret = eigen_translate_l2g_(ictr, y_nnod, inod)
        else
          ret = eigen_translate_l2g_(ictr, y_nnod, y_inod)
        end if
      case default
        ret = -1
      end select

      return

      end function eigen_translate_l2g_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_translate_g2l
!
! Purpose
! =======
!
! Translate the global index to the local index
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'W', 'X' or 'Y' (not case sensitive)
!        'W' : TRD_nnod, TRD_inod
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_translate_g2l_(ictr, nnod, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_translate_g2l_

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = (ictr-1)/nnod+1

      return

      end function eigen_translate_g2l_

      integer function eigen_translate_g2l_XY(ictr, pdir, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_translate_g2l_XY

      integer, intent(in)           ::  ictr
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      select case (pdir(1:1))
      case( 'W', 'w', 'T', 't' )
        if ( present(inod) ) then
          ret = eigen_translate_g2l_(ictr, TRD_nnod, inod)
        else
          ret = eigen_translate_g2l_(ictr, TRD_nnod, TRD_inod)
        end if
      case( 'X', 'x', 'R', 'r' )
        if ( present(inod) ) then
          ret = eigen_translate_g2l_(ictr, x_nnod, inod)
        else
          ret = eigen_translate_g2l_(ictr, x_nnod, x_inod)
        end if
      case( 'Y', 'y', 'C', 'c' )
        if ( present(inod) ) then
          ret = eigen_translate_g2l_(ictr, y_nnod, inod)
        else
          ret = eigen_translate_g2l_(ictr, y_nnod, y_inod)
        end if
      case default
        ret = -1
      end select

      return

      end function eigen_translate_g2l_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_owner_node
!
! Purpose
! =======
!
! Returns the owner ID of the specific global index
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'W', 'X' or 'Y' (not case sensitive)
!        'W' : TRD_nnod, TRD_inod
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_owner_node_(ictr, nnod, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_owner_node_

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = mod(ictr-1, nnod)+1

      return

      end function eigen_owner_node_

      integer function eigen_owner_node_XY(ictr, pdir, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_owner_node_XY

      integer, intent(in)           ::  ictr
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      select case (pdir(1:1))
      case( 'W', 'w', 'T', 't' )
        if ( present(inod) ) then
          ret = eigen_owner_node_(ictr, TRD_nnod, inod)
        else
          ret = eigen_owner_node_(ictr, TRD_nnod, TRD_inod)
        end if
      case( 'X', 'x', 'R', 'r' )
        if ( present(inod) ) then
          ret = eigen_owner_node_(ictr, x_nnod, inod)
        else
          ret = eigen_owner_node_(ictr, x_nnod, x_inod)
        end if
      case( 'Y', 'y', 'C', 'c' )
        if ( present(inod) ) then
          ret = eigen_owner_node_(ictr, y_nnod, inod)
        else
          ret = eigen_owner_node_(ictr, y_nnod, y_inod)
        end if
      case default
        ret = -1
      end select

      return

      end function eigen_owner_node_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_owner_index
!
! Purpose
! =======
!
! Returns the index if owner but -1 if not owner.
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
! OR
!
! pdir   (input) character
!        'W', 'X' or 'Y' (not case sensitive)
!        'W' : TRD_nnod, TRD_inod
!        'X' : x_nnod, x_inod
!        'Y' : y_nnod, y_inod
!
! inod   (input) integer, optional
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_owner_index_(ictr, nnod, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_owner_index_

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      integer  ::  j_2, j_3

      j_2 = eigen_loop_start_(ictr, nnod, inod)
      j_3 = eigen_loop_end_  (ictr, nnod, inod)
      if (j_2 == j_3) then
        ret = j_2
      else
        ret = -1
      end if

      return

      end function eigen_owner_index_

      integer function eigen_owner_index_XY(ictr, pdir, inod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_owner_index_XY

      integer, intent(in)           ::  ictr
      character(*), intent(in)      ::  pdir
      integer, intent(in), optional ::  inod

      select case (pdir(1:1))
      case( 'W', 'w', 'T', 't' )
        if ( present(inod) ) then
          ret = eigen_owner_index_(ictr, TRD_nnod, inod)
        else
          ret = eigen_owner_index_(ictr, TRD_nnod, TRD_inod)
        end if
      case( 'X', 'x', 'R', 'r' )
        if ( present(inod) ) then
          ret = eigen_owner_index_(ictr, x_nnod, inod)
        else
          ret = eigen_owner_index_(ictr, x_nnod, x_inod)
        end if
      case( 'Y', 'y', 'C', 'c' )
        if ( present(inod) ) then
          ret = eigen_owner_index_(ictr, y_nnod, inod)
        else
          ret = eigen_owner_index_(ictr, y_nnod, y_inod)
        endif
      case default
        ret = -1
      end select

      return

      end function eigen_owner_index_XY

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_convert_ID_xy2w
!
! Purpose
! =======
!
! Convert the 2D rank ID into the 1D rank ID.
!
!
! Arguments
! =========
!
! xinod  (input) integer
!        Local process rank ID in x_comm
!
! yinod  (input) integer
!        Local process rank ID in y_comm
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_convert_ID_xy2w(xinod, yinod)
     &     result(ret)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_convert_ID_xy2w

      integer, intent(in)    ::  xinod
      integer, intent(in)    ::  yinod

      character(1)           ::  order

      call eigen_get_grid_major( order )
      if ( order(1:1) == 'R' ) then
        ret = xinod * y_nnod + yinod
      else
        ret = yinod * x_nnod + xinod
      end if

      return

      end function eigen_convert_ID_xy2w

      subroutine eigen_convert_ID_w2xy(inod, xinod, yinod)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_convert_ID_w2xy

      integer, intent(in)    ::  inod
      integer, intent(out)   ::  xinod
      integer, intent(out)   ::  yinod

      character(1)           ::  order

      call eigen_get_grid_major( order )
      if ( order(1:1) == 'R' ) then
        xinod =    (inod-1)/y_nnod +1
        yinod = mod(inod-1, y_nnod)+1
      else
        xinod = mod(inod-1 ,x_nnod)+1
        yinod =    (inod-1)/x_nnod +1
      end if

      return

      end subroutine eigen_convert_ID_w2xy

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_vector_zeropad_x
!
! Purpose
! =======
!
! Fill zeros in padding area of the working buffer
!
!
! Arguments
! =========
!
! v_x    (input/output) real(8) array, dimension(*)
!        working buffer
!
! L      (input) integer
!        Length of global vector
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_vector_zeropad_x(v_x, L)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_vector_zeropad_x

      real(8), intent(inout) :: v_x(*)
      integer, intent(in)    :: L

      integer                :: j_3, n

      j_3 = eigen_loop_end_     (L, x_nnod, x_inod)
      n   = eigen_translate_g2l_(L, x_nnod, x_inod)
      if (j_3 < n) then
        v_x(j_3+1:n) = ZERO     ! in case
      end if

      end subroutine eigen_vector_zeropad_x

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_vector_zeropad_y
!
! Purpose
! =======
!
! Fill zeros in padding area of the working buffer
!
!
! Arguments
! =========
!
! v_y    (input/output) real(8) array, dimension(*)
!        working buffer
!
! L      (input) integer
!        Length of global vector
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_vector_zeropad_y(v_y, L)
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_vector_zeropad_y

      real(8), intent(inout) :: v_y(*)
      integer, intent(in)    :: L
      integer                :: j_3, n

      j_3 = eigen_loop_end_     (L, y_nnod, y_inod)
      n   = eigen_translate_g2l_(L, y_nnod, y_inod)
      if (j_3 < n) then
        v_y(j_3+1:n) = ZERO     ! in case
      end if

      end subroutine eigen_vector_zeropad_y

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_eps
!
! Purpose
! =======
!
! Constant number machine_epsilon
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_eps()
     &     result(r)
!DIR$ ATTRIBUTES FORCEINLINE :: get_constant_eps

      integer(8), target :: const_eps_i8
      real(8), pointer   :: const_eps_r8
      data const_eps_i8 / z'3CB0000000000000' /

      call c_f_pointer( c_loc(const_eps_i8), const_eps_r8 )
      r = const_eps_r8

      return

      end function get_constant_eps
   
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_nan
!
! Purpose
! =======
!
! Cautional number or signaling NaN(not a number)
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_nan()
     &     result(r)
!DIR$ ATTRIBUTES FORCEINLINE :: get_constant_nan

      integer(8), target :: const_nan_i8
      real(8), pointer   :: const_nan_r8
      data const_nan_i8 / z'7FFFFFFFFFFFFFFF' /

      call c_f_pointer( c_loc(const_nan_i8), const_nan_r8 )
      r = const_nan_r8

      return

      end function get_constant_nan
   
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_pai
!
! Purpose
! =======
!
! Constant number PAI
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_pai()
     &     result(r)
!DIR$ ATTRIBUTES FORCEINLINE :: get_constant_pai

      integer(8), target :: const_pai_i8
      real(8), pointer   :: const_pai_r8
      data const_pai_i8 / z'400921FB54442D18' /

      call c_f_pointer( c_loc(const_pai_i8), const_pai_r8 )
      r = const_pai_r8

      return

      end function get_constant_pai

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_2pai
!
! Purpose
! =======
!
! Constant number 2*PAI
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_2pai()
     &     result(r)
!DIR$ ATTRIBUTES FORCEINLINE :: get_constant_2pai

      integer(8), target :: const_2pai_i8
      real(8), pointer   :: const_2pai_r8
      data const_2pai_i8 / z'401921FB54442D18' /

      call c_f_pointer( c_loc(const_2pai_i8), const_2pai_r8 )
      r = const_2pai_r8

      return

      end function  get_constant_2pai
   
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_pai_2
!
! Purpose
! =======
!
! Constant number PAI/2
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_pai_2()
     &     result(r)
!DIR$ ATTRIBUTES FORCEINLINE :: get_constant_pai_2

      integer(8), target :: const_pai2_i8
      real(8), pointer   :: const_pai2_r8
      data const_pai2_i8 / z'3ff921FB54442D18' /

      call c_f_pointer( c_loc(const_pai2_i8), const_pai2_r8 )
      r = const_pai2_r8

      return

      end function get_constant_pai_2

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: 
!
! Purpose
! =======
!
! *
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_diag_loop_info ( lstart, lend,
     &           kstart, kend, istart, istep, jstart, jstep )
!DIR$ ATTRIBUTES FORCEINLINE :: eigen_diag_loop_info

      integer, intent(in)  :: lstart, lend
      integer, intent(out) :: kstart, kend
      integer, intent(out) :: istart, istep
      integer, intent(out) :: jstart, jstep

      integer :: iend, jend


      if ( diag_0 <= 0 .or. lstart > lend ) then
        istart = 0; istep = 1
        jstart = 0; jstep = 1
        kstart = 0; kend = -1
        return
      end if

      istart = eigen_loop_start(lstart, 'Y')
      iend   = eigen_loop_end  (lend,   'Y')
      istart = max( diag_0, istart )

      jstart = eigen_loop_start(lstart, 'X')
      jend   = eigen_loop_end  (lend,   'X')
      jstart = max( diag_1, jstart )

      if ( istart > iend .or. jstart > jend ) then
        istart = 0; istep = 1
        jstart = 0; jstep = 1
        kstart = 0; kend = -1
        return
      end if

      istep  = x_nnod / n_common
      jstep  = y_nnod / n_common

      kstart = 0
      kend   = (iend - istart) / istep

      end subroutine eigen_diag_loop_info

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: 
!
! Purpose
! =======
!
! * NaN (Not a Number) check functions
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_check_vec_NaN ( n, v, msg, info )
      integer       :: n
      real(8)       :: v(1:n)
      character(*)  :: msg
      integer       :: info, info_

      integer :: i_1,i_2,i_3

      info_ = info
      info = 0
      do i_1=1,n
        if ( ieee_is_nan( v(i_1) ) ) then
          if ( info_ /= 0 ) then
          print*,"NaN(",msg,") :",TRD_inod, i_1, v(i_1)
          end if
          info = i_1
        end if
      end do

      end subroutine eigen_check_vec_NaN

      subroutine eigen_check_mat_NaN ( n, nvec, z, ldz, msg, info )
      integer       :: n, nvec, ldz
      real(8)       :: z(ldz,*)
      character(*)  :: msg
      integer       :: info(2), info_

      integer :: i,i_1,i_2,i_3
      integer :: j,j_1,j_2,j_3

      i_2 = eigen_loop_start(1,    'Y')
      i_3 = eigen_loop_end  (nvec, 'Y')
      j_2 = eigen_loop_start(1,    'X')
      j_3 = eigen_loop_end  (n,    'X')

      info_ = info(1)
      info(1:2) = 0
      do i_1=i_2,i_3
        do j_1=j_2,j_3
          if ( ieee_is_nan( Z(j_1,i_1) ) ) then
            i = eigen_translate_l2g(i_1,'Y')
            j = eigen_translate_l2g(j_1,'X')
            if ( info_ /= 0 ) then
            print*,"NaN(",msg,") :",TRD_inod,
     &                  j_1,i_1,"/",j,i,"/",Z(j_1,i_1)
            end if
            info(1) = i_1
            info(2) = j_1
          end if
        end do
      end do

      end subroutine eigen_check_mat_NaN


      end module eigen_libs0_mod

