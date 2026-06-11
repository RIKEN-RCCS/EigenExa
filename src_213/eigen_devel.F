!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_devel_mod
!
! Purpose
! =======
!
! eigen_devel_mod manages the modules for development of the EigenExa
! library, that contains the minimal set to be loaded or linked.
! constant, common variales to coordinate the MPI process groups,
! and timer and profileing routines
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_devel_mod

      use mpi
!$    use omp_lib

      implicit none
      public

!     CONSTANT PARAMETERS
      real(8), parameter     :: ZERO   =  0.0D+00
      real(8), parameter     :: HALF   =  0.5D+00
      real(8), parameter     :: ONE    =  1.0D+00
      real(8), parameter     :: TWO    =  2.0D+00
      real(8), parameter     :: THREE  =  3.0D+00
      real(8), parameter     :: FOUR   =  4.0D+00
      real(8), parameter     :: FIVE   =  5.0D+00
      real(8), parameter     :: SIX    =  6.0D+00
      real(8), parameter     :: SEVEN  =  7.0D+00
      real(8), parameter     :: EIGHT  =  8.0D+00
      real(8), parameter     :: NINE   =  9.0D+00
      real(8), parameter     :: TEN    =  1.0D+01

!     PREFIX MINUS_ MEANS THE NEGATIVE SIGN
      real(8), parameter     :: MHALF  = -0.5D+00
      real(8), parameter     :: MONE   = -1.0D+00
      real(8), parameter     :: MTWO   = -2.0D+00

!     COMPLEX with postfix Z
      complex(8), parameter  :: ZEROZ  = (ZERO, ZERO)
      complex(8), parameter  :: ONEZ   = (ONE,  ZERO)
      complex(8), parameter  :: MONEZ  = (MONE, ZERO)
      complex(8), parameter  :: IONEZ  = (ZERO, ONE)
      complex(8), parameter  :: IMONEZ = (ZERO, MONE)

!$    integer                :: MPI_THREAD_MODE = MPI_THREAD_SINGLE
      character(8), private  :: Process_Grid_Major
      logical,      private  :: Eigen_initialized_flag = .false.

      integer                :: TRD_inod       = 0
      integer                :: TRD_nnod       = 0
      integer                :: TRD_COMM_WORLD = MPI_COMM_WORLD

      integer                :: x_inod, x_nnod, x_COMM_WORLD
      integer                :: y_inod, y_nnod, y_COMM_WORLD
      integer                :: z_inod, z_nnod, z_COMM_WORLD
      integer                :: w_inod, w_nnod, w_COMM_WORLD
      integer                :: n_common, diag_0, diag_1

      integer(8)             :: ERROR_INFO

      real(8)                :: Barrier_Overhead_x
      real(8)                :: Barrier_Overhead_y
      real(8)                :: Reduce_Overhead_x
      real(8)                :: Reduce_Overhead_y
      real(8)                :: Bcast_Overhead_x
      real(8)                :: Bcast_Overhead_y
      real(8)                :: Reduce_cont_Overhead_x
      real(8)                :: Reduce_cont_Overhead_y
      real(8)                :: Bcast_cont_Overhead_x
      real(8)                :: Bcast_cont_Overhead_y

      real(8)                :: comm_time_reduction
      real(8)                :: comm_time_dc
      real(8)                :: comm_time_backtrafo

      real(8)                :: timer_t1, timer_t2
      logical                :: flag_overlap
      logical                :: flag_oncache

!     These must be attributed with pointer for Fortran 200x later
      integer, pointer       :: p0_(:), q0_(:)

!     parameters for backward Householder transform
      integer, parameter     :: nsx = 480
      integer, parameter     :: nsm = 256
      integer, parameter     :: ns0 = nsm * nsm + 6
      integer, parameter     :: MBAND = 2

!$    integer(8) :: TRBK_MASK_FULL = 0

      integer                :: items_bcast, items_reduce
      integer                :: items_redist, items_gather

#if TIMER_PRINT
!     These must be attributed with pointer for Fortran 200x later
      real(8), pointer       :: rt_timer(:,:)
      integer                :: rt_ptr(0:16)

      real(8)                :: time_bcast, time_reduce
      real(8)                :: time_redist, time_gather
      real(8)                :: time_bcast_ (1:16)
      real(8)                :: time_reduce_(1:16)
      real(8)                :: time_redist_(1:16)
      real(8)                :: time_gather_(1:16)
      integer(8)             :: counter_bcast_ (1:16)
      integer(8)             :: counter_reduce_(1:16)
      integer(8)             :: counter_redist_(1:16)
      integer(8)             :: counter_gather_(1:16)
      integer(8)             :: messages_bcast_ (1:16)
      integer(8)             :: messages_reduce_(1:16)
      integer(8)             :: messages_redist_(1:16)
      integer(8)             :: messages_gather_(1:16)
#endif
      logical                :: repro_reduce

!     BLAS functions definitions
      real(8), external      :: ddot
      real(8), external      :: dnrm2
      real(8), external      :: dasum


      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!
! Subroutine: eigen_abort
!
! Purpose
! =======
!
! Abort with a message
!
! Arguments
! =========
!
! message  (input) character(*)
!          Error message
!
! code     (input) integer
!          Error code
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_abort(message, code)

      character(*), intent(in) :: message
      integer,       intent(in) :: code
      integer :: ierr

      call flush(6)
      print*,message
      call flush(6)
      call sleep(1)
      call flush(6)

      call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)

      return

      end subroutine eigen_abort

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_get_wtime
!
! Purpose
! =======
!
! The eigen_get_wtime function returns a double-precision floating point
! value equal to the elapsed wall clock time in seconds.
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function eigen_get_wtime()
     &     result(ret)

      real(8), external :: second


#if __FUJITSU
      ret = MPI_Wtime()
#else
#ifdef _OPENMP
!$    ret = omp_get_wtime()
#else
      ret = MPI_Wtime()
#endif
#endif

      return

      end function eigen_get_wtime

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_set_initialized
!
! Purpose
! =======
!
! eigen_set_initialized subroutine sets the flag TRUE that means 
! initialized.
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_set_initialized()

      Eigen_initialized_flag = .true.

      return

      end subroutine eigen_set_initialized

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_unset_initialized
!
! Purpose
! =======
!
! eigen_unset_initialized subroutine sets the flag FALSE that means 
! not initialized.
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_unset_initialized()

      Eigen_initialized_flag = .false.

      return

      end subroutine eigen_unset_initialized

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_initialized
!
! Purpose
! =======
!
! eigen_get_initialized subroutine gets the flag that means 
! initialized.
!
!
! Arguments
! =========
!
! flag    (output) logical
!         flag that means whether or not initialized.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_initialized(flag)

      logical, intent(out) :: flag

      flag = Eigen_initialized_flag

      return

      end subroutine eigen_get_initialized

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_set_grid_major
!
! Purpose
! =======
!
! eigen_set_grid_major sets up the grid major.
!
!
! Arguments
! =========
!
! Major   (input) character(*)
!         'R': row-major
!         'C': column-major
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_set_grid_major(Major)

      character(*), intent(in) :: Major

      Process_Grid_Major(1:1) = Major(1:1)

      return

      end subroutine eigen_set_grid_major

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_grid_major
!
! Purpose
! =======
!
! eigen_set_grid_major returns the grid major registered.
!
!
! Arguments
! =========
!
! Major   (output) character(*)
!         'R': row-major
!         'C': column-major
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_grid_major(Major)

      character(*), intent(out) :: Major

      Major(1:1) = Process_Grid_Major(1:1)

      return
      
      end subroutine eigen_get_grid_major

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_timer_reset
!
! Purpose
! =======
!
! eigen_timer_reset resets the timer function and the profiling tables.
!
!
! Arguments
! =========
!
! bcast  (input) integer
!        The number of checkpoints for broadcast
!
! reduce (input) integer
!        The number of checkpoints for allreduce
!
! redist (input) integer
!        The number of checkpoints for re-distribution function
!
! gather (input) integer
!        The number of checkpoints for allgather
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_timer_reset(bcast, reduce, redist, gather)

      integer, intent(in), optional  ::  bcast
      integer, intent(in), optional  ::  reduce
      integer, intent(in), optional  ::  redist
      integer, intent(in), optional  ::  gather

      if (present(bcast)) then
         items_bcast  = bcast
      else
         items_bcast   = 0
      end if

      if (present(reduce)) then
         items_reduce = reduce
      else
         items_reduce  = 0
      end if

      if (present(redist)) then
         items_redist = redist
      else
         items_redist  = 0
      end if

      if (present(gather)) then
         items_gather = gather
      else
         items_gather  = 0
      end if

#if TIMER_PRINT

      time_bcast  = ZERO
      time_reduce = ZERO
      time_redist = ZERO
      time_gather = ZERO

      time_bcast_ (1:items_bcast +items_redist) = ZERO
      time_reduce_(1:items_reduce             ) = ZERO
      time_redist_(1:items_redist             ) = ZERO
      time_gather_(1:items_gather+items_redist) = ZERO

      counter_bcast_ (1:items_bcast +items_redist) = 0
      counter_reduce_(1:items_reduce             ) = 0
      counter_redist_(1:items_redist             ) = 0
      counter_gather_(1:items_gather+items_redist) = 0

      messages_bcast_ (1:items_bcast +items_redist) = 0
      messages_reduce_(1:items_reduce             ) = 0
      messages_redist_(1:items_redist             ) = 0
      messages_gather_(1:items_gather+items_redist) = 0
#endif
      
      return

      end subroutine eigen_timer_reset

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_timer_print
!
! Purpose
! =======
!
! print the timer information with message string.
!
!
! Arguments
! =========
!
! message (input) character(*)
!         message string.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function eigen_timer_print(messages)
     &     result(ret)

      character(*), intent(in) :: messages

      real(8)       :: total_time, throughput
      integer       :: i
      integer(8)    :: mes

      total_time = 0.0d0

#if TIMER_PRINT
      if (TRD_inod == 1) then

        print*, "COMM_STAT / [ ",  messages, " ]"

        mes = sum(counter_bcast_(1:items_bcast+items_redist))*8
        if (mes > 0) then
          throughput = mes / time_bcast * 1D-9
        else
          throughput = 0d0
        end if
        print 20000, "   BCAST  :: ", time_bcast, throughput, "[GB/s]"
        do i=1,items_bcast
          print 30000, "          // ",
     &         time_bcast_(i), counter_bcast_(i), messages_bcast_(i)
        end do
        do i=items_bcast+1,items_bcast+items_redist
          print 30000, "          ;; ",
     &         time_bcast_(i), counter_bcast_(i), messages_bcast_(i)
        end do

        mes = sum(counter_reduce_(1:items_reduce))*8
        if (mes > 0) then
          throughput = mes / time_reduce * 1D-9
        else
          throughput = 0d0
        end if
        print 20000, "   REDUCE :: ", time_reduce, throughput, "[GB/s]"
        do i=1,items_reduce
          print 30000, "          // ",
     &         time_reduce_(i), counter_reduce_(i), messages_reduce_(i)
        end do

        mes = sum(counter_gather_(1:items_gather+items_redist))*8
        if (mes > 0) then
          throughput = mes / time_gather * 1D-9
        else
          throughput = 0d0
        end if
        print 20000, "   GATHER :: ", time_gather, throughput, "[GB/s]"
        do i=1,items_gather
          print 30000, "          // ",
     &         time_gather_(i), counter_gather_(i), messages_gather_(i)
        end do
        do i=items_gather+1,items_gather+items_redist
          print 30000, "          ;; ",
     &         time_gather_(i), counter_gather_(i), messages_gather_(i)
        end do

        print 10000, "   REDIST :: ", time_redist
        do i=1,items_redist
          print 30000, "          // ",
     &         time_redist_(i), counter_redist_(i), messages_redist_(i)
        end do

        total_time = time_bcast + time_reduce
     &       + time_redist
     &       + time_gather
        print 10000, "   Total  :: ", total_time

10000   format(X,A,E25.16e2)
20000   format(X,A,2E25.16e2,X,A)
30000   format(X,A,E25.16e2,I14,I10)

      end if
      items_bcast   = 0
      items_reduce  = 0
      items_redist  = 0
      items_gather  = 0

#endif
      ret = total_time

      return

      end function eigen_timer_print

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: sync_other_than_master_init
!
! Purpose
! =======
!
! Initialize the syncronization among threads other than the master.
!
!
! Arguments
! =========
!
! TRBK_lock  (input/output) integer(kind=omp_lock_kind)
!            mutex variable
!
! TRBK_mask  (input/output) integer(8)
!            internal working variable
!
!--------*---------*---------*---------*---------*---------*---------*-*
!$
!$    subroutine sync_other_than_master_init(TRBK_lock, TRBK_mask)
!$
!$    integer(kind=omp_lock_kind) :: TRBK_lock
!$    integer(8)                  :: TRBK_mask(1:2)
!$
!$    integer :: i
!$
!$    TRBK_MASK_FULL = 0
!$    if (omp_get_num_threads() > 64) then
!$      TRBK_MASK_FULL = omp_get_num_threads() - 1
!$    else
!$      do i=1,omp_get_num_threads()-1
!$        TRBK_MASK_FULL = ibset(TRBK_MASK_FULL, i)
!$      end do
!$    end if
!$
!$    call omp_init_lock(TRBK_lock)
!$
!$    call omp_set_lock(TRBK_lock)
!$    TRBK_mask(1) = 0
!$    TRBK_mask(2) = 0
!$    call omp_unset_lock(TRBK_lock)
!$
!$    end subroutine sync_other_than_master_init

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: sync_other_than_master_finalize
!
! Purpose
! =======
!
! Finalize the syncronization among threads other than the master.
!
!
! Arguments
! =========
!
! TRBK_lock  (input/output) integer(kind=omp_lock_kind)
!            mutex variable
!
!--------*---------*---------*---------*---------*---------*---------*-*
!$    subroutine sync_other_than_master_finalize(TRBK_lock)
!$
!$    integer(kind=omp_lock_kind) :: TRBK_lock
!$
!$    call omp_destroy_lock(TRBK_lock)
!$
!$    end subroutine sync_other_than_master_finalize
!$

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: sync_other_than_master
!
! Purpose
! =======
!
! Syncronization among threads other than the master.
!
!
! Arguments
! =========
!
! TRBK_lock  (input/output) integer(kind=omp_lock_kind)
!            mutex variable
!
! TRBK_mask  (input/output) integer(8)
!            internal working variable
!
!--------*---------*---------*---------*---------*---------*---------*-*
!$
!$    subroutine sync_other_than_master(TRBK_lock, TRBK_mask)
!$
!$    integer(kind=omp_lock_kind) :: TRBK_lock
!$    integer(8)                  :: TRBK_mask(1:2)
!$
!$    integer(8) :: T
!$    integer    :: local_size, local_rank
!$
!$    local_size = omp_get_num_threads()
!$    local_rank = omp_get_thread_num()
!$
!$    if (local_size == 1) return
!$    if (local_rank == 0) return
!$
!$    call omp_set_lock(TRBK_lock)
!$    if (local_size > 64) then
!$      T = TRBK_mask(1) + 1
!$    else
!$      T = ibset(TRBK_mask(1), local_rank)
!$    end if
!$    if (T == TRBK_MASK_FULL) TRBK_mask(2) = T
!$    TRBK_mask(1) = T
!$    call omp_unset_lock(TRBK_lock)
!$
!$    do
!$      call omp_set_lock(TRBK_lock)
!$      T = TRBK_mask(1)
!$      call omp_unset_lock(TRBK_lock)
!$      if (T == TRBK_MASK_FULL) exit
!$    end do
!$
!$    call omp_set_lock(TRBK_lock)
!$    if (local_size > 64) then
!$      T = TRBK_mask(2) - 1
!$    else
!$      T = ibclr(TRBK_mask(2), local_rank)
!$    end if
!$    if (T == 0) TRBK_mask(1) = T
!$    TRBK_mask(2) = T
!$    call omp_unset_lock(TRBK_lock)
!$
!$    do
!$      call omp_set_lock(TRBK_lock)
!$      T = TRBK_mask(2)
!$      call omp_unset_lock(TRBK_lock)
!$      if (T == 0) exit
!$    end do
!$
!$    end subroutine sync_other_than_master

      end module eigen_devel_mod
