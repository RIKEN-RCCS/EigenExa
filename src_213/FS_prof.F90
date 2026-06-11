!#define COUNT_CHECK
!>
!> @file   FS_prof.F90
!> @brief  module FS_prof_mod
!!>
!
!
!
!> module FS_prof_mod
!> @brief  @n
!> Purpose @n
!> ======= @n
!> FS_prof_mod include measurement of exection time control
!>
module FS_prof_mod
!$    use omp_lib
      implicit none

      !> max number of region for prof data
      integer,parameter   :: FS_max_region = 70

      !> type FS_prof
      !> @brief prof data struct
      type FS_prof
        character*32        :: region_name(FS_max_region)   !< region name
        real(8)              :: region_time(FS_max_region)   !< region time
        real(8)              :: region_start(FS_max_region)  !< region start time (temporary)
        integer             :: region_ecount(FS_max_region) !< region count of called FS_prof_end
#ifdef COUNT_CHECK
        integer             :: region_scount(FS_max_region) !< region count of called FS_prof_start
#endif
      end type FS_prof

contains

!--------*---------*---------*---------*---------*---------*---------*-*
! prof initialize
!--------*---------*---------*---------*---------*---------*---------*-*
      !> subroutine FS_prof_init
      !> @brief initialize prof information
      !> @param[in,out] prof  prof information
      subroutine FS_prof_init(prof)
      implicit none

      type(FS_prof), intent(inout) :: prof

                              !12345678901234567890123456789012
      prof%region_name(  1) = "total"
      prof%region_name( 10) = "FS_EDC"
      prof%region_name( 11) = "  DSTEDC(NPROCS=1)"
      prof%region_name( 20) = "  FS_PDLAED0"
      prof%region_name( 21) = "    FS_dividing"
      prof%region_name( 22) = "      FS_create_merge_comm"
      prof%region_name( 23) = "      barrier"
      prof%region_name( 28) = "    DSTEDC"
      prof%region_name( 30) = "    FS_PDLAED1"
      prof%region_name( 31) = "      FS_MERGE_D"
      prof%region_name( 40) = "      FS_PDLAEDZ"
      prof%region_name( 45) = "      FS_REDUCE_ZD"
      prof%region_name( 46) = "        barrier"
      prof%region_name( 47) = "        allreduce"
      prof%region_name( 50) = "      FS_PDLAED2"
      prof%region_name( 60) = "      FS_PDLAED3"
      prof%region_name( 61) = "        DLAED4(Z,D)"
      prof%region_name( 62) = "        allreduce x3"
      prof%region_name( 63) = "        COPY Q2"
      prof%region_name( 64) = "        WAIT/UNPACK"
      prof%region_name( 65) = "        PACK/ISEND/IRECV"
!     prof%region_name( 66) = "        DLAED4(DELTA)+DGEMM"
      prof%region_name( 66) = "        RECALC DELTA+DGEMM"
      prof%region_name( 67) = "          DGEMM"
      prof%region_name( 70) = "    FS_PDLASRT"

      prof%region_time  = 0.0d0
      prof%region_start = 0.0d0
      prof%region_ecount = 0
#ifdef COUNT_CHECK
      prof%region_scount = 0
#endif

      return
      end subroutine FS_prof_init

!--------*---------*---------*---------*---------*---------*---------*-*
! start prof region
!--------*---------*---------*---------*---------*---------*---------*-*
      !> subroutine FS_prof_start
      !> @brief start measurment of region [ID]
      !> @param[in,out] prof  prof information
      !> @param[in]     ID    region ID
      subroutine FS_prof_start(prof, ID)
        use mpi
        implicit none

        type(FS_prof), intent(inout) :: prof
        integer, intent(in) :: ID
        real(8) :: tt

#ifdef _OPENMP
!$      tt = omp_get_wtime()
#else
        tt = MPI_Wtime()
#endif
        prof%region_start(ID) = tt

#ifdef COUNT_CHECK
        prof%region_scount(ID) = prof%region_scount(ID) + 1
#endif

        return
      end subroutine FS_prof_start

!--------*---------*---------*---------*---------*---------*---------*-*
! end prof region
!--------*---------*---------*---------*---------*---------*---------*-*
      !> subroutine FS_prof_end
      !> @brief stop measurment of region [ID]
      !> @param[in,out] prof  prof information
      !> @param[in]     ID    region ID
      subroutine FS_prof_end(prof, ID)
        use mpi
        implicit none

        type(FS_prof), intent(inout) :: prof
        integer, intent(in) :: ID
        real(8) :: tt

#ifdef _OPENMP
!$      tt = omp_get_wtime()
#else
        tt = MPI_Wtime()
#endif
        prof%region_time(ID) = prof%region_time(ID) + (tt - prof%region_start(ID))

        prof%region_ecount(ID) = prof%region_ecount(ID) + 1

        return
      end subroutine FS_prof_end

!--------*---------*---------*---------*---------*---------*---------*-*
! prof finalize
!--------*---------*---------*---------*---------*---------*---------*-*
      !> subroutine FS_prof_finalize
      !> @brief finalize prof mesurement and print log
      !> @param[in,out] prof    prof information
      !> @param[in]     kout    file No. of log file (default:standard output)
      !> @param[in]     outall  output all procs flag (default:only rank0)
      subroutine FS_prof_finalize(prof, kout, outall)
        use mpi
        use FS_libs_mod
        implicit none

        type(FS_prof), intent(inout) :: prof
        integer, optional, intent(in) :: kout
        integer, optional, intent(in) :: outall
        integer :: kout_, i
        integer :: ierr, n, nprocs_out
        integer :: stat(MPI_STATUS_SIZE)

        real(8)   :: region_time(FS_max_region)
#ifdef COUNT_CHECK
        integer  :: region_scount(FS_max_region)
#endif
        integer  :: region_ecount(FS_max_region)

        ! set file No.
        if( present(kout) ) then
          kout_ = kout
        else
          kout_ = 6
        endif

        ! out procs
        call MPI_COMM_SIZE(FS_COMM_WORLD, nprocs_out, ierr)
        if( present(outall) ) then
          if( outall.eq.0 ) then
            nprocs_out = 1
          endif
        endif
#ifndef TIMER_ALLPROCS
        nprocs_out = 1
#endif

        ! 集計と出力
        if( FS_MYRANK.eq.0 ) then

          ! ランク0が出力
          write(kout_,'(a)') " ===================================================================="
#ifdef _BLOCKING_DGEMM
          write(kout_,'(a)') "  TIMING INFO (DGEMM BLOCKING)"
#else
          write(kout_,'(a)') "  TIMING INFO (DGEMM NON-BLOCKING)"
#endif
          write(kout_,'(a)') " ===================================================================="

          ! 各ランクから情報を受け取りながら出力
          do n=0, nprocs_out-1

            ! ランク0以外から情報を受信
            if( n.eq.0 ) then
              region_time   = prof%region_time
#ifdef COUNT_CHECK
              region_scount = prof%region_scount
#endif
              region_ecount = prof%region_ecount
            else
              call MPI_RECV(region_time  ,FS_max_region,MPI_DOUBLE_PRECISION  ,n,0,FS_COMM_WORLD,stat,ierr)
#ifdef COUNT_CHECK
              call MPI_RECV(region_scount,FS_max_region,MPI_INTEGER,n,0,FS_COMM_WORLD,stat,ierr)
#endif
              call MPI_RECV(region_ecount,FS_max_region,MPI_INTEGER,n,0,FS_COMM_WORLD,stat,ierr)
            endif

            write(kout_,'(a,i0)') "  RANK = ", n
            write(kout_,'(a)') " -ID-+----region name ----------------+----time [s] ---------+-count-"
            do i=1,FS_max_region
              if( region_ecount(i).gt.0 ) then
                write(kout_,'(1x,i5,1x,a32,1x,1pe22.15,1x,i6)') &
                  i, prof%region_name(i), region_time(i), region_ecount(i)
#ifdef COUNT_CHECK
                if( region_scount(i) .ne. region_ecount(i) ) then
                  write(kout_,'(1x,3a)') "  Warning : start/end count are different in [", trim(prof%region_name(i)), "]"
                endif
#endif
              endif
            enddo
            write(kout_,'(a)') " ===================================================================="
          enddo

        else ! FS_MYRANK.ne.0

          ! ランク0に情報を送信
#ifdef TIMER_ALLPROCS
          if( nprocs_out.gt.1 ) then
            call MPI_SEND(prof%region_time  ,FS_max_region,MPI_DOUBLE_PRECISION  ,0,0,FS_COMM_WORLD,ierr)
#ifdef COUNT_CHECK
            call MPI_SEND(prof%region_scount,FS_max_region,MPI_INTEGER,0,0,FS_COMM_WORLD,ierr)
#endif
            call MPI_SEND(prof%region_ecount,FS_max_region,MPI_INTEGER,0,0,FS_COMM_WORLD,ierr)
          endif
#endif

        endif

        return
      end subroutine FS_prof_finalize

!--------*---------*---------*---------*---------*---------*---------*-*
! prof add
!--------*---------*---------*---------*---------*---------*---------*-*
      subroutine FS_prof_add(prof, prof_add)
        use mpi
        implicit none

        type(FS_prof), intent(inout) :: prof
        type(FS_prof), intent(in) :: prof_add
        integer :: i

        do i=1, FS_max_region
          prof%region_time(i)   = prof%region_time(i)   + prof_add%region_time(i)
          prof%region_ecount(i) = prof%region_ecount(i) + prof_add%region_ecount(i)
#ifdef COUNT_CHECK
          prof%region_scount(i) = prof%region_scount(i) + prof_add%region_scount(i)
#endif
        end do

        return
      end subroutine FS_prof_add

end module FS_prof_mod
