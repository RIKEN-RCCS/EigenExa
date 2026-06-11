!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: trbakwy4_body_mod
!
! Purpose
! =======
!
! trbakwy4_body_mod manages the modules of the Householder back
! transformation routines.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module trbakwy4_body_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use eigen_trbak_mod
      use mpi
!$    use omp_lib

      implicit none
      private


      public  :: eigen_trbakwy_block_body
#if AT_BCAST_OVERLAP
!$    private :: trbk_decide_overlap_level
#endif
      private :: eigen_trbakwy_block_body1
      private :: eigen_trbakwy_block_body2

      contains

#if __FUJITSU
#define dgemm   dl_gemm
#define dtrsm   dl_trsm
#endif

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trbakwy_block_body
!
! Purpose
! =======
!
! eigen_trbakwy_block_body is the core routine of the Householder back
! transformation algorithm with the WY representation.
!
!
! Arguments
! =========
!
! local_nvec (input) integer
! < description >
!
! z       (input/output) real(8) array, dimension(nmz,*)
!         z contains the elements of the eigenvectors.
!
! nmz     (input) integer
!         The leading dimension of the array a.
!
! d       (input/output) real(8) array, dimension(*)
!         d contains auxiliary data.
!
! e       (input/output) real(8) array, dimension(*)
!         d contains auxiliary data.
!
! v       (input/output) real(8) array, dimension(nm,*)
!         v contains the block reflector.
!
! nm      (input) integer
!         The leading dimension of the array v.
!
! m       (input) integer
!         Block factor
!
! i       (input) integer
!         Global loop conter
!
! ss      (input/output) real(8) array, dimension(*)
!         ss contains the results of v^t*(Z,v).
!
! tt      (input/output) real(8) array, dimension(*)
!         working buffer
!
! nss     (input) integer
!         The leading dimension of the internal array ss.
!
! iblk    (input) integer
!         Auxiliary data for back-transformation (tri-/penta-/...).
!
! dcom    (input/output) real(8)
!         dcom contains the accumulation of the communication time.
!
! dx     (input/output) real(8)
!         dx contains the accumulation of the computaiontal time.
!
! dy     (input/output) real(8)
!         dy contains the accumulation of the computaiontal time.
!
! dz     (input/output) real(8)
!         dz contains the accumulation of the computaiontal time.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trbakwy_block_body(local_nvec, z, nmz,
     &     v, nm, m, i, ss, tt, nss, iblk,
     &     dcom, dx, dy, dz)
     
      integer, intent(in)    :: local_nvec
      real(8), intent(inout) :: z(nmz,*)
      integer, intent(in)    :: nmz
      real(8), intent(inout) :: v(nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: m
      integer, intent(in)    :: i
      real(8), intent(inout) :: ss(*)
      real(8), intent(inout) :: tt(*)
      integer, intent(in)    :: nss
      integer, intent(in)    :: iblk
      real(8), intent(inout) :: dcom
      real(8), intent(inout) :: dx
      real(8), intent(inout) :: dy
      real(8), intent(inout) :: dz

      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_2, j_3
      integer                :: i_0, m_0
      integer                :: ii, jj
      real(8)                :: fr
      real(8)                :: ds, de

      integer                :: local_rank, local_size

      local_size = 1
      local_rank = 0
!$    local_size = omp_get_num_threads()
!$    local_rank = omp_get_thread_num()

      i_2 = 1
      i_3 = local_nvec

      j_2 = eigen_loop_start(1,         'X')
      j_3 = eigen_loop_end  (i+m-1-iblk,'X')


#if AT_BCAST_OVERLAP
!$    if (DO_OVERLAP_BCAST_LEVEL == 2) then
!===  
!$      call sync_other_than_master(TRBK_LOCK, TRBK_MASK)
!===  
!$    end if
#endif
!$    if (
#if AT_BCAST_OVERLAP
!$   &     DO_OVERLAP_BCAST_LEVEL /= 2
#else
!$   &     .true.
#endif
!$   &     ) then
!===  
!$OMP     BARRIER
!===  
!$    end if
!-
      jj = i_3-i_2+1
#if AT_REDUCE_OVERLAP
      fr = TRBK_TIME_FR

      if (TRD_nnod > 1 .and.
     &     TRBK_TIME_REDUC_OVERHEAD_x
     &     < TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x) then

        ii = max(8,int(jj*fr))
        ii = ((ii-1)/13+1)*13
        ii = min(jj,ii)
        fr = (TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x)
     &       / (TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x + TRBK_TIME_REDUC)

      else

        fr = ZERO
        ii = jj

      end if

#if TIMER_PRINT>1
!$OMP  MASTER
      if (TRD_inod == 1) then
        write(6, '(A,X,I,X,I,X,4F)' )
     &       "Overlap[on,off] ", ii, (i_3-i_2+1)-ii,
     &       TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x,
     &       TRBK_TIME_REDUC,
     &       (TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x)*(jj-ii)/jj,
     &       TRBK_TIME_REDUC*ii/jj
      end if
!$OMP  END MASTER
#endif

#else
      fr = ZERO
      ii = i_3-i_2+1
#endif
!-
!     SM:= SM+lower(V^TV)
!     SS(1:m,1:n):= SS+V(j_2:j_3,1:m)^T*Z(j_2:j_3,1:n)
      if (local_size==1 .or. local_rank==1) then
         ds = eigen_get_wtime()
      end if
      call eigen_trbakwy_block_body1(z, nmz,
     &     v, nm, m, ss(1+ns0), ss(1),
     &     i_2,i_2+ii-1,j_2,j_3)
      if (local_size==1 .or. local_rank==1) then
        de = eigen_get_wtime()
!$      TRBK_TIME_TRBK1_ = de - ds
#if TIMER_PRINT
        dx = dx + (de - ds)
#endif
      end if

!-
!-----------------------------------------------------------------
!$OMP  BARRIER
!-----------------------------------------------------------------
!-
!$OMP  MASTER
      ds = eigen_get_wtime()
      ss(ns0-5) = TRBK_TIME_REDUC
      ss(ns0-4) = TRBK_TIME_BCAST
      ss(ns0-3) = TRBK_TIME_TRBK1_
      ss(ns0-2) = TRBK_TIME_TRBK1x_
      ss(ns0-1) = TRBK_TIME_TRBK2
      ss(ns0-0) = fr
      call reduce_dbl(ss(1), tt(1), ns0 + ii*m,
     &     3, x_COMM_WORLD)
!     call reduce_dbl(ss(1), tt(1), ns0,
!     &     3, x_COMM_WORLD)
      TRBK_TIME_REDUC  = ss(ns0-5) / x_nnod
      TRBK_TIME_BCAST  = ss(ns0-4) / x_nnod
      TRBK_TIME_TRBK1  = ss(ns0-3) / x_nnod
      TRBK_TIME_TRBK1x = ss(ns0-2) / x_nnod
      TRBK_TIME_TRBK2  = ss(ns0-1) / x_nnod
      TRBK_TIME_FR     = ss(ns0-0) / x_nnod
      de = eigen_get_wtime()
!$    TRBK_TIME_REDUC = (de-ds)
      dcom = dcom + (de - ds)
!$OMP  END MASTER
!-
#if AT_REDUCE_OVERLAP
      if (i_3 >= i_2 + ii) then
#if TIMER_PRINT
        if (local_size==1 .or. local_rank==1) then
          ds = eigen_get_wtime()
        end if
#endif
        call eigen_trbakwy_block_body1(z, nmz,
     &       v, nm, m, ss(1+ns0), ss(1),
     &       i_2+ii,i_3,j_2,j_3)
#if TIMER_PRINT
        if (local_size==1 .or. local_rank==1) then
          de = eigen_get_wtime()
          dx = dx + (de - ds)
        end if
#endif
      end if
#endif

!-----------------------------------------------------------------
!$OMP  BARRIER
!-----------------------------------------------------------------

#if AT_REDUCE_OVERLAP
      if (local_size==1 .or. local_rank==1) then
        if (i_3 >= i_2 + ii) then
!$        TRBK_TIME_TRBK1x_ = de - ds
        else
!$        TRBK_TIME_TRBK1x_ = ZERO
        end if
      end if
#endif

#if AT_REDUCE_OVERLAP
!$OMP  MASTER
!     if (i_3 >= i_2 + ii) then
      ds = eigen_get_wtime()
      call reduce_dbl(ss(1+ns0+ii*m), tt(1+ns0+ii*m),
     &     (i_3-i_2+1-ii)*m,
     &     4, x_COMM_WORLD)
      de = eigen_get_wtime()
!$    TRBK_TIME_REDUC = TRBK_TIME_REDUC + (de-ds)
      dcom = dcom + (de - ds)
!     end if
!$OMP  END MASTER
#endif
!$OMP  BARRIER

#if TIMER_PRINT
      ds = eigen_get_wtime()
#endif

!====================================
!     diag(SM):= diag(SM ? SM/2 : 1)
!$OMP  DO
      do m_0=1,m
         i_1 = m_0 + (m_0-1)*nsm
         if ( ss(i_1) == 0D0 ) then
           ss(i_1) = 1D0
         else
           ss(i_1) = ss(i_1) * 5D-1
         end if
      end do
!$OMP  ENDDO
!====================================

!$OMP  BARRIER
!$OMP  MASTER

#if AT_BCAST_OVERLAP
!$    if (TRBK_SWITCHED < 4) then
!$      call trbk_decide_overlap_level(i)
!$    end if
#endif

#if TIMER_PRINT
      de = eigen_get_wtime()
      dz = dz + (de - ds)
#endif
!$OMP  END MASTER

!-----------------------------------------------------------------
!$OMP  BARRIER
!-----------------------------------------------------------------

!     V:= V*(SM)^{-1}
!     Z(j_2:j_3,1:n):= Z + V(j_2:j_3,1:m)*SS(1:m,1:n)
      if (local_size==1 .or. local_rank==1) then
        ds = eigen_get_wtime()
      end if
      call eigen_trbakwy_block_body2(z,
     &     nmz, v, nm, m,
     &     ss(1+ns0), ss(1),
     &     i_2,i_3,j_2,j_3)

      if (local_size==1 .or. local_rank==1) then
        de = eigen_get_wtime()
!$      TRBK_TIME_TRBK2 = de - ds
#if TIMER_PRINT
        dy = dy + (de - ds)
#endif
      end if

#if AT_BCAST_OVERLAP
!$    if (TRBK_SWITCHED < 4) then
!$       if (DO_OVERLAP_BCAST_LEVEL == 0) then
!$OMP     BARRIER
!$       end if
!$    end if
#endif

      return

      end subroutine eigen_trbakwy_block_body

!------------------------------------------------------------------
!------------------------------------------------------------------

!===  
#if AT_BCAST_OVERLAP
!===  

!$    subroutine trbk_decide_overlap_level(i)
!$
!$    integer, intent(in)    :: i
!$
!$    integer                :: ll0, ll1, ll2, ll_(6)
!$    real(8)                :: BCAST_TIME(0:2)
!$    real(8)                :: f0, f1
!$    integer                :: local_rank, local_size
!$    integer                :: ierr
!-
!$    local_size = omp_get_num_threads()
!$    local_rank = omp_get_thread_num()
!$
!$    if (local_size > 1) then
!$
!$      TRBK_TIME_COUNTER = TRBK_TIME_COUNTER + 1
!$      if (TRBK_TIME_COUNTER >= TRBK_TIME_NEXT) then
!$        f0 = 1d0 * (local_size-1) / local_size
!$        f1 = 1d0 / f0
!$        if (DO_OVERLAP_BCAST_LEVEL == 2) then
!$          BCAST_TIME(0) = TRBK_TIME_BCAST
!$   &           + TRBK_TIME_TRBK2 * f0
!$   &           + TRBK_TIME_TRBK1 * f0
!$          BCAST_TIME(1) = max(TRBK_TIME_BCAST,
!$   &           TRBK_TIME_TRBK2)
!$   &           + TRBK_TIME_TRBK1 * f0
!$          BCAST_TIME(2) = max(TRBK_TIME_BCAST,
!$   &           TRBK_TIME_TRBK2
!$   &           + TRBK_TIME_TRBK1)
!$        end if
!$        if (DO_OVERLAP_BCAST_LEVEL == 1) then
!$          BCAST_TIME(0) = TRBK_TIME_BCAST
!$   &           + TRBK_TIME_TRBK2 * f0
!$   &           + TRBK_TIME_TRBK1
!$          BCAST_TIME(1) = max(TRBK_TIME_BCAST,
!$   &           TRBK_TIME_TRBK2)
!$   &           + TRBK_TIME_TRBK1
!$          BCAST_TIME(2) = max(TRBK_TIME_BCAST,
!$   &           TRBK_TIME_TRBK2
!$   &           + TRBK_TIME_TRBK1 * f1)
!$        end if
!$        if (DO_OVERLAP_BCAST_LEVEL == 0) then
!$          BCAST_TIME(0) = TRBK_TIME_BCAST
!$   &           + TRBK_TIME_TRBK2
!$   &           + TRBK_TIME_TRBK1
!$          BCAST_TIME(1) = max(TRBK_TIME_BCAST,
!$   &           TRBK_TIME_TRBK2 * f1)
!$   &           + TRBK_TIME_TRBK1
!$          BCAST_TIME(2) = max(TRBK_TIME_BCAST,
!$   &           TRBK_TIME_TRBK2 * f1
!$   &           + TRBK_TIME_TRBK1 * f1)
!$        end if
!$        ll0 = (DO_OVERLAP_BCAST_LEVEL)
!$        ll1 = 0
!$        do ll2 = 1, 2
!$          if (BCAST_TIME(ll1) >= BCAST_TIME(ll2)) then
!$            ll1 = ll2
!$          end if
!$        end do
!$        if (y_nnod >= 1024) then
!$          ll_(1:3)=0; ll_(ll1+1)=1
!$          call MPI_Allreduce(ll_(1), ll_(4),
!$   &           3, MPI_INTEGER, MPI_SUM, y_COMM_WORLD, ierr)
!$        else
!$          ll1 = ll_(1) + (ll_(2) + ll_(3) * 1024) * 1024
!$          ll_(1) = ll1
!$          call MPI_Allreduce(ll_(1), ll_(4),
!$   &           1, MPI_INTEGER, MPI_SUM, y_COMM_WORLD, ierr)
!$          ll1 = ll_(4)
!$          ll_(4) = mod(ll1,1024)
!$          ll1 = ll1 / 1024
!$          ll_(5) = mod(ll1,1024)
!$          ll_(6) = ll1 / 1024
!$        end if
#if OVERLAP_DECISION_TYPE==3
!     type 3 mean value
!$        ll1 = 0
!$        do ll2 = 1, 2
!$          if (ll_(4+ll2) > ll_(4+ll1)) then
!$            ll1 = ll2
!$          end if
!$        end do
#endif
#if OVERLAP_DECISION_TYPE==2
!     type 2 max value
!$        do ll2 = 0, 2
!$          if (ll_(4+ll2) > 0) then
!$            ll1 = ll2
!$          end if
!$        end do
#endif
#if OVERLAP_DECISION_TYPE==1
!     type 1 min value
!$        do ll2 = 2, 0, -1
!$          if (ll_(4+ll2) > 0) then
!$            ll1 = ll2
!$          end if
!$        end do
#endif
#if _DEBUG_
!$        if (TRD_inod == 1) then
!$          print*, i, "Overlap decision", ll0,"->",ll1,
!$   &           BCAST_TIME(0:2)
!$        end if
#endif
!$        if (ll1 == 0) then
!$          TRBK_SWITCHED = TRBK_SWITCHED + 1
!$        end if
!$        if (ll0 /= ll1) then
!$          if (TRBK_SWITCHED >= 8) then
!$            TRBK_TIME_INTERVAL = 32
!$          else
!$            TRBK_TIME_INTERVAL = 1
!$          end if
!$        else
!$          TRBK_TIME_INTERVAL = TRBK_TIME_INTERVAL*2
!$        end if
!$        TRBK_TIME_NEXT = TRBK_TIME_NEXT + TRBK_TIME_INTERVAL
!$        if (.true.) DO_OVERLAP_BCAST_LEVEL = ll1
!$      end if
!$    else
!$      if (.true.) DO_OVERLAP_BCAST_LEVEL = 0
!$    end if
!-
!$    end subroutine trbk_decide_overlap_level
!===  
#endif
!===

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine eigen_trbakwy_block_body1(z, nmz,
     &     v, nm, m, ss, sm, i_2,i_3,j_2,j_3)

      real(8), intent(inout) :: z(nmz,*)
      integer, intent(in)    :: nmz
      real(8), intent(inout) :: v(nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: m
      real(8), intent(inout) :: ss(m,*)
      real(8), intent(inout) :: sm(nsm,nsm)
      integer, intent(in)    :: i_2
      integer, intent(in)    :: i_3
      integer, intent(in)    :: j_2
      integer, intent(in)    :: j_3

      integer, parameter     :: BLAS_CHUNK1 = 16
      integer, parameter     :: BLAS_CHUNK2 = 16
      integer, parameter     :: BLAS_CHUNK3 = 64

      integer                :: i_1, i_4, i_5
      integer                :: j_1, j_5
      integer                :: ii_2, ii_3, jj_3

      integer                :: local_rank, local_size
      integer                :: ll_rank, ll_size

      integer                :: ii_step, blk_size1, blk_size2

     
      local_size = 1
      local_rank = 0
!$    local_size = omp_get_num_threads()
!$    local_rank = omp_get_thread_num()
!-
!$    if ((DO_OVERLAP_BCAST_LEVEL <= 1 .and. i_2 == 1) .or.
!$   &     local_size == 1 .or.
!$   &     local_rank >= 1) then
!$      if (
!$   &        (DO_OVERLAP_BCAST_LEVEL == 2 .and. local_size > 1) .or.
!$   &        (i_2 > 1 .and. local_size > 1)) then
!$        ll_size  = local_size - 1
!$        ll_rank  = local_rank - 1
!$      else
          ll_size  = local_size
          ll_rank  = local_rank
!$      end if

        if (i_2 == 1) then

          j_5 = j_3 - j_2 + 1
          if (j_5 > 0) then

            ii_step = 0
            do i_1 = 1, m, BLAS_CHUNK1
              do j_1 = 1, m, BLAS_CHUNK2

                ii_3 = min(m, i_1+BLAS_CHUNK1-1)
                jj_3 = min(m, j_1+BLAS_CHUNK2-1)

                blk_size1 = ii_3 - i_1 + 1
                blk_size2 = jj_3 - j_1 + 1

                if (blk_size1 > 0 .and. blk_size2 > 0) then
                  if (i_1+blk_size1 >= j_1+blk_size2) then

                    ii_2 = i_1

                    if (mod(ii_step, ll_size) == ll_rank) then

                      call dgemm('T','N',
     &                     blk_size1, blk_size2, j_5,
     &                     MONE, v(j_2, ii_2), nm,
     &                     v(j_2, j_1),  nm,
     &                     ZERO, sm(ii_2,j_1), nsm)

                    end if

                    ii_step = ii_step + 1

                  end if
                end if

              end do
            end do

          end if

        end if

#if 1 || defined(__INTEL_COMPILER)
        j_5 = j_3 - j_2 + 1
        ii_step = 0

        do ii_2 = i_2, i_3, BLAS_CHUNK3

          i_5 = min(BLAS_CHUNK3, i_3-ii_2+1)

          if (m > 0 .and. j_5 > 0 .and. i_5 > 0) then
            if (mod(ii_step, ll_size) == ll_rank) then

              call dgemm('T','N',
     &             m, i_5, j_5,
     &             ONE,  v (j_2 ,1       ), nm,
     &             z (j_2 ,ii_2), nmz,
     &             ZERO,  ss(1   ,ii_2), m)

            end if
            ii_step = ii_step + 1
          end if

        end do
#else
        j_5 = j_3 - j_2 + 1

        i_5 = (i_3-i_2) / ll_size + 1
        i_5 = i_5 + mod(i_5,2)
        i_4 = i_5 * ll_rank
        i_5 = min(i_5, i_3-(i_2+i_4)+1)

        if (m > 0 .and. j_5 > 0 .and. i_5 > 0) then

          call dgemm('T','N',
     &         m, i_5, j_5,
     &         ONE,  v (j_2 ,1       ), nm,
     &         z (j_2 ,i_2+i_4), nmz,
     &         ONE,  ss(1   ,i_2+i_4), m)

        end if
#endif
!$    end if

      return

      end subroutine eigen_trbakwy_block_body1

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine eigen_trbakwy_block_body2(z,
     &     nmz,
     &     v, nm, m, ss, sm, i_2,i_3, j_2,j_3 )

      real(8), intent(inout) :: z(nmz,*)
      integer, intent(in)    :: nmz
      real(8), intent(inout) :: v(nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: m
      real(8), intent(inout) :: ss(m,*)
      real(8), intent(inout) :: sm(nsm,*)
      integer, intent(in)    :: i_2
      integer, intent(in)    :: i_3
      integer, intent(in)    :: j_2
      integer, intent(in)    :: j_3

#if __FUJITSU
! should be multiple of 224
      integer, parameter     :: CHUNK_M = 96
! should be multiple of 40
      integer, parameter     :: CHUNK_N = 40
#else
      integer, parameter     :: CHUNK_M = 128
      integer, parameter     :: CHUNK_N = 48
#endif

      integer                :: local_rank, local_size
      integer                :: ll_rank, ll_size
      integer                :: ii_step
      integer                :: i_5, j_4, j_5
      integer                :: ii_0, ii_2, ii_3
      integer                :: jj_0, jj_2, jj_3

      local_size = 1
      local_rank = 0
!$    local_size = omp_get_num_threads()
!$    local_rank = omp_get_thread_num()

      j_5 = (j_3-j_2) / local_size + 1
      j_5 = ((j_5-1)/16+1)*16
      j_4 = j_5 * local_rank
      j_5 = min(j_5, j_3-(j_2+j_4)+1)
      
      if (j_5 > 0) then

        call dtrsm('R','L','N','N',
     &       j_5, m, ONE, sm(1,1), nsm, v(j_2+j_4,1), nm)

      end if

!$OMP BARRIER
!-
!$    if ((DO_OVERLAP_BCAST_LEVEL == 0) .or.
!$   &     local_size == 1 .or.
!$   &     local_rank >= 1) then
!     
!$      if (DO_OVERLAP_BCAST_LEVEL >= 1 .and. local_size > 1) then
!$        ll_size  = local_size - 1
!$        ll_rank  = local_rank - 1
!$      else
          ll_size  = local_size
          ll_rank  = local_rank
!$      end if
!     
        if (ll_rank >= 0) then

          ii_step = 0
          do ii_0=i_2,i_3,CHUNK_N
            ii_2=ii_0; ii_3=min(ii_0+CHUNK_N-1,i_3)

            do jj_0=j_2,j_3,CHUNK_M
              jj_2=jj_0; jj_3=min(jj_0+CHUNK_M-1,j_3)

              i_5 = ii_3 - ii_2 + 1
              j_5 = jj_3 - jj_2 + 1

              if (m > 0 .and. j_5 > 0 .and. i_5 > 0) then
                if (mod(ii_step, ll_size) == ll_rank) then

                  call dgemm('N','N',
     &                 j_5, i_5, m,
     &                 ONE,  v (jj_2 ,1   ), nm,
     &                 ss(1    ,ii_2), m,
     &                 ONE,  z (jj_2 ,ii_2), nmz)

                end if
                ii_step = ii_step + 1
              end if

            end do
          end do

        end if

!$    end if

!-
      return

      end subroutine  eigen_trbakwy_block_body2

      end module trbakwy4_body_mod
