       subroutine hrbakwy_body2(local_nvec, z, nmz,
     $            beta, v, nm, m, i, ss, tt, nss, iblk,
     $            dcom, dx, dy, dz)
!$     use OMP_LIB
       use MPI
       use eigen_devel_mod
       use eigen_libs_mod
       use eigen_hrbak_mod
       implicit NONE
*
       integer, intent(in)    :: local_nvec, nmz, nm, m, i, nss, iblk
       complex(8), intent(inout) :: z(nmz,*)
       complex(8), intent(inout) :: beta(*), v(nm,*)
       complex(8), intent(inout) :: ss(*), tt(*)
       real(8), intent(inout) :: dcom, dx, dy, dz

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: i_0, m_0
       integer                :: ii
       real(8)                :: fr
       real(8)                :: ds, de
*
       integer                :: local_rank, local_size

!$     integer                :: ierr
*
       local_size = 1
       local_rank = 0
!$     local_size = omp_get_num_threads()
!$     local_rank = omp_get_thread_num()
*-
          i_2 = 1
          i_3 = local_nvec
*
          j_2 = eigen_loop_start(1,         x_nnod,x_inod)
          j_3 = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
!
*-
#if AT_BCAST_OVERLAP
!$     if ( DO_OVERLAP_BCAST_LEVEL == 2 ) then
!===
!$        call sync_other_than_master( TRBK_LOCK, TRBK_MASK )
!===
!$     endif
#endif
!$     if (
#if AT_BCAST_OVERLAP
!$   &       DO_OVERLAP_BCAST_LEVEL /= 2
#else
!$   &       .TRUE.
#endif
!$   &       ) then
!===
!$OMP     BARRIER
!===
!$     endif
*-
#if AT_BCAST_OVERLAP
          fr = TRBK_TIME_FR
          ii = MAX(8,int((i_3-i_2+1)*fr))
          ii = MIN(i_3-i_2+1,ii)

          fr = (TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x)
     $       / (TRBK_TIME_TRBK1 + TRBK_TIME_TRBK1x + TRBK_TIME_REDUC)
#else
          fr = 0D0
          ii = i_3-i_2+1
#endif
*-
! SM:= SM+lower(V^TV)
! SS(1:m,1:n):= SS+V(j_2:j_3,1:m)^T*Z(j_2:j_3,1:n)
          if ( local_size==1 .OR. local_rank==1 ) then
             ds = eigen_get_wtime()
          end if
          call hrbk1(z, nmz,
     $               v, nm, m, ss(1+ns0), ss(1),
     $               i_2,i_2+ii-1,j_2,j_3)
          if ( local_size==1 .OR. local_rank==1 ) then
             de = eigen_get_wtime()
!$           TRBK_TIME_TRBK1 = de - ds
             dx = dx + (de - ds)
          end if
*-
*-----------------------------------------------------------------
!$OMP  BARRIER
*-----------------------------------------------------------------
*-
!$OMP  MASTER
          ds = eigen_get_wtime()
          ss(ns0-5) = TRBK_TIME_REDUC
          ss(ns0-4) = TRBK_TIME_BCAST
          ss(ns0-3) = TRBK_TIME_TRBK1
          ss(ns0-2) = TRBK_TIME_TRBK1y
          ss(ns0-1) = TRBK_TIME_TRBK2
          ss(ns0-0) = fr
          call reduce_cmplx(ss(1), tt(1), ii*m+ns0,
     $         1, x_COMM_WORLD)
          TRBK_TIME_REDUC  = DBLE(ss(ns0-5)) / x_nnod
          TRBK_TIME_BCAST  = DBLE(ss(ns0-4)) / x_nnod
          TRBK_TIME_TRBK1  = DBLE(ss(ns0-3)) / x_nnod
          TRBK_TIME_TRBK1x = DBLE(ss(ns0-2)) / x_nnod
          TRBK_TIME_TRBK2  = DBLE(ss(ns0-1)) / x_nnod
          TRBK_TIME_FR     = DBLE(ss(ns0-0)) / x_nnod
          de = eigen_get_wtime()
!$        TRBK_TIME_REDUC = (de-ds)
          dcom = dcom + (de - ds)
!$OMP  END MASTER
*-
#if AT_BCAST_OVERLAP
          if ( i_3 >= i_2 + ii ) then
          if ( local_size==1 .OR. local_rank==1 ) then
             ds = eigen_get_wtime()
          end if
          call hrbk1(z, nmz,
     $               v, nm, m, ss(1+ns0), ss(1),
     $               i_2+ii,i_3,j_2,j_3)
          if ( local_size==1 .OR. local_rank==1 ) then
             de = eigen_get_wtime()
             dx = dx + (de - ds)
          end if
          end if
#endif
*-
*-----------------------------------------------------------------
!$OMP  BARRIER
*-----------------------------------------------------------------
*-
#if AT_BCAST_OVERLAP
          if ( local_size==1 .OR. local_rank==1 ) then
          if ( i_3 >= i_2 + ii ) then
!$           TRBK_TIME_TRBK1y = de - ds
          else
!$           TRBK_TIME_TRBK1y = 0D0
          end if
          end if
#endif
*-
!$OMP  MASTER
#if AT_BCAST_OVERLAP
          if ( i_3 >= i_2 + ii ) then
          ds = eigen_get_wtime()
          call reduce_cmplx(ss(1+ns0+ii*m), tt(1), (i_3-i_2+1-ii)*m,
     $         1, x_COMM_WORLD)
          de = eigen_get_wtime()
!$        TRBK_TIME_REDUC = TRBK_TIME_REDUC + (de-ds)
          dcom = dcom + (de - ds)
          end if
#endif
*-
          ds = eigen_get_wtime()
!$OMP END MASTER
!$OMP BARRIER
!====================================
!     diag(SM):= diag(SM ? beta : 1)
!$OMP  DO
      do m_0=1,m
         i_1 = m_0 + (m_0-1)*nsm
         if ( ss(i_1) == ZEROZ ) then
           ss(i_1) = ONEZ
         else
           ss(i_1) = dconjg(beta(i+m_0-1))
         end if
      end do
!$OMP  ENDDO
!====================================
!$OMP MASTER
*-
#if AT_BCAST_OVERLAP
!$        if ( TRBK_SWITCHED < 4 ) then
!$           call trbk_decide_overlap_level( i )
!$        end if
#endif
*-
          de = eigen_get_wtime()
          dz = dz + (de - ds)
!$OMP  END MASTER
*-
*-----------------------------------------------------------------
!$OMP  BARRIER
*-----------------------------------------------------------------
*-
! V:= V*(I-SM)^{-1}
! Z(j_2:j_3,1:n):= Z + V(j_2:j_3,1:m)*SS(1:m,1:n)
       if ( local_size==1 .OR. local_rank==1 ) then
          ds = eigen_get_wtime()
       end if
          call hrbk2( z,
     $                nmz, v, nm, m,
     $                ss(1+ns0), ss(1),
     $                i_2,i_3,j_2,j_3 )

       if ( local_size==1 .OR. local_rank==1 ) then
          de = eigen_get_wtime()
!$        TRBK_TIME_TRBK2 = de - ds
          dy = dy + (de - ds)
       end if
*-
#if AT_BCAST_OVERLAP
!$     if ( TRBK_SWITCHED < 4 ) then
!$     if ( DO_OVERLAP_BCAST_LEVEL == 0 ) then
!$OMP     BARRIER
!$     end if
!$     end if
#endif
*-
       return
       end subroutine  hrbakwy_body2
!===
#if AT_BCAST_OVERLAP
!===
!$     subroutine trbk_decide_overlap_level( i )
!$     use OMP_LIB
!$     use MPI
!$     use eigen_devel
!$     use eigen_libs
!$     use eigen_hrbak
!$     implicit NONE
!$     integer, intent(in)    :: i
!$     integer                :: ll0, ll1, ll2, ll_(6)
!$     real(8)                :: BCAST_TIME(0:2)
!$     real(8)                :: f0, f1
!$     integer                :: local_rank, local_size
!$     integer                :: ierr
*-
!$        local_size = omp_get_num_threads()
!$        local_rank = omp_get_thread_num()
!$
!$        if ( local_size > 1 ) then
!$
!$        TRBK_TIME_COUNTER = TRBK_TIME_COUNTER + 1
!$        if ( TRBK_TIME_COUNTER >= TRBK_TIME_NEXT ) then
!$           f0 = 1D0 * (local_size-1) / local_size
!$           f1 = 1D0 / f0
!$           if ( DO_OVERLAP_BCAST_LEVEL == 2 ) then
!$              BCAST_TIME(0) = TRBK_TIME_BCAST 
!$   $                      + TRBK_TIME_TRBK2 * f0
!$   $                      + TRBK_TIME_TRBK1 * f0
!$              BCAST_TIME(1) = MAX( TRBK_TIME_BCAST,
!$   $                        TRBK_TIME_TRBK2 )
!$   $                      + TRBK_TIME_TRBK1 * f0
!$              BCAST_TIME(2) = MAX( TRBK_TIME_BCAST,
!$   $                        TRBK_TIME_TRBK2 
!$   $                      + TRBK_TIME_TRBK1 )
!$           endif
!$           if ( DO_OVERLAP_BCAST_LEVEL == 1 ) then
!$              BCAST_TIME(0) = TRBK_TIME_BCAST 
!$   $                      + TRBK_TIME_TRBK2 * f0
!$   $                      + TRBK_TIME_TRBK1
!$              BCAST_TIME(1) = MAX( TRBK_TIME_BCAST,
!$   $                        TRBK_TIME_TRBK2 )
!$   $                      + TRBK_TIME_TRBK1
!$              BCAST_TIME(2) = MAX( TRBK_TIME_BCAST,
!$   $                        TRBK_TIME_TRBK2
!$   $                      + TRBK_TIME_TRBK1 * f1)
!$           endif
!$           if ( DO_OVERLAP_BCAST_LEVEL == 0 ) then
!$              BCAST_TIME(0) = TRBK_TIME_BCAST
!$   $                      + TRBK_TIME_TRBK2
!$   $                      + TRBK_TIME_TRBK1
!$              BCAST_TIME(1) = MAX( TRBK_TIME_BCAST,
!$   $                        TRBK_TIME_TRBK2 * f1 )
!$   $                      + TRBK_TIME_TRBK1
!$              BCAST_TIME(2) = MAX( TRBK_TIME_BCAST,
!$   $                        TRBK_TIME_TRBK2 * f1
!$   $                      + TRBK_TIME_TRBK1 * f1 )
!$           end if
!$           ll0 = DO_OVERLAP_BCAST_LEVEL
!$           ll1 = 0; do ll2 = 1, 2
!$              if ( BCAST_TIME(ll1) >= BCAST_TIME(ll2) ) then
!$                 ll1 = ll2
!$              end if
!$           end do
!$           ll_(1:3)=0; ll_(ll1+1)=1
!$           call MPI_Allreduce( ll_(1), ll_(4),
!$   $                3, MPI_INTEGER, MPI_SUM, y_COMM_WORLD, ierr )
#if OVERLAP_DECISION_TYPE==3
!               type 3 mean value
!$           ll1 = 0; do ll2 = 1, 2
!$              if ( ll_(4+ll2) > ll_(4+ll1) ) then
!$                 ll1 = ll2
!$              end if
!$           end do
#endif
#if OVERLAP_DECISION_TYPE==2
!               type 2 max value
!$           do ll2 = 0, 2
!$              if ( ll_(4+ll2) > 0 ) then
!$                 ll1 = ll2
!$              end if
!$           end do
#endif
#if OVERLAP_DECISION_TYPE==1
!               type 1 min value
!$           do ll2 = 2, 0, -1
!$              if ( ll_(4+ll2) > 0 ) then
!$                 ll1 = ll2
!$              end if
!$           end do
#endif
#if _DEBUG_
!$           if ( TRD_inod == 1 ) then
!$                 print*, i, "Overlap decision", ll0,"->",ll1,
!$   $                    BCAST_TIME(0:2)
!$           end if
#endif
!$           if ( ll1 == 0 ) then
!$              TRBK_SWITCHED = TRBK_SWITCHED + 1
!$           end if
!$           if ( ll0 .NE. ll1 ) then
!$              if ( TRBK_SWITCHED >= 8 ) then
!$                 TRBK_TIME_INTERVAL = 32
!$              else
!$                 TRBK_TIME_INTERVAL = 1
!$              end if
!$           else
!$              TRBK_TIME_INTERVAL = TRBK_TIME_INTERVAL*2
!$           end if
!$           TRBK_TIME_NEXT = TRBK_TIME_NEXT + TRBK_TIME_INTERVAL
!$           DO_OVERLAP_BCAST_LEVEL = ll1
!$        end if
!$        else
!$           DO_OVERLAP_BCAST_LEVEL = 0
!$        end if
*-
!$     end subroutine trbk_decide_overlap_level
!===
#endif
*-
       subroutine hrbk1(z, nmz,
     $            v, nm, m, ss, sm, i_2,i_3,j_2,j_3)
!$     use OMP_LIB
       use MPI
       use eigen_devel_mod
       use eigen_libs_mod
       use eigen_hrbak_mod
       implicit NONE
*
       integer, intent(in)    :: nmz, nm, m
       integer, intent(in)    :: i_2, i_3, j_2, j_3
       complex(8), intent(inout) :: z(nmz,*),v(nm,*)
       complex(8), intent(inout) :: ss(m,*)
*
       complex(8), intent(inout) :: sm(nsm,nsm)
*
       integer                :: i_1, i_4, i_5
       integer                :: j_1, j_4, j_5
       integer                :: i_0, j_0, m_0
       integer                :: ii_2, ii_3, ii_4
       integer                :: jj_2, jj_3, jj_4
*
       real(8)                :: s0, s1, s2, s3
       real(8)                :: v0, v1, v2, v3
       real(8)                :: t0
*
       integer                :: local_rank, local_size
       integer                :: ll_rank, ll_size
*
       integer                :: ii_step, blk_size1, blk_size2
       integer                :: BLAS_CHUNK1 = 64
       integer                :: BLAS_CHUNK2 = 64
*

       local_size = 1
       local_rank = 0
!$     local_size = omp_get_num_threads()
!$     local_rank = omp_get_thread_num()
*-
!$     if ( (DO_OVERLAP_BCAST_LEVEL <= 1 .AND. i_2 == 1) .OR.
!$   &      local_size == 1 .OR.
!$   &      local_rank >= 1 ) then
!$        if (
!$   &       (DO_OVERLAP_BCAST_LEVEL == 2 .AND. local_size > 1) .OR.
!$   &       (i_2 > 1 .AND. local_size > 1) ) then
!$         ll_size  = local_size - 1
!$         ll_rank  = local_rank - 1
!$      else
           ll_size  = local_size
           ll_rank  = local_rank
!$      endif
*-
       if ( i_2 == 1 ) then
*-
          j_5 = j_3 - j_2 + 1
          if ( j_5 > 0 ) then

             ii_step = 0
             do i_1 = 1, m, BLAS_CHUNK1
             do j_1 = 1, m, BLAS_CHUNK2

                ii_3 = MIN(m, i_1+BLAS_CHUNK1-1)
                jj_3 = MIN(m, j_1+BLAS_CHUNK2-1)

                blk_size1 = ii_3 - i_1 + 1
                blk_size2 = jj_3 - j_1 + 1

                if ( blk_size1 > 0 .AND. blk_size2 > 0 ) then
                if ( i_1+blk_size1 >= j_1+blk_size2 ) then

                   ii_2 = i_1
!                   ii_2 = MIN(MAX(i_1,j_1),ii_3)
!                   blk_size1 = ii_3 - ii_2 + 1

                   if ( MOD(ii_step, ll_size) == ll_rank ) then

                      call ZGEMM('C','N',
     &                     blk_size1, blk_size2, j_5,
     &                     MONEZ, v(j_2, ii_2), nm,
     &                            v(j_2, j_1),  nm,
     &                     ZEROZ, sm(ii_2,j_1), nsm)

                   end if

                   ii_step = ii_step + 1

                end if
                end if

             end do
             end do

          end if
*-
        end if
*-
#if defined(__INTEL_COMPILER)
        j_5 = j_3 - j_2 + 1

        do ii_2 = i_2, i_3, (1024*ll_size)

          ii_3 = MIN(ii_2+(1024*ll_size)-1,i_3)
          ii_4 = ii_2-i_2+1

          i_0 = (ii_3-ii_2) / ll_size + 1
          i_4 = i_0 * ll_rank
          i_5 = MIN(i_0, ii_3-(ii_2+i_4)+1)

          if ( m > 0 .AND. j_5 > 0 .AND. i_5 > 0 ) then

             CALL ZGEMM('C','N',
     &               m, i_5, j_5,
     &               ONEZ, v (j_2    ,1       ), nm,
     &                     z (j_2    ,ii_2+i_4), nmz,
     &               ZEROZ, ss(1      ,ii_2+i_4), m)

          endif

        end do
#else
          j_5 = j_3 - j_2 + 1

          i_5 = (i_3-i_2) / ll_size + 1
          i_4 = i_5 * ll_rank
          i_5 = MIN(i_5, i_3-(i_2+i_4)+1)

          if ( m > 0 .AND. j_5 > 0 .AND. i_5 > 0 ) then

             CALL ZGEMM('C','N',
     &               m, i_5, j_5,
     &               ONEZ, v (j_2 ,1      ), nm,
     &                     z (j_2 ,i_2+i_4), nmz,
     &               ZEROZ, ss(1   ,i_2+i_4), m)

          endif
#endif
!$     end if
*-
       return
       end subroutine  hrbk1
*-
       subroutine hrbk2( z,
     $            nmz,
     $            v, nm, m, ss, sm, i_2,i_3, j_2,j_3 )
!$     use OMP_LIB
       use MPI
       use eigen_devel_mod
       use eigen_libs_mod
       use eigen_hrbak_mod
       implicit NONE
*
       integer, intent(in)    :: nmz, nm, m
       integer, intent(in)    :: i_2, i_3, j_2, j_3
       complex(8), intent(inout) :: z(nmz,*)
       complex(8), intent(inout) :: v(nm,*)
       complex(8), intent(inout) :: ss(m,*)
*
       complex(8), intent(inout) :: sm(nsm,*)

       integer                :: i_4, i_5
       integer                :: j_1, j_4, j_5, j_6, j_7, j_8
       integer                :: i_0, j_0, m_0
*
       real(8)                :: v0, v1, v2
       real(8)                :: u0, u1
*
       integer                :: local_rank, local_size
       integer                :: ll_rank, ll_size
*

       local_size = 1
       local_rank = 0
!$     local_size = omp_get_num_threads()
!$     local_rank = omp_get_thread_num()
*
*-
          j_5 = (j_3-j_2) / local_size + 1
          j_5 = ((j_5-1)/2+1)*2
          j_4 = j_5 * local_rank
          j_5 = MIN(j_5, j_3-(j_2+j_4)+1)
          j_6 = j_2+j_4
          j_7 = j_6+j_5-1
         
          if ( j_5 > 0 ) then
          call ZTRSM('R','L','N','N',
     &               j_5, m, ONEZ, sm(1,1), nsm, v(j_6,1), nm)
          end if

!$OMP BARRIER
*-
!$     if ( (DO_OVERLAP_BCAST_LEVEL == 0) .OR.
!$   &      local_size == 1 .OR.
!$   &      local_rank >= 1 ) then

!$        if ( DO_OVERLAP_BCAST_LEVEL >= 1 .AND. local_size > 1 ) then
!$           ll_size  = local_size - 1
!$           ll_rank  = local_rank - 1
!$        else
             ll_size  = local_size
             ll_rank  = local_rank
!$        endif

          j_5 = j_3 - j_2 + 1

          i_5 = (i_3-i_2) / ll_size + 1
          i_4 = i_5 * ll_rank
          i_5 = MIN(i_5, i_3-(i_2+i_4)+1)

          if ( m > 0 .AND. j_5 > 0 .AND. i_5 > 0 ) then

             CALL ZGEMM('N','N',
     &               j_5, i_5, m,
     &               ONEZ, v (j_2    ,1      ), nm,
     &                     ss(1      ,1  +i_4), m,
     &               ONEZ, z (j_2    ,i_2+i_4), nmz)

          end if

!$     end if
*-

       return
       end subroutine  hrbk2

